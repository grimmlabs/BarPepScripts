import os, shutil, warnings
import pandas as pd

configfile: "config/config.yaml"
input_directory = config["input_directory"]
output_directory = config["output_directory"]

INPUT_FILES = [
    os.path.join(input_directory, f) for f in os.listdir(input_directory)
    ]
INPUT_FILES_ENDING = ".".join(INPUT_FILES[0].split(".")[1:])
SAMPLES = [
    os.path.basename(s).split(".")[0] for s in INPUT_FILES
    ]

all_input = [
    expand("{output_directory}/logs/{sample}.recoveredReads.log", sample=SAMPLES, output_directory = output_directory)
]

# check necessary input for barcode analysis
if config["barcode_analysis"] and not config["annotation_file"]:
 raise Exception("Barcode analysis is not possible when no annotation file is given. Please adjust config.")
elif config["barcode_analysis"] and isinstance(config["annotation_file"], str):
    all_input.append(output_directory + "/BC_analysis/04.Bab.csv")

# check flanking region length
flanks = config["flanks"].split("...")
if len(flanks[0]) < 10 or len(flanks[1]) < 10:
    print("\033[35mWarning: Please consider longer flanking regions of at least 10 bps each.\033[0m")

rule all:
    input:
        all_input,
        output_directory + "/logs/used_config.yaml"


rule find_barcodes:
    input:
        input_directory + "/{sample}." + INPUT_FILES_ENDING
    output:
        temp(output_directory + "/{sample}.barcodes.fasta")
    params:
        flanks = config["flanks"],
        error_rate = config["error_rate"],
        cores = config["cutadapt_cores"]
    log:
        output_directory + "/logs/{sample}.cutadapt.log"
    shell:
        "cutadapt -a '{params.flanks}' --revcomp -e {params.error_rate} -j {params.cores} -o {output} {input} > {log}"


SEQKIT_FLAGS = "-sprgv" if config["reverse_complement_output"] else "-sgv"
rule count_barcodes:
    input:
        output_directory + "/{sample}.barcodes.fasta"
    output:
        temp(output_directory + "/{sample}.counts.csv")
    params:
        min_length = config["barcode_length_min"],
        max_length = config["barcode_length_max"]
    shell:
        "seqkit seq {SEQKIT_FLAGS} -t DNA -m {params.min_length} -M {params.max_length} {input} | sort | uniq -c | awk '{{print $1, \",\", $2}}' > {output}"


rule join_with_variants:
    input:
        counts=output_directory + "/{sample}.counts.csv"
    output:
        output_directory + "/{sample}.variantCounts.csv"
    params:
        variants=config["annotation_file"]
    run:
        cuad = pd.read_csv(input.counts, names=["freq", "barcode"])
        cuad["barcode"] = cuad["barcode"].str.strip().str.upper()
        if params.variants == False:
            cuad.to_csv(output[0], index=False)
        else:
            var = pd.read_csv(params.variants, names=["barcode", "variants"], delimiter = "\t")
            var["barcode"] = var["barcode"].str.strip().str.upper()
            joined = cuad.merge(var, how = "outer", on = "barcode").fillna(value = {"freq": 0, "variants": "Unknown"}).astype({'freq': 'int64'})
            joined.to_csv(output[0], index=False)

rule gather_input_stats:
    input:
        input_directory + "/{sample}." + INPUT_FILES_ENDING
    output:
        temp(output_directory + "/{sample}.stats.txt")
    shell:
        "seqkit stats {input} > {output}"


rule make_recoveredReads_log:
    input:
        stats=output_directory + "/{sample}.stats.txt",
        variantCounts=output_directory + "/{sample}.variantCounts.csv"
    output:
        output_directory + "/logs/{sample}.recoveredReads.log"
    run:
        # get number of reads, of variants
        with open(input.stats, "r") as fi: total_reads = int(fi.readlines()[1].split()[3].replace(",", ""))
        variantCounts = pd.read_csv(input.variantCounts)
        recoveredReads = variantCounts["freq"].sum()
        recoveredPerc = round(recoveredReads / total_reads*100, 2)
        negativeReads = total_reads - recoveredReads
        negativePerc = round(negativeReads / total_reads*100, 2)
        if config["annotation_file"] == False: # limited information when no annotation file is given
            with open(output[0], "w") as out:
                out.write(wildcards.sample + "\n\n")
                out.write(f"Total number of reads: {total_reads}\n")
                out.write(f"Reads recovered: {recoveredReads} ({recoveredPerc}%)\n\n")
                out.write(f"Reads with no detected flanking regions: {negativeReads} ({negativePerc}%)\n")
        else:
            # get number of unknown variants, and
            unknownVariants = variantCounts["variants"].value_counts()["Unknown"]
            unknownPerc = round(unknownVariants / total_reads*100, 2)
            with open(output[0], "w") as out:
                out.write(wildcards.sample + "\n\n")
                out.write(f"Total number of reads: {total_reads}\n")
                out.write(f"Reads recovered: {recoveredReads} ({recoveredPerc}%)\n\n")
                out.write(f"Reads with unknown variants: {unknownVariants} ({unknownPerc}%)\n")
                out.write(f"Reads with no detected flanking regions: {negativeReads} ({negativePerc}%)\n")
                out.write(f"Reads with contaminating variants: 0 (0.0%)\n") # TODO




rule BC_analysis:
    input:
        expand("{output_directory}/{sample}.variantCounts.csv", sample=SAMPLES, output_directory = output_directory)
    output:
        output_directory + "/BC_analysis/04.Bab.csv"
    params:
        tissue_annot = config["tissue_annotation"]
    run:
        # get input
        if not config["alternative_input"]:
            INPUT_SAMPLE = config["input_basename"]
            INPUT_variantCounts = pd.read_csv(f"{output_directory}/{INPUT_SAMPLE}.variantCounts.csv")
        elif not os.path.isfile(config["alternative_input"]):
            try:
                INPUT_variantCounts = pd.read_csv(config["alternative_input"])
            except:
                raise FileNotFoundError('Your supplied alternative input file was not found')
        i = INPUT_variantCounts.loc[INPUT_variantCounts.variants != "Unknown"].set_index("variants").drop("barcode", axis=1).rename({"freq": "input_norm"}, axis=1) + float(config["pseudo_count"])
        input_norm = i / i.sum(axis=0)

        # load data
        variantCounts = pd.concat(pd.read_csv(f).assign(Sample=os.path.basename(f).split(".")[0]) for f in [f"{output_directory}/{s}.variantCounts.csv" for s in SAMPLES]).reset_index(drop=True)
        tissue_annot = pd.read_csv(params.tissue_annot, delimiter = "\t")
        tissue_weights = pd.Series(tissue_annot.weight_variable.values, index=tissue_annot.Tissue).to_dict()
        master_df = pd.merge(variantCounts, tissue_annot, "outer", "Sample")
        # read counts
        tissue_ReadCounts = (master_df[["freq", "variants", "SampleType", "Animal", "Tissue"]]
                                .loc[master_df.variants != "Unknown"].dropna()
                                .pivot(index="variants", values="freq", columns=["SampleType", "Animal", "Tissue"]))
        tissue_ReadCounts.to_csv(f"{output_directory}/BC_analysis/01.readCounts.csv")
        # Pab
        tissue_Pab = tissue_ReadCounts / tissue_ReadCounts.sum(axis=0)
        tissue_Pab.to_csv(f"{output_directory}/BC_analysis/02.Pab.csv")

        # Pab*
        Pabs = pd.merge(tissue_ReadCounts.stack(level = ["SampleType", "Animal"]).reset_index().set_index("variants"), input_norm, how="left", left_index=True, right_index=True)
        tissue_list = [t for t in tissue_weights.keys()]
        Pabs[tissue_list] = Pabs[tissue_list] + float(config["pseudo_count"]) # add pseudocount
        Pabs[tissue_list] = Pabs[tissue_list].div(Pabs.input_norm, axis=0)
        Pabs.to_csv(f"{output_directory}/BC_analysis/03.Pabs.csv")
        # Bab
        Bab = Pabs.copy()
        for tissue, norm_val in tissue_weights.items():
            Bab[tissue] = Bab[tissue] * norm_val
        Bab.to_csv(f"{output_directory}/BC_analysis/04.Bab.csv")


rule save_config_provenance:
    output:
        output_directory + "/logs/used_config.yaml"
    run:
        # copy original
        shutil.copy(workflow.configfiles[0], output[0])
