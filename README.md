BarPep Scripts version 2.0 utilizes snakemake and cutadapt to allow for a faster and more accurate detection of barcode or peptide sequences from NGS data.

This repository provides an updated combined version of the [Barcode Detection Script and the Peptide Detection Extraction Script](https://github.com/grimmlabs/AAV_GrimmLab_JoVE2022) for the analysis of illumina data obtained from the screen of a diversified AAV capsid library.

# Table of Contents

[Installation](#installation)

[Quick Tutorial](#quick-tutorial)

[Configuration](#configuration)

[Barcode analysis output files explanation](#barcode-analysis-output-files-explanation)

# Installation

Easiest way is to install via the provided conda environment:

```
conda env create -f env/barpep.yml
```

# Quick Tutorial

This command will start the analysis of the provided example data.

```
conda activate barpep
snakemake --cores all
```

# Configuration

The settings are found in the config file under 'config/config.yaml' and should be adjusted for your analysis.

## general

**input directory**: path to the directory containing all fastq files.

**annotation file**: path to a tab-separated look-up table containing barcode sequences and associated names. Please adhere to the structure of example_data/seq_annotation.txt. Required for barcode analysis. Set to False if not required.

**output directory**: path to the desired output directory.

**reverse complement output**: dictates wether the found barcodes are put out as reverse complement, necessary for correct assignment using the barcode annotation table.

## cutadapt options

**flanks**: dictates the flanking sequences within the reads used to detect the barcode/peptide. Always in the format NNN...NNN, where NNN are the flanking nucleotides downstream and upstream of the barcode/peptide. The length of each flank is up to you, but 8 to 12 bp are recommended.

**error rate**: the allowed error rate for flanking region detection directly used by cutadapt.

**cutadapt cores**: the number of cores assigned to cutadapt for flank detection. Set to 0 for auto detection of maximum cores. If you do that, do not set the snakemake --cores argument too high.

**barcode length min**: minimal length of barcode/peptide that is accepted.

**barcode length max**: maximum length of barcode/peptide that is accepted.

## barcode analysis

**barcode_analysis**: Flag whether barcode analysis is performed.

**tissue_annotation**: Required tab-separated file, containing metadata on each input fastq file. Please use the scheme below. **The Sample should contain the basename of the input file without path or file extensions!**

| Sample   | SampleType| Animal     |Tissue        | weight_variable|
|----------|-----------|------------|--------------|----------------|
| Sample1  | cDNA      | M1         | Heart        | 0.00635        |
| Sample2  | gDNA      | M1         | Heart        | 0.00635        |
| Sample3  | cDNA      | M1         | Kidney       | 0.000293       |
| Sample4  | gDNA      | M1         | Kidney       | 0.000293       |
| Sample5  | cDNA      | M2         | Lung         | 0.00871        |
| Sample6  | gDNA      | M2         | Lung         | 0.00871        |

**input_basename**: Basename (no file extension or path) of the file in the input fastq directory that contains normalisation information (e.g. library before injection).

**alternative_input**: If required, you can supply the path to an alternative input variantCount.csv file. When there is a path here, it overwrites the input_basename option. Default=False.

**pseudo_count**: Pseudo count added dividing 0 values for input normalisation.

# Barcode analysis output files explanation

- **01.readCounts.csv**  
Table of raw read counts of each variant in each tissue for each unique combination of animal and sample type.

- **02.Pab.csv**  
Proportional read count values, or P<sub>αβ</sub> values. They are calculated by normalizing the read counts R of all variants α in tissue β to the sum of all variants α in β:

$$
\displaystyle
\ P_{αβ}= \frac{R_{αβ}}{\sum_{α} R_{αβ}}
$$  

- **03.Pabs.csv**  
Proportional count values normalized to the input library, or P*<sub>αβ</sub> values. They are calculated by normalizing P<sub>αβ</sub> to the proportion of each variant α in the initial library L<sub>α</sub>, thus correcting for the uneven composition in library:

$$
\displaystyle
\ P_{αβ}^*= \frac{P_{αβ}}{L_{α}}
$$

- **04.Bab.csv**  
P*<sub>αβ</sub> is weighted by the weight_variable (e.g. vg/dg or RQ values), termed G<sub>β</sub>, to allow a comparison of one variant α over all analyzed tissues β:

$$
\displaystyle
\ B_{αβ}= \frac{P_{αβ}}{L_{α}}*G_β
$$


## Currently not supported

- **V<sub>αβ</sub>**  
B<sub>αβ</sub> values are shown as proportions of the sum over all variants α of B<sub>αβ</sub>. These values can be useful to create bar plots which demonstrate the proportion of all variants α in one tissue β, exemplifying the efficiency of the individual vectors:  

$$
\displaystyle
\ V_{αβ}= \frac{B_{αβ}}{\sum_{α} B_{αβ}}
$$    

- **T<sub>αβ</sub>**  
B<sub>αβ</sub> values are shown as proportions of the sum over all tissues β of B<sub>αβ</sub>. These values can be useful to create bar plots which show the proportion of one variant α in all tissues β, allowing an analysis of the tissue specificity:

$$
\displaystyle
\ T_{αβ}= \frac{B_{αβ}}{\sum_{β} B_{αβ}}
$$  
