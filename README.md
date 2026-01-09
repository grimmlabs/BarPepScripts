<div style="text-align: justify">

BarPep Scripts version 2.0 utilizes snakemake and cutadapt to allow for a faster and more accurate detection of barcode or peptide sequences from NGS data.

This repository provides an updated combined version of the [Barcode Detection Script and the Peptide Detection Extraction Script](https://github.com/grimmlabs/AAV_GrimmLab_JoVE2022) for the analysis of illumina data obtained from the screen of a diversified AAV capsid library.
&emsp;

# Table of Contents

[Installation](#installation)

[Quick Tutorial](#quick-tutorial)


# Installation

```
conda create -f env/barpep.yaml
```

# Quick Tutorial

This command will start the analysis of the provided example data.

```
conda activate barpep
snakemake --cores 4
```

The settings are found in the config file under 'config/snakemake.config.yaml'.

_input directory_: path to the directory containing all fastq files

_variants file_: path to the file containig named barcodes

_output directory_: path to the desired output directory

_reverse complement output_: dictates wether the found barcodes put out as reverse complement

_flanks_: dictates the flanking sequences used to detect the barcode/peptide always in the format NNN...NNN. The length of each flank is up to you, but 8 to 12 bp are recommended.

_error rate_: the allowed error rate for flanking region detection
_cutadapt cores_: the number of cores assigned to cutadapt for flank detection. Set to 0 for auto detection of maximum cores. If you do that, do not set the snakemake --cores argument too high.
_barcode length min_: minimal length of barcode/peptide that is accepted
_barcode length max_: maximum length of barcode/peptide that is accepted
