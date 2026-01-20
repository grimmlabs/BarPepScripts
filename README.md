BarPep Scripts version 2.0 utilizes snakemake and cutadapt to allow for a faster and more accurate detection of barcode or peptide sequences from NGS data.

This repository provides an updated combined version of the [Barcode Detection Script and the Peptide Detection Extraction Script](https://github.com/grimmlabs/AAV_GrimmLab_JoVE2022) for the analysis of illumina data obtained from the screen of a diversified AAV capsid library.

# Table of Contents

[Installation](#installation)

[Quick Tutorial](#quick-tutorial)

[Configuration](#configuration)

# Installation

Easiest way is to install via the provided conda environment:

```
conda env create -f env/barpep.yaml
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

*input directory*: path to the directory containing all fastq files.

*barcode annotation file*: path to a tab-separated look-up table containing barcode sequences and associated names.

*output directory*: path to the desired output directory.

*reverse complement output*: dictates wether the found barcodes are put out as reverse complement, necessary for correct assignment using the barcode annotation table.

## cutadapt options

*flanks*: dictates the flanking sequences within the reads used to detect the barcode/peptide. Always in the format NNN...NNN, where NNN are the flanking nucleotides downstream and upstream of the barcode/peptide. The length of each flank is up to you, but 8 to 12 bp are recommended.

*error rate*: the allowed error rate for flanking region detection directly used by cutadapt.

*cutadapt cores*: the number of cores assigned to cutadapt for flank detection. Set to 0 for auto detection of maximum cores. If you do that, do not set the snakemake --cores argument too high.

*barcode length min*: minimal length of barcode/peptide that is accepted.

*barcode length max*: maximum length of barcode/peptide that is accepted.
