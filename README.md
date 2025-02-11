<div style="text-align: justify">

This repository provides an updated combined version of the [Barcode Detection Script and the Peptide Detection Extraction Script](https://github.com/grimmlabs/AAV_GrimmLab_JoVE2022) for the analysis of illumina data obtained from the screen of a diversified AAV capsid library. This README will give you a detailed instruction on how to use the updated script, i.e., how your input files should look like, what arguments you need to specifiy, and what kind of output files will be generated. 

&emsp;

# Table of Contents

[Quick Tutorial](#quick-tutorial)

[Requirements](#requirements)

A. [BarPepDetection](#barpepdetection)
1. [Prepare your Input Files](#prepare-your-input-files)
    - [For Barcode Detection](#for-barcode-detection)
    - [For Peptide Detection](#for-peptide-detection)
2. [Arguments BarPepDetection](#arguments-barpepdetection)
3. [Output Files BarPepDetection](#output-files-barpepdetection)
    - [For Barcode Detection](#for-barcode-detection-1)
    - [For Peptide Detection](#for-peptide-detection-1)

B. [BarPepAnalysis](#barpepanalysis)
1. [Prepare your Input Files](#prepare-your-input-files)
    - [For Barcode Analysis](#for-barcode-analysis)
    - [For Peptide Analysis](#for-peptide-analysis)
2. [Arguments BarPepAnalysis](#arguments-barpepanalysis)
3. [Output Files BarPepAnalysis](#output-files-barpepanalysis)
    - [For Barcode Analysis](#for-barcode-analysis-1)
    - [For Peptide Analysis](#for-peptide-analysis-1) 

&emsp;

# Quick Tutorial

These commands wil get you going in running the scripts on the provided barcoded example data set.

```
# Detection
mkdir example_output
python BarPepDetection.py -a BC -v example_data/variants.txt -d example_data/raw_reads/ -o example_output/ -l GGCCCA -r CCAGCC

# Analysis
python BarPepAnalysis.py -a BC -i example_data/assignment_file_for_analysis.csv -l example_output/m1_input_Variants.csv -d example_output/ -x
```
# Requirements

Python3 (version used to test this script: 3.12.3) with these modules is required for running:

- Biopython (1.83)
- Numpy (1.26.4)
- Pandas (2.2.2)
- Matplotlib (3.8.4)
- Seaborn (0.13.2)

**Note:** Especially for the pandas module, it is important that you use the indicated version! The script might not work with older or newer versions of pandas.

To install required modules :
```
# on linux
pip install -r requirements.txt

# on windows
py -m pip install -r requirements.txt
```
Read more about pip : https://pip.pypa.io/en/stable/installation/
&emsp;

# BarPepDetection

In principle the detection of barcode or peptide variants works by identifying user-supplied flanking regions within every input read in fastq format. The output consists of the count of each barcode or peptide variant for every input file.

## Prepare your Input Files
### _For Barcode Detection_
If you want to use the script to detect barcode sequences, you will need the following input files:
- A folder containing **only** the raw input files in fastq format (either gunzipped or not). It would be best if you already rename the files and give them a distinct name, e.g. _Sample1.txt.gz_, _Sample2.txt.gz_ or _M1_heart_cDNA.txt.gz_, _M2_liver_gDNA.txt.gz_ etc.
- A tab-delimited text file "Variants.txt" that includes your barcode sequences and the respective variants, e.g:

    > AGACTCGTTGTATAT&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;&nbsp;AAV1  
    TGGGCGGTCAGGGTC&emsp;&emsp;&emsp;&emsp;AAV2  
    TTGCCGTCCTTCGAG&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;AAV3  
    TTCAGCGGACGGGCC&emsp;&emsp;&emsp;&emsp;AAV4

- If you want to check your sequencing data for possible contaminations (i.e., other barcode sequnces that are used frequently in your lab), then you should also prepare **a tab-delimited text file "Contaminations.txt"** that includes the barcode sequences used in your lab and the respective variants, e.g.:

    > AGACTCGTTGTATAT&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;&nbsp;A2  
    TAGAGATTTAAACCG&emsp;&emsp;&emsp;&emsp;&nbsp;A3  
    CGTGACAGCGGATGG&emsp;&emsp;&emsp;&emsp;A4  
    TGGGCGGTCAGGGTC&emsp;&emsp;&emsp;&emsp;A5

 &emsp;

### _For Peptide Detection_
If you want to use the script to detect peptide insertions, you will need the following input:
- A folder containing **only** the raw input files in fastq format (either gunzipped or not). It would be best if you already rename the files and give them a distinct name, e.g. _Sample1.txt.gz_, _Sample2.txt.gz_ or _M1_heart_cDNA.txt.gz_, _M2_liver_gDNA.txt.gz_ etc.

 &emsp;

## Arguments BarPepDetection

 &emsp;

**<p style="font-size:15px;">Required arguments:</p>**
- **-a MODE:**  
Specify in which mode you want to run the script. If you want to use it for detection of barcodes, type "BC", if you want to use it for detection of peptide insertions, type "PV".

- **-d DIRECTORY:**  
Give the path to the directory where your sequencing data (gz files) are located.

- **-o OUTPUTDIR:**  
Give the path to the directory where the output files should be saved.

- **-l BCVLEFT:**  
Give the short flanking oligo sequence at 5' of the barcode/peptide location.  

- **-r BCVRIGHT:**  
Give the short flanking oligo sequence at 3' of the barcode/peptide location.

**Important:** The orientation in which you indicate the flanking regions is crucial in determining the correct orientation of the found barcode or peptide sequence. If you run into problems it is advisable to double check your input reads if the flanking sequences appear in the expected orientation.

&emsp;

**<p style="font-size:15px;">Optional arguments:</p>**

- **-n BCVLOC:**  
Give the position of the first expected barcode/peptide nucleotide if the read numbering starts with 0.
- **-m BCVMARGIN:**  
Give the number of nucleotides before and after BCV_loc to search for the barcode/peptide.
- **-k BCVLOCREVCOMP:**  
Give the position of the first expected barcode/peptide nucleotide on the reverse complement strand if the read numbering starts with 0.

- **-p PLOTS:**  
Set this flag if you want to generate plots showing the quality of your sequencing run which will be saved in a seperate pdf file for each sample. **Note:** The plots will only include a subsample of the data (100,000 reads).  

- **-w SILENCE:**  
Set this flag if you want to avoid printouts in the terminal.  

- **-z VERSION:**  
Set this flag if you want to show the script's version number and exit.

&emsp;

**<p style="font-size:15px;">Additional arguments for barcode detection:</p>**

- **-v VARIANTS:**  
Give the path to the tab-delimited text file that includes unique barcode sequences assigned to one of the cap variants. This argument is required for barcode detection. 

- **-c CONTAMINATIONS:**  
If you want to check your sequencing data for contaminations, give the path to the tab-delimited text file that includes unique barcode sequences assigned to contaminating cap variants.  

&emsp;

**<p style="font-size:15px;">Additional arguments for peptide detection:</p>**

- **-s BCVSIZE:**  
Give the length of the peptide sequence in nucleotides. This argument is required for peptide detection.  

&emsp;

By default, the script searches for the barcode or peptide sequence over the complete length of the read. If you want to restrict the search to a specific margin, you need to specifiy the arguments BCVLOC, BCVMARGIN, and BCVLOCREVCOMP. These arguments restrict the search for the flanking regions to a given area within all reads, which improves the performance in longer reads.

&emsp;

## Output Files BarPepDetection
### _For Barcode Detection_
Running the script in barcode detection mode will generate at least 3 different output files per gz input file.

- **Log file**:  
This text file summarizes all the important numbers: It includes the total number of reads as well as the number and proportion of recovered reads (i.e., reads in which the flanking constant regions were found), of reads with barcodes corresponding to your expected variants, of reads with barcodes corresponding to contaminating variants (if a "Contaminations.txt" file was given, otherwise this number will be 0), of reads with unknown variants, and of reads where the constant regions were not found. Additionally, you can find the mean sequence length and mean sequence quality here.

- **Variants file:**  
This CSV file can be loaded into Excel, giving you an alphabetically sorted table containing the expected variants with their barcode sequence and the frequency with which they were found in your sequencing data.

- **Unknown  Variants file**:  
This CSV file can be loaded into Excel, giving you a table containing the contaminating or unknown barcode sequences sorted by their frequency in the sequencing data. If you ran the script with a contamination file as input, the contaminating barcodes will be assigned to a variant name. Otherwise, they will all be labelled "Unknown".  

&emsp;

Optional:

- **Plots file**:  
This file will only be generated if you set the flag for the plot argument. The pdf file includes different plots that visualize the quality of your sequencing run: a box plot showing the per base sequence quality, a histogram showing the per sequence quality, a line plot showing the sequence content across all bases, and a histogram showing the sequence length distribution. One look at these plots should give you a rough overview over your sequencing quality.  

&emsp;

### _For Peptide Detection_
Running the script in peptide detection mode will generate at least 2 different output files per sample input file.

- **Log file:**  
Similarly to the barcode detection mode, a log file will be generated. This text file includes the total number of reads and of recovered reads, as well as the mean sequence length and quality.

- **PV file:**  
This CSV file can be loaded into Excel, giving you a table containing the detected peptide DNA sequences sorted by their frequency in the sequencing data.

Additionally, if you set the flag for the plot argument, a pdf file with plots showing the quality of your sequencing data, as described above for barcode detection, will be generated.

&emsp;
&emsp;

# BarPepAnalysis

This Script can be applied both for barcoded and peptide insertion screens, indicated by running the script in the respective _mode_. When running the script in barcode analysis mode, it performs three (or optionally five) normalization steps, creating an output file for each one. When running the script in peptide analysis mode, it translates the found DNA sequences of the peptides into amino acid sequences and ranks them by their frequency.

&emsp;

## Prepare your Input Files
### _For Barcode Analysis_

In this mode the script requires the output from the _BarPepDetection_ script as well as one input file containing a table with your filenames, animals, tissues, and a weight_varialble (vg/dg). This file is easily prepared in Excel or a similar program and then saved as a CSV file. **The separator in the CSV file has to be a ','.** 

&emsp;

The input file MUST follow this structure:
&emsp;

| Filename                | SampleType| Animal     |Tissue        | weight_variable|
|-------------------------|-----------|------------|--------------|----------------|
| Sample1_Variants.csv    | cDNA      | M1         | Heart        | 0.00635        |
| Sample2_Variants.csv    | gDNA      | M1         | Heart        | 0.00635        |
| Sample3_Variants.cs     | cDNA      | M1         | Kidney       | 0.000293       |
| Sample4_Variants.csv    | gDNA      | M1         | Kidney       | 0.000293       |
| Sample5_Variants.csv    | cDNA      | M2         | Lung         | 0.00871        |
| Sample6_Variants.csv    | gDNA      | M2         | Lung         | 0.00871        |

&emsp;

**IMPORTANT:** The names of the columns have to be written exactly how it is shown here! Do not use a '_' in your animal, sample type or tissue entries!

The weight variable is used for the calculation of the B<sub>αβ</sub>, V<sub>αβ</sub>, and T<sub>αβ</sub> values (see output files). Traditionally these are vg/dg measurements of every tissue, but can be replaced by a value of choice.
  
Other input files that you will need, but won't have to be specially prepared:

  - A folder with **the _Variant.csv output files** from the _BarPepDetection_ script containing the counts of the expected variants.  

  - **The _Variant.csv output file of your Input Library** from the _BarPepDetection_ script. The values for the normalization to the input library will be directly calculated from this.  

**Note**: If any column contains a zero value, a pseudocount of default 1x10<sup>-6</sup> will be added to every column. A warning will be printed in the terminal. 

&emsp;

### _For Peptide Analysis_
In this mode the script requires the output from the _BarPepDetection_ script as well as one input file containing a table with your filenames, animals, and tissues. This file is easily prepared in Excel or a similar program and then saved as a CSV file. **The separator in the CSV file has to be a ','.** 
&emsp;

The input file MUST follow this structure:
&emsp;

| Filename        | Animal     |Tissue        |       
|-----------------|------------|--------------|
| Sample1_PV.csv  | M1         | Heart        |
| Sample2_PV.csv  | M1         | Kidney       |
| Sample3_PV.csv  | M1         | Lung         |
| Sample4_PV.csv  | M2         | Kidney       |
| Sample5_PV.csv  | M2         | Lung         |
| Sample6_PV.csv  | M3         | Heart        |

&emsp;

**IMPORTANT:** The names of the columns have to be written exactly how it is shown here! Do not use a '_' in your animal or tissue entries!

Other input files that you will need, but won't have to be specially prepared:
  - A folder with **the _PV.csv output files from the _BarPepDetection_** containing the counts of the found peptides.

&emsp;

## Arguments BarPepAnalysis
If you are running the script for the first time, open the Terminal and call the script with 
```
python BarPepAnalysis.py -h
````

This will give you a list of the arguments that you need to specifiy for running the script successfully:  

&emsp;

**<p style="font-size:15px;">Required arguments:</p>**
- **-a MODE:**  
 Specify in which mode you want to run the script. If you want to use it for analysis of barcodes, type "BC", if you want to use it for analysis of peptide insertions, type "PV". 

- **-i INPUTFILE:**  
 Give the path to your input CSV file. For barcode analysis, this should contain the filenames, sample types, animals, tissues, and vg/dg or RQ values. For peptide analysis, it only includes the filenames, animals, and tissues.  

 - **-d DIRECTORY:**  
 Give the path to the directory containing the _Variants.csv  or _PV.csv output files from the _Detection Script_. 

&emsp;

 **<p style="font-size:15px;">Optional arguments:</p>**

- **-o OUTPUTDIR:**  
Give the path to the directory where the output files should be saved. If not specified, a new folder in the directory with the input data will be automatically created for the output files.

- **-w SILENCE:**  
Set this flag if you want to avoid printouts in the terminal.  

- **-z VERSION:**  
Set this flag if you want to show the script's version number and exit.

- **-p PSEUDO:**
Pseudocount to be added to columns containing zeroes. Default 1e-6.

&emsp;

**<p style="font-size:15px;">Additional arguments for barcode analysis:</p>**

- **-l LIBRARYNORM:**  
 Give the path to the _Variants.csv output file of the input library from the _Detection Script_. This argument is required for barcode analysis.

- **-x EXTRANORM:**  
 Set this flag if you additionally want to compute the V<sub>αβ</sub> and T<sub>αβ</sub> values.

&emsp;

**<p style="font-size:15px;">Additional arguments for peptide analysis:</p>**

- **-t TOPNUMBER:**  
You can specify the number of peptide sequences shown in the output files. The default is set to 100, i.e., the output will be limited to the top 100 ranking peptide sequences.

&emsp;

## Output Files BarPepAnalysis
### _For Barcode Analysis_
Running the script in barcode analysis mode will generate at least 4 output CSV files.

- **Read Counts File:**  
This CSV file contains a table with the raw read counts of each variant in each tissue for each unique combination of animal and sample type. Additionally, it shows the vg/dg or RQ values as well as the input library normalization values. Each output file will follow this table structure.  
  
- **P<sub>αβ</sub>  File:**  
This CSV file contains the proportional read count values, or P<sub>αβ</sub> values. They are calculated by normalizing the read counts R of all variants α in tissue β to the sum of all variants α in β:

$$
\displaystyle
\ P_{αβ}= \frac{R_{αβ}}{\sum_{α} R_{αβ}}
$$  

- __P*<sub>αβ</sub> File:__  
This CSV file contains the proportional count values normalized to the input library, or P*<sub>αβ</sub> values. They are calculated by normalizing P<sub>αβ</sub> to the proportion of each variant α in the initial library L<sub>α</sub>, thus correcting for the uneven composition in library:

$$
\displaystyle
\ P_{αβ}^*= \frac{P_{αβ}}{L_{α}}
$$

- **B<sub>αβ</sub> File:**  
In this CSV file,  P*<sub>αβ</sub> is weighted by the weight_variable (e.g. vg/dg or RQ values), termed G<sub>β</sub>, to allow a comparison of one variant α over all analyzed tissues β:

$$
\displaystyle
\ B_{αβ}= \frac{P_{αβ}}{L_{α}}*G_β
$$

&emsp;

Optional:

- **V<sub>αβ</sub> File:**  
In this CSV file,  B<sub>αβ</sub> values are shown as proportions of the sum over all variants α of B<sub>αβ</sub>. These values can be useful to create bar plots which demonstrate the proportion of all variants α in one tissue β, exemplifying the efficiency of the individual vectors:  

$$
\displaystyle
\ V_{αβ}= \frac{B_{αβ}}{\sum_{α} B_{αβ}}
$$    

- **T<sub>αβ</sub> File:**  
In this CSV file,  B<sub>αβ</sub> values are shown as proportions of the sum over all tissues β of B<sub>αβ</sub>. These values can be useful to create bar plots which show the proportion of one variant α in all tissues β, allowing an analysis of the tissue specificity:

$$
\displaystyle
\ T_{αβ}= \frac{B_{αβ}}{\sum_{β} B_{αβ}}
$$  

&emsp;

### _For Peptide Analysis_
Running the script in peptide analysis mode will generate 2 output CSV files.

- **Ranking DNAseq File:**  
This CSV file contains a table with the DNA sequences of the ranked top 100 (or otherwise specified number) peptides for each animal and tissue combination, their raw counts, and proportions.  

- **Ranking PeptideSeq File:**  
This CSV file contains a table with the amino acid sequences of the ranked top 100 (or otherwise specified number) peptides for each animal and tissue combination, their raw counts, and proportions.


</div>
