
#!/usr/bin/env python

#__==================================== Barcode & Peptide Analysis Script 2.0 ========================================__
#__=========================================== (written by E. Locke) =================================================__


# This script analyses output files from the BarPepDetection Script.
# For Barcode Analysis, it performs three (optionally five) normalization steps, creating output files for each step.
# For Peptide Analysis, it ranks the top DNA and amino acid sequences of detected peptides.



#_____________________________________________________NECESSARY IMPORTS______________________________________________________
#

import argparse
import pandas as pd
import os
from Bio.Seq import Seq
import numpy as np



#___________________________________________________DEFINITION OF ARGUMENTS____________________________________________________
#

# Set up argument parser
ap = argparse.ArgumentParser()

# Create argument groups
ap._action_groups.pop()
required = ap.add_argument_group("Required Arguments")
optional = ap.add_argument_group("Optional Arguments")
barcode = ap.add_argument_group("Additional Arguments for Barcode Analysis")
peptide = ap.add_argument_group("Additional Arguments for Peptide Analysis")

# Add arguments to argument groups
required.add_argument("-a", "--mode", required = True, choices=["BC", "PV"], help = "Which mode do you want to run the script in? Type 'BC' for barcode analysis or 'PV' for peptide analysis.")
required.add_argument("-i'", "--inputfile", required = True, help = "Path to your input CSV file containing the file names, animal, and tissue (+ for BC analysis sample type and weight_variable).")
required.add_argument("-d", "--directory", required = True, help = "Path to the directory containing the output CSV files from the Detection Script.")
optional.add_argument("-o", "--outputdir", required=False, help="Path to the directory where the output files should be saved. If not specified, a new folder in the directory with the input data will be automatically created for the output files.")
optional.add_argument("-w", "--silence", action="store_false", help="Set flag to avoid printouts in the terminal.")
optional.add_argument("-z" "--version", action="version", version="\n"*5+"====== Barcode & Peptide Analysis Script 2.0 (written by E. Locke) ======\n\n", help="Set flag if you want to print script's version number and exit.")
optional.add_argument("-p", "--pseudo", type=float, default=1e-6, help="Pseudocount to be added to columns containing zeroes.")
barcode.add_argument("-l", "--libraryNorm", required = False, help = "Path to the output CSV file of the input library from the Detection Script.")
barcode.add_argument("-x", "--extraNorm", action="store_true", help="Set flag to compute the V_αβ and T_αβ values.")
peptide.add_argument("-t", "--topNumber", required=False, type = int, default = 100, help="Number of peptide sequences you want to limit the output to. Default is set to 100.")
args = ap.parse_args()



#_____________________________________INITIALIZATION OF SILENCE VARIABLE TO AVOID PRINTOUTS____________________________________
#

silence = args.silence

if silence:
    print("\n"*5+"====== Barcode & Peptide Analysis Script 2.0 (written by E. Locke) ======\n\n")



#____________________________________________________NECESSARY FUNCTIONS_______________________________________________________
#

# Custom format function converts integer values to their string integer form
def custom_format(x):
    if pd.isna(x):
        return x
    if x == int(x):
        return str(int(x))
    else:
        return x


# Combining the DataFrames and saving them as one DataFrame in a CSV file
def save_combined_dataframe(df_dict, output_path, header, mode):
    combined_df_list = []

    # Mode Barcode Analysis
    if mode == "BC":
        for key, df in df_dict.items():
            animal, sample_type = key.split('_')
            df_copy = df.copy()  # Make a copy of the dataframe to avoid altering the original
            df_copy['Animal'] = animal
            df_copy['SampleType'] = sample_type
            combined_df_list.append(df_copy)

    # Mode Peptide Analysis
    if mode == "PV":
        for key, df in df_dict.items():
            animal, tissue = key.split('_')
            rank = df.index+1
            df_copy = df.copy()  # Make a copy of the dataframe to avoid altering the original
            df_copy['Animal'] = animal
            df_copy['Tissue'] = tissue
            df_copy['Rank'] = rank
            combined_df_list.append(df_copy)

    # Merge DataFrames
    combined_df = pd.concat(combined_df_list, ignore_index=True)

    # Save combined_df in a CSV file
    with open(output_path, 'w') as f:
        combined_df.to_csv(f, index=False)
    if silence:
        print(f"\n\n{header} DataFrames have been saved to {output_path}")



#__________________________________________INITIALIZATION OF VARIABLES WITH ARGUMENTS__________________________________________
#

my_dir = args.directory
inputlibrary = args.libraryNorm
input = args.inputfile

# Define the output directory
if args.outputdir:
    output_directory = args.outputdir
# If user did not specify an output directory, create a new folder in directory with input data
else:
    if args.mode == "BC":
        output_directory = os.path.join(my_dir, 'BarcodeAnalysis_Results')
    if args.mode == "PV":
        output_directory = os.path.join(my_dir, 'PeptideAnalysis_Results')
    os.makedirs(output_directory, exist_ok=True)



#___________________________________________________BARCODE ANALYSIS SCRIPT____________________________________________________
#

if args.mode == "BC":

    if silence:
        print("Barcode Analysis Script is running.\n")


    # Read in the CSV file containing the file names with sample type, animal, tissue, and DNA/RNA normalization values.
    inputfile = pd.read_csv(input, sep=',')

    # Count the number of zeros before replacement
    zero_count_before = (inputfile['weight_variable'] == 0).sum()
    # apply pseudocount
    inputfile['weight_variable'] = inputfile['weight_variable'] + args.pseudo
    # Print a warning message if any zero values were replaced
    if zero_count_before > 0:
        if silence:
            print(f"Warning: {zero_count_before} zero values found in weight_variable. {args.pseudo} was added to all values in column weight_variable.\n")

    # Print the inputfile DataFrame
    if silence:
        print("Given inputfile:\n")
        print(inputfile)


    # Read in the CSV file containing the variants and input library normalization values.
    variants_norm = pd.read_csv(inputlibrary, sep=',')

    # Select only the 'Variant' and 'Count' columns
    variants_norm = variants_norm[['Variant', 'Count']]
    # Count the number of zeros before replacement
    zero_count_before = (variants_norm['Count'] == 0).sum()
    # apply pseudocount
    variants_norm['Count'] = variants_norm['Count'] + args.pseudo
    # Print a warning message if any zero values were replaced
    if zero_count_before > 0:
        if silence:
            print(f"Warning: {zero_count_before} zero values found in counts. {args.pseudo} was added to all counts.")
    # Calculate the total count
    total_count = variants_norm['Count'].sum()
    # Calculate the proportions of the counts
    variants_norm['InputNorm'] = variants_norm['Count'] / total_count
    # Drop the 'Count' column
    variants_norm = variants_norm.drop(columns=['Count'])
    # Rename the columns to 'Variant' and 'InputNorm' explicitly
    variants_norm = variants_norm.rename(columns={'Variant': 'Variant', 'InputNorm': 'InputNorm'})
    # Print the resulting DataFrame
    if silence:
        print("\n\nComputed input library normalization values:\n")
        print(variants_norm)


    # Get unique values from the "Animal" and "SampleType" columns of the inputfile DataFrame
    unique_animals = inputfile["Animal"].unique()
    unique_sample_types = inputfile["SampleType"].unique()

    # Create an empty DataFrame for each combination of Animal and Sample Type and add the columns from variants_norm
    animal_sampletype_dfs = {}
    for animal in unique_animals:
        for sample_type in unique_sample_types:
            key = f"{animal}_{sample_type}"
            # Create a copy of variants_norm to avoid modifying the original
            animal_sampletype_dfs[key] = variants_norm.copy()

    # Add tissue-specific columns to the corresponding DataFrames
    for idx, row in inputfile.iterrows():
        animal = row["Animal"]
        sample_type = row["SampleType"]
        tissue = row["Tissue"]
        filename = row["Filename"]

        if os.path.exists(os.path.join(my_dir, filename)):
            # Read the CSV file for the current tissue
            tissue_df = pd.read_csv(my_dir+filename)

        # Print warning if file could not be found
        else:
            if silence:
                print(f"\n\nWarning: File {os.path.join(my_dir, filename)} does not exist. Skipping.\n")

        # Merge the tissue_df with the corresponding DataFrame on the "Variant" column
        key = f"{animal}_{sample_type}"
        animal_sampletype_dfs[key] = animal_sampletype_dfs[key].merge(tissue_df[["Variant", "Count"]].rename(columns={"Count": tissue}),
                                            on="Variant", how="left")

        # Shift the InputNorm column to the end of the dataframe
        new_columns = [col for col in animal_sampletype_dfs[key].columns if col != "InputNorm"] + ["InputNorm"]
        animal_sampletype_dfs[key] = animal_sampletype_dfs[key][new_columns]


    # Add weight_variable as the last row for each DataFrame
    for key in animal_sampletype_dfs.keys():

        # Initialize a dictionary for the new row
        new_row = {col: None for col in animal_sampletype_dfs[key].columns}
        new_row["Variant"] = "weight_variable"

        # Extract the corresponding rows from inputfile
        for idx, row in inputfile.iterrows():
            if f"{row["Animal"]}_{row["SampleType"]}" == key:
                new_row[row["Tissue"]] = row["weight_variable"]

        # Convert the new row to a DataFrame
        new_row_df = pd.DataFrame([new_row])

        # Append the new row DataFrame to the existing DataFrame using pd.concat
        animal_sampletype_dfs[key] = pd.concat([animal_sampletype_dfs[key], new_row_df], ignore_index=True)


    # Apply custom formatting, so that only numbers with non-zero decimal parts are shown with decimal places
    for key, df in animal_sampletype_dfs.items():
        for tissue in inputfile['Tissue'].unique():
            if tissue in df.columns:
                df[tissue] = df[tissue].apply(lambda x: custom_format(x) if pd.notnull(x) else x)



    #________________________________________________OUTPUTFILE 1: READ COUNTS______________________________________________________
    #

    # Save the read counts dataframes
    save_combined_dataframe(animal_sampletype_dfs, os.path.join(output_directory, 'BarcodeAnalysis_ReadCounts.csv'), "Read Counts", "BC")



    #_________________________________________OUTPUTFILE 2: PROPORTIONS (P_αβ VALUES)_________________________________________________
    #

    # Dictionary to hold new DataFrames with proportions
    proportion_dfs = {}

    # Calculate proportions
    for key, df in animal_sampletype_dfs.items():

        # Create a copy of the DataFrame excluding the last row ('weight_variable')
        temp_df = df.iloc[:-1].copy()

        # Exclude the 'Variant' and 'InputNorm' columns
        tissue_columns = temp_df.columns.difference(['Variant', 'InputNorm'])

        # Convert tissue/variant counts to floats
        temp_df[tissue_columns] = temp_df[tissue_columns].astype(float)

        # Calculate sums for each tissue column
        sums = temp_df[tissue_columns].sum()

        # Divide each value by the column sum to get the proportions
        for column in tissue_columns:
            # Check if the sum is zero to avoid division by zero
            if sums[column] != 0:
                temp_df[column] = temp_df[column] / sums[column]
            else:
                # If the sum is zero, set the proportion to zero
                temp_df[column] = np.where(temp_df[column] == 0, 0, np.inf)

        # Append the last row back to the DataFrame without performing any calculations on it
        proportion_df = pd.concat([temp_df, df.iloc[-1:].replace(np.inf, 0)])  # Use concat to append the last row

        # Add to the dictionary
        proportion_dfs[key] = proportion_df


    # Save the P_αβ dataframes
    save_combined_dataframe(proportion_dfs, os.path.join(output_directory, 'BarcodeAnalysis_P_ab.csv'), "P_αβ", "BC")



    #________________________________OUTPUTFILE 3: NORMALIZATION INPUT LIBRARY (P*_αβ VALUES)_______________________________________
    #

    # Dictionary to hold new DataFrames with divided proportions
    divided_proportion_dfs = {}

    # Divide proportions by InputNorm
    for key, df in proportion_dfs.items():

        # Create a copy of the DataFrame excluding the last row ('weight_variable')
        temp_df = df.iloc[:-1].copy()

        # Exclude the 'Variant' and 'InputNorm' columns
        tissue_columns = temp_df.columns.difference(['Variant', 'InputNorm'])

        # Divide each value by the corresponding InputNorm value
        temp_df[tissue_columns] = temp_df[tissue_columns].div(df['InputNorm'], axis=0)

        # Replace any resulting 'inf' values with 0
        temp_df.replace(np.inf, 0, inplace=True)

        # Append the last row back to the DataFrame without performing any calculations on it
        divided_proportion_df = pd.concat([temp_df, df.iloc[-1:]])

        # Add to the dictionary
        divided_proportion_dfs[key] = divided_proportion_df


    # Save the P*_αβ dataframes
    save_combined_dataframe(divided_proportion_dfs, os.path.join(output_directory, 'BarcodeAnalysis_Ps_ab.csv'), "P*_αβ", "BC")



    #________________________________OUTPUTFILE 4: NORMALIZATION TISSUES (B_αβ VALUES)____________________________________________
    #

    # Dictionary to hold new DataFrames with multiplied proportions
    multiplied_proportion_dfs = {}

    # Multiply proportions by weight_variable
    for key, df in divided_proportion_dfs.items():

        # Create a copy of the DataFrame excluding the last row ('weight_variable') and 'InputNorm' column
        temp_df = df.iloc[:-1].copy()

        # Exclude the 'Variant' and 'InputNorm' columns
        tissue_columns = temp_df.columns.difference(['Variant', 'InputNorm'])

        # Convert tissue/variant counts to floats
        temp_df[tissue_columns] = temp_df[tissue_columns].astype(float)

        # Get the weight_variable value from the last row based on the 'Variant' column
        vg_dg_value = df[tissue_columns].loc[df['Variant'] == 'weight_variable'].astype(float) if 'weight_variable' in df['Variant'].values else 1

        # Multiply each value by the weight_variable value
        temp_df[tissue_columns] = temp_df[tissue_columns].mul(vg_dg_value[tissue_columns].values)

        # Append the last row back to the DataFrame without performing any calculations on it
        multiplied_proportion_df = pd.concat([temp_df, df.iloc[[-1]]])

        # Add to the dictionary
        multiplied_proportion_dfs[key] = multiplied_proportion_df


    # Save the B_αβ dataframes
    save_combined_dataframe(multiplied_proportion_dfs, os.path.join(output_directory, 'BarcodeAnalysis_B_ab.csv'), "B_αβ", "BC")



    #_________________________________________________OPTIONAL NORMALIZATION STEPS___________________________________________________
    #

    if args.extraNorm:



        #_________________________________OUTPUTFILE 5: NORMALIZATION OVER VARIANTS (V_αβ VALUES)______________________________________
        #

        # Dictionary to hold new DataFrames with V_αβ values
        v_proportion_dfs = {}

        # Calculate V_αβ values
        for key, df in multiplied_proportion_dfs.items():

            # Create a copy of the DataFrame excluding the last row ('weight_variable')
            temp_df = df.iloc[:-1].copy()

            # Exclude the 'Variant' and 'InputNorm' columns
            tissue_columns = temp_df.columns.difference(['Variant', 'InputNorm'])

            # Convert tissue/variant counts to floats
            temp_df[tissue_columns] = temp_df[tissue_columns].astype(float)

           # Calculate the sum of B_αβ values for each column (excluding 'weight_variable' and 'InputNorm')
            column_sums = temp_df[tissue_columns].sum(axis=0)

            # Divide each B_αβ value by the corresponding column sum to get the V_αβ values
            for column in tissue_columns:
                temp_df[column] = temp_df[column] / column_sums[column]

            # Append the last row back to the DataFrame without performing any calculations on it
            v_proportion_df = pd.concat([temp_df, df.iloc[[-1]]])

            # Add to the dictionary
            v_proportion_dfs[key] = v_proportion_df

        # Save the V_αβ dataframes
        save_combined_dataframe(v_proportion_dfs, os.path.join(output_directory, 'BarcodeAnalysis_V_ab.csv'), "V_αβ", "BC")



        #_________________________________OUTPUTFILE 6: NORMALIZATION OVER TISSUES (T_αβ VALUES)____________________________________
        #

        # Dictionary to hold new DataFrames with T_αβ values
        t_proportion_dfs = {}

        # Calculate V_αβ values
        for key, df in multiplied_proportion_dfs.items():
            # Create a copy of the DataFrame excluding the last row ('weight_variable')
            temp_df = df.iloc[:-1].copy()

            # Exclude the 'Variant' and 'InputNorm' columns
            tissue_columns = temp_df.columns.difference(['Variant', 'InputNorm'])

            # Convert tissue/variant counts to floats
            temp_df[tissue_columns] = temp_df[tissue_columns].astype(float)

            # Calculate the sum of B_αβ values for each row (excluding 'Variant' and 'InputNorm')
            row_sums = temp_df[tissue_columns].sum(axis=1)

            # Divide each B_αβ value by the corresponding row sum to get the T_αβ values
            temp_df[tissue_columns] = temp_df[tissue_columns].div(row_sums, axis=0)

            # Append the last row back to the DataFrame without performing any calculations on it
            t_proportion_df = pd.concat([temp_df, df.iloc[[-1]]])

            # Add to the dictionary
            t_proportion_dfs[key] = t_proportion_df

        # Save the T_αβ dataframes
        save_combined_dataframe(t_proportion_dfs, os.path.join(output_directory, 'BarcodeAnalysis_T_ab.csv'), "T_αβ", "BC")



#_______________________________________________PEPTIDE ANALYSIS SCRIPT_______________________________________________________
#

if args.mode == "PV":

    if silence:
        print("Peptide Analysis Script is running.\n")



    #____________________________READING INPUT FILE & CREATING DATAFRAME WITH RAW READ COUNTS___________________________________
    #

    # Read the input file
    inputfile = pd.read_csv(input, sep=',')

    # Create a dictionary to store dataframes
    PV_dfs = {}

    # Loop through each row in the input file to read corresponding CSV files
    for idx, row in inputfile.iterrows():
        filename = row['Filename']
        animal = row['Animal']
        tissue = row['Tissue']

        # Construct the file path
        read_counts_file_path = os.path.join(my_dir, filename)

        # Read the CSV file into a DataFrame
        if os.path.exists(read_counts_file_path):
            df = pd.read_csv(read_counts_file_path).iloc[:, 1:]

            # Create a key for the dictionary
            key = f"{animal}_{tissue}"

            # Save the DataFrame in the dictionary
            PV_dfs[key] = df

        # Print warning if file could not be found
        else:
            if silence:
                print(f"\n\nWarning: File {read_counts_file_path} does not exist. Skipping.\n")


    # Print the keys of the dictionary to verify
    if silence:
        print("DataFrames loaded for the following animal_tissue combinations:")
    for key in PV_dfs.keys():
        print(key)



    #________________________________CALCULATING AND SAVING PROPORTIONS FOR DNA SEQUENCES___________________________________
    #

    # Create a dictionary to store the proportion DataFrames
    proportion_dfs = {}

    # Process each DataFrame to create proportion DataFrames
    for key, df in PV_dfs.items():

        # Calculate the sum of the counts
        count_sum = df['Count'].sum()

        # Create a new DataFrame for proportions
        proportion_df = df.copy()
        proportion_df['Proportion'] = proportion_df['Count'] / count_sum

        # Sort the DataFrame by proportions in descending order
        proportion_df = proportion_df.sort_values(by='Proportion', ascending=False)

        # Save the proportion DataFrame in the dictionary
        proportion_dfs[key] = proportion_df.head(args.topNumber)  # By default, keep only the top 100 sequences


    # Save the T_αβ dataframes
    save_combined_dataframe(proportion_dfs, os.path.join(output_directory, 'PeptideAnalysis_Ranking_DNAseq.csv'), "Ranking DNAseq", "PV")



    #__________________________________CALCULATING AND SAVING PROPORTIONS FOR AMINO ACID SEQUENCES_______________________________
    #

    # Create a dictionary to store the translated proportions DataFrames
    translated_dfs = {}

    # Translate DNA sequences to peptide sequences, recalculate counts and proportions
    for key, df in PV_dfs.items():
        df['PeptideSeq'] = df['Peptide'].apply(lambda x: str(Seq(x).translate()))
        df_grouped = df.groupby('PeptideSeq').agg({'Count': 'sum'}).reset_index()
        total_count = df_grouped['Count'].sum()
        df_grouped['Proportion'] = df_grouped['Count'] / total_count
        df_sorted = df_grouped.sort_values(by='Proportion', ascending=False).reset_index(drop=True)
        translated_dfs[key] = df_sorted.head(args.topNumber)  # Keep only the top 100 sequences


    # Save the T_αβ dataframes
    save_combined_dataframe(translated_dfs, os.path.join(output_directory, 'PeptideAnalysis_Ranking_PeptideSeq.csv'), "Ranking PeptideSeq", "PV")



# The script is completed!
if silence:
    print("\n\n======Script completed!======\n\n")
