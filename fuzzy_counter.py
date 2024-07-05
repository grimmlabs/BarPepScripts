from Bio import SeqIO, Align, Seq
import time
import numpy as np
import multiprocessing as mp
import argparse
from argparse import RawTextHelpFormatter
import os
import pdb
import pandas as pd


def align_barcodes(read_name, read_sequence, barcodes, aligner, threshhold):
    qualities = {}
    qualities[read_name] = []
    for BC_name, BC_fwd in barcodes.items():
        alignments = aligner.align(read_sequence, BC_fwd)
        for alignment in sorted(alignments):
            if alignment.score >= threshhold:
                qualities[read_name].append([alignment.score, alignment, BC_name, read_sequence])
    return(qualities)

def filter_alignments_score(x):
    """
    Keeps only the best alignment.
    Reads with multiple best alignments are put into the multibarcoded category (weird flag).
    Reads with no alignments into the no barcodes category.
    """
    weird_reads = []
    unbarcoded_reads = []
    new_dict = {}
    for read_name, aln_list in x.items():
        weird_flag = False
        aln_list = np.array(aln_list, dtype=object)
        if len(aln_list) == 0: # if there are no alignments
            unbarcoded_reads.append(read_name)
            continue
        score_list = np.array([x[0] for x in aln_list])
        max_score = int(np.max(score_list))
        # check if the maximum score appears more than once
        max_score_count = np.sum(score_list == max_score)
        if max_score_count > 1:
            # make a list of all the barcodes with the maximum score
            max_score_aln_list = [x[2] for x in aln_list[score_list == max_score]]
            # when the aligned, top-score barcodes are the same then do nothing
            # the first alignment is then chosen
            if len(np.unique(max_score_aln_list)) != 1:
                # the read has at least 2 different barcodes aligned with the same quality
                # set the weird flag
                weird_flag = True
        if weird_flag:
            weird_reads.append(read_name)
            continue
        new_dict[read_name] = [aln_list[np.argmax(score_list)]]
    return(new_dict, unbarcoded_reads, weird_reads)

def prepare_arguments(x, barcodes, aligner, threshold):
    # prepare arguments for mp
    async_readnames_list = list(x.keys())
    async_sequences_list = [x[j] for j in async_readnames_list]
    async_barcode_list = [barcodes for j in range(len(x))]
    async_aligner_list = [aligner for j in range(len(x))]
    async_thresh_list = [threshold for i in range(len(x))]

    starmap_args = []
    for i in range(len(x)):
        starmap_args.append((
            async_readnames_list[i],
            async_sequences_list[i],
            async_barcode_list[i],
            async_aligner_list[i],
            async_thresh_list[i]))
    return(starmap_args)


if __name__ == "__main__":

    # argument parser
    # Set up argument parser
    ap = argparse.ArgumentParser()

    # Create argument groups
    ap._action_groups.pop()
    required = ap.add_argument_group("Required Arguments")
    optional = ap.add_argument_group("Optional Arguments")

    required.add_argument("-i", '--input', type=str, required=True,  help='path to the single fastq file or directory containing multiple fastq files')
    required.add_argument("-v", '--variants', type=str, required=True,  help='path to the variants file in tsv format')
    required.add_argument("-l", '--barcode_length', type=str, required=False,  help='all barcode must be the same length. If this is not provided it is inferred from the last barcode in the variants file.')
    required.add_argument("-e", '--threshold', type=int, required=True,  help='This is the threshold for determining a correct alignment. E.g. with a threshold the same number as the length of the barcode and match_score, gap_open_score, and mismatch_score respectively 1, -1, -1 (default), then only perfectly matched barcodes are accepted. For more fuzzy matching one could consider a threshold 80% of the length of the barcode. Must be an integer.')
    required.add_argument("-t", "--threads", type=int, required=False, default=8)
    required.add_argument("-o", "--out_dir", type=str, required=True, help='Output directory for the demultiplexed file. Two files per barcode (fwd, rev) will be created, as well as a file with unbarcoded reads and multi-barcoded reads.')
    required.add_argument("-m", "--match", type=str, required=False, default="1")
    required.add_argument("-x", "--mismatch", type=str, required=False, default="-1")
    required.add_argument("-g", "--open_gap", type=str, required=False, default="-1")
    required.add_argument("-p", "--extend_gap", type=str, required=False, default="-10")
    optional.add_arguments("-rc", "--add_reverse", action="store_true", help="Set flag if you want to generate plots showing the quality of your sequencing run."))
    args = vars(parser.parse_args())

    # parameters
    ## aligner parameters
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.open_gap_score = int(args["open_gap"])
    aligner.extend_gap_score = int(args["extend_gap"])
    aligner.mismatch_score = int(args["mismatch"])
    aligner.match_score = int(args["match"])
    ## barcodes
    barcode_threshold = int(args["threshold"])

    pdb.set_trace()

    barcodes_FWD = pd.read_csv(args.variants, sep="\t", names = ["sequence", "name"]).set_index("name").to_dict()["sequence"]
    if add_reverse:
        barcodes_RC = {name+"RC":str(Seq.Seq(bc).reverse_complement()) for name, bc in barcodes_FWD.items()}
        barcodes = {**barcodes_FWD, **barcodes_RC}
    else:
        barcodes = barcodes_FWD

    if args["barcode_length"] == None:
        barcode_length = len(entry.seq)
    else:
        barcode_length = int(args["barcode_length"])


    ## fasta input
    fasta = {}
    for entry in SeqIO.parse(args["fasta"], "fasta"):
        fasta[entry.description] = str(entry.seq)
    ## out directory
    out_directory = args["out_directory"]
    if out_directory[-1] != "/":
        out_directory += "/"
    if not os.path.exists(out_directory):
        os.mkdir(out_directory)
    ## multiprocessing
    n_workers = args["threads"]


    # Main computation
    ## parallel workers
    p = mp.Pool(n_workers)
    start_time = time.time()
    ## do the calculation
    results = p.starmap(
        align_barcodes,
        prepare_arguments(fasta, barcodes, aligner, barcode_threshold)
    )

    ## combine the results
    qualities = {}
    for r in results:
        qualities.update(r)

    ## filter the alignments of every read, to only retain the best one
    filt_qualities, unbarcoded_reads, weird_reads = filter_alignments_score(qualities)

    # output
    ## make output dictionary
    out_dictionary = {}
    for barcode_name in list(barcodes.keys()):
        out_dictionary[barcode_name] = []

    ## populate output dict
    for read_name, aln_list in filt_qualities.items():
        out_dictionary[aln_list[0][2]].append([read_name, aln_list[0][3]])

    # format
    df = pd.DataFrame({key:len(val) for key, val in out_dictionary.items()}, index=[0])

    df.to_csv(out_directory + "results.csv")

    time_needed = (time.time() - start_time) / (len(fasta))
    time_per_1e6 = time_needed * 1e6
    print("Total time: {}s".format(round(time.time() - start_time, 3)))
    print("Time per million reads: {}s".format(round(time_per_1e6, 3)))
