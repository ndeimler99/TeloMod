#!/usr/bin/env python3

import argparse
import pysam
import pandas as pd


def rev_comp(seq):
    rev_dict = {"A":"T", "T":"A", "C":"G", "G":"C"}
    return ''.join([rev_dict[i] for i in seq[::-1]])


def main(args):

    # load cluster dict from out_file
    stats_df = pd.read_table(args.telo_stats, sep="\t")
    aln_file = pysam.AlignmentFile(args.telobam, "r", check_sq=False)

    aln_dict = {}
    telo_dict = {}

    for aln in aln_file:
        # if aln.query_name == "c48e0977-0e21-4d7b-a81c-c82c7aebc6cd":
        #     print(stats_df.loc[stats_df["read_id"]==aln.query_name, "strand"].item())
        #     if stats_df.loc[stats_df["read_id"]==aln.query_name, "strand"].item() == "C":
        #         print("C Identified")

        if stats_df.loc[stats_df["read_id"]==aln.query_name, "strand"].item() == "C":
            telo_dict[aln.query_name] = rev_comp(aln.query_sequence)
        else:
            telo_dict[aln.query_name] = aln.query_sequence

        aln_dict[aln.query_name] = aln

    for cluster in stats_df["Cluster"].dropna().unique():
        cluster_bam = pysam.AlignmentFile("cluster_{}.bam".format(int(cluster)), "wb", template=aln_file)
        with open("cluster_{}.fa".format(int(cluster)), "w") as fasta_fh:
            for read in stats_df.loc[stats_df["Cluster"]==cluster]["read_id"]:
                if len(telo_dict[read]) > 50000:
                    telo_dict[read] = telo_dict[read][-30000:]
                fasta_fh.write(">{}\n{}\n".format(read, telo_dict[read]))
                cluster_bam.write(aln_dict[read])
                
def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--telo_stats",required=True)
    parser.add_argument("--telobam", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)