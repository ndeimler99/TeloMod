#!/usr/bin/env python3

import argparse
import pysam
import gzip
import pandas as pd
import math

def rev_comp(seq):
    rev_dict = {"A":"T", "T":"A", "C":"G", "G":"C"}
    return ''.join([rev_dict[i] for i in seq[::-1]])

def get_telo_dict(telo_stats):
    
    telo_dict = {}
    with open(telo_stats,"r") as fh:
        linecount = 0
        for line in fh:
            if linecount == 0:
                linecount += 1
                continue
            line = line.strip().split()
            telo_dict[line[0]] = [line[1], int(line[3]), int(line[4])]
    return telo_dict


def cigar(full, trimmed):
    s = full.index(trimmed)
    e = len(full) - (s + len(trimmed))
    return f"{s}H{len(trimmed)}M{e}H"


def reorient_bam(bam_file, telo_stats, out_file, ref_bam):
    
    #telo_dict = get_telo_dict(telo_stats)
    
    telo_dict = pd.read_table(telo_stats, sep="\t")
    aln_file = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
    out_bam = pysam.AlignmentFile(out_file, "wb", template=aln_file)

    ref_dict = {}
    ref_file = pysam.AlignmentFile(ref_bam, "rb", check_sq=False)
    for aln in ref_file:
        ref_dict[aln.query_name] = aln 

    for aln in aln_file:
        if aln.query_name in list(telo_dict["read_id"]):
            if telo_dict.loc[telo_dict["read_id"]==aln.query_name, "strand"].item() == "C":
                qual = aln.query_qualities
                aln.query_sequence = rev_comp(aln.query_sequence)
                aln.query_qualities = qual[::-1]

            #aln.cigarstring = cigar(ref_dict[aln.query_name].query_sequence, aln.query_sequence)
            # tags = aln.get_tags()
            # tags = [(k, v) for k, v in tags if k not in {"XS", "XB"}]
            # print(tags)
            # aln.set_tags(tags)

            # print(aln.get_tags())

            out_bam.write(aln)
    
def main(args):

    reorient_bam(args.telo_bam, args.stats_file, args.out_bam, args.ref_bam)




def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats_file", required=True)
    parser.add_argument("--telo_bam", required=True)
    parser.add_argument('--out_bam', required=True)
    parser.add_argument("--ref_bam", required=True)

    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)