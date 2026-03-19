#!/usr/bin/env python3

import argparse
import pysam
import gzip
import cairo
import math

def rev_comp(seq):
    rev_dict = {"A":"T", "T":"A", "C":"G", "G":"C"}
    return ''.join([rev_dict[i] for i in seq[::-1]])

def extract_telomeric_modbam(telo_stats_fh, mod_bam, out_bam):
    
    # load telomeric stats 
    telo_stats = {}
    with open(telo_stats_fh, "r") as fh:
        linecount = 0
        for line in fh:
            if linecount == 0:
                linecount += 1
                continue
            line = line.strip().split()
            telo_stats[line[0]] = {"strand":line[1], "read_length":None, "telo_length":int(line[4]), "telo_start":int(line[3]), "telo_end":int(line[3])+int(line[4])}
    
    modbam = pysam.AlignmentFile(mod_bam, "rb", check_sq=False)
    out_bam = pysam.AlignmentFile(out_bam, "wb", template=modbam)

    for aln in modbam:
        if aln.query_name in telo_stats:
            out_bam.write(aln)

        
def main(args):

    extract_telomeric_modbam(args.telo_stats, args.mod_bam, args.out_file)




def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--telo_stats", required=True)
    parser.add_argument("--mod_bam", required=True)
    parser.add_argument('--out_file', required=True)

    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)