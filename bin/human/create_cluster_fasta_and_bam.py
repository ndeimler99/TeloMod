#!/usr/bin/env python3

import argparse
import pysam
import numpy as np
import math
import gzip
import cairo

def rev_comp(seq):
    rev_dict = {"A":"T", "T":"A", "C":"G", "G":"C"}
    return ''.join([rev_dict[i] for i in seq[::-1]])
                
def get_cluster_dict(cluster_fh):
    
    cluster_dict = {}

    with open(cluster_fh, "r") as clusters:
        linecount = 0
        for line in clusters:
            if linecount == 0:
                linecount += 1
                continue

            line = line.strip().split()
            if int(line[1]) not in cluster_dict:
                cluster_dict[int(line[1])] = []
            
            cluster_dict[int(line[1])].append(line[0])
    
    return cluster_dict


def get_telo_dict(mod_bam, telo_stats):

    aln_file = pysam.AlignmentFile(mod_bam, "rb", check_sq=False)
    aln_dict = {}
    telo_reads = {}

    for aln in aln_file:
        if aln.query_name in telo_stats and not aln.is_supplementary and not aln.is_secondary:
            
            if aln.is_reverse:
                aln.query_sequence = rev_comp(aln.query_sequence)
            #telo_orientated[aln.query_name] = aln.query_sequence
            if telo_stats[aln.query_name]["strand"] == "G":
                telo_reads[aln.query_name] = aln.query_sequence
            else:
                telo_reads[aln.query_name] = rev_comp(aln.query_sequence)
            
            aln_dict[aln.query_name] = aln
    
    return telo_reads, aln_dict

def get_stats_dict(telo_stats_fh):
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
                #telo_stats[line[0]] = {"strand":line[1], "read_length":len(telo_dict[line[0]]), "telo_length":int(line[4]), "telo_start":int(line[3]), "telo_end":int(line[3])+int(line[4])}

    return telo_stats



def main(args):

    # load cluster dict from out_file
    cluster_dict = get_cluster_dict(args.cluster_file)
    stats_dict = get_stats_dict(args.telo_stats) 
    telo_dict, aln_dict = get_telo_dict(args.mod_bam, stats_dict)


    aln_file = pysam.AlignmentFile(args.mod_bam, "rb", check_sq=False)

    for cluster in cluster_dict:
        out_aln = pysam.AlignmentFile("cluster_{}.bam".format(cluster), "wb", template=aln_file)
        with open("cluster_{}.fa".format(cluster), "w") as fasta_fh:
            for read in cluster_dict[cluster]:
                fasta_fh.write(">{}\n{}\n".format(read, telo_dict[read]))
                out_aln.write(aln_dict[read])


def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cluster_file",required=True)
    parser.add_argument("--mod_bam", required=True)
    parser.add_argument("--telo_stats", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)