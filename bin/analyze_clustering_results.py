#!/usr/bin/env python3

import argparse
import pysam
import numpy as np
        

def rev_comp(seq):
    rev_dict = {"A":"T", "T":"A", "C":"G", "G":"C"}
    return "".join([rev_dict[i] for i in seq[::-1]])

def analyze_cluster_results(cluster_file, cluster_out_fh, mod_bam, telo_stats_fh):
    
    cluster_sizes = {}
    cluster_dict = {}
    
    with open(cluster_file, "r") as cluster_fh:
        linecount = 0
        for line in cluster_fh:
            if linecount == 0:
                linecount += 1
                continue
            line = line.strip().split()
      
            cluster_sizes[int(line[3])] = len(line[-1].split(","))
            cluster_dict[int(line[3])] = [name.strip("@") for name in line[-1].split(",")]
    
    seq_count = sum([cluster_sizes[cluster] for cluster in cluster_sizes])

    cluster_perc = {}
    final_cluster_size = {}
    final_cluster_dict = {}
    for cluster in cluster_sizes:
        if cluster_sizes[cluster] > 0.002 * sum(cluster_sizes.values()):
            final_cluster_size[cluster] = cluster_sizes[cluster]
            final_cluster_dict[cluster] = cluster_dict[cluster]
            
    for cluster in final_cluster_size:     
        cluster_perc[cluster] = final_cluster_size[cluster] / sum(final_cluster_size.values()) * 100
        
    print("Number of Clusters: {}".format(len(final_cluster_size)))
    print("Percentage of Reads Clustered: {}".format(sum([final_cluster_size[i] for i in final_cluster_size])/seq_count * 100))
    print("Average Cluster Percentage: {}".format(np.mean(list(cluster_perc.values()))))
    
    telo_seqs = []
    with open(cluster_out_fh, "w") as cluster_fh:
        cluster_fh.write("read_id\tcluster\n")
        for cluster in final_cluster_dict:
            for read in final_cluster_dict[cluster]:
                telo_seqs.append(read)
                cluster_fh.write("{}\t{}\n".format(read, cluster))

    
    # create cluster specific fasta files - not needed yet 
    # # load telomeric stats 
    # telo_stats = {}
    # with open(telo_stats_fh, "r") as fh:
    #     linecount = 0
    #     for line in fh:
    #         if linecount == 0:
    #             linecount += 1
    #             continue
    #         line = line.strip().split()
    #         telo_stats[line[0]] = {"strand":line[1]}
    
    # telo_dict = {}
    # aln_file = pysam.AlignmentFile(mod_bam, "r", check_sq=False)
    # for aln in aln_file:
    #     if aln.query_name in telo_seqs and not aln.is_supplementary and not aln.is_secondary:
            
    #         if aln.is_reverse:
    #             aln.query_sequence = rev_comp(aln.query_sequence)
  
    #         if telo_stats[aln.query_name]["strand"] == "G":
    #             telo_dict[aln.query_name] = aln.query_sequence
    #         else:
    #             telo_dict[aln.query_name] = rev_comp(aln.query_sequence)

    # for cluster in final_cluster_dict:
    #     with open("cluster_{}.fa".format(cluster), "w") as cluster_fh:
    #         for read in final_cluster_dict[cluster]:
    #             cluster_fh.write(">{}\n{}\n".format(read, telo_dict[read]))
            
def main(args):

    analyze_cluster_results(args.cluster_file, args.cluster_out_fh, args.mod_bam, args.telo_stats)


def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cluster_file",required=True)
    parser.add_argument("--cluster_out_fh", required=True)
    parser.add_argument("--mod_bam", required=True)
    parser.add_argument("--telo_stats", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)