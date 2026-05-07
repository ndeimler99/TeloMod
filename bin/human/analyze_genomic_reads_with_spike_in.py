#!/usr/bin/env python3

import gzip
import argparse
import pysam
from Bio import SeqIO

def get_ref_dict(reference_fasta):
    
    seqio_dict = SeqIO.to_dict(SeqIO.parse(reference_fasta, "fasta"))
    ref_dict = {}
    for chrom in seqio_dict:
        ref_dict[chrom] = str(seqio_dict[chrom].seq)
    return ref_dict

def get_telo_read_ids(telo_stats_fh):
    
    telo_seqs = []
    with open(telo_stats_fh, "r") as fh:
        linecount = 0
        for line in fh:
            if linecount == 0:
                linecount += 1
                continue
            line = line.strip().split()
            telo_seqs.append(line[0])
    return telo_seqs

def process_genomic_alignment(alignment_file, reference_dict, minimum_length, removed_seqs=[]):

    removed_seqs = set(removed_seqs)
    aln_file = pysam.AlignmentFile(alignment_file, "r", check_sq=False)
    seqs = {}

    for aln in aln_file:
        if aln.is_unmapped:
            continue     
        if aln.reference_name == "MT" or aln.reference_name.startswith("KI") or aln.reference_name.startswith("GL"):
            continue       
        if aln.query_name in removed_seqs:
            continue
        if not aln.is_secondary and not aln.is_supplementary:
            if len(aln.query_sequence) > minimum_length and aln.mapping_quality > 20:
                if aln.query_name not in seqs:
                    if reference_dict is not None:
                        seqs[aln.query_name] = {"query_seq":aln.query_sequence,
                                            "reference_seq": reference_dict[aln.reference_name][aln.reference_start:aln.reference_end]}
                    else:
                        seqs[aln.query_name] = {"query_seq":aln.query_sequence}
    
    return seqs



def process_genomic_stats(genomic_seqs, mod_table, ref_dict, mod_nucl, modification, out_file):
    
    mod_dict = {}
    mod_possible = {}
    mod_eliminated = {}

    for read in genomic_seqs:
        mod_possible[read] = genomic_seqs[read]["reference_seq"].count(mod_nucl)

    with gzip.open(mod_table, "rt") as mod_fh:
        linecount = 0
        for line in mod_fh:
            if linecount == 0:
                linecount += 1
                line = line.strip().split()
                #print(line)
                continue

            line = line.strip().split()
            if line[0] not in genomic_seqs or line[19] == "true":
                continue
            
            if line[0] not in mod_dict:
                mod_dict[line[0]] = 0
                mod_eliminated[line[0]] = 0
            
            rev_dict = {"A":"T", "C":"G", "G":"C", "T":"A"}
            
            if line[13] == modification and int(line[2]) != -1:
                #if in region of read that is within mapping
                mod_dict[line[0]] += 1
    
    with open(out_file, "w") as out_fh:
        
        out_fh.write("read_id\tstrand\tread_length\tread_type\ttotal_bases\ttotal_mods\tmod_percentage\tmods_eliminated\n")
        for read in mod_dict:

            if mod_possible[read] == 0:
                print("Passed")
                continue

            out_fh.write("{}\tNA\t{}\tgenomic\t{}\t{}\t{}\t{}\n".format(read, len(genomic_seqs[read]["reference_seq"]),
                                                            mod_possible[read], mod_dict[read],
                                                            mod_dict[read]/mod_possible[read] * 100, mod_eliminated[read]))

def main(args):

    args.minimum_read_length = int(args.minimum_read_length)

    ref_dict = get_ref_dict(args.reference_fa)
    print(len(ref_dict))
    telo_seqs = get_telo_read_ids(args.telo_stats)
    print(len(telo_seqs))
    spike_in_dict = process_genomic_alignment(args.spike_in_aln, None, args.minimum_read_length)
    print(len(spike_in_dict))
    print(ref_dict.keys())
    genomic_seqs_for_analysis = process_genomic_alignment(args.reference_aln, ref_dict, \
                                                        args.minimum_read_length, \
                                                        telo_seqs + list(spike_in_dict.keys()))


    process_genomic_stats(genomic_seqs_for_analysis, args.mod_calls, ref_dict, args.modified_nucleotide, \
                         args.modification_code, args.out_file)


def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference_aln",required=True)
    parser.add_argument("--spike_in_aln",required=True)
    parser.add_argument("--telo_stats", required=True)
    parser.add_argument("--reference_fa", required=True)
    parser.add_argument("--minimum_read_length", required=True)
    parser.add_argument("--mod_calls", required=True)
    parser.add_argument('--modified_nucleotide', required=True)
    parser.add_argument('--modification_code', required=True)
    parser.add_argument("--out_file", required=True)
    
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)