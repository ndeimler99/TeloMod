#!/usr/bin/env python3

import argparse
import pysam
import numpy as np
import math
import gzip
import cairo
import pandas as pd

def rev_comp(seq):
    rev_dict = {"A":"T", "T":"A", "C":"G", "G":"C"}
    return ''.join([rev_dict[i] for i in seq[::-1]])

def plot_cluster_telo_seqs(telo_dict, telo_stats, mod_dict, file_name, c_strand_only, image_width=1500, image_height=1000):
    
    print("Script Started")
    max_telo_length = max([telo_stats[read]["vrr_telo_length"] for read in telo_stats])#list(telo_stats["vrr_telo_length"]))
    #print(max_telo_length)
    subtelo_stretch_max = 5000
    seq_length = subtelo_stretch_max + max_telo_length

    #sorted_telos = list(telo_stats.sort_values(by="vrr_telo_length", ascending=True)["read_id"])
    sorted_telos = sorted(telo_stats, key=lambda k: telo_stats[k]["vrr_telo_length"], reverse=False)

    strand_spacer = 50
    y_offset_top = 50
    y_offset_bottom = 50
    x_offset_left = 75
    x_offset_right = 10

    if not c_strand_only:
        seq_height = (image_height - strand_spacer - y_offset_top - y_offset_bottom)/len(sorted_telos)
    else:
        seq_height = (image_height - y_offset_top - y_offset_bottom)/len(sorted_telos)

    nucl_width = (image_width - x_offset_left - x_offset_right)/seq_length
    
    g_strands = len([telo_stats[read]["strand"] for read in telo_stats if telo_stats[read]["strand"] == "G"])
    c_strands = len([telo_stats[read]["strand"] for read in telo_stats if telo_stats[read]["strand"] == "C"])
    
    print(g_strands)
    print(seq_height)
    print(seq_height*g_strands)
  

    with cairo.PDFSurface(file_name, image_width, image_height) as surface:
    #if True:
#         surface = cairo.ImageSurface(cairo.FORMAT_RGB24,
#                                  image_width,
#                                  image_height)
        ctx = cairo.Context(surface)
        ctx.rectangle(0, 0, image_width, image_height)
        ctx.set_source_rgb(1,1,1)
        ctx.fill()
        
        ctx.set_font_size(9)
        ctx.select_font_face("Ariel",
                              cairo.FONT_SLANT_NORMAL,
                              cairo.FONT_WEIGHT_NORMAL)
        
        ### draw C telomeres 
        ctx.set_source_rgb(0,0,0)

        ctx.rectangle(x_offset_left, image_height-y_offset_bottom, image_width-x_offset_left - x_offset_right, 3)
        ctx.fill()

        # draw y axis
        if not c_strand_only:
            ctx.rectangle(x_offset_left - 3, y_offset_top+g_strands*seq_height+strand_spacer, 3, c_strands*seq_height + 3)
            
            ctx.save()
            a,b,width,height,c,d = ctx.text_extents('C Strand')
            ctx.translate(30, image_height - y_offset_bottom - (seq_height*c_strands)/2)
            ctx.rotate(-math.pi / 2)   
            ctx.move_to(-width/2, 0)
            ctx.text_path('C Strand')
            ctx.fill()
            ctx.restore()

            ctx.save()
            a,b,width,height,c,d = ctx.text_extents('G Strand')
            ctx.translate(30, image_height - y_offset_bottom - (seq_height*c_strands)-strand_spacer - (seq_height*g_strands/2))
            ctx.rotate(-math.pi / 2)   
            ctx.move_to(-width/2, 0)
            ctx.text_path('G Strand')
            ctx.fill()
            ctx.restore()
        else:
            ctx.rectangle(x_offset_left - 3, y_offset_top, 3, c_strands*seq_height + 3)

            ctx.save()
            a,b,width,height,c,d = ctx.text_extents('Number of Reads')
            ctx.translate(30, image_height - y_offset_bottom - (seq_height*c_strands)/2)
            ctx.rotate(-math.pi / 2)   
            ctx.move_to(-width/2, 0)
            ctx.text_path('Number of Reads')
            ctx.fill()
            ctx.restore()
        
        a,b,width,height,c,d = ctx.text_extents('Telomere Position')
        ctx.move_to(x_offset_left + (image_width - x_offset_left - x_offset_right)/2 - width/2, image_height - 10)
        ctx.text_path("Telomere Position")
        ctx.fill()

        # label x-axis positions
        ctx.set_source_rgb(0, 0, 0)
        for i in range(0, max_telo_length, 1000):
            if i % 3000 == 0:
                a,b,width,height,c,d = ctx.text_extents('{}'.format(i))
                #ctx.move_to(x_offset_left + i*nucl_width - width/2 + 1000*nucl_width - 2.5, image_height - y_offset_bottom + 25)
                ctx.move_to(x_offset_left + i*nucl_width - width/2 + subtelo_stretch_max*nucl_width, image_height - y_offset_bottom + 13)
                ctx.text_path('{}'.format(i))
                ctx.fill()
            ctx.rectangle(x_offset_left + i*nucl_width -1+ subtelo_stretch_max*nucl_width, image_height-y_offset_bottom, 2, 4)
            ctx.fill()
        
        for i in range(0, subtelo_stretch_max+1, 1000):
            if i == 0:
                continue
            if i % 1000 == 0:
                a,b,width,height,c,d = ctx.text_extents('-{}'.format(i))
                #ctx.move_to(x_offset_left + i*nucl_width - width/2 + 1000*nucl_width - 2.5, image_height - y_offset_bottom + 25)
                ctx.move_to(x_offset_left + subtelo_stretch_max*nucl_width - i*nucl_width - width/2, image_height - y_offset_bottom + 13)
                ctx.text_path('-{}'.format(i))
                ctx.fill()
            ctx.rectangle(x_offset_left + subtelo_stretch_max*nucl_width - i*nucl_width - 1, image_height-y_offset_bottom, 2, 4)
            ctx.fill()

        
        # label y-axis positions
        for i in range(0, c_strands, 20):
            if i % 100 == 0:
                a,b,width,height,c,d = ctx.text_extents('{}'.format(i))
                ctx.move_to(x_offset_left - width - 8, image_height - y_offset_bottom - i*seq_height + (height+1)/2)
                ctx.text_path('{}'.format(i))
            ctx.rectangle(x_offset_left-5, image_height - y_offset_bottom - i*seq_height, 4, 2)
            ctx.fill()
                
        # draw telomeres
        y_rect = image_height - y_offset_bottom - seq_height
        for telo in sorted_telos:
            if telo_stats[telo]["strand"] == "G":
                continue
    
            telo_seq = telo_dict[telo][telo_stats[telo]["vrr_start_pos"]:telo_stats[telo]["telo_end"]]
            
            subtelo = telo_dict[telo][0:telo_stats[telo]["vrr_start_pos"]][-5000:]
            x_rect = x_offset_left + subtelo_stretch_max*nucl_width
            i = 0
            while i < len(telo_seq):
                if telo_seq[i:i+6] == 'GGTTAG':
                    i += 6
                    ctx.set_source_rgba(0.118, 0.533, 0.898, 0.3)
                    ctx.rectangle(x_rect, y_rect, 6*nucl_width, seq_height)
                    ctx.fill()
                    x_rect += 6*nucl_width
                else:
                    ctx.set_source_rgba(1, 0.757, 0.027, 0.3)
                    ctx.rectangle(x_rect, y_rect, 1*nucl_width, seq_height)
                    ctx.fill()
                    x_rect += 1*nucl_width
                    i+=1
    
            # draw subtelomere
            i = 0
            x_rect = subtelo_stretch_max*nucl_width + x_offset_left
            while i < len(subtelo):
                if subtelo[-i:-i-6] == "GGTTAG":
                    i += 6
                    ctx.set_source_rgba(0.118, 0.533, 0.898, 0.3)
                    ctx.rectangle(x_rect-6*nucl_width, y_rect, 6*nucl_width, seq_height)
                    ctx.fill()
                    x_rect -= 6*nucl_width
                else:
                    i += 1
                    ctx.set_source_rgba(1, 0.757, 0.027, 0.3)
                    ctx.rectangle(x_rect-1*nucl_width, y_rect, 1*nucl_width, seq_height)
                    ctx.fill()
                    x_rect -= 1*nucl_width     
    
            # draw modifications
            
            if telo in mod_dict:
                telo_start = telo_stats[telo]["vrr_start_pos"]
                for pos in mod_dict[telo]["pos"]:
                    if pos < telo_start - subtelo_stretch_max:
                        continue
                    pos = pos - (telo_start - subtelo_stretch_max)
                    
                    ctx.set_source_rgb(0, 0.302, 0.251)
                    ctx.rectangle(x_offset_left+pos*nucl_width, y_rect, 1*nucl_width, seq_height)
                    ctx.fill()
                
            y_rect -= seq_height 


        if not c_strand_only:     
            # draw G telomeres
            ctx.set_source_rgb(0,0,0)
            # draw y axis
            ctx.rectangle(x_offset_left - 3, y_offset_top, 3, g_strands*seq_height + 3)
            ctx.fill()
            
            # label y-axis positions
            for i in range(0, g_strands, 20):
                if i % 100 == 0:
                    a,b,width,height,c,d = ctx.text_extents('{}'.format(i))
                    ctx.move_to(x_offset_left - width - 8, y_offset_top + g_strands*seq_height - i*seq_height + (height+1)/2)
                    ctx.text_path('{}'.format(i))
                ctx.rectangle(x_offset_left-5, y_offset_top + g_strands*seq_height - i*seq_height, 4, 2)
                ctx.fill()
            
            
            # draw telomeres
            y_rect = y_offset_top + g_strands*seq_height - seq_height
                
            for telo in sorted_telos:
                if telo_stats[telo]["strand"] == "C":
                    continue
                
                telo_seq = telo_dict[telo][telo_stats[telo]["vrr_start_pos"]:telo_stats[telo]["telo_end"]]
                subtelo = telo_dict[telo][0:telo_stats[telo]["vrr_start_pos"]][-5000:]
                x_rect = x_offset_left + subtelo_stretch_max*nucl_width
                i = 0
                while i < len(telo_seq):
                    if telo_seq[i:i+6] == 'GGTTAG':
                        i += 6
                        ctx.set_source_rgba(0.118, 0.533, 0.898, 0.3)
                        ctx.rectangle(x_rect, y_rect, 6*nucl_width, seq_height)
                        ctx.fill()
                        x_rect += 6*nucl_width
                    else:
                        ctx.set_source_rgba(1, 0.757, 0.027, 0.3)
                        ctx.rectangle(x_rect, y_rect, 1*nucl_width, seq_height)
                        ctx.fill()
                        x_rect += 1*nucl_width
                        i+=1
        
                # draw subtelomere
                i = 0
                x_rect = subtelo_stretch_max*nucl_width + x_offset_left
                while i < len(subtelo):
                    if subtelo[-i:-i-6] == "GGTTAG":
                        i += 6
                        ctx.set_source_rgba(0.118, 0.533, 0.898, 0.3)
                        ctx.rectangle(x_rect-6*nucl_width, y_rect, 6*nucl_width, seq_height)
                        ctx.fill()
                        x_rect -= 6*nucl_width
                    else:
                        i += 1
                        ctx.set_source_rgba(1, 0.757, 0.027, 0.3)
                        ctx.rectangle(x_rect-1*nucl_width, y_rect, 1*nucl_width, seq_height)
                        ctx.fill()
                        x_rect -= 1*nucl_width
                        
                if telo in mod_dict:
                    for pos in mod_dict[telo]["pos"]:
                        telo_start = telo_stats[telo]["vrr_start_pos"]#.loc[telo_stats["read_id"]==telo, "vrr_start_pos"].item()
                        if pos < telo_start - subtelo_stretch_max:
                            continue
                        pos = pos - (telo_start - subtelo_stretch_max)
                        
                        ctx.set_source_rgb(0, 0.302, 0.251)
                        ctx.rectangle(x_offset_left+pos*nucl_width, y_rect, 1*nucl_width, seq_height)
                        ctx.fill()
                
                
                y_rect -= seq_height            
            # draw modifications
        
    
#     # plot reads start with telomere and through telomere then backtrack and plot subtelomeric squences
#         # trim all reads to contain maximum of 10kb subtelomeric sequence
                
# def get_cluster_dict(cluster_fh):
    
#     cluster_dict = {}
#     telo_seqs = []

#     with open(cluster_fh, "r") as clusters:
#         linecount = 0
#         for line in clusters:
#             if linecount == 0:
#                 linecount += 1
#                 continue

#             line = line.strip().split()
#             telo_seqs.append(line[0])
#             if int(line[1]) not in cluster_dict:
#                 cluster_dict[int(line[1])] = []
            
#             cluster_dict[int(line[1])].append(line[0])
    
#     return cluster_dict, telo_seqs

def get_telo_dict(mod_bam, stats_dict):

    aln_file = pysam.AlignmentFile(mod_bam, "rb", check_sq=False)
    
    telo_reads = {}

    for aln in aln_file:
        if aln.query_name in stats_dict and not aln.is_supplementary and not aln.is_secondary:
            
            if aln.is_reverse:
                aln.query_sequence = rev_comp(aln.query_sequence)
            
            #telo_orientated[aln.query_name] = aln.query_sequence
            if stats_dict[aln.query_name]["strand"] == "G":
                telo_reads[aln.query_name] = aln.query_sequence
            else:
                telo_reads[aln.query_name] = rev_comp(aln.query_sequence)
    
    return telo_reads

def get_mod_dict(mod_table, modification, telo_stats):
      # load mod bam
    mod_dict = {}

    with gzip.open(mod_table, "rt") as fh:
        linecount = 0
        for line in fh:
            if linecount == 0:
                linecount += 1
                continue
            line = line.strip().split()
            if line[0] not in telo_stats:
                continue
            if line[0] not in mod_dict:
                mod_dict[line[0]] = {"full_read":0, "subtelo":0, "telo":0, "pos":[]}
            
            if line[19] == "true":
                continue
                
            if line[13] == modification:          
                if telo_stats[line[0]]["strand"] == "C":
                    mod_pos = telo_stats[line[0]]["read_len"] - int(line[1])
                else:
                    mod_pos = int(line[1])

                if mod_pos > telo_stats[line[0]]["telo_end"]:
                    continue
                elif mod_pos > telo_stats[line[0]]["vrr_start_pos"]:#telo_stats.loc[telo_stats["read_id"]==line[0], "vrr_start_pos"].item():
                    mod_dict[line[0]]["telo"] += 1
                else:
                    mod_dict[line[0]]["subtelo"] += 1

                mod_dict[line[0]]["full_read"] += 1
                mod_dict[line[0]]["pos"].append(mod_pos)
    return mod_dict


# def get_stats_dict(telo_stats_fh, seqs):
#     # load telomeric stats 
#     telo_stats = {}
#     with open(telo_stats_fh, "r") as fh:
#         linecount = 0
#         for line in fh:
#             if linecount == 0:
#                 linecount += 1
#                 continue
#             line = line.strip().split()

#             if line[0] in seqs:
#                 telo_stats[line[0]] = {"strand":line[1], "read_length":None, "telo_length":int(line[4]), "telo_start":int(line[3]), "telo_end":int(line[3])+int(line[4])}
#                 #telo_stats[line[0]] = {"strand":line[1], "read_length":len(telo_dict[line[0]]), "telo_length":int(line[4]), "telo_start":int(line[3]), "telo_end":int(line[3])+int(line[4])}

#     return telo_stats

def main(args):

    args.image_height = int(args.image_height)
    args.image_width = int(args.image_width)
    args.c_strand_only = args.c_strand_only == "true"
    # load cluster dict from out_file

    stats_table = pd.read_table(args.telo_stats, sep="\t")
    stats_table["telo_end"] = stats_table["vrr_start_pos"] + stats_table["vrr_telo_length"]

    telo_stats = {}
    for read in list(stats_table["read_id"]):
        telo_stats[read] = {"strand": stats_table.loc[stats_table["read_id"]==read, "strand"].item(),
                            "vrr_start_pos": stats_table.loc[stats_table["read_id"]==read, "vrr_start_pos"].item(),
                            "telo_end":  stats_table.loc[stats_table["read_id"]==read, "telo_end"].item(),
                            "read_len":stats_table.loc[stats_table["read_id"]==read, "read_len"].item(),
                            "vrr_telo_length":stats_table.loc[stats_table["read_id"]==read, "vrr_telo_length"].item(),
                            "Cluster":stats_table.loc[stats_table["read_id"]==read, "Cluster"].item()}

    telo_dict = get_telo_dict(args.mod_bam, telo_stats)
    print("Telo Dict Loaded")
    mod_dict = get_mod_dict(args.mod_table, args.modification, telo_stats)
    print("Mod Dict Loaded")

    for cluster in stats_table["Cluster"].dropna().unique():
        print(cluster)
        plot_cluster_telo_seqs(telo_dict, {k: v for k, v in telo_stats.items() if telo_stats[k]["Cluster"] == cluster}, mod_dict, "cluster_{}.pdf".format(int(cluster)), \
                            args.c_strand_only, image_width=args.image_width, image_height=args.image_height)


def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mod_bam", required=True)
    parser.add_argument("--telo_stats", required=True)
    parser.add_argument("--mod_table", required=True)
    parser.add_argument("--image_width", required=True)
    parser.add_argument("--image_height", required=True)
    parser.add_argument("--modification", required=True)
    parser.add_argument("--c_strand_only", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)