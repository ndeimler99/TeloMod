#!/usr/bin/env python3

import argparse
import pysam
import gzip
import cairo
import math

def rev_comp(seq):
    rev_dict = {"A":"T", "T":"A", "C":"G", "G":"C"}
    return ''.join([rev_dict[i] for i in seq[::-1]])

def plot_cluster_independent_telomeric_reads(mod_bam, telo_stats_fh, mod_table, mod_nucl, modification, file_name, summary_file, max_subtelo_stretch, image_width=1500, image_height=1000):
    
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
    
    # load telomeric reads
    mod_bam_fh = pysam.AlignmentFile(mod_bam, "rb", check_sq=False)
    telo_reads = {}
    telo_orientated = {}
    telo_possible = {}
    for aln in mod_bam_fh:
        if aln.query_name in telo_stats and not aln.is_supplementary and not aln.is_secondary:
            
            if aln.is_reverse:
                aln.query_sequence = rev_comp(aln.query_sequence)
                
            telo_reads[aln.query_name] = aln.query_sequence
            telo_stats[aln.query_name]["read_length"] = len(aln.query_sequence)
            
            #telo_orientated[aln.query_name] = aln.query_sequence
            if telo_stats[aln.query_name]["strand"] == "G":
                telo_orientated[aln.query_name] = aln.query_sequence
                telo_possible[aln.query_name] = {"full_read":telo_reads[aln.query_name].count(mod_nucl), 
                                      "subtelo":telo_reads[aln.query_name][0:telo_stats[aln.query_name]["telo_start"]].count(mod_nucl), 
                                      "telo":telo_reads[aln.query_name][telo_stats[aln.query_name]["telo_start"]:telo_stats[aln.query_name]["telo_end"]].count(mod_nucl)}
            else:
                telo_orientated[aln.query_name] = rev_comp(aln.query_sequence)
                telo_possible[aln.query_name] = {"full_read":telo_reads[aln.query_name].count(mod_nucl), 
                                      "subtelo":telo_reads[aln.query_name][telo_stats[aln.query_name]["read_length"]-telo_stats[aln.query_name]["telo_start"]:].count(mod_nucl), 
                                      "telo":telo_reads[aln.query_name][telo_stats[aln.query_name]["read_length"]-telo_stats[aln.query_name]["telo_end"]:telo_stats[aln.query_name]["read_length"]-telo_stats[aln.query_name]["telo_start"]].count(mod_nucl)}
            
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
                
            valid_mod = True
            if line[13] == modification:          
                # check strand and if a C strand telomere flip position to be G stranded
                # check position of mod, within telomere, subtelomere, or past end of telomere?
                if valid_mod:
                    if telo_stats[line[0]]["strand"] == "C":
                        mod_pos = telo_stats[line[0]]["read_length"] - int(line[1])
                    else:
                        mod_pos = int(line[1])

                    if mod_pos > telo_stats[line[0]]["telo_end"]:
                        continue
                    elif mod_pos > telo_stats[line[0]]["telo_start"]:
                        mod_dict[line[0]]["telo"] += 1
                    else:
                        mod_dict[line[0]]["subtelo"] += 1

                    mod_dict[line[0]]["full_read"] += 1
                    mod_dict[line[0]]["pos"].append(mod_pos)
                
    with open(summary_file, "w") as summary_fh:
        summary_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("read_id", "strand", "telo_length", "read_type", "possible_mods", "mods", "proportion"))
        for read in telo_stats:
            if read not in mod_dict:
                summary_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(read, telo_stats[read]["strand"], telo_stats[read]["telo_length"], "full_read", telo_possible[read]["full_read"], 0, 0))
                summary_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(read, telo_stats[read]["strand"], telo_stats[read]["telo_length"], "subtelo", telo_possible[read]["subtelo"], 0, 0))
                summary_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(read, telo_stats[read]["strand"], telo_stats[read]["telo_length"], "telo", telo_possible[read]["telo"], 0 ,0))
                continue
            if telo_possible[read]["telo"] == 0:
                summary_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(read, telo_stats[read]["strand"], telo_stats[read]["telo_length"], "telo", 0, 0, 0))
            else:
                summary_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(read, telo_stats[read]["strand"], telo_stats[read]["telo_length"], "telo", 
                                                                 telo_possible[read]["telo"], mod_dict[read]["telo"],
                                                                 mod_dict[read]["telo"]/telo_possible[read]["telo"] * 100))
            
            if telo_possible[read]["subtelo"] == 0:
                summary_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(read, telo_stats[read]["strand"], telo_stats[read]["telo_length"], "subtelo", 0, 0, 0))
            else:
                summary_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(read, telo_stats[read]["strand"], telo_stats[read]["telo_length"], "subtelo", 
                                                                   telo_possible[read]["subtelo"], mod_dict[read]["subtelo"],
                                                                 mod_dict[read]["subtelo"]/telo_possible[read]["subtelo"] * 100))
                
            if telo_possible[read]["full_read"] == 0:
                summary_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(read, telo_stats[read]["strand"], telo_stats[read]["telo_length"], "full_read", 0, 0, 0))
            else:
                summary_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(read, telo_stats[read]["strand"], telo_stats[read]["telo_length"], "full_read", 
                                                                 telo_possible[read]["full_read"], mod_dict[read]["full_read"],
                                                                 mod_dict[read]["full_read"]/telo_possible[read]["full_read"] * 100))

            
            
    max_telo_length = max([telo_stats[i]["telo_length"] for i in telo_stats])
    subtelo_stretch_max = max_subtelo_stretch
    seq_length = subtelo_stretch_max + max_telo_length

    sorted_telos = sorted(telo_stats, key=lambda k: telo_stats[k]["telo_length"], reverse=False)

    strand_spacer = 50
    y_offset_top = 50
    y_offset_bottom = 50
    x_offset_left = 75
    x_offset_right = 10
    seq_height = (image_height - strand_spacer - y_offset_top - y_offset_bottom)/len(sorted_telos)
    nucl_width = (image_width - x_offset_left - x_offset_right)/seq_length
    
    g_strands = len([read for read in telo_stats if telo_stats[read]["strand"]=="G"])
    c_strands = len([read for read in telo_stats if telo_stats[read]["strand"]=="C"])
    
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
        # draw y axis
        ctx.rectangle(x_offset_left - 3, y_offset_top+g_strands*seq_height+strand_spacer, 3, c_strands*seq_height + 3)
        # draw x axis
        ctx.rectangle(x_offset_left, image_height-y_offset_bottom, image_width-x_offset_left - x_offset_right, 3)
        ctx.fill()
        
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
        for i in range(0, c_strands, 200):
            if i % 1000 == 0:
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
    
            telo_seq = telo_orientated[telo][telo_stats[telo]["telo_start"]:telo_stats[telo]["telo_end"]]
            
            subtelo = telo_orientated[telo][0:telo_stats[telo]["telo_start"]][-5000:]
            x_rect = x_offset_left + subtelo_stretch_max*nucl_width
            i = 0
            while i < len(telo_seq):
                if telo_seq[i:i+6] == 'GGTTAG':
                    i += 6
                    ctx.set_source_rgba(0.118, 0.533, 0.898, 0.5)
                    ctx.rectangle(x_rect, y_rect, 6*nucl_width, seq_height)
                    ctx.fill()
                    x_rect += 6*nucl_width
                else:
                    ctx.set_source_rgba(1, 0.757, 0.027, 0.5)
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
                    ctx.set_source_rgba(0.118, 0.533, 0.898, 0.5)
                    ctx.rectangle(x_rect-6*nucl_width, y_rect, 6*nucl_width, seq_height)
                    ctx.fill()
                    x_rect -= 6*nucl_width
                else:
                    i += 1
                    ctx.set_source_rgba(1, 0.757, 0.027, 0.5)
                    ctx.rectangle(x_rect-1*nucl_width, y_rect, 1*nucl_width, seq_height)
                    ctx.fill()
                    x_rect -= 1*nucl_width     
    
            # draw modifications
            
            if telo in mod_dict:
                for pos in mod_dict[telo]["pos"]:
                    if pos < telo_stats[telo]["telo_start"] - subtelo_stretch_max:
                        continue
                    pos = pos - (telo_stats[telo]["telo_start"] - subtelo_stretch_max)
                    
                    ctx.set_source_rgb(0, 0.302, 0.251)
                    ctx.rectangle(x_offset_left+pos*nucl_width, y_rect, 1*nucl_width, seq_height)
                    ctx.fill()
                
            y_rect -= seq_height      
        # draw G telomeres
        ctx.set_source_rgb(0,0,0)
        # draw y axis
        ctx.rectangle(x_offset_left - 3, y_offset_top, 3, g_strands*seq_height + 3)
        ctx.fill()
        
        # label y-axis positions
        for i in range(0, g_strands, 200):
            if i % 1000 == 0:
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
            
            telo_seq = telo_orientated[telo][telo_stats[telo]["telo_start"]:telo_stats[telo]["telo_end"]]
            subtelo = telo_orientated[telo][0:telo_stats[telo]["telo_start"]][-5000:]
            x_rect = x_offset_left + subtelo_stretch_max*nucl_width
            i = 0
            while i < len(telo_seq):
                if telo_seq[i:i+6] == 'GGTTAG':
                    i += 6
                    ctx.set_source_rgba(0.118, 0.533, 0.898, 0.5)
                    ctx.rectangle(x_rect, y_rect, 6*nucl_width, seq_height)
                    ctx.fill()
                    x_rect += 6*nucl_width
                else:
                    ctx.set_source_rgba(1, 0.757, 0.027, 0.5)
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
                    ctx.set_source_rgba(0.118, 0.533, 0.898, 0.5)
                    ctx.rectangle(x_rect-6*nucl_width, y_rect, 6*nucl_width, seq_height)
                    ctx.fill()
                    x_rect -= 6*nucl_width
                else:
                    i += 1
                    ctx.set_source_rgba(1, 0.757, 0.027, 0.5)
                    ctx.rectangle(x_rect-1*nucl_width, y_rect, 1*nucl_width, seq_height)
                    ctx.fill()
                    x_rect -= 1*nucl_width
                    
            if telo in mod_dict:
                for pos in mod_dict[telo]["pos"]:
                    if pos < telo_stats[telo]["telo_start"] - subtelo_stretch_max:
                        continue
                    pos = pos - (telo_stats[telo]["telo_start"] - subtelo_stretch_max)
                    
                    ctx.set_source_rgb(0, 0.302, 0.251)
                    ctx.rectangle(x_offset_left+pos*nucl_width, y_rect, 1*nucl_width, seq_height)
                    ctx.fill()
            y_rect -= seq_height            

        
def main(args):

    args.max_subtelo_stretch = int(args.max_subtelo_stretch)
    args.image_width = int(args.image_width)
    args.image_height = int(args.image_height)

    plot_cluster_independent_telomeric_reads(args.reference_aln, \
                                            args.telo_stats, \
                                            args.mod_calls, \
                                            args.modified_nucleotide, \
                                            args.modification_code, \
                                            args.telomere_plot, \
                                            args.summary_file, \
                                            args.max_subtelo_stretch, \
                                            image_width=args.image_width, image_height=args.image_height)



def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference_aln",required=True)
    parser.add_argument("--telo_stats", required=True)
    parser.add_argument("--mod_calls", required=True)
    parser.add_argument('--modified_nucleotide', required=True)
    parser.add_argument('--modification_code', required=True)
    parser.add_argument("--telomere_plot", required=True)
    parser.add_argument("--summary_file", required=True)
    parser.add_argument("--max_subtelo_stretch", required=True)
    parser.add_argument("--image_width", required=True)
    parser.add_argument("--image_height", required=True)
    
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)