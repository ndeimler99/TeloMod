#!/usr/bin/env python3

import argparse
import pysam
import gzip
import cairo
import math
from Bio import SeqIO
import regex as re

def rev_comp(seq):
    rev_dict = {"A":"T", "T":"A", "C":"G", "G":"C"}
    return ''.join([rev_dict[i] for i in seq[::-1]])

def get_telo_start(read, repeat, sliding_window, sliding_window_interval, upper_threshold, lower_threshold, consecutive_threshold):
    # identifies start of telomeric read (read).  See manuscript for more details on how this is done and why it is done each way
    telo_found = False
    telo_start = None
    below_threshold = 0
    for i in range(0, len(read)-sliding_window, sliding_window_interval):
        telo_matches = list(re.finditer(r"(%s){s<=1}" % repeat, read[i:i+sliding_window]))
        telo_perc = len(telo_matches) * len(repeat) / sliding_window

        if telo_perc >= upper_threshold:
            below_threshold = 0
            if not telo_found:
                telo_found = True
                telo_start = i + telo_matches[0].span()[0]
        elif telo_perc < lower_threshold:
            if telo_found:
                below_threshold += 1
                if below_threshold >= consecutive_threshold:
                    telo_found = False
                    telo_start = None
    return telo_found, telo_start


def plot_pileup_results(consenus_fa, g_strand_pileup, c_strand_pileup, g_cov, c_cov, file_name, 
                        modification, modified_nucleotide, c_strand_only, image_height, image_width):
    
    consensus_seq = list(SeqIO.parse(consenus_fa, "fasta"))[0]
    
    nucleotide_frequency = {"G":{}, "C":{}}
    for i in range(0,len(consensus_seq),30):
        if not c_strand_only:
            nucleotide_frequency["G"][i] = (consensus_seq.seq[i:i+30].count(modified_nucleotide) * len(modified_nucleotide))/len(consensus_seq.seq[i:i+30])
        nucleotide_frequency["C"][i] = (rev_comp(consensus_seq.seq[i:i+30]).count(modified_nucleotide)*len(modified_nucleotide))/len(consensus_seq.seq[i:i+30])
        
        
    telo_start = get_telo_start(str(consensus_seq.seq), "GGTTAG", 100, 15, 0.6, 0.05, 15)
    # print(telo_start)
    # print(len(consensus_seq.seq))
    coverage_stats = {"G":{}, "C":{}}
    if not c_strand_only:
        with open(g_cov, "r") as fh:
            for line in fh:
                line = line.strip().split()
                coverage_stats["G"][int(line[1])-1] = int(line[2])

    with open(c_cov, "r") as fh:
        for line in fh:
            line = line.strip().split()
            coverage_stats["C"][int(line[1])-1] = int(line[2])
            

    # print(coverage_stats)            
    #print(len(telo_stats))
    pileup_stats = {"G":{}, "C":{}}
    # print(modification)
    if not c_strand_only:
        with open(g_strand_pileup, "r") as fh:
            for line in fh:
                line = line.strip().split()
                #print(line)
                if line[3] == modification:
                    print(line)
                    pileup_stats["G"][int(line[1])] = [float(line[11])/coverage_stats["G"][int(line[1])] * 100, int(line[11])]
    
    with open(c_strand_pileup, "r") as fh:
        for line in fh:
            line = line.strip().split()
            if line[3] == modification:
                pileup_stats["C"][int(line[1])] = [float(line[11])/coverage_stats["C"][int(line[1])]* 100, int(line[11])]

    # print(pileup_stats)        
    telo_length = len(consensus_seq.seq) - telo_start[1]
    #print(max_telo_length)
    seq_length = len(consensus_seq.seq)
    

    spacer = 8
    strand_spacer = 20
    y_offset_top = 30
    y_offset_bottom = 30
    x_offset_left = 75
    x_offset_right = 10
    if not c_strand_only:
        sub_height = image_height/12
        plot_height = (image_height - y_offset_top - y_offset_bottom - spacer*4 - strand_spacer - sub_height*4 - 25)/2

    else:
        sub_height = image_height/6
        plot_height=(image_height - y_offset_top - y_offset_bottom - spacer*2- strand_spacer - sub_height*2 - 25)
    
#     seq_height = (image_height - strand_spacer - y_offset_top - y_offset_bottom)/len(sorted_telos)
    nucl_width = (image_width - x_offset_left - x_offset_right)/seq_length
    
#     g_strands = len([read for read in telo_stats if telo_stats[read]["strand"]=="G"])
#     c_strands = len([read for read in telo_stats if telo_stats[read]["strand"]=="C"])
    
#     g_start = int(np.median([telo_stats[telo]["telo_start"] for telo in sorted_telos if telo_stats[telo]["strand"]=="G"]))  
#     c_start = int(np.median([telo_stats[telo]["telo_start"] for telo in sorted_telos if telo_stats[telo]["strand"]=="C"]))

    with cairo.PDFSurface(file_name, image_width, image_height) as surface:
#     #if True:
#         surface = cairo.ImageSurface(cairo.FORMAT_RGB24,
# #                                  image_width,
# #                                  image_height)
        ctx = cairo.Context(surface)
        ctx.rectangle(0, 0, image_width, image_height)
        ctx.set_source_rgb(1,1,1)
        ctx.fill()
        
        ctx.set_font_size(9)
        ctx.select_font_face("Ariel",
                              cairo.FONT_SLANT_NORMAL,
                              cairo.FONT_WEIGHT_NORMAL)
        
        #### x axis bar
        ctx.set_source_rgb(0,0,0)
        ctx.rectangle(x_offset_left, image_height-y_offset_bottom, image_width-x_offset_left - x_offset_right, 3)
        ctx.fill()
        
        ### C stand modification percentages
        ctx.rectangle(x_offset_left-3, image_height-y_offset_bottom-30-spacer-plot_height, 3, plot_height)
        ctx.fill()
        
        for i in range(0, 101, 25):
            a,b,width,height,c,d = ctx.text_extents('{}'.format(i))
            ctx.move_to(x_offset_left - width - 8, image_height - y_offset_bottom - 30 - spacer - i/100*plot_height + (height+1)/2)
            ctx.text_path('{}'.format(i))
            ctx.rectangle(x_offset_left-5, image_height - y_offset_bottom -30 - spacer - i/100*plot_height, 4, 2)
            ctx.fill()

        ctx.save()
        a,b,width,height,c,d = ctx.text_extents('Percentage')
        ctx.translate(x_offset_left-35, image_height - y_offset_bottom - 30 - spacer - plot_height/2 + width/2)
        ctx.rotate(-math.pi / 2)   
        ctx.move_to(0, 0)
        ctx.text_path('Percentage')
        ctx.fill()
        ctx.restore()
        ctx.save()
        
        #### C strand coverage
        ctx.rectangle(x_offset_left-3, image_height-y_offset_bottom-30-spacer*2-plot_height-sub_height, 3, sub_height)
        ctx.fill()
        
        c_max = max(coverage_stats["C"].values())
        for i in range(0, c_max + 1, int(c_max/4)+1):
            a,b,width,height,c,d = ctx.text_extents('{}'.format(i))
            ctx.move_to(x_offset_left - width - 8, image_height - y_offset_bottom - 30 - spacer*2 - plot_height - i/c_max*sub_height + (height+1)/2)
            ctx.text_path('{}'.format(i))
            ctx.rectangle(x_offset_left-5, image_height - y_offset_bottom -30 - spacer*2 - plot_height - i/c_max*sub_height, 4, 2)
            ctx.fill()
            
        ctx.save()
        a,b,width,height,c,d = ctx.text_extents('Coverage')
        ctx.translate(x_offset_left-35, image_height - y_offset_bottom - 30 - spacer*2 - plot_height - sub_height/2 + width/2)
        ctx.rotate(-math.pi / 2)   
        ctx.move_to(0, 0)
        ctx.text_path('Coverage')
        ctx.fill()
        ctx.restore()
        ctx.save()
        
        ctx.save()
        a,b,width,height,c,d = ctx.text_extents('C Strand')
        ctx.translate(x_offset_left-60, image_height - y_offset_bottom - 30 - (spacer*2 + plot_height + sub_height*2)/2 + width/2)
        ctx.rotate(-math.pi / 2)   
        ctx.move_to(0, 0)
        ctx.text_path('C Strand')
        ctx.fill()
        ctx.restore()
        ctx.save()
        
        #### C strand Density
        ctx.rectangle(x_offset_left-3, image_height-y_offset_bottom-30-spacer*3-plot_height-sub_height*2, 3, sub_height)
        ctx.fill()
        
        for i in range(0, 101, 50):
            a,b,width,height,c,d = ctx.text_extents('{}'.format(i))
            ctx.move_to(x_offset_left - width - 8, image_height - y_offset_bottom - 30 - spacer*3 - plot_height - sub_height - i/100*sub_height + (height+1)/2)
            ctx.text_path('{}'.format(i))
            ctx.rectangle(x_offset_left-5, image_height - y_offset_bottom -30 - spacer*3 - plot_height - sub_height - i/100*sub_height, 4, 2)
            ctx.fill()
        
        ctx.save()
        a,b,width,height,c,d = ctx.text_extents('Density')
        ctx.translate(x_offset_left-35, image_height - y_offset_bottom - 30 - spacer*3 - plot_height - sub_height - sub_height/2 + width/2)
        ctx.rotate(-math.pi / 2)   
        ctx.move_to(0, 0)
        ctx.text_path('Density')
        ctx.fill()
        ctx.restore()
        ctx.save()
        
        if not c_strand_only:


            ctx.rectangle(x_offset_left-3, image_height-y_offset_bottom-30-spacer*3-plot_height*2-sub_height*2-strand_spacer, 3, plot_height)
            ctx.fill()
        
            for i in range(0, 101, 25):
                a,b,width,height,c,d = ctx.text_extents('{}'.format(i))
                ctx.move_to(x_offset_left - width - 8, image_height - y_offset_bottom - strand_spacer - 30 - spacer*3 - plot_height - sub_height*2 - i/100*plot_height + (height+1)/2)
                ctx.text_path('{}'.format(i))
                ctx.rectangle(x_offset_left-5, image_height - y_offset_bottom - strand_spacer - 30 - spacer*3 - plot_height - sub_height*2 - i/100*plot_height, 4, 2)
                ctx.fill()
                
            ctx.save()
            a,b,width,height,c,d = ctx.text_extents('Percentage')
            ctx.translate(x_offset_left-35, image_height - y_offset_bottom - 30 - strand_spacer - sub_height*2 - spacer*3 - plot_height - plot_height/2 + width/2)
            ctx.rotate(-math.pi / 2)   
            ctx.move_to(0, 0)
            ctx.text_path('Percentage')
            ctx.fill()
            ctx.restore()
            ctx.save()
            
            ctx.rectangle(x_offset_left-3, image_height-y_offset_bottom-30-spacer*4-plot_height*2-sub_height*3-strand_spacer, 3, sub_height)
            ctx.fill()
        
            g_max = max(coverage_stats["G"].values())
            for i in range(0, g_max + 1, int(g_max/4)+1):
                a,b,width,height,c,d = ctx.text_extents('{}'.format(i))
                ctx.move_to(x_offset_left - width - 8, image_height - y_offset_bottom - strand_spacer - 30 - spacer*4 - plot_height*2 - sub_height*2 - i/g_max*sub_height + (height+1)/2)
                ctx.text_path('{}'.format(i))
                ctx.rectangle(x_offset_left-5, image_height - y_offset_bottom - strand_spacer - 30 - spacer*4 - plot_height*2 - sub_height*2 - i/g_max*sub_height, 4, 2)
                ctx.fill()
            
            ctx.save()
            a,b,width,height,c,d = ctx.text_extents('Coverage')
            ctx.translate(x_offset_left-35, image_height - y_offset_bottom - 30 - strand_spacer - sub_height*2 - spacer*4 - plot_height*2 - sub_height/2 + width/2)
            ctx.rotate(-math.pi / 2)   
            ctx.move_to(0, 0)
            ctx.text_path('Coverage')
            ctx.fill()
            ctx.restore()
            ctx.save()
        
        
            ctx.save()
            a,b,width,height,c,d = ctx.text_extents('G Strand')
            ctx.translate(x_offset_left-60, image_height - y_offset_bottom - 30 - strand_spacer - sub_height*2 - spacer*2 - plot_height - (plot_height+sub_height*2+spacer*2)/2 + width/2)
            ctx.rotate(-math.pi / 2)   
            ctx.move_to(0, 0)
            ctx.text_path('G Strand')
            ctx.fill()
            ctx.restore()
            ctx.save()
            
            ctx.rectangle(x_offset_left-3, image_height-y_offset_bottom-30-spacer*5-plot_height*2-sub_height*4-strand_spacer, 3, sub_height)
            ctx.fill()
            
            for i in range(0, 101, 50):
                a,b,width,height,c,d = ctx.text_extents('{}'.format(i))
                ctx.move_to(x_offset_left - width - 8, image_height - y_offset_bottom - strand_spacer - 30 - spacer*5 - plot_height*2 - sub_height*3 - i/100*sub_height + (height+1)/2)
                ctx.text_path('{}'.format(i))
                ctx.rectangle(x_offset_left-5, image_height - y_offset_bottom - strand_spacer - 30 - spacer*5 - plot_height*2 - sub_height*3 - i/100*sub_height, 4, 2)
                ctx.fill()
            
            ctx.save()
            a,b,width,height,c,d = ctx.text_extents('Density')
            ctx.translate(x_offset_left-35, image_height - y_offset_bottom - 30 - strand_spacer - sub_height*3 - spacer*5 - plot_height*2 - sub_height/2 + width/2)
            ctx.rotate(-math.pi / 2)   
            ctx.move_to(0, 0)
            ctx.text_path('Density')
            ctx.fill()
            ctx.restore()
            ctx.save()
        
        
        ctx.set_line_width(0.5)
        max_cov = max([coverage_stats["C"][pos] for pos in coverage_stats["C"]])
        for pos in pileup_stats["C"]:
            if pileup_stats["C"][pos][0] > 0:
                ctx.rectangle(x_offset_left + pos*nucl_width, image_height-y_offset_bottom-30-spacer-plot_height*(pileup_stats["C"][pos][0]/100), nucl_width, plot_height*pileup_stats["C"][pos][0]/100)
                ctx.set_source_rgb(1, 0, 0.365)
                ctx.fill()
                ctx.rectangle(x_offset_left + pos*nucl_width, image_height-y_offset_bottom-30-spacer-plot_height, nucl_width, plot_height*(1-pileup_stats["C"][pos][0]/100))
                ctx.set_source_rgba(0.5, 0.5, 0.5, 0.1)
                ctx.fill()
                
        
        positions = sorted(coverage_stats["C"].keys())
        ctx.set_source_rgb(0,0,0)
        for i,pos in enumerate(positions):
            x = pos*nucl_width + x_offset_left
            y = image_height-y_offset_bottom-30-spacer*2-plot_height-sub_height*(coverage_stats["C"][pos]/max_cov)
            if i == 0:
                ctx.move_to(x,y)
            else:
                ctx.line_to(x,y)
        ctx.stroke()
        
        
        positions = sorted(nucleotide_frequency["C"].keys())
        ctx.set_source_rgb(0,0,0)
        for i,pos in enumerate(positions):
            x = pos*nucl_width + x_offset_left
            y = image_height-y_offset_bottom-30-spacer*3-plot_height-sub_height-sub_height*nucleotide_frequency["C"][pos]
            if i == 0:
                ctx.move_to(x,y)
            else:
                ctx.line_to(x,y)
        ctx.stroke()
        
        
        if not c_strand_only:
            max_cov = max([coverage_stats["G"][pos] for pos in coverage_stats["G"]])
            for pos in pileup_stats["G"]:
                if pileup_stats["G"][pos][0] > 0:
                    ctx.rectangle(x_offset_left + pos*nucl_width, image_height-y_offset_bottom-30-spacer*3-plot_height-sub_height*2-strand_spacer-plot_height*(pileup_stats["G"][pos][0]/100), nucl_width, plot_height*pileup_stats["G"][pos][0]/100)
                    ctx.set_source_rgb(1, 0, 0.365)
                    ctx.fill()
                    ctx.rectangle(x_offset_left + pos*nucl_width, image_height-y_offset_bottom-30-spacer*3-plot_height*2-sub_height*2-strand_spacer, nucl_width, plot_height*(1-pileup_stats["G"][pos][0]/100))
                    ctx.set_source_rgba(0.5, 0.5, 0.5, 0.1)
                    ctx.fill()
                    
            
            positions = sorted(coverage_stats["G"].keys())
            ctx.set_source_rgb(0,0,0)
            for i,pos in enumerate(positions):
                x = pos*nucl_width + x_offset_left
                y = image_height-y_offset_bottom-30-spacer*4-plot_height*2-sub_height*2-strand_spacer-sub_height*(coverage_stats["G"][pos]/max_cov)
                if i == 0:
                    ctx.move_to(x,y)
                else:
                    ctx.line_to(x,y)
            ctx.stroke()
            
            
            positions = sorted(nucleotide_frequency["G"].keys())
            ctx.set_source_rgb(0,0,0)
            ctx.set_line_width(0.5)
            for i,pos in enumerate(positions):
                x = pos*nucl_width + x_offset_left
                y = image_height-y_offset_bottom-30-spacer*5-plot_height*2-sub_height*3-strand_spacer-sub_height*nucleotide_frequency["G"][pos]
                if i == 0:
                    ctx.move_to(x,y)
                else:
                    ctx.line_to(x,y)
            ctx.stroke()
            

        
        #label x-axis positions
        ctx.set_source_rgb(0, 0, 0)
        for i in range(0, telo_length, 1000):
            if i % 3000 == 0:
                a,b,width,height,c,d = ctx.text_extents('{}'.format(i))
                #ctx.move_to(x_offset_left + i*nucl_width - width/2 + 1000*nucl_width - 2.5, image_height - y_offset_bottom + 25)
                ctx.move_to(x_offset_left + i*nucl_width - width/2 + telo_start[1]*nucl_width, image_height - y_offset_bottom + 13)
                ctx.text_path('{}'.format(i))
                ctx.fill()
            ctx.rectangle(x_offset_left + i*nucl_width -1 + telo_start[1]*nucl_width, image_height-y_offset_bottom, 2, 4)
            ctx.fill()
        
        for i in range(0, telo_start[1]+1, 1000):
            if i == 0:
                continue
            if i % 1000 == 0:
                a,b,width,height,c,d = ctx.text_extents('-{}'.format(i))
                #ctx.move_to(x_offset_left + i*nucl_width - width/2 + 1000*nucl_width - 2.5, image_height - y_offset_bottom + 25)
                ctx.move_to(x_offset_left + telo_start[1]*nucl_width - i*nucl_width - width/2, image_height - y_offset_bottom + 13)
                ctx.text_path('-{}'.format(i))
                ctx.fill()
            ctx.rectangle(x_offset_left + telo_start[1]*nucl_width - i*nucl_width - 1, image_height-y_offset_bottom, 2, 4)
            ctx.fill()
            
        
        a,b,width,height,c,d = ctx.text_extents('Telomere Position')
        x_coords = (image_width-x_offset_left-x_offset_right)/2 + x_offset_left - width
        ctx.move_to(x_coords, image_height-y_offset_bottom/3)
        ctx.text_path('Telomere Position')
        ctx.fill()
        
        a,b,width,height,c,d = ctx.text_extents('Consensus')
        x_coords = (x_offset_left - width - 3)
        ctx.move_to(x_coords, image_height-y_offset_bottom-15 + height/2)
        ctx.text_path('Consensus')
        ctx.fill()
            
#         # draw telomeres
        y_rect = image_height - y_offset_bottom - 30
     
        telo_seq = str(consensus_seq.seq)[telo_start[1]:]
        subtelo = str(consensus_seq.seq)[0:telo_start[1]]
        
        #print(len(subtelo))
        x_rect = x_offset_left + len(subtelo)*nucl_width
        i = 0

        #print(telo_seq)
        while i < len(telo_seq):
            if telo_seq[i:i+6] == 'GGTTAG':
                i += 6
                ctx.set_source_rgba(0.118, 0.533, 0.898, 0.5)
                ctx.rectangle(x_rect, y_rect, 6*nucl_width, 30)
                ctx.fill()
                x_rect += 6*nucl_width
            elif telo_seq[i:i+1] == '-':
                i += 1
                ctx.set_source_rgba(0.5, 0.5, 0.5, 1)
                ctx.rectangle(x_rect, y_rect, 1*nucl_width, 30)
                ctx.fill()
                x_rect += 1*nucl_width
            else:
                ctx.set_source_rgba(1, 0.757, 0.027, 0.5)
                ctx.rectangle(x_rect, y_rect, 1*nucl_width, 30)
                ctx.fill()
                x_rect += 1*nucl_width
                i+=1
    
            # draw subtelomere
        i = 0
        x_rect = len(subtelo)*nucl_width + x_offset_left
        while i < len(subtelo):
            if subtelo[-i:-i-6] == "GGTTAG":
                i += 6
                ctx.set_source_rgba(0.118, 0.533, 0.898, 0.5)
                ctx.rectangle(x_rect-6*nucl_width, y_rect, 6*nucl_width, 30)
                ctx.fill()
                x_rect -= 6*nucl_width
            elif subtelo[-i-1:-i] == "-":
                i += 1
                ctx.set_source_rgba(0.5, 0.5, 0.5, 0.5)
                ctx.rectangle(x_rect-1*nucl_width, y_rect, 1*nucl_width, 30)
                ctx.fill()
                x_rect -= 1*nucl_width
            else:
                #print(subtelo[-i:-(i)-2])
                i += 1
                ctx.set_source_rgba(1, 0.757, 0.027, 0.5)
                ctx.rectangle(x_rect-1*nucl_width, y_rect, 1*nucl_width, 30)
                ctx.fill()
                x_rect -= 1*nucl_width     

        
def main(args):

    
    #args.max_subtelo_stretch = int(args.max_subtelo_stretch)
    args.image_width = int(args.image_width)
    args.image_height = int(args.image_height)
    args.c_strand_only = args.c_strand_only == "true"

    plot_pileup_results(args.consensus, args.g_strand_pileup, args.c_strand_pileup, args.g_strand_cov, args.c_strand_cov, \
                        args.out_file, args.modification, args.modified_nucleotide, args.c_strand_only, args.image_height, args.image_width)




def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--consensus",required=True)
    parser.add_argument("--g_strand_pileup", required=True)
    parser.add_argument("--c_strand_pileup", required=True)
    parser.add_argument('--g_strand_cov', required=True)
    parser.add_argument('--c_strand_cov', required=True)
    parser.add_argument("--out_file", required=True)
    parser.add_argument("--modification", required=True)
    parser.add_argument("--modified_nucleotide", required=True)
    parser.add_argument("--image_width", required=True)
    parser.add_argument("--image_height", required=True)
    parser.add_argument("--c_strand_only", required=True)
    
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)