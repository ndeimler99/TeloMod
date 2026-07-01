# TeloMod v1.1.2

TeloMod is a nextflow pipeline for the analysis of Telomeric Reads sequenced with Oxford Nanopore Technologies (ONT) sequencing devices and basecalled with a modified basecalling model.
TeloMod currently only works on human telomere sequencing data (and other vertebrates), but will in the future be expanded for other model species. The use of nextflow eliminates any dependency/version conflicts as the first time the pipeline is executed docker files containing the appropriate software will be downloaded. Please see the nextflow and docker documentation for more information on this.

## Input

TeloMod requires the following paramters as input: 

| Parameter      | Default | Description |  
| :---:        |    :----:   | :--- |
| --human | False | Is the data being analyzed human data? |
| --pombe| False | Is the data being analyzed S. pombe data? (currently does nothing) |
| --cerevisiae | False | Is the data being analyzed S. cerevisiae data? (currently does nothing) |
| --outdir | Required | The location of the directory for output files to be generated |
| --modbam | Required | The modified basecalling bam file |
| --telo_stats | Required | TARPON output (telo_stats.txt) file |
| --telo_bam | Required | TARPON output (telomeric.bam) file |
| --cluster_results | "" | Not required unless interested in cluster specific analyses. Telogator2 tlens_by_allele.tsv output file supplied by TARPON |
| --genomic_comparison | false | If set, will compare telomeric sequences to background genomic reads contained within the provided modbam by aligning to the provided reference genome |
| --reference | "" | Reference genome of species of interest to determine background levels of genomic modifications, required if genomic_comparison is set |
| --spike_in_reference | "" | Not required unless the telomeres were sequenced in parallel to another species. This allows for the removal of off-target non-telomeric background DNA used in the analysis|
| --mod_confidence | 0.9 | The minimum probability of the modified base calling | 
| --modification_code | m | The ONT modification code (5mC = m, 5hmC = h)|
| --modification_nucl | C | The nucleotide that is expected to be modified. For 5mC place C, for 5mCG place CG |
| --minimum_genomic_read_length | 5000 | The minimum read length required to be included for background information |
| --image_width | 1500 | Pixel width of telomere sequence plots generated |
| --image_height | 1000| Pixel height of telomere sequence plots generated |
| --max_subtelo_stretch | 5000 | The number of subtelomeric nucleotides to include in telomere sequence plots |
| --threads | 4 | Number of threads to use in analysis |
| --c_strand_only | false | If set, will only analyze C strand telomeric sequences. This parameter should only be set when telomeric sequences are collected through the splint-enriched protocols |


## Workflow

TeloMod currently only analyzes human (or other vertebrate) telomeric sequencing data and does so at a read by read level calculating the proportion of modified nucleotides per read as well as a pileup analysis if clustering results are provided. The TeloMod workflow can be found below.
![TeloMod Workflow](telomod_workflow.pdf)

To execute TeloMod basecalling with a modification specific model must be performed and TARPON must be executed on the results. Modifications are first added to the telomeric.bam file produced by TARPON via modkit repair. From there, modification information is extracted and telomeres are analyzed on a read by read basis segregating the subtelomere from the telomere. Modifications are required to pass the minimum threshold specified in the parameters (default 0.9). 

If --genomic_comparison is set all reads are additionally aligned to the provided reference genome and if provided the spike_in_reference. Primary alignments that do not align to the spike in reference, are not present within the TARPON telomeric.bam file, and are greater than --minimum_genomic_read_length are extracted and analyzed as a background level of modification.

If clustering information is provided, read-by-read modification statistics are determined in a cluster specific manner. Additionally, every telomeric sequence within a given cluster is aligned to each other via a SPOA local alignment. Telomeric sequences are then aligned back to the generated consensus sequence and modification information is plotted in a positional manner. The results of this can be found in the clusters subdirectory of the specified output folder.

A variety of plots are generated comparing genomic to telomeric reads as well as portions of telomeric reads to each other including clustery by cluster comparisons and correlations with telomere length.
Furthermore, both genomic and telomeric summary files are provided with read by read statistics on the number of modified bases identified.
