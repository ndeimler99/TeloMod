# TeloMod

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
| --cluster_results | "" | Not required unless interested in cluster specific analyses. Telogator2 tlens_by_allele.tsv output file supplied by TARPON |
| --reference | Required | Reference genome of species of interest to determine background levels of genomic modifications |
| --spike_in_reference | "" | Not required unless the telomeres were sequenced in parallel to another species. This allows for the removal of off-target non-telomeric background DNA used in the analysis|
| --mod_confidence | 0.9 | The minimum probability of the modified base calling | 
| --modification_code | m | The ONT modification code (5mC = m, 5hmC = h)|
| --modification_nucl | C | The nucleotide that is expected to be modified. For 5mC place C, for 5mCG place CG |
| --minimum_genomic_read_length | 5000 | The minimum read length required to be included for background information |
| --image_width | 1500 | Pixel width of telomere sequence plots generated |
| --image_height | 1000| Pixel height of telomere sequence plots generated |
| --max_subtelo_stretch | 5000 | The number of subtelomeric nucleotides to include in telomere sequence plots |
| --threads | 4 | Number of threads to use in analysis |

## Workflow

TeloMod currently only analyzes human (or other vertebrate) telomeric sequencing data and does so at a read by read level calculating the proportion of modified nucleotides per read. In the future, pileup analyses will also be offered.
Modifications are first extracted from the supplied modbam based on passing the minimum threshold specified (by default 0.9).
All sequenced reads are aligned to both the reference genome and the spike_in_reference if provided. Primary alignments to the reference that do not align to the spike in and are no telomeric are analyzed for background levels of modifications within the genome.
The proportion of nucleotides modified within telomere containing reads is calculated, as well as segregated by subtelomere and telomeric portions of the read, and further segregated by strand.
All telomere containing reads as well as their subtelomere stretches are plotted with the location of modifications labeled. This is additionally performed for each individual cluster if a clustering results file is provided.
A variety of plots are generated comparing genomic to telomeric reads as well as portions of telomeric reads to each other including clustery by cluster comparisons and correlations with telomere length.
Furthermore, both genomic and telomeric summary files are provided with read by read statistics on the number of modified bases identified.

