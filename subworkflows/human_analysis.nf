/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { EXTRACTING_MODCALLS } from "../bin/human/process.nf"
// include { EXTRACTING_MODCALLS as EXTRACT_TELO_MODS } from "../bin/human/process.nf"
include { CONVERT_BAM_TO_FASTQ } from "../bin/human/process.nf"
include { ALIGN_TO_REF as ALIGN_TO_REF } from "../bin/human/process.nf"
include { ALIGN_TO_REF as ALIGN_TO_SPIKE_IN } from "../bin/human/process.nf"
include { GENOMIC_READS_WITH_SPIKE } from "../bin/human/process.nf"
include { GENOMIC_READ_ANALYSIS } from "../bin/human/process.nf"
include { EXTRACT_TELO_MODS } from "../bin/human/process.nf"
include { EXTRACT_TELO_READS } from "../bin/human/process.nf"
include { TELO_READ_ANALYSIS } from "../bin/human/process.nf"
//include { TELO_CLUSTER_ANALYSIS } from "../bin/human/process.nf"
include { PLOT_RESULTS; PLOT_TELO_RESULTS } from "../bin/human/process.nf"
include { READ_BY_READ_ANALYSIS } from "../bin/human/process.nf"
include { CREATE_CLUSTER_FASTA } from "../bin/human/process.nf"
include { REORIENT_TELOBAM } from "../bin/human/process.nf"
include { REPAIR_TELOBAM } from "../bin/human/process.nf"
include { GENERATE_CONSENSUS_AND_ALIGN } from "../bin/human/process.nf"
include { REPAIR; INDEX; PILEUP; PLOT_PILEUP_RESULTS } from "../bin/human/process.nf"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow human_analysis {
    
    main:
        Pinguscript.ping_start(nextflow, workflow, params)

        if (params.genomic_comparison) {
            // modkit extract modifications calls
            modcalls = EXTRACTING_MODCALLS(params.modbam)

            // convert bam to fastq returning modcalls
            fastq = CONVERT_BAM_TO_FASTQ(params.modbam)

            // align fastq to reference
            aln_to_ref = ALIGN_TO_REF(fastq.mod_fastq, params.reference)

            //if spike in sample as well align to spike in sample and process reads
            if (params.spike_in_reference != ""){
                spike_aln = ALIGN_TO_SPIKE_IN(fastq.mod_fastq, params.spike_in_reference)
                genomic_results = GENOMIC_READS_WITH_SPIKE(aln_to_ref.alignment, spike_aln.alignment,
                                                            params.telo_stats, modcalls.modcalls)
            }
            // // or just process reads
            else {
                genomic_results = process_genomic_reads(aln_to_ref, modcalls)
            }
        }

        telo_modbam = EXTRACT_TELO_READS(params.modbam, params.telo_stats)

        // reorient tarpon telobam
        telo_bam = REORIENT_TELOBAM(params.telo_bam, params.telo_stats, telo_modbam.telo_modbam)

        // repair reoriented tarpon telo_bam using extracted telo reads
        telo_modbam_tarpon = REPAIR_TELOBAM(telo_modbam.telo_modbam, telo_bam.telo_bam)

        //extract telomeric relevant mods
        telo_modcalls = EXTRACT_TELO_MODS(telo_modbam_tarpon.telo_modbam)

        // process telomeric modifications
        telo_mods = TELO_READ_ANALYSIS(telo_modbam_tarpon.telo_modbam, params.telo_stats, telo_modcalls.modcalls)

        // process telomeric modifications in cluster specific manner
        if (params.clustering) {
            // cluster_mods = TELO_CLUSTER_ANALYSIS(telo_modbam_tarpon.telo_modbam,
            //                                     params.telo_stats, telo_modcalls.modcalls)
                                        
            read_by_read_results = READ_BY_READ_ANALYSIS(telo_modbam_tarpon.telo_modbam,
                                                params.telo_stats, telo_modcalls.modcalls, \
                                                telo_mods.telomeric_reads)

            cluster_fasta = CREATE_CLUSTER_FASTA(telo_modbam_tarpon.telo_modbam,//repair_results,
                                                params.telo_stats)

            cluster_fasta.fasta
                .flatten()
                .map { file -> tuple(file.baseName, file) }
                .set { cluster_fasta_ch }

            cluster_fasta.bam  
                .flatten()
                .map { file -> tuple(file.baseName, file) }
                .set { cluster_bam_ch }

    
            joined_channel = cluster_fasta_ch.join(cluster_bam_ch)
          
            consensus_alignment = GENERATE_CONSENSUS_AND_ALIGN(joined_channel)

            repaired_bam = REPAIR(consensus_alignment.cluster_files, telo_modbam.telo_modbam)
            indexed_bam = INDEX(repaired_bam.repaired)
            pileup = PILEUP(indexed_bam.indexed)

            PLOT_PILEUP_RESULTS(pileup.results)
        }

        if (params.genomic_comparison){
            PLOT_RESULTS(genomic_results.genomic_reads,telo_mods.telomeric_reads, params.telo_stats)
        }
        else {
            PLOT_TELO_RESULTS(telo_mods.telomeric_reads, params.telo_stats)
        }
}
