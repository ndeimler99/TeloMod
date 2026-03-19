/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { EXTRACTING_MODCALLS } from "../bin/process.nf"
// include { EXTRACTING_MODCALLS as EXTRACT_TELO_MODS } from "../bin/process.nf"
include { CONVERT_BAM_TO_FASTQ } from "../bin/process.nf"
include { ALIGN_TO_REF as ALIGN_TO_REF } from "../bin/process.nf"
include { ALIGN_TO_REF as ALIGN_TO_SPIKE_IN } from "../bin/process.nf"
include { GENOMIC_READS_WITH_SPIKE } from "../bin/process.nf"
include { GENOMIC_READ_ANALYSIS } from "../bin/process.nf"
include { EXTRACT_TELO_MODS } from "../bin/process.nf"
include { EXTRACT_TELO_READS } from "../bin/process.nf"
include { TELO_READ_ANALYSIS } from "../bin/process.nf"
include { TELO_CLUSTER_ANALYSIS } from "../bin/process.nf"
include { PLOT_RESULTS } from "../bin/process.nf"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow human_analysis {
    
    main:
        Pinguscript.ping_start(nextflow, workflow, params)

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

        telo_bam = EXTRACT_TELO_READS(params.modbam, params.telo_stats)

        //extract telomeric relevant mods
        telo_modcalls = EXTRACT_TELO_MODS(telo_bam.telo_modbam)


        // process telomeric modifications
        telo_mods = TELO_READ_ANALYSIS(telo_bam.telo_modbam, params.telo_stats, telo_modcalls.modcalls)

        // process telomeric modifications in cluster specific manner
        if (params.cluster_results != "") {
            cluster_mods = TELO_CLUSTER_ANALYSIS(params.cluster_results, telo_bam.telo_modbam,
                                                params.telo_stats, telo_modcalls.modcalls, \
                                                telo_mods.telomeric_reads)
        }

        PLOT_RESULTS(genomic_results.genomic_reads,telo_mods.telomeric_reads, params.telo_stats)
}
