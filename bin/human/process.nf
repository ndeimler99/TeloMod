import groovy.json.JsonOutput

process getParams {

    label "telomod"
    tag "Getting Parameters"

    output:
        path "params.json", emit: params

    script:
    json_str = JsonOutput.toJson(params)
    json_indented = JsonOutput.prettyPrint(json_str)

    """
    echo '${json_indented}' > "params.json"
    """
}

process getVersions {

    label "telomod"
    tag "Getting Versions"

    output:
        path "versions.txt", emit: versions

    script:
    """
    python --version | sed 's/ /,/' >> versions.txt
    python -c "import argparse; print(f'argparse,{argparse.__version__}')" >> versions.txt
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    python -c "import numpy; print(f'numpy,{numpy.__version__}')" >> versions.txt
    """
}

process getManifest {
    
    label 'telomod'
    tag "Collecting Manifest Data"

    output:
        path "manifest.json", emit:manifest

    script:
    json_str = JsonOutput.toJson(workflow.manifest)
    json_indented = JsonOutput.prettyPrint(json_str)
    """
    echo '${json_indented}' > "manifest.json"
    """
}


process EXTRACTING_MODCALLS {

    tag "Extracting Modification Calls from ModBam"
    label 'modkit'
    cpus Math.min(params.threads as int, Runtime.runtime.availableProcessors())

    input:
        path(modbam)
    
    output:
        path("all_modcalls.txt"), emit: modcalls

    //publishDir "${params.outdir}/", mode: 'copy', overwrite:true, pattern:"all_modcalls.txt"

    script:
    """
    modkit extract calls --reference ${params.reference} \
            --bgzf --pass-only \
            --filter-threshold ${params.mod_confidence} ${modbam} all_modcalls.txt
    """
}

process EXTRACT_TELO_MODS {

    tag "Extracting Modification Calls from Telo ModBam"
    label 'modkit'
    cpus Math.min(params.threads as int, Runtime.runtime.availableProcessors())

    input:
        path(modbam)
    
    output:
        path("telo_modcalls.txt.gz"), emit: modcalls

    publishDir "${params.outdir}/", mode: 'copy', overwrite:true, pattern:"telo_modcalls.txt.gz"

    script:
    """
    modkit extract calls --bgzf --pass-only \
            --filter-threshold ${params.mod_confidence} ${modbam} telo_modcalls.txt.gz
    """
}


process CONVERT_BAM_TO_FASTQ {

    label 'telomod'
    tag "Converting BAM to FASTQ"
    cpus Math.min(params.threads as int, Runtime.runtime.availableProcessors())


    input:
        path(modbam)
    
    output:
        path("mod_fastq.fastq.gz"), emit: mod_fastq

    script:
    """
    samtools fastq -@ ${task.cpus} -T MM,ML ${modbam} | gzip > mod_fastq.fastq.gz
    """
}

process ALIGN_TO_REF {

    label 'telomod'
    tag "Aligning FASTQ to Reference"
    cpus Math.min(params.threads as int, Runtime.runtime.availableProcessors())

    input:
        path(mod_fastq)
        path(reference)

    output:
        path("alignment.bam"), emit: alignment

    script:
    """
    minimap2 -x map-ont -a -t ${params.threads} -y --secondary=no ${reference} ${mod_fastq} | samtools view -b -@ ${params.threads} > alignment.bam
    """
}


process GENOMIC_READS_WITH_SPIKE {

    label 'telomod'
    tag 'Genomic Modification Analysis w/ Spike-In'

    input:
        path(ref_aln), stageAs: "reference.aln.bam"
        path(spike_aln), stageAs: "spikeIn.aln.bam"
        path(telo_stats)
        path(modcalls)
    
    output:
        path("genomic_summary.txt"), emit: genomic_reads

    publishDir "${params.outdir}/", mode: 'copy', overwrite:true, pattern:"genomic_summary.txt"

    script:
    """
    ${baseDir}/bin/human/analyze_genomic_reads_with_spike_in.py --reference_aln ${ref_aln} \
                                            --spike_in_aln ${spike_aln} \
                                            --telo_stats ${telo_stats} \
                                            --reference_fa ${params.reference} \
                                            --minimum_read_length ${params.minimum_genomic_read_length} \
                                            --mod_calls ${modcalls} \
                                            --modified_nucleotide ${params.modification_nucl} \
                                            --modification_code ${params.modification_code} \
                                            --out_file genomic_summary.txt
    """
}

process GENOMIC_READ_ANALYSIS {

    label 'telomod'
    tag 'Genomic Modification Analysis'

    input:
        path(ref_aln)
        path(modcalls)

    output:
        path("genomic_summary.txt"), emit: genomic_reads

    publishDir "${params.outdir}/", mode: 'copy', overwrite:true, pattern:"genomic_summary.txt"

    script:
    """
    ${baseDir}/bin/human/analyze_genomic_reads_with_spike_in.py --reference_aln ${ref_aln} \
                                            --telo_stats ${telo_stats} \
                                            --reference_fa ${params.reference} \
                                            --minimum_read_length ${params.minimum_genomic_read_length} \
                                            --mod_calls ${modcalls} \
                                            --modified_nucleotide ${params.modification_nucl} \
                                            --modification_code ${params.modification_code} \
                                            --out_file genomic_summary.txt
    """
}

process TELO_READ_ANALYSIS {

    label 'telomod'
    tag 'Telomeric Modification Analysis'

    input:
        path(mod_bam)
        path(telo_stats)
        path(modcalls)
    
    output:
        path("telomeric_mod_location.pdf")
        path("telomeric_mod_summary.txt"), emit: telomeric_reads

    publishDir "${params.outdir}/", mode: 'copy', overwrite:true, pattern:"telomeric_mod_location.pdf"
    publishDir "${params.outdir}/", mode: 'copy', overwrite:true, pattern:"telomeric_mod_summary.txt"

    script:
    """
    python -u ${baseDir}/bin/human/analyze_telomeric_reads.py --reference_aln ${mod_bam} \
                                --telo_stats ${telo_stats} \
                                --mod_calls ${modcalls} \
                                --modified_nucleotide ${params.modification_nucl} \
                                --modification_code ${params.modification_code} \
                                --telomere_plot telomeric_mod_location.pdf \
                                --summary_file telomeric_mod_summary.txt \
                                --max_subtelo_stretch ${params.max_subtelo_stretch} \
                                --c_strand_only ${params.c_strand_only} \
                                --image_width ${params.image_width}\
                                --image_height ${params.image_height}   
    """
}

process EXTRACT_TELO_READS {

    label 'telomod'
    tag 'Extracting Telomeric Reads from ModBam'

    input:
        path(mod_bam)
        path(telo_stats)

    output:
        path("telomeric_modbam.sorted.bam"), emit: telo_modbam

    script:
    """
    ${baseDir}/bin/human/extract_telomeric_modbam.py --mod_bam ${mod_bam} --telo_stats ${telo_stats} --out_file telomeric_modbam.bam
    samtools sort -n telomeric_modbam.bam > telomeric_modbam.sorted.bam
    """
}

process REORIENT_TELOBAM {

    label 'telomod'
    tag 'Reorienting TARPON BAM in Strand specificity'

    input:
        path(telo_bam)
        path(stats_file)
        path(ref_bam)

    output:
        path("telomeric.reorientated.sorted.bam"), emit: telo_bam

    script:
    """
    ${baseDir}/bin/human/reorient_telobam.py --telo_bam ${telo_bam} --stats_file ${stats_file} --out_bam telomeric.reorientated.bam --ref_bam ${ref_bam}
    samtools sort -n telomeric.reorientated.bam > telomeric.reorientated.sorted.bam
    """
}
process REPAIR_TELOBAM {

    label 'modkit'
    cpus Math.min(params.threads as int, Runtime.runtime.availableProcessors())
    tag 'Repairing TeloBam Modifications'
    errorStrategy 'retry'
    maxRetries 3

    input:
        path(donor_file)
        path(acceptor_file)
    
    output:
        path("telo_modbam.bam"), emit: telo_modbam

    script:
    """
    modkit repair --donor-bam ${donor_file} --acceptor-bam ${acceptor_file} --output-bam telo_modbam.bam --log-filepath repair.log.txt
    """
}

// process TELO_CLUSTER_ANALYSIS {

//     label 'telomod'
//     tag 'Cluster Specific Telomere Analysis'

//     input:
//         path(cluster_results)
//         path(mod_bam)
//         path(telo_stats)
//         path(modcalls)

//     output:
//         path("clustering_results.txt")
//         path("cluster_assignment.txt"), emit: cluster_assignment

//     publishDir "${params.outdir}/", mode: 'copy', overwrite:true, pattern:"clustering_results.txt"

//     // plots must include cluster size distribution
//     script:
//     """
//     ${baseDir}/bin/human/analyze_clustering_results.py --cluster_file ${cluster_results} \
//                             --cluster_out_fh cluster_assignment.txt \
//                             --mod_bam ${mod_bam} \
//                             --telo_stats ${telo_stats}

//     mv .command.out clustering_results.txt
//     """
// }

process READ_BY_READ_ANALYSIS {

    label 'telomod'
    tag 'Read-by-Read Analysis'

    input:
        path(mod_bam)
        path(telo_stats)
        path(modcalls)
        path(telo_summary)

    output:
        path("*.pdf")

    publishDir "${params.outdir}/CLUSTERS/", mode: 'copy', overwrite:true, pattern:"cluster*.pdf"
    publishDir "${params.outdir}/CLUSTERS/", mode: 'copy', overwrite:true, pattern:"*.pdf"

    script:
    """
    
    export PYTHONUNBUFFERED=1
    ${baseDir}/bin/human/cluster_specific_modification_analysis.py --mod_bam ${mod_bam} \
                                                --telo_stats ${telo_stats} \
                                                --mod_table ${modcalls} \
                                                --image_width ${params.image_width} \
                                                --image_height ${params.image_height} \
                                                --modification ${params.modification_code} \
                                                --c_strand_only ${params.c_strand_only}

    ${baseDir}/bin/human/cluster_plots.R ${telo_summary} ${telo_stats} ${params.c_strand_only}
    """
}


process PLOT_RESULTS {

    label 'telomod'
    tag 'Comparing Genomic to Telomeric Reads'

    input:
        path(genomic_summary)
        path(telomeric_summary)
        path(telo_stats)
    
    output:
        path("*.pdf")
    
    publishDir "${params.outdir}/", mode: 'copy', overwrite:true, pattern:"*.pdf"

    script:
    """
    ${baseDir}/bin/human/global_plots.R ${genomic_summary} ${telomeric_summary} ${telo_stats} ${params.c_strand_only}
    """
}

process PLOT_TELO_RESULTS {

    label 'telomod'
    tag 'Comparing Genomic to Telomeric Reads'

    input:
        path(telomeric_summary)
        path(telo_stats)
    
    output:
        path("*.pdf")
    
    publishDir "${params.outdir}/", mode: 'copy', overwrite:true, pattern:"*.pdf"

    script:
    """
    ${baseDir}/bin/human/global_telo_plots.R ${telomeric_summary} ${telo_stats} ${params.c_strand_only}
    """
}

process CREATE_CLUSTER_FASTA {

    label 'telomod'
    tag 'Cluster Cluster FASTA'

    input:
        path(mod_bam)
        path(telo_stats)

    output:
        path("*.fa"), emit: fasta
        path("*.bam"), emit: bam
       // tuple path("*sorted.bam"), path("*.consensus.fa"), emit: cluster_files

    script:
    """
    ${baseDir}/bin/human/create_cluster_fasta_and_bam_new.py --telobam ${mod_bam} \
                                                --telo_stats ${telo_stats} \
                                                --max_read_length ${params.max_read_length}
    """
}

process GENERATE_CONSENSUS_AND_ALIGN {
    maxForks 5
    label 'telomod'
    tag 'Generate Consensus and Align Reads'
    maxRetries 5


    input:
        tuple val(cluster), path(cluster_fasta), path(cluster_bam)

    output:
        tuple val(cluster), path("*.consensus.fa"), path("*.c_strand.sorted.bam"), path("*.g_strand.sorted.bam"), emit: cluster_files

    publishDir "${params.outdir}/CLUSTERS/", mode: 'copy', overwrite:true, pattern:"*.consensus.fa"

    script:
    def extra = params.c_strand_only ?
        """
        touch ${cluster}.asm10.g_strand.sorted.bam
        """
    
    :
        """
        samtools view -F 2068 -b ${cluster}.asm10.sam > ${cluster}.asm10.g_strand.bam
        samtools sort -n ${cluster}.asm10.g_strand.bam > ${cluster}.asm10.g_strand.sorted.bam
        """
    """
    spoa -l 0 -r 0 ${cluster_fasta} > ${cluster}.consensus.fa

    samtools fastq ${cluster_bam} > ${cluster}.fastq
    minimap2 -ax asm10 -N 1 -f 0.0001 ${cluster}.consensus.fa ${cluster}.fastq > ${cluster}.asm10.sam

    samtools view -f 16 -F 2048 -b ${cluster}.asm10.sam > ${cluster}.asm10.c_strand.bam
    samtools sort -n ${cluster}.asm10.c_strand.bam > ${cluster}.asm10.c_strand.sorted.bam
    ${extra}
    """
}


process REPAIR {
    maxForks 5
    label 'modkit'
    tag 'MODKIT Repair and PileUp'
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val(cluster), path(consensus), path(c_strand_aln), path(g_strand_aln)
        path(ref_telobam)

    output:
        tuple val(cluster), path(consensus), path("*c_strand_mod.bam"), path("*g_strand_mod.bam"), emit: repaired

    script:

    def extra = params.c_strand_only ?
        """
        touch ${cluster}.g_strand_mod.bam
        """
    :
        """
        modkit repair --donor ${ref_telobam} --acceptor ${g_strand_aln} --output-bam ${cluster}.g_strand_mod.bam
        """
    
    """
    modkit repair --donor ${ref_telobam} --acceptor ${c_strand_aln} --output-bam ${cluster}.c_strand_mod.bam
    ${extra}
    """
}

process INDEX {

    maxForks 5
    label 'telomod'

    tag 'Sorting and Indexing New BAM'

   input:
        tuple val(cluster), path(consensus), path(c_strand_aln), path(g_strand_aln)

    output:
        tuple val(cluster), path(consensus), path("*.c_strand_aln.sorted.bam"), path("*.g_strand_aln.sorted.bam"), path("*.bai"), emit:indexed

    script:
    def extra = params.c_strand_only ?
        """
        touch ${cluster}.g_strand_aln.sorted.bam
        touch ${cluster}.g_strand_aln.sorted.bam.bai
        """
    :
        """
        samtools sort ${g_strand_aln} > ${cluster}.g_strand_aln.sorted.bam
        samtools index ${cluster}.g_strand_aln.sorted.bam
        """
    """
    samtools sort ${c_strand_aln} > ${cluster}.c_strand_aln.sorted.bam
    samtools index ${cluster}.c_strand_aln.sorted.bam
    ${extra}
    """
}

process PILEUP {

    maxForks 5
    label 'modkit'
    tag 'MODKIT PileUp'

    input:
        tuple val(cluster), path(consensus), path(c_strand_aln), path(g_strand_aln), path(indexes)

    output:
        tuple val(cluster), path(consensus), path(c_strand_aln), path(g_strand_aln), path("*c_strand.bed"), path("*g_strand.bed"), emit: results
        path("*.bed")
        
    publishDir "${params.outdir}/CLUSTERS/", mode: 'copy', overwrite:true, pattern:"*.bed"

    script:
    def extra = params.c_strand_only ?
        """
        touch ${cluster}.g_strand.bed
        """
    :
        """
        modkit pileup --filter-threshold ${params.mod_confidence} ${g_strand_aln} ${cluster}.g_strand.bed
        """
    """
    modkit pileup --filter-threshold ${params.mod_confidence} ${c_strand_aln} ${cluster}.c_strand.bed
    ${extra}
    """
}


process PLOT_PILEUP_RESULTS {

    maxForks 8
    label 'telomod'
    tag 'Plotting PileUp Results'

    input:  
        tuple val(cluster), path(consensus), path(c_strand_aln), path(g_strand_aln), path(c_strand_bed), path(g_strand_bed)

    output:
        path("*.pdf")

    publishDir "${params.outdir}/CLUSTERS/", mode: 'copy', overwrite:true, pattern:"*.pdf"

    script:
    def extra = params.c_strand_only ?
        """
        touch g_strand.cov.txt
        """
    :
        """
        samtools depth ${g_strand_aln} > g_strand.cov.txt
        """
    """
    samtools depth ${c_strand_aln} > c_strand.cov.txt
    ${extra}   
    echo "depth calculated, starting plot"

    ${baseDir}/bin/human/plot_pileup_results.py --consensus ${consensus} --c_strand_pileup ${c_strand_bed} \
                            --g_strand_pileup ${g_strand_bed} \
                            --c_strand_cov c_strand.cov.txt \
                            --g_strand_cov g_strand.cov.txt \
                            --out_file ${cluster}.pileup.pdf \
                            --modification ${params.modification_code} \
                            --modified_nucleotide ${params.modification_nucl} \
                            --image_width ${params.image_width} --image_height ${params.image_height} \
                            --c_strand_only ${params.c_strand_only}
    """
}