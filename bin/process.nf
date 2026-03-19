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
            --filter-threshold 0.9 ${modbam} all_modcalls.txt
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
            --filter-threshold 0.9 ${modbam} telo_modcalls.txt.gz
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
    analyze_genomic_reads_with_spike_in.py --reference_aln ${ref_aln} \
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
    analyze_genomic_reads_with_spike_in.py --reference_aln ${ref_aln} \
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
    analyze_telomeric_reads.py --reference_aln ${mod_bam} \
                                --telo_stats ${telo_stats} \
                                --mod_calls ${modcalls} \
                                --modified_nucleotide ${params.modification_nucl} \
                                --modification_code ${params.modification_code} \
                                --telomere_plot telomeric_mod_location.pdf \
                                --summary_file telomeric_mod_summary.txt \
                                --max_subtelo_stretch ${params.max_subtelo_stretch} \
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
        path("telomeric_modbam.bam"), emit: telo_modbam

    script:
    """
    extract_telomeric_modbam.py --mod_bam ${mod_bam} --telo_stats ${telo_stats} --out_file telomeric_modbam.bam
    """
}

process TELO_CLUSTER_ANALYSIS {

    label 'telomod'
    tag 'Cluster Specific Telomere Analysis'

    input:
        path(cluster_results)
        path(mod_bam)
        path(telo_stats)
        path(modcalls)
        path(telo_summary)

    output:
        path("clustering_results.txt")
        path("cluster_assignment.txt"), emit: cluster_assignment
        path("*.pdf")

    publishDir "${params.outdir}/", mode: 'copy', overwrite:true, pattern:"clustering_results.txt"
    publishDir "${params.outdir}/CLUSTERS/", mode: 'copy', overwrite:true, pattern:"cluster*.pdf"
    publishDir "${params.outdir}/CLUSTERS/", mode: 'copy', overwrite:true, pattern:"*.pdf"

    // plots must include cluster size distribution
    script:
    """
    analyze_clustering_results.py --cluster_file ${params.cluster_results} \
                            --cluster_out_fh cluster_assignment.txt \
                            --mod_bam ${mod_bam} \
                            --telo_stats ${telo_stats}

    mv .command.out clustering_results.txt

    cluster_specific_modification_analysis.py --cluster_file cluster_assignment.txt \
                                                --mod_bam ${mod_bam} \
                                                --telo_stats ${telo_stats} \
                                                --mod_table ${modcalls} \
                                                --image_width ${params.image_width} \
                                                --image_height ${params.image_height} \
                                                --modification ${params.modification_code} 

    cluster_plots.R ${telo_summary} ${telo_stats} cluster_assignment.txt
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
    global_plots.R ${genomic_summary} ${telomeric_summary} ${telo_stats}
    """
}