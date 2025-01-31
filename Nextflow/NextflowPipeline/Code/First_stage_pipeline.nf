#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input = "/home/ilyass09/output_test/**/*.fastq" 
params.outdir = "pipeline_output_test"

log.info """\
    PIPELINE RIBO-SEQ - QC ET TRIMMING
    ==================================
    Input    : ${params.input}
    Outdir   : ${params.outdir}
    """
    .stripIndent()

Channel
    .fromPath(params.input)
    .map { file -> 
        def gse = file.parent.name
        tuple(gse, file.simpleName, file)
    }
    .set { input_files }

process FASTQC_PRE {
    tag "$sample_id"
    
    
    publishDir "${params.outdir}/${gse}/${sample_id}/fastqc_pre", mode: 'copy'
    

    input:
    tuple val(gse), val(sample_id), path(reads)

    output:
    tuple val(gse), val(sample_id), path("*_fastqc.html"), path("*_fastqc.zip"), emit: fastqc_output
    
    script:
    """/home/ilyass09/env.yml
    fastqc --quiet $reads
    """
}

process TRIM_GALORE {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${gse}/${sample_id}/trim", mode: 'copy'
    

    input:
    tuple val(gse), val(sample_id), path(reads)

    output:
    tuple val(gse), val(sample_id), path("*_trimmed.fq"), emit: trimmed_reads
    tuple val(gse), val(sample_id), path("*_trimming_report.txt"), emit: report
    
    script:
    """
    trim_galore \
        --fastqc \
        --trim-n \
        --length 25 \
        --quality 20 \
        --max_length 40 \
        $reads
    """/home/ilyass09/env.yml
}

process FASTQC_POST {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${gse}/${sample_id}/fastqc_post", mode: 'copy'
    

    input:
    tuple val(gse), val(sample_id), path(trimmed_reads)

    output:
    tuple val(gse), val(sample_id), path("*_fastqc.html"), path("*_fastqc.zip"), emit: fastqc_output
    
    script:
    """
    fastqc --quiet $trimmed_reads
    """
}

workflow {
    FASTQC_PRE(input_files)
    TRIM_GALORE(input_files)
    FASTQC_POST(TRIM_GALORE.out.trimmed_reads)
    

// 
}

workflow.onComplete {
    log.info """
    Pipeline terminé avec succès!
    maintenant execution de Multiqc 
    """

    
}

