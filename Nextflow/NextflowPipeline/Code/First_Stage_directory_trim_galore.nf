#!/usr/bin/env nextflow
// Le code marche bien il me reste de verifier si le MultiQc fonctionne comme c'est souhaiter 
nextflow.enable.dsl = 2

params.input = "/home/machinegun/output2/**/*.fastq"
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
    memory "1GB"

    input:
    tuple val(gse), val(sample_id), path(reads)

    output:
    tuple val(gse), val(sample_id), path("*_fastqc.html"), path("*_fastqc.zip"), emit: fastqc_output
    
    script:
    """
    fastqc --quiet $reads
    """
}

process TRIM_GALORE {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${gse}/${sample_id}/trim", mode: 'copy'
    memory "1GB"

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
    """
}

process FASTQC_POST {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${gse}/${sample_id}/fastqc_post", mode: 'copy'
    memory "1GB"

    input:
    tuple val(gse), val(sample_id), path(trimmed_reads)

    output:
    tuple val(gse), val(sample_id), path("*_fastqc.html"), path("*_fastqc.zip"), emit: fastqc_output
    
    script:
    """
    fastqc --quiet $trimmed_reads
    """
}

process MULTIQC {
    tag "$gse"
    publishDir "${params.outdir}/${gse}", mode: 'copy'
    memory "1GB"
    input:
    tuple val(gse), path(fastqc_pre), path(fastqc_post), path(trim_reports)

    output:
    path "multiqc/*"
    
    script:
    """
    mkdir -p fastqc_pre fastqc_post trim_reports
    cp -L $fastqc_pre fastqc_pre/
    cp -L $fastqc_post fastqc_post/
    cp -L $trim_reports trim_reports/
    
    mkdir -p multiqc
    multiqc . -o multiqc
    """
}
    // Combinaison des résultats pour MultiQC
    // fastqc_pre_grouped
    //     .join(fastqc_post_grouped)
    //     .join(trim_reports_grouped)
    //     .set { multiqc_input }
workflow {
    FASTQC_PRE(input_files)
    TRIM_GALORE(input_files)
    FASTQC_POST(TRIM_GALORE.out.trimmed_reads)

    // Regroupement des résultats par GSE
    FASTQC_PRE.out.fastqc_output
        .map { gse, sample_id, html, zip -> tuple(gse, html, zip) }
        .groupTuple(by: 0)
        .map { gse, html, zip -> tuple(gse, html.flatten() + zip.flatten()) }
        .set { fastqc_pre_grouped }

    FASTQC_POST.out.fastqc_output
        .map { gse, sample_id, html, zip -> tuple(gse, html, zip) }
        .groupTuple(by: 0)
        .map { gse, html, zip -> tuple(gse, html.flatten() + zip.flatten()) }
        .set { fastqc_post_grouped }

    TRIM_GALORE.out.report
        .map { gse, sample_id, report -> tuple(gse, report) }
        .groupTuple(by: 0)
        .set { trim_reports_grouped }

    // Combinaison des résultats pour MultiQC
    // fastqc_pre_grouped
    //     .join(fastqc_post_grouped)
    //     .join(trim_reports_grouped)
    //     .set { multiqc_input }

    // MULTIQC(multiqc_input)
}
