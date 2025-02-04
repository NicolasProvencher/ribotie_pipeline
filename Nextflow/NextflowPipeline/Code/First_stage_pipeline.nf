#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.input = "/home/machinegun/output/**/SRR*/*.fastq"
params.outdir = "pipeline_output"

def input_files = Channel
    .fromPath(params.input)
    .map { file ->
        def gse = file.parent.parent.name
        def srr = file.parent.name
        def srr_dir = file.parent
        def files = srr_dir.listFiles()
        def is_paired = files.size() == 2
        if (is_paired && !file.name.endsWith("_1.fastq")) {
            return null
        }
        tuple(gse, srr, is_paired, files)
    }
    .filter { it != null }
    .unique()
    .branch {
        view_tuples: true
            return it
        count_tuples: true
            return it
    }

input_files.view_tuples
           .filter {it[2]== true}
           .count()
           .view{"Paired-end files count: ${it}"}

input_files.view_tuples
          .filter{it[2]==false}
          .count()
          .view{"Single-end files count: ${it}"}

tuples_paires = input_files.view_tuples
    .filter { it[2] == true }
    .map { gse, srr, is_paired, files ->
        def read1 = file(files.find { it.name.contains("_1.fastq") }.toString())
        def read2 = file(files.find { it.name.contains("_2.fastq") }.toString())
        tuple(gse, srr, is_paired, [read1, read2])
    }

tuples_singles = input_files.view_tuples
    .filter { it[2] == false }
    .map { gse, srr, is_paired, files ->
        tuple(gse, srr, is_paired, files[0])
    }

process FASTQC_PRE_SINGLE {
    tag "$sample_id"
    publishDir "${params.outdir}/${gse}/${sample_id}/fastqc_pre", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(gse), val(sample_id), val(is_paired), path(reads)

    output:
    tuple val(gse), val(sample_id), val(is_paired), path("*.html"), path("*.zip"), emit: fastqc_pre_output
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc --quiet $reads
    """
}

process FASTQC_PRE_PAIR {
    tag "$sample_id"
    publishDir "${params.outdir}/${gse}/${sample_id}/fastqc_pre", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(gse), val(sample_id), val(is_paired), path(reads)

    output:
    tuple val(gse), val(sample_id), val(is_paired), path("*.html"), path("*.zip"), emit: fastqc_pre_output
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc --quiet ${reads[0]}
    fastqc --quiet ${reads[1]}
    """
}

process TRIM_GALORE_SINGLE {
    tag "$srr"
    publishDir "${params.outdir}/${gse}/${srr}/trimmed", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(gse), val(srr), val(is_paired), path(reads)

    output:
    tuple val(gse), val(srr), val(is_paired), path("*_trimmed.fq"), emit: trimmed
    path "*_trimming_report.txt"

    script:
    """
    trim_galore \
        --fastqc \
        --trim-n \
        --length 25 \
        --quality 20 \
        --max_length 40 \
        ${reads}
    """
}

process TRIM_GALORE_PAIRED {
    tag "$srr"
    publishDir "${params.outdir}/${gse}/${srr}/trimmed", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(gse), val(srr), val(is_paired), path(reads)

    output:
    tuple val(gse), val(srr), val(is_paired), path("*_val_{1,2}.fq"), emit: trimmed
    path "*_trimming_report.txt"

    script:
    """
    trim_galore \
        --paired \
        --fastqc \
        --trim-n \
        --length 25 \
        --quality 20 \
        ${reads[0]} ${reads[1]}
    """
}

process FASTQC_POST_SINGLE {
    tag "$sample_id"
    publishDir "${params.outdir}/${gse}/${sample_id}/fastqc_post", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(gse), val(sample_id), val(is_paired), path(trimmed_reads)

    output:
    tuple val(gse), val(sample_id), val(is_paired), path("*.html"), path("*.zip"), emit: fastqc_post_output
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc --quiet $trimmed_reads
    """
}

process FASTQC_POST_PAIR {
    tag "$sample_id"
    publishDir "${params.outdir}/${gse}/${sample_id}/fastqc_post", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(gse), val(sample_id), val(is_paired), path(trimmed_reads)

    output:
    tuple val(gse), val(sample_id), val(is_paired), path("*.html"), path("*.zip"), emit: fastqc_post_output
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc --quiet ${trimmed_reads[0]}
    fastqc --quiet ${trimmed_reads[1]}
    """
}

workflow {
    FASTQC_PRE_SINGLE(tuples_singles)
    FASTQC_PRE_PAIR(tuples_paires)
    
    TRIM_GALORE_SINGLE(tuples_singles)
    TRIM_GALORE_PAIRED(tuples_paires)
    
    FASTQC_POST_SINGLE(TRIM_GALORE_SINGLE.out.trimmed)
    FASTQC_POST_PAIR(TRIM_GALORE_PAIRED.out.trimmed)
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
}
