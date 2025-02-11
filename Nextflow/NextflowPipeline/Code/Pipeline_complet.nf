#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.max_retries = 3
params.outdir = '/project/def-xroucou/riboseq_pipeline/output_complet'

csvChannelPerLigne = Channel
    .fromPath('/home/ilyass09/Sample_Sheet_test.csv')
    .splitCsv(header: true)
    .map { row -> tuple(
        row.Study_accession,
        row.Sample_accession,
        row.Drug,
        row.Biological_type,
        row.Trim_arg,
        row.S_P_type
    )}
    

csvChannelPerGsm = csvChannelPerLigne
    .map { gse, gsm, drug, bio, trim, sp -> tuple(
        gse,
        gsm.split(';').collect{it.trim()},
        drug,
        bio,
        trim,
        sp
    )}
    .transpose(by: 1) 
    

process FASTQ_DUMP {
    // errorStrategy 'retry'
    maxRetries params.max_retries
    publishDir "${params.outdir}/${gse}_${drug}_${bio}/${gsm}", mode: 'copy'

    input:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp)
    
    output:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp), path("*.fastq"), emit: fastq_files

    script:
    """
    echo "[INFO] FASTQ_DUMP : Starting download of GSM: ${gsm}"
    module load sra-toolkit
    fasterq-dump ${gsm}
    echo "[SUCCESS] FASTQ_DUMP : Downloaded GSM: ${gsm}"
    """
    
}

process FASTQC_PRE {
    tag "$gsm"
    publishDir "${params.outdir}/${gse}_${drug}_${bio}/${gsm}/fastqc_pre", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp), path(reads)

    output:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp), path("*.html"), path("*.zip"), emit: fastqc_output
    
    script:
    if (sp.toLowerCase() == "paired") {
        """
        echo "[INFO] FASTQC_PRE : Starting download of GSM: ${gsm}"
        fastqc --quiet ${reads[0]}
        fastqc --quiet ${reads[1]}
        echo "[SUCCESS] FASTQC_PRE : Downloaded GSM: ${gsm}"

        """
    } else {
        """
        fastqc --quiet $reads
        """
    }
}

process TRIM_GALORE {
    tag "$gsm"
    publishDir "${params.outdir}/${gse}/${gsm}/trimmed", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp), path(reads)

    output:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp), path("*_trimmed.fq"), emit: trimmed
    path "*_trimming_report.txt"

    //s'il y a une valeur pour le trim on l'ajoute sinon le trim va etre juste une valeur vide " " donc pas e trim special .
    script:
    if (sp.toLowerCase() == "paired") {
       
        """
        echo "[INFO] TRIM_GALORE : Starting download of GSM: ${gsm}"
        trim_galore \
        --paired \
        --fastqc \
        --trim-n \
        --length 20 \
        --quality 25 \
        ${trim} \
        ${reads[0]} ${reads[1]}
        echo "[SUCCESS] TRIM_GALORE : Downloaded GSM: ${gsm}"        
        """
    } else {
        """
        echo "[INFO] TRIM_GALORE : Starting download of GSM: ${gsm}"        
         trim_galore \
        --fastqc \
        --trim-n \
        --length 20 \
        --quality 25 \
        --max_length 40 \
        ${trim} \
        ${reads}
        echo "[SUCCESS] TRIM_GALORE : Downloaded GSM: ${gsm}"           
        """
    }
}

process FASTQC_POST {
    tag "$gsm"
    publishDir "${params.outdir}/${gse}_${drug}_${bio}/${gsm}/fastqc_pre", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp), path(trimmed_reads)

    output:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp), path("*.html"), path("*.zip"), emit: fastqc_output
    
    script:
    if (sp.toLowerCase() == "paired") {
        """
        echo "[INFO] FASTQC_POST : Starting download of GSM: ${gsm}"     
        fastqc --quiet ${trimmed_reads[0]}
        fastqc --quiet ${trimmed_reads[1]}
        echo "[SUCCESS] FASTQC_POST : Downloaded GSM: ${gsm}"           

        """
    } else {
        """
        echo "[INFO] FASTQC_POST : Starting download of GSM: ${gsm}"   
        fastqc --quiet $trimmed_reads
        echo "[SUCCESS] FASTQC_POST : Downloaded GSM: ${gsm}"     
        """
    }
}

log.info "Starting pipeline..."

workflow {
    // Download des fichiers fastq
    FASTQ_DUMP(csvChannelPerGsm)
    // Monitoring des GSM téléchargés
    FASTQ_DUMP.out.fastq_files
        .map { gse, gsm, drug, bio, trim, sp, reads ->
            "[DOWNLOAD COMPLETED] GSM: $gsm (Study: $gse, Drug: $drug)"
        }
        .view()
    // Application du fastqc rapport
    FASTQC_PRE(FASTQ_DUMP.out.fastq_files)
    // Application du trim
    TRIM_GALORE(FASTQ_DUMP.out.fastq_files)
    // Application du fastqc rapport trimmed
    FASTQC_POST(TRIM_GALORE.out.trimmed)

}
workflow.onComplete {
    log.info """
    Pipeline execution summary
    -------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    """
}
