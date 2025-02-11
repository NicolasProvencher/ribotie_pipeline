#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.max_retries = 3
params.outdir = '/project/def-xroucou/riboseq_pipeline/output_complet'

// Fonctions de vérification pour chaque type de fichier
def checkFastqExists(gse, gsm, drug, bio, sp) {
    def outputPath = "${params.outdir}/${gse}_${drug}_${bio}/${gsm}"
    if (sp.toLowerCase() == "paired") {
        def file1 = new File("${outputPath}/${gsm}_1.fastq")
        def file2 = new File("${outputPath}/${gsm}_2.fastq")
        return file1.exists() && file2.exists()
    } else {
        def file = new File("${outputPath}/${gsm}.fastq")
        return file.exists()
    }
}

def checkFastqcPreExists(gse, gsm, drug, bio, sp) {
    def outputPath = "${params.outdir}/${gse}_${drug}_${bio}/${gsm}/fastqc_pre"
    if (sp.toLowerCase() == "paired") {
        return new File("${outputPath}/${gsm}_1_fastqc.html").exists() &&
               new File("${outputPath}/${gsm}_2_fastqc.html").exists() &&
               new File("${outputPath}/${gsm}_1_fastqc.zip").exists() &&
               new File("${outputPath}/${gsm}_2_fastqc.zip").exists()
    } else {
        return new File("${outputPath}/${gsm}_fastqc.html").exists() &&
               new File("${outputPath}/${gsm}_fastqc.zip").exists()
    }
}

def checkTrimExists(gse, gsm, drug, bio, sp) {
    def outputPath = "${params.outdir}/${gse}/${gsm}/trimmed"
    if (sp.toLowerCase() == "paired") {
        return new File("${outputPath}/${gsm}_1_val_1.fq").exists() &&
               new File("${outputPath}/${gsm}_2_val_2.fq").exists()
    } else {
        return new File("${outputPath}/${gsm}_trimmed.fq").exists()
    }
}

def checkFastqcPostExists(gse, gsm, drug, bio, sp) {
    def outputPath = "${params.outdir}/${gse}_${drug}_${bio}/${gsm}/fastqc_pre"
    if (sp.toLowerCase() == "paired") {
        return new File("${outputPath}/${gsm}_1_val_1_fastqc.html").exists() &&
               new File("${outputPath}/${gsm}_2_val_2_fastqc.html").exists() &&
               new File("${outputPath}/${gsm}_1_val_1_fastqc.zip").exists() &&
               new File("${outputPath}/${gsm}_2_val_2_fastqc.zip").exists()
    } else {
        return new File("${outputPath}/${gsm}_trimmed_fastqc.html").exists() &&
               new File("${outputPath}/${gsm}_trimmed_fastqc.zip").exists()
    }
}

// Définition des channels avec vérification

def create_initial_channels() {
    def csvChannelPerLigne = Channel
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

    return csvChannelPerLigne
        .map { gse, gsm, drug, bio, trim, sp -> tuple(
            gse,
            gsm.split(';').collect{it.trim()},
            drug,
            bio,
            trim,
            sp
        )}
        .transpose(by: 1)
}


process FASTQ_DUMP {
    maxRetries params.max_retries
    maxForks 10  // Limit parallel execution to 5 concurrent jobs

    publishDir "${params.outdir}/${gse}_${drug}_${bio}/${gsm}", mode: 'copy'

    input:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp)
    
    output:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp), path("*.fastq"), emit: fastq_files

    when:
    !checkFastqExists(gse, gsm, drug, bio, sp)
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
    // errorStrategy 'ignore'

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

// Process TRIM_GALORE_SINGLE pour les fichiers single-end
process TRIM_GALORE_SINGLE {
    tag "$gsm"
    publishDir "${params.outdir}/${gse}_${drug}_${bio}/${gsm}/trimmed", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp), path(reads)

    output:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp), path("*_trimmed.fq"), emit: trimmed
    path "*_trimming_report.txt"



    script:
    """
    echo "[INFO] TRIM_GALORE_SINGLE : Starting trim of GSM: ${gsm}"
    trim_galore \
    --fastqc \
    --trim-n \
    --length 20 \
    --quality 25 \
    --max_length 40 \
    ${trim} \
    ${reads}
    echo "[SUCCESS] TRIM_GALORE_SINGLE : Completed GSM: ${gsm}"
    """
}

// Process TRIM_GALORE_PAIRED pour les fichiers paired-end
process TRIM_GALORE_PAIRED {
    tag "$gsm"
    publishDir "${params.outdir}/${gse}_${drug}_${bio}/${gsm}/trimmed", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp), path(reads)

    output:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp), path("*_val_*.fq"), emit: trimmed
    path "*_trimming_report.txt"

    script:
    """
    echo "[INFO] TRIM_GALORE_PAIRED : Starting trim of GSM: ${gsm}"
    trim_galore \
    --paired \
    --fastqc \
    --trim-n \
    --length 20 \
    --quality 25 \
    ${trim} \
    ${reads[0]} ${reads[1]}
    echo "[SUCCESS] TRIM_GALORE_PAIRED : Completed GSM: ${gsm}"
    """
}

process FASTQC_POST {
    tag "$gsm"
    publishDir "${params.outdir}/${gse}_${drug}_${bio}/${gsm}/fastqc_post", mode: 'copy'
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
        fastqc --quiet ${trimmed_reads}
        echo "[SUCCESS] FASTQC_POST : Downloaded GSM: ${gsm}"     
        """
    }
}

log.info "Starting pipeline..."


workflow {
    // Créer le channel initial
    def csvChannelPerGsm = create_initial_channels()
    
    // Download des fichiers fastq
    FASTQ_DUMP(csvChannelPerGsm)
    
    // Création d'un channel pour les fichiers existants
    def existing_fastq_files = csvChannelPerGsm
        .filter { gse, gsm, drug, bio, trim, sp ->
            checkFastqExists(gse, gsm, drug, bio, sp)
        }
        .map { gse, gsm, drug, bio, trim, sp ->
            def outputPath = "${params.outdir}/${gse}_${drug}_${bio}/${gsm}"
            def fastqFiles = sp.toLowerCase() == "paired" ?
                [file("${outputPath}/${gsm}_1.fastq"), file("${outputPath}/${gsm}_2.fastq")] :
                file("${outputPath}/${gsm}.fastq")
            tuple(gse, gsm, drug, bio, trim, sp, fastqFiles)
        }
    

    def combined_fastq_files = FASTQ_DUMP.out.fastq_files.mix(existing_fastq_files)
    
    combined_fastq_files
        .map { gse, gsm, drug, bio, trim, sp, reads ->
            "[FILE READY] GSM: $gsm (Study: $gse, Drug: $drug)"
        }
        .view()
        
    FASTQC_PRE(combined_fastq_files)
    
    def (singleFiles, pairedFiles) = combined_fastq_files
        .branch {
            single: it[5].toLowerCase() == "single"
            paired: it[5].toLowerCase() == "paired"
        }

    TRIM_GALORE_SINGLE(singleFiles)
    TRIM_GALORE_PAIRED(pairedFiles)
    
    def trimmed_files = TRIM_GALORE_SINGLE.out.trimmed.mix(TRIM_GALORE_PAIRED.out.trimmed)
    
    FASTQC_POST(trimmed_files)
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
