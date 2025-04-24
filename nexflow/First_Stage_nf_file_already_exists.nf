#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
// In this part
// everything works well but I need to verify the MultiQC to ensure everything is good like the trim
// This code is customized for Nicolas' files which are in disorder

// Define parameters
params.outdir_file_exist_already = "/path/to/directory/containing/fastq"
params.outdir = "/path/to/trim/output/directory"
params.input_csv = "/path/to/samplesheet/directory/samplesheet.csv"
params.max_retries = 3

// Function to check if files exist in the fastq directory
def checkFastqExists(gsm, sp) {
    def outputPath = "${params.outdir_file_exist_already}"
    def found = false
    
    try {
        def dir = new File(outputPath)
        def files = dir.list()
        
        if (sp.toLowerCase() == "paired") {
            def file1Exists = files.any { it == "${gsm}_1.fastq" }
            def file2Exists = files.any { it == "${gsm}_2.fastq" }
            found = file1Exists && file2Exists
        } else {
            found = files.any { it == "${gsm}.fastq" }
        }
    } catch (Exception e) {
        log.error "Error while checking files for GSM: $gsm - ${e.message}"
    }
    
    return found
}

// Verification functions for post-download steps
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
    def outputPath = "${params.outdir}/${gse}_${drug}_${bio}/${gsm}/trimmed"
    if (sp.toLowerCase() == "paired") {
        return new File("${outputPath}/${gsm}_1_val_1.fq").exists() &&
               new File("${outputPath}/${gsm}_2_val_2.fq").exists()
    } else {
        return new File("${outputPath}/${gsm}_trimmed.fq").exists()
    }
}

def checkFastqcPostExists(gse, gsm, drug, bio, sp) {
    def outputPath = "${params.outdir}/${gse}_${drug}_${bio}/${gsm}/fastqc_post"
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

// Definition of the create_initial_channels function
def create_initial_channels() {
    def csvChannelPerLigne = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> tuple(
            row.Study_accession,
            row.Sample_accession,
            row.Drug,
            row.Biological_type,
            row.Trim_arg,
            row.S_P_type,
            row.Species
        )}

    return csvChannelPerLigne
        .map { gse, gsm, drug, bio, trim, sp, type_cells -> tuple(
            gse,
            gsm.split(';').collect{it.trim()},
            drug,
            bio,
            trim,
            sp,
            type_cells
        )}
        .transpose(by: 1)
}

// Process for FASTQC_PRE
process FASTQC_PRE {
    tag "$gsm"
    publishDir "${params.outdir}/${type_cells}/${gse}_${drug}_${bio}/${gsm}/fastqc_pre", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp),val(type_cells), path(reads)

    output:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp),val(type_cells), path("*.html"), path("*.zip"), emit: fastqc_output

    when:
    !checkFastqcPreExists(gse, gsm, drug, bio, sp)

    script:
    if (sp.toLowerCase() == "paired") {
        """
        echo "[INFO] FASTQC_PRE : Starting fastqc of GSM: ${gsm}"
        fastqc --quiet ${reads[0]}
        fastqc --quiet ${reads[1]}
        echo "[SUCCESS] FASTQC_PRE : fastqc GSM: ${gsm}"
        """
    } else {
        """
        echo "[INFO] FASTQC_PRE : Starting fastqc of GSM: ${gsm}"
        fastqc --quiet $reads
        echo "[SUCCESS] FASTQC_PRE : fastqc GSM: ${gsm}"
        """
    }
}

// Process TRIM_GALORE_SINGLE for single-end files
process TRIM_GALORE_SINGLE {
    tag "$gsm"
    publishDir "${params.outdir}/${type_cells}/${gse}_${drug}_${bio}/${gsm}/trimmed", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp),val(type_cells), path(reads)

    output:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp),val(type_cells), path("*_trimmed.fq"), emit: trimmed
    path "*_trimming_report.txt"

    when:
    !checkTrimExists(gse, gsm, drug, bio, sp)

    script:
    """
    echo "[INFO] TRIM_GALORE_SINGLE : Starting trim of GSM: ${gsm}"
    trim_galore \
    --fastqc \
    --trim-n \
    --length 18 \
    --quality 25 \
    --max_length 40 \
    ${trim} \
    ${reads}
    echo "[SUCCESS] TRIM_GALORE_SINGLE : Completed GSM: ${gsm}"
    """
}

// Process TRIM_GALORE_PAIRED for paired-end files
process TRIM_GALORE_PAIRED {
    tag "$gsm"
    publishDir "${params.outdir}/${type_cells}/${gse}_${drug}_${bio}/${gsm}/trimmed", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp),val(type_cells), path(reads)

    output:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp),val(type_cells), path("*_val_*.fq"), emit: trimmed
    path "*_trimming_report.txt"

    when:
    !checkTrimExists(gse, gsm, drug, bio, sp)

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
    publishDir "${params.outdir}/${type_cells}/${gse}_${drug}_${bio}/${gsm}/fastqc_post", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp),val(type_cells), path(trimmed_reads)

    output:
    tuple val(gse), val(gsm), val(drug), val(bio), val(trim), val(sp),val(type_cells), path("*.html"), path("*.zip"), emit: fastqc_output
    
    when:
    !checkFastqcPostExists(gse, gsm, drug, bio, sp)

    script:
    if (sp.toLowerCase() == "paired") {
        """
        echo "[INFO] FASTQC_POST : Starting fastqc of GSM: ${gsm}"     
        fastqc --quiet ${trimmed_reads[0]}
        fastqc --quiet ${trimmed_reads[1]}
        echo "[SUCCESS] FASTQC_POST : Downloaded GSM: ${gsm}"           
        """
    } else {
        """
        echo "[INFO] FASTQC_POST : Starting fastqc of GSM: ${gsm}"   
        fastqc --quiet ${trimmed_reads}
        echo "[SUCCESS] FASTQC_POST : Downloaded GSM: ${gsm}"     
        """
    }
}



// Main workflow
workflow {
    log.info "Starting pipeline for existing files only..."
    // Create initial channel
    def csvChannelPerGsm = create_initial_channels()
    
    // Filter only for existing files
    def existing_channel = csvChannelPerGsm.filter {
        gse, gsm, drug, bio, trim, sp, type_cells -> checkFastqExists(gsm, sp)
    }
    
    // // Log existing files
    // existing_channel.view { gse, gsm, drug, bio, trim, sp ->
    //     "Existing file: GSE=$gse, GSM=$gsm - will be processed"
    // }
    
    // Map existing files to include file paths
    def existing_fastq_files = existing_channel
        .map { gse, gsm, drug, bio, trim, sp, type_cells ->
            def fastqFiles
            if (sp.toLowerCase() == "paired") {
                fastqFiles = [
                    file("${params.outdir_file_exist_already}/${gsm}_1.fastq"),
                    file("${params.outdir_file_exist_already}/${gsm}_2.fastq")
                ]
            } else {
                fastqFiles = file("${params.outdir_file_exist_already}/${gsm}.fastq")
            }
            tuple(gse, gsm, drug, bio, trim, sp,type_cells, fastqFiles)
        }
    
    // Execute FASTQC_PRE for existing files
    FASTQC_PRE(existing_fastq_files)
    
    // Separate single and paired files
    def (singleFiles, pairedFiles) = existing_fastq_files
        .branch {
            single: it[5].toLowerCase() == "single"
            paired: it[5].toLowerCase() == "paired"
        }

    // Execute trimming
    TRIM_GALORE_SINGLE(singleFiles)
    TRIM_GALORE_PAIRED(pairedFiles)
    
    // Combine trimming results
    def trimmed_files = TRIM_GALORE_SINGLE.out.trimmed.mix(TRIM_GALORE_PAIRED.out.trimmed)
    
    // Execute FASTQC_POST on trimmed files
    FASTQC_POST(trimmed_files)
}

// workflow.onComplete {
//     log.info """
//     Pipeline execution summary
//     -------------------------
//     Completed at: ${workflow.complete}
//     Duration    : ${workflow.duration}
//     Success     : ${workflow.success}
//     workDir     : ${workflow.workDir}
//     exit status : ${workflow.exitStatus}
//     """
// }
