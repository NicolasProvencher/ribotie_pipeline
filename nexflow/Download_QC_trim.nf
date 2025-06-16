#!/usr/bin/env nextflow
// In this part
// everything works well but I need to verify the MultiQC to ensure everything is good like the trim

// Define parameters
params.fastq_dir = "/path/to/directory/containing/new/fastq"
params.pipeline_dir = "/path/to/trim/output/directory"
params.input_csv = "/path/to/samplesheet/directory/samplesheet.csv"
params.max_retries = 3

process FASTQ_DUMP {
    maxRetries params.max_retries
    // maxForks 1  // Limit parallel execution to 5 concurrent jobs
    tag "$meta.GSM"
    publishDir "${params.fastq_dir}", mode: 'copy'

    input:
    val (meta)
    
    output:
    val (meta), emit: meta
    path("${meta.GSM}*.fastq"), emit: fastq_files

    script:
    """
    module load sra-toolkit
    echo "[INFO] FASTQ_DUMP : Starting download of GSM: ${meta.GSM}"
    fasterq-dump ${meta.GSM}
    echo "[SUCCESS] FASTQ_DUMP : Downloaded GSM: ${meta.GSM}"
    """
}

// Process for FASTQC_PRE
process FASTQC_PRE {
    tag "$meta.GSM"
    publishDir "${params.pipeline_dir}/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/fastqc_pre", mode: 'copy'
    errorStrategy 'retry'
    maxRetries params.max_retries

    input:
    val (meta)

    output:
    val (meta), emit : meta
    tuple path("*.html") , path("*.zip"), emit : files

    script:
    def reads = meta.sp ? "${params.fastq_dir}/${meta.GSM}_1.fastq ${params.fastq_dir}/${meta.GSM}_2.fastq" : "${params.fastq_dir}/${meta.GSM}.fastq"
        """
        echo "[INFO] FASTQC_PRE : Starting fastqc of GSM: ${meta.GSM}"
        fastqc --quiet ${reads}
        echo "[SUCCESS] FASTQC_PRE : fastqc GSM: ${meta.GSM}"
        """
}


// TODO fix the output of the trim galore in the process
process TRIM_GALORE {
    tag "$meta.GSM"
    publishDir "${params.pipeline_dir}/trimmed", mode: 'copy'
    errorStrategy 'ignore'

    input:
    val(meta)

    output:
    val(meta), emit: meta
    path("${meta.GSM}*.trim.fq"), emit: trimmed
    path "${meta.GSM}*_trimming_report.txt"

    script:
    def files  = meta.sp ? "-- paired ${params.fastq_dir}/${meta.GSM}_1 ${params.fastq_dir}/${meta.GSM}_2" : "${params.fastq_dir}/${meta.GSM}"
    """
    echo "[INFO] TRIM_GALORE_PAIRED : Starting trim of GSM: ${meta.GSM}"
    trim_galore \
    --fastqc \
    --trim-n \
    --length 20 \
    --quality 25 \
    ${meta.trimming_arg} \
    ${files}
    echo "[SUCCESS] TRIM_GALORE_PAIRED : Completed GSM: ${meta.GSM}"
    """
}


process FASTQC_POST {
    tag "$meta.GSM"
    publishDir "${params.pipeline_dir}/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/fastqc_post", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(meta), val(trim)

    output:
    val (meta), emit : meta
    path("*.html"), path("*.zip"), emit: files

    script:
        """
        echo "[INFO] FASTQC_POST : Starting fastqc of GSM: ${meta.GSM}"   
        fastqc --quiet ${trim}
        echo "[SUCCESS] FASTQC_POST : Downloaded GSM: ${meta.GSM}"     
        """
}
// TODO add multiqc template to split pre and post fastqc
process MULTIQC {
    tag "multiqc"
    publishDir "${params.pipeline_dir}/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}/multiqc", mode: 'copy'
    errorStrategy 'ignore'

    input:
    val (meta)

    output:
    path("multiqc_report.html"), emit: multiqc_report

    script:
    """
    echo "[INFO] MULTIQC : Starting MultiQC report generation"
    multiqc -o ${params.pipeline_dir}/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type} ${params.pipeline_dir}/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}
    echo "[SUCCESS] MULTIQC : MultiQC report generated"
    """
}



// Main workflow
workflow {
    log.info "Starting pipeline for existing files only..."
    // Channel that check for previously downloaded fastqs
    //TODO add a 2nd path to check for the old fastq files 
    Channel
        .fromPath(["${params.fastq_dir}/*.f*"])
        .map { file -> 
            file.simpleName.replaceFirst(/_[12]$/, '')
        }
        .unique()
        .collect().set { existing_fastq }

    // Main channel containing the CSV data, branched into two sub-channels depending
    // on whether the sample is already downloaded or not for the FASTQ_DUMP process
    // remerged later for the TRIM_GALORE and FASTQC processes
    Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .flatMap { row -> 
            row.Sample_accession
                .split(';')
                .collect{it.trim()}
                .collect { GSM -> 
                    [
                        GSE: row.Study_accession,
                        GSM: GSM,
                        drug: row.Drug,
                        sample_type: row.Biological_type,
                        trimming_args: row.Trim_arg,
                        paired_end: row.paired_end.toBoolean(),
                        sp: row.Species
                    ]
                }
        }
        .branch {
            existing: existing_fastq.val.contains(it.GSM)
            to_download: true
        }.set { csv_channel }

    FASTQ_DUMP(csv_channel.to_download)

    all_samples = FASTQ_DUMP.out.meta.mix(csv_channel.existing)

    TRIM_GALORE(all_samples)
    FASTQC_PRE(all_samples)
    TRIM_GALORE.out.meta.view { "TRIM_GALORE output: $it" }
    FASTQC_POST(TRIM_GALORE.out.meta)

    // regrouping channels to run multiqc on a multitude of sample fomr the same study
    TRIM_GALORE.out.meta
        .map { item -> [item.GSE + "_" + item.drug + '_' + item.sample_type, item] }
        .groupTuple()
        .map { key, items -> [
            GSE: items[0].GSE,
            drug: items[0].drug,
            GSMs: items.collect { it.GSM },
            riboseq_type: items[0].riboseq_type,
            sample_type: items[0].sample_type,
            sp: items[0].sp
        ]}
        .set { collapsed_ch }

    MULTIQC(collapsed_ch)
}
// add a way to log what failed


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
