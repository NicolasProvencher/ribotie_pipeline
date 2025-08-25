process FASTQ_DUMP {
    tag "$meta.GSM"
    publishDir "${params.fastq_dir}", mode: 'link'

    input:
    val (meta)
    
    output:
    val (meta), emit: meta
    path("${meta.GSM}*.fastq"), emit: fastq_file

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
    publishDir "${params.path_pipeline_directory}/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/fastqc_pre", mode: 'link', overwrite: true
    conda "fastqc"
    input:
    val (meta)

    output:
    val (meta), emit : meta
    tuple path("*.html") , path("*.zip"), emit : files



    script:
    def reads = meta.paired_end ? "${params.fastq_dir}/${meta.GSM}_1.fastq ${params.fastq_dir}/${meta.GSM}_2.fastq" : "${params.fastq_dir}/${meta.GSM}.fastq"
        """
        echo "[INFO] FASTQC_PRE : Starting fastqc of GSM: ${meta.GSM}"
        pwd
        echo "fastqc -o \$(pwd) ${reads}"
        fastqc -o \$(pwd) ${reads}
        echo "[SUCCESS] FASTQC_PRE : fastqc GSM: ${meta.GSM}"
        """
}


// TODO fix the output of the trim galore in the process
process TRIM_GALORE {
    tag "$meta.GSM"
    publishDir "${params.trimmed_fastq_dir}", mode: 'link', overwrite: true

    publishDir "${params.path_pipeline_directory}/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/trimmed", mode: 'link', overwrite: true
    conda "trim-galore"
    input:
    val(meta)

    output:
    val(meta), emit: meta
    path("${meta.GSM}*.fq"), emit: trimmed
    path("${meta.GSM}*_trimming_report.txt")



    script:
    def files  = meta.paired_end ? "--paired ${params.fastq_dir}/${meta.GSM}_1.fastq ${params.fastq_dir}/${meta.GSM}_2.fastq" : "${params.fastq_dir}/${meta.GSM}.fastq"
    def trimming_arg = meta.trimming_args ?: "" 
    """
    echo "[INFO] TRIM_GALORE_PAIRED : Starting trim of GSM: ${meta.GSM}"
    echo "Trimming arguments: ${trimming_arg}, meta.trimming_args: ${meta.trimming_args}"
    trim_galore \
    --length 20 \
    --quality 25 \
    --max_length 39 \
    ${trimming_arg} \
    ${files}
    echo "[SUCCESS] TRIM_GALORE_PAIRED : Completed GSM: ${meta.GSM}"
    """
}


process FASTQC_POST {
    tag "$meta.GSM"
    publishDir "${params.path_pipeline_directory}/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/fastqc_post", mode: 'link', overwrite: true
    conda "fastqc"

    input:
    val(meta)
    val(trim)

    output:
    val (meta), emit : meta
    tuple path("*.html"), path("*.zip"), emit: files

    

    script:
    def files  = meta.paired_end ? "${trim[0]} ${trim[1]}" : "${trim}"
        """
        echo "[INFO] FASTQC_POST : Starting fastqc of GSM: ${meta.GSM}"   
        echo "fastqc -o \$(pwd) ${files}"
        fastqc -o \$(pwd) ${files}
        echo "[SUCCESS] FASTQC_POST : Downloaded GSM: ${meta.GSM}"     
        """
}
// TODO add multiqc template to split pre and post fastqc
process MULTIQC {
    tag "multiqc"
    publishDir "${params.path_pipeline_directory}/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}/multiqc", mode: 'link', overwrite: true
    conda "multiqc"
    input:
    val (meta)

    output:
    path("multiqc_report.html"), emit: multiqc_report



    script:
    """
    echo "[INFO] MULTIQC : Starting MultiQC report generation"
    multiqc -o \$(pwd) -c ${params.multiqc_config} ${params.path_pipeline_directory}/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}
    echo "[SUCCESS] MULTIQC : MultiQC report generated"
    """
}



// Main workflow
workflow {
    // Channel that check for previously downloaded fastqs
    //TODO add a 2nd path to check for the old fastq files 
    Channel
        .fromPath(["${params.fastq_dir}/*.f*"])
        .map { file -> 
            [file.simpleName.replaceFirst(/_[12]$/, ''), true]
        }
        .unique()
        .set { existing_fastq }

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
                .collect { GSM -> [GSM,
                    [
                        GSE: row.Study_accession,
                        GSM: GSM,
                        drug: row.Drug,
                        sample_type: row.Biological_type,
                        trimming_args: row.Trim_arg,
                        paired_end: row.paired_end.toBoolean(),
                        sp: row.Species
                    ]]
                }
        }.set { csv_channel }


    joined=csv_channel.join(existing_fastq, by : 0, remainder: true)

    no_download = joined
        .filter { it[1] != null && it[2] != null }
        .map { it[1] } // [META]

    to_download = joined
        .filter { it[1] != null && it[2] == null }
        .map { it[1] } // [GSM]


    FASTQ_DUMP(to_download)

    all_samples = FASTQ_DUMP.out.meta.mix(no_download)

    TRIM_GALORE(all_samples)
    FASTQC_PRE(all_samples)
    // TRIM_GALORE.out.meta.view { "TRIM_GALORE output: $it" }
    FASTQC_POST(TRIM_GALORE.out.meta, TRIM_GALORE.out.trimmed)

    // regrouping channels to run multiqc on a multitude of sample fomr the same study
    FASTQC_POST.out.meta
        .map { item -> [item.GSE + "_" + item.drug + '_' + item.sample_type, item] }
        .groupTuple()
        .map { _key, items -> [
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