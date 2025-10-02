

process BOWTIE_ALIGN {
    cache 'lenient'
    beforeScript 'module load bowtie2'  
    tag "${meta.GSM}"
    publishDir "${params.path_pipeline_directory}/../output/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/bowtie", mode: 'link', overwrite: true
       
    input:
    tuple val(meta), path(file)

    output:
    val(meta), emit: meta
    path("*_filtered*.fq.gz"), emit: unmapped
    path("${meta.GSM}_bowtie.log"), emit: log_file

    script:
    def files  = meta.paired_end ? "-1 ${params.trimmed_fastq_dir}/${meta.GSM}_1_val_1.fq.gz -2 ${params.trimmed_fastq_dir}/${meta.GSM}_2_val_2.fq.gz" : "-U ${params.trimmed_fastq_dir}/${meta.GSM}_trimmed.fq.gz"
    def out  = meta.paired_end ? "--un-conc-gz ${meta.GSM}_filtered_%.fq.gz" : "--un-gz ${meta.GSM}_filtered.fq.gz"
    
    """    
    # Run Bowtie2 with parameters for single-end and redirect to /dev/null
    # since we are only interested in unmapped reads
    bowtie2 \
        -p 4 \
        -t \
        -x ${params.bowtie_filter_index[meta.sp]} \
        ${files} \
        ${out} > /dev/null 2> ${meta.GSM}_bowtie.log
    """
}

process STAR_ALIGN {
    cache 'lenient'
    beforeScript 'module load star'
    tag "${meta.GSM}"
    publishDir "${params.path_pipeline_directory}/../output/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/star", mode: 'link', overwrite: true

    input:
    val(meta)


    output:
    val(meta), emit: meta
    path("*Aligned.toTranscriptome.out.bam"), emit: aligned_transcriptome
    path("*Log.final.out"), emit: log
    path("*")

    script:
    def filtered_fq = meta.paired_end ? "${params.path_pipeline_directory}/../output/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/bowtie/${meta.GSM}_filtered_1.fq.gz ${params.path_pipeline_directory}/../output/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/bowtie/${meta.GSM}_filtered_2.fq.gz" : "${params.path_pipeline_directory}/../output/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/bowtie/${meta.GSM}_filtered.fq.gz"
    """    
    STAR --runThreadN 5 \
         --genomeDir ${params.star_index[meta.sp]} \
         --genomeLoad NoSharedMemory \
         --readFilesIn ${filtered_fq} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${meta.GSM}_ \
         --outSAMtype BAM Unsorted \
         --quantMode TranscriptomeSAM \
         --outSAMattributes MD NH \
         --outFilterMultimapNmax 10 \
         --outMultimapperOrder Random \
         --outFilterMismatchNmax 2 \
         --seedSearchStartLmaxOverLread 0.5 \
         --alignEndsType EndToEnd \
         --outTmpDir \$SLURM_TMPDIR/star_tmp \

    """
}







workflow {
    trimmed_files_ch = Channel
        .fromPath("${params.trimmed_fastq_dir}/*.fq.gz")
        .map { file -> 
            def gsm = file.simpleName.split('_')[0]
            [gsm, file]  // Named tuple
        }

    metadata_ch = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .flatMap { row -> 
            row.Sample_accession
                .split(';')
                .collect { it.trim() }
                .collect { gsm -> 
                    [gsm, [
                        GSE: row.Study_accession,
                        GSM: gsm,
                        drug: row.Drug,
                        sample_type: row.Biological_type,
                        trimming_args: row.Trim_arg,
                        paired_end: row.paired_end.toBoolean(),
                        sp: row.Species
                    ]]
                }
        }

    joined_ch = metadata_ch.join(trimmed_files_ch, by: 0, remainder: true)

    // joined_ch.view()

    matched_ch = joined_ch
        .filter { it[2] != null && it[1] != null }
        .map { [ it[1], it[2]] } // [META, FILE]

    missing_files_ch = joined_ch
        .filter { it[2] == null && it[1] != null }
        .map { [it[0], it[1]] } // 



    missing_files_ch.view { a, b -> "missing_files_ch: ${a} - ${b}" }
    BOWTIE_ALIGN(matched_ch)
    STAR_ALIGN(BOWTIE_ALIGN.out.meta)
    
    // Use matched_ch for downstream processes
    // matched_ch emits: [GSM, file_path, metadata_map]
}