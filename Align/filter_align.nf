

process BOWTIE_ALIGN {

    beforeScript 'module load bowtie2'  
    tag "${meta.GSM}"
    publishDir "${params.path_pipeline_directory}/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/bowtie", mode: 'link', overwrite: true
       
    input:
    val(meta)
    path(file)

    output:
    val(meta), emit: meta
    path("*_filtered*.fq.gz"), emit: unmapped
    path("${meta.GSM}_bowtie.log"), emit: log_file

    script:
    def files  = meta.paired_end ? "-1 ${params.trimmed_fastq_dir}/${meta.GSM}_1_trimmed.fq.gz -2 ${params.trimmed_fastq_dir}/${meta.GSM}_2_trimmed.fq.gz" : "-U ${params.trimmed_fastq_dir}/${meta.GSM}_trimmed.fq.gz"
    def out  = meta.paired_end ? "--un-conc-gz ${meta.GSM}_filtered_%.fq.gz" : "--un-gz ${meta.GSM}_filtered.fq.gz"
    
    """    
    # Run Bowtie2 with parameters for single-end and redirect to /dev/null
    # since we are only interested in unmapped reads
    bowtie2 \
        -p 4 \
        -t \
        -x ${params.bowtie_index[meta.sp]} \
        ${files} \
        ${out}
    """
}

process STAR_ALIGN {
    beforeScript 'module load star'
    tag "${meta.GSM}"
    publishDir "${params.outdir_stage_stage}/STAR/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/star", mode: 'link', overwrite: true

    input:
    val(meta)


    output:
    val(meta), emit: meta
    path("*Aligned.toTranscriptome.out.bam"), emit: aligned_transcriptome
    path("*Log.final.out"), emit: log
    path("*")

    script:
    def filtered_fq = meta.paired_end ? "${meta.GSM}_filtered_1.fq ${meta.GSM}_filtered_2.fq" : "${meta.GSM}_filtered.fq"
    """    
    STAR --runThreadN 5 \
         --genomeDir ${params.star_index[meta.sp]} \
         --genomeLoad NoSharedMemory \
         --readFilesIn ${filtered_fq} \
         --outFileNamePrefix ${meta.GSM}_ \
         --outSAMtype BAM Unsorted \
         --quantMode TranscriptomeSAM \
         --outSAMattributes MD NH \
         --outFilterMultimapNmax 10 \
         --outMultimapperOrder Random \
         --outFilterMismatchNmax 2 \
         --seedSearchStartLmaxOverLread 0.5 \
         --alignEndsType EndToEnd \
         --outWigType bedGraph
    """
}







workflow {
    trimmed_files_ch = Channel
        .fromPath("${params.trimmed_fastq_dir}/*.fq.gz")
        .map { file -> 
            def gsm = file.simpleName.split('_')[0]
            [GSM_ID:gsm, FILE:file]  // Named tuple
        }

    metadata_ch = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .flatMap { row -> 
            row.Sample_accession
                .split(';')
                .collect { it.trim() }
                .collect { gsm -> 
                    [GSM_ID:gsm, META:[
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

    // Use the 'by' parameter to specify the join key field
    joined_ch = metadata_ch.join(trimmed_files_ch, by: 'GSM_ID', remainder: true)
    // Result: [GSM_ID: gsm_value, META: metadata_map, FILE: file_or_null]
    
    // Split into matched and missing file entries
    matched_ch = joined_ch
        .filter { entry -> entry.FILE != null }
        .map { entry -> [GSM_ID:entry.GSM_ID, FILE:entry.FILE, META:entry.META] }

    missing_files_ch = joined_ch
        .filter { entry -> entry.FILE == null }
        .map { entry -> [entry.GSM_ID, entry.META] }
    
    // Report missing files
    missing_files_ch.view { gsm, metadata ->
        "WARNING: Missing trimmed file for GSM ${gsm} from study ${metadata.GSE} (${metadata.sample_type}, ${metadata.drug})"
    }
    
    // Create summary of missing files
    missing_summary_ch = missing_files_ch
        .collect()
        .map { missing_list ->
            if (missing_list.size() > 0) {
                def count = missing_list.size()
                def gsm_list = missing_list.collect { gsm, _metadata -> gsm }.join(', ')
                "SUMMARY: ${count} GSM samples missing trimmed files: ${gsm_list}"
            } else {
                "All samples have corresponding files!"
            }
        }
    
    missing_summary_ch.view()
    BOWTIE_ALIGN(matched_ch.META, matched_ch.FILE)
    STAR_ALIGN(BOWTIE_ALIGN.out.meta)
    
    // Use matched_ch for downstream processes
    // matched_ch emits: [GSM, file_path, metadata_map]
}