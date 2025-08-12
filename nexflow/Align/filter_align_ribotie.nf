

process BOWTIE_ALIGN {

    beforeScript 'module load bowtie2'  
    // Maximum attempt of 3    
    // Dynamic strategy: retry up to the 3rd attempt, then ignore  
    tag "${meta.GSM}"
    publishDir "${params.path_pipeline_directory}/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/bowtie", mode: 'link', overwrite: true
       
    input:
    tuple val(meta), path(file)


    output:
    val(meta), emit: meta
    path("*_filtered.fq"), emit: unmapped
    path("${meta.GSM}_bowtie.log"), emit: log_file

    script:
    // Retrieve the index prefix (without extension)
    def files  = meta.paired_end ? "-1 ${params.trimmed_fastq_dir}/${meta.GSM}_1_trimmed.fq -2 ${params.trimmed_fastq_dir}/${meta.GSM}_2_trimmed.fq" : "-U ${params.trimmed_fastq_dir}/${meta.GSM}_trimmed.fq"
    def out  = meta.paired_end ? "--un-conc ${meta.GSM}_filtered_%.fq" : "--un ${meta.GSM}_filtered.fq"
    
    """    
    # Run Bowtie2 with parameters for single-end and redirect to /dev/null
    # since we are only interested in unmapped reads
    bowtie2 \
        -p 4 \
        -x ${params.bowtie_index[meta.sp]} \
        ${files} \
        ${out}
    """
}

process STAR_ALIGN {
    beforeScript '''module load star'''  
    tag "${gsm}"
    publishDir "${params.outdir_stage_stage}/STAR/${gse}_${drug}_${bio}/${gsm}", mode: 'link'
       
    input:
    val(meta)

    output:
    val(meta), emit: meta
    path("*Aligned.sortedByCoord.out.bam"), emit: aligned_bam
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
         --outSAMtype BAM SortedByCoordinate \
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

//TODO replace yml with actual arguments to simplify workflow
//TODO refer to the gtf and fasta using params.var[species]
//TODO arrange this for once per file
//TODO correclty output the h5 db

process RUN_RIBOTIE_DATA {
    tag "${gse}_${drug}_${bio}"
    
    // Use appropriate resources
    // cpus 8
    // memory '80 GB * ${task.attempt}'
    // time '24h'
    // maxRetries 3
    beforeScript 'module load python cuda cudnn arrow '
    // clusterOptions = '--account=rrg-xroucou'
    
    // Publish results
    publishDir "${params.outdir_ribotie}/${gse}_${drug}_${bio}/results_data", mode: 'copy'
    // errorStrategy { params.ignore_ribotie_errors ? 'ignore' : 'retry' }
    
    input:
    tuple val(meta), path(aligned_transcriptome)
    
    output:
    val(meta), emit: meta
    val(aligned_transcriptome), emit: aligned_transcriptome
    path("ribotie_data_log.txt")
    
    script:
    """
    # Prepare the environment
    export PATH="\$PATH:${params.ribotie_dir}/bin"
    virtualenv --no-download \$SLURM_TMPDIR/env
    source \$SLURM_TMPDIR/env/bin/activate
    pip install --no-index transcript_transformer
    
    # Run RiboTIE with the YAML file as template
    ribotie ${yaml_file} --data > ribotie_data_log.txt 2>&1

    """
}

// TODO correctly publish the csv results and the wieghts i guess, nah fuck the weights
process RUN_RIBOTIE {
    tag "${gse}_${drug}_${bio}"
     // Use appropriate resources
    // cpus 8
    // memory '60 GB'
    // time '24h'
    // maxRetries 2
    beforeScript 'module load python/3.9 cuda cudnn arrow'
    // clusterOptions = '--account=def-xroucou --gres=gpu:1'
    // Publish results
    publishDir "${params.outdir_ribotie}/${gse}_${drug}_${bio}/results_run", mode: 'copy'
    // errorStrategy { params.ignore_ribotie_errors ? 'ignore' : 'retry' }
     
    input:
    tuple val(gse), val(drug), val(bio), path(yaml_file)
    
    output:
    tuple val(gse), val(drug), val(bio), path("*.csv"), emit: ribotie_results 
    path("ribotie_log.txt")
    
    script:
    """
    # Prepare the environment
    export PATH="\$PATH:${params.ribotie_dir}/bin"
    virtualenv --no-download \$SLURM_TMPDIR/env
    source \$SLURM_TMPDIR/env/bin/activate
    pip install --no-index transcript_transformer
    
    # Run RiboTIE with the YAML file as template
    ribotie ${yaml_file}  > ribotie_log.txt 2>&1

    """
}






workflow {
    // Keep your named tuples as they are
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
        .map { entry -> [entry.GSM_ID, entry.FILE, entry.META] }
    
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
                def gsm_list = missing_list.collect { gsm, metadata -> gsm }.join(', ')
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