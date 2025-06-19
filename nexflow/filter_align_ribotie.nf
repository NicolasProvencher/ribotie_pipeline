


process BOWTIE_ALIGN {
    beforeScript 'module load bowtie2'  
    // Maximum attempt of 3    
    // Dynamic strategy: retry up to the 3rd attempt, then ignore  
    tag "${gsm}"
    publishDir "${params.outdir_stage_stage}/bowtie/${gse}_${drug}_${bio}/${gsm}", mode: 'link', overwrite: true
       
    input:
    val(meta)


    output:
    val(meta), emit: meta
    path("*_filtered.fq"), emit: unmapped
    path("${gsm}_bowtie.log"), emit: log_file

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
    beforeScript 'module load star'  
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

workflow {
    Channel
        .fromPath(["${params.trimmed_fastq_dir}/*.fq"])
        .map { file -> 
            file.simpleName.split('_')[0]
        }
        .unique()
        .collect().set { trimmed_fastq }
    // Main channel containing the CSV data, branched into two sub-channels depending
    // on whether the file is present in the trimmed folder or not
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
            trimmed_files: trimmed_fastq.val.contains(it.GSM)
            missing: true
        }.set { csv_channel }
    csv_channel.missing.view { "CSV channel missing files: $it" }

    BOWTIE_ALIGN(csv_channel.trimmed_files)
    STAR_ALIGN(BOWTIE_ALIGN.out.meta)




}