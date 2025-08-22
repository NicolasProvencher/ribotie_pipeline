//TODO replace yml with actual arguments to simplify workflow
//TODO refer to the gtf and fasta using params.var[species]
//TODO arrange this for once per file
//TODO correclty output the h5 db
process RUN_RIBOTIE_DATA {
    beforeScript 'module load python cuda cudnn arrow '
    publishDir "${params.outdir_ribotie}/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/ribotie", mode: 'link', overwrite: true
    cache 'lenient'
    tag "${meta.GSM}"

    input:
    tuple val(meta), path(aligned_transcriptome)
    
    output:
    val(meta), emit: meta
    val(aligned_transcriptome), emit: aligned_transcriptome
    path("ribotie_data_log.txt")
    
    script:
    """
    # Prepare the environment
    virtualenv --no-download \$SLURM_TMPDIR/env
    source \$SLURM_TMPDIR/env/bin/activate
    pip install --no-index transcript_transformer

    cp /home/noxatras/scratch/ribotie/reference/HS/stuff.h5 \$SLURM_TMPDIR/${meta.GSM}.h5
    # Run RiboTIE with the YAML file as template
    ribotie --data \
        --gtf_path ${params.gtf_path} \
        --fa_path ${params.fa_path} \
        --h5_path \$SLURM_TMPDIR/${meta.GSM}.h5 \
        --out_prefix ${meta.GSM} \
        --ribo_paths ${ribopath} \ 
        --samples ${samples}


    cp \$SLURM_TMPDIR/${meta.GSM}.h5 \$PWD

    """
    // fix bam path and sample here 
}

// TODO correctly publish the csv results and the wieghts i guess, nah fuck the weights
process RUN_RIBOTIE {
    tag "${meta.GSM}"
    beforeScript 'module load python/3.9 cuda cudnn arrow'
    publishDir "${params.outdir_ribotie}/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/ribotie", mode: 'link', overwrite: true
    cache 'lenient'

    input:
    val(meta)
    
    output:
    tuple val(gse), val(drug), val(bio), path("*.csv"), emit: ribotie_results 
    path("ribotie_log.txt")
    
    script:
    """
    # Prepare the environment
    virtualenv --no-download \$SLURM_TMPDIR/env
    source \$SLURM_TMPDIR/env/bin/activate
    pip install --no-index transcript_transformer
    
    # Run RiboTIE with the YAML file as template
    ribotie ${yaml_file}  > ribotie_log.txt 2>&1

    """
}
// TODO adapt to file output from alignment processes
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