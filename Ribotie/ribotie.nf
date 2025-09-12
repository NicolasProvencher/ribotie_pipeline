//TODO replace yml with actual arguments to simplify workflow
//TODO refer to the gtf and fasta using params.var[species]
//TODO arrange this for once per file
//TODO correclty output the h5 db
process RIBOTIE_DATA {
    beforeScript 'module load python/3.11 cuda cudnn arrow '
    publishDir "${projectDir}/../output/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/ribotie", mode: 'link', overwrite: true
    cache 'lenient'
    tag "${meta.GSM}"

    input:
    tuple val(meta), file(Transcriptome_Bam)
    
    output:
    val(meta), emit: meta
    path("${meta.GSM}.h5"), emit: h5_db
    path(Transcriptome_Bam), emit: transcriptome_bam
    
    script:
    """
    # Prepare the environment
    virtualenv --no-download \$SLURM_TMPDIR/env
    source \$SLURM_TMPDIR/env/bin/activate
    pip install --no-index ${params.ribotie_package_path}
    mkdir -p \$SLURM_TMPDIR/${meta.GSM}
    cp ${params.reference_files_directory}/HS/Homo_sapiens.GRCh38.114.h5 \$SLURM_TMPDIR/${meta.GSM}/${meta.GSM}.h5

    ribotie --data \
        --gtf_path ${params.annotation_GTF[meta.sp]} \
        --fa_path ${params.dna_assembly[meta.sp]} \
        --h5_path \$SLURM_TMPDIR/${meta.GSM}/${meta.GSM}.h5 \
        --out_prefix ${meta.GSM} \
        --ribo_paths '{"${meta.GSM}": "${Transcriptome_Bam}"}' \
        --samples ${meta.GSM} \

    ls \$SLURM_TMPDIR/${meta.GSM}
    cp \$SLURM_TMPDIR/${meta.GSM}/* \$PWD

    """
    // fix bam path and sample here 
}

// TODO correctly publish the csv results and the wieghts i guess, nah fuck the weights
process RIBOTIE_ML {
    tag "${meta.GSM}"
    beforeScript 'module load python/3.11 cuda cudnn arrow'
    publishDir "${projectDir}/../output/${meta.sp}/${meta.GSE}_${meta.drug}_${meta.sample_type}/${meta.GSM}/ribotie/results", mode: 'link', overwrite: true
    cache 'lenient'

    input:
    val(meta)
    path(h5_db)
    path(Transcriptome_Bam)


    output:
    val(meta), emit: meta
    path("*")
    
    script:
    """
    # Prepare the environment
    virtualenv --no-download \$SLURM_TMPDIR/env
    source \$SLURM_TMPDIR/env/bin/activate
    pip install --no-index ${params.ribotie_package_path}
    mkdir -p \$SLURM_TMPDIR/${meta.GSM}
    cp ${h5_db} \$SLURM_TMPDIR/${meta.GSM}/${meta.GSM}.h5

    ribotie \
        --gtf_path ${params.annotation_GTF[meta.sp]} \
        --fa_path ${params.dna_assembly[meta.sp]} \
        --h5_path \$SLURM_TMPDIR/${meta.GSM}/${meta.GSM}.h5 \
        --out_prefix ${meta.GSM} \
        --ribo_paths '{"${meta.GSM}": "${Transcriptome_Bam}"}' \
        --samples ${meta.GSM} \


    ls \$SLURM_TMPDIR/${meta.GSM}
    cp \$SLURM_TMPDIR/${meta.GSM}/* \$PWD
    """
}
// TODO adapt to file output from alignment processes
workflow {
    trimmed_files_ch = Channel
        .fromPath("${projectDir}/../output/*/*/*/star/*_Aligned.toTranscriptome.out.bam")
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

    matched_ch = joined_ch
        .filter { it[2] != null && it[1] != null }
        .map { [ it[1], it[2]] } // [META, FILE]

    missing_files_ch = joined_ch
        .filter { it[2] == null && it[1] != null }
        .map { [it[0], it[1]] } // [GSM, META]



    missing_files_ch.view { a, b -> "missing_files_ch: ${a} - ${b}" }
    RIBOTIE_DATA(matched_ch)
    RIBOTIE_ML(RIBOTIE_DATA.out)
}