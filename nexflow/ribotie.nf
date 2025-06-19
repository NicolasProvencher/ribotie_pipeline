#!/usr/bin/env nextflow

// DSL version
nextflow.enable.dsl = 2



// Process to run RiboTIE on each GSE YAML file
process RUN_RIBOTIE_DATA {
    tag "${gse}_${drug}_${bio}"
    
    // Use appropriate resources
    // cpus 8
    // memory '80 GB * ${task.attempt}'
    // time '24h'
    // maxRetries 3
    beforeScript 'module load python/3.9 cuda cudnn arrow'
    // clusterOptions = '--account=def-xroucou --gres=gpu:1'
    
    // Publish results
    publishDir "${params.outdir_ribotie}/${gse}_${drug}_${bio}/results_data", mode: 'copy'
    // errorStrategy { params.ignore_ribotie_errors ? 'ignore' : 'retry' }
    
    input:
    tuple val(gse), val(drug), val(bio), path(yaml_file)
    
    output:
    tuple val(gse), val(drug), val(bio) ,path(yaml_file), emit: ribotie_results_data 
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
process RUN_RIBOTIE {
    tag "${gse}_${drug}_${bio}"
     // Use appropriate resources
    // cpus 8
    // memory '60 GB'
    // time '24h'
    // maxRetries 2
    beforeScript 'module load python/3.9 cuda cudnn arrow'
    clusterOptions = '--account=def-xroucou --gres=gpu:1'
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
    
    // Display the results
    RUN_RIBOTIE_DATA.out.ribotie_results_data.view { gse, drug, bio, result_files ->
        return "RiboTIE completed for GSE: $gse, Drug: $drug, Bio: $bio"
    }
    RUN_RIBOTIE(csv_files.ribotie_results_data)
    // Display the results
    RUN_RIBOTIE.out.ribotie_results.view { gse, drug, bio, result_files ->
        return "RiboTIE completed for GSE: $gse, Drug: $drug, Bio: $bio"
    }
}
