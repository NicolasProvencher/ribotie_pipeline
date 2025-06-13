#!/usr/bin/env nextflow

// DSL version
nextflow.enable.dsl = 2

// Specific parameters for RiboTIE
// params.star_dir = '/path/to/star/pipeline/directory/STAR'
// params.input_csv = '/path/to/pipeline/Samples_sheet/test.csv'
// params.outdir_ribotie = '/path/to/ribotie/output/directory'
// params.ribotie_dir = '/path/to/ribotie/directory'

// params.ignore_ribotie_errors = true



// TODO: I don't think I will use CSV for this separate file
// Definition of function to create channels


// Function to create channel from STAR folders (simplified version)
def create_star_bam_channels() {
    def results = []
    
    // Get all GSE folders
    def starDir = new File(params.star_dir)
    if (!starDir.exists()) {
        log.error "STAR directory does not exist: ${params.star_dir}"
        return Channel.empty()
    }
    
    // Iterate through GSE folders
    starDir.eachDir { gseDir ->
        def gseDirName = gseDir.name
        // Extract GSE, drug, and bio_type from folder name
        def parts = gseDirName.split('_')
        if (parts.length >= 3) {
            def gse = parts[0]
            def drug = parts[1]
            def bio = parts.length > 2 ? parts[2..-1].join('_') : ""
            
            // Iterate through GSM folders
            gseDir.eachDir { gsmDir ->
                def gsm = gsmDir.name
                // Look for BAM file of the transcriptome
                def bamFile = new File(gsmDir, "${gsm}_Aligned.toTranscriptome.out.bam")
                if (bamFile.exists()) {
                    // Add each GSM entry directly to the result
                    results.add(tuple(gse, gsm, drug, bio, bamFile.absolutePath))
                }
            }
        }
    }
    
    // Create and return the channel
    return Channel.fromList(results)
}


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
    
    
    // Group BAMs by GSE/drug/bio
    
    // Create YAML files per GSE
    def gse_yaml_files = CREATE_GSE_YAML(grouped_bams)
    
    // Run RiboTIE on each GSE YAML file
    def csv_files = RUN_RIBOTIE_DATA(gse_yaml_files.gse_yaml_file)
    
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
