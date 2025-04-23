#!/usr/bin/env nextflow

// DSL version
nextflow.enable.dsl = 2

// Specific parameters for RiboTIE
// params.star_dir = '/home/ilyass09/scratch/riboseq_pipeline/Second_Stage_copy/HS/STAR'
// params.input_csv = '/home/ilyass09/scratch/riboseq_pipeline/Samples_sheet/test.csv'
// params.outdir_ribotie = '/home/ilyass09/scratch/riboseq_pipeline/Ribotie_complet'
// params.ribotie_dir = '/home/ilyass09/scratch/riboseq_pipeline/ribotie'

// params.ignore_ribotie_errors = true

// Helper functions to get file paths
def getFastaPath(type) {
    switch(type.toUpperCase()) {
        case 'HS': return params.path_fasta_HS
        case 'CE': return params.path_fasta_CE
        case 'DM': return params.path_fasta_DM
        case 'SC': return params.path_fasta_SC
        case 'DR': return params.path_fasta_DR
        case 'SM': return params.path_fasta_SM
        default: throw new RuntimeException("Unrecognized type: ${type}")
    }
}

def getGtfPath(type) {
    switch(type.toUpperCase()) {
        case 'HS': return params.path_fasta_HS_GTF
        case 'CE': return params.path_fasta_CE_GTF
        case 'DM': return params.path_fasta_DM_GTF
        case 'SC': return params.path_fasta_SC_GTF
        case 'DR': return params.path_fasta_DR_GTF
        case 'SM': return params.path_fasta_SM_GTF
        default: throw new RuntimeException("Unrecognized type: ${type}")
    }
}

// TODO: I don't think I will use CSV for this separate file
// Definition of function to create channels
def create_initial_channels() {
    def csvChannelPerLigne = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> tuple(
            row.Study_accession,
            row.Sample_accession,
            row.Drug,
            row.Biological_type,
            row.Trim_arg,
            row.S_P_type,
            row.Species
        )}

    return csvChannelPerLigne
        .map { gse, gsm, drug, bio, trim, sp, type -> tuple(
            gse,
            gsm.split(';').collect{it.trim()},
            drug,
            bio,
            trim,
            sp,
            type
        )}
        .transpose(by: 1)
}

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

// Process to create YAML file per GSE
process CREATE_GSE_YAML {
    tag "${gse}_${drug}_${bio}"
    cpus 1
    memory '1 GB'
    time '30m'
    
    // Try with 'copy' mode and 'create: true'
    publishDir "${params.outdir_ribotie}/${gse}_${drug}_${bio}", mode: 'copy', create: true
    
    input:
    tuple val(gse), val(gsm_list), val(drug), val(bio), val(bam_list), val(type)
    
    output:
    tuple val(gse), val(drug), val(bio), path("${gse}_${drug}_${bio}.yml"), emit: gse_yaml_file
    
    script:
    def gtf_path = getGtfPath(type)
    def fa_path = getFastaPath(type)
    def h5_dir = "${params.outdir_ribotie}/${gse}_${drug}_${bio}/ribotie_data_results"
    def ribo_paths_str = ""
    
    // Create entries for ribo_paths
    for (int i = 0; i < gsm_list.size(); i++) {
        ribo_paths_str += "  ${gsm_list[i]}: ${bam_list[i]}\n"
    }
    
    // Create sample groups with the specific format required
    def samples_str = "samples:\n"
    samples_str += "  - - ${gsm_list[0]}\n"
    
    for (int i = 1; i < gsm_list.size(); i++) {
        samples_str += "    - ${gsm_list[i]}\n"
    }
    
    """
    # Explicitly create output directories
    mkdir -p ${params.outdir_ribotie}/${gse}_${drug}_${bio}
    mkdir -p ${h5_dir}
    
    cat > ${gse}_${drug}_${bio}.yml << EOF
gtf_path: ${gtf_path}
fa_path: ${fa_path}

ribo_paths:
${ribo_paths_str}

${samples_str}
h5_path: ${h5_dir}/${gse}_${drug}_${bio}.h5

out_dir: ${params.outdir_ribotie}/${gse}_${drug}_${bio}/results
EOF
    
    # Copy file explicitly if needed
    cp ${gse}_${drug}_${bio}.yml ${params.outdir_ribotie}/${gse}_${drug}_${bio}/
    """
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
    // Create the channel from STAR folders
    def bam_channel = create_star_bam_channels()
    
    // Display the content of the channel
    bam_channel.view { gse, gsm, drug, bio, bam_path ->
        return "GSE: $gse, GSM: $gsm, Drug: $drug, Bio: $bio, BAM: $bam_path"
    }
    
    // Group BAMs by GSE/drug/bio
    def grouped_bams = bam_channel
        .map { gse, gsm, drug, bio, bam_path -> 
            def key = [gse, drug, bio]
            return tuple(key, tuple(gsm, bam_path))
        }
        .groupTuple()
        .map { key, values -> 
            def (gse, drug, bio) = key
            def gsm_list = values.collect { it[0] }
            def bam_list = values.collect { it[1] }
            // Add a default type (you can adapt it as needed)
            // TODO: We are working directly on HS because RiboTIE is adapted for HS, so we need to think about modifying the input path structure
            // instead of STAR/HS to STAR directly
            def type = "HS"  // Default to Human
            return tuple(gse, gsm_list, drug, bio, bam_list, type)
        }
    
    // Display the groups
    grouped_bams.view { gse, gsm_list, drug, bio, bam_list, type ->
        return "GSE: $gse, Drug: $drug, Bio: $bio, Type: $type\nGSMs: ${gsm_list.join(', ')}\nBAMs: ${bam_list.join(', ')}"
    }
    
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
