#!/usr/bin/env nextflow

params.outdir_first_stage = "/path/to/first_stage/output/directory"
params.multiqc_config = "/path/to/Ribotie_nextflow_pipeline/multiqc_custom_config.yaml"

log.info """\
=================================================================
MultiQC by GSE - Nextflow Pipeline (Folder Organization)
=================================================================
Input directory      : ${params.outdir_first_stage}
Output directory     : ${params.outdir_first_stage}
MultiQC Config       : ${params.multiqc_config}
=================================================================
"""

// Find all GSE directories
process FIND_GSE_DIRS {
    output:
    path "gse_dirs.txt", emit: gse_dirs

    script:
    """
    find ${params.outdir_first_stage} -mindepth 1 -maxdepth 1 -type d | sort > gse_dirs.txt
    """
}

// Run MultiQC directly in each GSE folder
process RUN_MULTIQC {
    publishDir "${params.outdir_first_stage}/${gse_name}/multiqc", mode: 'copy', overwrite: true
    
    conda "/path/to/conda/env_yaml/env.yml"
    memory '4 GB'
    
    input:
    val gse_dir
    
    output:
    path "multiqc_report.html", optional: true
    path "multiqc_data", optional: true
    
    script:
    gse_name = file(gse_dir).getName()
    """
    export CONDA_NO_PLUGINS=true
    
    # Move to the GSE directory and run MultiQC directly
    cd ${gse_dir}
    
    # Run MultiQC to analyze all FastQC files in the directory
    multiqc -v -f -c ${params.multiqc_config} \
        --filename multiqc_report.html \
        --ignore "*/work/*" \
        --ignore "*/multiqc/*" .
    
    # Check the result
    if [ ! -f "multiqc_report.html" ]; then
        touch multiqc_report.html
        mkdir -p multiqc_data
    fi
    """
}

workflow {
    // Find all GSE directories
    FIND_GSE_DIRS()
    
    // Create a channel from the paths of the GSE directories
    gse_dirs_ch = FIND_GSE_DIRS.out.gse_dirs
        .splitText()
        .map { it.trim() }
        .filter { it != "" }
    
    // Run MultiQC on each GSE
    RUN_MULTIQC(gse_dirs_ch)
}

// End notification and summary
workflow.onComplete {
    log.info """
    =================================================================
    Execution completed!
    
    Summary:
    - Status: ${workflow.success ? "SUCCESS" : "FAILURE"}
    - Output directory: ${params.outdir_first_stage}
    
    MultiQC reports have been organized in:
    [GSE]/multiqc/multiqc_report.html
    =================================================================
    """
}
