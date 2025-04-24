#!/usr/bin/env nextflow

process validate_outdir_first_stage {
    when:
    !params.outdir_first_stage

    script:
    """
    echo "Please specify the pipeline_output_complet directory path with --outdir_first_stage" >&2
    exit 1
    """
}
// Process to collect and rename MultiQC reports
process COLLECT_MULTIQC_REPORTS {
    publishDir "${params.outdir_first_stage}/multiqc_rapport_complet", mode: 'copy', overwrite: true
    
    input:
    tuple val(gse_name), path(multiqc_report)
    
    output:
    path "multiqc_report_${gse_name}.html"
    
    script:
    """
    cp ${multiqc_report} multiqc_report_${gse_name}.html
    """
}

workflow {
    validate_outdir_first_stage()
    // Création du canal pour les rapports MultiQC
    // Creation of the channel for MultiQC reports
    Channel
        .fromPath("${params.outdir_first_stage}/GSE*/multiqc/multiqc_report.html")
        .map { file -> 
            def gse_name = file.getParent().getParent().getName()
            tuple(gse_name, file)
        }
        .set { multiqc_reports }
    
    // Collecte et renommage des rapports
    // Collection and renaming of reports
    COLLECT_MULTIQC_REPORTS(multiqc_reports)
}

// workflow.onComplete {
//     log.info """
//     Pipeline completed successfully!
//     MultiQC reports have been copied to the multiqc_rapport_complet folder
//     """
// }

// workflow.onError {
//     log.error "An error occurred during pipeline execution"
//     log.error "Error: ${workflow.errorMessage}"
// }
