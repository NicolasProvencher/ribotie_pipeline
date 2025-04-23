#!/usr/bin/env nextflow

if (!params.outdir_first_stage) {
    log.error "Veuillez spécifier le chemin du dossier pipeline_output_complet avec --outdir_first_stage"
    log.error "Please specify the pipeline_output_complet directory path with --outdir_first_stage"
    exit 1
}

// Processus pour collecter et renommer les rapports MultiQC
// Process to collect and rename MultiQC reports
process COLLECT_MULTIQC_REPORTS {
    publishDir "${params.outdir_first_stage}/multiqc_rapport_complet", mode: 'copy', overwrite: true
    
    outdir_first_stage:
    tuple val(gse_name), path(multiqc_report)
    
    output:
    path "multiqc_report_${gse_name}.html"
    
    script:
    """
    cp ${multiqc_report} multiqc_report_${gse_name}.html
    """
}

workflow {
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

workflow.onComplete {
    log.info """
    Pipeline terminé avec succès!
    Les rapports MultiQC ont été copiés dans le dossier multiqc_rapport_complet
    Pipeline completed successfully!
    MultiQC reports have been copied to the multiqc_rapport_complet folder
    """
}

workflow.onError {
    log.error "Une erreur est survenue lors de l'exécution du pipeline"
    log.error "An error occurred during pipeline execution"
    log.error "Erreur: ${workflow.errorMessage}"
    log.error "Error: ${workflow.errorMessage}"
}
