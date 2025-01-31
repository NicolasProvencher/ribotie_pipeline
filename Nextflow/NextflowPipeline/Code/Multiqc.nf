#!/usr/bin/env nextflow

params.input = false
params.help = false

// Message d'aide
def helpMessage() {
    log.info"""
    Usage:
    nextflow run multiqc.nf --input /chemin/vers/pipeline_output
    
    Paramètres obligatoires:
        --input    Chemin vers le dossier pipeline_output contenant les résultats du pipeline
    """
}

if (params.help) {
    helpMessage()
    exit 0
}

if (!params.input) {
    log.error "Veuillez spécifier le chemin du dossier pipeline_output avec --input"
    exit 1
}

// Processus principal pour MultiQC
process MULTIQC_GSE {
    publishDir "${params.input}/${gse_dir.baseName}/multiqc", mode: 'copy', overwrite: true    tag "${gse_dir.baseName}"
    
    input:
    path gse_dir
    
    output:
    path "multiqc_report.html"
    path "multiqc_report_data"
    
    script:
    """

    echo "Traitement de ${gse_dir.baseName}..."
    
    # Collecter tous les rapports FastQC et logs de trimming du dossier GSE
    multiqc . \
        --filename multiqc_report.html \
        --force \
        --interactive \
        --dirs \
        --dirs-depth 3
    """
}

workflow {
    // Création du canal pour les dossiers GSE
    Channel
        .fromPath("${params.input}/GSE*", type: 'dir')
        .set { gse_dirs }
    
    // Exécution de MultiQC pour chaque dossier GSE
    MULTIQC_GSE(gse_dirs)
}

workflow.onComplete {
    log.info """
    Pipeline terminé avec succès!
    Rapports MultiQC générés dans chaque dossier GSE
    """
}

workflow.onError {
    log.error "Une erreur est survenue lors de l'exécution du pipeline"
    log.error "Erreur: ${workflow.errorMessage}"
}
