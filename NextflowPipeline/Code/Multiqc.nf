#!/usr/bin/env nextflow

params.outdir_first_stage = "/project/def-xroucou/riboseq_pipeline/Etape_final_GSE63570/HS"
params.multiqc_config = "/project/def-xroucou/riboseq_pipeline/config_multiqc.yaml"

log.info """\
=================================================================
MultiQC par GSE - Nextflow Pipeline (Organisation par dossier)
=================================================================
Répertoire d'entrée     : ${params.outdir_first_stage}
Répertoire de sortie    : ${params.outdir_first_stage}
Config MultiQC          : ${params.multiqc_config}
=================================================================
"""

// Trouver tous les répertoires GSE
process FIND_GSE_DIRS {
    output:
    path "gse_dirs.txt", emit: gse_dirs

    script:
    """
    find ${params.outdir_first_stage} -mindepth 1 -maxdepth 1 -type d | sort > gse_dirs.txt
    """
}

// Exécuter MultiQC directement dans chaque dossier GSE
process RUN_MULTIQC {
    publishDir "${params.outdir_first_stage}/${gse_name}/multiqc", mode: 'copy', overwrite: true
    
    conda "/home/ilyass09/project/pipeline/code/env.yml"
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
    
    # Se déplacer dans le répertoire GSE et exécuter MultiQC directement
    cd ${gse_dir}
    
    # Exécuter MultiQC pour analyser tous les fichiers FastQC dans le répertoire
    multiqc -v -f -c ${params.multiqc_config} \
        --filename multiqc_report.html \
        --ignore "*/work/*" \
        --ignore "*/multiqc/*" .
    
    # Vérifier le résultat
    if [ ! -f "multiqc_report.html" ]; then
        touch multiqc_report.html
        mkdir -p multiqc_data
    fi
    """
}

workflow {
    // Trouver tous les répertoires GSE
    FIND_GSE_DIRS()
    
    // Créer un canal à partir des chemins des répertoires GSE
    gse_dirs_ch = FIND_GSE_DIRS.out.gse_dirs
        .splitText()
        .map { it.trim() }
        .filter { it != "" }
    
    // Exécuter MultiQC sur chaque GSE
    RUN_MULTIQC(gse_dirs_ch)
}

// Notification de fin et récapitulatif
workflow.onComplete {
    log.info """
    =================================================================
    Exécution terminée!
    
    Résumé:
    - Statut: ${workflow.success ? "SUCCÈS" : "ÉCHEC"}
    - Répertoire de sortie: ${params.outdir_first_stage}
    
    Les rapports MultiQC ont été organisés dans:
    [GSE]/multiqc/multiqc_report.html
    =================================================================
    """
}
