#!/usr/bin/env nextflow


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

// Processus pour trouver tous les GSE dans le répertoire d'entrée
process FIND_GSE_DIRS {
    output:
    path "gse_dirs.txt", emit: gse_dirs

    script:
    """
    find ${params.outdir_first_stage} -mindepth 1 -maxdepth 1 -type d | sort > gse_dirs.txt
    echo "Répertoires GSE trouvés:" >&2
    cat gse_dirs.txt >&2
    """
}

// Processus pour générer la liste des fichiers FastQC pour chaque GSE
process GENERATE_FASTQC_LIST {
    input:
    path gse_dir

    output:
    tuple val(gse_name), path("${gse_name}_fastqc_files.txt")

    script:
    gse_name = gse_dir.getName()
    """
    # Trouver tous les fichiers FastQC HTML avec chemins absolus
    find ${gse_dir} -type f -name "*_fastqc.html" -o -name "*_trimmed_fastqc.html" | sort > ${gse_name}_fastqc_files.txt
    
    # Vérifier que des fichiers ont été trouvés
    if [ ! -s ${gse_name}_fastqc_files.txt ]; then
        echo "ATTENTION: Aucun fichier FastQC trouvé pour ${gse_name}" >&2
        echo "# Aucun fichier FastQC trouvé pour ${gse_name}" > ${gse_name}_fastqc_files.txt
    else
        NUM_FILES=\$(wc -l < ${gse_name}_fastqc_files.txt)
        echo "Fichiers FastQC trouvés pour ${gse_name}: \$NUM_FILES" >&2
    fi
    """
}

// Processus pour exécuter MultiQC sur chaque GSE et organiser par dossier
process RUN_MULTIQC {
    // Publier dans le dossier multiqc du GSE correspondant
    publishDir "${params.outdir_first_stage}/${gse_name}/multiqc", mode: 'copy', overwrite: true
    
    input:
    tuple val(gse_name), path(fastqc_list)
    
    output:
    tuple val(gse_name), path("multiqc_report.html"), optional: true
    tuple val(gse_name), path("multiqc_report_data"), optional: true
    
    script:
    """
    # Vérifier si la liste des fichiers est vide
    if [ ! -s ${fastqc_list} ] || grep -q "^#" ${fastqc_list}; then
        echo "Aucun fichier FastQC trouvé pour ${gse_name}, aucun rapport MultiQC n'est généré." >&2
        mkdir -p multiqc_report_data
        touch multiqc_report.html
    else
        echo "Exécution de MultiQC pour ${gse_name}" >&2
        
        # Exécuter MultiQC avec la liste des fichiers FastQC pour ce GSE
        multiqc -c ${params.multiqc_config} \\
            --filename multiqc_report.html \\
            -f \\
            --file-list ${fastqc_list}
        
        if [ ! -f "multiqc_report.html" ]; then
            echo "ERREUR: Échec de génération du rapport MultiQC pour ${gse_name}" >&2
            mkdir -p multiqc_report_data
            touch multiqc_report.html
        else
            echo "Rapport MultiQC créé avec succès pour ${gse_name}" >&2
            echo "Il sera publié dans: ${params.outdir_first_stage}/${gse_name}/multiqc/" >&2
        fi
    fi
    """
}

workflow {
    // Trouver tous les répertoires GSE
    FIND_GSE_DIRS()
    
    // Lire les répertoires GSE et les traiter comme un canal
    gse_dirs_ch = FIND_GSE_DIRS.out.gse_dirs
        .splitText()
        .map { it.trim() }
        .filter { it != "" }
    
    // Générer la liste des fichiers FastQC pour chaque GSE
    fastqc_lists_ch = GENERATE_FASTQC_LIST(gse_dirs_ch)
    
    // Exécuter MultiQC sur chaque GSE
    RUN_MULTIQC(fastqc_lists_ch)
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
