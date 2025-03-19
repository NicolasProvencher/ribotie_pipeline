// Pipeline pour l'analyse de données de séquençage
// 1. Indexation avec bowtie et STAR
// 2. Alignement avec bowtie pour filtrer les ARN
// 3. Alignement des lectures non mappées avec STAR

nextflow.enable.dsl = 2
params.max_retries = 3

// Fonctions de vérification pour chaque type de fichier
def checkFastqExists(gse, gsm, drug, bio, sp) {
    def outputPath = "${params.input_output_fastq_second_stage}/${gse}_${drug}_${bio}/${gsm}"
    if (sp.toLowerCase() == "paired") {
        def file1 = new File("${outputPath}/${gsm}_1.fastq")
        def file2 = new File("${outputPath}/${gsm}_2.fastq")
        return file1.exists() && file2.exists()
    } else {
        def file = new File("${outputPath}/${gsm}.fastq")
        return file.exists()
    }
}

// Définition de la fonction pour créer les channels
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

// Fonctions auxiliaires pour obtenir les chemins des fichiers
def getFastaPathBowtie(type) {
    switch(type.toUpperCase()) {
        case 'HS': return params.path_fasta_HS_B
        case 'CE': return params.path_fasta_CE_B
        case 'DM': return params.path_fasta_DM_B
        case 'SC': return params.path_fasta_SC_B
        case 'DR': return params.path_fasta_DR_B
        case 'SM': return params.path_fasta_SM_B
        default: throw new RuntimeException("Type non reconnu: ${type}")
    }
}

def getFastaPathStar(type) {
    switch(type.toUpperCase()) {
        case 'HS': return params.path_fasta_HS
        case 'CE': return params.path_fasta_CE
        case 'DM': return params.path_fasta_DM
        case 'SC': return params.path_fasta_SC
        case 'DR': return params.path_fasta_DR
        case 'SM': return params.path_fasta_SM
        default: throw new RuntimeException("Type non reconnu: ${type}")
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
        default: throw new RuntimeException("Type non reconnu: ${type}")
    }
}

// Ajouter ces fonctions avec les autres fonctions de vérification
def checkFastaExists(path) {
    if (!path) {
        println "[AVERTISSEMENT] Le chemin du fichier FASTA n'est pas défini"
        return false
    }
    def file = new File(path)
    if (!file.exists()) {
        println "[ERREUR] Fichier FASTA manquant: ${path}"
        return false
    }
    println "[INFO] Fichier FASTA trouvé: ${path}"
    return true
}

def checkGtfExists(path) {
    if (!path) {
        println "[AVERTISSEMENT] Le chemin du fichier GTF n'est pas défini"
        return false
    }
    def file = new File(path)
    if (!file.exists()) {
        println "[ERREUR] Fichier GTF manquant: ${path}"
        return false
    }
    println "[INFO] Fichier GTF trouvé: ${path}"
    return true
}

process BOWTIE_INDEX {
    cpus 5
    memory '5 GB'
    time '1h'    
    beforeScript 'module load bowtie2'  
    maxRetries params.max_retries
    tag "$type"
    publishDir "${params.outdir_stage_stage_parent}/index_bowtie/${type}", mode: 'copy'

    input:
    tuple val(type), path(fasta_file)

    output: 
    tuple val(type), path("*.bt2"), emit: bowtie_index_output
    tuple val(type), val("${type}"), emit: index_prefix

    script:
    """
    echo "Indexation Bowtie pour ${type} avec le fichier ${fasta_file}"
    bowtie2-build ${fasta_file} ${type}
    
    # Vérification des index générés
    echo "Index Bowtie générés:"
    ls -la *.bt2
    
    echo "Indexation Bowtie terminée avec préfixe: ${type}"
    """
}

process STAR_INDEX {
    cpus 30
    memory '40 GB'
    time '3h'
    beforeScript 'module load star'  
    maxRetries params.max_retries
    tag "$type"
    publishDir "${params.outdir_stage_stage_parent}/index_STAR", mode: 'copy'

    input:
    tuple val(type), val(gtf), path(fasta_file)

    output:
    tuple val(type), path("${type}_star_index/*"), emit: star_index_output

    script:
    """
    echo "Indexation STAR pour ${type} avec le fichier ${fasta_file} et GTF ${gtf}"
    mkdir -p ${type}_star_index
    STAR --runThreadN 30 \
     --genomeDir ${type}_star_index \
     --runMode genomeGenerate \
     --genomeFastaFiles ${fasta_file} \
     --sjdbGTFfile ${gtf} \
     --genomeSAindexNbases 8
    echo "Indexation STAR terminée"
    """
}

// Définition du process bowtie_single
process BOWTIE_SINGLE {
    cpus 10
    memory '20 GB'
    time '3h'
    beforeScript 'module load bowtie2'  
    maxRetries params.max_retries
    tag "${gsm}"
    publishDir "${params.outdir_stage_stage}/bowtie/${gse}_${drug}_${bio}/${gsm}", mode: 'copy'
       
    input:
    tuple val(type), val(gse), val(gsm), val(drug), val(bio), val(sp), path(fastq_file), path(bowtie_indexes)

    output:
    tuple val(type), val(gse), val(gsm), val(drug), val(bio), path("*_unmapped.fq"), emit: unmapped
    path("${gsm}_bowtie.log"), emit: log_file

    script:
    // Récupérer le préfixe des index (sans l'extension)
    def index_prefix = bowtie_indexes[0].toString() - ~/\.[0-9]\.bt2$/
    
    """
    # S'assurer que le fichier FASTQ est correctement identifié
    FASTQ=\$(ls -1 ${fastq_file})
    
    echo "Using index: ${index_prefix}"
    echo "Input file: \$FASTQ"
    echo "Index files: ${bowtie_indexes.join(', ')}"
    
    # Exécuter Bowtie2 avec les paramètres pour single-end et rediriger vers /dev/null
    # puisque nous ne sommes intéressés que par les lectures non mappées
    bowtie2 \
        -p 10 \
        -x ${index_prefix} \
        -U \$FASTQ \
        --un ${gsm}_unmapped.fq \
        2> ${gsm}_bowtie.log
    """
}

process BOWTIE_PAIRED {
    cpus 20
    memory '50 GB'
    time '3h'
    stageInMode 'copy'  // Forcer la copie des fichiers plutôt que des liens symboliques
    beforeScript 'module load bowtie2'  
    maxRetries params.max_retries
    tag "${gsm}"
    publishDir "${params.outdir_stage_stage}/bowtie/${gse}_${drug}_${bio}/${gsm}", mode: 'copy'
       
    input:
    tuple val(type), val(gse), val(gsm), val(drug), val(bio), val(sp), path(fastq_file_1), path(fastq_file_2), path(bowtie_indexes)

    output:
    tuple val(type), val(gse), val(gsm), val(drug), val(bio), path("*_unmapped_*.fq"), emit: unmapped
    path("${gsm}_bowtie.log"), emit: log_file

    script:
    // Récupérer le préfixe des index (sans l'extension)
    def index_prefix = bowtie_indexes[0].toString() - ~/\.[0-9]\.bt2$/

    """
    bowtie2 \
        -p 20 \
        -x ${index_prefix} \
        -1 ${fastq_file_1} \
        -2 ${fastq_file_2} \
        --un-conc ${gsm}_unmapped_%.fq \
        2> ${gsm}_bowtie.log > /dev/null
    """
}


// Mise à jour du process STAR_SINGLE
process STAR_SINGLE {
    cpus 10
    memory '40 GB'
    time '3h'
    beforeScript 'module load star'  
    maxRetries params.max_retries
    tag "${gsm}"
    publishDir "${params.outdir_stage_stage}/STAR/${gse}_${drug}_${bio}/${gsm}", mode: 'copy'
       
    input:
    tuple val(type), val(gse), val(gsm), val(drug), val(bio), path(unmapped_fq), path(star_indexes)

    output:
    tuple val(type), val(gse), val(gsm), val(drug), val(bio), path("*Aligned.sortedByCoord.out.bam"), emit: aligned_bam
    tuple val(type), val(gse), val(gsm), val(drug), val(bio), path("*Aligned.toTranscriptome.out.bam"), emit: aligned_transcriptome
    tuple val(type), val(gse), val(gsm), val(drug), val(bio), path("*Log.final.out"), emit: log
    path("*")

    script:
    """
    echo "Traitement STAR pour ${gsm} avec ${unmapped_fq}"

     # Déterminer le répertoire parent contenant les fichiers d'index
    INDEX_DIR=\$(dirname \$(readlink -f ${star_indexes[0]}))
    
    STAR --runThreadN 10 \
         --genomeDir \$INDEX_DIR \
         --genomeLoad NoSharedMemory \
         --readFilesIn ${unmapped_fq} \
         --outFileNamePrefix ${gsm}_ \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM \
         --outSAMattributes MD NH \
         --outFilterMultimapNmax 10 \
         --outMultimapperOrder Random \
         --outFilterMismatchNmax 2 \
         --seedSearchStartLmaxOverLread 0.5 \
         --alignEndsType EndToEnd \
         --outWigType bedGraph
         
    echo "Traitement STAR terminé pour ${gsm}"
    echo "Vérification des fichiers générés:"
    ls -la
    """
}

// Mise à jour du process STAR_PAIRED
process STAR_PAIRED {
    cpus 20
    memory '70 GB'
    time '6h' 
    beforeScript 'module load star'  
    maxRetries params.max_retries
    tag "${gsm}"
    publishDir "${params.outdir_stage_stage}/STAR/${gse}_${drug}_${bio}/${gsm}", mode: 'copy'
       
    input:
    tuple val(type), val(gse), val(gsm), val(drug), val(bio), path(unmapped_fq), path(star_indexes)

    output:
    tuple val(type), val(gse), val(gsm), val(drug), val(bio), path("*Aligned.sortedByCoord.out.bam"), emit: aligned_bam
    tuple val(type), val(gse), val(gsm), val(drug), val(bio), path("*Aligned.toTranscriptome.out.bam"), emit: aligned_transcriptome
    tuple val(type), val(gse), val(gsm), val(drug), val(bio), path("*Log.final.out"), emit: log
    path("*")

    script:
    // Détecter le pattern des fichiers non mappés (peut varier selon la version de Bowtie2)
    def read1 = unmapped_fq.find { it.toString().contains("_1.fq") || it.toString().contains("_unmapped_1.fq") }
    def read2 = unmapped_fq.find { it.toString().contains("_2.fq") || it.toString().contains("_unmapped_2.fq") }
    
    if (!read1 || !read2) {
        // Fallback si le pattern est différent
        def files = unmapped_fq.sort()
        read1 = files[0]
        read2 = files[1]
    }
    
    """
    echo "Traitement STAR pour ${gsm} avec ${read1} et ${read2}"
     
     # Déterminer le répertoire parent contenant les fichiers d'index
    INDEX_DIR=\$(dirname \$(readlink -f ${star_indexes[0]}))
    
    STAR --runThreadN 10 \
         --genomeDir \$INDEX_DIR \
         --genomeLoad NoSharedMemory \
         --readFilesIn ${read1} ${read2} \
         --outFileNamePrefix ${gsm}_ \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM \
         --outSAMattributes MD NH \
         --outFilterMultimapNmax 10 \
         --outMultimapperOrder Random \
         --outFilterMismatchNmax 2 \
         --seedSearchStartLmaxOverLread 0.5 \
         --alignEndsType EndToEnd \
         --outWigType bedGraph
         
    echo "Traitement STAR terminé pour ${gsm}"
    echo "Vérification des fichiers générés:"
    ls -la
    """
}

// Workflow principal

workflow {
    // Créer le channel initial
    def csvChannelPerGsm = create_initial_channels()

    // Création d'un channel pour les fichiers existants
    def existing_fastq_files = csvChannelPerGsm
        .filter { it -> 
            checkFastqExists(it[0], it[1], it[2], it[3], it[5])
        }
        .map { gse, gsm, drug, bio, trim, sp, type ->
            def outputPath = "${params.input_output_fastq_second_stage
        }/${gse}_${drug}_${bio}/${gsm}/trimmed"
            if (sp.toLowerCase() == "paired") {
                def fastq1 = file("${outputPath}/${gsm}_1_val_1.fq") // a modifier a /${gsm}_1_trimmed.fq donc il faut adaper le premier stage
                def fastq2 = file("${outputPath}/${gsm}_2_val_2.fq")
                tuple(type, gse, gsm, drug, bio, sp, fastq1, fastq2)
            } else {
                def fastq = file("${outputPath}/${gsm}_trimmed.fq")
                tuple(type, gse, gsm, drug, bio, sp, fastq)
            }
        }

    // Séparation en deux channels distincts selon le type de séquençage
    def (singleChannel, pairedChannel) = existing_fastq_files.branch {
        single: it[5].toLowerCase() == 'single'
        paired: it[5].toLowerCase() == 'paired'
    }

    // Debug des channels d'entrée
    singleChannel.count().subscribe { count ->
        println "\n[INFO] Nombre d'échantillons single-end à traiter: ${count}\n"
    }
    
    pairedChannel.count().subscribe { count ->
        println "\n[INFO] Nombre d'échantillons paired-end à traiter: ${count}\n"
    }
    
    // Extraire les types uniques puis la verification si les fichiers necessaires sont presents comme paths
    def species_types = existing_fastq_files
        .map { it -> it[0] }  // Extraire le type (species)
        .unique()
        .filter { type -> 
            def fasta_bowtie = getFastaPathBowtie(type)
            def fasta_star = getFastaPathStar(type)
            def gtf_star = getGtfPath(type)
            
            def filesExist = checkFastaExists(fasta_bowtie) && 
                          checkFastaExists(fasta_star) && 
                          checkGtfExists(gtf_star)
            
            if (!filesExist) {
                println "[AVERTISSEMENT] Fichiers manquants pour le type: ${type}"
            }
            return filesExist
        }

    // Créer les channels pour Bowtie et STAR séparément pour l'indexation
    def bowtie_input = species_types
        .map { type ->
            tuple(type, file(getFastaPathBowtie(type)))
        }

    def star_input = species_types
        .map { type ->
            tuple(
                type, 
                file(getGtfPath(type)),
                file(getFastaPathStar(type))
            )
        }

    // Lancer les processus d'indexation 
    def bowtie_indices = BOWTIE_INDEX(bowtie_input)
    def star_indices = STAR_INDEX(star_input)

    // Debug des informations d'index
    bowtie_indices.bowtie_index_output.view { type, files ->
        println "[INFO] Index Bowtie généré pour: $type (${files.size()} fichiers)"
    }

    star_indices.star_index_output.view { type, files ->
        println "[INFO] Index STAR généré pour: $type (${files.size()} fichiers)"
    }

    // Transformer le channel d'index pour avoir le bon format pour la combinaison
    def bowtie_indices_formatted = bowtie_indices.bowtie_index_output

    // SingleChnnel contient les lignes de csv qui sont des single    
    // Combiner les channels avec les index correspondants - FIXED
    def bowtie_single_with_index = singleChannel
        .map { type, gse, gsm, drug, bio, sp, fastq -> 
            return tuple(type, gse, gsm, drug, bio, sp, fastq)
        }
        .combine(bowtie_indices_formatted, by: 0) // La commande combine lie entre les lignes de fichier csv et les fichiers produits par l'indexation 
    
    def bowtie_paired_with_index = pairedChannel
        .map { type, gse, gsm, drug, bio, sp, fastq1, fastq2 -> 
            return tuple(type, gse, gsm, drug, bio, sp, fastq1, fastq2)
        }
        .combine(bowtie_indices_formatted, by: 0)
    
    // Exécuter Bowtie sur les single-end et paired-end avec leurs index respectifs
    def bowtie_single_output = BOWTIE_SINGLE(bowtie_single_with_index)
    def bowtie_paired_output = BOWTIE_PAIRED(bowtie_paired_with_index)

    // Log des fichiers de sortie Bowtie
    bowtie_single_output.log_file.view { log ->
        println "[INFO] Log Bowtie single-end généré: ${log}"
    }

    bowtie_paired_output.log_file.view { log ->
        println "[INFO] Log Bowtie paired-end généré: ${log}"
    }
    
    // FIX: Modifier également le format pour les index STAR
    def star_indices_formatted = star_indices.star_index_output

    
    // Combiner les sorties unmapped de Bowtie avec les index STAR - FIXED
    def star_single_input = bowtie_single_output.unmapped
        .combine(star_indices_formatted, by: 0)
    
    def star_paired_input = bowtie_paired_output.unmapped
        .combine(star_indices_formatted, by: 0)
    
    // Exécuter STAR sur les fichiers non mappés
    def star_single_output = STAR_SINGLE(star_single_input)
    def star_paired_output = STAR_PAIRED(star_paired_input)
    
    // Log des résultats STAR
    star_single_output.aligned_bam.view { type, gse, gsm, drug, bio, bam ->
        println "[SUCCESS] STAR single-end alignement terminé pour: $gsm"
    }

    star_paired_output.aligned_bam.view { type, gse, gsm, drug, bio, bam ->
        println "[SUCCESS] STAR paired-end alignement terminé pour: $gsm"
    }
    
    // Ajout pour le log des fichiers alignés au transcriptome
    star_single_output.aligned_transcriptome.view { type, gse, gsm, drug, bio, bam ->
        println "[SUCCESS] STAR single-end transcriptome alignement terminé pour: $gsm"
    }

    star_paired_output.aligned_transcriptome.view { type, gse, gsm, drug, bio, bam ->
        println "[SUCCESS] STAR paired-end transcriptome alignement terminé pour: $gsm"
    }
}


// Mise à jour de la section onComplete
workflow.onComplete {
    // Comptabiliser les résultats
    def bowtie_count = 0
    def star_count = 0
    def transcriptome_count = 0
    
    // Vérifier les fichiers de sortie Bowtie
    def bowtie_dir = new File("${params.outdir_stage_stage}/bowtie")
    if (bowtie_dir.exists()) {
        bowtie_count = bowtie_dir.listFiles().findAll { it.isDirectory() }.size()
    }
    
    // Vérifier les fichiers de sortie STAR
    def star_dir = new File("${params.outdir_stage_stage}/STAR")
    if (star_dir.exists()) {
        star_count = star_dir.listFiles().findAll { it.isDirectory() }.size()
        
        // Compter les fichiers transcriptome
        transcriptome_count = 0
        star_dir.eachFileRecurse { file ->
            if (file.name.endsWith("Aligned.toTranscriptome.out.bam")) {
                transcriptome_count++
            }
        }
    }
    
    log.info """
    Pipeline execution summary
    -------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    
    Results:
    --------
    Bowtie alignment completed: $bowtie_count échantillons
    STAR alignment completed: $star_count échantillons
    Transcriptome alignments: $transcriptome_count fichiers
    
    Output directory: ${params.outdir_stage_stage}
    """
}
