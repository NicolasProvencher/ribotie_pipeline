// Pipeline for sequencing data analysis
// 1. Indexing with bowtie and STAR
// 2. Alignment with bowtie to filter RNAs
// 3. Alignment of unmapped reads with STAR

nextflow.enable.dsl = 2
params.max_retries = 3
// The output_complet folder is the result of the first stage
params.input_output_fastq_second_stage = '/home/ilyass09/scratch/riboseq_pipeline/Etape_final_trimmed/HS'
params.input_csv_test = '/home/ilyass09/scratch/riboseq_pipeline/Samples_sheet/Samples_sheet_final.csv'
params.outdir_stage_stage = '/home/ilyass09/scratch/riboseq_pipeline/Second_Stage_HS_final_2'
params.outdir_stage_stage_parent = '/home/ilyass09/scratch/riboseq_pipeline'
// Verification functions for each file type
def checkFastqExists(gse, gsm, drug, bio, sp) {
    def outputPath = "${params.input_output_fastq_second_stage}/${gse}_${drug}_${bio}/${gsm}/trimmed"
    if (sp.toLowerCase() == "paired") {
        def file1 = new File("${outputPath}/${gsm}_1_val_1.fq")
        def file2 = new File("${outputPath}/${gsm}_2_val_2.fq")
        return file1.exists() && file2.exists()
    } else {
        def file = new File("${outputPath}/${gsm}_trimmed.fq")
        return file.exists()
    }
}

// Definition of the function to create channels
def create_initial_channels() {
    def csvChannelPerLigne = Channel
        .fromPath(params.input_csv_test)
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

// Helper function for Bowtie FASTA paths
def getFastaPathBowtie(type) {
    def upperType = type.toUpperCase()
    
    if (upperType == 'HS') {
        return params.path_fasta_HS_B
    } else if (upperType == 'CE') {
        return params.path_fasta_CE_B
    } else if (upperType == 'DM') {
        return params.path_fasta_DM_B
    } else if (upperType == 'SC') {
        return params.path_fasta_SC_B
    } else if (upperType == 'DR') {
        return params.path_fasta_DR_B
    } else if (upperType == 'SM') {
        return params.path_fasta_SM_B
    } else {
        throw new RuntimeException("Unrecognized type: ${type}")
    }
}

// Helper function for STAR FASTA paths
def getFastaPathStar(type) {
    def upperType = type.toUpperCase()
    
    if (upperType == 'HS') {
        return params.path_fasta_HS
    } else if (upperType == 'CE') {
        return params.path_fasta_CE
    } else if (upperType == 'DM') {
        return params.path_fasta_DM
    } else if (upperType == 'SC') {
        return params.path_fasta_SC
    } else if (upperType == 'DR') {
        return params.path_fasta_DR
    } else if (upperType == 'SM') {
        return params.path_fasta_SM
    } else {
        throw new RuntimeException("Unrecognized type: ${type}")
    }
}

// Helper function for GTF paths
def getGtfPath(type) {
    def upperType = type.toUpperCase()
    
    if (upperType == 'HS') {
        return params.path_fasta_HS_GTF
    } else if (upperType == 'CE') {
        return params.path_fasta_CE_GTF
    } else if (upperType == 'DM') {
        return params.path_fasta_DM_GTF
    } else if (upperType == 'SC') {
        return params.path_fasta_SC_GTF
    } else if (upperType == 'DR') {
        return params.path_fasta_DR_GTF
    } else if (upperType == 'SM') {
        return params.path_fasta_SM_GTF
    } else {
        throw new RuntimeException("Unrecognized type: ${type}")
    }
}

// Add these functions with other verification functions
def checkFastaExists(path) {
    if (!path) {
        println "[WARNING] FASTA file path is not defined"
        return false
    }
    def file = new File(path)
    if (!file.exists()) {
        println "[ERROR] Missing FASTA file: ${path}"
        return false
    }
    println "[INFO] FASTA file found: ${path}"
    return true
}

def checkGtfExists(path) {
    if (!path) {
        println "[WARNING] GTF file path is not defined"
        return false
    }
    def file = new File(path)
    if (!file.exists()) {
        println "[ERROR] Missing GTF file: ${path}"
        return false
    }
    println "[INFO] GTF file found: ${path}"
    return true
}

process BOWTIE_INDEX {
    cpus 1
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
    echo "Bowtie indexing for ${type} with file ${fasta_file}"
    bowtie2-build ${fasta_file} ${type}
    
    # Verify generated indices
    echo "Generated Bowtie indices:"
    ls -la *.bt2
    
    echo "Bowtie indexing completed with prefix: ${type}"
    """
}

process STAR_INDEX {
    cpus 5
    memory '40 GB'
    time '12h'
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
    echo "STAR indexing for ${type} with file ${fasta_file} and GTF ${gtf}"
    mkdir -p ${type}_star_index
    STAR --runThreadN 5\
     --genomeDir ${type}_star_index \
     --runMode genomeGenerate \
     --genomeFastaFiles ${fasta_file} \
     --sjdbGTFfile ${gtf} \
     --genomeSAindexNbases 8
    echo "STAR indexing completed"
    """
}

// Definition of the process bowtie_single
process BOWTIE_SINGLE {
    cpus 5
    memory '25 GB'
    time '12h'
    beforeScript 'module load bowtie2'  
    // Maximum attempt of 3    
    // Dynamic strategy: retry up to the 3rd attempt, then ignore
    errorStrategy 'ignore'     
    tag "${gsm}"
    publishDir "${params.outdir_stage_stage}/bowtie/${gse}_${drug}_${bio}/${gsm}", mode: 'copy'
       
    input:
    tuple val(type), val(gse), val(gsm), val(drug), val(bio), val(sp), path(fastq_file), path(bowtie_indexes)

    output:
    tuple val(type), val(gse), val(gsm), val(drug), val(bio), path("*_unmapped.fq"), emit: unmapped
    path("${gsm}_bowtie.log"), emit: log_file

    script:
    // Retrieve the index prefix (without extension)
    def index_prefix = bowtie_indexes[0].toString() - ~/\.[0-9]\.bt2$/
    
    """
    # Ensure the FASTQ file is correctly identified
    FASTQ=\$(ls -1 ${fastq_file})
    
    echo "Using index: ${index_prefix}"
    echo "Input file: \$FASTQ"
    echo "Index files: ${bowtie_indexes.join(', ')}"
    
    # Run Bowtie2 with parameters for single-end and redirect to /dev/null
    # since we are only interested in unmapped reads
    bowtie2 \
        -p 5 \
        -x ${index_prefix} \
        -U \$FASTQ \
        --un ${gsm}_unmapped.fq \
        2> ${gsm}_bowtie.log
    """
}

process BOWTIE_PAIRED {
    cpus 10
    memory '50 GB'
    time '3h'
    stageInMode 'copy'  // Force copying files instead of symbolic links
    beforeScript 'module load bowtie2'  
    // Maximum attempt of 3
    maxRetries 3
    
    // Dynamic strategy: retry up to the 3rd attempt, then ignore
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }    
    tag "${gsm}"
    publishDir "${params.outdir_stage_stage}/bowtie/${gse}_${drug}_${bio}/${gsm}", mode: 'copy'
       
    input:
    tuple val(type), val(gse), val(gsm), val(drug), val(bio), val(sp), path(fastq_file_1), path(fastq_file_2), path(bowtie_indexes)

    output:
    tuple val(type), val(gse), val(gsm), val(drug), val(bio), path("*_unmapped_*.fq"), emit: unmapped
    path("${gsm}_bowtie.log"), emit: log_file

    script:
    // Retrieve the index prefix (without extension)
    def index_prefix = bowtie_indexes[0].toString() - ~/\.[0-9]\.bt2$/

    """
    bowtie2 \
        -p 10 \
        -x ${index_prefix} \
        -1 ${fastq_file_1} \
        -2 ${fastq_file_2} \
        --un-conc ${gsm}_unmapped_%.fq \
        2> ${gsm}_bowtie.log > /dev/null
    """
}


// Update of the process STAR_SINGLE
process STAR_SINGLE {
    cpus 5
    memory '60 GB'
    time '12h'
    beforeScript 'module load star'  
    // Maximum attempt of 3
    errorStrategy 'ignore'     
 
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
    echo "STAR processing for ${gsm} with ${unmapped_fq}"

     # Determine the parent directory containing the index files
    INDEX_DIR=\$(dirname \$(readlink -f ${star_indexes[0]}))
    
    STAR --runThreadN 5 \
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
         
    echo "STAR processing completed for ${gsm}"
    echo "Verification of generated files:"
    ls -la
    """
}

// Update of the process STAR_PAIRED
process STAR_PAIRED {
    cpus 5
    memory '70 GB'
    time '6h' 
    beforeScript 'module load star'  
    // Maximum attempt of 3
    maxRetries 3
    
    // Dynamic strategy: retry up to the 3rd attempt, then ignore
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
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
    // Detect the pattern of unmapped files (may vary depending on Bowtie2 version)
    def read1 = unmapped_fq.find { it.toString().contains("_1.fq") || it.toString().contains("_unmapped_1.fq") }
    def read2 = unmapped_fq.find { it.toString().contains("_2.fq") || it.toString().contains("_unmapped_2.fq") }
    
    if (!read1 || !read2) {
        // Fallback if the pattern is different
        def files = unmapped_fq.sort()
        read1 = files[0]
        read2 = files[1]
    }
    
    """
    echo "STAR processing for ${gsm} with ${read1} and ${read2}"
     
     # Determine the parent directory containing the index files
    INDEX_DIR=\$(dirname \$(readlink -f ${star_indexes[0]}))
    
    STAR --runThreadN 5 \
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
         
    echo "STAR processing completed for ${gsm}"
    echo "Verification of generated files:"
    ls -la
    """
}

// Main workflow

workflow {
    // Create the initial channel
    def csvChannelPerGsm = create_initial_channels()

    // Create a channel for existing files
    def existing_fastq_files = csvChannelPerGsm
        .filter { it -> 
            checkFastqExists(it[0], it[1], it[2], it[3], it[5])
        }
        .map { gse, gsm, drug, bio, trim, sp, type ->
            def outputPath = "${params.input_output_fastq_second_stage}/${gse}_${drug}_${bio}/${gsm}/trimmed"
            if (sp.toLowerCase() == "paired") {
                def fastq1 = file("${outputPath}/${gsm}_1_val_1.fq") // to be modified to /${gsm}_1_trimmed.fq so the first stage needs to be adapted
                def fastq2 = file("${outputPath}/${gsm}_2_val_2.fq")
                tuple(type, gse, gsm, drug, bio, sp, fastq1, fastq2)
            } else {
                def fastq = file("${outputPath}/${gsm}_trimmed.fq")
                tuple(type, gse, gsm, drug, bio, sp, fastq)
            }
        }

    // Split into two distinct channels based on sequencing type
    def (singleChannel, pairedChannel) = existing_fastq_files.branch {
        single: it[5].toLowerCase() == 'single'
        paired: it[5].toLowerCase() == 'paired'
    }

    Channel.of(params.input_csv_test).view{ "CSV path: $it" }
    Channel.of(params.input_output_fastq_second_stage).view{ "FASTQ files path: $it" }

    // Debug input channels
    singleChannel.count().subscribe { count ->
        println "\n[INFO] Number of single-end samples to process: ${count}\n"
    }
    
    pairedChannel.count().subscribe { count ->
        println "\n[INFO] Number of paired-end samples to process: ${count}\n"
    }
    
    // Extract unique types and verify if necessary files are present as paths
    def species_types = existing_fastq_files
        .map { it -> it[0] }  // Extract the type (species)
        .unique()
        .filter { type -> 
            def fasta_bowtie = getFastaPathBowtie(type)
            def fasta_star = getFastaPathStar(type)
            def gtf_star = getGtfPath(type)
            
            def filesExist = checkFastaExists(fasta_bowtie) && 
                          checkFastaExists(fasta_star) && 
                          checkGtfExists(gtf_star)
            
            if (!filesExist) {
                println "[WARNING] Missing files for type: ${type}"
            }
            return filesExist
        }

    // Create channels for Bowtie and STAR separately for indexing
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

    // Launch indexing processes 
    def bowtie_indices = BOWTIE_INDEX(bowtie_input)
    def star_indices = STAR_INDEX(star_input)

    // Debug index information
    bowtie_indices.bowtie_index_output.view { type, files ->
        println "[INFO] Bowtie index generated for: $type (${files.size()} files)"
    }

    star_indices.star_index_output.view { type, files ->
        println "[INFO] STAR index generated for: $type (${files.size()} files)"
    }

    // Transform the index channel to have the correct format for combination
    def bowtie_indices_formatted = bowtie_indices.bowtie_index_output

    // SingleChannel contains the CSV lines that are single    
    // Combine channels with corresponding indices - FIXED
    def bowtie_single_with_index = singleChannel
        .map { type, gse, gsm, drug, bio, sp, fastq -> 
            return tuple(type, gse, gsm, drug, bio, sp, fastq)
        }
        .combine(bowtie_indices_formatted, by: 0) // The combine command links between the CSV file lines and the files produced by indexing 
    
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


// // Mise à jour de la section onComplete
// workflow.onComplete {
//     // Comptabiliser les résultats
//     def bowtie_count = 0
//     def star_count = 0
//     def transcriptome_count = 0
    
//     // Vérifier les fichiers de sortie Bowtie
//     def bowtie_dir = new File("${params.outdir_stage_stage}/bowtie")
//     if (bowtie_dir.exists()) {
//         bowtie_count = bowtie_dir.listFiles().findAll { it.isDirectory() }.size()
//     }
    
//     // Vérifier les fichiers de sortie STAR
//     def star_dir = new File("${params.outdir_stage_stage}/STAR")
//     if (star_dir.exists()) {
//         star_count = star_dir.listFiles().findAll { it.isDirectory() }.size()
        
//         // Compter les fichiers transcriptome
//         transcriptome_count = 0
//         star_dir.eachFileRecurse { file ->
//             if (file.name.endsWith("Aligned.toTranscriptome.out.bam")) {
//                 transcriptome_count++
//             }
//         }
//     }
    
//     log.info """
//     Pipeline execution summary
//     -------------------------
//     Completed at: ${workflow.complete}
//     Duration    : ${workflow.duration}
//     Success     : ${workflow.success}
//     workDir     : ${workflow.workDir}
//     exit status : ${workflow.exitStatus}
    
//     Results:
//     --------
//     Bowtie alignment completed: $bowtie_count échantillons
//     STAR alignment completed: $star_count échantillons
//     Transcriptome alignments: $transcriptome_count fichiers
    
//     Output directory: ${params.outdir_stage_stage}
//     """
// }
