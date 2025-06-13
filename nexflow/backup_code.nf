// Function to check if files exist in the fastq directory
def checkFastqExists(gsm, sp) {
    def outputPath = "${params.outdir_file_exist_already}"
    def found = false
    
    try {
        def dir = new File(outputPath)
        def files = dir.list()
        
        if (sp.toLowerCase() == "paired") {
            def file1Exists = files.any { it == "${gsm}_1.fastq" }
            def file2Exists = files.any { it == "${gsm}_2.fastq" }
            found = file1Exists && file2Exists
        } else {
            found = files.any { it == "${gsm}.fastq" }
        }
    } catch (Exception e) {
        log.error "Error while checking files for GSM: $gsm - ${e.message}"
    }
    
    return found
}

// Verification functions for post-download steps
def checkFastqcPreExists(gse, gsm, drug, bio, sp) {
    def outputPath = "${params.outdir}/${gse}_${drug}_${bio}/${gsm}/fastqc_pre"
    if (sp.toLowerCase() == "paired") {
        return new File("${outputPath}/${gsm}_1_fastqc.html").exists() &&
               new File("${outputPath}/${gsm}_2_fastqc.html").exists() &&
               new File("${outputPath}/${gsm}_1_fastqc.zip").exists() &&
               new File("${outputPath}/${gsm}_2_fastqc.zip").exists()
    } else {
        return new File("${outputPath}/${gsm}_fastqc.html").exists() &&
               new File("${outputPath}/${gsm}_fastqc.zip").exists()
    }
}

def checkTrimExists(gse, gsm, drug, bio, sp) {
    def outputPath = "${params.outdir}/${gse}_${drug}_${bio}/${gsm}/trimmed"
    if (sp.toLowerCase() == "paired") {
        return new File("${outputPath}/${gsm}_1_val_1.fq").exists() &&
               new File("${outputPath}/${gsm}_2_val_2.fq").exists()
    } else {
        return new File("${outputPath}/${gsm}_trimmed.fq").exists()
    }
}

def checkFastqcPostExists(gse, gsm, drug, bio, sp) {
    def outputPath = "${params.outdir}/${gse}_${drug}_${bio}/${gsm}/fastqc_post"
    if (sp.toLowerCase() == "paired") {
        return new File("${outputPath}/${gsm}_1_val_1_fastqc.html").exists() &&
               new File("${outputPath}/${gsm}_2_val_2_fastqc.html").exists() &&
               new File("${outputPath}/${gsm}_1_val_1_fastqc.zip").exists() &&
               new File("${outputPath}/${gsm}_2_val_2_fastqc.zip").exists()
    } else {
        return new File("${outputPath}/${gsm}_trimmed_fastqc.html").exists() &&
               new File("${outputPath}/${gsm}_trimmed_fastqc.zip").exists()
    }
}

workflow {
    log.info "Starting pipeline for existing files only..."
    // Create initial channel
    def csvChannelPerGsm = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> [
            GSE: row.Study_accession,
            GSM: row.Sample_accession.split(';').collect{it.trim()},
            drug: row.Drug,
            riboseq_type: row.Biological_type,
            trimming_args: row.Trim_arg,
            paired_end: row.paired_end.toBoolean(),
            sp: row.Species
        ]}.transpose(by: 1).view()
    
    // Filter only for existing files
    def existing_channel = csvChannelPerGsm.filter {
        v -> checkFastqExists(v.gsm, v.sp)
    }
    
    // // Log existing files
    // existing_channel.view { gse, gsm, drug, bio, trim, sp ->
    //     "Existing file: GSE=$gse, GSM=$gsm - will be processed"
    // }
    
    // Map existing files to include file paths
    def existing_fastq_files = existing_channel
        .map { gse, gsm, drug, bio, trim, sp, type_cells ->
            def fastqFiles
            if (sp.toLowerCase() == "paired") {
                fastqFiles = [
                    file("${params.outdir_file_exist_already}/${gsm}_1.fastq"),
                    file("${params.outdir_file_exist_already}/${gsm}_2.fastq")
                ]
            } else {
                fastqFiles = file("${params.outdir_file_exist_already}/${gsm}.fastq")
            }
            tuple(gse, gsm, drug, bio, trim, sp,type_cells, fastqFiles)
        }
    
    // Execute FASTQC_PRE for existing files
    FASTQC_PRE(
        input: existing_fastq_files
    )


    // Map existing files to include file paths
    def existing_fastq_files = existing_channel
        .map { gse, gsm, drug, bio, trim, sp, type_cells ->
            def fastqFiles
            if (sp.toLowerCase() == "paired") {
                fastqFiles = [
                    file("${params.outdir_file_exist_already}/${gsm}_1.fastq"),
                    file("${params.outdir_file_exist_already}/${gsm}_2.fastq")
                ]
            } else {
                fastqFiles = file("${params.outdir_file_exist_already}/${gsm}.fastq")
            }
            tuple(gse, gsm, drug, bio, trim, sp,type_cells, fastqFiles)
        }