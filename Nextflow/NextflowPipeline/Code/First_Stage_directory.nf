// Paramètres
params.input_dir = null
params.outdir = "${baseDir}/Files"

process fastqc_1 {
    publishDir "${params.outdir}/${fastq_file.simpleName}/Premier_Fastqc", mode: 'copy'
    memory '4 GB'


    input:
    tuple val(gse_dir), path(fastq_file)
    
    output:
    path "*_fastqc.{zip,html}", emit: reports
    
    script:
    """
    fastqc ${fastq_file}
    """
}

process first_fastp {
    publishDir "${params.outdir}/${reads.simpleName}/Premier_Fastp", mode: 'copy'
    memory '4 GB'


    input:
    tuple val(gse_dir), path(reads)
    
    output:
    tuple val(gse_dir), path("trimmed_*.fastq.gz"), emit: trimmed_reads
    path "fastp.json"
    path "fastp.html"
    
    script:
    """
    fastp -i ${reads} \
          -o trimmed_${reads.simpleName}.fastq.gz \
          -h fastp.html \
          -j fastp.json
    """
}

process sec_fastp {
    publishDir "${params.outdir}/${reads.baseName.replaceFirst(/^trimmed_/, '').replaceFirst(/\.fastq$/, '')}/Deuxieme_Fastp", mode: 'copy'
    memory '4 GB'


    input:
    tuple val(gse_dir), path(reads)
    
    output:
    tuple val(gse_dir), path("trimmed_*.fastq.gz"), emit: trimmed_reads
    path "fastp.json"
    path "fastp.html"
    
    script:
    """
    fastp -i ${reads} \
          -o trimmed_${reads.simpleName}.fastq.gz \
          -h fastp.html \
          -j fastp.json \
          -b 40 \
          -M 25 \
          -3 \
          -W 2 \
          --length_limit 40 \
          --length_required 20
    """
}

process fastqc_2 {
    publishDir "${params.outdir}/${fastq_file.baseName.replaceFirst(/^trimmed_/, '').replaceFirst(/^trimmed_/, '').replaceFirst(/\.fastq$/, '').replaceFirst(/\.fastq$/, '')}/Dernier_Fastqc", mode: 'copy'
    memory '4 GB'

    
    input:
    tuple val(gse_dir), path(fastq_file)
    
    output:
    path "*_fastqc.{zip,html}"
    
    script:
    """
    fastqc ${fastq_file}
    """
}

workflow {
    Channel
        .fromPath("${params.input_dir}/GSE*/SRR*.fastq")
        .map { file -> 
            def gse_dir = file.parent.name
            tuple(gse_dir, file)
        }
        .set { input_ch }
    
    fastqc_1(input_ch)
    first_fastp(input_ch)
    sec_fastp(first_fastp.out.trimmed_reads)
    fastqc_2(sec_fastp.out.trimmed_reads)
}