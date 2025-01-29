// Paramètres
params.input = null
params.outdir = "${baseDir}/Files"

process fastqc_1 {
    publishDir "${params.outdir}/${fastq_file.baseName}/Premier_Fastqc", mode: 'copy'
    
    input:
    path fastq_file
    
    output:
    path "*_fastqc.{zip,html}", emit: reports
    
    script:
    """
    fastqc ${fastq_file}
    """
}

process first_fastp {
    publishDir "${params.outdir}/${reads.baseName}/Premier_Fastp", mode: 'copy'
    
    input:
    path reads
    
    output:
    path "trimmed_*.fastq.gz", emit: trimmed_reads
    path "fastp.json"
    path "fastp.html"
    
    script:
    """
    fastp -i ${reads} \
          -o trimmed_${reads.baseName}.fastq.gz \
          -h fastp.html \
          -j fastp.json
    """
}

process sec_fastp {
    publishDir "${params.outdir}/${reads.baseName.replaceFirst(/^trimmed_/, '').replaceFirst(/\.fastq$/, '')}/Deuxieme_Fastp", mode: 'copy'
    
    input:
    path reads
    
    output:
    path "trimmed_*.fastq.gz", emit: trimmed_reads
    path "fastp.json"
    path "fastp.html"
    
    script:
    """
    fastp -i ${reads} \
          -o trimmed_${reads.baseName}.fastq.gz \
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
    
    input:
    path fastq_file
    
    output:
    path "*_fastqc.{zip,html}"
    
    script:
    """
    fastqc ${fastq_file}
    """
}

workflow {
    input_ch = Channel.fromPath(params.input)
    
    fastqc_1(input_ch)
    first_fastp(input_ch)
    sec_fastp(first_fastp.out.trimmed_reads)
    fastqc_2(sec_fastp.out.trimmed_reads)
}