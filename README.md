# Analysis Pipeline for Ribosome Profiling Data

This repository contains a set of Nextflow workflows designed for large scale riboseq reanalysis from the download step to detected TIS identification step. 
Packages used
Download:
   sra-tools
QC and filtering:
   multiqc
   fastqc
   trim galore
   bowtie (rRNA reads filtering)
Alignement:
   STAR
ORF calling
   RiboTIE

## HPC Environment Configuration

This pipeline utilizes a HPC configuration as a way to run the RiboTIE step which requires GPU hardware acceleration

## Pipeline Structure

The pipeline consists of three main stages:

- **First and Second Stage**: Using conda for dependency management
  ```bash
  # Activate conda environment
  # After creating it
  conda activate nextflow_env
  ```

1. **First Stage**: Data Acquisition and Preprocessing
   - Download FASTQ files from GSM/GSE
   - Quality control with FastQC (pre-processing)
   - Read cleaning with Trim Galore
   - Quality control after cleaning (post-processing)

2. **Second Stage**: Sequence Alignment
   - Genome indexing with Bowtie and STAR
   - RNA filtering with Bowtie
   - Unmapped reads alignment with STAR
   - Generation of sorted BAM files for further analysis

3. **RiboTIE Analysis**: Translation Initiation Site Identification
   - Creation of experiment-specific YAML files
   - Data preparation with `ribotie --data`
   - Complete analysis with RiboTIE algorithm

## Input CSV File Format

The pipeline uses a structured CSV file containing sample metadata for analysis. Here are the required columns:

| Column | Description |
|--------|-------------|
| Name | Study or author name |
| Species | Biological species (e.g., HS for Homo sapiens) |
| Study_accession | GSE study identifier |
| Project_link | Link to GEO study page |
| PMID | PubMed ID of associated publication |
| Treatment_type | Treatment type (pre-lysis, post-lysis) |
| Drug | Drug used (e.g., cycloheximide) |
| Sample_accession | GSM sample identifiers, separated by semicolons |
| Biological_type | Sample biological type (e.g., Cells_A549, HAEC) |
| Ribo_type | Ribosome type (monosome, polysome) |
| Trim_arg | Additional arguments for Trim Galore (optional) |
| S_P_type | Sequencing type (single or paired) |

## Installation and Configuration

1. Clone this repository:
   ```bash
   git clone https://github.com/ilokinn/Nextflow.git
   cd Nextflow
   ```

2. Before running the pipeline, you need to:
   - Modify paths in the configuration file (`nextflow.config`) to match your environment

## Paths to Modify in Configuration File

You need to modify the following paths in the `nextflow.config` file to match your environment:

1. **Reference File Paths**:
   - `params.path_fasta_HS_B`: FASTA file for Bowtie filtering (human non-coding RNAs)
   - `params.path_fasta_HS`: Human reference genome (FASTA)
   - `params.path_fasta_HS_GTF`: Human genome annotation (GTF)

2. **Output Directory Paths**:
   - `params.outdir_first_stage`: Output directory for first stage
   - `params.input_output_fastq_second_stage`: Directory containing processed FASTQ files
   - `params.outdir_stage_stage`: Output directory for second stage
   - `params.outdir_stage_stage_parent`: Parent directory for indexes
   - `params.outdir_ribotie`: Output directory for RiboTIE

3. **Input File Paths**:
   - `params.input_csv`: CSV file containing sample metadata
   - `params.multiqc_config`: MultiQC configuration file
   - `params.star_dir`: Directory containing STAR results
   - `params.ribotie_dir`: RiboTIE installation directory

4. **Cluster Options**:
   - `process.clusterOptions`: Modify `--account=rrg-xroucou` according to your cluster account

## Usage

To run the complete pipeline, follow these steps:

1. Run first stage (download and preprocessing):
   ```bash
   nextflow run First_stage_pipeline.nf --input_csv samples.csv
   ```

2. Generate MultiQC quality reports (using config file):
   ```bash
   nextflow run MultiQc.nf -c my_config.yaml
   ```

3. (Optional) Centralize MultiQC reports:
   ```bash
   nextflow run MultiQc_Rapport_Management.nf
   ```

4. Run second stage (alignment):
   ```bash
   nextflow run Second_Stage.3.1.nf -profile beluga
   ```

5. Run RiboTIE analysis:
   ```bash
   nextflow run Third_Stage.2.0.nf
   ```

## Prerequisites

- Nextflow
- Bowtie2
- STAR
- FastQC
- Trim Galore
- MultiQC
- Python 3.9+ (for RiboTIE)
- CUDA (for RiboTIE GPU acceleration)

## HPC Environment Configuration

The pipeline is optimized to run on HPC clusters using Slurm as a resource manager. Pre-configured profiles are available in `nextflow.config`.
