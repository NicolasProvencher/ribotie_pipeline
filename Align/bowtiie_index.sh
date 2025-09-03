#!/bin/bash
#SBATCH --account=def-xroucou
#SBATCH --mem=30G
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=nicolas.provencher@usherbrooke.ca
#SBATCH --mail-type=ALL

module load bowtie2

bowtie2-build reference/HS/snRNA_snoRNA_rRNA_HS_ensembl113.fasta reference/HS/snRNA_snoRNA_rRNA_HS_ensembl113_index
