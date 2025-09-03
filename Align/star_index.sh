#!/bin/bash
#SBATCH --account=def-xroucou
#SBATCH --mem=50G
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=nicolas.provencher@usherbrooke.ca
#SBATCH --mail-type=ALL

module load star

STAR --runMode genomeGenerate \
    --runThreadN 8 \
    --genomeDir reference/HS/HS_GRCh38_ensembl_114_star_index \
    --genomeFastaFiles reference/HS/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --genomeSAindexNbases 8 \
    --sjdbGTFfile reference/HS/Homo_sapiens.GRCh38.114.gtf
