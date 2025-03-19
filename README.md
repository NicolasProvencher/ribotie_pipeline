# Pipeline d'analyse pour données de ribosome profiling

Ce dépôt contient un ensemble de workflows Nextflow conçus pour l'analyse complète des données de séquençage de ribosome profiling, de l'acquisition des données brutes jusqu'aux analyses d'expression différentielle.

## Structure du pipeline

Le pipeline se compose de trois étapes principales :

1. **Première étape** : Acquisition et prétraitement des données
   - Téléchargement des fichiers FASTQ à partir de GSM/GSE
   - Contrôle qualité avec FastQC (pré-traitement)
   - Nettoyage des lectures avec Trim Galore
   - Contrôle qualité après nettoyage (post-traitement)

2. **Deuxième étape** : Alignement des séquences
   - Indexation du génome avec Bowtie et STAR
   - Filtrage des ARN avec Bowtie
   - Alignement des lectures non mappées avec STAR
   - Génération de fichiers BAM triés pour l'analyse ultérieure

3. **Analyse RiboTIE** : Identification des sites d'initiation de la traduction
   - Création de fichiers YAML spécifiques à chaque expérience
   - Préparation des données avec `ribotie --data`
   - Analyse complète avec l'algorithme RiboTIE

## Caractéristiques

- Support pour plusieurs espèces (humain, souris, drosophile, C. elegans, etc.)
- Traitement parallélisé pour une exécution efficace
- Vérification automatique des fichiers intermédiaires pour éviter les redondances
- Compatibilité avec les données single-end et paired-end
- Intégration optimisée pour les systèmes HPC avec Slurm

## Utilisation

Le pipeline est configuré pour fonctionner avec un fichier CSV d'entrée contenant les métadonnées des échantillons (accessions GEO, type de traitement, paramètres de nettoyage, etc.).

Pour exécuter le pipeline complet :

```bash
nextflow run main.nf -profile slurm --input_csv samples.csv
```

## Prérequis

- Nextflow
- Bowtie2
- STAR
- FastQC
- Trim Galore
- Python 3.9+ (pour RiboTIE)
- CUDA (pour l'accélération GPU de RiboTIE)