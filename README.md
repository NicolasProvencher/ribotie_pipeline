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

## Installation et configuration

1. Clonez ce dépôt :
   ```bash
   git clone https://github.com/ilokinn/Nextflow.git
   cd Nextflow
   ```

2. Avant d'exécuter le pipeline, vous devez :
   - Faire un `git pull` pour obtenir les dernières mises à jour
   - Modifier les chemins dans le fichier de configuration (`nextflow.config`) pour les adapter à votre environnement
   - Configurer les chemins vers vos fichiers de référence (FASTA, GTF) et répertoires de sortie

## Structure des fichiers

Le pipeline est organisé en plusieurs scripts Nextflow :
- `First_stage_pipeline.nf` : Acquisition et prétraitement des données
- `MultiQc.nf` : Génération des rapports qualité par GSE
- `MultiQc_Rapport_Management.nf` : Centralisation des rapports MultiQC (optionnel)
- `Second_Stage.3.1.nf` : Alignement et filtrage des séquences
- `Third_Stage.2.0.nf` : Analyse RiboTIE pour l'identification des sites d'initiation

## Utilisation

Le pipeline est configuré pour fonctionner avec un fichier CSV d'entrée contenant les métadonnées des échantillons (accessions GEO, type de traitement, paramètres de nettoyage, etc.).

Pour exécuter le pipeline complet, suivez ces étapes :

1. Exécution de la première étape (téléchargement et prétraitement) :
   ```bash
   nextflow run First_stage_pipeline.nf --input_csv samples.csv
   ```

2. Génération des rapports qualité MultiQC :
   ```bash
   nextflow run MultiQc.nf
   ```

3. (Optionnel) Centralisation des rapports MultiQC :
   ```bash
   nextflow run MultiQc_Rapport_Management.nf
   ```

4. Exécution de la deuxième étape (alignement) :
   ```bash
   nextflow run Second_Stage.3.1.nf -profile slurm
   ```

5. Exécution de l'analyse RiboTIE :
   ```bash
   nextflow run Third_Stage.2.0.nf
   ```

## Prérequis

- Nextflow
- Bowtie2
- STAR
- FastQC
- Trim Galore
- MultiQC
- Python 3.9+ (pour RiboTIE)
- CUDA (pour l'accélération GPU de RiboTIE)

## Configuration pour environnements HPC

Le pipeline est optimisé pour fonctionner sur des clusters HPC utilisant Slurm comme gestionnaire de ressources. Des profils préconfigurés sont disponibles dans `nextflow.config`.