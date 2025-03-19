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

## Format du fichier CSV d'entrée

Le pipeline utilise un fichier CSV structuré contenant les métadonnées des échantillons à analyser. Voici les colonnes requises :

| Colonne | Description |
|---------|-------------|
| Name | Nom de l'étude ou de l'auteur |
| Species | Espèce biologique (ex: HS pour Homo sapiens) |
| Study_accession | Identifiant GSE de l'étude |
| Project_link | Lien vers la page GEO de l'étude |
| PMID | Identifiant PubMed de la publication associée |
| Treatment_type | Type de traitement (pre-lysis, post-lysis) |
| Drug | Médicament utilisé (ex: cycloheximide) |
| Sample_accession | Identifiants GSM des échantillons, séparés par des points-virgules |
| Biological_type | Type biologique de l'échantillon (ex: Cells_A549, HAEC) |
| Ribo_type | Type de ribosome (monosome, polysome) |
| Trim_arg | Arguments additionnels pour Trim Galore (optionnel) |
| S_P_type | Type de séquençage (single ou paired) |

## Installation et configuration

1. Clonez ce dépôt :
   ```bash
   git clone https://github.com/ilokinn/Nextflow.git
   cd Nextflow
   ```

2. Avant d'exécuter le pipeline, vous devez :
   - Faire un `git pull` pour obtenir les dernières mises à jour
   - Modifier les chemins dans le fichier de configuration (`nextflow.config`) pour les adapter à votre environnement

## Chemins à modifier dans le fichier de configuration

Vous devez modifier les chemins suivants dans le fichier `nextflow.config` pour les adapter à votre environnement :

1. **Chemins des fichiers de référence** :
   - `params.path_fasta_HS_B` : Fichier FASTA pour filtrage Bowtie (ARNs non-codants humains)
   - `params.path_fasta_HS` : Génome de référence humain (FASTA)
   - `params.path_fasta_HS_GTF` : Annotation du génome humain (GTF)

2. **Chemins des répertoires de sortie** :
   - `params.outdir_first_stage` : Répertoire de sortie pour la première étape
   - `params.input_output_fastq_second_stage` : Répertoire contenant les fichiers FASTQ traités
   - `params.outdir_stage_stage` : Répertoire de sortie pour la deuxième étape
   - `params.outdir_stage_stage_parent` : Répertoire parent pour les index
   - `params.outdir_ribotie` : Répertoire de sortie pour RiboTIE

3. **Chemins des fichiers d'entrée** :
   - `params.input_csv` : Fichier CSV contenant les métadonnées des échantillons
   - `params.multiqc_config` : Fichier de configuration pour MultiQC
   - `params.star_dir` : Répertoire contenant les résultats STAR
   - `params.ribotie_dir` : Répertoire d'installation de RiboTIE

4. **Options de cluster** :
   - `process.clusterOptions` : Modifier l'option `--account=rrg-xroucou` selon votre compte sur le cluster

## Utilisation

Pour exécuter le pipeline complet, suivez ces étapes :

1. Exécution de la première étape (téléchargement et prétraitement) :
   ```bash
   nextflow run First_stage_pipeline.nf --input_csv samples.csv
   ```

2. Génération des rapports qualité MultiQC (en utilisant le fichier de config) :
   ```bash
   nextflow run MultiQc.nf -c my_config.yaml
   ```

3. (Optionnel) Centralisation des rapports MultiQC :
   ```bash
   nextflow run MultiQc_Rapport_Management.nf
   ```

4. Exécution de la deuxième étape (alignement) :
   ```bash
   nextflow run Second_Stage.3.1.nf -profile beluga
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