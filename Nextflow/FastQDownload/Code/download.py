# Description: Script Python pour télécharger des données SRA à partir de GEO.
# Auteur: Ilyass Jenjare
# Date: 2025-01-16
# Résumé: Ce script permet de télécharger des données SRA à partir de GEO en utilisant des identifiants GSE et GSM.
# Entrées: Fichier CSV contenant les identifiants GSE et GSM.
# Sorties: Fichiers FASTQ téléchargés pour chaque GSM.
# les coloonnes du fichier csv sont : Study_accession ; Drug ; Biological_type ; Sample_accession ...


import pandas as pd
import os
import shutil
from pysradb import SRAweb
import argparse
import subprocess
import sys
import time
import code


# --- Étape 1 : Importer les données ---
def import_data(file_path):
    return pd.read_csv(file_path)


# --- Étape 2 : Traiter les données ---
def extract_gsms_until_nan(df):
    gsm_lists = []
    for index, row in df.iterrows():
        if pd.isna(row['Sample_accession']):
            print(f"Arrêt à l'index {index} en raison d'une valeur NaN.")
            break
        gsm_list = row['Sample_accession'].split("; ")
        gsm_lists.append(gsm_list)
    return gsm_lists


def process_data(data_frame):
    gse_list = data_frame['Study_accession'].tolist()
    drug_list = data_frame['Drug'].tolist()
    biological_type_list = data_frame['Biological_type'].tolist()
    gsm_list = extract_gsms_until_nan(data_frame)
    return gse_list, drug_list, biological_type_list, gsm_list


# --- Étape 3 : Conversion GSM en SRR ---
def convert_gsms_to_srrs_with_db(gsm_lists):
    """
    Convertit une liste de listes de GSMs en une liste de listes de SRRs en utilisant pysradb.

    :param gsm_lists: Liste de listes des GSMs.
    :return: Liste de listes des SRRs correspondant aux GSMs.
    """
    db = SRAweb()  # Initialiser l'instance SRAweb
    srr_lists = []  # Liste principale pour stocker les SRRs

    # Parcourir chaque liste de GSMs
    for gsm_list in gsm_lists:
        srr_list = []  # Liste pour stocker les SRRs pour une liste de GSMs
        for gsm in gsm_list:
            try:
                # Utiliser pysradb pour convertir GSM en SRR
                srr_df = db.gsm_to_srr(gsm)
                print(f"GSM {gsm} converti en SRRs:")
                print(srr_df)

                if not srr_df.empty:
                    # Récupérer les SRRs (colonne 'Run_accession') sous forme de liste
                    srrs = srr_df['run_accession'].tolist()
                    srr_list.extend(srrs)  # Ajouter tous les SRRs pour le GSM
                else:
                    print(f"Aucun SRR trouvé pour le GSM {gsm}")
                    srr_list.append(None)  # Aucun SRR trouvé
            except Exception as e:
                print(f"Erreur avec le GSM {gsm}: {e}")
                srr_list.append(None)  # Erreur lors de la conversion
            time.sleep(1)

        srr_lists.append(srr_list)

    return srr_lists


def download_fastq(acc, gse, btype, drug):
    """
    Télécharge directement les fichiers FASTQ à partir de l'accession SRA en utilisant fasterq-dump.
    """
    print(f"Processing accession: {acc}...")

    # Définir le répertoire de sortie
    output_dir = f'output/{gse}_{btype}_{drug}'
    os.makedirs(output_dir, exist_ok=True)

    # Vérifier si les fichiers sont déjà téléchargés
    fastq_files = [
        os.path.join(output_dir, f"{acc}_1.fastqsanger"),
        os.path.join(output_dir, f"{acc}_2.fastqsanger"),
        os.path.join(output_dir, f"{acc}_single.fastqsanger"),
    ]
    if any(os.path.exists(f) for f in fastq_files):
        print(f"Accession {acc} already processed. Skipping download.")
        return

    # Télécharger directement avec fasterq-dump
    fasterq_command = [
        "fasterq-dump", "--split-files", "--outdir", output_dir, acc
    ]
    try:
        subprocess.run(fasterq_command, check=True)  # Exécute la commande
    except subprocess.CalledProcessError as e:
        print(f"Error during fasterq-dump for {acc}: {e}")
        return

    # Vérifier les fichiers générés
    fastq_files_generated = [
        os.path.join(output_dir, f) for f in os.listdir(output_dir)
        if f.startswith(acc) and f.endswith('.fastq')
    ]
    if len(fastq_files_generated) == 2:
        print(f"Fichiers générés : {fastq_files_generated[0]} (forward), {fastq_files_generated[1]} (reverse)")
    elif len(fastq_files_generated) == 1:
        print(f"Fichier généré : {fastq_files_generated[0]} (single)")
    else:
        print(f"Aucun fichier FASTQ généré pour {acc}.")
        return

    print(f"Completed processing for {acc}")


# --- Étape 5 : Processer les échantillons ---
def process_samples(combined, loaded_list):
    for (gse, cell_type, drug), srr_list in zip(combined, loaded_list):
        print(f"Processing GSE: {gse}, Cell Type: {cell_type}, Drug: {drug}")
        for srr in srr_list:
            if srr:
                download_fastq(srr, gse, cell_type, drug)


# --- Point d'entrée principal ---
if __name__ == "__main__":
    # Charger les données
    df = import_data("/home/machinegun/copy_df.csv")
    gse_list, drug_list, bio_type_list, gsm_list = process_data(df)

    print("GSE:", gse_list)
    print("Drug:", drug_list)
    print("Biological Types:", bio_type_list)
    print("GSMs:", gsm_list)

    # Créer la liste combinée
    combined = [[gse, btype, drug] for gse, btype, drug in zip(gse_list, bio_type_list, drug_list)]

    print("Combined:", combined)

    # Conversion GSM -> SRR
    # list_sra = [['SRR31074191', 'SRR31074192', 'SRR31074193', 'SRR31074194', 'SRR31074195', 'SRR31074196', 'SRR31074197', 'SRR31074198', 'SRR31074199', 'SRR31074200', 'SRR31074201', 'SRR31074202', 'SRR31074203', 'SRR31074204', 'SRR31074205', 'SRR31074206', 'SRR31074207', 'SRR31074208', 'SRR31074209', 'SRR31074210'], ['SRR27352753', 'SRR27352754', 'SRR27352745', 'SRR27352746', 'SRR27352747', 'SRR27352748', 'SRR27352749', 'SRR27352750', 'SRR27352751', 'SRR27352752', 'SRR27352755', 'SRR27352756', 'SRR27352757', 'SRR27352758', 'SRR27352759', 'SRR27352760', 'SRR27352761', 'SRR27352762'], ['SRR12318118', 'SRR12318117', 'SRR12318116', 'SRR12318115', 'SRR12318114', 'SRR12318113', 'SRR12318112', 'SRR12318111', 'SRR12318110', 'SRR12318109', 'SRR12318108', 'SRR12318107', 'SRR12318106', 'SRR12318105'], ['SRR14189736', 'SRR14189735', 'SRR14189734', 'SRR14189733', 'SRR14189732', 'SRR14189731', 'SRR14189730', 'SRR14189729', 'SRR14189728', 'SRR14189727', 'SRR14189726', 'SRR14189725'], ['SRR26543959', 'SRR26543960', 'SRR26543961', 'SRR26543962', 'SRR26543963', 'SRR26543964', 'SRR26543965', 'SRR26543966', 'SRR26543967', 'SRR26543968', 'SRR26543969', 'SRR26543970', 'SRR26543971', 'SRR26543972', 'SRR26543973', 'SRR26543974', 'SRR26543975', 'SRR26543976', 'SRR26543977', 'SRR26543978', 'SRR26543979', 'SRR26543980', 'SRR26543981', 'SRR26543982', 'SRR26543983', 'SRR26543984', 'SRR26543985', 'SRR26543986'], ['SRR29790646', 'SRR29790647', 'SRR29790648', 'SRR29790649', 'SRR29790650', 'SRR29790651', 'SRR29790658', 'SRR29790659', 'SRR29790660', 'SRR29790661', 'SRR29790662', 'SRR29790663'], ['SRR30163335', 'SRR30163336', 'SRR30163337', 'SRR30163338', 'SRR30163339', 'SRR30163340'], ['SRR27534371', 'SRR27534372', 'SRR27534373', 'SRR27534374', 'SRR27534375', 'SRR27534376'], ['SRR23190729', 'SRR23190730', 'SRR23190731', 'SRR23190732', 'SRR23190733', 'SRR23190734', 'SRR23190735', 'SRR23190736', 'SRR23190737', 'SRR23190738', 'SRR23190739', 'SRR23190740', 'SRR23190741', 'SRR23190742', 'SRR23190743', 'SRR23190744'], ['SRR17429210', 'SRR17429211', 'SRR17429212', 'SRR17429213', 'SRR17429214', 'SRR17429215', 'SRR17429216', 'SRR17429217', 'SRR17429218', 'SRR17429219', 'SRR17429220', 'SRR17429221', 'SRR17429222', 'SRR17429223', 'SRR17429224', 'SRR17429225', 'SRR17429226', 'SRR17429227'], ['SRR15513226', 'SRR15513225', 'SRR15513224', 'SRR15513223', 'SRR15513222', 'SRR15513221', 'SRR15513220', 'SRR15513219', 'SRR15513218', 'SRR15513217', 'SRR15513216'], ['SRR15513215', 'SRR15513214', 'SRR15513213'], ['SRR15513212', 'SRR15513211', 'SRR15513210', 'SRR15513209', 'SRR15513208'], ['SRR15513207', 'SRR15513206', 'SRR15513205', 'SRR15513204', 'SRR15513203'], ['SRR15513202', 'SRR15513201', 'SRR15513200', 'SRR15513199', 'SRR15513198', 'SRR15513197'], ['SRR15513226', 'SRR15513225', 'SRR15513224', 'SRR15513223', 'SRR15513222', 'SRR15513221', 'SRR15513220', 'SRR15513219', 'SRR15513218', 'SRR15513217', 'SRR15513216'], ['SRR15513196', 'SRR15513195', 'SRR15513194', 'SRR15513193', 'SRR15513192', 'SRR15513191', 'SRR15513190', 'SRR15513189', 'SRR15513188', 'SRR15513187', 'SRR15513186', 'SRR15513185', 'SRR15513184', 'SRR15513183', 'SRR15513182', 'SRR15513181', 'SRR15513180', 'SRR15513179', 'SRR15513178', 'SRR15513177', 'SRR15513176', 'SRR15513175', 'SRR15513174', 'SRR15513173', 'SRR15513172', 'SRR15513171', 'SRR15513170', 'SRR15513169', 'SRR15513168', 'SRR15513167', 'SRR15513166', 'SRR15513165'], ['SRR15513164', 'SRR15513163', 'SRR15513162', 'SRR15513161', 'SRR15513160', 'SRR15513159'], ['SRR15513158', 'SRR15513157', 'SRR15513156', 'SRR15513155', 'SRR15513154', 'SRR15513153'], ['SRR15513152', 'SRR15513151', 'SRR15513150', 'SRR15513149', 'SRR15513148'], ['SRR13664946', 'SRR13664945', 'SRR13664944', 'SRR13664943', 'SRR13664942'], ['SRR11432008', 'SRR11431992', 'SRR11431993', 'SRR11431998', 'SRR11431997', 'SRR11431996', 'SRR11431995', 'SRR11431994'], ['SRR19387516'], ['SRR29874569', 'SRR29874567', 'SRR29874568', 'SRR29874570', 'SRR29874571'], ['SRR28214280', 'SRR28214281', 'SRR28214282', 'SRR28214283', 'SRR28214284', 'SRR28214285', 'SRR28214286', 'SRR28214287', 'SRR28214288', 'SRR28214289', 'SRR28214290', 'SRR28214291', 'SRR28214292', 'SRR28214293', 'SRR28214294', 'SRR28214295', 'SRR28214296', 'SRR28214297', 'SRR28214298', 'SRR28214299', 'SRR28214300', 'SRR28214301', 'SRR28214302', 'SRR28214303', 'SRR28214280', 'SRR28214281', 'SRR28214282', 'SRR28214283', 'SRR28214284', 'SRR28214285', 'SRR28214286', 'SRR28214287', 'SRR28214288', 'SRR28214289', 'SRR28214290', 'SRR28214291'], ['SRR25764045', 'SRR25764046', 'SRR25764045', 'SRR25764046'], ['SRR23760656', 'SRR23760657', 'SRR23760658', 'SRR23760659', 'SRR23760656', 'SRR23760657', 'SRR23760658', 'SRR23760659'], ['SRR27731317', 'SRR27731318', 'SRR27731319', 'SRR27731320', 'SRR27731321', 'SRR27731322'], ['SRR25421838', 'SRR25421839', 'SRR25421840', 'SRR25421841', 'SRR25421842', 'SRR25421843', 'SRR25421844', 'SRR25421845', 'SRR25421846', 'SRR25421847', 'SRR25421848', 'SRR25421849', 'SRR25421850', 'SRR25421851', 'SRR25421852', 'SRR25421853', 'SRR25421854', 'SRR25421855', 'SRR25421847', 'SRR25421848', 'SRR25421849', 'SRR25421850', 'SRR25421851', 'SRR25421852', 'SRR25421853', 'SRR25421854', 'SRR25421855'], ['SRR25764046', 'SRR25764045'], ['SRR25764045', 'SRR25764046']]

    list_sra = convert_gsms_to_srrs_with_db(gsm_list)
    print("Converted SRRs:", list_sra)

    process_samples(combined,list_sra)




