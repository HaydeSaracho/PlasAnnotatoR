"""
Script para generar lista de genomas a descargar de RefSeq
1 genoma por especie presente en PLSDB
"""
import pandas as pd
from pathlib import Path

# Rutas
BASE = Path("/home/bionfo/Escritorio/Hayde/PlasAnnotatoR/data")
TAXONOMY = BASE / "plsdb/meta/taxonomy.csv"
ASSEMBLY = BASE / "assembly_summary_refseq.txt"
OUTPUT_TSV = BASE / "genomes_to_download.tsv"
OUTPUT_SH = BASE / "download_chromosomes.sh"
OUTPUT_URLS = BASE / "urls.txt"
CHROMOSOMES_DIR = BASE / "chromosomes"

# Cargar especies de PLSDB
print("Cargando taxonomy de PLSDB...")
tax = pd.read_csv(TAXONOMY)
species = set(tax['TAXONOMY_species'].dropna().str.replace('_', ' '))
print("Total especies en PLSDB: {}".format(len(species)))

# Cargar assembly summary
print("Cargando assembly summary de RefSeq...")
asm = pd.read_csv(ASSEMBLY, sep='\t', skiprows=1, low_memory=False)

# Filtrar Complete Genomes de especies en PLSDB
complete = asm[asm['assembly_level'] == 'Complete Genome']
match = complete[complete['organism_name'].isin(species)]

# Un genoma por especie (el mas reciente)
one_per_species = match.sort_values('seq_rel_date', ascending=False).drop_duplicates('organism_name')
one_per_species = one_per_species[one_per_species['ftp_path'] != 'na']

print("Genomas seleccionados: {}".format(len(one_per_species)))
print("Tamaño estimado: ~{} GB".format(round(len(one_per_species) * 4 / 1000, 1)))

# Guardar lista TSV
one_per_species[['organism_name', 'ftp_path']].to_csv(OUTPUT_TSV, sep='\t', index=False)
print("Lista guardada en {}".format(OUTPUT_TSV))

# Generar script de descarga bash
lines = ['#!/bin/bash', 'mkdir -p {}'.format(CHROMOSOMES_DIR), '']
for _, row in one_per_species.iterrows():
    ftp = row['ftp_path'].rstrip('/')
    name = ftp.split('/')[-1]
    url = '{}/{}_genomic.fna.gz'.format(ftp, name)
    lines.append('wget -q -P {} {}'.format(CHROMOSOMES_DIR, url))

with open(OUTPUT_SH, 'w') as f:
    f.write('\n'.join(lines))
print("Script de descarga guardado en {}".format(OUTPUT_SH))

# Guardar lista de URLs para wget -i
urls = []
for _, row in one_per_species.iterrows():
    ftp = row['ftp_path'].rstrip('/')
    name = ftp.split('/')[-1]
    url = '{}/{}_genomic.fna.gz'.format(ftp, name)
    urls.append(url)

with open(OUTPUT_URLS, 'w') as f:
    f.write('\n'.join(urls))
print("URLs guardadas: {}".format(len(urls)))
print("Para descargar ejecutar:")
print("nohup wget -q -P data/chromosomes/ --continue --tries=3 --timeout=30 -i data/urls.txt >> data/download_log.txt 2>&1 &")