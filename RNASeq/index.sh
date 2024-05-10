#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --account=def-chantalg
#SBATCH --mem=35G

# Charger les modules importants pour le fonctionnement de Kallisto
module load StdEnv/2020 gcc/9.3.0 intel/2020.1.217 kallisto/0.46.1

# Indexation du transcriptome de référence
kallisto index --index=Homo_sapiens.GRCh38.cdna.all.index Homo_sapiens.GRCh38.cdna.all.fa.gz 1>>log.txt 2>>err.txt
