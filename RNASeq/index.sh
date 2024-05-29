#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --account=def-chantalg
#SBATCH --mem=35G


<<index
 The present code is to indexed the transcriptome for further use in the quantification process (kallisto)
 Before head, you must fetch to correct file corresponding to the transcriptome.
 To do so, you could run the following command line:

            wget https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

 Visit the ENSEMBL website for other files (or any other website)
index


# Load modules for dependencies of kallisto and kallisto itself
module load StdEnv/2020 gcc/9.3.0 intel/2020.1.217 kallisto/0.46.1


# Indexation of the reference transcriptome
kallisto index --index=Homo_sapiens.GRCh38.cdna.all.index Homo_sapiens.GRCh38.cdna.all.fa.gz 1>>log.txt 2>>err.txt
