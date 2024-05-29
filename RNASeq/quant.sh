#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --account=def-chantalg
#SBATCH --mem=5G
#SBATCH --array=1-72
#SBATCH --cpus-per-task=2


# Charger les modules nécessaires à l'utilisation de kallisto
module load StdEnv/2020 gcc/9.3.0 intel/2020.1.217 kallisto/0.46.1

DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" working_dir.txt)

# Générer le nom de chacun des reads d'un même échantillon
PRE=$( grep $DIR working_files.txt )
file1="$PRE"_R1.fastq.gz
file2="$PRE"_R2.fastq.gz

# Faire la quantification
kallisto quant \
	--threads 8 \
	--index=../index/Homo_sapiens.GRCh38.cdna.all.index \
	--output-dir=$DIR \
	"$file1" "$file2" \
	1>$DIR/log.txt 2>$DIR/err.txt

