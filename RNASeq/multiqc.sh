#!/bin/bash
#SBATCH --time=00:00:30
#SBATCH --account=def-chantalg
#SBATCH --mem=0.5G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-12

# Charger l'emplacement virtuel
module load python/3.11.5
virtualen --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r requirements.txt

DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" multiQC_dir.txt)
cd $DIR

multiqc .  -o /home/josdes34/projects/def-chantalg/josdes34/rna_seq_analysis/QC/fastqc/before/multiqc/
