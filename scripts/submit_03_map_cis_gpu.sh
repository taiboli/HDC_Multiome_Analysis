#!/bin/bash
#SBATCH --job-name=HDC_tensorqtl
#SBATCH --partition=gpu
#SBATCH --account=pi-gilad
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --array=1-22
#SBATCH -o logs/map_cis_%A_%a.log

# One task per cell type x modality (rows of tasks.tsv); loops chr1-22, skipping finished chromosomes.
# Optional: sbatch --export=ALL,CHRS=22 to restrict chromosomes (pilot).
read CT MOD < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${HDC_CODE}/tasks.tsv)
${PYTHON} ${HDC_CODE}/03_map_cis.py "$CT" "$MOD" ${CHRS:-}
