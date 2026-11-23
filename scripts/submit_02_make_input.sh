#!/bin/bash
#SBATCH --job-name=HDC_qtl_input
#SBATCH --partition=gilad-hm
#SBATCH --account=pi-gilad
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --array=1-22
#SBATCH -o logs/make_input_%A_%a.log

module load R/4.2.0
read CT MOD < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${HDC_CODE}/tasks.tsv)
Rscript ${HDC_CODE}/02_make_tensorqtl_input.R "$CT" "$MOD"
