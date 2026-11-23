#!/bin/bash
#SBATCH --job-name=gwas_harm
#SBATCH --account=pi-gilad
#SBATCH --partition=gilad-hm
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --output=logs/harmonize_%A_%a.log
# array index = row of GWAS/gwas_manifest.tsv (0-based); default 0-189
${PYTHON} ${HDC_CODE}/06b_harmonize_gwas.py run $SLURM_ARRAY_TASK_ID
