#!/bin/bash
#SBATCH --job-name=hdc_coloc
#SBATCH --account=pi-gilad
#SBATCH --partition=gilad-hm
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=logs/coloc_%A_%a.log
# array index = line (1-based) of GWAS/coloc_tasks.tsv: <CellType> <RNA|ATAC> <chr>; 1-440
read CT MOD CHR <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" ${HDC_DIR}/GWAS/coloc_tasks.tsv)
${PYTHON} ${HDC_CODE}/07_coloc_smr.py $CT $MOD $CHR
