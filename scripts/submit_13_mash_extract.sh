#!/bin/bash
#SBATCH --job-name=mash_extract
#SBATCH --account=pi-gilad
#SBATCH --partition=gilad-hm
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=4:00:00
#SBATCH --array=1-44%12
#SBATCH --output=logs/mash_extract_%A_%a.log
read MOD CHR < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${HDC_CODE}/mash_tasks.tsv)
${PYTHON} ${HDC_CODE}/13_mash_extract.py $MOD $CHR
