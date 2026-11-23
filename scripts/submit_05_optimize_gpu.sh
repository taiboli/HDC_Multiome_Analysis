#!/bin/bash
#SBATCH --job-name=HDC_covopt
#SBATCH --partition=gpu
#SBATCH --account=pi-gilad
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --array=1-10
#SBATCH -o logs/covopt_%A_%a.log

module load R/4.2.0
export R_HOME=$(R RHOME)
CT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${HDC_CODE}/celltypes_mapped.txt)
${PYTHON} ${HDC_CODE}/05_optimize_covariates.py "$CT" ATAC
