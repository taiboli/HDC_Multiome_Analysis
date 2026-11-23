#!/bin/bash
#SBATCH --job-name=HDC_covopt_RNA
#SBATCH --partition=gilad-hm
#SBATCH --account=pi-gilad
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --array=1-10
#SBATCH -o logs/covopt_RNA_%A_%a.log

module load R/4.2.0
export R_HOME=$(R RHOME)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
CT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${HDC_CODE}/celltypes_mapped.txt)
${PYTHON} ${HDC_CODE}/05_optimize_covariates.py "$CT" RNA
