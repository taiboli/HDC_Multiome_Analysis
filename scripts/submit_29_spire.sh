#!/bin/bash
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --job-name=spire
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=logs/29_spire_%A_%a.log
# SPIRE runs (29_spire.py); RUNS = space-separated run names, one per array task; B = bootstraps (default 50)
R=($RUNS)
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
${PYTHON} -I ${HDC_CODE}/29_spire.py ${R[$SLURM_ARRAY_TASK_ID]} ${B:-50}
