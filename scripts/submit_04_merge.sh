#!/bin/bash
#SBATCH --job-name=HDC_qtl_merge
#SBATCH --partition=gilad-hm
#SBATCH --account=pi-gilad
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH -o logs/merge_%j.log

module load R/4.2.0
export R_HOME=$(R RHOME)
COVSET=${COVSET:-} ${PYTHON} ${HDC_CODE}/04_merge_qvalues.py
