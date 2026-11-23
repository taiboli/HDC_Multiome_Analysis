#!/bin/bash
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --job-name=scmultimap
#SBATCH --array=0-13
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/23_scmultimap_%A_%a.log
# scMultiMap peak-gene links: 11 cell types, AllCells pooled, and 2 within-donor RNA-permuted nulls
T=(Ectoderm_1 Ectoderm_2 Ectoderm_3 Endoderm Mesoderm Neural_crest Neuronal Retinal Unknown Vascular_endothelial iPSC
   AllCells "Neuronal --perm" "iPSC --perm")
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
${PYTHON} -I ${HDC_CODE}/23_scmultimap.py ${T[$SLURM_ARRAY_TASK_ID]}
