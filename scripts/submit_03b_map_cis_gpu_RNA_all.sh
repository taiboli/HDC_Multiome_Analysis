#!/bin/bash
#SBATCH --job-name=HDC_tensorqtl_RNA
#SBATCH --partition=gpu
#SBATCH --account=pi-gilad
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH -o logs/map_cis_RNA_all_%j.log

# All RNA (eQTL) tasks sequentially in one GPU job (gpu QOS allows max 12 submitted jobs per user).
for CT in $(cat ${HDC_CODE}/celltypes.txt); do
  ${PYTHON} ${HDC_CODE}/03_map_cis.py "$CT" RNA
done
