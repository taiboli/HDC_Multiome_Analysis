#!/bin/bash
#SBATCH --job-name=HDC_pseudobulk
#SBATCH --partition=gilad-hm
#SBATCH --account=pi-gilad
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH -o logs/pseudobulk_%j.log

${PYTHON_SCVI} -I ${HDC_CODE}/01_pseudobulk.py
