#!/bin/bash
#SBATCH --job-name=hdc_universe
#SBATCH --account=pi-gilad
#SBATCH --partition=gilad-hm
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=logs/universe_%j.log
module load R/4.2.0
Rscript ${HDC_CODE}/06a_variant_universe.R
