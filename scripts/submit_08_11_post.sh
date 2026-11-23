#!/bin/bash
#SBATCH --job-name=hdc_gwas_post
#SBATCH --account=pi-gilad
#SBATCH --partition=gilad-hm
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=logs/gwas_post_%j.log
# summary / SMR (08), lead-variant enrichment (09), peak-gene links + GTEx (10), cell-type peak enrichment (11)
P=${PYTHON}
S=${HDC_CODE}

$P $S/08_summarize_coloc_smr.py
$P $S/09_enrichment_qtl_leads.py
$P $S/10_multiome_link_gtex.py
$P $S/11_celltype_peak_enrichment.py
