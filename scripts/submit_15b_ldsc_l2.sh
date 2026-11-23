#!/bin/bash
#SBATCH --job-name=ldsc_l2
#SBATCH --account=pi-gilad
#SBATCH --partition=gilad-hm
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --array=1-22
#SBATCH --output=logs/ldsc_l2_%x_%A_%a.log
# usage: sbatch --export=ALL,BATCH=spec,POP=EUR submit_15b_ldsc_l2.sh
set -e
CHR=$SLURM_ARRAY_TASK_ID
S=${HDC_CODE}
REF=${LDSC_REF_DIR}
A=${HDC_DIR}/GWAS/ldsc/annot/$BATCH/$POP
PY=${PYTHON}
if [ "$POP" = EUR ]; then
  BF=$REF/1000G_Phase3_plinkfiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.$CHR
  BL=$REF/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.$CHR.l2.ldscore.gz
else
  BF=$(ls $REF/1000G_Phase3_EAS_plinkfiles/*/1000G.EAS.QC.$CHR.bim | sed 's/.bim$//')
  BL=$(ls $REF/1000G_Phase3_EAS_baselineLD_v2.2_ldscores/*$CHR.l2.ldscore.gz | grep "\.$CHR\.l2")
fi
$PY $S/15b_ldsc_annot.py make $BATCH $POP $CHR
zcat $BL | tail -n +2 | cut -f2 > $A/printsnps.$CHR.txt
${LDSC_PYTHON} ${LDSC_SCRIPT} \
  --l2 --bfile $BF --ld-wind-cm 1 --annot $A/all.$CHR.annot.gz --print-snps $A/printsnps.$CHR.txt --out $A/all.$CHR
$PY $S/15b_ldsc_annot.py split $BATCH $POP $CHR
echo done
