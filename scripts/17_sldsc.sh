#!/bin/bash
#SBATCH --job-name=sldsc
#SBATCH --account=pi-gilad
#SBATCH --partition=gilad-hm
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/sldsc_%x_%A_%a.log
# S-LDSC (stratified LD score regression) for one GWAS (row SLURM_ARRAY_TASK_ID of GWAS/ldsc/sumstats/${TRAITS:-traits.txt})
# and one annotation batch: each test set is fitted on top of baselineLD v2.2 + the batch's control set:
#   spec : peakspec_<CT> | HDC_peaks ; genespec_<CT> | HDC_genes ; HDC_peaks and HDC_genes alone (baselineLD only)
#   qtl  : caQTLpeak_<X> | HDC_peaks_tested ; eGene_<X> | HDC_genes_tested ; the two controls alone
#   link : linked_<CT> | HDC_peaks ; speclinked_<CT> | HDC_peaks + peakspec_<CT>   (scMultiMap links)
#   tf   : tfprog_<CT> | HDC_peaks + peakspec_<CT> ; tfdyn_<lin> | HDC_peaks   (SPIRE TF programs)
# Reference LD: EAS for BBJ, EUR otherwise. No AFR reference in Zenodo 10515792, so PAGE / UKB_AFR are run with
# EUR LD scores (LD mismatch: interpret as exploratory).
# usage: sbatch --array=1-N%40 --export=ALL,BATCH=spec 17_sldsc.sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
REF=${LDSC_REF_DIR}
L=${HDC_DIR}/GWAS/ldsc
LDSC="${LDSC_PYTHON} ${LDSC_SCRIPT}"
T=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $L/sumstats/${TRAITS:-traits.txt})
if [[ $T == BBJ__* ]]; then
  POP=EAS; BL=$REF/1000G_Phase3_EAS_baselineLD_v2.2_ldscores/baselineLD.
  W=$(ls -d $REF/1000G_Phase3_EAS_weights_hm3_no_MHC/*/)weights.EAS.hm3_noMHC.
  FRQ="--not-M-5-50"
else
  POP=EUR; BL=$REF/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.
  W=$(ls -d $REF/1000G_Phase3_weights_hm3_no_MHC/*/)weights.hm3_noMHC.
  FRQ="--frqfile-chr $(ls -d $REF/1000G_Phase3_frq/*/)1000G.EUR.QC."
fi
O=$L/results/$BATCH/$T; mkdir -p $O
LS=$L/ldscores/$POP
run () {  # name, ref-ld prefixes
  [ -s $O/$1.results ] && return
  $LDSC --h2 $L/sumstats/$T.sumstats.gz --ref-ld-chr $2 --w-ld-chr $W --overlap-annot $FRQ \
        --print-coefficients --out $O/$1 > /dev/null 2>&1 || echo "FAILED $T $1"
}
if [ "$BATCH" = spec ]; then
  run baseline $BL
  for C in HDC_peaks HDC_genes; do run $C $BL,$LS/$C/$C.; done
  for X in $(ls $LS | grep "^peakspec_"); do run $X $BL,$LS/HDC_peaks/HDC_peaks.,$LS/$X/$X.; done
  for X in $(ls $LS | grep "^genespec_"); do run $X $BL,$LS/HDC_genes/HDC_genes.,$LS/$X/$X.; done
elif [ "$BATCH" = tf ]; then  # SPIRE TF programs (31_tfnet_analysis.py): beyond all peaks + the cell type's specific peaks
  for X in $(ls $LS | grep "^tfprog_"); do C=${X#tfprog_}; run $X $BL,$LS/HDC_peaks/HDC_peaks.,$LS/peakspec_$C/peakspec_$C.,$LS/$X/$X.; done
  for X in $(ls $LS | grep "^tfdyn_"); do run $X $BL,$LS/HDC_peaks/HDC_peaks.,$LS/$X/$X.; done
elif [ "$BATCH" = link ]; then
  for X in $(ls $LS | grep "^linked_"); do run $X $BL,$LS/HDC_peaks/HDC_peaks.,$LS/$X/$X.; done
  for X in $(ls $LS | grep "^speclinked_"); do C=${X#speclinked_}; run $X $BL,$LS/HDC_peaks/HDC_peaks.,$LS/peakspec_$C/peakspec_$C.,$LS/$X/$X.; done
else
  for C in HDC_peaks_tested HDC_genes_tested; do run $C $BL,$LS/$C/$C.; done
  for X in $(ls $LS | grep "^caQTLpeak_"); do run $X $BL,$LS/HDC_peaks_tested/HDC_peaks_tested.,$LS/$X/$X.; done
  for X in $(ls $LS | grep "^eGene_"); do run $X $BL,$LS/HDC_genes_tested/HDC_genes_tested.,$LS/$X/$X.; done
fi
echo done
