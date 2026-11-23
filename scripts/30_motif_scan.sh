#!/bin/bash
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --job-name=motifscan
#SBATCH --array=0-23
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=logs/30_motif_%A_%a.log
# HOMER known vertebrate motif instances (default log-odds thresholds) in the central 200 bp of every HDC peak
# (24 chunks of PeakAnnot/sc/features_ATAC.tsv.gz). Output TFNet/motif/find/c<NN>.txt (PositionID, offset, seq, motif, strand, score)
E=${TOOLS_BIN}
D=${HDC_DIR}/TFNet/motif
FA=${HDC_DIR}/PeakAnnot/ref/hg38.fa
M=${HOMER_MOTIFS}
C=$(printf "c%02d" $SLURM_ARRAY_TASK_ID)
export PATH=$E:$PATH
mkdir -p $D/find $D/pre_$C
findMotifsGenome.pl $D/chunks/$C.bed $FA $D/tmp_$C -size 200 -find $M -p 8 -preparsedDir $D/pre_$C > $D/find/$C.txt 2> $D/find/$C.err
rm -rf $D/pre_$C $D/tmp_$C
