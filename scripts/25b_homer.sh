#!/bin/bash
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --job-name=homer
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=6:00:00
#SBATCH --output=logs/25b_homer_%A_%a.log
# HOMER known motifs (vertebrate known library, GC-normalized vs the given background), one set per array task.
# Sets = PeakAnnot/homer/input/*.target.bed (sorted); usage: sbatch --array=0-<n-1> 25b_homer.sh
E=${TOOLS_BIN}
D=${HDC_DIR}/PeakAnnot/homer
FA=${HDC_DIR}/PeakAnnot/ref/hg38.fa
SETS=($(ls $D/input/*.target.bed | xargs -n1 basename | sed 's/.target.bed//' | sort))
S=${SETS[$SLURM_ARRAY_TASK_ID]}
mkdir -p $D/out $D/preparsed_$S
export PATH=$E:$PATH
findMotifsGenome.pl $D/input/$S.target.bed $FA $D/out/$S -size 200 -bg $D/input/$S.bg.bed -nomotif -p 8 \
  -preparsedDir $D/preparsed_$S
rm -rf $D/preparsed_$S
