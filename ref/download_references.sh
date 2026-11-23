set -e
R=${HDC_DIR}/PeakAnnot/ref
cd $R
B=https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final
for f in $(cat chromhmm18/list.txt); do [ -s chromhmm18/$f ] || curl -s --retry 3 -o chromhmm18/$f $B/$f; done
[ -s GRCh38-cCREs.V4.bed ] || curl -s --retry 3 -o GRCh38-cCREs.V4.bed https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed
[ -s phastConsElements100way.txt.gz ] || curl -s --retry 3 -o phastConsElements100way.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/phastConsElements100way.txt.gz
[ -s gencode.v32.basic.annotation.gtf.gz ] || curl -s --retry 3 -o gencode.v32.basic.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.basic.annotation.gtf.gz
echo DONE; ls chromhmm18 | wc -l; du -sh .
B15=https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final
for f in $(cat chromhmm15/list.txt); do [ -s chromhmm15/$f ] || curl -s --retry 3 -o chromhmm15/$f $B15/$f; done
echo DONE15; ls chromhmm15 | wc -l
