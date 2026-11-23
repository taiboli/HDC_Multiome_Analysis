## Lift all HDC multiome features (524k peaks, 36.6k genes; hg38) to hg19 for S-LDSC (1000G Phase 3 reference is hg19).
## Features whose interval splits or maps to another chromosome are dropped; the hg19 interval is the range of the
## lifted pieces when they stay on one chromosome within 2x the original length.
## Output: GWAS/ldsc/features_hg19_{ATAC,RNA}.tsv.gz (id, chr, start, end; chr without "chr")
suppressMessages({library(data.table); library(GenomicRanges); library(rtracklayer)})
base = paste0(Sys.getenv("HDC_DIR"), "/")
ch = import.chain(paste0(Sys.getenv("GWAS_DIR"), "/hg38ToHg19.over.chain"))
dir.create(paste0(base, "GWAS/ldsc"), showWarnings = FALSE)
for (mod in c("ATAC", "RNA")) {
  f = fread(paste0(base, "PseudoBulk/features_", mod, ".tsv.gz"))
  f = f[chr %in% paste0("chr", 1:22)]
  gr = GRanges(f$chr, IRanges(f$start + 1, f$end), id = f$id)  # features files are 0-based starts (BED-like)
  lo = unlist(liftOver(gr, ch))
  d = as.data.table(lo)[, .(chr = as.character(seqnames[1]), nchr = uniqueN(seqnames), start = min(start), end = max(end)), by = id]
  len0 = setNames(width(gr), gr$id)
  d = d[nchr == 1 & (end - start + 1) <= 2 * len0[id]]
  d[, chr := sub("chr", "", chr)]
  fwrite(d[, .(id, chr, start, end)], paste0(base, "GWAS/ldsc/features_hg19_", mod, ".tsv.gz"), sep = "\t")
  cat(mod, nrow(f), "->", nrow(d), "\n")
}
