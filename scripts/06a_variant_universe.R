# HDC genotype variant universe (hg38, rsID) with hg19 positions, for matching GWAS summary statistics.
# Output: HDC_Filtered/GWAS/variant_universe.parquet
#   vidx (0-based global index), chr, pos38, id (rsID = tensorQTL variant_id), ref, alt, pos19 (NA if unmapped)
# Variants with ID "." (308 genome-wide) are dropped: they cannot be matched to tensorQTL output.
suppressMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(arrow)
})

ch = import.chain(paste0(Sys.getenv("GWAS_DIR"), "/hg38ToHg19.over.chain"))
out = list()
for (chr in 1:22) {
  pvar = fread(cmd = paste0("grep -v '^##' ", Sys.getenv("GENOTYPE_DIR"), "/Genotype/plink/chr", chr, ".pvar"),
               select = 1:5, col.names = c("chr", "pos38", "id", "ref", "alt"))
  pvar = pvar[id != "."]
  gr = GRanges(paste0("chr", chr), IRanges(pvar$pos38, width = 1))
  lo = liftOver(gr, ch)
  ok = lengths(lo) == 1
  lo1 = unlist(lo[ok])
  pvar$pos19 = NA_integer_
  same.chr = as.character(seqnames(lo1)) == paste0("chr", chr)
  pvar$pos19[which(ok)[same.chr]] = start(lo1)[same.chr]
  out[[chr]] = pvar
  cat("chr", chr, nrow(pvar), "variants;", sum(!is.na(pvar$pos19)), "lifted to hg19\n")
}
u = rbindlist(out)
u[, chr := as.integer(chr)]
u[, vidx := 0:(.N - 1)]
setcolorder(u, c("vidx", "chr", "pos38", "id", "ref", "alt", "pos19"))
dir.create(paste0(Sys.getenv("HDC_DIR"), "/GWAS"), showWarnings = FALSE)
write_parquet(u, paste0(Sys.getenv("HDC_DIR"), "/GWAS/variant_universe.parquet"))
cat("total", nrow(u), "\n")
