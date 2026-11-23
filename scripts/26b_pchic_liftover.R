## Promoter-capture Hi-C significant promoter-other interactions (Jung et al. 2019 Nat Genet, GSE86189_all_interaction.po,
## 27 cell types / tissues incl. H1 ESC and H1-derived mesendoderm (ME), MSC, NPC, trophoblast (TB); hg19 HindIII fragments).
## Target fragments lifted to hg38 (range of lifted pieces, same chromosome, <= 2x original length); gene names split on ",".
## Output: PeakAnnot/ref/pchic/pchic_po_hg38.tsv.gz (cell, gene, chr, start, end, dist, dist_pvalue)
suppressMessages({library(data.table); library(GenomicRanges); library(rtracklayer)})
d = fread(paste0(Sys.getenv("HDC_DIR"), "/PeakAnnot/ref/pchic/GSE86189_all_interaction.po.txt.gz"), skip = 1, header = FALSE,
          col.names = c("cell", "chr", "gene", "target", "raw", "dist", "all_capture_res", "dist_res", "dist_pvalue"))  # header is space-delimited
t = tstrsplit(d$target, ".", fixed = TRUE)
gr = GRanges(t[[1]], IRanges(as.integer(t[[2]]) + 1, as.integer(t[[3]])), row = seq_len(nrow(d)))
ch = import.chain(paste0(Sys.getenv("HDC_DIR"), "/PeakAnnot/ref/hg19ToHg38.over.chain"))
lo = unlist(liftOver(unique(granges(gr)) |> (\(u) { u$key = paste(seqnames(u), start(u), end(u)); u })(), ch))
lo = as.data.table(lo)[, .(chr = unique(seqnames)[1], nchr = uniqueN(seqnames), start = min(start) - 1L, end = max(end)), by = key]
w = as.data.table(gr)[, key := paste(seqnames, start, end)][, .(row, key, w0 = width)]
m = merge(w, lo[nchr == 1], by = "key")[(end - start) <= 2 * w0]
out = cbind(d[m$row, .(cell, gene, dist, dist_pvalue)], m[, .(chr, start, end)])
out = out[, .(gene = unlist(strsplit(gene, ","))), by = .(cell, chr, start, end, dist, dist_pvalue)][gene != ""]
fwrite(out, paste0(Sys.getenv("HDC_DIR"), "/PeakAnnot/ref/pchic/pchic_po_hg38.tsv.gz"), sep = "\t")
cat(nrow(d), "interactions;", nrow(out), "gene-level hg38 rows;", length(unique(out$cell)), "cells\n")
