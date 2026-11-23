# companion of 07v_validate_coloc.py: rerun coloc.abf on exported pairs and compare with the Python values
args = commandArgs(trailingOnly = TRUE)
od = args[1]
suppressMessages(library(coloc))
meta = read.delim(file.path(od, "meta.tsv"))
for (i in seq_len(nrow(meta))) {
  d = read.delim(file.path(od, paste0("pair", meta$pair[i], ".tsv")))
  D1 = list(beta = d$slope, varbeta = d$varbeta_qtl, type = "quant", sdY = 1, snp = d$variant_id)
  D2 = list(beta = d$beta, varbeta = d$varbeta_gwas, type = meta$type[i], snp = d$variant_id)
  if (meta$type[i] == "quant") D2$sdY = meta$sdY[i] else D2$s = 0.5
  r = suppressMessages(coloc.abf(D1, D2))$summary
  cat(sprintf("%-22s %-28s nsnps py=%d R=%d  H3 py=%.4f R=%.4f  H4 py=%.4f R=%.4f\n", meta$trait[i], meta$feature[i],
              meta$nsnps[i], as.integer(r[["nsnps"]]), meta$H3[i], r[["PP.H3.abf"]], meta$H4[i], r[["PP.H4.abf"]]))
}
