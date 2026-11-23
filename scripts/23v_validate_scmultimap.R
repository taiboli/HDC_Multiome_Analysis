# Run the original scMultiMap R functions on the subset written by 23v_validate_scmultimap.py write
V <- paste0(Sys.getenv("HDC_DIR"), "/PeakAnnot/scMultiMap/validate/")
for (f in list.files(paste0(Sys.getenv("HDC_DIR"), "/PeakAnnot/scMultiMap/Rsrc/"), full.names = TRUE)) source(f)
cells <- read.delim(file.path(V, "cells.tsv"))
G <- as.matrix(Matrix::readMM(file.path(V, "gene.mtx"))); colnames(G) <- readLines(file.path(V, "genes.txt"))
P <- as.matrix(Matrix::readMM(file.path(V, "peak.mtx"))); colnames(P) <- readLines(file.path(V, "peaks.txt"))
ig <- scMultiMap_IRLS(G, cells$sd_g, bsample = cells$donor, irls = TRUE)
ip <- scMultiMap_IRLS(P, cells$sd_p, bsample = cells$donor, irls = TRUE)
pairs <- expand.grid(peak = colnames(P), gene = colnames(G), stringsAsFactors = FALSE)[, c("gene", "peak")]
res <- scMultiMap_WLS(list(gene = G, peak = P), list(cells$sd_g, cells$sd_p), pairs, list(ig, ip))
write.table(res, file.path(V, "R_wls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sigma_sq = ig$sigma_sq), file.path(V, "R_irls_gene.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(as.character(unique(ig$theta)), file.path(V, "R_theta.txt"))
