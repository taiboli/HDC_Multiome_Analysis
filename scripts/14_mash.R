## mashr across the 10 HDC cell types (one modality), on z-scores from 13_mash_extract.py
##   Bhat = z, Shat = 1 (z-scale mash); feature not tested in a cell type -> Bhat = 0, Shat = 1e3 (no information),
##   and its lfsr / posterior are set to NA in the output.
##   null correlation (sample overlap: same 17 donors in every cell type): estimate_null_correlation_simple on
##     complete random pairs
##   covariances: canonical + ED-refined PCA (5) from strong pairs with max |z| > 4 (top 5,000 if fewer)
##   mixture weights fitted on complete-case random pairs (tested in all 10 cell types; <= 200k), posteriors on all
##   strong pairs (top variant per feature; missing cell types uninformative)
## Output (mash/): <MOD>_mash_strong.parquet (phenotype_id, variant_id, lfsr_<CT>, pm_<CT>, z_<CT>),
##                 <MOD>_mash_summary.tsv (per cell type: tested, FDR10 from tensorQTL, lfsr<0.05/0.1),
##                 <MOD>_mash_sharing.tsv (pairwise sharing by sign and magnitude within 2x), <MOD>_mash_model.rds
## Usage: Rscript 14_mash.R <RNA|ATAC>
suppressMessages({library(mashr); library(data.table); library(arrow)})
set.seed(2025)
mod = commandArgs(TRUE)[1]
base = paste0(Sys.getenv("HDC_DIR"), "/mash/")
cts = c("Ectoderm_1", "Ectoderm_3", "Endoderm", "iPSC", "Mesoderm", "Neural_crest", "Neuronal", "Retinal",
        "Unknown", "Vascular_endothelial")

load_set = function(kind) {
  fs = Sys.glob(sprintf("%s/extract/%s_chr*_%s.parquet", base, mod, kind))
  d = rbindlist(lapply(fs, function(f) as.data.table(read_parquet(f))))
  # an rsID mapping to >1 site gives duplicated feature x variant x cell type rows: ambiguous, drop them
  d[, n := .N, by = .(phenotype_id, variant_id, celltype)]
  cat(kind, "dropped duplicated rows:", sum(d$n > 1), "\n")
  d = d[n == 1][, n := NULL]
  w = dcast(d, phenotype_id + variant_id ~ celltype, value.var = "z", fun.aggregate = NULL)
  for (ct in setdiff(cts, names(w))) w[[ct]] = NA_real_
  w
}
to_mash = function(w, V) {
  Z = as.matrix(w[, ..cts])
  S = matrix(1, nrow(Z), ncol(Z), dimnames = dimnames(Z))
  miss = is.na(Z)
  Z[miss] = 0; S[miss] = 1e3
  list(data = mash_set_data(Z, S, V = V), miss = miss)
}

strong = load_set("strong"); random = load_set("random")
cat(mod, "strong pairs:", nrow(strong), " random pairs:", nrow(random), "\n")
if (nrow(random) > 2e5) random = random[sample(.N, 2e5)]

rc = random[complete.cases(random[, ..cts])]
cat("complete random pairs for null correlation:", nrow(rc), "\n")
Vhat = estimate_null_correlation_simple(mash_set_data(as.matrix(rc[, ..cts]), matrix(1, nrow(rc), length(cts))))
dimnames(Vhat) = list(cts, cts)
print(round(Vhat, 2))

maxz = apply(abs(as.matrix(strong[, ..cts])), 1, max, na.rm = TRUE)
top = which(maxz > 4)
if (length(top) < 5000) top = order(-maxz)[1:min(5000, nrow(strong))]
cat("strong subset for data-driven covariances:", length(top), "\n")
d.top = to_mash(strong[top], Vhat)$data
U.pca = cov_pca(d.top, 5)
U.ed = cov_ed(d.top, U.pca)
U.c = cov_canonical(d.top)

## mixture weights from complete-case random pairs: with partially observed pairs (most ATAC pairs are tested in
## only a few cell types) the fit put ~50% weight on small cell-type-specific components (null weight 0.3-0.5 vs
## 0.84 on complete cases; mash/debug/dbg.log), i.e. missingness rather than signal
d.rand = to_mash(rc, Vhat)$data
m = mash(d.rand, Ulist = c(U.ed, U.c), outputlevel = 1)
g = get_fitted_g(m)
w = get_estimated_pi(m)
print(round(sort(w[w > 0.01], decreasing = TRUE), 3))

## posteriors on all strong pairs in chunks
res = list()
idx = split(seq_len(nrow(strong)), ceiling(seq_len(nrow(strong)) / 50000))
for (i in seq_along(idx)) {
  x = to_mash(strong[idx[[i]]], Vhat)
  mm = mash(x$data, g = g, fixg = TRUE, outputlevel = 2)
  lf = get_lfsr(mm); pm = get_pm(mm)
  lf[x$miss] = NA; pm[x$miss] = NA
  colnames(lf) = paste0("lfsr_", cts); colnames(pm) = paste0("pm_", cts)
  res[[i]] = data.table(strong[idx[[i]], .(phenotype_id, variant_id)], lf, pm)
  cat("chunk", i, "/", length(idx), "\n")
}
res = rbindlist(res)
zz = strong[, ..cts]; setnames(zz, paste0("z_", cts))
res = cbind(res, zz)
write_parquet(res, sprintf("%s/%s_mash_strong.parquet", base, mod))
fwrite(res, sprintf("%s/%s_mash_strong.tsv.gz", base, mod), sep = "\t")  # readable by older pyarrow
saveRDS(list(Vhat = Vhat, g = g, pi = w, U.ed = U.ed), sprintf("%s/%s_mash_model.rds", base, mod))

## per cell type summary (tensorQTL FDR from the permutation q-values for comparison)
qd = if (mod == "RNA") "eQTL" else "caQTL"
summ = rbindlist(lapply(cts, function(ct) {
  q = fread(sprintf(paste0(Sys.getenv("HDC_DIR"), "/%s_Output_optCov/%s/allchr_cis_noBootstrap_TMM_invnorm.cis_qtl.qval.txt.gz"), qd, ct))
  lf = res[[paste0("lfsr_", ct)]]
  data.table(celltype = ct, modality = mod, n_tested = sum(!is.na(lf)), tensorqtl_FDR10 = sum(q$qval < 0.1, na.rm = TRUE),
             mash_lfsr05 = sum(lf < 0.05, na.rm = TRUE), mash_lfsr10 = sum(lf < 0.1, na.rm = TRUE))
}))
L = as.matrix(res[, paste0("lfsr_", cts), with = FALSE])
summ = rbind(summ, data.table(celltype = "any", modality = mod, n_tested = nrow(res), tensorqtl_FDR10 = NA,
                              mash_lfsr05 = sum(apply(L, 1, min, na.rm = TRUE) < 0.05),
                              mash_lfsr10 = sum(apply(L, 1, min, na.rm = TRUE) < 0.1)))
fwrite(summ, sprintf("%s/%s_mash_summary.tsv", base, mod), sep = "\t")
print(summ)

## pairwise sharing (same sign, within 2-fold) among pairs significant in either cell type
PM = as.matrix(res[, paste0("pm_", cts), with = FALSE])
sh = matrix(NA, length(cts), length(cts), dimnames = list(cts, cts))
for (i in seq_along(cts)) for (j in seq_along(cts)) {
  k = which((L[, i] < 0.05 | L[, j] < 0.05) & !is.na(L[, i]) & !is.na(L[, j]))
  if (length(k) < 10) next
  r = PM[k, i] / PM[k, j]
  sh[i, j] = mean(r > 0.5 & r < 2)
}
fwrite(data.table(celltype = cts, sh), sprintf("%s/%s_mash_sharing.tsv", base, mod), sep = "\t")
print(round(sh, 2))
cat("done\n")
