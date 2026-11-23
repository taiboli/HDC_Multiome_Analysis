## Make tensorQTL phenotype BED + covariates for one HDC cell type and modality.
## Same logic as EB_CellTypeData (EBQTL/RNA/EB_RNA_tensorQTL_TMM_MakeInput.R and
## EB_CellTypeData/EB_Celltypedata_Assign.R):
##   donors with >= 20 cells; skip cell type if < 10 donors
##   features detected (>0) in >= 5 cells AND >= 1% of cells in >= 80% of donors
##   pseudobulk sum -> TMM log-CPM -> per-feature inverse normal transform
##   covariates: Sex + PhenoPC (princomp loadings) + 5 GenoPC
## Only change: number of PhenoPCs is reduced when there are few donors so that
## tensorQTL residual dof (n - 2 - n_cov) stays >= 3 (17 donors -> 5 PCs, as before).
##
## Usage: Rscript 02_make_tensorqtl_input.R <CellType> <RNA|ATAC>

suppressMessages({
  library(data.table)
  library(dplyr)
  library(edgeR)
})

args = commandArgs(trailingOnly = TRUE)
this.label = args[1]
this.mod = args[2]

base.dir = paste0(Sys.getenv("HDC_DIR"), "/")
pb.dir = paste0(base.dir, "PseudoBulk/")

if(this.mod == "RNA"){
  this.input.dir = paste0(base.dir, "eQTL/", this.label, "/")
  expr.file = "tensorqtl.geneExpr_noBootstrap_TMM_invnorm.bed"
  cov.file = "tensorqtl.geneCov_noBootstrap_TMM_invnorm.txt"
  keep.chr = paste0("chr", 1:22)
} else {
  this.input.dir = paste0(base.dir, "caQTL/", this.label, "/")
  expr.file = "tensorqtl.atacExpr_noBootstrap_TMM_invnorm.bed"
  cov.file = "tensorqtl.atacCov_noBootstrap_TMM_invnorm.txt"
  keep.chr = paste0("chr", c(1:22, "X"))
}

invnorm_transform <- function(x) {
  qqnorm(x, plot.it=F)$x
}

ncell = fread(paste0(pb.dir, this.label, "_NCell.tsv"))
to.keep.donors = ncell$donor[which(ncell$ncell >= 20)]
if(length(to.keep.donors) < 10){
  cat(this.label, this.mod, "SKIPPED: only", length(to.keep.donors), "donors with >= 20 cells\n")
  quit(save = "no", status = 0)
}

inds.ordered = scan(paste0(Sys.getenv("GENOTYPE_DIR"), "/EB_Multiome_highpass_individuals.txt"), what = character(), quiet = T)
inds.ordered = intersect(inds.ordered, to.keep.donors)
n.ind = length(inds.ordered)

features = fread(paste0(pb.dir, "features_", this.mod, ".tsv.gz"))
sum.data = fread(paste0(pb.dir, this.label, "_", this.mod, "_SumCount.tsv.gz"))
nz.data = fread(paste0(pb.dir, this.label, "_", this.mod, "_NCellExpr.tsv.gz"))
stopifnot(identical(sum.data$id, features$id), identical(nz.data$id, features$id))

chr.keep = which(features$chr %in% keep.chr)
features = features[chr.keep]
sum.data = sum.data[chr.keep]
nz.data = nz.data[chr.keep]

cells.count = as.matrix(nz.data[, ..inds.ordered])
cells.perc = sweep(cells.count, 2, ncell$ncell[match(inds.ordered, ncell$donor)], "/")

tmp = rowSums(cells.count >= 5)
tmp2 = rowSums(cells.perc >= 0.01)
feat.keep = which(tmp >= 0.8 * n.ind & tmp2 >= 0.8 * n.ind)

counts = as.matrix(sum.data[feat.keep, ..inds.ordered])
features = features[feat.keep]

# expression normalization
# 1. TMM normalization to get log TMM-normalized counts
expression.mat <- DGEList(counts=counts) %>% calcNormFactors(method="TMM") %>% cpm(log=T)
# 2. normalize each feature to have zero mean, unit variance across individuals
expression.mat <- apply(expression.mat, 1, invnorm_transform) %>% t
rownames(expression.mat) <- features$id
colnames(expression.mat) <- inds.ordered

all.data = data.table(expression.mat)
all.data$phenotype_id = features$id
all.data[["#chr"]] = gsub("chr", "", features$chr)
all.data$start = as.integer(features$start)
all.data$end = as.integer(features$end)

all.data$`#chr` = factor(all.data[["#chr"]], levels = c(1:22, "X"))
all.data = setorderv(all.data, c("#chr", "start", "end"))

dir.create(this.input.dir, showWarnings = F, recursive = T)
to.choose <- c("#chr", "start", "end", "phenotype_id", inds.ordered)
fwrite(all.data[, ..to.choose], file = paste0(this.input.dir, expr.file),
       col.names = T, row.names = F, sep = "\t", quote = F)
system(paste0("gzip -f ", this.input.dir, expr.file))


###################
## Cov
###################

n.pheno.pc = max(0, min(5, n.ind - 2 - 3 - 1 - 5))

this.cov = NULL
yri.cov = fread(paste0(Sys.getenv("GENOTYPE_DIR"), "/igsr-yri.tsv.tsv"))
yri.cov$Gender = ifelse(yri.cov$Sex == "female", 0, 1)
yri.cov = subset(yri.cov, `Sample name` %in% inds.ordered)
this.cov = as.data.frame(t(yri.cov[, c("Gender")]))
names(this.cov) = yri.cov$`Sample name`
this.cov = this.cov[, inds.ordered]
this.cov$id = c("Sex")

if(n.pheno.pc > 0){
  this.pca <- princomp(all.data[, ..inds.ordered])
  this.pca <- as.data.frame(t(this.pca$loadings[, 1:n.pheno.pc, drop = F]))
  this.pca = this.pca[, inds.ordered]
  this.pca$id = paste0("PhenoPC_", 1:n.pheno.pc)
  this.cov = rbind(this.cov, this.pca)
}

genotype.pca <- fread(Sys.getenv("GENOTYPE_PCS"))
genotype.pca = subset(genotype.pca, V2 %in% inds.ordered)
geno.inds <- genotype.pca$V2
genotype.pca <- as.data.frame(t(genotype.pca[,3:ncol(genotype.pca)]))
names(genotype.pca) <- geno.inds
genotype.pca = genotype.pca[, inds.ordered]
genotype.pca$id <- paste0("GenoPC_", 1:nrow(genotype.pca))
this.cov = rbind(this.cov, genotype.pca[1:5, ])

to.choose = c("id", inds.ordered)
setDT(this.cov)
fwrite(this.cov[, ..to.choose], file = paste0(this.input.dir, cov.file),
       col.names = T, row.names = F, sep = "\t", quote = F)

fwrite(data.table(CellType = this.label, Modality = this.mod, nDonor = n.ind,
                  nFeature = nrow(all.data), nPhenoPC = n.pheno.pc,
                  nCov = nrow(this.cov), Donors = paste(inds.ordered, collapse = ",")),
       file = paste0(this.input.dir, "input_summary.tsv"), sep = "\t")
cat(this.label, this.mod, "donors:", n.ind, "features:", nrow(all.data), "PhenoPCs:", n.pheno.pc, "\n")
