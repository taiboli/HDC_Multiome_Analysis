# Figure legends

Files: `Main/Figure<n>.{pdf,png}` and `Supplementary/FigureS<n>.{pdf,png}`. The main legends match the manuscript doc. The supplementary legends are new and follow the figure list in the doc.

HDC: human iPSC-derived heterogeneous differentiating cultures from 17 YRI lines (128,588 paired multiome cells). CT: cell type. "Vasc_endo" = Vascular_endothelial.

## Main figures

### Figure 1. Multiome profiling defines and links cell-type-specific developmental cis-regulatory elements.

- **a** ENCODE cCRE v4 class of 524,022 HDC peaks vs width-matched shuffled regions.
- **b** Enhancer-state enrichment (Roadmap 15-state EnhG + Enh) of cell-type-specific peaks across matched epigenomes.
- **c** Top HOMER known motifs in cell-type-specific peaks vs other HDC peaks.
- **d** Positive scMultiMap links per cell type (FDR 10%); title gives links in within-donor permuted data.
- **e** Fraction of linked vs unlinked pairs with promoter-capture Hi-C contact (H1 and H1-derived cells); Mantel–Haenszel OR above bars.
- **f** At caQTL lead variants of linked peaks: nominal eQTL of the linked gene and allelic direction concordance.
- **g** Recovery of per-cell-type links by the pooled analysis, by number of cell types in which a link is significant.
- **h** GAP43 locus: pseudobulk accessibility per cell type and 19 neuronal links (solid arcs, PCHi-C contact).

### Figure 2. Developmental QTLs and regulatory elements in adult traits.

- **a** Genome-wide π₁ of permutation p-values per caQTL and eQTL map.
- **b** mashR QTLs per cell type (lfsr < 0.05); totals across cell types at top right.
- **c** Replication in adult brain neuron and glia caQTLs: π₁ and sign concordance for mashR caQTLs vs null peaks; one dot per HDC cell type, bars are medians.
- **d** S-LDSC of cell-type-specific peaks conditioned on baselineLD and all HDC peaks: −log10 p(τ) for traits with any FDR < 10% result (dots).
- **e** Independent colocalized regions (PP.H4 > 0.5) found by per-cell-type maps, the pooled map, or both.
- **f** 11q13: angina GWAS (top) and neural crest caQTL (bottom) for peak chr11:65,586,716–65,587,229 (shaded). The diamond marks shared lead variant rs6591182. A second, broader GWAS signal around 65.7 Mb is not shared with the caQTL.
- **g** GTEx v8 eGene status of genes with an HDC eQTL colocalization vs genes at non-colocalized GWAS leads.

### Figure 3. TF regulatory networks of the HDC (SPIRE).

- **a** Edge Jaccard of each paired network with an independent replicate (new pseudocells and seed; grey) and with the unpaired network (RNA and ATAC from different cells of the same state; orange). Percentages are unpaired / paired edge counts.
- **b** Cognate HOMER motif enrichment among positive targets per TF in the joint network, TFs with ≥ 50 targets; pink = binomial FDR < 5%.
- **c** Fraction of caQTL peaks whose allelic direction is predicted by motif disruption (motif Δ × edge sign), for SPIRE-assigned vs other expressed TFs; Fisher p above bars.
- **d** Top 3 key TFs per cell type (joint network): colour = expression z across cell types; circle = accessibility z of the TF's targets.
- **e** Spearman correlation of TF regulatory scope between within-cell-type networks (Vascular_endothelial excluded); replicate ceilings in title.
- **f** Drivers along the neuronal and neural crest lineages, ordered by expression timing: expression z (left) and target accessibility z (right) over 2 iPSC + 8 pseudotime bins.
- **g** Drivers used by two lineages: target Jaccard vs mean expression timing; orange = early (< 0.35), purple = late.

## Supplementary figures

### Figure S1. Genomic and epigenomic features of HDC peaks.

- **a** Genomic context of the 524,022 HDC peaks (promoter, exon, intron, distal) vs width- and chromosome-matched shuffled regions.
- **b** Fraction of HDC peaks in active Roadmap 18-state chromatin per epigenome (dots), grouped by tissue; boxes show the median and quartiles.
- **c** Peak features by cell-type specificity decile (τ index; 10 = most specific): promoter overlap, active or bivalent/Polycomb state in ESC (H1/H9), and phastCons element overlap.
- **d** HOMER known-motif enrichment in the top 3,000 CT-specific peaks of each cell type vs 30,000 other HDC peaks (200 bp, GC-normalized). Colour is −log10 p, capped at 200. Full version of Fig. 1c.

### Figure S2. Properties of peak–gene links.

- **a** Positive scMultiMap links per 1,000 tested peak–TSS pairs, by distance bin and cell type.
- **b** Enrichment of annotations among linked vs unlinked tested peaks. Odds ratios come from logistic regression adjusted for accessibility, GC content, width, TSS distance and number of candidate genes; bars are 95% CIs.
- **c** Link sharing across cell types: π₁ of the discovery cell type's links (rows) among p-values in each replication cell type (columns) and in the pooled analysis (AllCells).
- **d** Fraction of linked and tested genes and peaks that are cell-type specific.
- **e** S-LDSC of linked CT-specific peaks conditioned on baselineLD, all HDC peaks and the cell type's specific peaks. Bars show the fraction of primary traits with τ p < 0.05 (dashed line = 5% expected); labels give the number of traits at FDR < 10% (0 in every cell type).

### Figure S3. Enhancer-state enrichment of CT-specific peaks across all Roadmap epigenomes.

- **a** Roadmap 18-state model (98 epigenomes): log2 enrichment of enhancer states (EnhG, EnhA, EnhWk) at each cell type's specific peaks relative to all HDC peaks.
- **b** Roadmap 15-state model (127 epigenomes): the same for EnhG + Enh.

The colour strip under each heatmap marks the epigenome's tissue group. Fig. 1b shows the matched epigenomes from panel b.

### Figure S4. Calibration of scMultiMap detection thresholds by within-donor permutation.

- **a, b** iPSC (a) and Neuronal (b). Each cell is one combination of gene and peak detection thresholds (minimum fraction of the cell type's cells). Colour and the top number give the ratio of permuted to real links at FDR 10%; the bottom number is the count of real links. The box marks the thresholds used (gene ≥ 10%, peak ≥ 2%), where the ratio is 0.05 (iPSC) and 0.01 (Neuronal).
- **c** Real vs permuted link counts across all threshold combinations (circles iPSC, squares Neuronal; colour = gene threshold). The dashed line is y = x.

### Figure S5. Additional peak–gene link loci.

- **a** TDGF1 (iPSC): 16 positive links.
- **b** GATA6: 13 positive links (Endoderm 11, Mesoderm 2).

Tracks show pseudobulk accessibility per cell type. Arcs connect linked peaks (triangles) to the gene TSS: solid arcs have promoter-capture Hi-C contact, dashed arcs do not. Arc colour gives the cell type in which the link was found.

### Figure S6. Covariate optimization for QTL mapping.

π₁ of permutation p-values is computed on four held-out chromosomes for each cell type × modality.

- **a, b** Number of phenotype PCs with no other covariates, for RNA (a) and ATAC (b). Dashed lines mark the selected values (RNA 2, ATAC 3).
- **c** Change in π₁ from adding sex and/or genotype PCs (gPC), relative to the best phenotype-PC-only model of each set.
- **d** Cis-window size, relative to the 1 Mb window, with no MAF filter.
- **e** MAF filter, relative to no filter, with a 1 Mb window.

In d and e each line is one set at its own covariate configuration (blue RNA, red ATAC). Ectoderm_2 is absent because it was not QTL-mapped: only 9 donors had ≥ 20 cells.

### Figure S7. Sharing of QTL effects across cell types (mashR).

Pairwise sharing of caQTLs (**a**) and eQTLs (**b**): among features significant (lfsr < 0.05) in either cell type, the fraction whose posterior effects have the same sign and are within 2-fold of each other.

### Figure S8. Replication of HDC caQTLs in adult brain, per cell type.

Each HDC caQTL is matched to adult neuron or glia caQTL summary statistics at its best shared SNP.

- **a, b** π₁ among mashR caQTLs (purple) and mashR null peaks (lfsr > 0.5 in all cell types; grey), in neurons (a) and glia (b). Bars are 95% bootstrap intervals (100 bootstraps).
- **c, d** Sign concordance of adult effects at p < 0.05 (dashed line = 0.5 expected by chance).
- **e, f** Odds ratio that an HDC caQTL peak overlaps an adult brain caQTL peak, for mashR caQTLs (purple) and permutation p < 0.01 caQTLs (green). Odds ratios come from logistic regression adjusted for accessibility (log CPM) and peak width; bars are 95% CIs.

CT-specific mashR caQTLs (1–27 per cell type) are too few for π₁ and are not shown.

### Figure S9. Full S-LDSC results for peak, gene, QTL and link annotations.

Primary traits have h² z > 5 (109 per annotation). Each annotation is fitted with baselineLD v2.2 plus a control:
- CT-specific peaks: all HDC peaks
- CT-specific genes: all HDC genes
- caQTL peaks: tested peaks
- eGenes: tested genes
- linked peaks: all HDC peaks
- linked CT-specific peaks: all HDC peaks + that cell type's specific peaks

Panels:
- **a** Fraction of traits with τ p < 0.05 per annotation (dashed line = 5% expected). Numbers above bars are traits at FDR < 10%.
- **b, c** τ z for the 25 traits with the lowest p across CT-specific peak and gene annotations (b), and across caQTL peak and eGene annotations (c). Where a trait is present in several cohorts, the largest τ z is shown. * = FDR < 10%.

Only CT-specific peak and gene annotations reach FDR < 10%. QTL and link annotations stay at or below the 5% nominal rate.

### Figure S10. Audit of the 155 colocalizations (PP.H4 > 0.5).

- **a** PP.H4 from this study's implementation vs R coloc 5.2.3 (coloc.abf, same priors). 153 of 155 agree to four decimals; the largest difference is 0.009.
- **b** Sensitivity to the colocalization prior: PP.H4 at p12 = 1e-5 (used) vs 5e-6; 59 of 155 stay above 0.5 (purple per-cell-type maps, orange pooled map).
- **c** Feature-level QTL permutation p (beta-approximated) at colocalized features; 100 of 155 have p < 0.05 (dashed line).
- **d** Colocalizations per QTL map and modality.

The full table, with the R check, the p12 sensitivity and the permutation p of each colocalization, is `GWAS/results/coloc_audit/hits_audited.tsv`.

### Figure S11. SPIRE drivers along all five differentiation lineages.

Neuronal (**a**), neural crest (**b**), retinal (**c**), mesoderm (**d**) and endoderm (**e**) slingshot lineages. Each panel shows the top 20 drivers by number of targets, ordered by expression timing (centre of mass). Left: TF expression z. Right: accessibility z of the TF's targets. Columns are 2 iPSC bins followed by 8 equal-count pseudotime bins.

Drivers have ≥ 50 targets, an expression range ≥ 1 log2 unit, and r(expression, target activity) ≥ 0.5. Panels a and b are shown in Fig. 3f.

### Figure S12. SPIRE network quality control.

- **a** Network sparsity: stable TFs per tested peak in each design:
  - all: joint network over 11 cell types
  - nn_: within-cell-type MultiVI neighbourhood pseudocells
  - traj_: lineage pseudotime bins
  - ct_iPSC: random pooling negative control
- **b** Consensus edges with paired vs unpaired pseudocells (RNA and ATAC from different cells of the same state).
- **c** Reproducibility: Spearman correlation of TF regulatory scope with an independent replicate (new pseudocells and seed).
- **d** Cognate-motif enrichment among positive edges in each network (bars), against 100 TF-label permutations (mean ± 2 s.d.).

The Vascular_endothelial within-cell-type network is dense (1.3 TFs per peak) with no motif enrichment (1.05×), so it is excluded from cross-cell-type comparisons.

### Figure S13. S-LDSC of TF programs and lineage-driver targets.

TF programs are the positive targets of each cell type's top 10 key TFs, fitted on baselineLD + all HDC peaks + that cell type's specific peaks. Lineage-driver sets ("dyn") are the targets of each lineage's top 10 drivers, fitted on baselineLD + all HDC peaks.

- **a** τ z for the 12 traits with the lowest p across the 16 sets. No result reaches FDR < 10% (lowest FDR 0.73).
- **b** Number of the 109 primary traits with τ p < 0.05 per set; the dashed line is the 5% expectation. In red: the only set with an excess at binomial p < 0.05, the endothelial program (11 of 109, p = 0.02, led by angina).

Over all sets, 67 of 1,744 tests (3.8%) are nominal. The per-set summary is Table 1 and `TFNet/results/sldsc_tf_summary.tsv`.
