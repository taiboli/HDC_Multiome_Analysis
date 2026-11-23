# HDC per-cell-type eQTL / caQTL (tensorQTL)

Reuses the EB_CellTypeData pipeline (EBQTL/RNA/EB_RNA_tensorQTL_TMM_MakeInput.R,
EB_CellTypeData/EB_Celltypedata_Assign.R, {eQTL,caQTL}/submit_*_PerCellType_TMM.sh) on the
filtered HDC multiome object.

Input: `MultiVI/HDC_multiome_RemoveHighRibo_posttraining_adata_germlayer.h5ad`
(raw counts, 36,601 genes + 524,022 peaks); cells = labeled `paired` cells (128,588),
cell type = `highres_celltype`, donor = `donor` (17 YRI lines).

| Step | Script | SLURM | What |
|---|---|---|---|
| 1 | 01_pseudobulk.py | submit_01_pseudobulk.sh (gilad-hm) | one pass over h5ad: per celltype x donor summed counts + # cells with count>0 -> `PseudoBulk/` |
| 2 | 02_make_tensorqtl_input.R | submit_02_make_input.sh (gilad-hm, array 1-22 = tasks.tsv) | donor >=20 cells (skip celltype if <10 donors); feature in >=5 cells & >=1% cells in >=80% donors; TMM logCPM -> invnorm; cov = Sex + PhenoPCs + 5 GenoPCs -> `eQTL/<CT>/`, `caQTL/<CT>/` |
| 3 | 03_map_cis.py | ATAC: submit_03_map_cis_gpu.sh (gpu, array over ATAC rows of tasks.tsv); RNA: submit_03b_map_cis_gpu_RNA_all.sh (one gpu job, all cell types; gpu QOS max 12 jobs/user) | tensorQTL cis (10k permutations, seed 2025) + cis_nominal (all pairs, eigenMT) per chr1-22 |
| 4 | 04_merge_qvalues.py | submit_04_merge.sh | genome-wide q-values per celltype x modality, completeness check -> `QTL_summary{_optCov}.tsv` |
| 5 | 05_optimize_covariates.py | submit_05_optimize_gpu.sh (ATAC), submit_05_optimize_cpu_RNA.sh (RNA) | covariate / maf / window search on chr2,10,17,20 (1k perms) -> `CovOptimization/` (`all_summary.tsv`) |

Differences from the EB run:
- PhenoPCs: 5 when >= 16 donors (as before); fewer when fewer donors, so tensorQTL residual dof >= 3.
- Permutation q-values are computed genome-wide (step 4) instead of per chromosome.
- BED sorted by chr, start, end (tensorQTL 1.0.10 requires end order for tied starts).
- Variants monomorphic in the included donors are dropped before mapping (applied from Neuronal ATAC chr2 onward; they only produce NaN tests).
- Genotypes for a chromosome are loaded once for both permutation and nominal passes; GPU.

## Final covariate set: `optCov` (run steps 3-4 with `COVSET=optCov`)
Chosen by step 5 (maximize pi1; hit counts too small to rank): RNA = 2 PhenoPCs, ATAC = 3 PhenoPCs,
no Sex / GenoPCs (always lowered pi1), maf 0, +/-1 Mb. Files `tensorqtl.{gene,atac}Cov_optCov.txt`.
Genome-wide pi1 rose from ~0-0.08 (11-cov) to ~0.1-0.33 in 18/20 maps (not Vascular_endothelial, 12 donors).
Hits FDR 10%: Retinal RNA 4, Neural_crest ATAC 2 (11-cov run: 0 everywhere). Old 11-cov outputs deleted
2026-10-14 (inputs `tensorqtl.*Cov_noBootstrap_TMM_invnorm.txt` and old `QTL_summary.tsv` kept).

Outputs per cell type in `eQTL_Output_optCov/<CT>/` (RNA) and `caQTL_Output_optCov/<CT>/` (ATAC):
- `chr{N}_cis_noBootstrap_TMM_invnorm.cis_qtl.txt.gz` — permutation, lead variant per phenotype
- `chr{N}_cis_nominal_noBootstrap_TMM_invnorm.cis_qtl_pairs.{N}.parquet` — all cis pairs (for coloc / enrichment)
- `allchr_cis_noBootstrap_TMM_invnorm.cis_qtl.qval.txt.gz` — all chromosomes with genome-wide q-values
- per-run `.log` files; inputs + `input_summary.tsv` (donors, # features, # PCs) in `eQTL/<CT>/`, `caQTL/<CT>/`

Not run: Ectoderm_2 (only 9 donors with >= 20 cells). Vascular_endothelial: 12 donors, 1 PhenoPC (low power).

## Trait relevance: GWAS harmonization, coloc, SMR, enrichment (steps 06-12; outputs in `../GWAS/`)
GWAS resources reused from the EB projects (read-only): FinnGen R10 (77 traits, hg38), Pan-UKB 21 traits x
{meta_hq, EUR, AFR} (hg19), BBJ 22 traits (hg19), PAGE 28 traits (multi-ethnic, AFR-rich, hg38); GTEx v8 from
$GTEX_DIR.

| step | script (submit) | what | output |
|---|---|---|---|
| 06a | `06a_variant_universe.R` (`submit_06a_universe.sh`) | HDC genotype variants, hg38 + hg19 liftOver | `GWAS/variant_universe.parquet` |
| 06b | `06b_harmonize_gwas.py manifest` / `run <row>` (`submit_06b_harmonize.sh`, array 0-189) | align each GWAS to HDC ALT allele; QC (match rate, lambda, sdY) | `GWAS/gwas_manifest.tsv`, `GWAS/harmonized/<trait>.parquet`, `GWAS/qc/` |
| 07 | `07_coloc_smr.py <CT> <MOD> <chr>` (`submit_07_coloc.sh`, tasks `GWAS/coloc_tasks.tsv`, 440) | vectorized coloc.abf (p1=p2=1e-4, p12=1e-5; validated against R coloc 5.2.3 with `07v_validate_coloc.*`) for every feature whose window has GWAS p<1e-6; lead-variant table with GWAS z for all traits | `GWAS/coloc/<MOD>/<CT>/chr<N>.{coloc,lead}.parquet` |
| 08 | `08_summarize_coloc_smr.py` | collect coloc; SMR at GWS lead variants (QTL z from pval_beta, avoids winner's curse); BH per trait and global | `GWAS/results/coloc_all.tsv.gz`, `coloc_H4gt0.5.tsv`, `coloc_summary.tsv`, `smr_gwsLeads.tsv.gz` |
| 09 | `09_enrichment_qtl_leads.py` | locus-free: GWAS chi2 at strong (pval_beta<0.01) vs null (>0.5) QTL leads, MAF x distance standardized, block jackknife; pooled ALL_<MOD> | `enrichment_qtl_leads.tsv` |
| 10 | `10_multiome_link_gtex.py` | caQTL<->eQTL coloc at trait-coloc loci; GTEx eGene / eQTL annotation | `peak_gene_links.tsv`, `coloc_H4gt0.5_annotated.tsv` |
| 11 | `11_celltype_peak_enrichment.py` | donor-independent: GWAS chi2 in cell-type-specific peaks vs all HDC peaks (11 CTs incl. Ectoderm_2) | `GWAS/celltype_specific_peaks.tsv.gz`, `celltype_peak_enrichment.tsv` |
| 12 | `12_gtex_background.py` | GTEx overlap of coloc SNPs vs GWAS-lead SNPs at non-colocalized HDC features | `gtex_novelty_vs_background.tsv` |

`submit_08_11_post.sh` runs 08-11 sequentially. Caveats: n = 17 donors, so coloc PP.H4 is sensitive to p12
(27/78 PP.H4>0.5 survive p12=5e-6); SMR without HEIDI (no LD reference); chi2 enrichments do not model LD
(no S-LDSC yet); "Unknown" cell type = stress / integrated-stress-response cells (HSPA1A, DDIT3, ATF3, GDF15).

## Pooled map, mash, S-LDSC, brain replication (steps 01b, 05b, 13-21)

| step | script (submit) | what | output |
|---|---|---|---|
| 01b | `01b_pseudobulk_allcells.py` | pseudobulk over all cells per donor ("AllCells"), mapped with steps 02-04 like a cell type | `PseudoBulk/AllCells_*`, `{eQTL,caQTL}_Output_optCov/AllCells/` |
| 05b | `05b_write_optcov.py` | writes the chosen covariate files (`tensorqtl.*Cov_optCov.txt`; AllCells RNA Sex + GenoPC1 + 2 PhenoPCs, ATAC 2 PhenoPCs) | |
| 13-14 | `13_mash_extract.py` (`submit_13_mash_extract.sh`), `14_mash.R` (mashr env) | strong (per-feature best) and random subsets across cell types; mash fit on complete cases | `mash/<MOD>_mash_strong.tsv.gz` (lfsr / posterior per CT) |
| 15a/15/15b | `15a_liftover_features.R`, `15_ldsc_sets.py <spec\|qtl\|link>`, `15b_ldsc_annot.py` (`submit_15b_ldsc_l2.sh`) | hg19 feature sets -> binary annotations + LD scores (1000G Phase 3 EUR / EAS) | `GWAS/ldsc/sets_*.tsv.gz`, `GWAS/ldsc/ldscores/<POP>/<set>/` |
| 16-18 | `16_ldsc_sumstats.py`, `17_sldsc.sh` (array over traits; `BATCH=spec\|qtl\|link`), `18_sldsc_summarize.py` | S-LDSC: each set on top of baselineLD v2.2 + its control; tau*, one-sided p_tau, BH per family over primary traits (h2 z > 5, LD-matched) | `GWAS/ldsc/sldsc_<batch>.tsv` |
| 19 | `19_mash_trait_enrichment.py` | GWAS chi2 at mash-significant vs null features | `GWAS/results/mash_trait_enrichment.tsv` |
| 20 | `20_brain_caqtl_replication.py` | external replication in adult brain neuron / glia caQTLs (data/BraincaQTL; AMP-AD peak beds): pi1, sign concordance, peak-overlap enrichment | `BrainRep/` |
| 21 | `21_covopt_brainrep.py` | covariate grid scored by brain-replicated discoveries (selected on odd, evaluated on even chromosomes); no configuration robustly beat optCov -> no remapping | `BrainRep/covopt/` |

LD reference caveat: the Zenodo 10515792 release has EUR and EAS LD scores only. BBJ uses EAS, other GWAS EUR; PAGE and
UKB AFR (AFR-rich) are run with EUR LD scores and are exploratory (`ld_match = False`). Annotations are genomic
intervals and do not depend on the YRI donors' LD.

## Peak characterization and cell-type-specific peak-gene links (steps 22-27; outputs in `../PeakAnnot/`)

References (`PeakAnnot/ref/`, `ref/download_references.sh`): GENCODE v32 basic (= cellranger-arc 2020-A), ENCODE SCREEN cCRE v4,
Roadmap ChromHMM 18-state core_K27ac (98 epigenomes) and 15-state core (127) hg38 liftovers, UCSC phastConsElements100way,
UCSC gaps + ENCODE blacklist v2 (shuffle exclusion), hg38.fa ($GENOME_FA), promoter-capture
Hi-C (Jung et al. 2019 GSE86189, lifted hg19 -> hg38), GTEx v8 signif pairs ($GTEX_DIR). Tools env:
`$TOOLS_ENV` (bedtools, HOMER). Python: tensorqtl env (pyarrow); 22 uses scvi-gpu (h5py).

| step | script | what | output |
|---|---|---|---|
| 22 | `22_extract_ct_matrices.py` (bigmem, 3 min) | per-cell-type single-cell RNA / ATAC count matrices (same cells as 01) | `PeakAnnot/sc/<CT>/{RNA,ATAC}.npz, cells.tsv.gz` |
| 23 | `23_scmultimap.py <CT\|AllCells> [--perm]` (`submit_23_scmultimap.sh`, caslake array) + `scmultimap_py.py` | scMultiMap (Su et al. 2025) IRLS + WLS ported to NumPy (identical to the R package to ~1e-11, `23v_validate_scmultimap.{py,R}`); bsample = donor; pairs = peak center within 500 kb of TSS; `--perm` permutes RNA cells within donor | `PeakAnnot/scMultiMap/<CT>.pairs.parquet` (all tested pairs) |
| 24 | `24_peak_annotation.py [--shuffle]` | per-peak genomic context, GC, cCRE class, phastCons, ChromHMM state at center (98 epigenomes), CPM / tau / CT-specific flags; `--shuffle` = width/chromosome-matched null | `PeakAnnot/peak_annotation{,_shuffled}.parquet`, `peak_chromhmm18*.parquet` |
| 25 / 25c | `25_peak_characterization.py`, `25c_chromhmm15.py` | HDC vs shuffled; ChromHMM enrichment per epigenome; specificity deciles; CT-specific peaks x epigenome enhancer enrichment; ESC poised states; novel (not in registry) peaks | `PeakAnnot/results/*.tsv`, `PeakAnnot/figures/fig_peaks_*.pdf` |
| 25b | `25b_homer_prep.py`, `25b_homer.sh` (array) | HOMER known motifs: top 3000 CT-specific peaks vs 30000 other HDC peaks, 200 bp, GC-normalized | `PeakAnnot/homer/out/<set>/knownResults.txt` |
| 26a / 26b / 26c | `26a_gtex_peak_index.py`, `26b_pchic_liftover.R`, `26c_scmultimap_calibration.py` | GTEx eVariant-in-peak index; PCHi-C hg38; real vs permuted links over detection thresholds | `PeakAnnot/ref/gtex_peak_gene.parquet`, `ref/pchic/pchic_po_hg38.tsv.gz`, `results/scmultimap_calibration.tsv` |
| 26 | `26_link_analysis.py` | links = BH FDR < 10% after the calibrated filter (gene detected in >= 10%, peak in >= 2% of the cell type's cells, both >= 100 cells); summaries, distance decay, linked-peak annotation enrichment (logistic, adjusted for accessibility, GC, width, TSS distance, # candidate genes), PCHi-C / GTEx support (Mantel-Haenszel over distance bins), caQTL->eQTL concordance, cross-CT sharing, AllCells recovery, CT-specificity | `PeakAnnot/results/links_fdr0.1.tsv.gz`, `link_*.tsv` |
| 27 / 27b | `27_figures.py`, `27b_locus.py GENE ...` (`COMPACT=1` for a wide version; `draw_locus()` importable) | main Figure 1 (cCRE, matched-epigenome enhancers, motifs, links + permuted null, PCHi-C, caQTL->eQTL, pooled recovery, GAP43 locus) + supplement (other peak / link panels, S-LDSC of linked peaks) + per-panel PDFs; locus views | `Figures/Main/Figure1.{pdf,png}`, `PeakAnnot/figures/FigureS_peaks_links.pdf`, `figures/panels/` |
| 08c | `08c_coloc_audit.py` (caslake, 48G), then `07v_validate_coloc.R GWAS/results/coloc_audit/R` (R/4.2.0 + coloc 5.2.3) | audit of coloc: integrity, loci at trait / phenotype / region level (pooled vs per-CT), p12 = 5e-6 sensitivity, QTL permutation p at hits, export of all PP.H4 > 0.5 pairs for R coloc. 2026-10-16: 155 pairs, 62 regions (34 per-CT only, 24 pooled only, 4 both); all 155 reproduce in R (max dPP.H4 0.009); 59/155 keep PP.H4 > 0.5 at p12 = 5e-6; 55/155 features have permutation p >= 0.05 | `GWAS/results/coloc_audit/` |
| 28 | `28_figure2.py` | main Figure 2: QTL pi1, mashR counts, brain replication, S-LDSC of CT-specific peaks, colocalized regions, 11q13 neural crest caQTL x angina locus, GTEx novelty (12 rerun 2026-10-16 on all 155 hits: 80% vs 92% eGenes, p = 0.045) | `Figures/Main/Figure2.{pdf,png}`, `QTL_pi1_optCov.tsv` |
| 29 | `29_spire.py <run> [B]` (`submit_29_spire.sh`, env `RUNS`, `B`; caslake 96G) + `spire_omp_fast.py` | SPIRE TF -> peak networks (`$SPIRE_DIR`): per-peak OMP (L = 10) of TMM log-CPM peak accessibility on centred, unit-norm TF expression, bootstrap stability (B = 50, tau = 0.5) with pseudocells rebuilt every replicate (`resample="cells"`). Designs: `all` (11 CT x 60 pseudocells x 50 cells), `nn_<CT>` (300 seeds + 49 nearest neighbours in the MultiVI latent `TFNet/multivi_latent.tsv.gz`), `traj_<lin>` (slingshot pseudotime; 2 iPSC + 8 pseudotime bins x 40 pseudocells x 30 cells), `_unpaired` (RNA from different cells of the same state), `_s2` (seed 2027), `ct_<CT>` (random pooling, negative control). `omp_fast` = identical supports / coefficients to `spire.omp.omp`, ~3x faster | `TFNet/runs/<run>/` (edges, scope, freq, strata, log) |
| 30 / 30b | `30_motif_scan.sh` (array 0-23), `30b_motif_matrix.py` | HOMER known motifs (central 200 bp) in all 524,022 peaks; motif -> TF symbol map (425 of 471 motifs, 382 TFs) | `TFNet/motif/peak_tf_motif.npz`, `tf_motif.tsv` |
| 31 / 31b | `31_tfnet_analysis.py`, `31b_caqtl_motif.py` (caslake 64G) | motif concordance of edges (TF-label permutation null), key TFs (expression z x target activity), cross-CT scope / edge concordance, trajectory drivers and timing, paired vs unpaired, S-LDSC sets (`tfprog_<CT>`, `tfdyn_<lin>`, `BATCH=tf` in 15b / 17 / 18); caQTL motif-disruption x SPIRE edge sign concordance, TFs at colocalized peaks | `TFNet/results/` |
| 32 | `32_figure3.py` | main Figure 3 (TF networks, panels a-g) + supplement (S11 lineage cascades, S12 QC, S13 S-LDSC) | `Figures/Main/Figure3.{pdf,png}`, `PeakAnnot/figures/FigureS_tfnet*.{pdf,png}` |
| 33 | `33_sldsc_tf_summary.py` (after `BATCH=tf` 15b / 17 / 18) | per TF set: traits with tau p < 0.05 vs 5% (binomial), mean tau z, median enrichment, best trait / FDR. 2026-10-16: no set at FDR 10% (min 0.73); 67 / 1,744 tests p < 0.05 | `TFNet/results/sldsc_tf_summary.tsv` |
| 34 | `34_supp_figures.py` (run after 27, 27b, 28, 32; needs `gs`) | supplementary S1-S13: draws S4 (scMultiMap calibration), S6 (covariates), S7 (mashR sharing), S8 (brain replication per CT), S9 (full S-LDSC), S10 (coloc audit) from tables; assembles S1-S2, S3, S5, S11-S13 from earlier outputs. Legends: `Figures/Figure_legends.md` | `Figures/Supplementary/FigureS<n>.{pdf,png}` |

Calibration note: with the raw scMultiMap test, sparse features (detected in a few % of cells) give a positive
heavy-tailed null (within-donor permutation: iPSC 8,399 vs 15,264 real FDR-10% pairs). The bulk of the null is calibrated
(lambda ~0.95); after the detection filter the permuted / real ratio is 1.5% (Neuronal) and 4.8% (iPSC).
Vascular_endothelial (919 cells) has little power (44 links). AllCells was run with the filter applied up front
(`MINFRAC_GENE=0.10 MINFRAC_PEAK=0.02 MIN_CELLS=100`) to keep IRLS tractable on 128k cells.
