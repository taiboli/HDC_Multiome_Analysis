# HDC Multiome Analysis

Analysis code for the paper on single-cell multiome (10x RNA + ATAC) profiling of human iPSC-derived heterogeneous differentiating cultures (HDCs). The cultures come from 17 Yoruba (YRI) iPSC lines: 128,588 paired cells, 36,601 genes and 524,022 peaks across 11 cell types.

What the code covers:
- cell-type-specific caQTL and eQTL mapping (tensorQTL, mashR)
- colocalization and S-LDSC with 187 GWAS
- replication in adult brain caQTLs
- peak annotation and peak–gene links (scMultiMap)
- TF regulatory networks (SPIRE / orthogonal matching pursuit)
- all main and supplementary figures

## Layout

| Path | Contents |
|---|---|
| `scripts/` | The numbered pipeline (`01_` … `34_`), SLURM submit scripts and small config tables (`celltypes.txt`, `tasks.tsv`, `mash_tasks.tsv`). |
| `scripts/README.md` | Step-by-step description: inputs, filters, parameters, outputs and key results of every script. |
| `upstream/HDC_MultiVI_Cleanup.ipynb` | Adds cell-type labels to the MultiVI-integrated object that the pipeline starts from. Outputs are stripped. |
| `ref/download_references.sh` | Downloads the public annotation references (Roadmap ChromHMM 18/15-state, ENCODE cCRE v4, phastCons, GENCODE v32). |
| `config/paths.example.env` | Template for the input, output and software paths (environment variables) used by every script. |
| `Figure_legends.md` | Legends for Figures 1–3 and S1–S13. |
| `LICENSE` | MIT license. |

## Pipeline overview

The numbers are the script prefixes; `scripts/README.md` gives the details.

| Steps | Analysis |
|---|---|
| 01–05 | Pseudobulk per cell type × donor, tensorQTL input, cis mapping, q-values, covariate optimization |
| 06–12 | GWAS harmonization, coloc / SMR, coloc audit, QTL-lead enrichment, multiome-link and GTEx comparisons |
| 13–14, 19 | mashR across cell types and trait enrichment of mashR QTLs |
| 15–18 | S-LDSC annotations, LD scores, regression and summaries |
| 20–21 | Replication in adult brain neuron and glia caQTLs |
| 22–26 | Cell-type matrices, scMultiMap peak–gene links (NumPy reimplementation, validated against R), peak annotation, HOMER, ChromHMM, PCHi-C and GTEx validation |
| 27–28 | Figures 1–2 |
| 29–33 | SPIRE TF networks, motif scans, network analysis, caQTL motif test, Figure 3, S-LDSC of TF programs |
| 34 | Supplementary figures S1–S13 |

## Running

All file locations are read from environment variables, so no paths are hard-coded. To set them:

1. Copy `config/paths.example.env` to `config/paths.local.env`.
2. Fill in your own paths. The `.local.env` file is ignored by git.
3. Load the file into your shell before running anything:

```bash
set -a; source config/paths.local.env; set +a
python scripts/27_figures.py           # run a script directly
sbatch scripts/submit_03_map_cis_gpu.sh  # or submit a SLURM job
```

`HDC_DIR` is the analysis directory: it holds the input objects and receives all outputs. The other variables point to external inputs (genotypes, GWAS, GTEx, brain caQTLs, LD-score references, pseudotime, TF list, genome) and to software environments.

SLURM jobs inherit the variables from the shell they were submitted from. Submit them from a directory that contains a `logs/` folder. The `#SBATCH --account` and `--partition` lines are set for the original cluster; change them for yours.

Large inputs are not in this repository: the h5ad objects, genotypes, GWAS summary statistics, LD scores and QTL outputs.

## Dependencies

**Python 3.11:**
- numpy, pandas, scipy
- pyarrow, h5py
- matplotlib, seaborn
- tensorqtl 1.0.10 with torch (GPU)
- anndata / scanpy / scvi-tools (upstream notebook)

**R 4.2:**
- data.table, dplyr
- edgeR
- mashr
- coloc 5.2.3
- GenomicRanges, rtracklayer
- arrow

**External tools:**
- ldsc (Python 2.7 environment) with baselineLD v2.2
- HOMER (`findMotifsGenome.pl`, known vertebrate motifs)
- ghostscript (`34_supp_figures.py`)

**Third-party methods used but not included:**
- **SPIRE**: `29_spire.py` imports the package; `spire_omp_fast.py` is a vectorized OMP that gives identical supports and coefficients.
- **scMultiMap** (Su et al. 2025): `scmultimap_py.py` is a NumPy reimplementation. `23v_validate_scmultimap.R` compares it with the original R code, which must be installed separately.

## Not included

The MultiVI model training and cell-type annotation that produced the input object are not in this repository; only the labelling notebook in `upstream/` is.
