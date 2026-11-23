"""
Annotate every HDC ATAC peak (524k; PseudoBulk/features_ATAC.tsv.gz, hg38) for the peak-characterization figure.

References (PeakAnnot/ref/, downloaded by PeakAnnot/ref/download.sh):
  GENCODE v32 basic (= cellranger-arc 2020-A), ENCODE SCREEN cCRE v4 (GRCh38), Roadmap 18-state ChromHMM
  (core_K27ac joint model, hg38 liftover, 98 epigenomes), UCSC phastConsElements100way, hg38.fa ($GENOME_FA).

Per peak:
  width, GC fraction, genomic category (promoter = within 1 kb of a protein_coding/lncRNA transcript TSS >
  exon > intron > distal), nearest gene TSS (protein_coding/lncRNA) and signed distance,
  cCRE v4 class with maximal bp overlap (none = not in ENCODE registry), phastCons element bp overlap,
  pseudobulk CPM per cell type (all donors summed), tau specificity (log2 CPM+1, 10 annotated cell types),
  top cell type, CT-specific flags (top 10% specificity; same definition as 11 / S-LDSC peakspec).
  ChromHMM state at the peak center in each of 98 epigenomes (separate int8 matrix).
Genome-wide bp coverage of each ChromHMM state and cCRE class (for fold enrichment = peak share / genome share).

--shuffle: the same annotation for a null set (bedtools shuffle of every peak within its chromosome, seed 1, avoiding
UCSC gaps + ENCODE blacklist v2; no cell-type columns) -> *_shuffled outputs.

Output: PeakAnnot/peak_annotation.parquet, PeakAnnot/peak_chromhmm18.parquet, PeakAnnot/chromhmm18_genome_bp.tsv,
        PeakAnnot/ccre_genome_bp.tsv
"""
import glob
import os
import subprocess
import sys

import numpy as np
import pandas as pd

BASE = os.environ["HDC_DIR"] + "/"
REF = f"{BASE}/PeakAnnot/ref/"
OUT = f"{BASE}/PeakAnnot/"
TMP = f"{OUT}/tmp/"
BT = os.environ["TOOLS_BIN"] + "/bedtools"
FA = os.environ["GENOME_FA"]
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
SHUF = "--shuffle" in sys.argv
SUF = "_shuffled" if SHUF else ""
TMP = f"{TMP}/shuffled/" if SHUF else TMP
os.makedirs(TMP, exist_ok=True)


def bt(args, out):
    with open(out, "w") as fh:
        subprocess.run([BT] + args, stdout=fh, check=True)


def write_bed(df, path, cols=("chr", "start", "end", "peak")):
    df[list(cols)].sort_values(["chr", "start"]).to_csv(path, sep="\t", header=False, index=False)


pk = pd.read_csv(f"{BASE}/PseudoBulk/features_ATAC.tsv.gz", sep="\t").rename(columns={"id": "peak"})
pk = pk[pk.chr.isin(CHROMS)].reset_index(drop=True)
pk["width"] = pk.end - pk.start
pk["center"] = (pk.start + pk.end) // 2
if SHUF:
    fai = pd.read_csv(f"{FA}.fai", sep="\t", header=None, usecols=[0, 1])
    fai[fai[0].isin(CHROMS)].to_csv(f"{TMP}/genome.txt", sep="\t", header=False, index=False)
    gap = pd.read_csv(f"{REF}/hg38.gap.txt.gz", sep="\t", header=None, usecols=[1, 2, 3])
    bl = pd.read_csv(f"{REF}/hg38-blacklist.v2.bed.gz", sep="\t", header=None, usecols=[0, 1, 2])
    pd.concat([gap.set_axis([0, 1, 2], axis=1), bl]).to_csv(f"{TMP}/excl.bed", sep="\t", header=False, index=False)
    write_bed(pk, f"{TMP}/orig.bed")
    bt(["shuffle", "-i", f"{TMP}/orig.bed", "-g", f"{TMP}/genome.txt", "-chrom", "-seed", "1", "-noOverlapping",
        "-excl", f"{TMP}/excl.bed"], f"{TMP}/shuf.bed")
    pk = pd.read_csv(f"{TMP}/shuf.bed", sep="\t", header=None, names=["chr", "start", "end", "peak"])
    pk["peak"] = "shuf_" + pk.peak
    pk["width"] = pk.end - pk.start
    pk["center"] = (pk.start + pk.end) // 2
write_bed(pk, f"{TMP}/peaks.bed")
print(f"{len(pk)} peaks", flush=True)

# ---------------- gene annotation ----------------
g = pd.read_csv(f"{REF}/gencode.v32.basic.annotation.gtf.gz", sep="\t", comment="#", header=None,
                usecols=[0, 2, 3, 4, 6, 8], names=["chr", "type", "start", "end", "strand", "attr"])
g = g[g.chr.isin(CHROMS)]
a = g.attr.str.extract(r'gene_id "([^".]+)[^"]*";.*gene_type "([^"]+)";.*gene_name "([^"]+)"')
g = g.assign(gene_id=a[0], gene_type=a[1], gene_name=a[2]).drop(columns="attr")
g = g[g.gene_type.isin(["protein_coding", "lncRNA"])]
tx = g[g.type == "transcript"]
tss = np.where(tx.strand == "+", tx.start - 1, tx.end - 1)
prom = pd.DataFrame({"chr": tx.chr, "start": np.maximum(tss - 1000, 0), "end": tss + 1000, "peak": tx.gene_name})
write_bed(prom, f"{TMP}/prom.bed")
write_bed(g[g.type == "exon"].assign(start=lambda d: d.start - 1, peak=lambda d: d.gene_name), f"{TMP}/exon.bed")
write_bed(g[g.type == "gene"].assign(start=lambda d: d.start - 1, peak=lambda d: d.gene_name), f"{TMP}/gene.bed")
cat = pd.Series("distal", index=pk.peak)
for name in ("gene", "exon", "prom"):  # later assignments win: promoter > exon > intron > distal
    bt(["intersect", "-u", "-a", f"{TMP}/peaks.bed", "-b", f"{TMP}/{name}.bed"], f"{TMP}/hit_{name}.bed")
    hit = pd.read_csv(f"{TMP}/hit_{name}.bed", sep="\t", header=None)[3]
    cat[hit] = {"gene": "intron", "exon": "exon", "prom": "promoter"}[name]
pk["category"] = cat.reindex(pk.peak).to_numpy()

gene_tss = g[g.type == "gene"].copy()
gene_tss["tss"] = np.where(gene_tss.strand == "+", gene_tss.start - 1, gene_tss.end - 1)
pk["nearest_gene"], pk["dist_tss"] = None, 0
for c, d in pk.groupby("chr"):
    t = gene_tss[gene_tss.chr == c].sort_values("tss")
    tp = t.tss.to_numpy()
    j = np.searchsorted(tp, d.center.to_numpy())
    jl, jr = np.clip(j - 1, 0, len(tp) - 1), np.clip(j, 0, len(tp) - 1)
    jj = np.where(np.abs(d.center.to_numpy() - tp[jl]) <= np.abs(tp[jr] - d.center.to_numpy()), jl, jr)
    sgn = np.where(t.strand.to_numpy()[jj] == "+", 1, -1)
    pk.loc[d.index, "nearest_gene"] = t.gene_name.to_numpy()[jj]
    pk.loc[d.index, "dist_tss"] = (d.center.to_numpy() - tp[jj]) * sgn
pk["dist_tss"] = pk.dist_tss.astype(int)
print("genes done", pk.category.value_counts().to_dict(), flush=True)

# ---------------- GC ----------------
bt(["nuc", "-fi", FA, "-bed", f"{TMP}/peaks.bed"], f"{TMP}/nuc.txt")
nuc = pd.read_csv(f"{TMP}/nuc.txt", sep="\t")
pk["gc"] = pk.peak.map(dict(zip(nuc.iloc[:, 3], nuc.iloc[:, 5])))

# ---------------- ENCODE cCRE v4 ----------------
cc = pd.read_csv(f"{REF}/GRCh38-cCREs.V4.bed", sep="\t", header=None, names=["chr", "start", "end", "d", "e", "cls"])
cc = cc[cc.chr.isin(CHROMS)]
cc[["chr", "start", "end", "cls"]].sort_values(["chr", "start"]).to_csv(f"{TMP}/ccre.bed", sep="\t", header=False, index=False)
bt(["intersect", "-wo", "-a", f"{TMP}/peaks.bed", "-b", f"{TMP}/ccre.bed"], f"{TMP}/ccre_ov.txt")
ov = pd.read_csv(f"{TMP}/ccre_ov.txt", sep="\t", header=None, usecols=[3, 7, 8], names=["peak", "cls", "bp"])
ov = ov.sort_values("bp", ascending=False).drop_duplicates("peak")
pk["ccre"] = pk.peak.map(dict(zip(ov.peak, ov.cls))).fillna("none")
if not SHUF:
    cc.assign(bp=cc.end - cc.start).groupby("cls").bp.sum().to_csv(f"{OUT}/ccre_genome_bp.tsv", sep="\t")
print("cCRE done", pk.ccre.value_counts().to_dict(), flush=True)

# ---------------- phastCons elements ----------------
ph = pd.read_csv(f"{REF}/phastConsElements100way.txt.gz", sep="\t", header=None, usecols=[1, 2, 3], names=["chr", "start", "end"])
ph = ph[ph.chr.isin(CHROMS)].sort_values(["chr", "start"])
ph.to_csv(f"{TMP}/phast.bed", sep="\t", header=False, index=False)
bt(["intersect", "-wao", "-a", f"{TMP}/peaks.bed", "-b", f"{TMP}/phast.bed"], f"{TMP}/phast_ov.txt")
po = pd.read_csv(f"{TMP}/phast_ov.txt", sep="\t", header=None, usecols=[3, 7], names=["peak", "bp"]).groupby("peak").bp.sum()
pk["phastcons_bp"] = pk.peak.map(po).fillna(0).astype(int)
pd.Series({"phastcons_bp": (ph.end - ph.start).sum()}).to_csv(f"{OUT}/phastcons_genome_bp.tsv", sep="\t", header=False)

# ---------------- cell-type accessibility / specificity ----------------
cts = [] if SHUF else sorted(f.split("/")[-1].replace("_ATAC_SumCount.tsv.gz", "") for f in glob.glob(f"{BASE}/PseudoBulk/*_ATAC_SumCount.tsv.gz")
             if not f.split("/")[-1].startswith("AllCells"))
cpm = {}
for ct in cts if not SHUF else []:
    x = pd.read_csv(f"{BASE}/PseudoBulk/{ct}_ATAC_SumCount.tsv.gz", sep="\t", index_col=0).sum(axis=1)
    cpm[ct] = x / x.sum() * 1e6
if not SHUF:
    cpm = pd.DataFrame(cpm).fillna(0).reindex(pk.peak)
    for ct in cts:
        pk[f"cpm_{ct}"] = cpm[ct].to_numpy(np.float32)
    named = [c for c in cts if c != "Unknown"]
    L = np.log2(cpm[named].to_numpy() + 1)
    mx = L.max(1)
    pk["tau"] = np.where(mx > 0, (1 - L / np.where(mx > 0, mx, 1)[:, None]).sum(1) / (len(named) - 1), np.nan)
    pk["top_ct"] = np.array(named)[L.argmax(1)]
    spec = pd.read_csv(f"{BASE}/GWAS/celltype_specific_peaks.tsv.gz", sep="\t").set_index("peak")
    for ct in cts:
        pk[f"spec_{ct}"] = spec[ct].reindex(pk.peak).fillna(False).to_numpy(bool)

# ---------------- ChromHMM 18-state, 98 epigenomes ----------------
files = sorted(glob.glob(f"{REF}/chromhmm18/E*_18_core_K27ac_hg38lift_dense.bed.gz"))
states = {}
gbp = {}
for f in files:
    eid = os.path.basename(f).split("_")[0]
    s = pd.read_csv(f, sep="\t", header=None, skiprows=1, usecols=[0, 1, 2, 3], names=["chr", "start", "end", "state"])
    s = s[s.chr.isin(CHROMS)].sort_values(["chr", "start"])
    s["code"] = s.state.str.split("_").str[0].astype(np.int8)
    gbp[eid] = (s.end - s.start).groupby(s.code).sum()
    out = np.zeros(len(pk), np.int8)
    for c, d in pk.groupby("chr"):
        sc = s[s.chr == c]
        j = np.searchsorted(sc.start.to_numpy(), d.center.to_numpy(), side="right") - 1
        ok = (j >= 0) & (d.center.to_numpy() < sc.end.to_numpy()[np.maximum(j, 0)])
        out[d.index[ok]] = sc.code.to_numpy()[j[ok]]
    states[eid] = out
    print(eid, flush=True)
st = pd.DataFrame(states)
st.insert(0, "peak", pk.peak.to_numpy())
st.to_parquet(f"{OUT}/peak_chromhmm18{SUF}.parquet", index=False)
if not SHUF:
    pd.DataFrame(gbp).fillna(0).astype(np.int64).to_csv(f"{OUT}/chromhmm18_genome_bp.tsv", sep="\t")
pk.to_parquet(f"{OUT}/peak_annotation{SUF}.parquet", index=False)
print("done", flush=True)
