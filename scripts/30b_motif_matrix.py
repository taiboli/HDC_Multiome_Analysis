"""
Peak x TF cognate-motif matrix from the HOMER scan (30_motif_scan.sh).

Each HOMER known vertebrate motif is assigned to TF symbols by its name (text before "(", case-insensitive, "-" removed,
plus a small alias table for names that are not HGNC symbols, e.g. Oct4 -> POU5F1). Composite motifs ("A:B") are assigned
to every named partner. A peak "has the TF's motif" if any of the TF's motifs has an instance in the central 200 bp.
Output TFNet/motif/: peak_motif.npz (sparse peaks x motifs, bool; row order = PeakAnnot/sc/features_ATAC.tsv.gz),
motifs.tsv (motif, symbols), peak_tf_motif.npz (sparse peaks x TFs with >= 1 motif), tf_motif.tsv
"""
import os
import glob
import re

import numpy as np
import pandas as pd
import scipy.sparse as sp

BASE = os.environ["HDC_DIR"] + "/"
D = f"{BASE}/TFNet/motif/"
ALIAS = {"OCT4": ["POU5F1"], "OCT2": ["POU2F2"], "OCT6": ["POU3F1"], "OCT11": ["POU2F3"], "BRN1": ["POU3F3"],
         "BRN2": ["POU3F2"], "OCT": ["POU5F1", "POU2F1", "POU2F2", "POU3F1", "POU3F2", "POU3F3"],
         "AP2GAMMA": ["TFAP2C"], "AP2ALPHA": ["TFAP2A"], "AP4": ["TFAP4"], "AP1": ["FOS", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND"],
         "P53": ["TP53"], "P63": ["TP63"], "P73": ["TP73"], "NF1": ["NFIA", "NFIB", "NFIC", "NFIX"],
         "NFKBP65": ["RELA"], "NFKBP50": ["NFKB1"], "NFKBP65REL": ["RELA", "REL"], "EKLF": ["KLF1"], "EWS": [],
         "PU.1": ["SPI1"], "SPIB": ["SPIB"], "EBF": ["EBF1"], "MYOD": ["MYOD1"], "NGN2": ["NEUROG2"], "NEUROD1": ["NEUROD1"],
         "OLIG2": ["OLIG2"], "TCF4": ["TCF7L2"], "LEF1": ["LEF1"], "TCFL2": ["TCF7L2"], "TCF7": ["TCF7"], "TCF3": ["TCF7L1"],
         "HNF4A": ["HNF4A"], "HNF1B": ["HNF1B"], "HNF6": ["ONECUT1"], "COUPTFII": ["NR2F2"], "EAR2": ["NR2F6"],
         "ERRA": ["ESRRA"], "ERRG": ["ESRRG"], "ERE": ["ESR1"], "GRE": ["NR3C1"], "PGR": ["PGR"], "ARE": ["AR"],
         "RAR:RXR": ["RARA", "RXRA"], "THRB": ["THRB"], "TR4": ["NR2C2"], "NUR77": ["NR4A1"], "SF1": ["NR5A1"],
         "LRH1": ["NR5A2"], "RORGT": ["RORC"], "REVERBA": ["NR1D1"], "PPARE": ["PPARG"], "TLX?": ["NR2E1"],
         "CHOP": ["DDIT3"], "ATF1": ["ATF1"], "HRE": ["HSF1"], "HIF2A": ["EPAS1"], "HIF1A": ["HIF1A"], "HIF1B": ["ARNT"],
         "CLOCK": ["CLOCK"], "NPAS": ["NPAS2"], "BMAL1": ["ARNTL"], "CMYC": ["MYC"], "NMYC": ["MYCN"], "MNT": ["MNT"],
         "USF1": ["USF1"], "USF2": ["USF2"], "TFE3": ["TFE3"], "TFEB": ["TFEB"], "MITF": ["MITF"], "SREBP1A": ["SREBF1"],
         "SREBP2": ["SREBF2"], "ZNF143|STAF": ["ZNF143"], "STAT4": ["STAT4"], "CEBP:AP1": ["CEBPB"], "CEBP": ["CEBPA", "CEBPB"],
         "NFE2L2": ["NFE2L2"], "NRF2": ["NFE2L2"], "NRF1": ["NRF1"], "NRF": ["NRF1"], "GABPA": ["GABPA"], "YY1": ["YY1"],
         "CTCF": ["CTCF"], "BORIS": ["CTCFL"], "REST-NRSF": ["REST"], "TEAD": ["TEAD1", "TEAD2", "TEAD3", "TEAD4"],
         "ETS": [], "E2A": ["TCF3"], "HEB": ["TCF12"], "PTF1A": ["PTF1A"], "ZFX": ["ZFX"], "ZIC": ["ZIC1", "ZIC2", "ZIC3"],
         "SIX1": ["SIX1"], "SIX2": ["SIX2"], "PAX3:FKHR": ["PAX3"], "PAX7": ["PAX7"], "PAX5": ["PAX5"], "PAX6": ["PAX6"],
         "PAX8": ["PAX8"], "NKX2.1": ["NKX2-1"], "NKX2.2": ["NKX2-2"], "NKX2.5": ["NKX2-5"], "NKX3.1": ["NKX3-1"],
         "NKX6.1": ["NKX6-1"], "HOXA9": ["HOXA9"], "PIT1": ["POU1F1"], "LHX1": ["LHX1"], "ISL1": ["ISL1"],
         "SMAD2": ["SMAD2"], "SMAD3": ["SMAD3"], "SMAD4": ["SMAD4"], "ZNF41": ["ZNF41"], "GFY": [], "GFX": [],
         "RFX": ["RFX1", "RFX2", "RFX3"], "X-BOX": ["RFX1", "RFX2", "RFX3"], "NFY": ["NFYA", "NFYB", "NFYC"],
         "SP1": ["SP1"], "BRACHYURY": ["TBXT"], "BAPX1": ["NKX3-2"], "E2F": ["E2F1", "E2F2", "E2F3"], "FRA1": ["FOSL1"],
         "FRA2": ["FOSL2"], "HNF1": ["HNF1A", "HNF1B"], "HNF6B": ["ONECUT2"], "NFAT": ["NFATC1", "NFATC2", "NFATC3", "NFATC4"],
         "OCT4SOX2TCFNANOG": ["POU5F1", "SOX2", "NANOG"], "RUNX": ["RUNX1", "RUNX2", "RUNX3"], "RUNXAML": ["RUNX1"],
         "SLUG": ["SNAI2"], "SNAIL1": ["SNAI1"], "TWIST": ["TWIST1", "TWIST2"], "TBET": ["TBX21"], "SCL": ["TAL1"],
         "TCFCP2L1": ["TFCP2L1"], "ZFP281": ["ZNF281"], "RORG": ["RORC"], "JUNAP1": ["JUN"], "CJUNCRE": ["JUN"],
         "STAT5": ["STAT5A", "STAT5B"], "LRF": ["ZBTB7A"], "PR": ["PGR"], "ZAC1": ["PLAGL1"], "RBPJ1": ["RBPJ"],
         "ISRE": ["IRF1", "IRF2", "IRF9"], "GATA": ["GATA1", "GATA2", "GATA3", "GATA4", "GATA5", "GATA6"],
         "CRE": ["CREB1", "ATF1"], "NF1HALFSITE": ["NFIA", "NFIB", "NFIC", "NFIX"], "TATABOX": ["TBP"], "SP5": ["SP5"], "KLF": [], "MEF2A": ["MEF2A"], "MEF2B": ["MEF2B"], "MEF2C": ["MEF2C"], "MEF2D": ["MEF2D"]}


def symbols(name, valid):
    pre = name.split("(")[0].split("/")[0]
    out = []
    for part in pre.split(":"):
        p = part.upper().replace("-", "")
        if p in ALIAS:
            out += ALIAS[p]
        elif part.upper() in valid:
            out.append(part.upper())
        elif p in valid:
            out.append(p)
    if pre.upper() in ALIAS:
        out = ALIAS[pre.upper()]
    return sorted(set(s for s in out if s in valid))


peaks = pd.read_csv(f"{BASE}/PeakAnnot/sc/features_ATAC.tsv.gz", sep="\t", usecols=["id"]).id
pidx = pd.Series(np.arange(len(peaks)), index=peaks)
rows, mot = [], []
for f in sorted(glob.glob(f"{D}/find/c*.txt")):
    d = pd.read_csv(f, sep="\t", usecols=[0, 3])
    d.columns = ["peak", "motif"]
    d = d.drop_duplicates()
    rows.append(pidx.reindex(d.peak).to_numpy())
    mot.append(d.motif.to_numpy())
    print(f, len(d), flush=True)
r = np.concatenate(rows)
m = np.concatenate(mot)
ok = ~np.isnan(r)
names, mc = np.unique(m[ok], return_inverse=True)
PM = sp.csr_matrix((np.ones(ok.sum(), bool), (r[ok].astype(int), mc)), shape=(len(peaks), len(names)))
sp.save_npz(f"{D}/peak_motif.npz", PM)
rna = pd.read_csv(f"{BASE}/PeakAnnot/sc/features_RNA.tsv.gz", sep="\t", usecols=["id"]).id
valid = set(rna.str.upper()) & set(pd.read_csv(os.environ["TF_LIST"], sep="\t").Symbol)
mt = pd.DataFrame({"motif": names, "symbols": [",".join(symbols(n, valid)) for n in names]})
mt.to_csv(f"{D}/motifs.tsv", sep="\t", index=False)
tfm = mt.assign(tf=mt.symbols.str.split(",")).explode("tf").query("tf != ''").dropna()
tfs = sorted(tfm.tf.unique())
cols = []
for t in tfs:
    j = np.where(np.isin(names, tfm.motif[tfm.tf == t]))[0]
    cols.append(sp.csr_matrix(PM[:, j].sum(1) > 0))
PT = sp.hstack(cols).tocsr()
sp.save_npz(f"{D}/peak_tf_motif.npz", PT)
tfm.groupby("tf").motif.apply(";".join).reset_index().assign(
    frac_peaks=np.asarray(PT.mean(0)).ravel()).to_csv(f"{D}/tf_motif.tsv", sep="\t", index=False)
print("motifs", len(names), "with TF", (mt.symbols != "").sum(), "| TFs with motif", len(tfs))
print(mt[mt.symbols == ""].motif.str.split("/").str[0].tolist())
