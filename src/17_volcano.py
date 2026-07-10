#!/usr/bin/env python3
import os
import re
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 600

# ---- environment ---------------------------------------------------------
REQUIRED_VARS = ["SHARED_DATA_DIR", "WILCOXON_OUT"]
for var in REQUIRED_VARS:
    if not os.environ.get(var):
        raise EnvironmentError(
            f"'{var}' not set; source 00_vars.sh before running this script."
        )

SHARED_DATA_DIR = Path(os.environ["SHARED_DATA_DIR"])
DE_DIR = Path(os.environ["WILCOXON_OUT"])

# ---- config --------------------------------------------------------------
COND2 = os.environ.get("COND2", "DCM")     # disease group in this comparison
COMPARISON = f"NF_vs_{COND2}"
LFC_THR = 0.585
FDR_THR = 0.05
TOP_LABELS = 0
POINT_SIZE = 8

samples_file = SHARED_DATA_DIR / "samples.tsv"
de_analysis_file = DE_DIR / f"{COMPARISON}.wilcoxon.csv"
output_pdf_file = DE_DIR / f"{COMPARISON}.L2FC.{LFC_THR}.volcano_plot.pdf"
output_png_file = DE_DIR / f"{COMPARISON}.L2FC.{LFC_THR}.volcano_plot.png"

# ---- read ----------------------------------------------------------------
df = pd.read_csv(de_analysis_file)
samples_tsv = pd.read_csv(samples_file, sep="\t")

id_col = "sample" if "sample" in samples_tsv.columns else "sample_id"
manifest_ids = set(samples_tsv[id_col].astype(str))
conditions = sorted(samples_tsv["condition"].astype(str).unique(),
                    # longest-first (ICM_DCM before DCM)
                    key=len, reverse=True)

# DE CSV names per-sample columns as <sample_id>_<condition>; strip the
# trailing _<condition> so they match manifest sample IDs.
cond_pat = re.compile(r"_(" + "|".join(map(re.escape, conditions)) + r")$")
sample_cols = [c for c in df.columns if cond_pat.sub("", c) in manifest_ids]

# light expression filter (> 5 CPM in >= 2 samples), mirrors MA script
# if len(sample_cols) >= 2:
#     df = df[(df[sample_cols] > 5).sum(axis=1) >= 2].copy()

# ---- columns -------------------------------------------------------------
gene_col = "name" if "name" in df.columns else df.columns[0]
if "log2FoldChange" in df.columns:
    lfc_col = "log2FoldChange"
elif "logFC" in df.columns:
    lfc_col = "logFC"
elif "foldChange" in df.columns:
    df["__log2FC__"] = np.log2(pd.to_numeric(
        df["foldChange"], errors="coerce"))
    lfc_col = "__log2FC__"
else:
    raise ValueError("No log2FoldChange/logFC/foldChange column found.")

padj_col = next((c for c in ["FDR", "padj", "PAdj"] if c in df.columns), None)
p_col = next(
    (c for c in ["PValue", "pvalue", "pval", "p"] if c in df.columns), None)
if p_col is None and padj_col is None:
    raise ValueError("Need a p-value or FDR column for a volcano.")

# ---- vectors -------------------------------------------------------------
genes = df[gene_col].astype(str).to_numpy()
lfc = pd.to_numeric(df[lfc_col], errors="coerce").to_numpy()
padj = (pd.to_numeric(df[padj_col], errors="coerce").to_numpy()
        if padj_col else np.full(lfc.shape, np.nan))
pval = (pd.to_numeric(df[p_col], errors="coerce").to_numpy()
        if p_col else padj)   # fall back to FDR for the y-axis if no raw p

# y-axis: -log10(p), floored so exact zeros don't become inf
pval = np.clip(pval, 1e-300, 1.0)
neglog10p = -np.log10(pval)

valid = np.isfinite(lfc) & np.isfinite(neglog10p)
lfc, neglog10p, padj, genes = lfc[valid], neglog10p[valid], padj[valid], genes[valid]
# keep aligned copy for cutoff calc
pval_v = pval[valid]

# ---- masks (significance by FDR, direction by LFC) -----------------------
sig = np.isfinite(padj) & (padj < FDR_THR)
up = sig & (lfc >= LFC_THR)
dn = sig & (lfc <= -LFC_THR)
ns = ~(up | dn)

# horizontal cutoff = -log10(p) of the least-significant gene still clearing
# FDR < 0.05 (where the FDR threshold lands on the raw-p axis)
if sig.any():
    p_at_fdr = pval_v[sig].max()
    y_cut = -np.log10(min(p_at_fdr, 1.0))
else:
    y_cut = None

# ---- plot ----------------------------------------------------------------


def remove_spines(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


cond2_label = COND2.replace("_", " + ")   # ICM_DCM -> "ICM + DCM"

fig, ax = plt.subplots(figsize=(6.5, 6.5))
ax.scatter(lfc[ns], neglog10p[ns], s=POINT_SIZE,
           alpha=0.5, linewidths=0, color="#9a9a9a")
ax.scatter(lfc[up], neglog10p[up], s=POINT_SIZE + 4, alpha=0.9, linewidths=0,
           label=f"Up in {cond2_label} ({int(up.sum())})", color="#d7191c")
ax.scatter(lfc[dn], neglog10p[dn], s=POINT_SIZE + 4, alpha=0.9, linewidths=0,
           label=f"Down in {cond2_label} ({int(dn.sum())})", color="#2c7bb6")

# guides
ax.axvline(LFC_THR,  ls="--", lw=1, color="#2E1760")
ax.axvline(-LFC_THR, ls="--", lw=1, color="#2E1760")
if y_cut is not None:
    ax.axhline(y_cut, ls="--", lw=1, color="#2E1760")
    ax.text(ax.get_xlim()[1], y_cut, f"  FDR {FDR_THR}", va="bottom", ha="right",
            fontsize=8, color="#2E1760")

ax.set_xlabel(r"${\rm Log}_{2}$ Fold Change (" + f"{cond2_label} / NF)")
ax.set_ylabel(r"$-{\rm Log}_{10}$ P-value")
ax.legend(frameon=True, fontsize=9, loc="upper right")

# label top genes (low FDR, large |LFC|), split up/down


def label_top(mask, k):
    idx = np.where(mask)[0]
    if idx.size == 0 or k <= 0:
        return
    order = np.lexsort((-np.abs(lfc[idx]), padj[idx]))
    for i in idx[order[:min(k, idx.size)]]:
        ax.annotate(genes[i], (lfc[i], neglog10p[i]),
                    xytext=(3, 1), textcoords="offset points", fontsize=7)


k_each = TOP_LABELS // 2
label_top(up, k_each + (TOP_LABELS % 2))
label_top(dn, k_each)

remove_spines(ax)
plt.tight_layout()
plt.savefig(output_png_file, dpi=300)
plt.savefig(output_pdf_file)
print(f"Up: {int(up.sum())}  Down: {int(dn.sum())}  NS: {int(ns.sum())}")
print(f"Wrote {output_png_file.name} / {output_pdf_file.name}")
