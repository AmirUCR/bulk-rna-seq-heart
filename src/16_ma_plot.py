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
DE_DIR          = Path(os.environ["WILCOXON_OUT"])

# ---- config --------------------------------------------------------------
COND2       = os.environ.get("COND2", "DCM")     # disease group in this comparison
COMPARISON  = f"NF_vs_{COND2}"
LFC_THR     = 0.585
FDR_THR     = 0.05
TOP_LABELS  = 0
POINT_SIZE  = 8

# Files
samples_file     = SHARED_DATA_DIR / "samples.tsv"
de_analysis_file = DE_DIR / f"{COMPARISON}.wilcoxon.csv"
output_csv_file  = DE_DIR / f"{COMPARISON}.L2FC.{LFC_THR}.filtered.csv"
output_pdf_file  = DE_DIR / f"{COMPARISON}.L2FC.{LFC_THR}.ma_plot.pdf"
output_png_file  = DE_DIR / f"{COMPARISON}.L2FC.{LFC_THR}.ma_plot.png"

# ---- read ----------------------------------------------------------------
df          = pd.read_csv(de_analysis_file)
samples_tsv = pd.read_csv(samples_file, sep="\t")

id_col = "sample" if "sample" in samples_tsv.columns else "sample_id"
manifest_ids = set(samples_tsv[id_col].astype(str))
conditions   = sorted(samples_tsv["condition"].astype(str).unique(),
                      key=len, reverse=True)          # longest-first (ICM_DCM before DCM)

# DE CSV names per-sample columns as <sample_id>_<condition>. Strip the
# trailing _<condition> so they can be matched against manifest sample IDs.
cond_pat = re.compile(r"_(" + "|".join(map(re.escape, conditions)) + r")$")

def strip_cond(col):
    return cond_pat.sub("", col)

# map each DE column to its cleaned id; keep only those that are real samples
sample_cols = [c for c in df.columns if strip_cond(c) in manifest_ids]
if len(sample_cols) < 2:
    raise ValueError("Fewer than 2 sample columns from samples.tsv found in the results table.")

# expression filter (optional): > 5 (CPM) in at least 2 samples
# df = df[(df[sample_cols] > 5).sum(axis=1) >= 2].copy()
df.to_csv(output_csv_file, index=False)

# ---- columns -------------------------------------------------------------
gene_col = "name" if "name" in df.columns else df.columns[0]
if "log2FoldChange" in df.columns:
    lfc_col = "log2FoldChange"
elif "logFC" in df.columns:
    lfc_col = "logFC"
elif "foldChange" in df.columns:
    df["__log2FC__"] = np.log2(pd.to_numeric(df["foldChange"], errors="coerce"))
    lfc_col = "__log2FC__"
else:
    raise ValueError("No log2FoldChange/logFC/foldChange column found.")

padj_col = next((c for c in ["FDR", "padj", "PAdj"] if c in df.columns), None)

# ---- vectors -------------------------------------------------------------
genes     = df[gene_col].astype(str).to_numpy()
lfc       = pd.to_numeric(df[lfc_col], errors="coerce").to_numpy()
mean_norm = df[sample_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1).to_numpy()
padj      = (pd.to_numeric(df[padj_col], errors="coerce").to_numpy()
             if padj_col else np.full(lfc.shape, np.nan))

valid = np.isfinite(mean_norm) & np.isfinite(lfc)
if np.isfinite(padj).any():
    valid &= np.isfinite(padj)
mean_norm, lfc, padj, genes = mean_norm[valid], lfc[valid], padj[valid], genes[valid]

# ---- masks ---------------------------------------------------------------
sig = np.isfinite(padj) & (padj < FDR_THR)
up  = sig & (lfc >=  LFC_THR)
dn  = sig & (lfc <= -LFC_THR)
ns  = ~(up | dn)

# ---- plot ----------------------------------------------------------------
def remove_spines(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

cond2_label = COND2.replace("_", " + ")   # ICM_DCM -> "ICM + DCM"

fig, ax = plt.subplots(figsize=(7.5, 6))
ax.scatter(mean_norm[ns], lfc[ns], s=POINT_SIZE, alpha=0.6, linewidths=0, color="#131313")
ax.scatter(mean_norm[up], lfc[up], s=POINT_SIZE + 4, marker="^", alpha=0.9, linewidths=0,
           label=f"Up in {cond2_label} ({int(up.sum())})", color="#d7191c")
ax.scatter(mean_norm[dn], lfc[dn], s=POINT_SIZE + 4, marker="v", alpha=0.9, linewidths=0,
           label=f"Down in {cond2_label} ({int(dn.sum())})", color="#2c7bb6")

ax.set_xscale("log")
ax.axhline(0,        ls="-",  lw=1, color="#DE316B4D")
ax.axhline(LFC_THR,  ls="--", lw=1, color="#2E1760")
ax.axhline(-LFC_THR, ls="--", lw=1, color="#2E1760")
if np.isfinite(padj).any():
    ax.text(0.02, 0.98, f"FDR < {FDR_THR}", transform=ax.transAxes,
            va="top", ha="left", fontsize=9)

ax.set_xlabel("Mean of Normalized Counts (CPM)")
ax.set_ylabel(r"${\rm Log}_{2}$ Fold Change (" + f"{cond2_label} / NF)")
ax.legend(frameon=False, fontsize=9, loc="upper right")

def label_top(mask, k):
    idx = np.where(mask)[0]
    if idx.size == 0 or k <= 0:
        return
    order = np.lexsort((-np.abs(lfc[idx]), padj[idx]))   # low FDR, then large |LFC|
    for i in idx[order[:min(k, idx.size)]]:
        ax.annotate(genes[i], (mean_norm[i], lfc[i]),
                    xytext=(3, 1), textcoords="offset points", fontsize=8)

k_each = TOP_LABELS // 2
label_top(up, k_each + (TOP_LABELS % 2))
label_top(dn, k_each)

remove_spines(ax)
plt.tight_layout()
plt.savefig(output_png_file, dpi=300)
plt.savefig(output_pdf_file)
print(f"Up: {int(up.sum())}  Down: {int(dn.sum())}  NS: {int(ns.sum())}")
print(f"Wrote {output_png_file.name} / {output_pdf_file.name} and {output_csv_file.name}")