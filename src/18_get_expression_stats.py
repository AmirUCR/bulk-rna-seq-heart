#!/usr/bin/env python3
import os
import pandas as pd
from pathlib import Path

# ---- environment ---------------------------------------------------------
REQUIRED_VARS = ["WILCOXON_OUT"]
for var in REQUIRED_VARS:
    if not os.environ.get(var):
        raise EnvironmentError(
            f"'{var}' not set; source 00_vars.sh before running this script."
        )

DE_DIR = Path(os.environ["WILCOXON_OUT"])

conditions = ["ICM", "DCM", "ICM_DCM"]
lfc_thrs   = [0.585, 1]
fdr_thrs   = [0.05]

rows = []
for cond2 in conditions:
    comparison = f"NF_vs_{cond2}"
    csv_path = DE_DIR / f"{comparison}.wilcoxon.csv"
    if not csv_path.exists():
        print(f"WARNING: missing {csv_path}, skipping")
        continue
    df = pd.read_csv(csv_path)

    lfc = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    fdr = pd.to_numeric(df["FDR"], errors="coerce")

    for lfc_thr in lfc_thrs:
        for fdr_thr in fdr_thrs:
            sig  = fdr < fdr_thr                       # significant = FDR below threshold
            up   = int((sig & (lfc >=  lfc_thr)).sum())
            down = int((sig & (lfc <= -lfc_thr)).sum())
            rows.append({
                "comparison":            comparison,
                "fdr_thresh":            fdr_thr,
                "log2FoldChange_thresh": lfc_thr,
                "num_up_in_disease":     up,
                "num_down_in_disease":   down,
                "num_total_sig":         up + down,
            })

out = pd.DataFrame(rows)
out_path = DE_DIR / "de_summary_stats.csv"
out.to_csv(out_path, index=False)
print(out.to_string(index=False))
print(f"\nWrote {out_path}")
