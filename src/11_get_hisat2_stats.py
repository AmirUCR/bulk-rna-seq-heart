#!/usr/bin/env python3
"""
Summarize HISAT2 alignment logs across samples.

Expects logs at:  <SHARED_BAM_DIR>/<SAMPLE_ID>/trimmed.hisat2.log
(single-end HISAT2 stderr format)

Writes a per-sample CSV and prints aggregate statistics.
"""
import os
import re
import csv
import sys
from pathlib import Path

# Define the required environment variables
REQUIRED_VARS = ["LOCAL_RESULTS_DIR", "SHARED_BAM_DIR"]

# Check if each variable is set
for var in REQUIRED_VARS:
    if not os.environ.get(var):
        raise EnvironmentError(
            f"'{var}' not set; source 00_vars.sh before running this script."
        )

SHARED_BAM_DIR = Path(os.environ.get("SHARED_BAM_DIR"))
LOCAL_RESULTS_DIR = Path(os.environ.get("LOCAL_RESULTS_DIR"))

# --- regexes matched against the HISAT2 stderr block ------------------------
RE_TOTAL = re.compile(r'^\s*([\d,]+)\s+reads;\s+of these:')
RE_ZERO = re.compile(r'([\d,]+)\s+\(([\d.]+)%\)\s+aligned 0 times')
RE_UNIQUE = re.compile(r'([\d,]+)\s+\(([\d.]+)%\)\s+aligned exactly 1 time')
RE_MULTI = re.compile(r'([\d,]+)\s+\(([\d.]+)%\)\s+aligned >1 times')
RE_OVERALL = re.compile(r'([\d.]+)%\s+overall alignment rate')


def _num(s):
    """'40,752,470' -> 40752470"""
    return int(s.replace(",", ""))


def parse_log(path):
    """Return a dict of metrics for one HISAT2 log, or None if unparseable."""
    text = path.read_text()

    def grab(rx, group=1, cast=_num):
        m = rx.search(text)
        return cast(m.group(group)) if m else None

    total = grab(RE_TOTAL)
    if total is None:
        return None  # not a recognizable HISAT2 log

    rec = {
        "total_reads":      total,
        "unaligned":        grab(RE_ZERO),
        "unaligned_pct":    grab(RE_ZERO, 2, float),
        "unique":           grab(RE_UNIQUE),
        "unique_pct":       grab(RE_UNIQUE, 2, float),
        "multi":            grab(RE_MULTI),
        "multi_pct":        grab(RE_MULTI, 2, float),
        "overall_pct":      grab(RE_OVERALL, 1, float),
    }
    return rec

# ---

out = LOCAL_RESULTS_DIR / "hisat2_summary.csv"
logname = "trimmed.hisat2.log"

root = Path(SHARED_BAM_DIR)
logs = sorted(root.glob(f"*/{logname}"))
if not logs:
    sys.exit(f"No logs matching */{logname} under {root}")

rows = []
for log in logs:
    sample = log.parent.name          # <sample_id>/trimmed.hisat2.log -> sample_id
    rec = parse_log(log)
    if rec is None:
        print(f"WARNING: could not parse {log}", file=sys.stderr)
        continue
    rec = {"sample": sample, **rec}
    rows.append(rec)

if not rows:
    sys.exit("No parseable logs found.")

# --- write per-sample CSV ----------------------------------------------
fields = ["sample", "total_reads", "unaligned", "unaligned_pct",
          "unique", "unique_pct", "multi", "multi_pct", "overall_pct"]
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)

# --- aggregate stats ----------------------------------------------------


def stats(key):
    """Return (mean, sample_sd) for a metric, or (None, None) if empty."""
    vals = [r[key] for r in rows if r[key] is not None]
    if not vals:
        return (None, None)
    mean = sum(vals) / len(vals)
    # sample standard deviation (n-1 denominator); 0 if only one value
    if len(vals) > 1:
        var = sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)
        sd = var ** 0.5
    else:
        sd = 0.0
    return (mean, sd)


n = len(rows)
tot = sum(r["total_reads"] for r in rows)
ov_mean, ov_sd = stats("overall_pct")
uq_mean, uq_sd = stats("unique_pct")
mu_mean, mu_sd = stats("multi_pct")

print(f"Samples parsed:        {n}")
print(f"Total reads (all):     {tot:,}")
print(f"Reads/sample (mean):   {tot // n:,}")
print(f"Overall align rate %:  mean {ov_mean:.2f}  SD {ov_sd:.2f}")
print(f"Uniquely aligned  %:   mean {uq_mean:.2f}  SD {uq_sd:.2f}")
print(f"Multi-aligned     %:   mean {mu_mean:.2f}  SD {mu_sd:.2f}")

# flag any low-alignment outliers worth inspecting
low = [r for r in rows if r["overall_pct"]
       is not None and r["overall_pct"] < 80]
if low:
    print("\nSamples below 80% overall alignment:")
    for r in sorted(low, key=lambda x: x["overall_pct"]):
        print(f"  {r['sample']}: {r['overall_pct']:.2f}%")

print(f"\nWrote {out}")
