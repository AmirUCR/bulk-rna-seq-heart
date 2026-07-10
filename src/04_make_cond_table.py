#!/usr/bin/env python3
import os
import pandas as pd
from pathlib import Path

exclude_rows = ["SRR7426823"]  # This one fails the download

# Define the required environment variables
REQUIRED_VARS = ["SHARED_DATA_DIR"]

# Check if each variable is set
for var in REQUIRED_VARS:
    if not os.environ.get(var):
        raise EnvironmentError(
            f"'{var}' not set; source 00_vars.sh before running this script."
        )

SHARED_DATA_DIR = Path(os.environ.get("SHARED_DATA_DIR"))
df = pd.read_csv(f"{SHARED_DATA_DIR}/SraRunTable.csv")

cols_to_keep = ["Run", "disease"]
df = df[cols_to_keep]
df = df[~df['Run'].isin(exclude_rows)]
df = df.rename(columns={"Run": "sample_id", "disease": "condition"})

# df['r1'] = df['sample_id'] + "_1.fastq.gz"
df['r1'] = df['sample_id'] + ".fastq.gz"
# df['r2'] = df['sample_id'] + "_2.fastq.gz"


def rename(name):
    if name == "non-failing":
        return "NF"
    elif name == "ischemic cardiomyopathy":
        return "ICM"
    else:
        return "DCM"


df['condition'] = df['condition'].apply(lambda x: rename(x))

df.to_csv(f'{SHARED_DATA_DIR}/samples.tsv', sep='\t', index=False)


def rename_icm_dcm(name):
    if name == "NF":
        return "NF"
    else:
        return "ICM_DCM"


df['condition'] = df['condition'].apply(lambda x: rename_icm_dcm(x))
df.to_csv(f'{SHARED_DATA_DIR}/disease_samples.tsv', sep='\t', index=False)
