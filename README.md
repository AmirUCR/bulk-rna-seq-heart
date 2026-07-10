# Bulk RNA-Seq for Public Heart Data
Do not run `00_vars.sh` and `01_common.sh` manually.

# Script order and execution
Install the conda environment: 
> `conda env create -f env.yaml`

Then give execution permission:
1. `cd src`
1. `chmod +x *`

Inspect `data/accessions.txt` which contains the files that will be downloaded.
The file `data/SraRunTable.csv` contains more information about the accessions. Inspect the samples files to see for which conditions differential expression analyses will be performed.

Inspect `src/00_vars.sh` and set the number of THREADS you want your machine to use. You may also update the download URLs for the reference genome and annotation. The other paths should remain as is.

1. To begin reference genome downloads, run
    > `./02_download_ref.sh`

1. To begin FASTQ downloads, run
    > `./03_download_sra.sh`

1. You may skip `./04_make_cond_table.py` as the `data/samples.tsv` and `data/disease_samples.tsv` are included in the repo.

1. When the `.fastq.gz` files are downloaded under `data/reads/untrimmed`, run
    > `./05_fastqc.sh`

    When finished, inspect `src/results/reads/untrimmed/multiqc_report.html` and check if you want to do any trimming.

1. Open `src/06_trim.sh` and inspect the `fastp` parameters. If you need to trim your reads, add the options here. For the heart dataset, these are single-end layout so the `run_fastp_single()` will run. Then run
    > `./06_trim.sh`

1. Run
    > `./07_fastqc_again.sh`
    
    to check post-trim results under `src/reads/trimmed/multiqc_report.html`

1. Run
    > `./08_hisat2_index.sh`

    to generate a reference genome index.

1. Continue with
    > `./09_find_strandedness.sh`
    
    > `./10_hisat2.sh`

    > `./11_get_hisat2_stats.sh`

    Check `src/results/hisat2_summary.csv` to get a count of library depth and alignment rates. Individual library stats are under `results/bam/SRRXXX` flagstat and hisat2.log files.

1. Run
    > `./12_featurecounts.sh`

    > `source 01_common.sh && ./13_parse_featurecounts.r`

    to generate a counts matrix under `src/results/counts`

1. Run
    > `source 01_common.sh && ./14_wilcoxon.r` 
    
    to run this analysis for NF vs DCM and generate results under `src/results/wilcoxon_de`
    
    Then, open `src/14_wilcoxon.r` manually. You will need to edit the source. Specifically, this block:

    ```
    ref_grp  <- "NF"     # group A / baseline
    trt_grp  <- "DCM"    # group B
    fdrThres <- 0.05

    COUNTS_OUT <- Sys.getenv("COUNTS_OUT")
    counts_file <- file.path(COUNTS_OUT, "counts_se.txt")

    SHARED_DATA_DIR <- Sys.getenv("SHARED_DATA_DIR")
    samples_file <- file.path(SHARED_DATA_DIR, "samples.tsv")
    ```

    Change `trt_grp  <- "DCM"` to `trt_grp  <- "ICM"` and rerun the script for NF vs ICM.

    Then change 
    > `trt_grp  <- "ICM"` 
    
    to 
    
    > `trt_grp  <- "ICM_DCM"`
    
    and also 
    
    > `samples_file <- file.path(SHARED_DATA_DIR, "samples.tsv")` 
    
    to
    
    > `samples_file <- file.path(SHARED_DATA_DIR, "disease_samples.tsv")`

    and rerun the script for NF vs (ICM + DCM)

1. Run
    > `source 01_common.sh && ./15_create_heatmap.r`

    to create a heatmap for NF vs DCM under `src/results/wilcoxon_de` with MIN_ABS_LFC < 0.585. Make the necessary adjustments for other combinations, and __REMEMBER to change `samples.tsv` to `disease_samples.tsv` for ICM+DCM.__

1. Run
    > `source 01_common.sh && ./16_ma_plot.py`

    to create an MA plot for NF vs DCM with MIN_ABS_LFC < 0.585. Make the necessary adjustments for other combinations, and __REMEMBER to change `samples.tsv` to `disease_samples.tsv` for ICM+DCM.__

1. Run
    > `source 01_common.sh && ./17_volcano.py`

    to create an volcano plot for NF vs DCM with MIN_ABS_LFC < 0.585. Make the necessary adjustments for other combinations, and __REMEMBER to change `samples.tsv` to `disease_samples.tsv` for ICM+DCM.__

# Support
Create a GitHub issue if you run into trouble and we will help you out.
