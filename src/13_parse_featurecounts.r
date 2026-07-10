#!/usr/bin/env Rscript

# Parse featureCounts output into a simple counts matrix.

# Project structure
for (v in c("COUNTS_OUT","SHARED_DATA_DIR")) {
  if (Sys.getenv(v) == "") stop(sprintf("%s not set; source 00_vars.sh before Rscript", v))
}

# Input files
COUNTS_OUT <- Sys.getenv("COUNTS_OUT")
counts_file <- file.path(COUNTS_OUT, "counts_se.txt")

SHARED_DATA_DIR <- Sys.getenv("SHARED_DATA_DIR")
samples_file <- file.path(SHARED_DATA_DIR, "samples.tsv")

# Output file
output_file <- file.path(COUNTS_OUT, "counts_se.csv")

# Inform the user
cat("# Tool: Parse featureCounts\n")
cat("# Samples:", samples_file, "\n")
cat("# Input:", counts_file, "\n")

# Read sample sheet
sample_data <- read.table(
    samples_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

sample_data <- sample_data[, c("sample_id", "condition")]
names(sample_data) <- c("sample", "condition")

if (!all(c("sample", "condition") %in% names(sample_data))) {
    stop("samples.tsv must contain columns named 'sample_id' and 'condition'")
}

sample_data$condition <- factor(sample_data$condition)
sample_data$condition <- stats::relevel(
    sample_data$condition,
    ref = as.character(sample_data$condition[1])
)

# Read featureCounts output
df <- read.table(
    counts_file,
    header = TRUE,
    sep = "\t",
    comment.char = "#",
    check.names = FALSE
)

# featureCounts columns:
# 1 Geneid
# 2 Chr
# 3 Start
# 4 End
# 5 Strand
# 6 Length
# 7 gene_name
# 8+ sample count columns

# counts <- df[, c(1, 8:ncol(df))]

# Combine Geneid (col 1) and gene_name (col 7) into one identifier column
counts <- data.frame(
    name = paste(df[[1]], df[[7]], sep = "|"),
    df[, 8:ncol(df), drop = FALSE],
    check.names = FALSE,
    stringsAsFactors = FALSE
)

# Check sample count consistency
if ((ncol(counts) - 1) != nrow(sample_data)) {
    stop("Number of samples in counts.txt does not match number of rows in samples.tsv")
}

# Rename columns using sample sheet sample names
colnames(counts) <- c("name", sample_data$sample)

# Write output
write.csv(counts, file = output_file, row.names = FALSE, quote = FALSE)

# Inform the user
cat("# Output:", output_file, "\n")
