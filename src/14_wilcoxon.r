#!/usr/bin/env Rscript
# Wilcoxon rank-sum DE: NF vs ICM/DCM, output formatted to match the old edgeR table.
# Test/normalization logic follows Li et al. 2022 (Genome Biol 23:79).

suppressWarnings(suppressMessages(library(edgeR)))

for (v in c("COUNTS_OUT","SHARED_DATA_DIR","WILCOXON_OUT")) {
  if (Sys.getenv(v) == "") stop(sprintf("%s not set; source 00_vars.sh before Rscript", v))
}

## ---- inputs (edit paths) -------------------------------------------------
ref_grp  <- "NF"     # group A / baseline
trt_grp  <- "DCM"    # group B
fdrThres <- 0.05

COUNTS_OUT <- Sys.getenv("COUNTS_OUT")
counts_file <- file.path(COUNTS_OUT, "counts_se.txt")

SHARED_DATA_DIR <- Sys.getenv("SHARED_DATA_DIR")
samples_file <- file.path(SHARED_DATA_DIR, "samples.tsv")

full_count_file <- file.path(COUNTS_OUT, "counts_se.csv")

WILCOXON_OUT <- Sys.getenv("WILCOXON_OUT")
out_file_format_name <- sprintf("%s_vs_%s.wilcoxon.csv", ref_grp, trt_grp)
output_file     <- file.path(WILCOXON_OUT, out_file_format_name)

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

## ---- read ----------------------------------------------------------------
readCount <- read.table(full_count_file, header = TRUE, row.names = 1, sep = ",",
                        stringsAsFactors = FALSE, check.names = FALSE)
meta <- read.table(samples_file, header = TRUE, stringsAsFactors = FALSE)
if (!"sample" %in% colnames(meta) && "sample_id" %in% colnames(meta)) {
  meta$sample <- meta$sample_id
}
stopifnot(all(c("sample","condition") %in% colnames(meta)))

cat("Full count table n rows:", nrow(readCount), "\n")

## ---- subset columns to the two groups; align conditions BY NAME ----------
keep_samples <- meta$sample[meta$condition %in% c(ref_grp, trt_grp)]
readCount <- readCount[, colnames(readCount) %in% keep_samples, drop = FALSE]
conditions <- meta$condition[match(colnames(readCount), meta$sample)]
stopifnot(!any(is.na(conditions)))
conditions <- relevel(factor(conditions), ref = ref_grp)   # NF = level 1 (group A)
cat("Samples:", ncol(readCount), "|",
    paste(levels(conditions), table(conditions)[levels(conditions)],
          collapse = "  "), "\n")

## ---- normalize on FULL transcriptome -------------------------------------
y    <- DGEList(counts = readCount, group = conditions)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
count_norm <- as.data.frame(cpm(y))
aveLog     <- aveLogCPM(y); names(aveLog) <- rownames(y)   # edgeR-style logCPM

stopifnot(identical(colnames(count_norm), colnames(readCount)))

## ---- Wilcoxon rank-sum per gene ----------------------------------------
pvalues <- sapply(seq_len(nrow(count_norm)), function(i) {
  d <- cbind.data.frame(gene = as.numeric(t(count_norm[i, ])), conditions)
  wilcox.test(gene ~ conditions, data = d, exact = FALSE)$p.value
})

## ---- assemble edgeR-style columns ----------------------------------------
lvl <- levels(conditions)                                  # c("NF","DCM")
A <- colnames(count_norm)[conditions == lvl[1]]            # NF sample ids
B <- colnames(count_norm)[conditions == lvl[2]]            # DCM sample ids

baseMeanA <- rowMeans(count_norm[, A, drop = FALSE])
baseMeanB <- rowMeans(count_norm[, B, drop = FALSE])
baseMean  <- rowMeans(count_norm[, c(A, B), drop = FALSE])
log2FoldChange <- log2((baseMeanB + 1) / (baseMeanA + 1)) # B over A (DCM over NF)
foldChange <- 2 ^ log2FoldChange
logCPM <- aveLog[rownames(count_norm)]
FDR      <- p.adjust(pvalues, method = "BH")              # adjusted p (significance)
PAdj_hoc <- p.adjust(pvalues, method = "hochberg")        # FWER, conservative

out <- data.frame(
  name           = rownames(count_norm),
  baseMean       = baseMean,
  baseMeanA      = baseMeanA,
  baseMeanB      = baseMeanB,
  foldChange     = foldChange,
  log2FoldChange = log2FoldChange,
  logCPM         = logCPM,
  PValue         = pvalues,
  PAdj_hoc       = PAdj_hoc,
  FDR            = FDR,
  check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE
)

## sort by p-value (then larger fold change), add cumulative expected FPs
out <- out[order(out$PValue, -out$foldChange), ]
out$falsePos <- seq_len(nrow(out)) * out$FDR

## append per-sample normalized columns (group A then group B), aligned by name
samp <- count_norm[out$name, c(A, B), drop = FALSE]
out  <- cbind(out, samp)

## final column order (matches the edgeR script's layout)
new_cols <- c("name", "baseMean", "baseMeanA", "baseMeanB",
              "foldChange", "log2FoldChange", "logCPM",
              "PValue", "PAdj_hoc", "FDR", "falsePos", A, B)
out <- out[, new_cols, drop = FALSE]

## ---- rename per-sample columns to <sample>_<condition> for output --------
# build sample_id -> condition map from meta
id2cond <- setNames(meta$condition, meta$sample)

# the per-sample columns in `out` are exactly c(A, B); rename those
sample_cols <- c(A, B)

if (any(is.na(id2cond[sample_cols]))) {
  stop("Some sample columns have no condition in meta: ",
       paste(sample_cols[is.na(id2cond[sample_cols])], collapse = ", "))
}

new_names    <- paste0(sample_cols, "_", id2cond[sample_cols])

# rename only those columns, leave the stat columns untouched
idx <- match(sample_cols, colnames(out))
colnames(out)[idx] <- new_names

write.csv(out, file = output_file, row.names = FALSE, quote = FALSE)
cat("# Tool: Wilcoxon rank-sum (Li et al. 2022)\n")
cat("Tested:", nrow(out),
    "| Significant (FDR <", fdrThres, "):", sum(out$FDR < fdrThres, na.rm = TRUE), "\n")
cat("Wrote:", output_file, "\n")
