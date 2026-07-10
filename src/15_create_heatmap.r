#!/usr/bin/env Rscript
#
# Heatmap from a Wilcoxon differential-expression table (edgeR-style columns).
# Rows = significant genes (FDR + |log2FC| cutoffs), ordered by fold change
# (up-in-disease at top) so layout is identical across LFC thresholds.
# Columns = samples, grouped by condition (not clustered). Values are row
# z-scored CPMs.
#
suppressPackageStartupMessages(library(gplots))

## ---- environment ---------------------------------------------------------
for (v in c("COUNTS_OUT", "SHARED_DATA_DIR", "WILCOXON_OUT")) {
  if (Sys.getenv(v) == "") stop(sprintf("%s not set; source 00_vars.sh before Rscript", v))
}
SHARED_DATA_DIR <- Sys.getenv("SHARED_DATA_DIR")
WILCOXON_OUT    <- Sys.getenv("WILCOXON_OUT")

## ---- config --------------------------------------------------------------
COMPARISON  <- Sys.getenv("COMPARISON", unset = "NF_vs_DCM")  # e.g. NF_vs_ICM, NF_vs_ICM_DCM
MIN_FDR     <- 0.05
MIN_ABS_LFC <- 0.585   # require |log2FC| >= this for display (0 to disable)
TOP_N       <- 0       # show top N up + top N down by FDR (0 = show all)

WIDTH   <- 12
HEIGHT  <- 14
MARGINS <- c(14, 12)   # heatmap.2 margins: c(bottom, right)

samples_file <- file.path(SHARED_DATA_DIR, "samples.tsv")
counts_file  <- file.path(WILCOXON_OUT, paste0(COMPARISON, ".wilcoxon.csv"))

lfc_tag     <- sprintf("L2FC.%s", format(MIN_ABS_LFC, trim = TRUE, drop0trailing = TRUE))
output_file <- file.path(WILCOXON_OUT, paste0(COMPARISON, ".", lfc_tag, ".heatmap.pdf"))
png_file    <- file.path(WILCOXON_OUT, paste0(COMPARISON, ".", lfc_tag, ".heatmap.png"))

cat("# Tool: Create Heatmap (Wilcoxon)\n")
cat("# Input:", counts_file, "\n")
cat("# Output:", output_file, "\n")

## ---- load + validate -----------------------------------------------------
data <- read.csv(counts_file, header = TRUE, as.is = TRUE, check.names = FALSE)
sample_data <- read.table(samples_file, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, check.names = FALSE)

# manifest uses 'sample'; fall back to 'sample_id'
id_col <- if ("sample" %in% colnames(sample_data)) "sample" else "sample_id"
stopifnot(id_col %in% colnames(sample_data),
          "condition" %in% colnames(sample_data),
          all(c("FDR", "name", "log2FoldChange") %in% colnames(data)))

## ---- select significant genes --------------------------------------------
data$FDR            <- as.numeric(data$FDR)
data$log2FoldChange <- as.numeric(data$log2FoldChange)
data <- subset(data, !is.na(FDR) & FDR <= MIN_FDR)
if (MIN_ABS_LFC > 0) {
  data <- subset(data, abs(log2FoldChange) >= MIN_ABS_LFC)
}
if (nrow(data) == 0) stop("No rows pass the cutoffs.")

if (TOP_N > 0) {
  up   <- head(data[data$log2FoldChange > 0, , drop = FALSE][
                 order(data$FDR[data$log2FoldChange > 0]), ], TOP_N)
  down <- head(data[data$log2FoldChange < 0, , drop = FALSE][
                 order(data$FDR[data$log2FoldChange < 0]), ], TOP_N)
  data <- rbind(up, down)
  cat("# Top-N selection:", nrow(up), "up +", nrow(down), "down =", nrow(data), "rows\n")
} else {
  cat("# Genes shown:", nrow(data), "\n")
}

# Fixed row order: highest positive LFC (up-in-disease) at top. This makes the
# quadrant pattern identical regardless of MIN_ABS_LFC or gene set.
data <- data[order(data$log2FoldChange, decreasing = TRUE), ]

## ---- resolve sample columns ----------------------------------------------
# DE CSV names per-sample columns as <sample_id>_<condition>; strip the
# trailing _<condition> so they match manifest sample IDs. Sort conditions
# longest-first so e.g. ICM_DCM is matched before DCM.
conds     <- unique(sample_data$condition)
conds     <- conds[order(nchar(conds), decreasing = TRUE)]
conds_pat <- paste0("_(", paste(conds, collapse = "|"), ")$")
data_ids  <- sub(conds_pat, "", colnames(data))

id_to_col   <- setNames(colnames(data), data_ids)          # cleaned id -> original column
sample_ids  <- intersect(data_ids, sample_data[[id_col]])  # CSV order
sample_cols <- id_to_col[sample_ids]                       # original (suffixed) names
stopifnot(!anyNA(sample_cols))
if (length(sample_cols) == 0) stop("No sample columns from samples.tsv found in the results table.")
cat("# Samples in heatmap:", length(sample_cols), "\n")

## ---- build value matrix + row z-scores -----------------------------------
values <- as.matrix(sapply(data[, sample_cols, drop = FALSE], as.numeric))
stopifnot(identical(unname(sample_cols), colnames(values)))
if (anyNA(values)) cat("# WARNING:", sum(is.na(values)), "NA cells in value matrix\n")
rownames(values) <- data$name
colnames(values) <- sample_ids                              # clean IDs on the axis

row_means <- rowMeans(values, na.rm = TRUE)
row_sds   <- apply(values, 1, sd, na.rm = TRUE)
row_sds[row_sds == 0 | is.na(row_sds)] <- 1                 # flat rows -> z = 0 (mid palette)
zscores   <- sweep(sweep(values, 1, row_means, "-"), 1, row_sds, "/")
stopifnot(identical(rownames(zscores), data$name))

## ---- adaptive sizing -----------------------------------------------------
cexRow      <- max(0.15, min(0.95, 22 / nrow(zscores)))
plot_height <- min(48, max(HEIGHT, nrow(zscores) * 0.13))
png_res     <- round(min(600, 8000 / plot_height))

# Keep the key (~1") and condition bar (~0.2") a fixed physical size as the
# figure grows: make the heatmap row proportional to total height.
LHEI <- c(1, plot_height)   # heatmap.2 inserts 0.2 for the ColSideColors bar

## ---- condition annotation bar --------------------------------------------
condition    <- sample_data$condition[match(sample_ids, sample_data[[id_col]])]
cond_palette <- setNames(c("#2c7bb6", "#d7191c", "#fdae61", "#abd9e9")[seq_along(unique(condition))],
                         unique(condition))
col_colors   <- cond_palette[condition]

## ---- draw ----------------------------------------------------------------
hr  <- hclust(dist(zscores))
dd  <- reorder(as.dendrogram(hr),
               wts = as.numeric(data$log2FoldChange), agglo.FUN = mean)

draw_heatmap <- function() {
  op <- par(oma = c(0, 0, 3, 0))   # reserve top strip on the device for the legend
  on.exit(par(op))

  heatmap.2(
    zscores,
    col           = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(75),
    breaks        = seq(-3, 3, length.out = 76),  # symmetric; beyond +/-3 saturates
    ColSideColors = col_colors,
    density.info  = "none",
    Rowv          = dd,          # fixed row order (by LFC), don't cluster genes
    Colv          = FALSE,          # keep samples grouped by condition
    dendrogram    = "row",
    trace         = "none",
    margins       = MARGINS,
    lhei          = LHEI,
    labRow        = rownames(zscores),
    cexRow        = cexRow,
    srtCol        = 45,
    cexCol        = 0.8,
    key.title     = NA,
    key.xlab      = "row z-score",
    keysize       = 1.0,
    lwid          = c(2, 8),
    key.par       = list(mar = c(3, 1, 2, 1)),
    key.xtickfun  = function() {
      # force an integer tick at every break value -3..3 (prevents -2 dropping)
      list(at = seq(0, 1, length.out = 7),
           labels = c("-3", "-2", "-1", "0", "1", "2", "3"))
    }
  )

  # condition legend in the reserved outer top margin (won't collide with rows)
  par(xpd = NA)
  legend("top", inset = c(0, -0.02), horiz = TRUE,
         legend = gsub("_", " + ", names(cond_palette)),
         fill = cond_palette, border = NA, bty = "n",
         cex = 0.9, title = "Condition")
}

pdf(output_file, width = WIDTH, height = plot_height); draw_heatmap(); invisible(dev.off())
png(png_file, width = WIDTH, height = plot_height, units = "in", res = png_res); draw_heatmap(); invisible(dev.off())
cat("# Wrote:", output_file, "and", png_file, "\n")
