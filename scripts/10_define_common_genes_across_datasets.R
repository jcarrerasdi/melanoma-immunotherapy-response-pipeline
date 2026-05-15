# ============================================================
# 10_define_common_genes_across_datasets.R
# Cleaned/renamed version of: 09b_define_common_genes_across_datasets.r
# Purpose: reproducible TFM pipeline while preserving the output filenames used in the memoria.
# ============================================================

################################
# 09b_define_common_genes_across_datasets.R
################################

suppressPackageStartupMessages({
  library(dplyr)
})

# Project root: run scripts either from the project root or from the scripts/ folder.
base_dir <- normalizePath(getwd(), mustWork = TRUE)
if (basename(base_dir) == "scripts") base_dir <- dirname(base_dir)
data_processed_dir <- file.path(base_dir, "data_processed")

normalize_gene_ids <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.\\d+$", "", x)
  x <- gsub("\\s+", "", x)
  toupper(x)
}

train <- readRDS(file.path(data_processed_dir, "GSE160638_raw_counts_aligned.rds"))
gse91061 <- readRDS(file.path(data_processed_dir, "GSE91061_expression_aligned.rds"))
gse78220 <- readRDS(file.path(data_processed_dir, "GSE78220_expression_aligned.rds"))

train <- as.matrix(train)
gse91061 <- as.matrix(gse91061)
gse78220 <- as.matrix(gse78220)

rownames(train) <- normalize_gene_ids(rownames(train))
rownames(gse91061) <- normalize_gene_ids(rownames(gse91061))
rownames(gse78220) <- normalize_gene_ids(rownames(gse78220))

train <- train[!duplicated(rownames(train)), ]
gse91061 <- gse91061[!duplicated(rownames(gse91061)), ]
gse78220 <- gse78220[!duplicated(rownames(gse78220)), ]

common_genes <- Reduce(intersect, list(
  rownames(train),
  rownames(gse91061),
  rownames(gse78220)
))

cat("Número genes comunes:", length(common_genes), "\n")

if (length(common_genes) == 0) {
  stop("ERROR: 0 genes comunes → revisar IDs")
}

train_common <- train[common_genes, ]
gse91061_common <- gse91061[common_genes, ]
gse78220_common <- gse78220[common_genes, ]

saveRDS(train_common, file.path(data_processed_dir, "GSE160638_common_genes_aligned.rds"))
saveRDS(gse91061_common, file.path(data_processed_dir, "GSE91061_common_genes_aligned.rds"))
saveRDS(gse78220_common, file.path(data_processed_dir, "GSE78220_common_genes_aligned.rds"))

cat("OK: genes comunes definidos\n")

