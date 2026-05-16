################################
# 04_exploratory_analysis.R
################################

suppressPackageStartupMessages({
  library(ggplot2)
})

# Project root: run scripts either from the project root or from the scripts/ folder.
base_dir <- normalizePath(getwd(), mustWork = TRUE)
if (basename(base_dir) == "scripts") base_dir <- dirname(base_dir)

data_processed_dir <- file.path(base_dir, "data_processed")
figures_dir <- file.path(base_dir, "figures")

dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

expr <- readRDS(
  file.path(data_processed_dir, "GSE160638_expression_matrix_exploration.rds")
)

metadata <- read.csv(
  file.path(data_processed_dir, "metadata_GSE160638_aligned.csv"),
  stringsAsFactors = FALSE
)

colnames(expr) <- trimws(colnames(expr))
metadata$sample_name <- trimws(metadata$sample_name)

idx <- match(colnames(expr), metadata$sample_name)

cat("Nombre de mostres a expr:", ncol(expr), "\n")
cat("Nombre de mostres a metadata:", nrow(metadata), "\n")
cat("Nombre de mostres no alineades:", sum(is.na(idx)), "\n")

if (any(is.na(idx))) {
  cat("\nMostres d'expr sense coincidència a metadata:\n")
  print(colnames(expr)[is.na(idx)])
  stop("No s'ha pogut alinear metadata$sample_name amb colnames(expr)")
}

metadata <- metadata[idx, , drop = FALSE]

stopifnot(all(metadata$sample_name == colnames(expr)))

pca <- prcomp(t(expr), scale. = FALSE)

pca_df <- data.frame(
  sample_name = colnames(expr),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  response = metadata$response
)

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = response)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA of GSE160638 expression data",
    x = "PC1",
    y = "PC2",
    color = "Response"
  )

print(p)

ggsave(
  filename = file.path(figures_dir, "PCA_GSE160638.png"),
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

