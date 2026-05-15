# ============================================================
# 05_differential_expression_exploratory.R
# Cleaned/renamed version of: 04_differential_expression.R
# Purpose: reproducible TFM pipeline while preserving the output filenames used in the memoria.
# ============================================================

################################
# 04_differential_expression.R
################################

suppressPackageStartupMessages({
  library(edgeR)
})

# Project root: run scripts either from the project root or from the scripts/ folder.
base_dir <- normalizePath(getwd(), mustWork = TRUE)
if (basename(base_dir) == "scripts") base_dir <- dirname(base_dir)

data_processed_dir <- file.path(base_dir, "data_processed")
results_dir <- file.path(base_dir, "results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

counts <- readRDS(
  file.path(data_processed_dir, "GSE160638_raw_counts_aligned.rds")
)

metadata <- read.csv(
  file.path(data_processed_dir, "metadata_GSE160638_aligned.csv"),
  stringsAsFactors = FALSE
)

colnames(counts) <- trimws(colnames(counts))
metadata$sample_name <- trimws(metadata$sample_name)

idx <- match(colnames(counts), metadata$sample_name)

cat("Número de muestras no alineadas:", sum(is.na(idx)), "\n")

if (any(is.na(idx))) {
  cat("\nMuestras de counts sin match en metadata:\n")
  print(colnames(counts)[is.na(idx)])
  stop("No se ha podido alinear metadata$sample_name con colnames(counts)")
}

metadata2 <- metadata[idx, , drop = FALSE]
stopifnot(all(metadata2$sample_name == colnames(counts)))

group <- factor(metadata2$response, levels = c("NonResponder", "Responder"))

cat("\nTabla response después del filtrado:\n")
print(table(group))

if (length(levels(group)) < 2) {
  stop("El factor 'group' tiene menos de 2 niveles. Revisa metadata$response.")
}

dge <- DGEList(counts = counts, group = group)

keep_genes <- filterByExpr(dge, group = group)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

design <- model.matrix(~ group)
dge <- estimateDisp(dge, design)

fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = 2)

de <- topTags(qlf, n = Inf)$table
de$gene <- rownames(de)

write.csv(
  de,
  file.path(results_dir, "differential_expression_GSE160638.csv"),
  row.names = FALSE
)

cat("\nAnálisis DE completado correctamente.\n")
cat("Número de genes analizados:", nrow(de), "\n")
cat("\nNOTA: este análisis DE es descriptivo/biológico.\n")
cat("No debe utilizarse para seleccionar genes antes de la validación cruzada del modelo.\n")

