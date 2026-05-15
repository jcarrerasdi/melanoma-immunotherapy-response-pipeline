# ============================================================
# 03_preprocessing_train.R
# Cleaned/renamed version of: 02_preprocessing.R
# Purpose: reproducible TFM pipeline while preserving the output filenames used in the memoria.
# ============================================================

################################
# 02_preprocessing.R
################################

suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
})

# Project root: run scripts either from the project root or from the scripts/ folder.
base_dir <- normalizePath(getwd(), mustWork = TRUE)
if (basename(base_dir) == "scripts") base_dir <- dirname(base_dir)

data_processed_dir <- file.path(base_dir, "data_processed")

counts <- readRDS(
  file.path(data_processed_dir, "GSE160638_raw_counts.rds")
)

metadata <- read.csv(
  file.path(data_processed_dir, "metadata_GSE160638.csv"),
  stringsAsFactors = FALSE
)

# ===========================
# 0. Limpieza básica
# ===========================

colnames(counts) <- trimws(colnames(counts))

metadata[] <- lapply(metadata, function(x) {
  if (is.character(x)) trimws(x) else x
})

if (!"sample_name" %in% colnames(metadata)) {
  stop("La columna 'sample_name' no existe en metadata")
}

metadata$sample_name <- trimws(metadata$sample_name)

# ===========================
# 1. Alineación counts-metadata
# ===========================

idx <- match(colnames(counts), metadata$sample_name)

cat("Número de muestras en counts:", ncol(counts), "\n")
cat("Número de muestras en metadata:", nrow(metadata), "\n")
cat("Número de muestras no alineadas:", sum(is.na(idx)), "\n")

if (any(is.na(idx))) {
  cat("\nMuestras de counts sin match en metadata:\n")
  print(colnames(counts)[is.na(idx)])
  stop("No se ha podido alinear metadata con counts en preprocessing")
}

metadata <- metadata[idx, , drop = FALSE]

stopifnot(all(metadata$sample_name == colnames(counts)))

# ===========================
# 2. Control explícito de pre-treatment
# ===========================
# Este bloque es crítico para evitar mezclar muestras baseline con
# muestras on-treatment o post-treatment.
#
# Si existe una variable temporal, se verifica y se filtra de forma explícita.
# Si no existe, el script informa de ello y mantiene las muestras, asumiendo
# que corresponden al diseño basal del dataset.

timepoint_candidates <- c(
  "timepoint",
  "Timepoint",
  "time_point",
  "sample_timepoint",
  "treatment_phase",
  "treatment_time",
  "collection_timepoint"
)

timepoint_col <- intersect(timepoint_candidates, colnames(metadata))

if (length(timepoint_col) > 0) {
  timepoint_col <- timepoint_col[1]
  
  cat("\nColumna temporal detectada:", timepoint_col, "\n")
  cat("Distribución original de timepoint:\n")
  print(table(metadata[[timepoint_col]], useNA = "ifany"))
  
  timepoint_values <- tolower(trimws(as.character(metadata[[timepoint_col]])))
  
  pre_labels <- c(
    "pre", "pretreatment", "pre-treatment", "pre treatment",
    "baseline", "before treatment", "prior to treatment",
    "screening", "week0", "week 0", "day0", "day 0"
  )
  
  keep_pre <- !is.na(timepoint_values) & timepoint_values %in% pre_labels
  
  cat("\nNúmero de muestras pre-treatment detectadas:", sum(keep_pre), "\n")
  cat("Número de muestras excluidas por no ser pre-treatment:", sum(!keep_pre), "\n")
  
  counts <- counts[, keep_pre, drop = FALSE]
  metadata <- metadata[keep_pre, , drop = FALSE]
  
  cat("\nDistribución de timepoint tras el filtrado:\n")
  print(table(metadata[[timepoint_col]], useNA = "ifany"))
  
} else {
  cat("\nNo se ha detectado columna temporal explícita en metadata.\n")
  cat("Se asume que las muestras corresponden al diseño basal/pre-treatment del dataset.\n")
}

# Revalidar alineación tras posible filtrado temporal
stopifnot(all(metadata$sample_name == colnames(counts)))

# ===========================
# 3. Filtrado por respuesta clínica válida
# ===========================

if (!"response" %in% colnames(metadata)) {
  stop("La columna 'response' no existe en metadata")
}

keep_samples <- !is.na(metadata$response) & metadata$response != ""

counts <- counts[, keep_samples, drop = FALSE]
metadata <- metadata[keep_samples, , drop = FALSE]

group <- factor(metadata$response, levels = c("NonResponder", "Responder"))

cat("\nTabla response:\n")
print(table(group, useNA = "ifany"))

if (nlevels(droplevels(group)) < 2) {
  stop("Después del filtrado no quedan dos clases válidas en 'response'")
}

# ===========================
# 4. Crear objeto DGE
# ===========================

dge <- DGEList(counts = counts, group = group)

# ===========================
# 5. Filtrar genes poco expresados
# ===========================

keep <- filterByExpr(dge, group = group)
dge <- dge[keep, , keep.lib.sizes = FALSE]

cat("\nDimensiones después del filtrado de genes:", dim(dge), "\n")

# ===========================
# 6. Normalización TMM
# ===========================

dge <- calcNormFactors(dge)

# ===========================
# 7. Expresión logCPM
# ===========================

expr <- cpm(dge, log = TRUE, prior.count = 1)

cat("Dimensiones expr:", dim(expr), "\n")
cat("Primeros genes:\n")
print(head(rownames(expr)))

# ===========================
# 8. Guardado de resultados
# ===========================

# Matriz para PCA / exploración descriptiva
saveRDS(
  expr,
  file.path(data_processed_dir, "GSE160638_expression_matrix_exploration.rds")
)

# Counts alineados y filtrados para modelización sin fugas
saveRDS(
  counts,
  file.path(data_processed_dir, "GSE160638_raw_counts_aligned.rds")
)

write.csv(
  metadata,
  file.path(data_processed_dir, "metadata_GSE160638_aligned.csv"),
  row.names = FALSE
)

cat("\nPreprocessament exploratori completat.\n")
cat("NOTA: aquesta matriu d'expressió és per PCA i anàlisi descriptiva,\n")
cat("no per a la validació del model predictiu.\n")

