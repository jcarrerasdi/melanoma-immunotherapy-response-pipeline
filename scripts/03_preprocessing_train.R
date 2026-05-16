################################
# 03_preprocessing.R
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
# 0. Neteja bàsica
# ===========================

colnames(counts) <- trimws(colnames(counts))

metadata[] <- lapply(metadata, function(x) {
  if (is.character(x)) trimws(x) else x
})

if (!"sample_name" %in% colnames(metadata)) {
  stop("La columna 'sample_name' no existeix a metadata")
}

metadata$sample_name <- trimws(metadata$sample_name)

# ===========================
# 1. Alineació counts-metadata
# ===========================

idx <- match(colnames(counts), metadata$sample_name)

cat("Nombre de mostres en counts:", ncol(counts), "\n")
cat("Nombre de mostres en metadata:", nrow(metadata), "\n")
cat("Nombre de mostres no alineades:", sum(is.na(idx)), "\n")

if (any(is.na(idx))) {
  cat("\nMostres de counts sense coincidència en metadata:\n")
  print(colnames(counts)[is.na(idx)])
  stop("No s'ha pogut alinear metadata amb counts en el preprocessament")
}

metadata <- metadata[idx, , drop = FALSE]

stopifnot(all(metadata$sample_name == colnames(counts)))

# ===========================
# 2. Control explícit de pre-treatment
# ===========================
# Aquest bloc és crític per evitar barrejar mostres baseline amb
# mostres on-treatment o post-treatment.
#
# Si existeix una variable temporal, es verifica i es filtra de forma explícita.
# Si no existeix, el script n'informa i manté les mostres, assumint
# que corresponen al disseny basal del dataset.

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
  cat("Distribució original de timepoint:\n")
  print(table(metadata[[timepoint_col]], useNA = "ifany"))
  
  timepoint_values <- tolower(trimws(as.character(metadata[[timepoint_col]])))
  
  pre_labels <- c(
    "pre", "pretreatment", "pre-treatment", "pre treatment",
    "baseline", "before treatment", "prior to treatment",
    "screening", "week0", "week 0", "day0", "day 0"
  )
  
  keep_pre <- !is.na(timepoint_values) & timepoint_values %in% pre_labels
  
  cat("\nNombre de mostres pre-treatment detectades:", sum(keep_pre), "\n")
  cat("Nombre de mostres excloses per no ser pre-treatment:", sum(!keep_pre), "\n")
  
  counts <- counts[, keep_pre, drop = FALSE]
  metadata <- metadata[keep_pre, , drop = FALSE]
  
  cat("\nDistribució de timepoint després del filtratge:\n")
  print(table(metadata[[timepoint_col]], useNA = "ifany"))
  
} else {
  cat("\nNo s'ha detectat cap columna temporal explícita a metadata.\n")
  cat("S'assumeix que les mostres corresponen al disseny basal/pre-treatment del dataset.\n")
}

# Revalidar l'alineació després del possible filtratge temporal
stopifnot(all(metadata$sample_name == colnames(counts)))

# ===========================
# 3. Filtratge per resposta clínica vàlida
# ===========================

if (!"response" %in% colnames(metadata)) {
  stop("La columna 'response' no existeix a metadata")
}

keep_samples <- !is.na(metadata$response) & metadata$response != ""

counts <- counts[, keep_samples, drop = FALSE]
metadata <- metadata[keep_samples, , drop = FALSE]

group <- factor(metadata$response, levels = c("NonResponder", "Responder"))

cat("\nTaula response:\n")
print(table(group, useNA = "ifany"))

if (nlevels(droplevels(group)) < 2) {
  stop("Després del filtratge no queden dues classes vàlides a 'response'")
}

# ===========================
# 4. Crear objecte DGE
# ===========================

dge <- DGEList(counts = counts, group = group)

# ===========================
# 5. Filtrar gens poc expressats
# ===========================

keep <- filterByExpr(dge, group = group)
dge <- dge[keep, , keep.lib.sizes = FALSE]

cat("\nDimensions després del filtratge de gens:", dim(dge), "\n")

# ===========================
# 6. Normalització TMM
# ===========================

dge <- calcNormFactors(dge)

# ===========================
# 7. Expressió logCPM
# ===========================

expr <- cpm(dge, log = TRUE, prior.count = 1)

cat("Dimensions expr:", dim(expr), "\n")
cat("Primers gens:\n")
print(head(rownames(expr)))

# ===========================
# 8. Desament de resultats
# ===========================

# Matriu per a PCA / exploració descriptiva
saveRDS(
  expr,
  file.path(data_processed_dir, "GSE160638_expression_matrix_exploration.rds")
)

# Counts alineats i filtrats per a modelització sense fugues
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

