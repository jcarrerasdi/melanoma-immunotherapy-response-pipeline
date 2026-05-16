################################
# 01_load_data.R
################################

suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
})

# Project root: run scripts either from the project root or from the scripts/ folder.
base_dir <- normalizePath(getwd(), mustWork = TRUE)
if (basename(base_dir) == "scripts") base_dir <- dirname(base_dir)

data_raw_dir <- file.path(base_dir, "data_raw")
data_processed_dir <- file.path(base_dir, "data_processed")

dir.create(data_processed_dir, showWarnings = FALSE, recursive = TRUE)

# ===========================
# 1. Carregar counts
# ===========================

counts_path <- file.path(data_raw_dir, "GSE160638_combined_raw_counts.csv.gz")

counts <- read.csv(
  counts_path,
  row.names = 1,
  check.names = FALSE
)

colnames(counts) <- trimws(colnames(counts))
rownames(counts) <- trimws(rownames(counts))

cat("Dimensions counts:", dim(counts), "\n")
cat("Primeres files:\n")
print(head(rownames(counts)))
cat("Primeres mostres:\n")
print(head(colnames(counts)))

saveRDS(
  counts,
  file.path(data_processed_dir, "GSE160638_raw_counts.rds")
)

# ===========================
# 2. Descarregar metadata de GEO
# ===========================

gse <- getGEO("GSE160638", GSEMatrix = TRUE)
eset <- gse[[1]]

pheno <- pData(eset) |> as.data.frame()

write.csv(
  pheno,
  file.path(data_processed_dir, "GSE160638_pheno_FULL.csv"),
  row.names = FALSE
)

# ===========================
# 3. Construir metadata mínima
# ===========================

char_cols <- grep("^characteristics_ch1", colnames(pheno), value = TRUE)

pheno$char_all <- apply(
  pheno[, char_cols, drop = FALSE],
  1,
  paste,
  collapse = " | "
)

pheno$response_raw <- NA_character_

pheno$response_raw[
  grepl(
    "responder|response[:= ]?R\\b|benefit[:= ]?yes|clinical benefit[:= ]?yes|CB[:= ]?yes|CR|PR",
    pheno$char_all,
    ignore.case = TRUE
  )
] <- "Responder"

pheno$response_raw[
  grepl(
    "non[- ]?responder|no response|benefit[:= ]?no|no clinical benefit|PD|SD",
    pheno$char_all,
    ignore.case = TRUE
  )
] <- "NonResponder"

# Extreure identificador entre parèntesis final
sample_name_clean <- sub("^.*\\(([^()]*)\\)\\s*$", "\\1", trimws(pheno$title))

meta <- pheno |>
  transmute(
    sample_id = geo_accession,
    sample_name = sample_name_clean,
    title = trimws(title),
    source_name = source_name_ch1,
    characteristics = char_all,
    response = response_raw
  )

write.csv(
  meta,
  file.path(data_processed_dir, "metadata_GSE160638.csv"),
  row.names = FALSE
)

cat("\nMetadata generada correctament.\n")
cat("Taula response:\n")
print(table(meta$response, useNA = "ifany"))

cat("\nPrimeres sample_name netes:\n")
print(head(meta$sample_name))

cat("\nMostres de counts no presents en metadata:\n")
print(setdiff(colnames(counts), meta$sample_name))

cat("\nMostres de metadata no presents en counts:\n")
print(setdiff(meta$sample_name, colnames(counts)))

