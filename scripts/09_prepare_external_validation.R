################################
# 09_prepare_external_validation.R
################################

suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# Project root: run scripts either from the project root or from the scripts/ folder.
base_dir <- normalizePath(getwd(), mustWork = TRUE)
if (basename(base_dir) == "scripts") base_dir <- dirname(base_dir)
data_raw_dir <- file.path(base_dir, "data_raw")
data_processed_dir <- file.path(base_dir, "data_processed")

dir.create(data_raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_processed_dir, showWarnings = FALSE, recursive = TRUE)

clean_text <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "N/A", "na")] <- NA_character_
  x
}

clean_gene_ids <- function(x) {
  x <- clean_text(x)
  x <- sub("\\.\\d+$", "", x)
  x <- gsub("\\s+", "", x)
  toupper(x)
}

# =========================
# GSE91061 (CORRECCIÓ CRÍTICA)
# =========================
read_gse91061_expression <- function(dest_dir) {
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz"
  dest <- file.path(dest_dir, basename(url))
  
  if (!file.exists(dest)) {
    download.file(url, destfile = dest, mode = "wb")
  }
  
  x <- read_csv(dest, show_col_types = FALSE)
  x <- as.data.frame(x)
  
  gene_ids <- clean_text(x[[1]])
  x[[1]] <- NULL
  
  expr <- as.matrix(x)
  mode(expr) <- "numeric"
  
  # CONVERSIÓ A SYMBOL
  gene_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_ids,
    columns = "SYMBOL",
    keytype = "ENTREZID"
  )
  
  gene_map <- gene_map[!is.na(gene_map$SYMBOL), ]
  
  idx <- match(gene_ids, gene_map$ENTREZID)
  gene_symbols <- gene_map$SYMBOL[idx]
  
  keep <- !is.na(gene_symbols)
  expr <- expr[keep, ]
  gene_symbols <- gene_symbols[keep]
  
  rownames(expr) <- clean_gene_ids(gene_symbols)
  expr <- expr[!duplicated(rownames(expr)), ]
  
  expr
}

# =========================
# GSE78220
# =========================
read_gse78220_expression <- function(dest_dir) {
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/GSE78220_PatientFPKM.xlsx"
  dest <- file.path(dest_dir, basename(url))
  
  if (!file.exists(dest)) {
    download.file(url, destfile = dest, mode = "wb")
  }
  
  x <- readxl::read_xlsx(dest)
  x <- as.data.frame(x)
  
  gene_ids <- clean_gene_ids(x[[1]])
  x[[1]] <- NULL
  
  expr <- as.matrix(x)
  mode(expr) <- "numeric"
  
  rownames(expr) <- gene_ids
  expr <- expr[!duplicated(rownames(expr)), ]
  
  expr
}

# =========================
# DESAMENT FINAL
# =========================
datasets <- c("GSE91061", "GSE78220")

for (ds in datasets) {
  cat("\nProcessant:", ds, "\n")
  
  expr <- if (ds == "GSE91061") {
    read_gse91061_expression(data_raw_dir)
  } else {
    read_gse78220_expression(data_raw_dir)
  }
  
  saveRDS(
    expr,
    file.path(data_processed_dir, paste0(ds, "_expression_aligned.rds"))
  )
  
  cat("Desat:", ds, "\n")
}

