################################
# 02_clinical_harmonization.R
################################

suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

# --------------------------------------------------
# RUTES
# --------------------------------------------------

# Project root: run scripts either from the project root or from the scripts/ folder.
base_dir <- normalizePath(getwd(), mustWork = TRUE)
if (basename(base_dir) == "scripts") base_dir <- dirname(base_dir)
data_processed_dir <- file.path(base_dir, "data_processed")
results_dir <- file.path(base_dir, "results")

dir.create(data_processed_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

datasets <- c("GSE160638", "GSE78220", "GSE91061")

# --------------------------------------------------
# FUNCIONS AUXILIARS
# --------------------------------------------------

clean_text <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "N/A", "na", "n/a", "NULL", "null")] <- NA_character_
  x
}

collapse_characteristics <- function(pheno) {
  char_cols <- grep("^characteristics_ch1", colnames(pheno), value = TRUE)
  if (length(char_cols) == 0) return(rep(NA_character_, nrow(pheno)))
  
  apply(pheno[, char_cols, drop = FALSE], 1, function(z) {
    z <- clean_text(z)
    z <- z[!is.na(z)]
    paste(z, collapse = " | ")
  })
}

extract_sample_name_gse160638 <- function(title_vec) {
  out <- sub("^.*\\(([^()]*)\\)\\s*$", "\\1", clean_text(title_vec))
  out <- ifelse(is.na(out), clean_text(title_vec), out)
  trimws(out)
}

extract_sample_name_default <- function(pheno, dataset_id) {
  if (dataset_id == "GSE160638") {
    return(extract_sample_name_gse160638(pheno$title))
  }
  
  if (dataset_id %in% c("GSE78220", "GSE91061") && "title" %in% colnames(pheno)) {
    return(clean_text(pheno$title))
  }
  
  if ("geo_accession" %in% colnames(pheno)) {
    return(clean_text(pheno$geo_accession))
  }
  
  stop(paste0("[", dataset_id, "] no s'ha pogut inferir sample_name"))
}

# --------------------------------------------------
# RESPOSTA CLÍNICA
# --------------------------------------------------

extract_response_from_text <- function(x) {
  x0 <- tolower(clean_text(x))
  out <- rep(NA_character_, length(x0))
  
  out[grepl("\\bcomplete response\\b|\\bpartial response\\b|\\bcr\\b|\\bpr\\b", x0)] <- "Responder"
  out[grepl("\\bprogressive disease\\b|\\bstable disease\\b|\\bpd\\b|\\bsd\\b", x0)] <- "NonResponder"
  
  out[grepl("\\bresponder\\b", x0) & !grepl("non[- ]?responder", x0)] <- "Responder"
  out[grepl("non[- ]?responder|\\bnr\\b", x0)] <- "NonResponder"
  
  out[grepl("clinical benefit[:= ]?yes|benefit[:= ]?yes|cb[:= ]?yes", x0)] <- "Responder"
  out[grepl("clinical benefit[:= ]?no|benefit[:= ]?no|cb[:= ]?no", x0)] <- "NonResponder"
  
  out
}

extract_exact_response_value <- function(x) {
  x <- clean_text(x)
  out <- rep(NA_character_, length(x))
  
  pat <- "(?i)(?:anti-pd-1 response|response)\\s*[:=]\\s*([^|;]+)"
  has_pat <- grepl(pat, x, perl = TRUE)
  out[has_pat] <- stringr::str_match(x[has_pat], pat)[, 2]
  
  clean_text(out)
}

map_response_value <- function(x) {
  x0 <- tolower(clean_text(x))
  x0 <- gsub("\\s+", " ", x0)
  
  out <- rep(NA_character_, length(x0))
  
  out[x0 %in% c(
    "complete response", "partial response", "cr", "pr", "prcr"
  )] <- "Responder"
  
  out[x0 %in% c(
    "progressive disease", "stable disease", "pd", "sd", "nr",
    "non-response", "non response"
  )] <- "NonResponder"
  
  out
}

# --------------------------------------------------
# TEMPS DE TRACTAMENT
# --------------------------------------------------

infer_treatment_time <- function(dataset_id, title_vec, source_vec, char_all, sample_name) {
  text_all <- paste(
    clean_text(title_vec),
    clean_text(source_vec),
    clean_text(char_all),
    clean_text(sample_name),
    sep = " | "
  )
  
  text_all <- tolower(text_all)
  out <- rep(NA_character_, length(text_all))
  
  # regles generals
  out[grepl("pre[- ]?treatment|pretreatment|baseline|before treatment", text_all)] <- "PreTreatment"
  out[grepl("on[- ]?treatment|post[- ]?treatment|after treatment|progression", text_all)] <- "PostOrOnTreatment"
  
  # regles específiques
  if (dataset_id == "GSE160638") {
    out[is.na(out)] <- "PreTreatment"
  }
  
  if (dataset_id == "GSE78220") {
    out[is.na(out)] <- "PreTreatment"
  }
  
  if (dataset_id == "GSE91061") {
    out[grepl("_pre_", text_all)] <- "PreTreatment"
    out[grepl("visit (pre or on treatment): pre", text_all, fixed = TRUE)] <- "PreTreatment"
    
    out[grepl("_on_", text_all)] <- "PostOrOnTreatment"
    out[grepl("visit (pre or on treatment): on", text_all, fixed = TRUE)] <- "PostOrOnTreatment"
  }
  
  out
}

# --------------------------------------------------
# PROCESSAMENT
# --------------------------------------------------

all_metadata <- list()

for (ds in datasets) {
  
  cat("\nProcessant:", ds, "\n")
  
  gse <- getGEO(ds, GSEMatrix = TRUE)
  eset <- gse[[1]]
  pheno <- pData(eset) |> as.data.frame()
  
  pheno$char_all <- collapse_characteristics(pheno)
  sample_name <- extract_sample_name_default(pheno, ds)
  
  if (ds == "GSE160638") {
    meta_train <- read.csv(
      file.path(data_processed_dir, "metadata_GSE160638.csv"),
      stringsAsFactors = FALSE
    )
    
    meta_train$sample_name <- trimws(meta_train$sample_name)
    sample_name <- trimws(sample_name)
    
    m_idx <- match(sample_name, meta_train$sample_name)
    
    response_raw <- meta_train$response[m_idx]
    response <- response_raw
  } else {
    response_raw <- extract_exact_response_value(pheno$char_all)
    idx_na <- is.na(response_raw)
    
    response_raw[idx_na] <- extract_response_from_text(pheno$char_all)[idx_na]
    response <- map_response_value(response_raw)
    
    idx_na2 <- is.na(response)
    response[idx_na2] <- extract_response_from_text(pheno$char_all)[idx_na2]
  }
  
  treatment_time <- infer_treatment_time(
    ds,
    pheno$title,
    pheno$source_name_ch1,
    pheno$char_all,
    sample_name
  )
  
  meta <- tibble(
    dataset = ds,
    sample_id = clean_text(pheno$geo_accession),
    sample_name = sample_name,
    response_raw = response_raw,
    response = response,
    treatment_time = treatment_time
  )
  
  meta <- meta |>
    mutate(
      keep = !is.na(response) & treatment_time == "PreTreatment",
      reason_excluded = case_when(
        is.na(response) ~ "No response info",
        treatment_time != "PreTreatment" ~ "Not pre-treatment",
        TRUE ~ NA_character_
      )
    )
  
  write_csv(
    meta,
    file.path(data_processed_dir, paste0("metadata_", ds, "_harmonized_full.csv"))
  )
  
  write_csv(
    meta |> filter(keep),
    file.path(data_processed_dir, paste0("metadata_", ds, "_harmonized_kept.csv"))
  )
  
  write_csv(
    meta |> filter(!keep),
    file.path(data_processed_dir, paste0("metadata_", ds, "_harmonized_excluded.csv"))
  )
  
  all_metadata[[ds]] <- meta
}

# --------------------------------------------------
# RESUM
# --------------------------------------------------

all_df <- bind_rows(all_metadata)

summary_df <- all_df |>
  group_by(dataset) |>
  summarise(
    n_total = n(),
    n_kept = sum(keep),
    n_excluded = sum(!keep),
    .groups = "drop"
  )

write_csv(
  summary_df,
  file.path(results_dir, "samples_kept_summary_step1.csv")
)

write_csv(
  all_df,
  file.path(results_dir, "clinical_metadata_harmonization_audit_step1.csv")
)

write_csv(
  all_df |>
    count(dataset, response, treatment_time, keep, sort = TRUE),
  file.path(results_dir, "clinical_metadata_harmonization_summary_step1.csv")
)

write_csv(
  all_df |>
    filter(!keep) |>
    count(dataset, reason_excluded, sort = TRUE),
  file.path(results_dir, "clinical_metadata_exclusion_reasons_step1.csv")
)

cat("\nPAS 1 COMPLETAT\n")
print(summary_df)

