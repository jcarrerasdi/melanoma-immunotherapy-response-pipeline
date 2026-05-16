################################
# 13_gene_stability_across_datasets.R
################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# Project root: run scripts either from the project root or from the scripts/ folder.
base_dir <- normalizePath(getwd(), mustWork = TRUE)
if (basename(base_dir) == "scripts") base_dir <- dirname(base_dir)
results_dir <- file.path(base_dir, "results")
data_processed_dir <- file.path(base_dir, "data_processed")
figures_dir <- file.path(base_dir, "figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

normalize_gene_ids <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.\\d+$", "", x)
  x <- gsub("\\s+", "", x)
  toupper(x)
}

convert_entrez_to_symbol_if_needed <- function(expr_mat, dataset_id) {
  expr_mat <- as.matrix(expr_mat)
  rn <- rownames(expr_mat)
  numeric_like <- grepl("^[0-9]+$", rn)
  prop_numeric <- mean(numeric_like, na.rm = TRUE)
  
  message("[", dataset_id, "] Proporció d'IDs numèrics: ", round(prop_numeric, 3))
  
  # Si ja són SYMBOL, no convertir
  if (prop_numeric <= 0.5) {
    message("[", dataset_id, "] Els IDs ja semblen SYMBOL. No s'aplica la conversió Entrez -> SYMBOL.")
    rownames(expr_mat) <- normalize_gene_ids(rownames(expr_mat))
    expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]
    return(expr_mat)
  }
  
  message("[", dataset_id, "] Entrez -> SYMBOL")
  
  gene_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = rn[numeric_like],
    columns = "SYMBOL",
    keytype = "ENTREZID"
  )
  
  gene_map <- gene_map[!is.na(gene_map$SYMBOL), , drop = FALSE]
  gene_map <- gene_map[!duplicated(gene_map$ENTREZID), , drop = FALSE]
  
  mapped <- match(rn, gene_map$ENTREZID)
  symbols <- gene_map$SYMBOL[mapped]
  
  keep <- !is.na(symbols)
  
  if (sum(keep) == 0) {
    stop("[", dataset_id, "] No s'ha pogut convertir cap gen Entrez a SYMBOL.")
  }
  
  expr_mat <- expr_mat[keep, , drop = FALSE]
  rownames(expr_mat) <- normalize_gene_ids(symbols[keep])
  
  expr_df <- as.data.frame(expr_mat, check.names = FALSE)
  expr_df$SYMBOL <- rownames(expr_mat)
  
  expr_df <- expr_df %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(
      dplyr::across(where(is.numeric), mean, na.rm = TRUE),
      .groups = "drop"
    )
  
  mat <- as.matrix(expr_df[, setdiff(colnames(expr_df), "SYMBOL"), drop = FALSE])
  rownames(mat) <- expr_df$SYMBOL
  
  return(mat)
}

# -----------------------------
# 1. Gens estables del model (≥4 folds)
# -----------------------------
stable_genes <- read_csv(
  file.path(results_dir, "stable_genes_4plus_final_100.csv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(gene_symbol = normalize_gene_ids(gene))

stable_list <- unique(stable_genes$gene_symbol)

# -----------------------------
# 2. Gens presents en datasets externs
# -----------------------------
expr_78220 <- readRDS(file.path(data_processed_dir, "GSE78220_expression_aligned.rds"))
expr_91061 <- readRDS(file.path(data_processed_dir, "GSE91061_expression_aligned.rds"))

expr_78220 <- convert_entrez_to_symbol_if_needed(expr_78220, "GSE78220")
expr_91061 <- convert_entrez_to_symbol_if_needed(expr_91061, "GSE91061")

genes_78220 <- normalize_gene_ids(rownames(expr_78220))
genes_91061 <- normalize_gene_ids(rownames(expr_91061))

# -----------------------------
# 3. Construir taula detallada
# -----------------------------
stability_table <- data.frame(
  gene = stable_list,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    in_GSE78220 = gene %in% genes_78220,
    in_GSE91061 = gene %in% genes_91061,
    present_in_both = in_GSE78220 & in_GSE91061
  ) %>%
  dplyr::arrange(dplyr::desc(present_in_both), dplyr::desc(in_GSE78220), dplyr::desc(in_GSE91061), gene)

# -----------------------------
# 4. Taula resum
# -----------------------------
summary_table <- data.frame(
  category = c(
    "Stable genes (>=4 folds)",
    "Present in GSE78220",
    "Present in GSE91061",
    "Present in both datasets"
  ),
  n_genes = c(
    nrow(stability_table),
    sum(stability_table$in_GSE78220),
    sum(stability_table$in_GSE91061),
    sum(stability_table$present_in_both)
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 5. Desar taules
# -----------------------------
write_csv(
  stability_table,
  file.path(results_dir, "gene_stability_across_datasets_final_100.csv")
)

write_csv(
  summary_table,
  file.path(results_dir, "gene_stability_across_datasets_summary_final_100.csv")
)

# -----------------------------
# 6. Gràfic resum de presència gènica
# -----------------------------
p_presence <- ggplot(summary_table, aes(x = reorder(category, n_genes), y = n_genes)) +
  geom_col() +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Estabilitat dels gens entre datasets externs",
    x = "Categoria",
    y = "Nombre de gens"
  )

ggsave(
  filename = file.path(figures_dir, "gene_stability_across_datasets_final_100.png"),
  plot = p_presence,
  width = 8,
  height = 5,
  dpi = 300
)

# -----------------------------
# 7. NOU: figura tipus paper per a 5.7
# -----------------------------
metrics_file <- file.path(results_dir, "stable_signature_external_validation_metrics.csv")

if (file.exists(metrics_file)) {
  
  auc_df <- read.csv(metrics_file, stringsAsFactors = FALSE)
  
  auc_df <- auc_df %>%
    dplyr::mutate(
      signature_label = dplyr::case_when(
        signature == "current_100" ~ "Current (100)",
        signature == "A" ~ "Stable A (≥3 folds)",
        signature == "B" ~ "Stable B (≥4 folds)",
        signature == "C" ~ "Stable C (ranking)",
        TRUE ~ signature
      )
    )
  
  auc_df$signature_label <- factor(
    auc_df$signature_label,
    levels = c(
      "Current (100)",
      "Stable A (≥3 folds)",
      "Stable B (≥4 folds)",
      "Stable C (ranking)"
    )
  )
  
  # taula neta per a l'informe
  auc_pretty <- auc_df %>%
    dplyr::select(dataset, signature_label, n_samples, n_genes_used, AUC) %>%
    dplyr::arrange(dataset, signature_label)
  
  write.csv(
    auc_pretty,
    file.path(results_dir, "stable_signature_external_validation_metrics_pretty.csv"),
    row.names = FALSE
  )
  
  # figura tipus paper
  p_auc <- ggplot(auc_df, aes(x = signature_label, y = AUC, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    theme_minimal(base_size = 12) +
    labs(
      title = "External validation performance by gene signature",
      x = "Gene signature",
      y = "AUC",
      fill = "Dataset"
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  
  ggsave(
    filename = file.path(figures_dir, "stable_signature_external_auc.png"),
    plot = p_auc,
    width = 8.5,
    height = 5.2,
    dpi = 300
  )
  
  cat("\nFigura tipus paper generada: stable_signature_external_auc.png\n")
  
} else {
  cat("\nNo s'ha trobat stable_signature_external_validation_metrics.csv\n")
  cat("Executa abans el script 10B_external_validation_stable_signatures.R\n")
}

# -----------------------------
# 8. Mostrar resum a la consola
# -----------------------------
cat("\nTaula resum d'estabilitat entre datasets:\n")
print(summary_table)

cat("\nScript completat correctament.\n")

