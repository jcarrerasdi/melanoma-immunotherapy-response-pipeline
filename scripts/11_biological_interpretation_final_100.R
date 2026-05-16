################################
# 11_biological_interpretation_final_100.R
################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(pheatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

has_reactome <- requireNamespace("ReactomePA", quietly = TRUE)
if (has_reactome) {
  library(ReactomePA)
}

# Project root: run scripts either from the project root or from the scripts/ folder.
base_dir <- normalizePath(getwd(), mustWork = TRUE)
if (basename(base_dir) == "scripts") base_dir <- dirname(base_dir)
data_processed_dir <- file.path(base_dir, "data_processed")
results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(base_dir, "figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

expr_matrix <- readRDS(file.path(data_processed_dir, "GSE160638_expression_matrix_exploration.rds"))
metadata <- read_csv(file.path(data_processed_dir, "metadata_GSE160638_aligned.csv"), show_col_types = FALSE)

importance_clean <- read_csv(
  file.path(results_dir, "rf_variable_importance_clean_final_100_exploratory.csv"),
  show_col_types = FALSE
)

stable_4plus <- read.csv(
  file.path(results_dir, "stable_genes_4plus_final_100.csv"),
  stringsAsFactors = FALSE
)

stable_5of5 <- read.csv(
  file.path(results_dir, "stable_genes_5of5_final_100.csv"),
  stringsAsFactors = FALSE
)

rownames(expr_matrix) <- gsub("\\.", "-", rownames(expr_matrix))

importance_clean <- importance_clean %>%
  mutate(gene_symbol = gsub("\\.", "-", gene))

stable_4plus <- stable_4plus %>%
  mutate(gene_symbol = gsub("\\.", "-", gene))

stable_5of5 <- stable_5of5 %>%
  mutate(gene_symbol = gsub("\\.", "-", gene))

n_genes_match_rows <- sum(importance_clean$gene_symbol %in% rownames(expr_matrix))
n_genes_match_cols <- sum(importance_clean$gene_symbol %in% colnames(expr_matrix))

if (n_genes_match_cols > n_genes_match_rows) {
  expr_matrix <- t(expr_matrix)
  rownames(expr_matrix) <- gsub("\\.", "-", rownames(expr_matrix))
}

candidate_cols <- colnames(metadata)
match_counts <- sapply(candidate_cols, function(cn) {
  vals <- as.character(metadata[[cn]])
  sum(vals %in% colnames(expr_matrix), na.rm = TRUE)
})

best_col <- names(which.max(match_counts))
best_n <- max(match_counts)

if (best_n == 0) {
  stop("No s'ha trobat cap columna de metadata que coincideixi amb les mostres de expr_matrix.")
}

metadata <- metadata %>%
  mutate(sample_id = as.character(.data[[best_col]]))

if (!"response" %in% colnames(metadata)) {
  stop("La columna 'response' no existeix a metadata.")
}

metadata <- metadata %>%
  filter(sample_id %in% colnames(expr_matrix)) %>%
  distinct(sample_id, .keep_all = TRUE)

expr_matrix <- expr_matrix[, metadata$sample_id, drop = FALSE]

top20_genes <- importance_clean %>%
  slice_max(order_by = importance, n = 20) %>%
  pull(gene_symbol)

top20_genes_present <- top20_genes[top20_genes %in% rownames(expr_matrix)]

if (length(top20_genes_present) < 2) {
  stop("No hi ha prou gens del top20 presents a la matriu per generar el heatmap.")
}

expr_top20 <- expr_matrix[top20_genes_present, metadata$sample_id, drop = FALSE]
expr_top20 <- expr_top20[match(top20_genes_present, rownames(expr_top20)), , drop = FALSE]

expr_top20 <- as.matrix(expr_top20)
mode(expr_top20) <- "numeric"

expr_top20 <- expr_top20[
  apply(expr_top20, 1, function(x) all(is.finite(x))),
  ,
  drop = FALSE
]

expr_top20 <- expr_top20[
  apply(expr_top20, 1, function(x) sd(x, na.rm = TRUE) > 0),
  ,
  drop = FALSE
]

expr_top20_scaled <- t(scale(t(expr_top20)))
expr_top20_scaled[!is.finite(expr_top20_scaled)] <- 0

annotation_col <- data.frame(Response = as.character(metadata$response))
rownames(annotation_col) <- metadata$sample_id

annotation_col$Response <- trimws(tolower(annotation_col$Response))
annotation_col$Response[annotation_col$Response %in% c("r", "responder", "response", "yes", "y", "1")] <- "Responder"
annotation_col$Response[annotation_col$Response %in% c("nr", "non-responder", "nonresponder", "non responder", "no", "n", "0")] <- "Non-responder"
annotation_col$Response <- factor(annotation_col$Response, levels = c("Responder", "Non-responder"))

ann_colors <- list(
  Response = c("Responder" = "#1b9e77", "Non-responder" = "#d95f02")
)

png(
  filename = file.path(figures_dir, "heatmap_top20_genes_final_100.png"),
  width = 2200, height = 1400, res = 200
)

pheatmap(
  expr_top20_scaled,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  scale = "none",
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  fontsize_row = 11,
  fontsize_col = 7,
  show_colnames = FALSE,
  border_color = NA,
  main = "Heatmap dels 20 gens més importants del model final (100 gens)"
)

dev.off()

genes_for_enrichment <- stable_4plus$gene_symbol %>% unique()
genes_for_enrichment <- genes_for_enrichment[!is.na(genes_for_enrichment)]

gene_df <- bitr(
  genes_for_enrichment,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

if (is.null(gene_df) || nrow(gene_df) == 0) {
  stop("No s'han pogut convertir gens a ENTREZID per a l'anàlisi d'enriquiment.")
}

genes_entrez <- unique(gene_df$ENTREZID)

ego_bp <- enrichGO(
  gene = genes_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

ego_mf <- enrichGO(
  gene = genes_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

ekegg <- enrichKEGG(
  gene = genes_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05
)

ereact <- NULL
if (has_reactome) {
  ereact <- ReactomePA::enrichPathway(
    gene = genes_entrez,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = TRUE
  )
}

write_csv(as.data.frame(ego_bp), file.path(results_dir, "GO_BP_enrichment_stable_4plus_final_100.csv"))
write_csv(as.data.frame(ego_mf), file.path(results_dir, "GO_MF_enrichment_stable_4plus_final_100.csv"))
write_csv(as.data.frame(ekegg), file.path(results_dir, "KEGG_enrichment_stable_4plus_final_100.csv"))

if (!is.null(ereact)) {
  write_csv(as.data.frame(ereact), file.path(results_dir, "Reactome_enrichment_stable_4plus_final_100.csv"))
}

if (nrow(as.data.frame(ego_bp)) > 0) {
  png(file.path(figures_dir, "dotplot_GO_BP_stable_4plus_final_100.png"), width = 2200, height = 1400, res = 200)
  print(dotplot(ego_bp, showCategory = 12, font.size = 12, title = "GO Biological Process"))
  dev.off()
}

if (nrow(as.data.frame(ego_mf)) > 0) {
  png(file.path(figures_dir, "dotplot_GO_MF_stable_4plus_final_100.png"), width = 2200, height = 1400, res = 200)
  print(dotplot(ego_mf, showCategory = 12, font.size = 12, title = "GO Molecular Function"))
  dev.off()
}

if (nrow(as.data.frame(ekegg)) > 0) {
  png(file.path(figures_dir, "dotplot_KEGG_stable_4plus_final_100.png"), width = 2200, height = 1400, res = 200)
  print(dotplot(ekegg, showCategory = 12, font.size = 12, title = "KEGG pathways"))
  dev.off()
}

if (!is.null(ereact) && nrow(as.data.frame(ereact)) > 0) {
  png(file.path(figures_dir, "dotplot_Reactome_stable_4plus_final_100.png"), width = 2200, height = 1400, res = 200)
  print(dotplot(ereact, showCategory = 12, font.size = 12, title = "Reactome pathways"))
  dev.off()
}

top20_table <- importance_clean %>%
  dplyr::slice_max(order_by = importance, n = 20) %>%
  dplyr::mutate(rank = dplyr::row_number()) %>%
  dplyr::select(rank, gene, gene_symbol, importance)

readr::write_csv(
  top20_table,
  file.path(results_dir, "top20_genes_final_100_for_report.csv")
)

cat("\nAnàlisi biològica final_100 completada correctament.\n")

