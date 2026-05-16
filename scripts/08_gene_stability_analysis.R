################################
# 08_gene_stability_analysis.R
################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# Project root: run scripts either from the project root or from the scripts/ folder.
base_dir <- normalizePath(getwd(), mustWork = TRUE)
if (basename(base_dir) == "scripts") base_dir <- dirname(base_dir)

results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(base_dir, "figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ===========================
# 1. Llegir fitxers
# ===========================

gene_freq <- read.csv(
  file.path(results_dir, "rf_gene_selection_frequency_final_100_no_leakage.csv"),
  stringsAsFactors = FALSE
)

imp <- read.csv(
  file.path(results_dir, "rf_variable_importance_clean_final_100_exploratory.csv"),
  stringsAsFactors = FALSE
)

# ===========================
# 2. Resum global
# ===========================

stability_summary <- gene_freq |>
  count(n_folds_selected, name = "n_genes") |>
  arrange(n_folds_selected)

write.csv(
  stability_summary,
  file.path(results_dir, "gene_stability_summary.csv"),
  row.names = FALSE
)

cat("\nResum d'estabilitat:\n")
print(stability_summary)

# ===========================
# 3. Genes estables
# ===========================

genes_5of5 <- gene_freq |>
  filter(n_folds_selected == 5) |>
  arrange(gene)

genes_4plus <- gene_freq |>
  filter(n_folds_selected >= 4) |>
  arrange(desc(n_folds_selected), gene)

write.csv(
  genes_5of5,
  file.path(results_dir, "genes_selected_in_5_of_5_folds.csv"),
  row.names = FALSE
)

write.csv(
  genes_4plus,
  file.path(results_dir, "genes_selected_in_4plus_folds.csv"),
  row.names = FALSE
)

cat("\nNombre de gens en 5/5 folds:", nrow(genes_5of5), "\n")
cat("Nombre de gens en >=4 folds:", nrow(genes_4plus), "\n")

# ===========================
# 4. Creuar estabilitat + importància
# ===========================

stable_imp <- gene_freq |>
  inner_join(imp, by = "gene") |>
  arrange(desc(n_folds_selected), desc(importance))

write.csv(
  stable_imp,
  file.path(results_dir, "stable_genes_with_importance.csv"),
  row.names = FALSE
)

stable_imp_5of5 <- stable_imp |>
  filter(n_folds_selected == 5)

stable_imp_4plus <- stable_imp |>
  filter(n_folds_selected >= 4)

write.csv(
  stable_imp_5of5,
  file.path(results_dir, "stable_genes_5of5_with_importance.csv"),
  row.names = FALSE
)

write.csv(
  stable_imp_4plus,
  file.path(results_dir, "stable_genes_4plus_with_importance.csv"),
  row.names = FALSE
)

cat("\nTop gens estables + importants:\n")
print(head(stable_imp_4plus, 20))

# ===========================
# 5. Figura: distribució d’estabilitat
# ===========================

p1 <- ggplot(stability_summary, aes(x = factor(n_folds_selected), y = n_genes)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Distribució d’estabilitat gènica entre folds",
    x = "Nombre de folds en què va ser seleccionat",
    y = "Nombre de gens"
  )

ggsave(
  filename = file.path(figures_dir, "gene_stability_distribution.png"),
  plot = p1,
  width = 8,
  height = 6,
  dpi = 300
)

# ===========================
# 6. Figura: top gens estables + importants
# ===========================

top_stable_imp <- stable_imp_4plus |>
  head(20)

p2 <- ggplot(top_stable_imp, aes(x = reorder(gene, importance), y = importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top gens estables (>=4 folds) amb major importància exploratòria",
    x = "Gen",
    y = "Mean Decrease Accuracy"
  )

ggsave(
  filename = file.path(figures_dir, "top_stable_important_genes.png"),
  plot = p2,
  width = 9,
  height = 7,
  dpi = 300
)

cat("\nAnàlisi d’estabilitat completada.\n")

