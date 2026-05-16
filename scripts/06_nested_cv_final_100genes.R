################################
# 06_nested_cv_final_100genes.R
################################

suppressPackageStartupMessages({
  library(edgeR)
  library(randomForest)
  library(pROC)
  library(ggplot2)
  library(caret)
  library(dplyr)
})

# Project root: run scripts either from the project root or from the scripts/ folder.
base_dir <- normalizePath(getwd(), mustWork = TRUE)
if (basename(base_dir) == "scripts") base_dir <- dirname(base_dir)

data_processed_dir <- file.path(base_dir, "data_processed")
results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(base_dir, "figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

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

if (any(is.na(idx))) {
  cat("\nMostres de counts sense coincidència a metadata:\n")
  print(colnames(counts)[is.na(idx)])
  stop("La metadata no està alineada correctament amb colnames(counts)")
}

metadata <- metadata[idx, , drop = FALSE]
stopifnot(all(metadata$sample_name == colnames(counts)))

y <- factor(metadata$response, levels = c("NonResponder", "Responder"))

cat("\nTaula de resposta:\n")
print(table(y, useNA = "ifany"))

set.seed(123)

top_n_genes <- 100
outer_folds <- 5

folds <- createFolds(y, k = outer_folds, returnTrain = FALSE)

all_predictions <- data.frame(
  sample_name = character(),
  observed = character(),
  predicted = character(),
  prob_responder = numeric(),
  fold = integer(),
  stringsAsFactors = FALSE
)

selected_genes_by_fold <- list()

best_mtry_by_fold <- data.frame(
  fold = integer(),
  best_mtry = numeric(),
  stringsAsFactors = FALSE
)

fold_metrics_all <- data.frame()

get_top_genes_from_train <- function(counts_train, y_train, top_n = 100) {
  dge_train <- DGEList(counts = counts_train, group = y_train)
  keep_genes <- filterByExpr(dge_train, group = y_train)
  dge_train <- dge_train[keep_genes, , keep.lib.sizes = FALSE]
  dge_train <- calcNormFactors(dge_train)
  
  design_train <- model.matrix(~ y_train)
  dge_train <- estimateDisp(dge_train, design_train)
  fit_train <- glmQLFit(dge_train, design_train)
  qlf_train <- glmQLFTest(fit_train, coef = 2)
  
  de_table <- topTags(qlf_train, n = Inf)$table
  top_genes <- rownames(de_table)[1:min(top_n, nrow(de_table))]
  
  list(
    top_genes = top_genes,
    de_table = de_table,
    kept_genes = rownames(dge_train)
  )
}

for (i in seq_along(folds)) {
  cat("\n===========================\n")
  cat("OUTER FOLD", i, "\n")
  cat("===========================\n")
  
  test_idx <- folds[[i]]
  train_idx <- setdiff(seq_len(ncol(counts)), test_idx)
  
  counts_train <- counts[, train_idx, drop = FALSE]
  counts_test  <- counts[, test_idx, drop = FALSE]
  
  y_train <- y[train_idx]
  y_test  <- y[test_idx]
  
  sample_test_names <- colnames(counts_test)
  
  sel <- get_top_genes_from_train(
    counts_train = counts_train,
    y_train = y_train,
    top_n = top_n_genes
  )
  
  top_genes <- sel$top_genes
  selected_genes_by_fold[[paste0("Fold_", i)]] <- top_genes
  
  cat("Gens seleccionats en train:", length(top_genes), "\n")
  
  dge_train_full <- DGEList(counts = counts_train)
  dge_train_full <- calcNormFactors(dge_train_full)
  expr_train <- cpm(dge_train_full, log = TRUE, prior.count = 1)
  
  dge_test_full <- DGEList(counts = counts_test)
  dge_test_full <- calcNormFactors(dge_test_full)
  expr_test <- cpm(dge_test_full, log = TRUE, prior.count = 1)
  
  common_genes <- intersect(top_genes, intersect(rownames(expr_train), rownames(expr_test)))
  
  expr_train_model <- expr_train[common_genes, , drop = FALSE]
  expr_test_model  <- expr_test[common_genes, , drop = FALSE]
  
  X_train <- as.data.frame(t(expr_train_model))
  X_test  <- as.data.frame(t(expr_test_model))
  
  colnames(X_train) <- make.names(colnames(X_train), unique = TRUE)
  colnames(X_test)  <- make.names(colnames(X_test), unique = TRUE)
  
  X_test <- X_test[, colnames(X_train), drop = FALSE]
  
  inner_ctrl <- trainControl(
    method = "cv",
    number = 3
  )
  
  max_p <- ncol(X_train)
  candidate_mtry <- unique(sort(pmax(1, c(2, floor(sqrt(max_p)), floor(max_p / 4)))))
  tunegrid <- expand.grid(mtry = candidate_mtry)
  
  train_df <- X_train
  train_df$response <- y_train
  
  set.seed(100 + i)
  
  rf_tuned <- train(
    response ~ .,
    data = train_df,
    method = "rf",
    metric = "Accuracy",
    trControl = inner_ctrl,
    tuneGrid = tunegrid,
    ntree = 500
  )
  
  best_mtry <- rf_tuned$bestTune$mtry
  
  best_mtry_by_fold <- rbind(
    best_mtry_by_fold,
    data.frame(fold = i, best_mtry = best_mtry)
  )
  
  cat("Best mtry in fold", i, "=", best_mtry, "\n")
  
  set.seed(200 + i)
  
  rf_final <- randomForest(
    x = X_train,
    y = y_train,
    ntree = 500,
    mtry = best_mtry,
    importance = TRUE
  )
  
  pred_class <- predict(rf_final, newdata = X_test, type = "response")
  pred_prob  <- predict(rf_final, newdata = X_test, type = "prob")[, "Responder"]
  
  fold_predictions <- data.frame(
    sample_name = sample_test_names,
    observed = as.character(y_test),
    predicted = as.character(pred_class),
    prob_responder = as.numeric(pred_prob),
    fold = i,
    stringsAsFactors = FALSE
  )
  
  all_predictions <- rbind(all_predictions, fold_predictions)
  
  fold_cm <- confusionMatrix(
    data = factor(pred_class, levels = c("NonResponder", "Responder")),
    reference = factor(y_test, levels = c("NonResponder", "Responder")),
    positive = "Responder"
  )
  
  fold_roc <- roc(
    response = factor(y_test, levels = c("NonResponder", "Responder")),
    predictor = pred_prob,
    levels = c("NonResponder", "Responder"),
    direction = "<",
    quiet = TRUE
  )
  
  fold_metrics <- data.frame(
    fold = i,
    Accuracy = unname(fold_cm$overall["Accuracy"]),
    Kappa = unname(fold_cm$overall["Kappa"]),
    Sensitivity = unname(fold_cm$byClass["Sensitivity"]),
    Specificity = unname(fold_cm$byClass["Specificity"]),
    Balanced_Accuracy = unname(fold_cm$byClass["Balanced Accuracy"]),
    AUC = as.numeric(auc(fold_roc))
  )
  
  fold_metrics_all <- rbind(fold_metrics_all, fold_metrics)
}

all_predictions$observed <- factor(
  all_predictions$observed,
  levels = c("NonResponder", "Responder")
)

all_predictions$predicted <- factor(
  all_predictions$predicted,
  levels = c("NonResponder", "Responder")
)

conf_mat <- confusionMatrix(
  data = all_predictions$predicted,
  reference = all_predictions$observed,
  positive = "Responder"
)

roc_obj <- roc(
  response = all_predictions$observed,
  predictor = all_predictions$prob_responder,
  levels = c("NonResponder", "Responder"),
  direction = "<",
  quiet = TRUE
)

auc_value <- as.numeric(auc(roc_obj))

cat("\n===========================\n")
cat("RESULTATS CV FINAL (100 GENS)\n")
cat("===========================\n")
print(conf_mat)
cat("\nAUC global CV:", round(auc_value, 4), "\n")

gene_frequency <- sort(table(unlist(selected_genes_by_fold)), decreasing = TRUE)

gene_frequency_df <- data.frame(
  gene = names(gene_frequency),
  n_folds_selected = as.integer(gene_frequency),
  stringsAsFactors = FALSE
)

stability_summary <- gene_frequency_df |>
  count(n_folds_selected, name = "n_genes") |>
  arrange(n_folds_selected)

stable_genes_5of5 <- gene_frequency_df |>
  filter(n_folds_selected == 5) |>
  arrange(gene)

stable_genes_4plus <- gene_frequency_df |>
  filter(n_folds_selected >= 4) |>
  arrange(desc(n_folds_selected), gene)

write.csv(
  all_predictions,
  file.path(results_dir, "rf_outer_cv_predictions_final_100_no_leakage.csv"),
  row.names = FALSE
)

write.csv(
  best_mtry_by_fold,
  file.path(results_dir, "rf_best_mtry_by_fold_final_100_no_leakage.csv"),
  row.names = FALSE
)

write.csv(
  fold_metrics_all,
  file.path(results_dir, "rf_fold_metrics_final_100_no_leakage.csv"),
  row.names = FALSE
)

metrics_summary <- data.frame(
  Metric = c(
    "Accuracy",
    "Kappa",
    "Sensitivity",
    "Specificity",
    "Balanced_Accuracy",
    "AUC",
    "Mean_best_mtry",
    "Unique_genes_selected",
    "Genes_selected_4plus_folds",
    "Genes_selected_5of5_folds"
  ),
  Value = c(
    unname(conf_mat$overall["Accuracy"]),
    unname(conf_mat$overall["Kappa"]),
    unname(conf_mat$byClass["Sensitivity"]),
    unname(conf_mat$byClass["Specificity"]),
    unname(conf_mat$byClass["Balanced Accuracy"]),
    auc_value,
    mean(best_mtry_by_fold$best_mtry),
    nrow(gene_frequency_df),
    sum(gene_frequency_df$n_folds_selected >= 4),
    sum(gene_frequency_df$n_folds_selected == 5)
  )
)

write.csv(
  metrics_summary,
  file.path(results_dir, "rf_cv_metrics_summary_final_100_no_leakage.csv"),
  row.names = FALSE
)

write.csv(
  gene_frequency_df,
  file.path(results_dir, "rf_gene_selection_frequency_final_100_no_leakage.csv"),
  row.names = FALSE
)

write.csv(
  stability_summary,
  file.path(results_dir, "gene_stability_summary_final_100.csv"),
  row.names = FALSE
)

write.csv(
  stable_genes_5of5,
  file.path(results_dir, "stable_genes_5of5_final_100.csv"),
  row.names = FALSE
)

write.csv(
  stable_genes_4plus,
  file.path(results_dir, "stable_genes_4plus_final_100.csv"),
  row.names = FALSE
)

png(
  file.path(figures_dir, "ROC_GSE160638_random_forest_final_100_no_leakage.png"),
  width = 1800, height = 1400, res = 220
)
plot(
  roc_obj,
  main = paste("ROC curve - Random Forest final 100-gene signature (AUC =", round(auc_value, 3), ")"),
  col = "darkblue",
  lwd = 3
)
dev.off()

p_stability <- ggplot(stability_summary, aes(x = factor(n_folds_selected), y = n_genes)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Distribució d'estabilitat gènica - firma final de 100 gens",
    x = "Nombre de folds en els que va ser seleccionat",
    y = "Nombre de gens"
  )

ggsave(
  filename = file.path(figures_dir, "gene_stability_distribution_final_100.png"),
  plot = p_stability,
  width = 8,
  height = 6,
  dpi = 300
)

sel_full <- get_top_genes_from_train(
  counts_train = counts,
  y_train = y,
  top_n = top_n_genes
)

dge_full <- DGEList(counts = counts)
dge_full <- calcNormFactors(dge_full)
expr_full <- cpm(dge_full, log = TRUE, prior.count = 1)

common_genes_full <- intersect(sel_full$top_genes, rownames(expr_full))
expr_model_full <- expr_full[common_genes_full, , drop = FALSE]

X_full <- as.data.frame(t(expr_model_full))
colnames(X_full) <- make.names(colnames(X_full), unique = TRUE)

set.seed(999)

rf_full <- randomForest(
  x = X_full,
  y = y,
  ntree = 500,
  importance = TRUE
)

imp_mat <- importance(rf_full)

imp <- data.frame(
  gene = rownames(imp_mat),
  importance = imp_mat[, "MeanDecreaseAccuracy"],
  stringsAsFactors = FALSE
)

imp <- imp[order(-imp$importance), ]

bad_pattern <- "^(RPL|RPS|MIR|SNORD|SNORA|RNU|RN7SK)"
imp_clean <- imp[!grepl(bad_pattern, imp$gene), , drop = FALSE]

write.csv(
  imp,
  file.path(results_dir, "rf_variable_importance_full_dataset_final_100_exploratory.csv"),
  row.names = FALSE
)

write.csv(
  imp_clean,
  file.path(results_dir, "rf_variable_importance_clean_final_100_exploratory.csv"),
  row.names = FALSE
)

top20_clean <- head(imp_clean, 20)

p_imp <- ggplot(
  top20_clean,
  aes(x = reorder(gene, importance), y = importance)
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top 20 gens per importància en RF - firma final 100 gens",
    x = "Gen",
    y = "Mean Decrease Accuracy"
  )

ggsave(
  filename = file.path(figures_dir, "Top20_gene_importance_RF_final_100_exploratory.png"),
  plot = p_imp,
  width = 9,
  height = 7,
  dpi = 300
)

cat("\n===========================\n")
cat("RESUM FINAL - FIRMA 100 GENS\n")
cat("===========================\n")
cat("Nombre de mostres:", ncol(counts), "\n")
cat("Nombre de folds:", outer_folds, "\n")
cat("Top gens per fold:", top_n_genes, "\n")
cat("Accuracy global CV:", round(unname(conf_mat$overall["Accuracy"]), 4), "\n")
cat("Balanced Accuracy global CV:", round(unname(conf_mat$byClass["Balanced Accuracy"]), 4), "\n")
cat("AUC global CV:", round(auc_value, 4), "\n")
cat("Gens únics seleccionats:", nrow(gene_frequency_df), "\n")
cat("Gens seleccionats en >=4 folds:", sum(gene_frequency_df$n_folds_selected >= 4), "\n")
cat("Gens seleccionats en 5/5 folds:", sum(gene_frequency_df$n_folds_selected == 5), "\n")
cat("===========================\n")

