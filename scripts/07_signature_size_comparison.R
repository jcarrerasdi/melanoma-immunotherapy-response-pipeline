################################
# 07_signature_size_comparison.R
################################

suppressPackageStartupMessages({
  library(edgeR)
  library(randomForest)
  library(pROC)
  library(caret)
  library(dplyr)
  library(ggplot2)
})

# Project root: run scripts either from the project root or from the scripts/ folder.
base_dir <- normalizePath(getwd(), mustWork = TRUE)
if (basename(base_dir) == "scripts") base_dir <- dirname(base_dir)

data_processed_dir <- file.path(base_dir, "data_processed")
results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(base_dir, "figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ===========================
# 1. Carregar dades
# ===========================

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

# ===========================
# 2. Paràmetres generals
# ===========================

set.seed(123)

signature_sizes <- c(50, 100, 250, 500)
outer_folds <- 5

# mateixos folds per a totes les comparacions
folds <- createFolds(y, k = outer_folds, returnTrain = FALSE)

# ===========================
# 3. Funció auxiliar:
#    selecció de gens NOMÉS en train
# ===========================

get_top_genes_from_train <- function(counts_train, y_train, top_n = 500) {
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

# ===========================
# 4. Objectes globals
# ===========================

all_metrics <- data.frame()
all_fold_metrics <- data.frame()
all_stability_summary <- data.frame()
all_top_stable_genes <- data.frame()

# ===========================
# 5. Bucle principal per mida de signatura
# ===========================

for (top_n_genes in signature_sizes) {
  
  cat("\n============================================\n")
  cat("INICIANT ANÀLISI PER A top_n_genes =", top_n_genes, "\n")
  cat("============================================\n")
  
  all_predictions <- data.frame(
    sample_name = character(),
    observed = character(),
    predicted = character(),
    prob_responder = numeric(),
    fold = integer(),
    signature_size = integer(),
    stringsAsFactors = FALSE
  )
  
  selected_genes_by_fold <- list()
  best_mtry_by_fold <- data.frame(
    fold = integer(),
    best_mtry = numeric(),
    signature_size = integer(),
    stringsAsFactors = FALSE
  )
  
  # ---------------------------
  # 5.1 CV externa
  # ---------------------------
  for (i in seq_along(folds)) {
    cat("\n---------------------------\n")
    cat("Signature size:", top_n_genes, "| OUTER FOLD", i, "\n")
    cat("---------------------------\n")
    
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_len(ncol(counts)), test_idx)
    
    counts_train <- counts[, train_idx, drop = FALSE]
    counts_test  <- counts[, test_idx, drop = FALSE]
    
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    sample_test_names <- colnames(counts_test)
    
    # Selecció NOMÉS amb train
    sel <- get_top_genes_from_train(
      counts_train = counts_train,
      y_train = y_train,
      top_n = top_n_genes
    )
    
    top_genes <- sel$top_genes
    selected_genes_by_fold[[paste0("Fold_", i)]] <- top_genes
    
    cat("Gens seleccionats en train:", length(top_genes), "\n")
    
    # Transformació train/test
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
    
    # Tuning NOMÉS amb train
    inner_ctrl <- trainControl(
      method = "cv",
      number = 3
    )
    
    max_p <- ncol(X_train)
    candidate_mtry <- unique(sort(pmax(1, c(2, floor(sqrt(max_p)), floor(max_p / 4)))))
    tunegrid <- expand.grid(mtry = candidate_mtry)
    
    train_df <- X_train
    train_df$response <- y_train
    
    set.seed(1000 + top_n_genes + i)
    
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
      data.frame(
        fold = i,
        best_mtry = best_mtry,
        signature_size = top_n_genes
      )
    )
    
    cat("Best mtry in fold", i, "=", best_mtry, "\n")
    
    # Model final del fold
    set.seed(2000 + top_n_genes + i)
    
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
      signature_size = top_n_genes,
      stringsAsFactors = FALSE
    )
    
    all_predictions <- rbind(all_predictions, fold_predictions)
    
    # Mètriques per fold
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
      signature_size = top_n_genes,
      fold = i,
      Accuracy = unname(fold_cm$overall["Accuracy"]),
      Kappa = unname(fold_cm$overall["Kappa"]),
      Sensitivity = unname(fold_cm$byClass["Sensitivity"]),
      Specificity = unname(fold_cm$byClass["Specificity"]),
      Balanced_Accuracy = unname(fold_cm$byClass["Balanced Accuracy"]),
      AUC = as.numeric(auc(fold_roc))
    )
    
    all_fold_metrics <- rbind(all_fold_metrics, fold_metrics)
  }
  
  # ---------------------------
  # 5.2 Mètriques globals
  # ---------------------------
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
  
  # ---------------------------
  # 5.3 Estabilitat gènica
  # ---------------------------
  gene_frequency <- sort(table(unlist(selected_genes_by_fold)), decreasing = TRUE)
  
  gene_frequency_df <- data.frame(
    gene = names(gene_frequency),
    n_folds_selected = as.integer(gene_frequency),
    signature_size = top_n_genes,
    stringsAsFactors = FALSE
  )
  
  stability_summary <- gene_frequency_df |>
    count(signature_size, n_folds_selected, name = "n_genes") |>
    arrange(signature_size, n_folds_selected)
  
  n_unique_genes <- nrow(gene_frequency_df)
  n_genes_5of5 <- sum(gene_frequency_df$n_folds_selected == outer_folds)
  n_genes_4plus <- sum(gene_frequency_df$n_folds_selected >= 4)
  
  # top gens més estables
  top_stable_genes <- gene_frequency_df |>
    arrange(desc(n_folds_selected), gene) |>
    head(30)
  
  all_top_stable_genes <- rbind(all_top_stable_genes, top_stable_genes)
  all_stability_summary <- rbind(all_stability_summary, stability_summary)
  
  # ---------------------------
  # 5.4 Desar fitxers per mida
  # ---------------------------
  write.csv(
    all_predictions,
    file.path(
      results_dir,
      paste0("rf_outer_cv_predictions_signature_", top_n_genes, "_no_leakage.csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    best_mtry_by_fold,
    file.path(
      results_dir,
      paste0("rf_best_mtry_by_fold_signature_", top_n_genes, "_no_leakage.csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    gene_frequency_df,
    file.path(
      results_dir,
      paste0("rf_gene_selection_frequency_signature_", top_n_genes, "_no_leakage.csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    stability_summary,
    file.path(
      results_dir,
      paste0("gene_stability_summary_signature_", top_n_genes, ".csv")
    ),
    row.names = FALSE
  )
  
  # ---------------------------
  # 5.5 Resum global de mètriques
  # ---------------------------
  metrics_summary <- data.frame(
    signature_size = top_n_genes,
    Accuracy = unname(conf_mat$overall["Accuracy"]),
    Kappa = unname(conf_mat$overall["Kappa"]),
    Sensitivity = unname(conf_mat$byClass["Sensitivity"]),
    Specificity = unname(conf_mat$byClass["Specificity"]),
    Balanced_Accuracy = unname(conf_mat$byClass["Balanced Accuracy"]),
    AUC = auc_value,
    n_unique_genes_selected = n_unique_genes,
    n_genes_selected_4plus_folds = n_genes_4plus,
    n_genes_selected_5of5_folds = n_genes_5of5,
    mean_best_mtry = mean(best_mtry_by_fold$best_mtry)
  )
  
  all_metrics <- rbind(all_metrics, metrics_summary)
  
  write.csv(
    metrics_summary,
    file.path(
      results_dir,
      paste0("rf_cv_metrics_summary_signature_", top_n_genes, "_no_leakage.csv")
    ),
    row.names = FALSE
  )
  
  cat("\n============================================\n")
  cat("RESULTATS PER A SIGNATURE SIZE =", top_n_genes, "\n")
  cat("============================================\n")
  print(metrics_summary)
}

# ===========================
# 6. Desar resultats globals comparatius
# ===========================

write.csv(
  all_metrics,
  file.path(results_dir, "signature_size_comparison_metrics.csv"),
  row.names = FALSE
)

write.csv(
  all_fold_metrics,
  file.path(results_dir, "signature_size_comparison_fold_metrics.csv"),
  row.names = FALSE
)

write.csv(
  all_stability_summary,
  file.path(results_dir, "signature_size_comparison_stability_summary.csv"),
  row.names = FALSE
)

write.csv(
  all_top_stable_genes,
  file.path(results_dir, "signature_size_comparison_top_stable_genes.csv"),
  row.names = FALSE
)

# ===========================
# 7. Figures comparatives
# ===========================

# 7.1 AUC per mida de signatura
p_auc <- ggplot(all_metrics, aes(x = factor(signature_size), y = AUC, group = 1)) +
  geom_point(size = 3) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Comparison of AUC by signatue size",
    x = "Number of genes selected per fold",
    y = "Global CV AUC"
  )

ggsave(
  filename = file.path(figures_dir, "signature_size_comparison_auc.png"),
  plot = p_auc,
  width = 8,
  height = 6,
  dpi = 300
)

# 7.2 Accuracy per mida de signatura
p_acc <- ggplot(all_metrics, aes(x = factor(signature_size), y = Accuracy, group = 1)) +
  geom_point(size = 3) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Comparació d'accuracy segons la mida de la signatura",
    x = "Nombre de gens seleccionats per fold",
    y = "Accuracy global CV"
  )

ggsave(
  filename = file.path(figures_dir, "signature_size_comparison_accuracy.png"),
  plot = p_acc,
  width = 8,
  height = 6,
  dpi = 300
)

# 7.3 Balanced Accuracy per mida de signatura
p_bacc <- ggplot(all_metrics, aes(x = factor(signature_size), y = Balanced_Accuracy, group = 1)) +
  geom_point(size = 3) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Comparació de balanced accuracy segons la mida de la signatura",
    x = "Nombre de gens seleccionats per fold",
    y = "Balanced accuracy global CV"
  )

ggsave(
  filename = file.path(figures_dir, "signature_size_comparison_balanced_accuracy.png"),
  plot = p_bacc,
  width = 8,
  height = 6,
  dpi = 300
)

# 7.4 Nombre de gens únics seleccionats
p_unique <- ggplot(all_metrics, aes(x = factor(signature_size), y = n_unique_genes_selected, group = 1)) +
  geom_point(size = 3) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Nombre de gens únics seleccionats segons la mida de la signatura",
    x = "Nombre de gens seleccionats per fold",
    y = "Gens únics seleccionats"
  )

ggsave(
  filename = file.path(figures_dir, "signature_size_comparison_unique_genes.png"),
  plot = p_unique,
  width = 8,
  height = 6,
  dpi = 300
)

# 7.5 Nombre de gens seleccionats en 5/5 folds
p_stable5 <- ggplot(all_metrics, aes(x = factor(signature_size), y = n_genes_selected_5of5_folds, group = 1)) +
  geom_point(size = 3) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Gens seleccionats en 5/5 folds segons la mida de la signatura",
    x = "Nombre de gens seleccionats per fold",
    y = "Nombre de gens en 5/5 folds"
  )

ggsave(
  filename = file.path(figures_dir, "signature_size_comparison_genes_5of5.png"),
  plot = p_stable5,
  width = 8,
  height = 6,
  dpi = 300
)

# 7.6 Boxplot AUC per fold
p_auc_fold <- ggplot(all_fold_metrics, aes(x = factor(signature_size), y = AUC)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Distribució de l'AUC per fold segons la mida de la signatura",
    x = "Nombre de gens seleccionats per fold",
    y = "AUC por fold"
  )

ggsave(
  filename = file.path(figures_dir, "signature_size_comparison_auc_by_fold.png"),
  plot = p_auc_fold,
  width = 8,
  height = 6,
  dpi = 300
)

# 7.7 Boxplot Accuracy per fold
p_acc_fold <- ggplot(all_fold_metrics, aes(x = factor(signature_size), y = Accuracy)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Distribució de l'accuracy per fold segons la mida de la signatura",
    x = "Nombre de gens seleccionats per fold",
    y = "Accuracy por fold"
  )

ggsave(
  filename = file.path(figures_dir, "signature_size_comparison_accuracy_by_fold.png"),
  plot = p_acc_fold,
  width = 8,
  height = 6,
  dpi = 300
)

cat("\n============================================\n")
cat("ANÀLISI DE COMPARACIÓ DE MIDA DE SIGNATURA COMPLETADA\n")
cat("============================================\n")

cat("\nResum comparatiu final:\n")
print(all_metrics[order(-all_metrics$AUC), ])

