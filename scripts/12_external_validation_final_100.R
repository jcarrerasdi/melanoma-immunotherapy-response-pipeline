# ============================================================
# 12_external_validation_final_100.R
# Cleaned/renamed version of: 10_external_validation_final_100.R
# Purpose: reproducible TFM pipeline while preserving the output filenames used in the memoria.
# ============================================================

################################
# 10_external_validation_final_100.R
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

set.seed(999)

top_n_genes <- 100
external_datasets <- c("GSE91061", "GSE78220")

# =========================
# Funciones auxiliares
# =========================

clean_ids <- function(x) {
  trimws(as.character(x))
}

normalize_gene_ids <- function(x) {
  x <- clean_ids(x)
  x <- sub("\\.\\d+$", "", x)
  x <- gsub("\\s+", "", x)
  toupper(x)
}

normalize_sample_key <- function(x) {
  x <- clean_ids(x)
  x <- toupper(x)
  gsub("[^A-Z0-9]", "", x)
}

normalize_sample_key_by_dataset <- function(x, dataset_id, is_expression = FALSE) {
  x <- clean_ids(x)
  x <- toupper(x)
  
  if (dataset_id == "GSE78220") {
    x <- gsub("\\.BASELINE$", "", x)
    x <- gsub("_BASELINE$", "", x)
    x <- gsub("-BASELINE$", "", x)
    x <- gsub("BASELINE$", "", x)
  }
  
  gsub("[^A-Z0-9]", "", x)
}

is_raw_counts_matrix <- function(mat) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  
  vals <- as.numeric(mat)
  vals <- vals[is.finite(vals)]
  
  if (length(vals) == 0) return(FALSE)
  
  prop_integer_like <- mean(abs(vals - round(vals)) < 1e-6)
  max_val <- max(vals, na.rm = TRUE)
  median_val <- median(vals, na.rm = TRUE)
  
  prop_integer_like > 0.9 && max_val > 50 && median_val >= 1
}

normalize_external_expression <- function(expr_mat) {
  expr_mat <- as.matrix(expr_mat)
  mode(expr_mat) <- "numeric"
  
  if (is_raw_counts_matrix(expr_mat)) {
    dge <- DGEList(counts = expr_mat)
    dge <- calcNormFactors(dge)
    expr_norm <- cpm(dge, log = TRUE, prior.count = 1)
    method_used <- "TMM_logCPM"
  } else {
    if (max(expr_mat, na.rm = TRUE) > 50) {
      expr_norm <- log2(expr_mat + 1)
      method_used <- "log2(x+1)"
    } else {
      expr_norm <- expr_mat
      method_used <- "already_log_like"
    }
  }
  
  list(expr = expr_norm, method = method_used)
}

compute_binary_metrics <- function(observed, predicted, prob_responder, dataset_name) {
  observed <- factor(observed, levels = c("NonResponder", "Responder"))
  predicted <- factor(predicted, levels = c("NonResponder", "Responder"))
  
  cm <- confusionMatrix(
    data = predicted,
    reference = observed,
    positive = "Responder"
  )
  
  roc_obj <- roc(
    response = observed,
    predictor = prob_responder,
    levels = c("NonResponder", "Responder"),
    direction = "<",
    quiet = TRUE
  )
  
  data.frame(
    dataset = dataset_name,
    n_samples = length(observed),
    Accuracy = unname(cm$overall["Accuracy"]),
    Kappa = unname(cm$overall["Kappa"]),
    Sensitivity = unname(cm$byClass["Sensitivity"]),
    Specificity = unname(cm$byClass["Specificity"]),
    Balanced_Accuracy = unname(cm$byClass["Balanced Accuracy"]),
    AUC = as.numeric(auc(roc_obj)),
    stringsAsFactors = FALSE
  )
}

compute_metrics_at_threshold <- function(observed, prob_responder, threshold, dataset_name, threshold_type) {
  observed <- factor(observed, levels = c("NonResponder", "Responder"))
  
  predicted <- ifelse(prob_responder >= threshold, "Responder", "NonResponder")
  predicted <- factor(predicted, levels = c("NonResponder", "Responder"))
  
  cm <- confusionMatrix(
    data = predicted,
    reference = observed,
    positive = "Responder"
  )
  
  roc_obj <- roc(
    response = observed,
    predictor = prob_responder,
    levels = c("NonResponder", "Responder"),
    direction = "<",
    quiet = TRUE
  )
  
  data.frame(
    dataset = dataset_name,
    threshold_type = threshold_type,
    threshold_value = as.numeric(threshold),
    n_samples = length(observed),
    Accuracy = unname(cm$overall["Accuracy"]),
    Kappa = unname(cm$overall["Kappa"]),
    Sensitivity = unname(cm$byClass["Sensitivity"]),
    Specificity = unname(cm$byClass["Specificity"]),
    Balanced_Accuracy = unname(cm$byClass["Balanced Accuracy"]),
    AUC = as.numeric(auc(roc_obj)),
    stringsAsFactors = FALSE
  )
}

select_threshold_youden <- function(observed, prob_responder) {
  observed <- factor(observed, levels = c("NonResponder", "Responder"))
  
  roc_obj <- roc(
    response = observed,
    predictor = prob_responder,
    levels = c("NonResponder", "Responder"),
    direction = "<",
    quiet = TRUE
  )
  
  best_coords <- coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  data.frame(
    threshold_type = "youden",
    threshold_value = as.numeric(best_coords["threshold"]),
    sensitivity_cv = as.numeric(best_coords["sensitivity"]),
    specificity_cv = as.numeric(best_coords["specificity"]),
    stringsAsFactors = FALSE
  )
}

select_threshold_sensitivity <- function(observed, prob_responder, min_specificity = 0.40) {
  observed <- factor(observed, levels = c("NonResponder", "Responder"))
  
  roc_obj <- roc(
    response = observed,
    predictor = prob_responder,
    levels = c("NonResponder", "Responder"),
    direction = "<",
    quiet = TRUE
  )
  
  coords_all <- coords(
    roc_obj,
    x = "all",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  coords_df <- as.data.frame(coords_all)
  coords_df <- coords_df[is.finite(coords_df$threshold), , drop = FALSE]
  
  candidate <- coords_df %>%
    filter(specificity >= min_specificity) %>%
    arrange(desc(sensitivity), desc(specificity), threshold)
  
  if (nrow(candidate) == 0) {
    candidate <- coords_df %>%
      arrange(desc(sensitivity), desc(specificity), threshold)
  }
  
  candidate <- candidate[1, , drop = FALSE]
  
  data.frame(
    threshold_type = "sensitivity_oriented",
    threshold_value = as.numeric(candidate$threshold),
    sensitivity_cv = as.numeric(candidate$sensitivity),
    specificity_cv = as.numeric(candidate$specificity),
    stringsAsFactors = FALSE
  )
}

align_external_expression_and_metadata <- function(expr_ext, meta_ext, dataset_id) {
  expr_colnames_original <- colnames(expr_ext)
  
  if (dataset_id == "GSE78220") {
    expr_names_clean <- toupper(clean_ids(expr_colnames_original))
    keep_expr <- !grepl("ONTX", expr_names_clean)
    expr_ext <- expr_ext[, keep_expr, drop = FALSE]
    expr_colnames_original <- colnames(expr_ext)
  }
  
  expr_keys <- normalize_sample_key_by_dataset(
    expr_colnames_original,
    dataset_id = dataset_id,
    is_expression = TRUE
  )
  
  meta_ext$sample_name <- clean_ids(meta_ext$sample_name)
  meta_ext$sample_key <- normalize_sample_key_by_dataset(
    meta_ext$sample_name,
    dataset_id = dataset_id,
    is_expression = FALSE
  )
  
  meta_ext <- meta_ext %>%
    filter(!is.na(sample_name), sample_name != "") %>%
    filter(!is.na(sample_key), sample_key != "") %>%
    distinct(sample_key, .keep_all = TRUE)
  
  common_keys <- intersect(expr_keys, meta_ext$sample_key)
  
  cat("\n[", dataset_id, "] muestras en expresión:", ncol(expr_ext), "\n", sep = "")
  cat("[", dataset_id, "] muestras en metadata:", nrow(meta_ext), "\n", sep = "")
  cat("[", dataset_id, "] muestras comunes tras normalización:", length(common_keys), "\n", sep = "")
  
  if (length(common_keys) < 2) {
    cat("\n[", dataset_id, "] Primeras columnas de expresión:\n", sep = "")
    print(head(expr_colnames_original, 20))
    cat("\n[", dataset_id, "] Primeros sample_name de metadata:\n", sep = "")
    print(head(meta_ext$sample_name, 20))
    stop(paste0("[", dataset_id, "] Menos de 2 muestras comunes entre expresión y metadata"))
  }
  
  expr_idx <- match(common_keys, expr_keys)
  meta_idx <- match(common_keys, meta_ext$sample_key)
  
  expr_ext <- expr_ext[, expr_idx, drop = FALSE]
  meta_ext <- meta_ext[meta_idx, , drop = FALSE]
  
  meta_ext$sample_name <- colnames(expr_ext)
  
  stopifnot(ncol(expr_ext) == nrow(meta_ext))
  stopifnot(all(meta_ext$sample_name == colnames(expr_ext)))
  
  list(expr = expr_ext, meta = meta_ext)
}

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

# =========================
# 1. Cargar training común
# =========================

counts <- readRDS(
  file.path(data_processed_dir, "GSE160638_common_genes_aligned.rds")
)

metadata <- read.csv(
  file.path(data_processed_dir, "metadata_GSE160638_aligned.csv"),
  stringsAsFactors = FALSE
)

counts <- as.matrix(counts)
mode(counts) <- "numeric"

rownames(counts) <- normalize_gene_ids(rownames(counts))
colnames(counts) <- clean_ids(colnames(counts))

metadata[] <- lapply(metadata, function(x) {
  if (is.character(x)) clean_ids(x) else x
})

if (!"sample_name" %in% colnames(metadata)) {
  stop("La columna 'sample_name' no existe en metadata_GSE160638_aligned.csv")
}

if (!"response" %in% colnames(metadata)) {
  stop("La columna 'response' no existe en metadata_GSE160638_aligned.csv")
}

idx <- match(colnames(counts), metadata$sample_name)

if (any(is.na(idx))) {
  cat("\nMuestras de counts sin match en metadata:\n")
  print(colnames(counts)[is.na(idx)])
  stop("La metadata no está alineada correctamente con colnames(counts)")
}

metadata <- metadata[idx, , drop = FALSE]
stopifnot(all(metadata$sample_name == colnames(counts)))

y <- factor(metadata$response, levels = c("NonResponder", "Responder"))

cat("\n===========================\n")
cat("TRAINING SET\n")
cat("===========================\n")
cat("Dimensiones counts:", dim(counts), "\n")
cat("Distribución de respuesta:\n")
print(table(y, useNA = "ifany"))

# =========================
# 2. Selección final de genes
# =========================

sel_full <- get_top_genes_from_train(
  counts_train = counts,
  y_train = y,
  top_n = top_n_genes
)

selected_genes <- sel_full$top_genes
selected_genes <- intersect(selected_genes, rownames(counts))
selected_genes <- unique(selected_genes)

cat("\nGenes seleccionados inicialmente:", length(sel_full$top_genes), "\n")
cat("Genes seleccionados tras restricción a espacio común:", length(selected_genes), "\n")

if (length(selected_genes) < 10) {
  stop("Quedan muy pocos genes tras restringir al espacio común. Revisa la preparación previa.")
}

write.csv(
  data.frame(gene = selected_genes),
  file.path(results_dir, "final_selected_genes_external_validation_final_100.csv"),
  row.names = FALSE
)

# =========================
# 3. Normalización training
# =========================

dge_full <- DGEList(counts = counts)
dge_full <- calcNormFactors(dge_full)
expr_train <- cpm(dge_full, log = TRUE, prior.count = 1)

expr_train <- expr_train[selected_genes, , drop = FALSE]

X_train <- as.data.frame(t(expr_train))
colnames(X_train) <- make.names(colnames(X_train), unique = TRUE)

# =========================
# 4. Tuning final en training
# =========================

inner_ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  savePredictions = "final"
)

max_p <- ncol(X_train)
candidate_mtry <- unique(sort(pmax(1, c(2, floor(sqrt(max_p)), floor(max_p / 4)))))
tunegrid <- expand.grid(mtry = candidate_mtry)

train_df <- X_train
train_df$response <- y

cat("\n===========================\n")
cat("TUNING FINAL MODEL\n")
cat("===========================\n")
cat("Candidate mtry:", paste(candidate_mtry, collapse = ", "), "\n")

set.seed(2026)

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

cat("Best mtry final =", best_mtry, "\n")

write.csv(
  rf_tuned$results,
  file.path(results_dir, "rf_final_model_tuning_results_external_validation_final_100.csv"),
  row.names = FALSE
)

# =========================
# 5. Modelo final
# =========================

set.seed(3030)

rf_final <- randomForest(
  x = X_train,
  y = y,
  ntree = 500,
  mtry = best_mtry,
  importance = TRUE
)

imp_mat <- importance(rf_final)
imp <- data.frame(
  gene = rownames(imp_mat),
  importance = imp_mat[, "MeanDecreaseAccuracy"],
  stringsAsFactors = FALSE
)
imp <- imp[order(-imp$importance), ]

write.csv(
  imp,
  file.path(results_dir, "rf_variable_importance_external_validation_final_100.csv"),
  row.names = FALSE
)

# =========================
# 5B. Threshold analysis
#     Umbrales definidos SOLO
#     con training / CV
# =========================

cv_pred <- rf_tuned$pred

if (!is.null(cv_pred) && nrow(cv_pred) > 0) {
  cv_pred_best <- cv_pred %>%
    filter(mtry == best_mtry)
  
  if (!"Responder" %in% colnames(cv_pred_best)) {
    stop("No se encontró la columna de probabilidad 'Responder' en rf_tuned$pred")
  }
  
  observed_cv <- factor(cv_pred_best$obs, levels = c("NonResponder", "Responder"))
  prob_cv <- as.numeric(cv_pred_best$Responder)
  
} else {
  warning("rf_tuned$pred no está disponible. Se usarán probabilidades aparentes del training para definir thresholds.")
  observed_cv <- y
  prob_cv <- as.numeric(predict(rf_final, newdata = X_train, type = "prob")[, "Responder"])
}

threshold_fixed <- data.frame(
  threshold_type = "fixed_0.5",
  threshold_value = 0.5,
  sensitivity_cv = NA_real_,
  specificity_cv = NA_real_,
  stringsAsFactors = FALSE
)

threshold_youden <- select_threshold_youden(
  observed = observed_cv,
  prob_responder = prob_cv
)

threshold_sens <- select_threshold_sensitivity(
  observed = observed_cv,
  prob_responder = prob_cv,
  min_specificity = 0.40
)

threshold_table <- bind_rows(
  threshold_fixed,
  threshold_youden,
  threshold_sens
)

write.csv(
  threshold_table,
  file.path(results_dir, "thresholds_external_validation_final_100.csv"),
  row.names = FALSE
)

cat("\n===========================\n")
cat("THRESHOLD ANALYSIS\n")
cat("===========================\n")
print(threshold_table)

# =========================
# 6. Validación externa
# =========================

all_external_predictions <- data.frame()
all_external_metrics <- data.frame()
normalization_audit <- data.frame()

for (ds in external_datasets) {
  cat("\n===========================\n")
  cat("VALIDANDO EN", ds, "\n")
  cat("===========================\n")
  
  expr_ext <- readRDS(
    file.path(data_processed_dir, paste0(ds, "_common_genes_aligned.rds"))
  )
  
  meta_ext <- read.csv(
    file.path(data_processed_dir, paste0("metadata_", ds, "_harmonized_kept.csv")),
    stringsAsFactors = FALSE
  )
  
  expr_ext <- as.matrix(expr_ext)
  mode(expr_ext) <- "numeric"
  
  rownames(expr_ext) <- normalize_gene_ids(rownames(expr_ext))
  colnames(expr_ext) <- clean_ids(colnames(expr_ext))
  
  meta_ext[] <- lapply(meta_ext, function(x) {
    if (is.character(x)) clean_ids(x) else x
  })
  
  if (!"sample_name" %in% colnames(meta_ext)) {
    stop(paste0("[", ds, "] falta la columna 'sample_name' en metadata alineada"))
  }
  if (!"response" %in% colnames(meta_ext)) {
    stop(paste0("[", ds, "] falta la columna 'response' en metadata alineada"))
  }
  
  aligned_ext <- align_external_expression_and_metadata(
    expr_ext = expr_ext,
    meta_ext = meta_ext,
    dataset_id = ds
  )
  
  expr_ext <- aligned_ext$expr
  meta_ext <- aligned_ext$meta
  
  y_ext <- factor(meta_ext$response, levels = c("NonResponder", "Responder"))
  
  norm_res <- normalize_external_expression(expr_ext)
  expr_ext_norm <- norm_res$expr
  method_used <- norm_res$method
  
  cat("Método de normalización para", ds, ":", method_used, "\n")
  cat("Dimensiones expr_ext después de alinear muestras:", dim(expr_ext_norm), "\n")
  
  normalization_audit <- rbind(
    normalization_audit,
    data.frame(
      dataset = ds,
      normalization_method = method_used,
      n_samples = ncol(expr_ext_norm),
      n_genes = nrow(expr_ext_norm),
      stringsAsFactors = FALSE
    )
  )
  
  genes_available <- intersect(selected_genes, rownames(expr_ext_norm))
  
  cat("Genes de la firma presentes en", ds, ":", length(genes_available), "de", length(selected_genes), "\n")
  
  if (length(genes_available) < 10) {
    warning(paste0("[", ds, "] Hay menos de 10 genes disponibles de la firma final."))
  }
  
  X_ext_mat <- matrix(
    0,
    nrow = ncol(expr_ext_norm),
    ncol = length(selected_genes),
    dimnames = list(colnames(expr_ext_norm), selected_genes)
  )
  
  X_ext_mat[, genes_available] <- t(expr_ext_norm[genes_available, , drop = FALSE])
  
  X_ext <- as.data.frame(X_ext_mat)
  colnames(X_ext) <- make.names(colnames(X_ext), unique = TRUE)
  
  X_ext <- X_ext[, colnames(X_train), drop = FALSE]
  
  pred_prob <- predict(rf_final, newdata = X_ext, type = "prob")[, "Responder"]
  
  for (j in seq_len(nrow(threshold_table))) {
    th_type <- threshold_table$threshold_type[j]
    th_value <- threshold_table$threshold_value[j]
    
    pred_class <- ifelse(pred_prob >= th_value, "Responder", "NonResponder")
    
    pred_df <- data.frame(
      dataset = ds,
      threshold_type = th_type,
      threshold_value = th_value,
      sample_name = rownames(X_ext),
      observed = as.character(y_ext),
      predicted = as.character(pred_class),
      prob_responder = as.numeric(pred_prob),
      genes_available_from_signature = length(genes_available),
      signature_size_total = length(selected_genes),
      stringsAsFactors = FALSE
    )
    
    all_external_predictions <- rbind(all_external_predictions, pred_df)
    
    metrics_df <- compute_metrics_at_threshold(
      observed = y_ext,
      prob_responder = pred_prob,
      threshold = th_value,
      dataset_name = ds,
      threshold_type = th_type
    )
    
    metrics_df$genes_available_from_signature <- length(genes_available)
    metrics_df$signature_size_total <- length(selected_genes)
    
    all_external_metrics <- rbind(all_external_metrics, metrics_df)
    
    cat("\n---", ds, "|", th_type, "| threshold =", round(th_value, 4), "---\n")
    print(metrics_df)
  }
  
  roc_obj_ds <- roc(
    response = factor(y_ext, levels = c("NonResponder", "Responder")),
    predictor = pred_prob,
    levels = c("NonResponder", "Responder"),
    direction = "<",
    quiet = TRUE
  )
  
  auc_ds <- as.numeric(auc(roc_obj_ds))
  
  png(
    file.path(figures_dir, paste0("ROC_", ds, "_external_validation_final_100.png")),
    width = 1800, height = 1400, res = 220
  )
  plot(
    roc_obj_ds,
    main = paste("ROC curve -", ds, "(AUC =", round(auc_ds, 3), ")"),
    col = "darkblue",
    lwd = 3
  )
  dev.off()
}

# =========================
# 7. Métricas globales combinadas
# =========================

combined_metrics <- data.frame()

for (j in seq_len(nrow(threshold_table))) {
  th_type <- threshold_table$threshold_type[j]
  th_value <- threshold_table$threshold_value[j]
  
  pred_subset <- all_external_predictions %>%
    filter(threshold_type == th_type)
  
  metrics_comb <- compute_metrics_at_threshold(
    observed = factor(pred_subset$observed, levels = c("NonResponder", "Responder")),
    prob_responder = pred_subset$prob_responder,
    threshold = th_value,
    dataset_name = "Combined_external",
    threshold_type = th_type
  )
  
  combined_metrics <- rbind(combined_metrics, metrics_comb)
}

cat("\n===========================\n")
cat("RESULTADOS EXTERNOS COMBINADOS\n")
cat("===========================\n")
print(combined_metrics)

roc_combined <- roc(
  response = factor(
    all_external_predictions %>%
      filter(threshold_type == "fixed_0.5") %>%
      pull(observed),
    levels = c("NonResponder", "Responder")
  ),
  predictor = all_external_predictions %>%
    filter(threshold_type == "fixed_0.5") %>%
    pull(prob_responder),
  levels = c("NonResponder", "Responder"),
  direction = "<",
  quiet = TRUE
)

auc_combined <- as.numeric(auc(roc_combined))

png(
  file.path(figures_dir, "ROC_external_validation_combined_final_100.png"),
  width = 1800, height = 1400, res = 220
)
plot(
  roc_combined,
  main = paste("ROC curve - Combined external validation (AUC =", round(auc_combined, 3), ")"),
  col = "darkblue",
  lwd = 3
)
dev.off()

# =========================
# 8. Guardado de resultados
# =========================

write.csv(
  all_external_predictions,
  file.path(results_dir, "external_validation_predictions_with_thresholds_final_100.csv"),
  row.names = FALSE
)

write.csv(
  all_external_metrics,
  file.path(results_dir, "external_validation_metrics_by_dataset_with_thresholds_final_100.csv"),
  row.names = FALSE
)

write.csv(
  combined_metrics,
  file.path(results_dir, "external_validation_metrics_combined_with_thresholds_final_100.csv"),
  row.names = FALSE
)

write.csv(
  normalization_audit,
  file.path(results_dir, "external_validation_normalization_audit_final_100.csv"),
  row.names = FALSE
)

# Compatibilidad con tus nombres previos:
metrics_fixed_05 <- all_external_metrics %>%
  filter(threshold_type == "fixed_0.5")

combined_fixed_05 <- combined_metrics %>%
  filter(threshold_type == "fixed_0.5")

pred_fixed_05 <- all_external_predictions %>%
  filter(threshold_type == "fixed_0.5")

write.csv(
  pred_fixed_05,
  file.path(results_dir, "external_validation_predictions_final_100.csv"),
  row.names = FALSE
)

write.csv(
  metrics_fixed_05,
  file.path(results_dir, "external_validation_metrics_by_dataset_final_100.csv"),
  row.names = FALSE
)

write.csv(
  combined_fixed_05,
  file.path(results_dir, "external_validation_metrics_combined_final_100.csv"),
  row.names = FALSE
)

# =========================
# 9. Figuras resumen
# =========================

plot_metrics <- bind_rows(
  all_external_metrics %>%
    mutate(dataset_plot = paste(dataset, threshold_type, sep = " | ")) %>%
    dplyr::select(dataset_plot, AUC),
  combined_metrics %>%
    mutate(dataset_plot = paste(dataset, threshold_type, sep = " | ")) %>%
    dplyr::select(dataset_plot, AUC)
)

p_auc <- ggplot(plot_metrics, aes(x = dataset_plot, y = AUC)) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "AUC en validación externa según threshold",
    x = "Dataset | threshold",
    y = "AUC"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(figures_dir, "external_validation_auc_summary_with_thresholds_final_100.png"),
  plot = p_auc,
  width = 10,
  height = 6,
  dpi = 300
)

plot_bacc <- bind_rows(
  all_external_metrics %>%
    mutate(dataset_plot = paste(dataset, threshold_type, sep = " | ")) %>%
    dplyr::select(dataset_plot, Balanced_Accuracy),
  combined_metrics %>%
    mutate(dataset_plot = paste(dataset, threshold_type, sep = " | ")) %>%
    dplyr::select(dataset_plot, Balanced_Accuracy)
)

p_bacc <- ggplot(plot_bacc, aes(x = dataset_plot, y = Balanced_Accuracy)) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "Balanced accuracy en validación externa según threshold",
    x = "Dataset | threshold",
    y = "Balanced accuracy"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(figures_dir, "external_validation_balanced_accuracy_summary_with_thresholds_final_100.png"),
  plot = p_bacc,
  width = 10,
  height = 6,
  dpi = 300
)

# =========================
# 10. Resumen final
# =========================

cat("\n===========================\n")
cat("RESUMEN FINAL VALIDACIÓN EXTERNA\n")
cat("===========================\n")
cat("Tamaño firma final:", length(selected_genes), "\n")
cat("Best mtry final:", best_mtry, "\n")
cat("\nThresholds definidos en training/CV:\n")
print(threshold_table)
cat("\nMétricas por dataset y threshold:\n")
print(all_external_metrics)
cat("\nMétricas combinadas:\n")
print(combined_metrics)
cat("===========================\n")
