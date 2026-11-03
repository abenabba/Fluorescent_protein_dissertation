# =========================================================
# Raw inferred brightness benchmark
# File: brighter_flowseq_inferred_brightness_updated.csv
#
# Compares:
#   1) aa_length >= 150
#   2) aa_length >= 150 and barcode clusters >= 3
#   3) aa_length >= 150 and barcode clusters >= 5
#
# Models:
#   null_mean
#   support_only_lm
#   ridge_sequence_features
#   support_residual_plus_ridge
#   xgboost_sequence_only
#   xgboost_support_plus_sequence
# =========================================================

library(tidyverse)
library(glmnet)

has_xgb <- requireNamespace("xgboost", quietly = TRUE)

# -----------------------------
# Find raw brightness file
# -----------------------------
search_root <- "/Users/anissab/Desktop/illumina_shuffled_multi-bin_barcode_analysis"
target_file <- "brighter_flowseq_inferred_brightness_updated.csv"

matches <- list.files(
  search_root,
  pattern = target_file,
  recursive = TRUE,
  full.names = TRUE
)

print(matches)

if (length(matches) == 0) {
  stop("Could not find brighter_flowseq_inferred_brightness_updated.csv")
}

csv_path <- matches[1]
cat("Using brightness file:\n", csv_path, "\n")

df <- read_csv(csv_path, show_col_types = FALSE)

# -----------------------------
# Merge barcode-cluster support
# -----------------------------
support_path <- "/Users/anissab/Desktop/illumina_shuffled_multi-bin_barcode_analysis/barcodes_and_mutants/Bright_bc3_bc4_barcode_cluster_support.csv"

support_df <- read_csv(support_path, show_col_types = FALSE) %>%
  select(
    aa_seq,
    n_barcode_clusters_total
  )

df <- df %>%
  left_join(support_df, by = "aa_seq")

support_col <- "n_barcode_clusters_total"

cat("Merged support column:", support_col, "\n")
cat("Rows missing support after merge:", sum(is.na(df[[support_col]])), "\n")

# -----------------------------
# Identify columns / target
# -----------------------------
seq_col <- "aa_seq"

if (!seq_col %in% names(df)) {
  stop("Could not find aa_seq column.")
}

if ("log10_MEFL" %in% names(df)) {
  df <- df %>% mutate(target_log10_MEFL = log10_MEFL)
} else if ("median_MEFL" %in% names(df)) {
  df <- df %>% mutate(target_log10_MEFL = log10(median_MEFL))
} else {
  stop("Could not find log10_MEFL or median_MEFL.")
}

# -----------------------------
# Build sequence-derived features
# -----------------------------
make_seq_features <- function(df) {
  df %>%
    mutate(
      aa_seq = toupper(.data[[seq_col]]),
      aa_length = nchar(aa_seq),
      aromatic_frac = str_count(aa_seq, "[FWY]") / aa_length,
      charged_frac  = str_count(aa_seq, "[DEKRH]") / aa_length,
      polar_frac    = str_count(aa_seq, "[STNQCY]") / aa_length,
      gly_pro_frac  = str_count(aa_seq, "[GP]") / aa_length,
      hydrophobic_frac = str_count(aa_seq, "[AILMFWV]") / aa_length,
      has_gxy = str_detect(aa_seq, "G.Y"),
      has_syg = str_detect(aa_seq, "SYG"),
      has_tyg = str_detect(aa_seq, "TYG"),
      has_cyg = str_detect(aa_seq, "CYG"),
      has_fyg = str_detect(aa_seq, "FYG"),
      has_wyg = str_detect(aa_seq, "WYG")
    )
}

df_feat <- make_seq_features(df)

feature_cols <- c(
  "aa_length",
  "aromatic_frac",
  "charged_frac",
  "polar_frac",
  "gly_pro_frac",
  "hydrophobic_frac",
  "has_gxy",
  "has_syg",
  "has_tyg",
  "has_cyg",
  "has_fyg",
  "has_wyg"
)

# -----------------------------
# Metrics
# -----------------------------
calc_metrics <- function(actual, pred) {
  tibble(
    n_test = length(actual),
    r2 = 1 - sum((actual - pred)^2) / sum((actual - mean(actual))^2),
    rmse = sqrt(mean((actual - pred)^2)),
    pearson = suppressWarnings(cor(actual, pred, method = "pearson")),
    spearman = suppressWarnings(cor(actual, pred, method = "spearman"))
  )
}

# -----------------------------
# Benchmark function
# -----------------------------
run_benchmark <- function(df, filter_name, min_len = 150, min_support = NULL) {
  
  sub <- df %>%
    filter(
      aa_length >= min_len,
      !is.na(target_log10_MEFL),
      is.finite(target_log10_MEFL)
    )
  
  if (!is.null(min_support)) {
    sub <- sub %>%
      filter(.data[[support_col]] >= min_support)
  }
  
  sub <- sub %>%
    select(target_log10_MEFL, all_of(support_col), all_of(feature_cols)) %>%
    drop_na()
  
  cat("\n", filter_name, "\n")
  cat("Rows:", nrow(sub), "\n")
  
  set.seed(123)
  train_idx <- sample(seq_len(nrow(sub)), size = floor(0.8 * nrow(sub)))
  train <- sub[train_idx, ]
  test  <- sub[-train_idx, ]
  
  y_train <- train$target_log10_MEFL
  y_test  <- test$target_log10_MEFL
  
  # null mean
  pred_null <- rep(mean(y_train), length(y_test))
  null_res <- calc_metrics(y_test, pred_null) %>%
    mutate(filter = filter_name, model = "null_mean")
  
  # support-only linear model
  fit_support <- lm(
    as.formula(paste0("target_log10_MEFL ~ log10(", support_col, " + 1)")),
    data = train
  )
  pred_support <- predict(fit_support, newdata = test)
  
  support_res <- calc_metrics(y_test, pred_support) %>%
    mutate(filter = filter_name, model = "support_only_lm")
  
  # ridge on sequence-derived features
  x_train <- as.matrix(train %>% select(all_of(feature_cols)))
  x_test  <- as.matrix(test %>% select(all_of(feature_cols)))
  
  ridge_fit <- cv.glmnet(
    x_train,
    y_train,
    alpha = 0,
    nfolds = 5,
    standardize = TRUE
  )
  
  pred_ridge <- as.numeric(predict(ridge_fit, newx = x_test, s = "lambda.min"))
  
  ridge_res <- calc_metrics(y_test, pred_ridge) %>%
    mutate(filter = filter_name, model = "ridge_sequence_features")
  
  # support residual + ridge
  train$resid_support <- residuals(fit_support)
  test$expected_support <- predict(fit_support, newdata = test)
  test$resid_support_actual <- test$target_log10_MEFL - test$expected_support
  
  resid_fit <- cv.glmnet(
    x_train,
    train$resid_support,
    alpha = 0,
    nfolds = 5,
    standardize = TRUE
  )
  
  pred_resid <- as.numeric(predict(resid_fit, newx = x_test, s = "lambda.min"))
  pred_total <- test$expected_support + pred_resid
  
  resid_res <- calc_metrics(y_test, pred_total) %>%
    mutate(filter = filter_name, model = "support_residual_plus_ridge")
  
  # XGBoost sequence-only
  if (has_xgb) {
    xgb_train_seq <- as.matrix(train %>% select(all_of(feature_cols)))
    xgb_test_seq  <- as.matrix(test %>% select(all_of(feature_cols)))
    
    dtrain_seq <- xgboost::xgb.DMatrix(data = xgb_train_seq, label = y_train)
    dtest_seq  <- xgboost::xgb.DMatrix(data = xgb_test_seq, label = y_test)
    
    xgb_seq_fit <- xgboost::xgb.train(
      params = list(
        objective = "reg:squarederror",
        max_depth = 3,
        eta = 0.03,
        subsample = 0.8,
        colsample_bytree = 0.8
      ),
      data = dtrain_seq,
      nrounds = 100,
      verbose = 0
    )
    
    pred_xgb_seq <- predict(xgb_seq_fit, dtest_seq)
    
    xgb_seq_res <- calc_metrics(y_test, pred_xgb_seq) %>%
      mutate(filter = filter_name, model = "xgboost_sequence_only")
  } else {
    xgb_seq_res <- tibble(
      n_test = nrow(test),
      r2 = NA_real_,
      rmse = NA_real_,
      pearson = NA_real_,
      spearman = NA_real_,
      filter = filter_name,
      model = "xgboost_sequence_only_not_installed"
    )
  }
  
  # XGBoost support + sequence
  if (has_xgb) {
    xgb_train_all <- as.matrix(train %>% select(all_of(c(support_col, feature_cols))))
    xgb_test_all  <- as.matrix(test %>% select(all_of(c(support_col, feature_cols))))
    
    dtrain_all <- xgboost::xgb.DMatrix(data = xgb_train_all, label = y_train)
    dtest_all  <- xgboost::xgb.DMatrix(data = xgb_test_all, label = y_test)
    
    xgb_all_fit <- xgboost::xgb.train(
      params = list(
        objective = "reg:squarederror",
        max_depth = 3,
        eta = 0.03,
        subsample = 0.8,
        colsample_bytree = 0.8
      ),
      data = dtrain_all,
      nrounds = 100,
      verbose = 0
    )
    
    pred_xgb_all <- predict(xgb_all_fit, dtest_all)
    
    xgb_all_res <- calc_metrics(y_test, pred_xgb_all) %>%
      mutate(filter = filter_name, model = "xgboost_support_plus_sequence")
  } else {
    xgb_all_res <- tibble(
      n_test = nrow(test),
      r2 = NA_real_,
      rmse = NA_real_,
      pearson = NA_real_,
      spearman = NA_real_,
      filter = filter_name,
      model = "xgboost_support_plus_sequence_not_installed"
    )
  }
  
  bind_rows(
    null_res,
    support_res,
    ridge_res,
    resid_res,
    xgb_seq_res,
    xgb_all_res
  )
}

# -----------------------------
# Run all filters
# -----------------------------
summary_tbl <- bind_rows(
  run_benchmark(
    df_feat,
    filter_name = "aa_length >= 150",
    min_len = 150,
    min_support = NULL
  ),
  run_benchmark(
    df_feat,
    filter_name = "aa_length >= 150 and support >= 3",
    min_len = 150,
    min_support = 3
  ),
  run_benchmark(
    df_feat,
    filter_name = "aa_length >= 150 and support >= 5",
    min_len = 150,
    min_support = 5
  )
) %>%
  select(filter, model, n_test, r2, rmse, pearson, spearman) %>%
  mutate(
    r2 = round(r2, 3),
    rmse = round(rmse, 3),
    pearson = round(pearson, 3),
    spearman = round(spearman, 3)
  ) %>%
  arrange(filter, desc(r2))

print(summary_tbl, n = Inf)

# -----------------------------
# Plot R2 by model/filter
# -----------------------------
p <- summary_tbl %>%
  filter(model != "null_mean") %>%
  ggplot(aes(x = model, y = r2, fill = filter)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "Raw inferred brightness benchmark",
    subtitle = "Test-set R² across filtering strategies",
    x = "Model",
    y = "Test R²"
  ) +
  theme_bw(base_size = 12)

print(p)

ggsave(
  filename = file.path(dirname(csv_path), "raw_brightness_model_benchmark_r2.png"),
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)