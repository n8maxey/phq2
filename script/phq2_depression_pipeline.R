# ============================================================
# PHQ-2 DEPRESSION: REGULAR vs EMBEDDINGS vs COMBINED MODELS
# CALIBRATION CURVES & ROC PLOTS (BLACK & WHITE, 3x3 GRIDS)
# ============================================================
# SURGICAL FIXES APPLIED (v9):
#   FIX 20 — Brier scores + 95% bootstrap CI added to calibration plots
#             and saved as CSVs (test + train)
# ============================================================

# ============================================================
# SET PERMANENT LIBRARY PATH (CRITICAL FOR HPC)
# ============================================================
user_lib <- "/N/scratch/natmaxey/R_libs"
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(user_lib, .libPaths()))
cat("R library paths:\n")
cat("  ", paste(.libPaths(), collapse = "\n   "), "\n\n")

local({
  r <- getOption("repos")
  r["CRAN"] <- "https://cloud.r-project.org"
  options(repos = r)
})
options(repos = c(CRAN = "https://cloud.r-project.org"))

packages <- c(
  "readr", "dplyr", "stringr", "purrr", "tidyr",
  "parsnip", "workflows", "recipes", "rsample", "yardstick",
  "tune", "dials", "ranger", "xgboost", "glmnet",
  "ggplot2", "tibble", "doParallel", "gtsummary", "flextable", "officer",
  "patchwork", "pROC", "ggpubr"
)

cat("Loading packages...\n")
for (pkg in packages) {
  if (!suppressWarnings(require(pkg, character.only = TRUE, lib.loc = user_lib))) {
    cat("Installing ", pkg, " to ", user_lib, "...\n")
    install.packages(pkg, lib = user_lib, repos = "https://cloud.r-project.org", quiet = FALSE)
    suppressWarnings(require(pkg, character.only = TRUE, lib.loc = user_lib))
  }
}
cat("All packages loaded\n\n")

# ============================================================
# PATHS
# ============================================================
is_quartz <- Sys.getenv("QUARTZ_JOB", unset = "FALSE") == "TRUE"
if (is_quartz) {
  base_dir <- "/N/u/natmaxey/Quartz/phq2_depression_ml"
  cat("Running on Quartz HPC\n\n")
} else {
  base_dir <- "/N/u/natmaxey/Quartz/phq2_depression_ml"
  cat("Running on local PC\n\n")
}

data_dir    <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")
extract_dir <- file.path(base_dir, "extracted_csvs")

dir.create(data_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

phq2_clean_rds_path <- file.path(data_dir, "hps_cleaned_depression.rds")
phq2_clean_csv_path <- file.path(data_dir, "hps_cleaned_depression.csv")

cat("Base Dir:    ", base_dir,    "\n")
cat("Data Dir:    ", data_dir,    "\n")
cat("Results Dir: ", results_dir, "\n")
cat("Extract Dir: ", extract_dir, "\n\n")

y_var  <- "depression"
x_vars <- c("age", "sex", "race", "hispanic", "education", "marital", "numkids")

# Checkpoints
checkpoint_table1          <- file.path(results_dir, "Table1_Univariate_Statistics.docx")
checkpoint_reliability     <- file.path(results_dir, "phq2_reliability_by_cycle.csv")
checkpoint_ml_df           <- file.path(results_dir, "ml_dataset_processed.rds")
checkpoint_train_test      <- file.path(results_dir, "train_test_split.rds")
checkpoint_tfidf_train     <- file.path(results_dir, "tfidf_train_features.rds")
checkpoint_tfidf_test      <- file.path(results_dir, "tfidf_test_features.rds")
checkpoint_cal_plots       <- file.path(results_dir, "checkpoint_calibration_plots.txt")
checkpoint_roc_plots       <- file.path(results_dir, "checkpoint_roc_plots.txt")
checkpoint_shap            <- file.path(results_dir, "checkpoint_shap_plots.txt")

# ============================================================
# PHQ-2 DETECTION & PROCESS FUNCTIONS
# ============================================================
find_adult_phq2_items <- function(df_names) {
  item1_var <- if ("INTEREST" %in% df_names) "INTEREST" else NULL
  item2_var <- if ("DOWN"     %in% df_names) "DOWN"     else NULL
  list(item1 = item1_var, item2 = item2_var)
}

process_hps_data <- function(df_raw) {
  phq2_detection <- find_adult_phq2_items(names(df_raw))
  if (is.null(phq2_detection$item1) || is.null(phq2_detection$item2)) return(NULL)
  
  age_col     <- names(df_raw)[names(df_raw) %in% c("TBIRTH_YEAR","TAGE1","TAGE")][1]
  sex_col     <- names(df_raw)[names(df_raw) %in% c("EGENID_BIRTH","ESEX1","EGENDER")][1]
  educ_col    <- names(df_raw)[names(df_raw) %in% c("EEDUC","REDUC1")][1]
  marital_col <- names(df_raw)[names(df_raw) %in% c("MS","MARITAL1")][1]
  numkids_col <- names(df_raw)[names(df_raw) %in% c("THHLD_NUMKID","AHHLD_NUMKID")][1]
  weight_col  <- names(df_raw)[names(df_raw) %in% c("PWEIGHT","WEIGHT")][1]
  
  if (is.na(sex_col)) return(NULL)
  df <- df_raw %>% mutate(across(everything(), as.character))
  
  df <- df %>%
    mutate(
      phq2_item1 = suppressWarnings(as.numeric(.data[[phq2_detection$item1]])) - 1,
      phq2_item2 = suppressWarnings(as.numeric(.data[[phq2_detection$item2]])) - 1,
      phq2_item1 = pmax(0, pmin(3, phq2_item1)),
      phq2_item2 = pmax(0, pmin(3, phq2_item2)),
      phq2_score = phq2_item1 + phq2_item2,
      depression = if_else(phq2_score >= 3, 1, 0, missing = NA_integer_)
    ) %>%
    filter(!is.na(depression))
  
  df <- df %>%
    mutate(
      age = if (!is.na(age_col)) {
        if (age_col == "TBIRTH_YEAR") {
          2025 - suppressWarnings(as.numeric(.data[[age_col]]))
        } else {
          suppressWarnings(as.numeric(.data[[age_col]]))
        }
      } else NA_real_,
      
      sex = if (!is.na(sex_col))
        factor(case_when(
          .data[[sex_col]] == "1" ~ "Male",
          .data[[sex_col]] == "2" ~ "Female",
          TRUE ~ NA_character_
        ), levels = c("Male","Female"))
      else NA_character_,
      
      education = if (!is.na(educ_col))
        factor(case_when(
          .data[[educ_col]] %in% c("1","2","3") ~ "HS or less",
          .data[[educ_col]] %in% c("4","5")     ~ "Some college",
          .data[[educ_col]] %in% c("6","7")     ~ "Bachelor+",
          TRUE ~ NA_character_
        ), levels = c("Bachelor+","Some college","HS or less"))
      else NA_character_,
      
      marital = if (!is.na(marital_col))
        factor(case_when(
          .data[[marital_col]] == "1" ~ "Married",
          .data[[marital_col]] == "2" ~ "Widowed",
          .data[[marital_col]] == "3" ~ "Divorced",
          .data[[marital_col]] == "4" ~ "Separated",
          .data[[marital_col]] == "5" ~ "Never married",
          TRUE ~ NA_character_
        ), levels = c("Married","Widowed","Divorced","Separated","Never married"))
      else NA_character_,
      
      race = case_when(
        RRACE == "1" ~ "White",
        RRACE == "2" ~ "Black",
        RRACE == "3" ~ "Asian",
        RRACE == "4" ~ "Other",
        TRUE         ~ "Other"
      ),
      
      hispanic = case_when(
        RHISPANIC == "2" ~ "Yes",
        RHISPANIC == "1" ~ "No",
        TRUE             ~ NA_character_
      ),
      
      numkids = if (!is.na(numkids_col)) suppressWarnings(as.numeric(.data[[numkids_col]])) else NA_real_,
      weight  = if (!is.na(weight_col))  suppressWarnings(as.numeric(.data[[weight_col]]))  else NA_real_
    ) %>%
    select(any_of(c("phq2_item1","phq2_item2","phq2_score","depression",
                    "age","sex","race","hispanic","education","marital","numkids","weight")))
  return(df)
}

load_hps_csvs <- function(extract_dir) {
  csvs <- list.files(extract_dir, pattern = "hps_04_.*puf\\.csv$",
                     full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
  csvs <- csvs[!grepl("repwgt", csvs, ignore.case = TRUE)]
  if (length(csvs) == 0) { cat("No CSVs found in: ", extract_dir, "\n"); return(NULL) }
  cat("Found ", length(csvs), " CSV(s)\n\n")
  
  cycle_map <- list(
    "hps_04_00_01" = 1, "hps_04_00_02" = 2, "hps_04_00_03" = 3,
    "hps_04_01_04" = 4, "hps_04_01_05" = 5, "hps_04_01_06" = 6,
    "hps_04_01_07" = 7, "hps_04_02_08" = 8, "hps_04_02_09" = 9
  )
  out <- list()
  
  for (i in seq_along(csvs)) {
    f <- csvs[i]; basename_f <- basename(f)
    message("[", i, "/", length(csvs), "] ", basename_f)
    
    df_raw <- tryCatch(
      readr::read_csv(f, col_types = readr::cols(.default = "c"), show_col_types = FALSE),
      error = function(e) { message("  Read error"); NULL }
    )
    if (is.null(df_raw) || nrow(df_raw) == 0) { message("  Empty"); next }
    
    clean <- tryCatch(process_hps_data(df_raw),
                      error = function(e) { message("  Process error"); NULL })
    if (!is.null(clean) && nrow(clean) > 0) {
      cycle <- NA_integer_
      for (key in names(cycle_map)) if (grepl(key, basename_f)) { cycle <- cycle_map[[key]]; break }
      clean$WEEK <- cycle
      out[[length(out)+1]] <- clean
      message("  ", nrow(clean), " rows (CYCLE=", cycle, ")")
    }
  }
  if (length(out) == 0) { cat("No data extracted\n"); return(NULL) }
  cat("\nCombining ", length(out), " files...\n")
  dplyr::bind_rows(out)
}

# ============================================================
# MAIN: LOAD DATA
# ============================================================
cat("Loading HPS data...\n\n")
hps_cleaned <- load_hps_csvs(extract_dir)
if (is.null(hps_cleaned) || nrow(hps_cleaned) == 0) stop("No data loaded")
cat("\nLoaded ", nrow(hps_cleaned), " rows\n")
cat("Depression: ", sum(hps_cleaned$depression),
    " (", round(100 * mean(hps_cleaned$depression), 1), "%)\n\n")
saveRDS(hps_cleaned, phq2_clean_rds_path)
readr::write_csv(hps_cleaned, phq2_clean_csv_path)

# ============================================================
# RELIABILITY BY CYCLE
# ============================================================
cat("Computing reliability by cycle...\n\n")
reliability_by_cycle <- hps_cleaned %>%
  filter(!is.na(WEEK)) %>%
  group_by(WEEK) %>%
  summarise(
    n     = sum(!is.na(phq2_item1) & !is.na(phq2_item2)),
    r12   = suppressWarnings(cor(phq2_item1, phq2_item2, use = "complete.obs")),
    omega = 2 * r12 / (1 + r12),
    .groups = "drop"
  )
print(reliability_by_cycle)
readr::write_csv(reliability_by_cycle, file.path(results_dir, "phq2_reliability_by_cycle.csv"))
cat("Saved reliability\n\n")

# ============================================================
# MISSINGNESS TABLE
# ============================================================
cat("Computing missingness by variable...\n\n")
checkpoint_missingness <- file.path(results_dir, "missingness_by_variable.csv")
if (!file.exists(checkpoint_missingness)) {
  missingness_summary <- hps_cleaned %>%
    dplyr::summarise(dplyr::across(dplyr::everything(),
                                   ~sum(is.na(.)) / dplyr::n() * 100)) %>%
    tidyr::pivot_longer(cols = dplyr::everything(),
                        names_to  = "variable",
                        values_to = "percent_missing") %>%
    dplyr::arrange(dplyr::desc(percent_missing)) %>%
    dplyr::mutate(
      n_total         = nrow(hps_cleaned),
      n_missing       = round(percent_missing * n_total / 100),
      n_complete      = n_total - n_missing,
      percent_missing = round(percent_missing, 2)
    ) %>%
    dplyr::select(variable, n_total, n_missing, n_complete, percent_missing)
  readr::write_csv(missingness_summary, checkpoint_missingness)
  cat("Saved missingness table\n\n")
} else {
  cat("Missingness table already exists\n\n")
}

# ============================================================
# ML DATASET
# ============================================================
if (file.exists(checkpoint_ml_df)) {
  cat("Loading ML dataset from checkpoint...\n\n")
  ml_df <- readRDS(checkpoint_ml_df)
  if (!"phq2_score" %in% names(ml_df)) {
    cat("  Checkpoint missing phq2_score - rebuilding...\n\n")
    ml_df <- hps_cleaned %>%
      select(all_of(c("WEEK", y_var, "phq2_score", x_vars, "weight"))) %>%
      filter(!is.na(.data[[y_var]])) %>%
      filter(!if_any(all_of(c(y_var, x_vars)), is.na))
    saveRDS(ml_df, checkpoint_ml_df)
  }
} else {
  cat("Creating ML dataset...\n\n")
  ml_df <- hps_cleaned %>%
    select(all_of(c("WEEK", y_var, "phq2_score", x_vars, "weight"))) %>%
    filter(!is.na(.data[[y_var]]))
  ml_df_cca <- ml_df %>%
    filter(!if_any(all_of(c(y_var, x_vars)), is.na))
  cat("ML dataset: ", nrow(ml_df), " -> after CCA: ", nrow(ml_df_cca), " rows\n\n")
  ml_df <- ml_df_cca
  saveRDS(ml_df, checkpoint_ml_df)
  cat("Saved ML dataset checkpoint\n\n")
}

# ============================================================
# TRAIN/TEST SPLIT
# ============================================================
cat("Creating train/test split...\n")
set.seed(67)
split    <- rsample::initial_split(ml_df, prop = 0.80, strata = depression)
train_df <- rsample::training(split)
test_df  <- rsample::testing(split)
saveRDS(list(train = train_df, test = test_df), file.path(results_dir, "train_test_split.rds"))
train_df[[y_var]] <- factor(train_df[[y_var]], levels = c(0, 1))
test_df[[y_var]]  <- factor(test_df[[y_var]],  levels = c(0, 1))
cat("Train: ", nrow(train_df), " | Test: ", nrow(test_df), "\n\n")

# ============================================================
# TF-IDF EMBEDDINGS
# ============================================================
cat("Setting up TF-IDF function...\n\n")
create_tfidf_features <- function(df) {
  cat("  Input df: ", nrow(df), " rows\n")
  df_with_doc <- df %>% dplyr::mutate(doc_id = dplyr::row_number())
  cat("  After adding doc_id: ", nrow(df_with_doc), " rows\n")
  
  terms_long <- df_with_doc %>%
    dplyr::select(doc_id, sex, race, hispanic, education, marital) %>%
    tidyr::pivot_longer(cols = -doc_id, names_to = "feature",
                        values_to = "term", values_drop_na = TRUE) %>%
    dplyr::mutate(term_name = paste0(feature, "_", term)) %>%
    dplyr::select(doc_id, term_name) %>%
    dplyr::distinct()
  cat("  Terms created\n")
  
  total_docs <- dplyr::n_distinct(df_with_doc$doc_id)
  idf_vals <- terms_long %>%
    dplyr::group_by(term_name) %>%
    dplyr::summarise(doc_freq = dplyr::n(),
                     idf = log(total_docs / (1 + doc_freq)), .groups = "drop")
  
  tfidf_long <- terms_long %>%
    dplyr::group_by(doc_id) %>%
    dplyr::mutate(tf = 1 / dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(idf_vals, by = "term_name") %>%
    dplyr::mutate(tf_idf = tf * idf)
  
  top_features <- tfidf_long %>%
    dplyr::group_by(term_name) %>%
    dplyr::summarise(mean_tfidf = mean(tf_idf, na.rm = TRUE), .groups = "drop") %>%
    dplyr::slice_max(order_by = mean_tfidf, n = 10) %>%
    dplyr::pull(term_name)
  cat("  Top 10 features selected\n")
  
  embedding_wide <- tfidf_long %>%
    dplyr::filter(term_name %in% top_features) %>%
    tidyr::pivot_wider(id_cols = doc_id, names_from = term_name,
                       values_from = tf_idf, values_fill = 0) %>%
    dplyr::arrange(doc_id)
  
  all_rows <- tibble::tibble(doc_id = 1:total_docs)
  embedding_wide <- all_rows %>%
    dplyr::left_join(embedding_wide, by = "doc_id") %>%
    dplyr::mutate(dplyr::across(-doc_id, ~tidyr::replace_na(., 0))) %>%
    dplyr::select(-doc_id)
  cat("  Final: ", nrow(embedding_wide), " rows\n")
  embedding_wide
}

cat("Setting up TF-IDF embeddings...\n\n")
if (file.exists(checkpoint_tfidf_train)) {
  cat("Loading train TF-IDF features from checkpoint...\n")
  tfidf_train <- readRDS(checkpoint_tfidf_train)
} else {
  cat("Creating train TF-IDF features...\n")
  tfidf_train <- create_tfidf_features(train_df)
  saveRDS(tfidf_train, checkpoint_tfidf_train)
}
cat("  Train TF-IDF: ", nrow(tfidf_train), " rows, ", ncol(tfidf_train), " features\n\n")

if (file.exists(checkpoint_tfidf_test)) {
  cat("Loading test TF-IDF features from checkpoint...\n")
  tfidf_test <- readRDS(checkpoint_tfidf_test)
} else {
  cat("Creating test TF-IDF features...\n")
  tfidf_test <- create_tfidf_features(test_df)
  saveRDS(tfidf_test, checkpoint_tfidf_test)
}
cat("  Test TF-IDF: ", nrow(tfidf_test), " rows, ", ncol(tfidf_test), " features\n\n")

tfidf_cols <- names(tfidf_train)
for (col in tfidf_cols) {
  if (!col %in% names(tfidf_test)) tfidf_test[[col]] <- 0
}
tfidf_test <- tfidf_test %>% dplyr::select(dplyr::all_of(tfidf_cols))
cat("TF-IDF features ready\n\n")

# ============================================================
# BUILD TF-IDF / COMBINED DATASETS
# ============================================================
train_df_tfidf <- dplyr::bind_cols(
  train_df %>% dplyr::select(all_of(c(y_var, "phq2_score", "age", "numkids", "weight"))),
  tfidf_train
)
test_df_tfidf <- dplyr::bind_cols(
  test_df %>% dplyr::select(all_of(c(y_var, "phq2_score", "age", "numkids", "weight"))),
  tfidf_test
)
train_df_combined <- train_df_tfidf
test_df_combined  <- test_df_tfidf

cat("TF-IDF Dataset Check:\n")
cat("  train_df rows:       ", nrow(train_df),       "\n")
cat("  train_df_tfidf rows: ", nrow(train_df_tfidf), "\n")
cat("  test_df rows:        ", nrow(test_df),        "\n")
cat("  test_df_tfidf rows:  ", nrow(test_df_tfidf),  "\n\n")

# ============================================================
# CLASS IMBALANCE ANALYSIS
# ============================================================
cat("Analyzing class imbalance...\n\n")
class_dist      <- table(train_df$depression)
n_neg           <- as.numeric(class_dist["0"])
n_pos           <- as.numeric(class_dist["1"])
total_n         <- n_neg + n_pos
prevalence      <- n_pos / total_n
imbalance_ratio <- n_neg / n_pos

cat("Class Distribution in Training Set:\n")
cat("  No Depression (0): ", n_neg, " (", round((n_neg/total_n)*100, 1), "%)\n", sep="")
cat("  Depression (1):    ", n_pos, " (", round((n_pos/total_n)*100, 1), "%)\n", sep="")
cat("  Imbalance Ratio:   ", round(imbalance_ratio, 2), ":1\n\n", sep="")

weight_neg_balanced <- 1.0
weight_pos_balanced <- imbalance_ratio
scale_pos_weight    <- imbalance_ratio

cat("Cost-Sensitive Weights (RF & XGB only):\n")
cat("  Class 0 weight: ",           round(weight_neg_balanced, 3), "\n", sep="")
cat("  Class 1 weight: ",           round(weight_pos_balanced, 3), "\n", sep="")
cat("  XGBoost scale_pos_weight: ", round(scale_pos_weight, 3),    "\n")
cat("  NOTE: Logistic regression uses NO case weights (glmnet incompatible)\n\n")

class_imbalance_summary <- tibble::tibble(
  class                = c("No Depression","Depression","Total"),
  n                    = c(n_neg, n_pos, total_n),
  percent              = c(round((n_neg/total_n)*100,2), round((n_pos/total_n)*100,2), 100.00),
  imbalance_ratio      = c(NA, imbalance_ratio, NA),
  balanced_weight      = c(weight_neg_balanced, weight_pos_balanced, NA),
  xgb_scale_pos_weight = c(NA, scale_pos_weight, NA)
)
readr::write_csv(class_imbalance_summary, file.path(results_dir, "class_imbalance_analysis.csv"))
cat("Saved class imbalance analysis\n\n")

# ============================================================
# MODEL SPECS
# ============================================================
cat("=== SETTING UP MODEL SPECS ===\n\n")
n_cores <- 20L
cat("Using ", n_cores, " cores (fixed)\n\n")
doParallel::registerDoParallel(cores = n_cores)

rf_spec <- parsnip::rand_forest(
  mtry  = tune::tune(),
  min_n = tune::tune(),
  trees = 500
) %>%
  parsnip::set_engine("ranger",
                      num.threads   = n_cores,
                      importance    = "impurity",
                      class.weights = c("0" = weight_neg_balanced,
                                        "1" = weight_pos_balanced)) %>%
  parsnip::set_mode("classification")

xgb_spec <- parsnip::boost_tree(
  trees          = tune::tune(),
  tree_depth     = tune::tune(),
  min_n          = tune::tune(),
  learn_rate     = tune::tune(),
  loss_reduction = tune::tune()
) %>%
  parsnip::set_engine("xgboost",
                      nthread          = n_cores,
                      scale_pos_weight = scale_pos_weight) %>%
  parsnip::set_mode("classification")

logit_spec <- parsnip::logistic_reg(
  penalty = tune::tune(),
  mixture = tune::tune()
) %>%
  parsnip::set_engine("glmnet") %>%
  parsnip::set_mode("classification")

cat("Model specs defined (RF, XGB, Logistic Regression)\n\n")

# ============================================================
# CLASSIFICATION TRAINING HELPER
# ============================================================
grid_size_for <- function(wf) {
  spec_class <- class(workflows::extract_spec_parsnip(wf))[1]
  if (spec_class == "boost_tree")  return(50L)
  if (spec_class == "rand_forest") return(25L)
  return(25L)
}

train_class_quick <- function(wf, train_data, name) {
  cat("Training ", name, "...\n")
  
  param_set <- workflows::extract_parameter_set_dials(wf)
  param_set <- dials::finalize(param_set, train_data)
  
  grid_size <- grid_size_for(wf)
  set.seed(42)
  lhc_grid <- dials::grid_space_filling(
    param_set,
    size = grid_size,
    type = "latin_hypercube"
  )
  cat("  grid_space_filling size = ", grid_size,
      " | params = ", nrow(param_set), "\n", sep = "")
  
  folds <- rsample::vfold_cv(train_data, v = 5, strata = depression)
  
  tune_result <- tune::tune_grid(
    wf,
    resamples = folds,
    grid      = lhc_grid,
    metrics   = yardstick::metric_set(yardstick::roc_auc),
    control   = tune::control_grid(
      verbose       = FALSE,
      parallel_over = "everything"
    )
  )
  
  best <- tune::select_best(tune_result, metric = "roc_auc")
  tune::finalize_workflow(wf, best) %>% parsnip::fit(train_data)
}

# ============================================================
# RECIPES
# ============================================================
rec <- recipes::recipe(
  depression ~ age + sex + race + hispanic + education + marital + numkids,
  data = train_df
) %>%
  recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = TRUE) %>%
  recipes::step_zv(recipes::all_predictors()) %>%
  recipes::step_normalize(recipes::all_numeric_predictors())

tfidf_term_str <- paste(sprintf("`%s`", tfidf_cols), collapse = " + ")
rec_tfidf <- recipes::recipe(
  as.formula(paste("depression ~ age + numkids +", tfidf_term_str)),
  data = train_df_tfidf
) %>%
  recipes::step_zv(recipes::all_predictors()) %>%
  recipes::step_normalize(recipes::all_numeric_predictors())

rec_combined <- recipes::recipe(
  as.formula(paste("depression ~ age + numkids + weight +", tfidf_term_str)),
  data = train_df_combined
) %>%
  recipes::step_zv(recipes::all_predictors()) %>%
  recipes::step_normalize(recipes::all_numeric_predictors())

# ============================================================
# BUILD WORKFLOWS
# ============================================================
wf_rf        <- workflows::workflow() %>% workflows::add_model(rf_spec)    %>% workflows::add_recipe(rec)
wf_xgb       <- workflows::workflow() %>% workflows::add_model(xgb_spec)   %>% workflows::add_recipe(rec)
wf_logit_reg <- workflows::workflow() %>% workflows::add_model(logit_spec) %>% workflows::add_recipe(rec)

wf_rf_tfidf    <- workflows::workflow() %>% workflows::add_model(rf_spec)    %>% workflows::add_recipe(rec_tfidf)
wf_xgb_tfidf   <- workflows::workflow() %>% workflows::add_model(xgb_spec)   %>% workflows::add_recipe(rec_tfidf)
wf_logit_tfidf <- workflows::workflow() %>% workflows::add_model(logit_spec) %>% workflows::add_recipe(rec_tfidf)

wf_rf_combined    <- workflows::workflow() %>% workflows::add_model(rf_spec)    %>% workflows::add_recipe(rec_combined)
wf_xgb_combined   <- workflows::workflow() %>% workflows::add_model(xgb_spec)   %>% workflows::add_recipe(rec_combined)
wf_logit_combined <- workflows::workflow() %>% workflows::add_model(logit_spec) %>% workflows::add_recipe(rec_combined)

# ============================================================
# TRAIN 9 CLASSIFICATION MODELS
# ============================================================
cat("=== TRAINING CLASSIFICATION MODELS ===\n\n")

final_rf_class_reg    <- train_class_quick(wf_rf,        train_df, "RF (Regular)")
final_xgb_class_reg   <- train_class_quick(wf_xgb,       train_df, "XGB (Regular)")
final_logit_class_reg <- train_class_quick(wf_logit_reg, train_df, "Logit (Regular)")

final_rf_class_tfidf    <- train_class_quick(wf_rf_tfidf,    train_df_tfidf, "RF (TF-IDF)")
final_xgb_class_tfidf   <- train_class_quick(wf_xgb_tfidf,   train_df_tfidf, "XGB (TF-IDF)")
final_logit_class_tfidf <- train_class_quick(wf_logit_tfidf, train_df_tfidf, "Logit (TF-IDF)")

final_rf_class_combined    <- train_class_quick(wf_rf_combined,    train_df_combined, "RF (Combined)")
final_xgb_class_combined   <- train_class_quick(wf_xgb_combined,   train_df_combined, "XGB (Combined)")
final_logit_class_combined <- train_class_quick(wf_logit_combined, train_df_combined, "Logit (Combined)")

cat("All 9 classification models trained (RF + XGB + Logistic Regression x 3 feature sets)\n\n")

# ============================================================
# EVALUATION HELPERS
# ============================================================
cat("============================================================\n")
cat("EVALUATING ALL 9 CLASSIFICATION MODELS\n")
cat("============================================================\n\n")

compute_boot_ci <- function(pred_prob, pred_binary, truth, B = 500, seed = 42) {
  set.seed(seed)
  n <- length(truth)
  
  boot_one <- function(idx) {
    pb  <- pred_prob[idx]; bb  <- pred_binary[idx]; tb  <- truth[idx]
    pos <- which(tb == 1); neg <- which(tb == 0)
    if (length(pos) == 0 || length(neg) == 0) return(rep(NA_real_, 6))
    auc_b <- mean(outer(pb[pos], pb[neg], ">")) +
      0.5 * mean(outer(pb[pos], pb[neg], "=="))
    tp <- sum(bb == 1 & tb == 1); tn <- sum(bb == 0 & tb == 0)
    fp <- sum(bb == 1 & tb == 0); fn <- sum(bb == 0 & tb == 1)
    c(auc_b,
      (tp + tn) / n,
      if ((tp+fn)>0) tp/(tp+fn) else NA_real_,
      if ((tn+fp)>0) tn/(tn+fp) else NA_real_,
      if ((tp+fp)>0) tp/(tp+fp) else NA_real_,
      if ((tn+fn)>0) tn/(tn+fn) else NA_real_)
  }
  mat <- replicate(B, boot_one(sample(n, n, replace = TRUE)))
  list(
    lo = apply(mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    hi = apply(mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
}

# ============================================================
# BOOTSTRAP CI FOR YOUDEN THRESHOLD (NEW)
# ============================================================
bootstrap_youden_ci <- function(prob, truth, B = 1000, seed = 42) {
  set.seed(seed)
  n <- length(truth)
  
  boot_thresh <- replicate(B, {
    idx <- sample(n, n, replace = TRUE)
    
    roc_b <- tryCatch(
      pROC::roc(response = truth[idx], predictor = prob[idx], quiet = TRUE),
      error = function(e) NULL
    )
    
    if (is.null(roc_b)) return(NA_real_)
    
    t <- tryCatch(
      pROC::coords(roc_b, x = "best", ret = "threshold", best.method = "youden")$threshold[1],
      error = function(e) NA_real_
    )
    
    if (!is.finite(t)) NA_real_ else t
  })
  
  list(
    lower = as.numeric(quantile(boot_thresh, 0.025, na.rm = TRUE)),
    upper = as.numeric(quantile(boot_thresh, 0.975, na.rm = TRUE))
  )
}

eval_class <- function(model, train_data, test_data, model_name, B = 1000) {
  train_prob  <- predict(model, train_data, type = "prob")$.pred_1
  train_truth <- as.numeric(as.character(train_data[[y_var]]))
  
  roc_train  <- pROC::roc(response = train_truth, predictor = train_prob, quiet = TRUE)
  thresh_raw <- pROC::coords(roc_train, x = "best", ret = "threshold",
                             best.method = "youden")$threshold
  opt_thresh <- thresh_raw[1]
  if (!is.finite(opt_thresh)) opt_thresh <- 0.5
  # --- NEW: bootstrap CI for Youden threshold ---
  thresh_ci <- bootstrap_youden_ci(train_prob, train_truth, B = 1000)
  truth       <- factor(test_data[[y_var]], levels = c(0, 1))
  pred_prob   <- predict(model, test_data, type = "prob")$.pred_1
  pred_binary <- factor(ifelse(pred_prob >= opt_thresh, 1, 0), levels = c(0, 1))
  
  auc_val <- yardstick::roc_auc(
    tibble::tibble(truth = truth, .pred_1 = pred_prob),
    truth, .pred_1, event_level = "second"
  )$.estimate
  
  tp <- sum(pred_binary == 1 & truth == 1)
  tn <- sum(pred_binary == 0 & truth == 0)
  fp <- sum(pred_binary == 1 & truth == 0)
  fn <- sum(pred_binary == 0 & truth == 1)
  n  <- length(truth)
  
  acc_val  <- (tp + tn) / n
  sens_val <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  spec_val <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  ppv_val  <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  npv_val  <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
  
  cat("  Computing bootstrap CIs for", model_name, "(B=1000)...\n")
  ci <- compute_boot_ci(pred_prob, pred_binary, truth, B = B)
  
  cat(sprintf(paste0(
    "  %s (thresh=%.3f [%.3f, %.3f], from train):\n",
    "    AUC  = %.3f [%.3f, %.3f]\n",
    "    Acc  = %.3f [%.3f, %.3f]\n",
    "    Sens = %.3f [%.3f, %.3f]\n",
    "    Spec = %.3f [%.3f, %.3f]\n",
    "    PPV  = %.3f [%.3f, %.3f]\n",
    "    NPV  = %.3f [%.3f, %.3f]\n\n"),
    model_name, opt_thresh, thresh_ci$lower, thresh_ci$upper,
    auc_val,  ci$lo[1], ci$hi[1],
    acc_val,  ci$lo[2], ci$hi[2],
    replace(sens_val, is.na(sens_val), NaN), ci$lo[3], ci$hi[3],
    replace(spec_val, is.na(spec_val), NaN), ci$lo[4], ci$hi[4],
    replace(ppv_val,  is.na(ppv_val),  NaN), ci$lo[5], ci$hi[5],
    replace(npv_val,  is.na(npv_val),  NaN), ci$lo[6], ci$hi[6]
  ))
  
  tibble::tibble(
    model             = model_name,
    optimal_threshold = opt_thresh,
    thresh_lo         = thresh_ci$lower,
    thresh_hi         = thresh_ci$upper,
    threshold_source  = "train_youden",
    auc               = auc_val,   accuracy    = acc_val,
    sensitivity       = sens_val,  specificity = spec_val,
    ppv               = ppv_val,   npv         = npv_val,
    auc_lo            = ci$lo[1],  acc_lo      = ci$lo[2],
    sens_lo           = ci$lo[3],  spec_lo     = ci$lo[4],
    ppv_lo            = ci$lo[5],  npv_lo      = ci$lo[6],
    auc_hi            = ci$hi[1],  acc_hi      = ci$hi[2],
    sens_hi           = ci$hi[3],  spec_hi     = ci$hi[4],
    ppv_hi            = ci$hi[5],  npv_hi      = ci$hi[6],
    type              = "classification"
  )
}

cat("CLASSIFICATION METRICS with Bootstrapped 95% CI (B=1000):\n\n")
class_results <- dplyr::bind_rows(
  eval_class(final_rf_class_reg,         train_df,          test_df,          "RF (Regular)"),
  eval_class(final_xgb_class_reg,        train_df,          test_df,          "XGB (Regular)"),
  eval_class(final_logit_class_reg,      train_df,          test_df,          "Logit (Regular)"),
  eval_class(final_rf_class_tfidf,       train_df_tfidf,    test_df_tfidf,    "RF (TF-IDF)"),
  eval_class(final_xgb_class_tfidf,      train_df_tfidf,    test_df_tfidf,    "XGB (TF-IDF)"),
  eval_class(final_logit_class_tfidf,    train_df_tfidf,    test_df_tfidf,    "Logit (TF-IDF)"),
  eval_class(final_rf_class_combined,    train_df_combined, test_df_combined, "RF (Combined)"),
  eval_class(final_xgb_class_combined,   train_df_combined, test_df_combined, "XGB (Combined)"),
  eval_class(final_logit_class_combined, train_df_combined, test_df_combined, "Logit (Combined)")
)
readr::write_csv(class_results,
                 file.path(results_dir, "classification_metrics_9models_with_CI.csv"))
cat("Saved: classification_metrics_9models_with_CI.csv\n\n")

# ============================================================
# COLLECT PREDICTED PROBABILITIES
# ============================================================
cat("Collecting predicted probabilities for plots...\n\n")
model_list <- list(
  "RF (Regular)"     = list(model = final_rf_class_reg,         train = train_df,          test = test_df),
  "XGB (Regular)"    = list(model = final_xgb_class_reg,        train = train_df,          test = test_df),
  "Logit (Regular)"  = list(model = final_logit_class_reg,      train = train_df,          test = test_df),
  "RF (TF-IDF)"      = list(model = final_rf_class_tfidf,       train = train_df_tfidf,    test = test_df_tfidf),
  "XGB (TF-IDF)"     = list(model = final_xgb_class_tfidf,      train = train_df_tfidf,    test = test_df_tfidf),
  "Logit (TF-IDF)"   = list(model = final_logit_class_tfidf,    train = train_df_tfidf,    test = test_df_tfidf),
  "RF (Combined)"    = list(model = final_rf_class_combined,    train = train_df_combined, test = test_df_combined),
  "XGB (Combined)"   = list(model = final_xgb_class_combined,   train = train_df_combined, test = test_df_combined),
  "Logit (Combined)" = list(model = final_logit_class_combined, train = train_df_combined, test = test_df_combined)
)

all_probs_test <- purrr::imap_dfr(model_list, function(entry, nm) {
  probs <- predict(entry$model, entry$test, type = "prob")$.pred_1
  tibble::tibble(
    model = nm,
    prob  = probs,
    truth = as.integer(as.character(factor(entry$test[[y_var]], levels = c(0, 1))))
  )
})

all_probs_train <- purrr::imap_dfr(model_list, function(entry, nm) {
  probs <- predict(entry$model, entry$train, type = "prob")$.pred_1
  tibble::tibble(
    model = nm,
    prob  = probs,
    truth = as.integer(as.character(factor(entry$train[[y_var]], levels = c(0, 1))))
  )
})

model_levels <- c(
  "RF (Regular)",    "XGB (Regular)",    "Logit (Regular)",
  "RF (TF-IDF)",     "XGB (TF-IDF)",     "Logit (TF-IDF)",
  "RF (Combined)",   "XGB (Combined)",   "Logit (Combined)"
)

all_probs_test$model  <- factor(all_probs_test$model,  levels = model_levels)
all_probs_train$model <- factor(all_probs_train$model, levels = model_levels)

# ============================================================
# FIX 20: BRIER SCORES + 95% BOOTSTRAP CI
# ============================================================
cat("Computing Brier scores with 95% bootstrap CI...\n\n")

brier_score <- function(prob, truth) mean((prob - truth)^2)

bootstrap_brier_ci <- function(prob, truth, R = 1000, alpha = 0.05) {
  set.seed(42)
  n       <- length(prob)
  boot_bs <- replicate(R, {
    idx <- sample(n, n, replace = TRUE)
    brier_score(prob[idx], truth[idx])
  })
  list(
    estimate = brier_score(prob, truth),
    lower    = quantile(boot_bs, alpha / 2),
    upper    = quantile(boot_bs, 1 - alpha / 2)
  )
}

compute_brier_table <- function(probs_df) {
  probs_df %>%
    dplyr::group_by(model) %>%
    dplyr::group_map(~ {
      ci <- bootstrap_brier_ci(.x$prob, .x$truth)
      tibble::tibble(
        model          = .y$model,
        brier_estimate = round(ci$estimate, 4),
        ci_lower_95    = round(ci$lower,    4),
        ci_upper_95    = round(ci$upper,    4)
      )
    }) %>%
    dplyr::bind_rows()
}

brier_stats_test  <- compute_brier_table(all_probs_test)
brier_stats_train <- compute_brier_table(all_probs_train)

readr::write_csv(brier_stats_test,
                 file.path(results_dir, "Brier_Scores_Test_95CI.csv"))
readr::write_csv(brier_stats_train,
                 file.path(results_dir, "Brier_Scores_Train_95CI.csv"))
cat("Saved Brier_Scores_Test_95CI.csv\n")
cat("Saved Brier_Scores_Train_95CI.csv\n\n")

# ============================================================
# ROC CURVE BUILDER
# ============================================================
build_roc_df <- function(prob, truth_int) {
  ord    <- order(prob, decreasing = TRUE)
  tp_vec <- cumsum(truth_int[ord])
  fp_vec <- cumsum(1 - truth_int[ord])
  n_pos  <- sum(truth_int)
  n_neg  <- length(truth_int) - n_pos
  tibble::tibble(
    fpr = c(0, fp_vec / n_neg, 1),
    tpr = c(0, tp_vec / n_pos, 1)
  )
}

# ============================================================
# ROC — 3x3 BW GRIDS
# ============================================================
cat("Generating ROC plots (3x3 grids)...\n\n")

make_roc_panel_one <- function(nm, probs_df, auc_override = NULL) {
  d <- dplyr::filter(probs_df, model == nm)
  
  if (!is.null(auc_override)) {
    auc_v <- auc_override
  } else {
    pos <- d$prob[d$truth == 1]; neg <- d$prob[d$truth == 0]
    auc_v <- mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "=="))
  }
  
  roc_df <- build_roc_df(d$prob, d$truth)
  
  ggplot2::ggplot(roc_df, ggplot2::aes(x = fpr, y = tpr)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dotdash", colour = "black", linewidth = 0.4) +
    ggplot2::geom_line(colour = "black", linewidth = 0.7) +
    ggplot2::annotate("text", x = 0.65, y = 0.08,
                      label = sprintf("AUC = %.3f", auc_v), size = 3.2) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    ggplot2::labs(title = nm, x = "1 - Specificity", y = "Sensitivity") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 8),
      axis.title       = ggplot2::element_text(size = 7),
      axis.text        = ggplot2::element_text(size = 6),
      panel.grid.minor = ggplot2::element_blank()
    )
}

roc_test_panels <- lapply(model_levels, function(nm) {
  auc_v <- class_results$auc[class_results$model == nm]
  make_roc_panel_one(nm, all_probs_test, auc_override = auc_v)
})

roc_train_panels <- lapply(model_levels, function(nm) {
  make_roc_panel_one(nm, all_probs_train)
})

roc_test_grid <- ggpubr::ggarrange(plotlist = roc_test_panels, ncol = 3, nrow = 3)
roc_test_grid <- ggpubr::annotate_figure(
  roc_test_grid,
  top = ggpubr::text_grob(
    "ROC Curves - Test Set (All 9 Models) | PHQ-2 Depression",
    face = "bold", size = 11
  )
)

roc_train_grid <- ggpubr::ggarrange(plotlist = roc_train_panels, ncol = 3, nrow = 3)
roc_train_grid <- ggpubr::annotate_figure(
  roc_train_grid,
  top = ggpubr::text_grob(
    "ROC Curves - Training Set (All 9 Models) | PHQ-2 Depression",
    face = "bold", size = 11
  )
)

ggplot2::ggsave(file.path(results_dir, "ROC_Test_9Grid_BW.png"),
                roc_test_grid,  width = 10, height = 10, dpi = 300)
ggplot2::ggsave(file.path(results_dir, "ROC_Train_9Grid_BW.png"),
                roc_train_grid, width = 10, height = 10, dpi = 300)
writeLines("done", checkpoint_roc_plots)
cat("Saved ROC_Test_9Grid_BW.png\n")
cat("Saved ROC_Train_9Grid_BW.png\n\n")

# ============================================================
# CALIBRATION — 3x3 BW GRIDS WITH BRIER ANNOTATION
# ============================================================
cat("Generating calibration plots (3x3 grids)...\n\n")

build_cal_df_adaptive <- function(prob, truth_int, n_bins = 10) {
  p_min  <- min(prob)
  p_max  <- max(prob)
  breaks <- seq(p_min, p_max, length.out = n_bins + 1)
  bin_id <- cut(prob, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  tibble::tibble(prob = prob, truth = truth_int, bin = bin_id) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      mean_pred = mean(prob),
      mean_obs  = mean(truth),
      n         = dplyr::n(),
      se        = sqrt(mean_obs * (1 - mean_obs) / n),
      .groups   = "drop"
    ) %>%
    dplyr::filter(!is.na(mean_pred))
}

make_cal_panel_one <- function(nm, probs_df, brier_df) {
  d      <- dplyr::filter(probs_df, model == nm)
  cal_df <- build_cal_df_adaptive(d$prob, d$truth)
  
  x_min <- max(0, min(cal_df$mean_pred) - 0.02)
  x_max <- min(1, max(cal_df$mean_pred) + 0.02)
  y_min <- max(0, min(cal_df$mean_obs - 1.96 * cal_df$se) - 0.02)
  y_max <- min(1, max(cal_df$mean_obs + 1.96 * cal_df$se) + 0.02)
  
  bs <- dplyr::filter(brier_df, model == nm)
  brier_label <- sprintf(
    "Brier: %.3f\n95%% CI [%.3f, %.3f]",
    bs$brier_estimate, bs$ci_lower_95, bs$ci_upper_95
  )
  
  ggplot2::ggplot(cal_df, ggplot2::aes(x = mean_pred, y = mean_obs)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "black", linewidth = 0.4) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_obs - 1.96 * se,
                   ymax = mean_obs + 1.96 * se),
      width = 0.01, linewidth = 0.3, colour = "black"
    ) +
    ggplot2::geom_line(colour  = "black", linewidth = 0.6) +
    ggplot2::geom_point(colour = "black", size = 2, shape = 16) +
    ggplot2::annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.05, vjust = 1.3,
      label = brier_label,
      size  = 2.2, colour = "black"
    ) +
    ggplot2::scale_x_continuous(limits = c(x_min, x_max)) +
    ggplot2::scale_y_continuous(limits = c(y_min, y_max)) +
    ggplot2::labs(title = nm,
                  x = "Mean Predicted Prob.",
                  y = "Observed Event Rate") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 8),
      axis.title       = ggplot2::element_text(size = 7),
      axis.text        = ggplot2::element_text(size = 6),
      panel.grid.minor = ggplot2::element_blank()
    )
}

cal_test_panels  <- lapply(model_levels, make_cal_panel_one,
                           probs_df = all_probs_test,  brier_df = brier_stats_test)
cal_train_panels <- lapply(model_levels, make_cal_panel_one,
                           probs_df = all_probs_train, brier_df = brier_stats_train)

cal_test_grid <- ggpubr::ggarrange(plotlist = cal_test_panels, ncol = 3, nrow = 3)
cal_test_grid <- ggpubr::annotate_figure(
  cal_test_grid,
  top = ggpubr::text_grob(
    "Calibration Curves - Test Set (All 9 Models) | PHQ-2 Depression\nDashed = perfect calibration  |  Error bars = +/-1.96 SE",
    face = "bold", size = 10
  )
)

cal_train_grid <- ggpubr::ggarrange(plotlist = cal_train_panels, ncol = 3, nrow = 3)
cal_train_grid <- ggpubr::annotate_figure(
  cal_train_grid,
  top = ggpubr::text_grob(
    "Calibration Curves - Training Set (All 9 Models) | PHQ-2 Depression\nDashed = perfect calibration  |  Error bars = +/-1.96 SE",
    face = "bold", size = 10
  )
)

ggplot2::ggsave(file.path(results_dir, "Calibration_Test_9Grid_BW.png"),
                cal_test_grid,  width = 10, height = 10, dpi = 300)
ggplot2::ggsave(file.path(results_dir, "Calibration_Train_9Grid_BW.png"),
                cal_train_grid, width = 10, height = 10, dpi = 300)
writeLines("done", checkpoint_cal_plots)
cat("Saved Calibration_Test_9Grid_BW.png\n")
cat("Saved Calibration_Train_9Grid_BW.png\n\n")

# ============================================================
# FINAL SUMMARY
# ============================================================
cat("============================================================\n")
cat("PIPELINE COMPLETE\n")
cat("============================================================\n")
cat("Results saved to: ", results_dir, "\n\n")
cat("Output files:\n")
cat("  classification_metrics_9models_with_CI.csv\n")
cat("  Brier_Scores_Test_95CI.csv\n")
cat("  Brier_Scores_Train_95CI.csv\n")
cat("  ROC_Test_9Grid_BW.png\n")
cat("  ROC_Train_9Grid_BW.png\n")
cat("  Calibration_Test_9Grid_BW.png  (Brier + 95% CI annotated per panel)\n")
cat("  Calibration_Train_9Grid_BW.png (Brier + 95% CI annotated per panel)\n\n")