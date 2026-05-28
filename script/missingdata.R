# =============================================================================
# MISSING DATA EXCLUSIONS — PHQ-2 ANALYTIC SAMPLE
# Household Pulse Survey (HPS) 2021–2024
#
# Analytic sample logic follows Maxey et al. (manuscript_v8):
#   - Start: all eligible adults across waves
#   - Step 2: exclude rows missing PHQ-2 outcome variables (INTEREST, DOWN)
#   - Step 3: exclude rows missing any predictor/confounder variable
#   - Final target: ~605,877 (complete-case analytic sample)
#
# NOTE: Missingness is calculated on the FULLY COMBINED dataset, not per-file.
# =============================================================================

# ---- PACKAGES ----------------------------------------------------------------
library(readr)
library(dplyr)
cat("✓ Packages loaded\n\n")

# ---- PATHS (unchanged) -------------------------------------------------------
base_dir    <- "C:/Users/natmaxey/OneDrive - Indiana University/Desktop/ml phq2 project/"
results_dir <- file.path(base_dir, "results")
extract_dir <- file.path(base_dir, "extracted_csvs")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# STEP 1 — IMPORT AND PROCESS DATA
# Matches pipeline: read each CSV with tryCatch, process per-file, then combine.
# No explicit na_if pass — missing values surface through as.numeric() coercion
# exactly as in the pipeline's process_hps_data().
# =============================================================================

# PHQ-2 item detection — identical to pipeline
find_adult_phq2_items <- function(df_names) {
  item1_var <- if ("INTEREST" %in% df_names) "INTEREST" else NULL
  item2_var <- if ("DOWN"     %in% df_names) "DOWN"     else NULL
  list(item1 = item1_var, item2 = item2_var)
}

# Per-file processing:
#   - coerce all cols to character
#   - recode -99/-88 on raw PHQ-2 items to NA_real_ before any clamping
#   - derive phq2_item1, phq2_item2, phq2_score, depression
#   - derive age, sex, race, hispanic, education, marital, numkids
# filter(!is.na(depression)) is NOT applied here — rows with missing
# PHQ-2 are kept so Step 2 can count and exclude them.
process_hps_data <- function(df_raw) {
  phq <- find_adult_phq2_items(names(df_raw))
  if (is.null(phq$item1) || is.null(phq$item2)) return(NULL)
  
  age_col     <- names(df_raw)[names(df_raw) %in% c("TBIRTH_YEAR","TAGE1","TAGE")][1]
  sex_col     <- names(df_raw)[names(df_raw) %in% c("EGENID_BIRTH","ESEX1","EGENDER")][1]
  educ_col    <- names(df_raw)[names(df_raw) %in% c("EEDUC","REDUC1")][1]
  marital_col <- names(df_raw)[names(df_raw) %in% c("MS","MARITAL1")][1]
  numkids_col <- names(df_raw)[names(df_raw) %in% c("THHLD_NUMKID")][1]
  
  if (is.na(sex_col)) return(NULL)
  
  df <- df_raw %>% mutate(across(everything(), as.character))
  
  
  df <- df %>%
    mutate(
      phq2_item1 = suppressWarnings(as.numeric(.data[[phq$item1]])) - 1,
      phq2_item2 = suppressWarnings(as.numeric(.data[[phq$item2]])) - 1,
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
      
      numkids = if (!is.na(numkids_col))
        suppressWarnings(as.numeric(.data[[numkids_col]])) else NA_real_
    )
  
  df
}

# Load and process all CSVs — mirrors pipeline load_hps_csvs()
cat("Loading HPS data...\n\n")
csvs <- list.files(extract_dir, pattern = "hps_04_.*puf\\.csv$",
                   full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
csvs <- csvs[!grepl("repwgt", csvs, ignore.case = TRUE)]
cat("Found", length(csvs), "CSV(s)\n\n")

out <- list()
for (i in seq_along(csvs)) {
  f <- csvs[i]
  message("[", i, "/", length(csvs), "] ", basename(f))
  df_raw <- tryCatch(
    readr::read_csv(f, col_types = readr::cols(.default = "c"), show_col_types = FALSE),
    error = function(e) { message("  Read error"); NULL }
  )
  if (is.null(df_raw) || nrow(df_raw) == 0) { message("  Empty"); next }
  clean <- tryCatch(process_hps_data(df_raw),
                    error = function(e) { message("  Process error"); NULL })
  if (!is.null(clean) && nrow(clean) > 0) out[[length(out) + 1]] <- clean
}
combined_data <- dplyr::bind_rows(out)
n_original <- nrow(combined_data)
cat("\nCombined dataset:", format(n_original, big.mark = ","), "rows\n\n")

# =============================================================================
# RAW MISSINGNESS RATES (before any exclusions)
# Denominator = n_original (610,755) — used for reporting in Results
# =============================================================================
raw_vars <- c("phq2_item1", "phq2_item2", "marital",
              "age", "sex", "race", "hispanic", "education", "numkids")
raw_vars_present <- intersect(raw_vars, names(combined_data))
cat("--- Raw missingness rates (% of", n_original, "eligible respondents) ---\n")
for (v in raw_vars_present) {
  n_miss   <- sum(is.na(combined_data[[v]]))
  pct_miss <- round(n_miss / n_original * 100, 2)
  cat(sprintf("  %-15s  %s missing  (%.2f%%)\n",
              v, format(n_miss, big.mark = ","), pct_miss))
}
cat("\n")

# =============================================================================
# STEP 2 — OUTCOME MISSINGNESS
# -99 and -88 were recoded to NA before clamping, so depression is NA for
# any row where INTEREST or DOWN was missing/not reported.
# =============================================================================

# Excluding those with missing outcomes
analytic_v1 <- combined_data %>% filter(!is.na(depression))
exclude1 <- n_original - nrow(analytic_v1)
cat("INTEREST and DOWN:", exclude1, "missing,",
    round((exclude1 / n_original) * 100, 1), "%\n")

# =============================================================================
# STEP 3 — PREDICTOR / CONFOUNDER MISSINGNESS
# Matches pipeline ml_df_cca logic:
#   filter(!if_any(all_of(x_vars), is.na))
# where x_vars = c("age","sex","race","hispanic","education","marital","numkids")
# depression is already complete after Step 2, so it is excluded here.
# =============================================================================
x_vars <- c("age", "sex", "race", "hispanic", "education", "marital", "numkids")

# Excluding those with missing predictors/confounders
analytic_v2 <- analytic_v1 %>%
  filter(!if_any(all_of(x_vars), is.na))
exclude2 <- nrow(analytic_v1) - nrow(analytic_v2)
cat("Predictors/confounders:", exclude2, "missing,",
    round((exclude2 / n_original) * 100, 1), "%\n")

cat("Final analytic sample:", nrow(analytic_v2), "\n")

# Save the final analytic sample
analytic_v2_path <- file.path(results_dir, "analytic_sample_complete_cases.csv")
write_csv(analytic_v2, analytic_v2_path)
cat("Saved to:", analytic_v2_path, "\n")