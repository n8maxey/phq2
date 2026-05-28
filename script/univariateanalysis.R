# PHQ-2 DEPRESSION: UNIVARIATE STATISTICS TABLE
# Load packages
library(readr)
library(dplyr)
library(gtsummary)
library(flextable)
library(officer)
library(lubridate)
cat("✓ All packages loaded\n\n")
# PATHS
base_dir <- "C:/Users/natmaxey/OneDrive - Indiana University/Desktop/ml phq2 project/"
data_dir <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")
extract_dir <- file.path(base_dir, "extracted_csvs")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
phq2_clean_rds_path <- file.path(data_dir, "hps_cleaned_depression.rds")
phq2_clean_csv_path <- file.path(data_dir, "hps_cleaned_depression.csv")
checkpoint_table1 <- file.path(results_dir, "Table1_Univariate_Statistics.docx")
y_var <- "depression"
x_vars <- c("age", "sex", "race", "hispanic", "education", "marital", "numkids")
# ---- PHQ-2 DETECTION ----
find_adult_phq2_items <- function(df_names) {
  item1_var <- if ("INTEREST" %in% df_names) "INTEREST" else NULL
  item2_var <- if ("DOWN"     %in% df_names) "DOWN"     else NULL
  list(item1 = item1_var, item2 = item2_var)
}
# ---- PROCESS DATA ----
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
    filter(!is.na(depression)) %>%
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
      
      numkids = if (!is.na(numkids_col)) suppressWarnings(as.numeric(.data[[numkids_col]])) else NA_real_
    ) %>%
    select(any_of(c("phq2_item1","phq2_item2","phq2_score","depression","age","sex","race","hispanic","education","marital","numkids")))
  
  df
}
# ---- LOAD CSVs ----
load_hps_csvs <- function(extract_dir) {
  csvs <- list.files(extract_dir, pattern = "hps_04_.*puf\\.csv$",
                     full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
  csvs <- csvs[!grepl("repwgt", csvs, ignore.case = TRUE)]
  if (length(csvs) == 0) { cat("No CSVs found in:", extract_dir, "\n"); return(NULL) }
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
  if (length(out) == 0) { cat("No data extracted\n"); return(NULL) }
  dplyr::bind_rows(out)
}
# ---- MAIN ----
cat("Loading HPS data...\n")
hps_cleaned <- load_hps_csvs(extract_dir)

# Ensure numkids exists even if not found in any CSV
if (!"numkids" %in% names(hps_cleaned)) {
  hps_cleaned$numkids <- NA_real_
}

saveRDS(hps_cleaned, phq2_clean_rds_path)
write_csv(hps_cleaned, phq2_clean_csv_path)
# ---- ANALYTIC SAMPLE (complete-case, matches manuscript N = 605,877) ----
# Apply the same CCA logic as the pipeline:
#   filter(!if_any(all_of(c(y_var, x_vars)), is.na))
analytic_sample <- hps_cleaned %>%
  filter(!if_any(all_of(c(y_var, x_vars, "phq2_score")), is.na))
cat("Analytic sample N =", nrow(analytic_sample), "(target: 605,877)\n\n")

# ---- TABLE 1 ----
# CREATE UNIVARIATE TABLE
cat("Creating univariate statistics table...\n\n")
table1 <- analytic_sample %>%
  mutate(numkids_cat = factor(numkids)) %>%
  select(
    depression,
    age,
    sex,
    race,
    hispanic,
    education,
    marital,
    numkids_cat,
    phq2_score
  ) %>%
  tbl_summary(
    by = depression,
    type = list(
      age         ~ "continuous",
      phq2_score  ~ "continuous",
      numkids_cat ~ "categorical"
    ),
    statistic = list(
      age        ~ "{mean} ({sd})",
      phq2_score ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      age         ~ "Age (years), mean (SD)",
      sex         ~ "Sex",
      race        ~ "Race",
      hispanic    ~ "Hispanic/Latino",
      education   ~ "Education Level",
      marital     ~ "Marital Status",
      numkids_cat ~ "Number of Children",
      phq2_score  ~ "PHQ-2 Score, mean (SD)"
    )
  ) %>%
  add_overall(last = TRUE) %>%
  add_p() %>%
  modify_header(label ~ "**Characteristic**")
cat("✓ Table created\n\n")

flextable::save_as_docx(as_flex_table(table1), path = checkpoint_table1)
cat("✓ Table 1 saved successfully\n")