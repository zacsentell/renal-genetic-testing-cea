# 01_import_validate.R
# Purpose: Import and strict validation of Excel data parameters
# Author: Renal Genetics CEA Team
# Date: 2025-12-10

# Set CRAN mirror for non-interactive sessions
local({
    r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org"
    options(repos = r)
})

if (!require("readxl")) install.packages("readxl")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("purrr")) install.packages("purrr")
if (!require("janitor")) install.packages("janitor")

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
if (!require("assertthat")) install.packages("assertthat")
library(assertthat)

# IO Paths
INPUT_COHORT <- "data/raw/cohort_curation.xlsx"
INPUT_VARIANT <- "data/raw/variant_curation.xlsx"
OUTPUT_RDS <- "data/intermediate/01_imported_data.rds"

# Ensure output directory exists
if (!dir.exists(dirname(OUTPUT_RDS))) dir.create(dirname(OUTPUT_RDS), recursive = TRUE)

cat("--- Starting Import & Validation ---\n")

# 1. Import Cohort Data
cat("Reading", INPUT_COHORT, "...\n")
# Assuming sheet 1 is the main data
cohort_raw <- read_excel(INPUT_COHORT) %>%
    clean_names() %>%
    rename(
        phenotype = phenotype_category,
        n_total = n_probands_in_category,
        n_monogenic = n_solved_in_category
    )

cat("Cohort Columns:", paste(colnames(cohort_raw), collapse = ", "), "\n")

# 2. Validation: Cohort
# Expecting: phenotype, count (or similar).
# Methods demands: Phenotype categories, counts, and proportions.
# We need to standardize column names if they vary.
# Let's inspect what we have and fail if critical columns are missing.

required_cols_cohort <- c("phenotype", "n_total", "n_monogenic") # Minimal expectation
missing_cols <- setdiff(required_cols_cohort, colnames(cohort_raw))
if (length(missing_cols) > 0) {
    stop("Cohort data missing columns: ", paste(missing_cols, collapse = ", "))
}

# Validate Types and Logic
# if (!exists("assert_that")) library(assertthat) # Ensure loaded if installed above (but I removed install logic for simplicity, assuming it's there or base is enough)
assert_that(
    is.numeric(cohort_raw$n_total),
    all(cohort_raw$n_total >= 0),
    all(cohort_raw$n_monogenic <= cohort_raw$n_total, na.rm = TRUE)
)

# 3. Import Variant Data
cat("Reading", INPUT_VARIANT, "...\n")
variant_raw <- read_excel(INPUT_VARIANT) %>%
    clean_names() %>%
    rename(
        variant_class = modeled_variant_class
    )

cat("Variant Columns:", paste(colnames(variant_raw), collapse = ", "), "\n")

# 4. Validation: Variant
# Methods: PMID, patient ID, phenotype, gene symbol...
required_cols_variant <- c("gene", "phenotype", "variant_class") # Minimal assumption
missing_cols_var <- setdiff(required_cols_variant, colnames(variant_raw))

# Strict validation: stop if required columns are missing
if (length(missing_cols_var) > 0) {
    stop("Variant data missing required columns: ", paste(missing_cols_var, collapse = ", "))
}

# 5. Centralized Phenotype Standardization & Filtering
# Standardize names (tubulopathy -> tubulopathies) and filter to allowed sets
cat("Applying phenotype standardization and filtering (Target: 5 allowed groups)...\n")

allowed_groups <- c("cystic", "glomerular", "tubulointerstitial", "tubulopathies", "CKDu")

# Standardize
# Maps:
# - tubulopathy -> tubulopathies (Standardization)
# - nephrolithiasis -> tubulopathies (User Request)
# - other_mixed -> CKDu (User Request)

cohort_raw <- cohort_raw %>%
    mutate(phenotype = case_when(
        phenotype == "tubulopathy" ~ "tubulopathies",
        phenotype == "nephrolithiasis" ~ "tubulopathies",
        phenotype == "other_mixed" ~ "CKDu",
        TRUE ~ phenotype
    ))

variant_raw <- variant_raw %>%
    mutate(phenotype = case_when(
        phenotype == "tubulopathy" ~ "tubulopathies",
        phenotype == "nephrolithiasis" ~ "tubulopathies",
        phenotype == "other_mixed" ~ "CKDu",
        TRUE ~ phenotype
    )) %>%
    # Standardize Variant Class (User Request: Unify CNV/SV under "CNV")
    mutate(variant_class = case_when(
        variant_class %in% c("Exon_CNV", "Structural_SV") ~ "CNV",
        grepl("CNV|Structural", variant_class, ignore.case = TRUE) ~ "CNV",
        TRUE ~ variant_class
    ))

# Filter
cohort_clean <- cohort_raw %>% filter(phenotype %in% allowed_groups)
variant_clean <- variant_raw %>% filter(phenotype %in% allowed_groups)

# Report drops
dropped_cohort_rows <- nrow(cohort_raw) - nrow(cohort_clean)
if (dropped_cohort_rows > 0) {
    dropped_types <- setdiff(unique(cohort_raw$phenotype), allowed_groups)
    cat(sprintf(
        "  - Dropped %d cohort rows. Excluded phenotypes: %s\n",
        dropped_cohort_rows, paste(dropped_types, collapse = ", ")
    ))
}

dropped_variant_rows <- nrow(variant_raw) - nrow(variant_clean)
if (dropped_variant_rows > 0) {
    dropped_types <- setdiff(unique(variant_raw$phenotype), allowed_groups)
    cat(sprintf(
        "  - Dropped %d variant rows. Excluded phenotypes: %s\n",
        dropped_variant_rows, paste(dropped_types, collapse = ", ")
    ))
}

# Overwrite with clean data
cohort_raw <- cohort_clean
variant_raw <- variant_clean

# Final check
if (nrow(cohort_raw) == 0) warning("WARNING: Cohort data is empty after filtering!")


# 6. Export
output_list <- list(
    cohort = cohort_raw,
    variant = variant_raw,
    metadata = list(
        timestamp = Sys.time(),
        source_files = c(INPUT_COHORT, INPUT_VARIANT),
        phenotype_remapping = "complement/aHUS merged into other_mixed"
    )
)

saveRDS(output_list, OUTPUT_RDS)
cat("SUCCESS: Data imported and saved to", OUTPUT_RDS, "\n")

# 7. Audit Exports (For Transparency)
AUDIT_DIR <- "data/audit"
if (!dir.exists(AUDIT_DIR)) dir.create(AUDIT_DIR, recursive = TRUE)

cat("Exporting audit CSVs to", AUDIT_DIR, "...\n")
if (!require("readr")) install.packages("readr")
library(readr)

write_csv(cohort_raw, file.path(AUDIT_DIR, "01_cohort_imported.csv"))
write_csv(variant_raw, file.path(AUDIT_DIR, "01_variant_imported.csv"))
cat("  - Saved: 01_cohort_imported.csv\n")
cat("  - Saved: 01_variant_imported.csv\n")
