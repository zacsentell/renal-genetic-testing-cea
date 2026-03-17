# 01c_load_vus.R
# Purpose: Load VUS probability parameters (Beta distributions) by test modality
# Author: Renal Genetics CEA Team
# Date: 2025-12-24

# Set CRAN mirror for non-interactive sessions
local({
    r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org"
    options(repos = r)
})

if (!require("readxl")) install.packages("readxl")
if (!require("dplyr")) install.packages("dplyr")
if (!require("stringr")) install.packages("stringr")
if (!require("purrr")) install.packages("purrr")

library(readxl)
library(dplyr)
library(stringr)
library(purrr)

# IO Paths
INPUT_VUS <- "data/raw/vus_curation.xlsx"
OUTPUT_VUS_RDS <- "data/params/vus_betas.rds"

# Ensure output directory exists
if (!dir.exists(dirname(OUTPUT_VUS_RDS))) dir.create(dirname(OUTPUT_VUS_RDS), recursive = TRUE)

cat("--- Starting VUS Parameter Loading ---\n")
cat("Reading", INPUT_VUS, "...\n")

# Read Data
vus_raw <- read_excel(INPUT_VUS)

# Inspect columns to handle different column naming conventions if needed
# We expect: "Test modality", "alpha", "beta"
if ("Test modality" %in% names(vus_raw)) {
    vus_raw <- vus_raw %>% rename(test_modality = `Test modality`)
} else if ("modality" %in% names(vus_raw)) {
    vus_raw <- vus_raw %>% rename(test_modality = modality)
}

cat("Available Modalities:\n")
print(unique(vus_raw$test_modality))

# Helper to extract alphas/betas
get_params <- function(pattern) {
    row <- vus_raw %>%
        filter(str_detect(test_modality, regex(pattern, ignore_case = TRUE)))

    if (nrow(row) != 1) {
        stop(paste0("Could not uniquely identify row for pattern: '", pattern, "' - Found ", nrow(row), " matches."))
    }

    list(
        shape1 = row$alpha,
        shape2 = row$beta,
        mean_sens = row$alpha / (row$alpha + row$beta) # Expected value of Beta
    )
}

# -------------------------------------------------------------------------
# Map Modalities to Parameters
# -------------------------------------------------------------------------

# 1. Panels (All Bins)
# We load all available gene count bins to allow downstream logic to pick the best fit.
# Expected bins: "11-25", "26-50", "51-100", "101-200", ">200"

panel_bins <- list(
    bin_11_25   = get_params("11-25"),
    bin_26_50   = get_params("26-50"),
    bin_51_100  = get_params("51-100"),
    bin_101_200 = get_params("101-200"),
    bin_gt_200  = get_params(">200")
)

# 2. ES and GS
# We map these directly to their specific rows
es_params <- get_params("Exome sequencing")
gs_params <- get_params("Genome sequencing")


# -------------------------------------------------------------------------
# Construct Output Object
# -------------------------------------------------------------------------

vus_params <- list(
    panels = panel_bins,
    es = es_params,
    gs = gs_params,
    metadata = list(
        timestamp = Sys.time(),
        source_file = INPUT_VUS
    )
)

print(str(vus_params))

saveRDS(vus_params, OUTPUT_VUS_RDS)
cat("\nSUCCESS: VUS parameters saved to", OUTPUT_VUS_RDS, "\n")
