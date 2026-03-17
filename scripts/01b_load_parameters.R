# 01b_load_parameters.R
# Purpose: Load genetic test costs and gene panels into a structured parameter object
# Author: Renal Genetics CEA Team
# Date: 2025-12-24

# Set CRAN mirror for non-interactive sessions
local({
    r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org"
    options(repos = r)
})

if (!require("readr")) install.packages("readr")
if (!require("dplyr")) install.packages("dplyr")
if (!require("purrr")) install.packages("purrr")
if (!requireNamespace("config", quietly = TRUE)) {
    stop("Package 'config' is required. Install it with install.packages('config').")
}

library(readr)
library(dplyr)
library(purrr)

# IO Paths
INPUT_COSTS <- "data/raw/genetic_test_costs.csv"
INPUT_PANELS_DIR <- "data/raw/panels"
OUTPUT_RDS <- "data/intermediate/01b_parameters.rds"

# Ensure output directory exists
if (!dir.exists(dirname(OUTPUT_RDS))) dir.create(dirname(OUTPUT_RDS), recursive = TRUE)

cat("--- Starting Parameter Loading ---\n")

# -------------------------------------------------------------------------
# 1. Load Costs
# -------------------------------------------------------------------------
cat("Loading costs from", INPUT_COSTS, "...\n")

raw_costs <- read_csv(INPUT_COSTS, show_col_types = FALSE)

required_cost_cols <- c(
    "cost_item_id", "cost_item_label", "applied_when",
    "mean_cad", "min_cad", "max_cad", "source"
)
missing_cost_cols <- setdiff(required_cost_cols, names(raw_costs))
if (length(missing_cost_cols) > 0) {
    stop("Cost table missing required columns: ", paste(missing_cost_cols, collapse = ", "))
}

required_ids <- c(
    "panel_multigene",
    "clinical_exome",
    "clinical_genome",
    "targeted_single_variant",
    "consult_pretest",
    "consult_posttest"
)

if (anyDuplicated(raw_costs$cost_item_id) > 0) {
    dupes <- unique(raw_costs$cost_item_id[duplicated(raw_costs$cost_item_id)])
    stop("Duplicate cost_item_id values found: ", paste(dupes, collapse = ", "))
}

missing_ids <- setdiff(required_ids, raw_costs$cost_item_id)
unexpected_ids <- setdiff(raw_costs$cost_item_id, required_ids)
if (length(missing_ids) > 0) {
    stop("Cost table missing required cost_item_id values: ", paste(missing_ids, collapse = ", "))
}
if (length(unexpected_ids) > 0) {
    stop("Cost table has unexpected cost_item_id values: ", paste(unexpected_ids, collapse = ", "))
}

num_cols <- c("mean_cad", "min_cad", "max_cad")
non_numeric <- num_cols[sapply(num_cols, function(col) !is.numeric(raw_costs[[col]]))]
if (length(non_numeric) > 0) {
    stop("Cost table has non-numeric columns where numeric required: ", paste(non_numeric, collapse = ", "))
}

if (anyNA(raw_costs[num_cols])) {
    stop("Cost table contains NA in numeric cost columns: ", paste(num_cols, collapse = ", "))
}

bad_range <- raw_costs %>%
    filter(!(min_cad <= mean_cad & mean_cad <= max_cad))
if (nrow(bad_range) > 0) {
    stop(
        "Invalid min/mean/max ordering for cost_item_id values: ",
        paste(bad_range$cost_item_id, collapse = ", ")
    )
}

find_cost_by_id <- function(item_id, data) {
    row <- data %>% filter(cost_item_id == item_id)
    if (nrow(row) != 1) {
        stop("Could not uniquely identify cost for cost_item_id: ", item_id)
    }
    as.numeric(row$mean_cad[[1]])
}

costs_list <- list(
    tests = list(
        panel    = find_cost_by_id("panel_multigene", raw_costs),
        es       = find_cost_by_id("clinical_exome", raw_costs),
        gs       = find_cost_by_id("clinical_genome", raw_costs),
        familial = find_cost_by_id("targeted_single_variant", raw_costs)
    ),
    consultations = list(
        pretest  = find_cost_by_id("consult_pretest", raw_costs),
        posttest = find_cost_by_id("consult_posttest", raw_costs)
    )
)

print(costs_list)

# -------------------------------------------------------------------------
# Cascade Parameters (loaded from config.yml, saved separately)
# These are epidemiological count parameters, NOT unit costs.
# They must not be co-located with costs_list to prevent accidental
# Gamma-sampling in the probabilistic sensitivity analysis.
# -------------------------------------------------------------------------
cfg <- config::get()

cascade_base <- cfg$cascade$eligible_relatives_base
cascade_min <- cfg$cascade$eligible_relatives_min
cascade_max <- cfg$cascade$eligible_relatives_max

if (!is.numeric(cascade_base) || cascade_base < 0) {
    stop("cascade.eligible_relatives_base must be a non-negative number in config.yml")
}
if (!is.numeric(cascade_min) || cascade_min < 0) {
    stop("cascade.eligible_relatives_min must be a non-negative number in config.yml")
}
if (!is.numeric(cascade_max) || cascade_max < cascade_base) {
    stop("cascade.eligible_relatives_max must be >= eligible_relatives_base in config.yml")
}

cascade_params <- list(
    eligible_relatives_base = cascade_base,
    eligible_relatives_min  = cascade_min,
    eligible_relatives_max  = cascade_max
)

cat(sprintf(
    "\nCascade parameters (from config.yml): base = %g, range = %g–%g relatives\n",
    cascade_base, cascade_min, cascade_max
))

# -------------------------------------------------------------------------
# 2. Load and Process Panels
# -------------------------------------------------------------------------
cat("\nLoading panels from", INPUT_PANELS_DIR, "...\n")

# Define file mappings
# Mappings:
# Cystic: cystic_kidney_disease_v4.3.tsv
# Glomerular: glomerulopathy_v5.2.tsv
# Tubulointerstitial: tubulointerstitial_kidney_disease_v3.6.tsv
# Tubulopathies: tubulopathies_v5.6.tsv AND nephrocalcinosis_5.6.tsv
# CKDu: unexplained_kf_v1.120.tsv

# Helper to read a TSV and extract gene symbols
read_panel_genes <- function(filename) {
    path <- file.path(INPUT_PANELS_DIR, filename)
    if (!file.exists(path)) stop(paste("Panel file not found:", path))

    df <- read_tsv(path, show_col_types = FALSE)

    # Gene column is usually "Gene Symbol" or "Entity Name" depending on PanelApp version
    # Based on investigation, it's "Entity Name" or "Gene Symbol".
    # Let's check for "Gene Symbol" first, then "Entity Name".

    if ("Gene Symbol" %in% names(df)) {
        genes <- df$`Gene Symbol`
    } else if ("Entity Name" %in% names(df)) {
        genes <- df$`Entity Name`
    } else {
        stop(paste("Could not find gene column in", filename))
    }

    # Filter out empty or NA
    unique(na.omit(genes))
}

# 2.1 Load specific phenotype panels
panel_cystic <- read_panel_genes("cystic_kidney_disease_v4.3.tsv")
panel_glomerular <- read_panel_genes("glomerulopathy_v5.2.tsv")
panel_tubulointerstitial <- read_panel_genes("tubulointerstitial_kidney_disease_v3.6.tsv")
panel_ckdu <- read_panel_genes("unexplained_kf_v1.120.tsv")

# Tubulopathies (Union of Tubulopathies + Nephrocalcinosis)
p_tub_main <- read_panel_genes("tubulopathies_v5.6.tsv")
p_nephro <- read_panel_genes("nephrocalcinosis_5.6.tsv")
panel_tubulopathies <- unique(c(p_tub_main, p_nephro))

# 2.2 Construct Comprehensive Panel (Union of ALL files in directory)
all_panel_files <- list.files(INPUT_PANELS_DIR, pattern = "\\.tsv$", full.names = FALSE)
cat("Merging the following", length(all_panel_files), "panel files for Comprehensive list:\n")
print(all_panel_files)

all_genes_list <- lapply(all_panel_files, read_panel_genes)
panel_comprehensive <- unique(unlist(all_genes_list))

cat("Comprehensive Panel Size:", length(panel_comprehensive), "unique genes\n")

# Store panels in list
# NOTE: CKDu uses the comprehensive panel (union of all) since patients with
# unknown etiology are tested with the broadest panel to maximize yield.
panels_list <- list(
    cystic = panel_cystic,
    glomerular = panel_glomerular,
    tubulointerstitial = panel_tubulointerstitial,
    tubulopathies = panel_tubulopathies,
    CKDu = panel_comprehensive, # Use comprehensive panel for unknown etiology
    comprehensive = panel_comprehensive
)

# Print sizes
cat("\nPanel Sizes (Gene Counts):\n")
print(map_int(panels_list, length))

# -------------------------------------------------------------------------
# 3. Final Object Construction and Export
# -------------------------------------------------------------------------

parameters <- list(
    panels = panels_list,
    costs = costs_list,
    metadata = list(
        timestamp = Sys.time(),
        source_files = list(
            costs = INPUT_COSTS,
            panels_dir = INPUT_PANELS_DIR
        )
    )
)

saveRDS(parameters, OUTPUT_RDS)
cat("\nSUCCESS: Parameters saved to", OUTPUT_RDS, "\n")

# -------------------------------------------------------------------------
# 4. Export to data/params/ for Simulation Engine (Step 4)
# -------------------------------------------------------------------------
PARAMS_DIR <- "data/params"
if (!dir.exists(PARAMS_DIR)) dir.create(PARAMS_DIR, recursive = TRUE)

saveRDS(costs_list, file.path(PARAMS_DIR, "costs.rds"))
saveRDS(panels_list, file.path(PARAMS_DIR, "panels.rds"))
saveRDS(cascade_params, file.path(PARAMS_DIR, "cascade_params.rds"))

cat("Exported to data/params/:\n")
cat("  - costs.rds\n")
cat("  - panels.rds\n")
cat("  - cascade_params.rds\n")
