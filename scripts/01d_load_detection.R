# 01d_load_detection.R
# Purpose: Load analytic detection parameters (Beta distributions) by modality and variant class
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
if (!require("tidyr")) install.packages("tidyr")

library(readxl)
library(dplyr)
library(tidyr)

# ==============================================================================
# IO Paths
# ==============================================================================
INPUT_FILE <- "data/raw/variant_performance_curation.xlsx"
INPUT_SHEET <- "tab2_pooled"
OUTPUT_RDS <- "data/params/detection_betas.rds"
OUTPUT_SUMMARY_CSV <- "outputs/parameters/detection_performance_matrix.csv"

# Ensure output directory exists
if (!dir.exists(dirname(OUTPUT_RDS))) dir.create(dirname(OUTPUT_RDS), recursive = TRUE)
if (!dir.exists(dirname(OUTPUT_SUMMARY_CSV))) dir.create(dirname(OUTPUT_SUMMARY_CSV), recursive = TRUE)

cat("--- Starting Analytic Detection Parameter Loading ---\n")
cat("Reading", INPUT_FILE, "sheet:", INPUT_SHEET, "...\n")

# ==============================================================================
# 1. Load and Validate Data
# ==============================================================================
raw_data <- read_excel(INPUT_FILE, sheet = INPUT_SHEET)

# Required columns
required_cols <- c("modality_sim", "variant_class_sim", "alpha", "beta", "mean_sensitivity")
missing_cols <- setdiff(required_cols, names(raw_data))
if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

cat("Loaded", nrow(raw_data), "rows from curated detection data.\n\n")

# ==============================================================================
# 2. Harmonize Modality and Variant Class Names
# ==============================================================================
# Modality mapping to model conventions:
#   PanelExome -> applies to both Panel and ES (shared detection characteristics)
#   Genome     -> GS
#   MUC1_assay -> specialized assay for MUC1 VNTR
#   PKD1_assay -> specialized assay for PKD1 pseudogene region

# Variant class mapping to model conventions:
#   SNV_INDEL       -> SNV_Indel (standard variants)
#   CNV_SV          -> CNV_SV (copy number and structural variants)
#   PKD1_pseudogene -> Difficult_Locus (PKD1-specific)
#   MUC1_VNTR       -> Difficult_Locus (MUC1-specific)

detection_clean <- raw_data %>%
    select(modality_sim, variant_class_sim, alpha, beta, mean_sensitivity) %>%
    mutate(
        # Harmonize variant class
        variant_class_model = case_when(
            variant_class_sim == "SNV_INDEL" ~ "SNV_Indel",
            variant_class_sim == "CNV_SV" ~ "CNV_SV",
            variant_class_sim %in% c("PKD1_pseudogene", "MUC1_VNTR") ~ "Difficult_Locus",
            TRUE ~ variant_class_sim
        ),
        # Keep difficult locus subtype for conditional logic
        difficult_locus_type = case_when(
            variant_class_sim == "PKD1_pseudogene" ~ "PKD1",
            variant_class_sim == "MUC1_VNTR" ~ "MUC1",
            TRUE ~ NA_character_
        )
    )

cat("Harmonized Modality/Variant Class Mapping:\n")
print(detection_clean %>% select(modality_sim, variant_class_sim, variant_class_model, difficult_locus_type))
cat("\n")

# ==============================================================================
# 3. Build Detection Matrix Structure
# ==============================================================================
# Create a helper to extract Beta parameters as a list
get_beta_params <- function(data, modality, variant_class_orig) {
    row <- data %>%
        filter(modality_sim == modality, variant_class_sim == variant_class_orig)

    if (nrow(row) == 0) {
        # Not applicable - return zero sensitivity
        return(list(alpha = 0, beta = 1, mean_sensitivity = 0, applicable = FALSE))
    }

    if (is.na(row$alpha) || is.na(row$beta)) {
        # NA indicates cannot detect this variant class
        return(list(alpha = 0, beta = 1, mean_sensitivity = 0, applicable = FALSE))
    }

    list(
        alpha = row$alpha,
        beta = row$beta,
        mean_sensitivity = row$mean_sensitivity,
        applicable = TRUE
    )
}

# Build the nested detection matrix
# Structure: detection[[modality]][[variant_class]] = list(alpha, beta, mean_sensitivity)

# ==============================================================================
# RAW DETECTION MATRIX (intrinsic modality capabilities)
# ==============================================================================
# This captures the intrinsic detection capability of each modality alone,
# WITHOUT any bundled specialized assays.

raw_detection_matrix <- list(
    # Panel and ES share PanelExome detection characteristics (from methods.md)
    Panel = list(
        SNV_Indel = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = list(
            PKD1 = get_beta_params(detection_clean, "PanelExome", "PKD1_pseudogene"),
            MUC1 = get_beta_params(detection_clean, "PanelExome", "MUC1_VNTR")
        )
    ),
    ES = list(
        SNV_Indel = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = list(
            PKD1 = get_beta_params(detection_clean, "PanelExome", "PKD1_pseudogene"),
            MUC1 = get_beta_params(detection_clean, "PanelExome", "MUC1_VNTR")
        )
    ),
    GS = list(
        SNV_Indel = get_beta_params(detection_clean, "Genome", "SNV_INDEL"),
        CNV_SV = get_beta_params(detection_clean, "Genome", "CNV_SV"),
        Difficult_Locus = list(
            PKD1 = get_beta_params(detection_clean, "Genome", "PKD1_pseudogene"),
            MUC1 = get_beta_params(detection_clean, "Genome", "MUC1_VNTR")
        )
    ),
    # Specialized assays for difficult loci (standalone reference)
    MUC1_assay = list(
        SNV_Indel = list(alpha = 0, beta = 1, mean_sensitivity = 0, applicable = FALSE),
        CNV_SV = list(alpha = 0, beta = 1, mean_sensitivity = 0, applicable = FALSE),
        Difficult_Locus = list(
            PKD1 = list(alpha = 0, beta = 1, mean_sensitivity = 0, applicable = FALSE),
            MUC1 = get_beta_params(detection_clean, "MUC1_assay", "MUC1_VNTR")
        )
    ),
    PKD1_assay = list(
        SNV_Indel = list(alpha = 0, beta = 1, mean_sensitivity = 0, applicable = FALSE),
        CNV_SV = list(alpha = 0, beta = 1, mean_sensitivity = 0, applicable = FALSE),
        Difficult_Locus = list(
            PKD1 = get_beta_params(detection_clean, "PKD1_assay", "PKD1_pseudogene"),
            MUC1 = list(alpha = 0, beta = 1, mean_sensitivity = 0, applicable = FALSE)
        )
    )
)

# ==============================================================================
# PHENOTYPE-SPECIFIC PANEL DETECTION MATRIX (for CEA model use)
# ==============================================================================
# MODELING ASSUMPTION (per user specification):
#
# For phenotype-specific panel testing modalities (but NOT the comprehensive
# panel for CKDu), we assume the inclusion of a targeted PKD1 or MUC1 assay:
#
#   - CYSTIC panel: includes PKD1 assay (but NOT MUC1)
#   - TUBULOINTERSTITIAL (ADTKD) panel: includes MUC1 assay (but NOT PKD1)
#   - GLOMERULAR panel: no bundled difficult locus assays
#   - TUBULOPATHIES panel: no bundled difficult locus assays
#   - COMPREHENSIVE panel (CKDu): no bundled difficult locus assays
#
# ES and GS strategies do NOT include these specialized assays regardless of
# phenotype - they rely only on their intrinsic detection capabilities.

# Helper: create a standard difficult locus entry with optional assay override
make_difficult_locus <- function(
    pkd1_use_assay = FALSE,
    muc1_use_assay = FALSE,
    detection_data = detection_clean) {
    list(
        PKD1 = if (pkd1_use_assay) {
            get_beta_params(detection_data, "PKD1_assay", "PKD1_pseudogene")
        } else {
            get_beta_params(detection_data, "PanelExome", "PKD1_pseudogene")
        },
        MUC1 = if (muc1_use_assay) {
            get_beta_params(detection_data, "MUC1_assay", "MUC1_VNTR")
        } else {
            get_beta_params(detection_data, "PanelExome", "MUC1_VNTR")
        }
    )
}

# Phenotype-specific panel detection parameters
panel_by_phenotype <- list(
    # Cystic panel: includes PKD1 assay (for ADPKD detection)
    cystic = list(
        SNV_Indel = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = make_difficult_locus(pkd1_use_assay = TRUE, muc1_use_assay = FALSE),
        bundled_assays = c("PKD1")
    ),
    # Tubulointerstitial (ADTKD) panel: includes MUC1 assay
    tubulointerstitial = list(
        SNV_Indel = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = make_difficult_locus(pkd1_use_assay = FALSE, muc1_use_assay = TRUE),
        bundled_assays = c("MUC1")
    ),
    # Glomerular panel: no bundled difficult locus assays
    glomerular = list(
        SNV_Indel = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = make_difficult_locus(pkd1_use_assay = FALSE, muc1_use_assay = FALSE),
        bundled_assays = character(0)
    ),
    # Tubulopathies panel: no bundled difficult locus assays
    tubulopathies = list(
        SNV_Indel = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = make_difficult_locus(pkd1_use_assay = FALSE, muc1_use_assay = FALSE),
        bundled_assays = character(0)
    ),
    # Comprehensive panel (CKDu): no bundled difficult locus assays
    CKDu = list(
        SNV_Indel = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = make_difficult_locus(pkd1_use_assay = FALSE, muc1_use_assay = FALSE),
        bundled_assays = character(0)
    )
)

# ==============================================================================
# ES AUGMENTED: phenotype-keyed detection with bundled difficult-locus assays
# ==============================================================================
# PKD1 assay applied universally (all phenotypes) — PKD1 pseudogene variants can
# underlie atypical presentations across phenotype categories and a co-ordered
# locus-specific assay would not be clinically restricted by phenotype.
# MUC1 assay applied to tubulointerstitial only — ADTKD-MUC1 has a disease-
# specific indication not present in other phenotype categories.
# SNV/Indel and CNV/SV detection identical to base-case ES (PanelExome).
# Bundled assay costs assumed absorbed within ES unit price (scenario assumption).
es_augmented_by_phenotype <- list(
    cystic = list(
        SNV_Indel       = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV          = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = make_difficult_locus(pkd1_use_assay = TRUE, muc1_use_assay = FALSE),
        bundled_assays  = c("PKD1")
    ),
    tubulointerstitial = list(
        SNV_Indel       = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV          = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = make_difficult_locus(pkd1_use_assay = TRUE, muc1_use_assay = TRUE),
        bundled_assays  = c("PKD1", "MUC1")
    ),
    glomerular = list(
        SNV_Indel       = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV          = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = make_difficult_locus(pkd1_use_assay = TRUE, muc1_use_assay = FALSE),
        bundled_assays  = c("PKD1")
    ),
    tubulopathies = list(
        SNV_Indel       = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV          = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = make_difficult_locus(pkd1_use_assay = TRUE, muc1_use_assay = FALSE),
        bundled_assays  = c("PKD1")
    ),
    CKDu = list(
        SNV_Indel       = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV          = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = make_difficult_locus(pkd1_use_assay = TRUE, muc1_use_assay = FALSE),
        bundled_assays  = c("PKD1")
    )
)

# ES detection (same for all phenotypes, no bundled assays)
es_detection <- list(
    SNV_Indel = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
    CNV_SV = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
    Difficult_Locus = make_difficult_locus(pkd1_use_assay = FALSE, muc1_use_assay = FALSE),
    bundled_assays = character(0)
)

# GS detection (same for all phenotypes, no bundled assays, uses Genome intrinsic)
gs_detection <- list(
    SNV_Indel = get_beta_params(detection_clean, "Genome", "SNV_INDEL"),
    CNV_SV = get_beta_params(detection_clean, "Genome", "CNV_SV"),
    Difficult_Locus = list(
        PKD1 = get_beta_params(detection_clean, "Genome", "PKD1_pseudogene"),
        MUC1 = get_beta_params(detection_clean, "Genome", "MUC1_VNTR")
    ),
    bundled_assays = character(0)
)

# Construct the full strategy detection matrix
strategy_detection_matrix <- list(
    # Phenotype-specific panels
    Panel = panel_by_phenotype,
    # ES (phenotype-agnostic, no bundled assays)
    ES = es_detection,
    # GS (phenotype-agnostic, no bundled assays)
    GS = gs_detection,
    # ES Augmented: phenotype-keyed with PKD1 assay universal, MUC1 assay for tubulointerstitial
    ES_Augmented = es_augmented_by_phenotype
)

# For backward compatibility, set detection_matrix to the strategy version
detection_matrix <- strategy_detection_matrix

# ==============================================================================
# 4. Summary Table for Verification
# ==============================================================================
cat("=== Detection Matrix Summary ===\n\n")

# Build a flat summary table for the phenotype-specific structure
summary_rows <- list()

# Helper to add rows for a detection entry
add_detection_rows <- function(modality_name, det_entry, phenotype = NA) {
    rows <- list()
    display_name <- if (is.na(phenotype)) modality_name else paste0("Panel (", phenotype, ")")
    bundled <- if (!is.null(det_entry$bundled_assays)) paste(det_entry$bundled_assays, collapse = ",") else ""

    for (variant_class in c("SNV_Indel", "CNV_SV")) {
        params <- det_entry[[variant_class]]
        rows[[length(rows) + 1]] <- data.frame(
            Modality = display_name,
            Variant_Class = variant_class,
            Alpha = params$alpha,
            Beta = params$beta,
            Mean_Sensitivity = sprintf("%.3f", params$mean_sensitivity),
            Applicable = params$applicable,
            Bundled_Assays = bundled,
            stringsAsFactors = FALSE
        )
    }
    # Difficult locus subtypes
    for (locus in c("PKD1", "MUC1")) {
        params <- det_entry[["Difficult_Locus"]][[locus]]
        rows[[length(rows) + 1]] <- data.frame(
            Modality = display_name,
            Variant_Class = paste0("Difficult_Locus (", locus, ")"),
            Alpha = params$alpha,
            Beta = params$beta,
            Mean_Sensitivity = sprintf("%.3f", params$mean_sensitivity),
            Applicable = params$applicable,
            Bundled_Assays = bundled,
            stringsAsFactors = FALSE
        )
    }
    rows
}

# Panel (phenotype-specific)
for (pheno in names(panel_by_phenotype)) {
    summary_rows <- c(summary_rows, add_detection_rows("Panel", panel_by_phenotype[[pheno]], pheno))
}

# ES
summary_rows <- c(summary_rows, add_detection_rows("ES", es_detection))

# GS
summary_rows <- c(summary_rows, add_detection_rows("GS", gs_detection))

# ES_Augmented (phenotype-specific)
for (pheno in names(es_augmented_by_phenotype)) {
    summary_rows <- c(summary_rows, add_detection_rows("ES_Augmented", es_augmented_by_phenotype[[pheno]], pheno))
}

summary_df <- bind_rows(summary_rows)
print(as.data.frame(summary_df))

summary_export <- summary_df %>%
    transmute(
        modality = Modality,
        variant_class = Variant_Class,
        alpha = Alpha,
        beta = Beta,
        mean_sensitivity = as.numeric(Mean_Sensitivity),
        applicable = Applicable,
        bundled_assays = Bundled_Assays
    )

write.csv(summary_export, OUTPUT_SUMMARY_CSV, row.names = FALSE)

# ==============================================================================
# 5. Validation Checks
# ==============================================================================
cat("\n=== Validation Checks ===\n")

# Check 1: All alpha/beta should be >= 0
check_positive <- function(lst, path = "") {
    for (name in names(lst)) {
        current_path <- paste0(path, "/", name)
        if (is.list(lst[[name]]) && !is.null(lst[[name]]$alpha)) {
            if (lst[[name]]$alpha < 0 || lst[[name]]$beta < 0) {
                warning("Negative alpha/beta at: ", current_path)
            }
        } else if (is.list(lst[[name]])) {
            check_positive(lst[[name]], current_path)
        }
    }
}
check_positive(detection_matrix)
cat("CHECK 1: All alpha/beta values are non-negative.\n")

# Check 2: Mean sensitivity should be in [0,1]
check_mean <- function(lst, path = "") {
    for (name in names(lst)) {
        current_path <- paste0(path, "/", name)
        if (is.list(lst[[name]]) && !is.null(lst[[name]]$mean_sensitivity)) {
            m <- lst[[name]]$mean_sensitivity
            if (!is.na(m) && (m < 0 || m > 1)) {
                warning("Mean sensitivity out of range at: ", current_path, " (", m, ")")
            }
        } else if (is.list(lst[[name]])) {
            check_mean(lst[[name]], current_path)
        }
    }
}
check_mean(detection_matrix)
cat("CHECK 2: All mean sensitivities are in [0, 1].\n")

# ==============================================================================
# 6. Save Output
# ==============================================================================
output_obj <- list(
    strategy_detection_matrix = strategy_detection_matrix,
    raw_data = detection_clean,
    metadata = list(
        timestamp = Sys.time(),
        source_file = INPUT_FILE,
        source_sheet = INPUT_SHEET,
        modeling_note = paste(
            "Phenotype-specific panel testing includes targeted assays as follows:",
            "CYSTIC panel includes PKD1 assay (but NOT MUC1).",
            "TUBULOINTERSTITIAL (ADTKD) panel includes MUC1 assay (but NOT PKD1).",
            "GLOMERULAR, TUBULOPATHIES, and COMPREHENSIVE (CKDu) panels have no bundled assays.",
            "ES and GS strategies use only intrinsic detection capabilities (no bundled assays).",
            "ES_Augmented (scenario): PKD1 assay applied universally (all phenotypes);",
            "MUC1 assay applied to tubulointerstitial only. Assay costs absorbed in ES unit price."
        )
    )
)

saveRDS(output_obj, OUTPUT_RDS)

cat("\n=== SUCCESS ===\n")
cat("Detection parameters saved to:", OUTPUT_RDS, "\n")
cat("Detection performance matrix exported to:", OUTPUT_SUMMARY_CSV, "\n")
