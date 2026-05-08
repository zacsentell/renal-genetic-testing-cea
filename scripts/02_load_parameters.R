# scripts/02_load_parameters.R
# Purpose: Load costs, panels, uptake, VUS, and analytic detection parameters into data/params/.
# Author: Zachary Sentell

library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# ==============================================================================
# IO Paths
# ==============================================================================
INPUT_COSTS         <- "data/raw/genetic_test_costs.csv"
INPUT_UPTAKE        <- "data/raw/uptake_parameters.csv"
INPUT_PANELS_DIR    <- "data/raw/panels"
INPUT_VUS           <- "data/raw/vus_curation.xlsx"
INPUT_DETECTION     <- "data/raw/variant_performance_curation.xlsx"
INPUT_DETECTION_TAB <- "tab2_pooled"

OUTPUT_INTERMEDIATE <- "data/intermediate/parameters.rds"
PARAMS_DIR          <- "data/params"
OUTPUT_DETECTION_CSV <- "outputs/results/parameters/detection_performance_matrix.csv"

if (!dir.exists(dirname(OUTPUT_INTERMEDIATE))) dir.create(dirname(OUTPUT_INTERMEDIATE), recursive = TRUE)
if (!dir.exists(PARAMS_DIR)) dir.create(PARAMS_DIR, recursive = TRUE)
if (!dir.exists(dirname(OUTPUT_DETECTION_CSV))) dir.create(dirname(OUTPUT_DETECTION_CSV), recursive = TRUE)

cat("--- Loading model parameters ---\n")

# ==============================================================================
# 1. Genetic test costs
# ==============================================================================
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

# ==============================================================================
# 2. Cascade parameters (epidemiological counts, loaded from config.yml)
# ==============================================================================
# Kept separate from costs_list to prevent accidental Gamma-sampling in PSA.
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

# ==============================================================================
# 3. Uptake parameters
# ==============================================================================
cat("\nLoading uptake parameters from", INPUT_UPTAKE, "...\n")

raw_uptake <- read_csv(INPUT_UPTAKE, show_col_types = FALSE)

required_uptake_cols <- c(
    "parameter", "base_case", "beta_shape1", "beta_shape2",
    "dsa_min", "dsa_max", "source"
)
missing_uptake_cols <- setdiff(required_uptake_cols, names(raw_uptake))
if (length(missing_uptake_cols) > 0) {
    stop("Uptake table missing required columns: ", paste(missing_uptake_cols, collapse = ", "))
}

required_uptake_ids <- c("reflex_uptake", "cascade_uptake")
missing_uptake_ids <- setdiff(required_uptake_ids, raw_uptake$parameter)
if (length(missing_uptake_ids) > 0) {
    stop("Uptake table missing required parameters: ", paste(missing_uptake_ids, collapse = ", "))
}

for (idx in seq_len(nrow(raw_uptake))) {
    pid <- raw_uptake$parameter[idx]
    bc <- raw_uptake$base_case[idx]
    if (!is.numeric(bc) || bc < 0 || bc > 1) {
        stop(sprintf("Uptake parameter '%s': base_case must be in [0, 1], got %s", pid, bc))
    }
    s1 <- raw_uptake$beta_shape1[idx]
    s2 <- raw_uptake$beta_shape2[idx]
    if (!is.na(s1) && (!is.numeric(s1) || s1 <= 0)) {
        stop(sprintf("Uptake parameter '%s': beta_shape1 must be > 0 when specified", pid))
    }
    if (!is.na(s2) && (!is.numeric(s2) || s2 <= 0)) {
        stop(sprintf("Uptake parameter '%s': beta_shape2 must be > 0 when specified", pid))
    }
    dmin <- raw_uptake$dsa_min[idx]
    dmax <- raw_uptake$dsa_max[idx]
    if (!is.na(dmin) && !is.na(dmax)) {
        if (dmin < 0 || dmin > 1 || dmax < 0 || dmax > 1) {
            stop(sprintf("Uptake parameter '%s': dsa_min and dsa_max must be in [0, 1]", pid))
        }
        if (dmin > bc || bc > dmax) {
            stop(sprintf("Uptake parameter '%s': must have dsa_min <= base_case <= dsa_max", pid))
        }
    }
}

uptake_params <- list()
for (idx in seq_len(nrow(raw_uptake))) {
    pid <- raw_uptake$parameter[idx]
    key <- sub("_uptake$", "", pid)
    uptake_params[[key]] <- list(
        base_case   = raw_uptake$base_case[idx],
        beta_shape1 = raw_uptake$beta_shape1[idx],
        beta_shape2 = raw_uptake$beta_shape2[idx],
        dsa_min     = raw_uptake$dsa_min[idx],
        dsa_max     = raw_uptake$dsa_max[idx]
    )
}

cat("Uptake parameters loaded:\n")
for (key in names(uptake_params)) {
    up <- uptake_params[[key]]
    if (!is.na(up$beta_shape1)) {
        cat(sprintf("  %s: base = %.2f, Beta(%.0f, %.0f), DSA range = %.2f–%.2f\n",
            key, up$base_case, up$beta_shape1, up$beta_shape2,
            up$dsa_min, up$dsa_max))
    } else {
        cat(sprintf("  %s: base = %.2f (fixed, not sampled in PSA)\n",
            key, up$base_case))
    }
}

# ==============================================================================
# 4. Gene panels
# ==============================================================================
cat("\nLoading panels from", INPUT_PANELS_DIR, "...\n")

read_panel_genes <- function(filename) {
    path <- file.path(INPUT_PANELS_DIR, filename)
    if (!file.exists(path)) stop(paste("Panel file not found:", path))

    df <- read_tsv(path, show_col_types = FALSE)

    if ("Gene Symbol" %in% names(df)) {
        genes <- df$`Gene Symbol`
    } else if ("Entity Name" %in% names(df)) {
        genes <- df$`Entity Name`
    } else {
        stop(paste("Could not find gene column in", filename))
    }

    unique(na.omit(genes))
}

panel_cystic <- read_panel_genes("cystic_kidney_disease_v4.3.tsv")
panel_glomerular <- read_panel_genes("glomerulopathy_v5.2.tsv")
panel_tubulointerstitial <- read_panel_genes("tubulointerstitial_kidney_disease_v3.6.tsv")
panel_ckdu <- read_panel_genes("unexplained_kf_v1.120.tsv")

p_tub_main <- read_panel_genes("tubulopathies_v5.6.tsv")
p_nephro <- read_panel_genes("nephrocalcinosis_5.6.tsv")
panel_tubulopathies <- unique(c(p_tub_main, p_nephro))

all_panel_files <- list.files(INPUT_PANELS_DIR, pattern = "\\.tsv$", full.names = FALSE)
cat("Merging the following", length(all_panel_files), "panel files for Comprehensive list:\n")
print(all_panel_files)

all_genes_list <- lapply(all_panel_files, read_panel_genes)
panel_comprehensive <- unique(unlist(all_genes_list))

cat("Comprehensive Panel Size:", length(panel_comprehensive), "unique genes\n")

# CKDu uses the comprehensive panel since unknown etiology requires the broadest yield.
panels_list <- list(
    cystic = panel_cystic,
    glomerular = panel_glomerular,
    tubulointerstitial = panel_tubulointerstitial,
    tubulopathies = panel_tubulopathies,
    CKDu = panel_comprehensive,
    comprehensive = panel_comprehensive
)

cat("\nPanel Sizes (Gene Counts):\n")
print(map_int(panels_list, length))

# ==============================================================================
# 5. VUS probability parameters (Beta distributions by modality)
# ==============================================================================
cat("\nLoading VUS parameters from", INPUT_VUS, "...\n")

vus_raw <- read_excel(INPUT_VUS)

if ("Test modality" %in% names(vus_raw)) {
    vus_raw <- vus_raw %>% rename(test_modality = `Test modality`)
} else if ("modality" %in% names(vus_raw)) {
    vus_raw <- vus_raw %>% rename(test_modality = modality)
}

cat("Available VUS Modalities:\n")
print(unique(vus_raw$test_modality))

get_vus_params <- function(pattern) {
    row <- vus_raw %>%
        filter(str_detect(test_modality, regex(pattern, ignore_case = TRUE)))
    if (nrow(row) != 1) {
        stop(paste0("Could not uniquely identify row for pattern: '", pattern, "' - Found ", nrow(row), " matches."))
    }
    list(
        shape1 = row$alpha,
        shape2 = row$beta,
        mean_sens = row$alpha / (row$alpha + row$beta)
    )
}

panel_bins <- list(
    bin_11_25   = get_vus_params("11-25"),
    bin_26_50   = get_vus_params("26-50"),
    bin_51_100  = get_vus_params("51-100"),
    bin_101_200 = get_vus_params("101-200"),
    bin_gt_200  = get_vus_params(">200")
)

es_vus <- get_vus_params("Exome sequencing")
gs_vus <- get_vus_params("Genome sequencing")

vus_params <- list(
    panels = panel_bins,
    es = es_vus,
    gs = gs_vus,
    metadata = list(
        timestamp = Sys.time(),
        source_file = INPUT_VUS
    )
)

# ==============================================================================
# 6. Analytic detection parameters (Beta distributions by modality and variant class)
# ==============================================================================
cat("\nLoading detection parameters from", INPUT_DETECTION, "sheet:", INPUT_DETECTION_TAB, "...\n")

raw_detection <- read_excel(INPUT_DETECTION, sheet = INPUT_DETECTION_TAB)

required_det_cols <- c("modality_sim", "variant_class_sim", "alpha", "beta", "mean_sensitivity")
missing_det_cols <- setdiff(required_det_cols, names(raw_detection))
if (length(missing_det_cols) > 0) {
    stop("Detection table missing required columns: ", paste(missing_det_cols, collapse = ", "))
}

cat("Loaded", nrow(raw_detection), "rows from curated detection data.\n\n")

# Modality mapping:
#   PanelExome -> Panel and ES (shared detection characteristics)
#   Genome     -> GS
#   MUC1_assay -> specialized assay for MUC1 VNTR
#   PKD1_assay -> specialized assay for PKD1 pseudogene region
# Variant class mapping:
#   SNV_INDEL       -> SNV_Indel
#   CNV_SV          -> CNV_SV
#   PKD1_pseudogene -> Difficult_Locus (PKD1)
#   MUC1_VNTR       -> Difficult_Locus (MUC1)

detection_clean <- raw_detection %>%
    select(modality_sim, variant_class_sim, alpha, beta, mean_sensitivity) %>%
    mutate(
        variant_class_model = case_when(
            variant_class_sim == "SNV_INDEL" ~ "SNV_Indel",
            variant_class_sim == "CNV_SV" ~ "CNV_SV",
            variant_class_sim %in% c("PKD1_pseudogene", "MUC1_VNTR") ~ "Difficult_Locus",
            TRUE ~ variant_class_sim
        ),
        difficult_locus_type = case_when(
            variant_class_sim == "PKD1_pseudogene" ~ "PKD1",
            variant_class_sim == "MUC1_VNTR" ~ "MUC1",
            TRUE ~ NA_character_
        )
    )

get_beta_params <- function(data, modality, variant_class_orig) {
    row <- data %>%
        filter(modality_sim == modality, variant_class_sim == variant_class_orig)

    if (nrow(row) == 0) {
        return(list(alpha = 0, beta = 1, mean_sensitivity = 0, applicable = FALSE))
    }
    if (is.na(row$alpha) || is.na(row$beta)) {
        return(list(alpha = 0, beta = 1, mean_sensitivity = 0, applicable = FALSE))
    }

    list(
        alpha = row$alpha,
        beta = row$beta,
        mean_sensitivity = row$mean_sensitivity,
        applicable = TRUE
    )
}

# Phenotype-specific panel detection assumptions:
#   CYSTIC panel: includes PKD1 assay (not MUC1)
#   TUBULOINTERSTITIAL (ADTKD) panel: includes MUC1 assay (not PKD1)
#   GLOMERULAR, TUBULOPATHIES, COMPREHENSIVE (CKDu) panels: no bundled assays
# ES and GS strategies use only intrinsic detection capabilities.

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

panel_by_phenotype <- list(
    cystic = list(
        SNV_Indel = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = make_difficult_locus(pkd1_use_assay = TRUE, muc1_use_assay = FALSE),
        bundled_assays = c("PKD1")
    ),
    tubulointerstitial = list(
        SNV_Indel = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = make_difficult_locus(pkd1_use_assay = FALSE, muc1_use_assay = TRUE),
        bundled_assays = c("MUC1")
    ),
    glomerular = list(
        SNV_Indel = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = make_difficult_locus(pkd1_use_assay = FALSE, muc1_use_assay = FALSE),
        bundled_assays = character(0)
    ),
    tubulopathies = list(
        SNV_Indel = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = make_difficult_locus(pkd1_use_assay = FALSE, muc1_use_assay = FALSE),
        bundled_assays = character(0)
    ),
    CKDu = list(
        SNV_Indel = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
        CNV_SV = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
        Difficult_Locus = make_difficult_locus(pkd1_use_assay = FALSE, muc1_use_assay = FALSE),
        bundled_assays = character(0)
    )
)

# ES Augmented (scenario): PKD1 assay applied universally; MUC1 assay only for tubulointerstitial.
# Bundled assay costs are absorbed within the ES unit price (scenario assumption).
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

es_detection <- list(
    SNV_Indel = get_beta_params(detection_clean, "PanelExome", "SNV_INDEL"),
    CNV_SV = get_beta_params(detection_clean, "PanelExome", "CNV_SV"),
    Difficult_Locus = make_difficult_locus(pkd1_use_assay = FALSE, muc1_use_assay = FALSE),
    bundled_assays = character(0)
)

gs_detection <- list(
    SNV_Indel = get_beta_params(detection_clean, "Genome", "SNV_INDEL"),
    CNV_SV = get_beta_params(detection_clean, "Genome", "CNV_SV"),
    Difficult_Locus = list(
        PKD1 = get_beta_params(detection_clean, "Genome", "PKD1_pseudogene"),
        MUC1 = get_beta_params(detection_clean, "Genome", "MUC1_VNTR")
    ),
    bundled_assays = character(0)
)

strategy_detection_matrix <- list(
    Panel        = panel_by_phenotype,
    ES           = es_detection,
    GS           = gs_detection,
    ES_Augmented = es_augmented_by_phenotype
)

detection_matrix <- strategy_detection_matrix

summary_rows <- list()
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

for (pheno in names(panel_by_phenotype)) {
    summary_rows <- c(summary_rows, add_detection_rows("Panel", panel_by_phenotype[[pheno]], pheno))
}
summary_rows <- c(summary_rows, add_detection_rows("ES", es_detection))
summary_rows <- c(summary_rows, add_detection_rows("GS", gs_detection))
for (pheno in names(es_augmented_by_phenotype)) {
    summary_rows <- c(summary_rows, add_detection_rows("ES_Augmented", es_augmented_by_phenotype[[pheno]], pheno))
}

summary_df <- bind_rows(summary_rows)
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

write.csv(summary_export, OUTPUT_DETECTION_CSV, row.names = FALSE)

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

detection_obj <- list(
    strategy_detection_matrix = strategy_detection_matrix,
    raw_data = detection_clean,
    metadata = list(
        timestamp = Sys.time(),
        source_file = INPUT_DETECTION,
        source_sheet = INPUT_DETECTION_TAB,
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

# ==============================================================================
# 7. Save consolidated parameter objects
# ==============================================================================
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

saveRDS(parameters,     OUTPUT_INTERMEDIATE)
saveRDS(costs_list,     file.path(PARAMS_DIR, "costs.rds"))
saveRDS(panels_list,    file.path(PARAMS_DIR, "panels.rds"))
saveRDS(cascade_params, file.path(PARAMS_DIR, "cascade_params.rds"))
saveRDS(uptake_params,  file.path(PARAMS_DIR, "uptake_params.rds"))
saveRDS(vus_params,     file.path(PARAMS_DIR, "vus_betas.rds"))
saveRDS(detection_obj,  file.path(PARAMS_DIR, "detection_betas.rds"))

cat("\nParameters saved to:\n")
cat("  -", OUTPUT_INTERMEDIATE, "\n")
cat("  -", file.path(PARAMS_DIR, "costs.rds"), "\n")
cat("  -", file.path(PARAMS_DIR, "panels.rds"), "\n")
cat("  -", file.path(PARAMS_DIR, "cascade_params.rds"), "\n")
cat("  -", file.path(PARAMS_DIR, "uptake_params.rds"), "\n")
cat("  -", file.path(PARAMS_DIR, "vus_betas.rds"), "\n")
cat("  -", file.path(PARAMS_DIR, "detection_betas.rds"), "\n")
cat("  -", OUTPUT_DETECTION_CSV, "\n")
