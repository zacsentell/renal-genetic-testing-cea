# 07_global_sensitivity_analysis.R
# Purpose: Perform Global Sensitivity Analysis (PRCC) for Section 8.2

if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("sensitivity")) install.packages("sensitivity")
if (!require("readr")) install.packages("readr")

library(dplyr)
library(tidyr)
library(sensitivity)
library(readr)

# ==============================================================================
# 1. Setup and Inputs
# ==============================================================================
INPUT_BASE_OUTCOMES <- "outputs/results/base_case/iteration_level/strategy_iteration_outcomes.csv"
INPUT_PARAM_TRACE <- "outputs/results/uncertainty_sensitivity/global_sensitivity_analysis_prcc/iteration_parameter_trace.csv"
INPUT_GS_SCENARIO_ITER <- "outputs/results/scenario_analysis/gs_uplift/gs_scenario_iteration_outcomes.csv"
OUTPUT_DIR <- "outputs/results/uncertainty_sensitivity/global_sensitivity_analysis_prcc"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

if (!file.exists(INPUT_BASE_OUTCOMES)) stop("Iteration outcomes not found: ", INPUT_BASE_OUTCOMES)
if (!file.exists(INPUT_PARAM_TRACE)) stop("Parameter trace not found: ", INPUT_PARAM_TRACE)
if (!file.exists(INPUT_GS_SCENARIO_ITER)) stop("GS scenario iteration outcomes not found: ", INPUT_GS_SCENARIO_ITER)


# ==============================================================================
# 2. Load and Prepare Data
# ==============================================================================
outcomes_df <- read_csv(INPUT_BASE_OUTCOMES, show_col_types = FALSE)
params_df <- read_csv(INPUT_PARAM_TRACE, show_col_types = FALSE)
gs_iter_df <- read_csv(INPUT_GS_SCENARIO_ITER, show_col_types = FALSE)

params_wide <- params_df %>%
    pivot_wider(names_from = parameter_name, values_from = parameter_value)

base_wide <- outcomes_df %>%
    select(iteration_id, strategy_label, total_cost_per_proband_cad, diagnoses_per_proband) %>%
    pivot_wider(names_from = strategy_label, values_from = c(total_cost_per_proband_cad, diagnoses_per_proband))

gs_scenario_only <- gs_iter_df %>%
    filter(grepl("^uplift_", scenario_label)) %>%
    transmute(
        iteration_id,
        total_cost_per_proband_cad_GS_scenario = total_cost_per_proband_cad,
        diagnoses_per_proband_GS_scenario = diagnoses_per_proband
    )

comparison_frame <- base_wide %>%
    left_join(gs_scenario_only, by = "iteration_id")

comparisons <- list(
    list(
        focal_cost_col = "total_cost_per_proband_cad_Panel_Reflex_ES",
        comp_cost_col  = "total_cost_per_proband_cad_Panel",
        focal_diag_col = "diagnoses_per_proband_Panel_Reflex_ES",
        comp_diag_col  = "diagnoses_per_proband_Panel",
        label          = "reflex_vs_panel",
        title_focal    = "Reflex (Panel→ES)",
        title_comp     = "Panel"
    ),
    list(
        focal_cost_col = "total_cost_per_proband_cad_ES",
        comp_cost_col  = "total_cost_per_proband_cad_Panel",
        focal_diag_col = "diagnoses_per_proband_ES",
        comp_diag_col  = "diagnoses_per_proband_Panel",
        label          = "es_vs_panel",
        title_focal    = "ES",
        title_comp     = "Panel"
    ),
    list(
        focal_cost_col = "total_cost_per_proband_cad_GS_scenario",
        comp_cost_col  = "total_cost_per_proband_cad_Panel_Reflex_ES",
        focal_diag_col = "diagnoses_per_proband_GS_scenario",
        comp_diag_col  = "diagnoses_per_proband_Panel_Reflex_ES",
        label          = "genome_vs_reflex",
        title_focal    = "GS (+10% yield)",
        title_comp     = "Reflex (Panel→ES)"
    )
)

# ==============================================================================
# 3. Run PRCC
# ==============================================================================
for (comp in comparisons) {
    req <- c(comp$focal_cost_col, comp$comp_cost_col, comp$focal_diag_col, comp$comp_diag_col)
    missing <- setdiff(req, names(comparison_frame))
    if (length(missing) > 0) {
        warning("Skipping ", comp$label, " due to missing columns: ", paste(missing, collapse = ", "))
        next
    }

    analysis_df <- comparison_frame %>%
        mutate(
            delta_cost_cad = .data[[comp$focal_cost_col]] - .data[[comp$comp_cost_col]],
            delta_diagnoses = .data[[comp$focal_diag_col]] - .data[[comp$comp_diag_col]]
        ) %>%
        select(iteration_id, delta_cost_cad, delta_diagnoses) %>%
        inner_join(params_wide, by = "iteration_id")

    X <- analysis_df %>%
        select(-iteration_id, -delta_cost_cad, -delta_diagnoses) %>%
        select(where(is.numeric))

    names(X) <- make.names(names(X), unique = TRUE)

    col_vars <- apply(X, 2, function(x) var(x, na.rm = TRUE))
    X <- X[, !(is.na(col_vars) | col_vars == 0), drop = FALSE]

    if (ncol(X) < 2) {
        warning("Skipping ", comp$label, ": not enough non-constant parameters")
        next
    }

    nboot <- 500
    cat("Running PRCC for", comp$label, "(nboot=", nboot, ")\n", sep = "")
    pcc_cost <- pcc(X, analysis_df$delta_cost_cad, rank = TRUE, nboot = nboot, conf = 0.95)
    pcc_diag <- pcc(X, analysis_df$delta_diagnoses, rank = TRUE, nboot = nboot, conf = 0.95)

    prcc_cost <- pcc_cost$PRCC
    prcc_diag <- pcc_diag$PRCC

    prcc_table <- data.frame(
        parameter_name = rownames(prcc_cost),
        prcc_delta_cost_cad = prcc_cost$original,
        prcc_delta_cost_cad_ci_low = prcc_cost$`min. c.i.`,
        prcc_delta_cost_cad_ci_high = prcc_cost$`max. c.i.`,
        prcc_delta_diagnoses = prcc_diag$original,
        prcc_delta_diagnoses_ci_low = prcc_diag$`min. c.i.`,
        prcc_delta_diagnoses_ci_high = prcc_diag$`max. c.i.`,
        stringsAsFactors = FALSE
    ) %>%
        arrange(desc(abs(prcc_delta_cost_cad)))

    out_table <- file.path(OUTPUT_DIR, paste0("prcc_results_", comp$label, ".csv"))
    write_csv(prcc_table, out_table)

    cat("Saved PRCC results for", comp$label, "\n")
}

# ==============================================================================
# 4. Note: composite figure
# ==============================================================================
# The composite PRCC figure (forest plot + heatmap) is built by 07b_prcc_composite_figure.R,
# which reads the CSV outputs produced above. Run 07b after this script.

cat("Global sensitivity analysis complete. Run 07b_prcc_composite_figure.R to build the composite figure.\n")
