# 07_global_sensitivity_analysis.R
# Purpose: Perform Global Sensitivity Analysis (PRCC) for Section 8.2

if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("sensitivity")) install.packages("sensitivity")
if (!require("readr")) install.packages("readr")
if (!require("patchwork")) install.packages("patchwork")

library(dplyr)
library(tidyr)
library(ggplot2)
library(sensitivity)
library(readr)
library(patchwork)

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
# 2. Human-Readable Parameter Labels
# ==============================================================================
param_labels <- c(
    "prev_CKDu" = "CKDu prevalence",
    "prev_cystic" = "Cystic prevalence",
    "prev_glomerular" = "Glomerular prevalence",
    "prev_tubulointerstitial" = "Tubulointerstitial prevalence",
    "prev_tubulopathies" = "Tubulopathies prevalence",
    "yield_CKDu" = "CKDu yield",
    "yield_cystic" = "Cystic yield",
    "yield_glomerular" = "Glomerular yield",
    "yield_tubulointerstitial" = "Tubulointerstitial yield",
    "yield_tubulopathies" = "Tubulopathies yield",
    "prob_vus_panel_mean" = "VUS probability (Panel)",
    "prob_vus_es" = "VUS probability (ES)",
    "prob_vus_gs" = "VUS probability (GS)",
    "det_snv_es" = "SNV detection (ES)",
    "det_cnv_es" = "CNV detection (ES)",
    "det_snv_gs" = "SNV detection (GS)",
    "det_cnv_gs" = "CNV detection (GS)",
    "det_snv_panel_weighted" = "SNV detection (Panel)",
    "det_cnv_panel_weighted" = "CNV detection (Panel)",
    "det_pkd1_es" = "PKD1 detection (ES)",
    "det_pkd1_gs" = "PKD1 detection (GS)",
    "det_muc1_es" = "MUC1 detection (ES)",
    "det_muc1_gs" = "MUC1 detection (GS)",
    "det_pkd1_panel" = "PKD1 detection (Cystic panel)",
    "det_muc1_panel" = "MUC1 detection (ADTKD panel)",
    "cost_panel" = "Panel Cost",
    "cost_es" = "Exome Sequencing Cost",
    "cost_gs" = "Genome Sequencing Cost",
    "cost_familial" = "Targeted Familial Test Cost",
    "cost_consult_pre" = "Pre-test Consultation Cost",
    "cost_consult_post" = "Post-test Consultation Cost"
)

get_label <- function(param_name) {
    clean_name <- gsub("\\.", "_", param_name)
    if (clean_name %in% names(param_labels)) return(param_labels[[clean_name]])
    tools::toTitleCase(gsub("_", " ", clean_name))
}

assign_category <- function(param_name) {
    if (grepl("^prev_", param_name)) return("Phenotype Prevalence")
    if (grepl("^cost_", param_name)) return("Unit Cost")
    if (grepl("det_pkd1|det_muc1", param_name)) return("Difficult Locus Detection")
    if (grepl("^det_", param_name)) return("Detection Sensitivity")
    if (grepl("^yield_", param_name)) return("Diagnostic Yield")
    if (grepl("^prob_vus", param_name)) return("VUS Probability")
    "Other"
}

plot_tornado <- function(df, outcome_col, outcome_label, focal_label, comp_label, top_n = 20) {
    ci_low_col <- paste0(outcome_col, "_ci_low")
    ci_high_col <- paste0(outcome_col, "_ci_high")

    plot_data <- df %>%
        mutate(
            value = .data[[outcome_col]],
            ci_low = .data[[ci_low_col]],
            ci_high = .data[[ci_high_col]],
            is_significant = sign(ci_low) == sign(ci_high),
            display_label = sapply(parameter_name, get_label),
            category = sapply(parameter_name, assign_category)
        ) %>%
        arrange(desc(abs(value))) %>%
        slice_head(n = top_n) %>%
        mutate(display_label = factor(display_label, levels = rev(display_label)))

    max_abs <- max(abs(c(plot_data$ci_low, plot_data$ci_high, plot_data$value)), na.rm = TRUE)
    lim <- max(0.1, ceiling(max_abs * 10) / 10)

    ggplot(plot_data, aes(x = display_label, y = value, fill = category)) +
        geom_bar(stat = "identity", width = 0.7, aes(alpha = is_significant)) +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        geom_errorbar(
            aes(ymin = ci_low, ymax = ci_high, alpha = is_significant),
            width = 0.2, color = "gray20", linewidth = 0.4
        ) +
        coord_flip() +
        scale_fill_brewer(palette = "Set2", name = "Parameter Category") +
        scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1), guide = "none") +
        scale_y_continuous(limits = c(-lim, lim), breaks = scales::pretty_breaks(n = 5)) +
        labs(
            title = paste0("Global Sensitivity: ", outcome_label),
            subtitle = paste0(focal_label, " vs ", comp_label),
            x = NULL,
            y = "Partial Rank Correlation Coefficient (PRCC)"
        ) +
        theme_minimal(base_size = 11) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(color = "black", size = 10),
            axis.title.x = element_text(color = "black", size = 11, face = "bold"),
            plot.title = element_text(size = 13, face = "bold", hjust = 0),
            plot.subtitle = element_text(size = 10, color = "gray30", hjust = 0),
            legend.position = "bottom"
        )
}

# ==============================================================================
# 3. Load and Prepare Data
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
        comp_cost_col = "total_cost_per_proband_cad_Panel",
        focal_diag_col = "diagnoses_per_proband_Panel_Reflex_ES",
        comp_diag_col = "diagnoses_per_proband_Panel",
        label = "reflex_vs_panel",
        title_focal = "Reflex to ES",
        title_comp = "Panel"
    ),
    list(
        focal_cost_col = "total_cost_per_proband_cad_GS_scenario",
        comp_cost_col = "total_cost_per_proband_cad_Panel_Reflex_ES",
        focal_diag_col = "diagnoses_per_proband_GS_scenario",
        comp_diag_col = "diagnoses_per_proband_Panel_Reflex_ES",
        label = "genome_vs_reflex",
        title_focal = "GS (+10% yield)",
        title_comp = "Reflex to ES"
    )
)

# ==============================================================================
# 4. Run PRCC
# ==============================================================================
# Collect plot objects across iterations for the composite figure
prcc_plots <- list()   # keys: "cost_A", "diag_A", "cost_B", "diag_B" in order

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

    min_thr <- 0.02
    plot_cost <- prcc_table %>% filter(abs(prcc_delta_cost_cad) >= min_thr)
    plot_diag <- prcc_table %>%
        mutate(category = sapply(parameter_name, assign_category)) %>%
        filter(category != "Unit Cost", abs(prcc_delta_diagnoses) >= min_thr)

    p_cost <- plot_tornado(plot_cost, "prcc_delta_cost_cad", "Incremental Cost", comp$title_focal, comp$title_comp) +
        labs(caption = paste0("Parameters with |PRCC| < ", min_thr, " are excluded."))

    p_diag <- plot_tornado(plot_diag, "prcc_delta_diagnoses", "Incremental Diagnoses", comp$title_focal, comp$title_comp) +
        labs(caption = paste0("Parameters with |PRCC| < ", min_thr, " are excluded. Unit costs removed for diagnosis outcome."))

    # Save individual plots
    ggsave(file.path(OUTPUT_DIR, paste0("prcc_tornado_", comp$label, "_delta_cost.png")), p_cost,
        width = 7, height = 5.5, dpi = 300, bg = "white"
    )
    ggsave(file.path(OUTPUT_DIR, paste0("prcc_tornado_", comp$label, "_delta_cost.svg")), p_cost,
        width = 7, height = 5.5, bg = "white"
    )
    ggsave(file.path(OUTPUT_DIR, paste0("prcc_tornado_", comp$label, "_delta_diagnoses.png")), p_diag,
        width = 7, height = 5.5, dpi = 300, bg = "white"
    )
    ggsave(file.path(OUTPUT_DIR, paste0("prcc_tornado_", comp$label, "_delta_diagnoses.svg")), p_diag,
        width = 7, height = 5.5, bg = "white"
    )

    # Add concise panel titles for composite; strip captions
    p_cost_panel <- p_cost + labs(
        title    = paste0(comp$title_focal, " vs ", comp$title_comp),
        subtitle = "Outcome: Incremental Cost",
        caption  = NULL
    )
    p_diag_panel <- p_diag + labs(
        title    = paste0(comp$title_focal, " vs ", comp$title_comp),
        subtitle = "Outcome: Incremental Diagnoses",
        caption  = NULL
    )

    prcc_plots[[paste0("cost_", comp$label)]] <- p_cost_panel
    prcc_plots[[paste0("diag_", comp$label)]] <- p_diag_panel

    cat("Saved PRCC outputs for", comp$label, "\n")
}

# ==============================================================================
# 5. Composite 2x2 figure
# ==============================================================================
# Layout: row 1 = Reflex vs Panel (cost | yield), row 2 = GS vs Reflex (cost | yield)
if (length(prcc_plots) == 4) {
    p_A_cost <- prcc_plots[["cost_reflex_vs_panel"]]
    p_A_diag <- prcc_plots[["diag_reflex_vs_panel"]]
    p_B_cost <- prcc_plots[["cost_genome_vs_reflex"]]
    p_B_diag <- prcc_plots[["diag_genome_vs_reflex"]]

    composite <- (p_A_cost | p_A_diag) / (p_B_cost | p_B_diag)
    composite <- composite +
        plot_layout(guides = "collect") +
        plot_annotation(
            tag_levels = "A",
            theme = theme(plot.tag = element_text(size = 13, face = "bold"))
        )

    ggsave(
        file.path(OUTPUT_DIR, "prcc_tornado_composite.png"),
        composite,
        width = 14, height = 10, dpi = 300, bg = "white"
    )
    ggsave(
        file.path(OUTPUT_DIR, "prcc_tornado_composite.svg"),
        composite,
        width = 14, height = 10, bg = "white"
    )
    cat("Saved composite PRCC figure.\n")
} else {
    warning("Expected 4 PRCC plots for composite; found ", length(prcc_plots), ". Composite not saved.")
}

cat("Global sensitivity analysis complete.\n")
