# scripts/08_base_case.R
# Purpose: Base-case summary, cost composition, incremental analysis, and efficiency frontier.
# Author: Zachary Sentell

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggrepel)
source("scripts/utils_schema.R")
source("scripts/utils_cea.R")

# ==============================================================================
# IO Paths
# ==============================================================================
INPUT_ITER_OUTCOMES <- "outputs/results/base_case/iteration_level/strategy_iteration_outcomes.csv"
INPUT_ITER_COMP     <- "outputs/results/base_case/cost_composition/iteration_cost_components.csv"

OUTPUT_SUMMARY      <- "outputs/results/base_case/summary_tables/base_case_outcomes_by_strategy.csv"
OUTPUT_COMP_TABLE   <- "outputs/results/base_case/cost_composition/base_case_cost_components_by_strategy.csv"
COMP_FIG_BASE       <- "outputs/results/base_case/cost_composition/mean_cost_composition_by_strategy_base_case"

INCREMENTAL_DIR     <- "outputs/results/incremental_analysis/base_case"
OUTPUT_INCREMENTAL  <- file.path(INCREMENTAL_DIR, "incremental_table_base_case.csv")
OUTPUT_VS_PANEL     <- file.path(INCREMENTAL_DIR, "incremental_vs_panel.csv")
FRONTIER_FIG_BASE   <- file.path(INCREMENTAL_DIR, "efficiency_frontier_base_case")

if (!dir.exists(INCREMENTAL_DIR)) dir.create(INCREMENTAL_DIR, recursive = TRUE)

# ==============================================================================
# 1. Summary statistics by strategy
# ==============================================================================
cat("Loading iteration-level data from:", INPUT_ITER_OUTCOMES, "\n")
iter_data <- read_csv(INPUT_ITER_OUTCOMES, show_col_types = FALSE)

iter_required <- c(
    "iteration_id", "strategy_id", "strategy_label",
    "total_cost_per_proband_cad", "diagnoses_per_proband",
    "cost_per_diagnosis_cad", "pr_at_least_one_vus"
)
assert_required_columns(iter_data, iter_required, "strategy_iteration_outcomes")

summary_stats <- iter_data %>%
    group_by(strategy_id, strategy_label) %>%
    summarise(
        total_cost_per_proband_cad_mean = mean(total_cost_per_proband_cad),
        total_cost_per_proband_cad_ui_low = quantile(total_cost_per_proband_cad, 0.025),
        total_cost_per_proband_cad_ui_high = quantile(total_cost_per_proband_cad, 0.975),
        diagnoses_per_proband_mean = mean(diagnoses_per_proband),
        diagnoses_per_proband_ui_low = quantile(diagnoses_per_proband, 0.025),
        diagnoses_per_proband_ui_high = quantile(diagnoses_per_proband, 0.975),
        cost_per_diagnosis_cad_mean = mean(cost_per_diagnosis_cad, na.rm = TRUE),
        cost_per_diagnosis_cad_ui_low = quantile(cost_per_diagnosis_cad, 0.025, na.rm = TRUE),
        cost_per_diagnosis_cad_ui_high = quantile(cost_per_diagnosis_cad, 0.975, na.rm = TRUE),
        pr_at_least_one_vus_mean = mean(pr_at_least_one_vus),
        pr_at_least_one_vus_ui_low = quantile(pr_at_least_one_vus, 0.025),
        pr_at_least_one_vus_ui_high = quantile(pr_at_least_one_vus, 0.975),
        .groups = "drop"
    ) %>%
    arrange(total_cost_per_proband_cad_mean)

summary_required <- c(
    "strategy_id", "strategy_label",
    "total_cost_per_proband_cad_mean", "total_cost_per_proband_cad_ui_low", "total_cost_per_proband_cad_ui_high",
    "diagnoses_per_proband_mean", "diagnoses_per_proband_ui_low", "diagnoses_per_proband_ui_high",
    "cost_per_diagnosis_cad_mean", "cost_per_diagnosis_cad_ui_low", "cost_per_diagnosis_cad_ui_high",
    "pr_at_least_one_vus_mean", "pr_at_least_one_vus_ui_low", "pr_at_least_one_vus_ui_high"
)
assert_required_columns(summary_stats, summary_required, "base_case_outcomes_by_strategy", exact = TRUE)
assert_no_na(summary_stats, c("strategy_id", "strategy_label"), "base_case_outcomes_by_strategy")

write_csv_validated(summary_stats, OUTPUT_SUMMARY, "base_case_outcomes_by_strategy")

# ==============================================================================
# 2. Cost composition by strategy
# ==============================================================================
cat("Loading cost components from:", INPUT_ITER_COMP, "\n")
comp_data <- read_csv(INPUT_ITER_COMP, show_col_types = FALSE)

component_iter_required <- c(
    "iteration_id", "strategy_id", "strategy_label",
    "cost_consultation_cad", "cost_testing_cad",
    "cost_cascade_cad", "cost_vus_followup_cad", "total_cost_per_proband_cad"
)
assert_required_columns(comp_data, component_iter_required, "iteration_cost_components")

comp_summary_wide <- comp_data %>%
    group_by(strategy_id, strategy_label) %>%
    summarise(
        cost_consultation_mean = mean(cost_consultation_cad),
        cost_testing_mean = mean(cost_testing_cad),
        cost_cascade_mean = mean(cost_cascade_cad),
        cost_vus_followup_mean = mean(cost_vus_followup_cad),
        total_cost_per_proband_cad_mean = mean(total_cost_per_proband_cad),
        total_cost_per_proband_cad_ui_low = quantile(total_cost_per_proband_cad, 0.025),
        total_cost_per_proband_cad_ui_high = quantile(total_cost_per_proband_cad, 0.975),
        .groups = "drop"
    ) %>%
    arrange(total_cost_per_proband_cad_mean)

component_label_map <- c(
    cost_consultation_mean = "Genetics Visits",
    cost_testing_mean = "Laboratory Testing",
    cost_cascade_mean = "Cascade Testing",
    cost_vus_followup_mean = "VUS Follow-up"
)

comp_summary_long <- comp_summary_wide %>%
    pivot_longer(
        cols = c(cost_consultation_mean, cost_testing_mean, cost_cascade_mean, cost_vus_followup_mean),
        names_to = "component_name",
        values_to = "component_cost_per_proband_cad_mean"
    ) %>%
    mutate(component_name = unname(component_label_map[component_name])) %>%
    select(
        strategy_id,
        strategy_label,
        component_name,
        component_cost_per_proband_cad_mean,
        total_cost_per_proband_cad_mean,
        total_cost_per_proband_cad_ui_low,
        total_cost_per_proband_cad_ui_high
    )

comp_required <- c(
    "strategy_id", "strategy_label", "component_name", "component_cost_per_proband_cad_mean",
    "total_cost_per_proband_cad_mean", "total_cost_per_proband_cad_ui_low", "total_cost_per_proband_cad_ui_high"
)
assert_required_columns(comp_summary_long, comp_required, "base_case_cost_components_by_strategy", exact = TRUE)
assert_no_na(comp_summary_long, c("strategy_id", "strategy_label", "component_name"), "base_case_cost_components_by_strategy")

write_csv_validated(comp_summary_long, OUTPUT_COMP_TABLE, "base_case_cost_components_by_strategy")

strategy_display_order <- strategy_display_label(comp_summary_wide$strategy_label)

plot_data <- comp_summary_long %>%
    mutate(
        component_name = factor(component_name,
            levels = c("Genetics Visits", "Laboratory Testing", "Cascade Testing", "VUS Follow-up")
        ),
        strategy_display = factor(strategy_display_label(strategy_label),
            levels = strategy_display_order
        )
    )

errorbars <- comp_summary_wide %>%
    mutate(strategy_display = factor(strategy_display_label(strategy_label),
        levels = strategy_display_order)) %>%
    select(strategy_display, total_cost_per_proband_cad_mean, total_cost_per_proband_cad_ui_low, total_cost_per_proband_cad_ui_high)

p_comp <- ggplot(plot_data, aes(x = strategy_display, y = component_cost_per_proband_cad_mean, fill = component_name)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    geom_errorbar(
        data = errorbars,
        aes(
            x = strategy_display,
            y = total_cost_per_proband_cad_mean,
            ymin = total_cost_per_proband_cad_ui_low,
            ymax = total_cost_per_proband_cad_ui_high
        ),
        inherit.aes = FALSE,
        width = 0.2,
        linewidth = 0.6
    ) +
    scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.85, direction = -1) +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.05)),
        labels = scales::dollar_format(prefix = "$", suffix = "")
    ) +
    labs(
        title = "Mean Cost Composition by Strategy (Base Case)",
        x = NULL,
        y = "Cost per Proband (CAD)",
        fill = "Cost Component"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        legend.position = "right",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )

ggsave(paste0(COMP_FIG_BASE, ".png"), p_comp, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(paste0(COMP_FIG_BASE, ".svg"), p_comp, width = 8, height = 6, bg = "white")

# ==============================================================================
# 3. Incremental analysis and efficiency frontier
# ==============================================================================
cat("\n--- Performing Incremental Analysis ---\n")

analysis_data <- summary_stats %>%
    select(
        strategy_id, strategy_label,
        cost = total_cost_per_proband_cad_mean,
        effect = diagnoses_per_proband_mean,
        cost_ui_low = total_cost_per_proband_cad_ui_low,
        cost_ui_high = total_cost_per_proband_cad_ui_high,
        effect_ui_low = diagnoses_per_proband_ui_low,
        effect_ui_high = diagnoses_per_proband_ui_high
    )

vus_data <- summary_stats %>%
    select(strategy_id, vus_prob = pr_at_least_one_vus_mean)

results <- perform_incremental_analysis(analysis_data) %>%
    left_join(vus_data, by = "strategy_id")

output_table <- results %>%
    arrange(cost) %>%
    mutate(
        nnt_probands_per_additional_diagnosis_mean = ifelse(
            !is.na(incremental_effect) & incremental_effect > 0,
            1 / incremental_effect,
            NA_real_
        )
    ) %>%
    select(
        strategy_id,
        strategy_label,
        total_cost_per_proband_cad_mean = cost,
        total_cost_per_proband_cad_ui_low = cost_ui_low,
        total_cost_per_proband_cad_ui_high = cost_ui_high,
        diagnoses_per_proband_mean = effect,
        diagnoses_per_proband_ui_low = effect_ui_low,
        diagnoses_per_proband_ui_high = effect_ui_high,
        incremental_cost_vs_prev_cad_mean = incremental_cost,
        incremental_diagnoses_vs_prev_mean = incremental_effect,
        icpd_cad_per_additional_diagnosis_mean = icer,
        nnt_probands_per_additional_diagnosis_mean,
        dominance_status
    )

incremental_required <- c(
    "strategy_id", "strategy_label",
    "total_cost_per_proband_cad_mean", "total_cost_per_proband_cad_ui_low", "total_cost_per_proband_cad_ui_high",
    "diagnoses_per_proband_mean", "diagnoses_per_proband_ui_low", "diagnoses_per_proband_ui_high",
    "incremental_cost_vs_prev_cad_mean", "incremental_diagnoses_vs_prev_mean",
    "icpd_cad_per_additional_diagnosis_mean", "nnt_probands_per_additional_diagnosis_mean",
    "dominance_status"
)
assert_required_columns(output_table, incremental_required, "incremental_table_base_case", exact = TRUE)
assert_no_na(output_table, c("strategy_id", "strategy_label", "dominance_status"), "incremental_table_base_case")
assert_frontier_first_incremental_na(output_table, table_name = "incremental_table_base_case")

write_csv_validated(output_table, OUTPUT_INCREMENTAL, "incremental_table_base_case")
cat("\nIncremental analysis table exported to:", OUTPUT_INCREMENTAL, "\n")

# ==============================================================================
# 3b. Incremental outcomes vs Panel (reference) for every other strategy
# ==============================================================================
# Unlike the frontier table (incremental vs the next non-dominated strategy), this
# table reports every alternate strategy's incremental cost, yield, and cost per
# additional diagnosis relative to the Panel reference, for direct report display.
panel_ref <- summary_stats %>% filter(strategy_label == "Panel")
if (nrow(panel_ref) != 1) {
    stop("Expected exactly one Panel strategy row for incremental-vs-Panel computation")
}

vs_panel_table <- summary_stats %>%
    filter(strategy_label != "Panel") %>%
    transmute(
        strategy_id,
        strategy_label,
        strategy_display = strategy_display_label(strategy_label),
        total_cost_per_proband_cad_mean,
        diagnoses_per_proband_mean,
        cost_per_diagnosis_cad_mean,
        incremental_cost_vs_panel_cad_mean =
            total_cost_per_proband_cad_mean - panel_ref$total_cost_per_proband_cad_mean,
        incremental_diagnoses_vs_panel_mean =
            diagnoses_per_proband_mean - panel_ref$diagnoses_per_proband_mean
    ) %>%
    mutate(
        icpd_vs_panel_cad_per_additional_diagnosis_mean = ifelse(
            incremental_diagnoses_vs_panel_mean > 0,
            incremental_cost_vs_panel_cad_mean / incremental_diagnoses_vs_panel_mean,
            NA_real_
        )
    ) %>%
    arrange(strategy_id)

write_csv_validated(vs_panel_table, OUTPUT_VS_PANEL, "incremental_vs_panel")
cat("\nIncremental-vs-Panel table exported to:", OUTPUT_VS_PANEL, "\n")
print(as.data.frame(vs_panel_table %>%
    select(strategy_display, incremental_cost_vs_panel_cad_mean,
           incremental_diagnoses_vs_panel_mean,
           icpd_vs_panel_cad_per_additional_diagnosis_mean)), row.names = FALSE)

cat("\n=== Incremental Analysis Results ===\n\n")
print(
    as.data.frame(output_table %>%
        select(
            strategy_label, total_cost_per_proband_cad_mean,
            diagnoses_per_proband_mean, icpd_cad_per_additional_diagnosis_mean,
            dominance_status
        )),
    row.names = FALSE
)

# ==============================================================================
# 4. Efficiency frontier plot with VUS-burden bubble sizing
# ==============================================================================
cat("\n--- Generating Efficiency Frontier Plot ---\n")

results <- results %>%
    mutate(
        strategy_display = strategy_display_label(strategy_label),
        strategy_label_full = paste0(
            strategy_display, "\n(VUS: ", round(vus_prob * 100, 0), "%)"
        )
    )

plot_data <- results %>%
    mutate(
        is_dominated = dominance_status != "non_dominated",
        alpha_val = ifelse(is_dominated, 0.5, 1.0),
        bubble_size = 3 + (vus_prob * 8)
    )

frontier_data <- results %>%
    filter(dominance_status == "non_dominated") %>%
    arrange(effect)

p_frontier <- ggplot(plot_data, aes(x = effect, y = cost)) +
    geom_errorbar(
        aes(ymin = cost_ui_low, ymax = cost_ui_high, alpha = alpha_val),
        width = 0.005, color = "gray60", linewidth = 0.5
    ) +
    geom_errorbarh(
        aes(xmin = effect_ui_low, xmax = effect_ui_high, alpha = alpha_val),
        height = 50, color = "gray60", linewidth = 0.5
    ) +
    geom_line(
        data = frontier_data,
        aes(x = effect, y = cost),
        color = "black", linewidth = 1.0, linetype = "solid"
    ) +
    geom_point(
        aes(
            color = strategy_display,
            alpha = alpha_val,
            size = bubble_size
        )
    ) +
    geom_label_repel(
        aes(label = strategy_label_full, alpha = alpha_val),
        size = 4.5, family = "sans", fontface = "bold",
        box.padding = 1.0, point.padding = 0.6,
        min.segment.length = 0.1, segment.size = 0.4,
        segment.color = "gray40",
        max.overlaps = 20,
        force = 3,
        force_pull = 0.4,
        fill = "white",
        label.size = 0.2,
        label.padding = unit(0.25, "lines")
    ) +
    scale_color_viridis_d(option = "D", begin = 0.15, end = 0.85, direction = 1) +
    scale_alpha_identity() +
    scale_size_identity() +
    scale_x_continuous(
        labels = scales::percent_format(accuracy = 1),
        expand = expansion(mult = 0.18)
    ) +
    scale_y_continuous(
        labels = scales::dollar_format(prefix = "$", suffix = ""),
        expand = expansion(mult = c(0.08, 0.12))
    ) +
    labs(
        title = "Efficiency Frontier: Cost vs Diagnostic Yield",
        subtitle = "Bubble size represents VUS burden (larger = higher probability of VUS).\nError bars show 95% UI.",
        x = "Diagnostic Yield (Proportion of Probands Diagnosed)",
        y = "Total Cost per Proband (CAD)",
        color = NULL
    ) +
    theme_minimal(base_size = 13, base_family = "sans") +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.6),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 13, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 11, color = "gray30", hjust = 0),
        legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(15, 15, 10, 10)
    )

ggsave(paste0(FRONTIER_FIG_BASE, ".png"), p_frontier, width = 7, height = 5.5, dpi = 300, bg = "white")
ggsave(paste0(FRONTIER_FIG_BASE, ".svg"), p_frontier, width = 7, height = 5.5, bg = "white")

# Verify ICER monotonicity along the frontier
frontier_check <- output_table %>%
    filter(dominance_status == "non_dominated") %>%
    arrange(total_cost_per_proband_cad_mean)

if (nrow(frontier_check) > 1) {
    icers <- na.omit(frontier_check$icpd_cad_per_additional_diagnosis_mean)
    if (length(icers) > 1) {
        is_monotonic <- all(diff(icers) >= 0)
        if (is_monotonic) {
            cat("PASS: Frontier ICERs are monotonically increasing\n")
        } else {
            warning("Frontier ICERs are NOT monotonically increasing!")
            print(frontier_check %>% select(strategy_label, icpd_cad_per_additional_diagnosis_mean))
        }
    }
}

cat("\nBase-case analysis complete.\n")
cat("Outputs:\n")
cat("  -", OUTPUT_SUMMARY, "\n")
cat("  -", OUTPUT_COMP_TABLE, "\n")
cat("  -", paste0(COMP_FIG_BASE, ".{png,svg}"), "\n")
cat("  -", OUTPUT_INCREMENTAL, "\n")
cat("  -", OUTPUT_VS_PANEL, "\n")
cat("  -", paste0(FRONTIER_FIG_BASE, ".{png,svg}"), "\n")
