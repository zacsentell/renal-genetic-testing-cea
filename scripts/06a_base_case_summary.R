# 06a_base_case_summary.R
# Purpose: Generate base-case outcomes summary table (§7.3) and cost composition (§7.4)

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
source("scripts/utils_schema.R")

# ==============================================================================
# 1. Configuration
# ==============================================================================
INPUT_FILE <- "outputs/results/base_case/iteration_level/strategy_iteration_outcomes.csv"
OUTPUT_FILE <- "outputs/results/base_case/summary_tables/base_case_outcomes_by_strategy.csv"

COMP_INPUT_FILE <- "outputs/results/base_case/cost_composition/iteration_cost_components.csv"
COMP_OUTPUT_TABLE <- "outputs/results/base_case/cost_composition/base_case_cost_components_by_strategy.csv"
FIG_BASE <- "outputs/results/base_case/cost_composition/mean_cost_composition_by_strategy_base_case"

# ==============================================================================
# 2. Load Data
# ==============================================================================
cat("Loading iteration-level data from:", INPUT_FILE, "\n")
iter_data <- read_csv(INPUT_FILE, show_col_types = FALSE)

iter_required <- c(
    "iteration_id", "strategy_id", "strategy_label",
    "total_cost_per_proband_cad", "diagnoses_per_proband",
    "cost_per_diagnosis_cad", "pr_at_least_one_vus"
)
assert_required_columns(iter_data, iter_required, "strategy_iteration_outcomes")

# ==============================================================================
# 3. Compute Summary Statistics (§7.3)
# ==============================================================================
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

write_csv_validated(summary_stats, OUTPUT_FILE, "base_case_outcomes_by_strategy")

# ==============================================================================
# 4. Cost Composition (§7.4)
# ==============================================================================
cat("Loading cost components from:", COMP_INPUT_FILE, "\n")
comp_data <- read_csv(COMP_INPUT_FILE, show_col_types = FALSE)

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

write_csv_validated(comp_summary_long, COMP_OUTPUT_TABLE, "base_case_cost_components_by_strategy")

# Plot from saved long table-style structure
plot_data <- comp_summary_long %>%
    mutate(
        component_name = factor(component_name,
            levels = c("Genetics Visits", "Laboratory Testing", "Cascade Testing", "VUS Follow-up")
        ),
        strategy_label = factor(strategy_label,
            levels = comp_summary_wide$strategy_label
        )
    )

errorbars <- comp_summary_wide %>%
    select(strategy_label, total_cost_per_proband_cad_mean, total_cost_per_proband_cad_ui_low, total_cost_per_proband_cad_ui_high)

p <- ggplot(plot_data, aes(x = strategy_label, y = component_cost_per_proband_cad_mean, fill = component_name)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    geom_errorbar(
        data = errorbars,
        aes(
            x = strategy_label,
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

ggsave(paste0(FIG_BASE, ".png"), p, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(paste0(FIG_BASE, ".svg"), p, width = 8, height = 6, bg = "white")

cat("Exported:\n")
cat(" -", OUTPUT_FILE, "\n")
cat(" -", COMP_OUTPUT_TABLE, "\n")
cat(" -", paste0(FIG_BASE, ".png"), "\n")
cat(" -", paste0(FIG_BASE, ".svg"), "\n")
