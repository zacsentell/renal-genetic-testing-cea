# 06c_cystic_subgroup.R
# Purpose: Cystic phenotype subgroup analysis (§7.6)

library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(tidyr)
source("scripts/utils_cea.R")
source("scripts/utils_schema.R")

# ==============================================================================
# Configuration & IO
# ==============================================================================
INPUT_FILE <- "outputs/results/supplement/phenotype_stratified/base_case/phenotype_iteration_outcomes.csv"
OUT_DIR <- "outputs/results/supplement/cystic_subgroup"

SUMMARY_TABLE <- file.path(OUT_DIR, "cystic_outcomes_by_strategy.csv")
COMP_TABLE <- file.path(OUT_DIR, "cystic_cost_components_by_strategy.csv")
COMP_FIG_BASE <- file.path(OUT_DIR, "mean_cost_composition_by_strategy_cystic")
INC_TABLE <- file.path(OUT_DIR, "incremental_table_cystic.csv")
FRONTIER_FIG_BASE <- file.path(OUT_DIR, "efficiency_frontier_cystic")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ==============================================================================
# 1. Load and Filter Data
# ==============================================================================
cat("Loading phenotype outcomes from:", INPUT_FILE, "\n")
pheno_data <- read_csv(INPUT_FILE, show_col_types = FALSE)

required_cols <- c(
    "iteration_id", "phenotype_category", "strategy_id", "strategy_label",
    "total_cost_per_proband_cad", "diagnoses_per_proband", "cost_per_diagnosis_cad",
    "cost_consultation_cad", "cost_testing_cad", "cost_cascade_cad", "cost_vus_followup_cad"
)
assert_required_columns(pheno_data, required_cols, "phenotype_iteration_outcomes")

cystic_data <- pheno_data %>% filter(phenotype_category == "cystic")
if (nrow(cystic_data) == 0) stop("No cystic rows found in phenotype iteration outcomes")

# ==============================================================================
# 2. Summary Statistics (schema aligned with §7.3)
# ==============================================================================
summary_stats <- cystic_data %>%
    group_by(strategy_id, strategy_label) %>%
    summarise(
        total_cost_per_proband_cad_mean = mean(total_cost_per_proband_cad, na.rm = TRUE),
        total_cost_per_proband_cad_ui_low = quantile(total_cost_per_proband_cad, 0.025, na.rm = TRUE),
        total_cost_per_proband_cad_ui_high = quantile(total_cost_per_proband_cad, 0.975, na.rm = TRUE),
        diagnoses_per_proband_mean = mean(diagnoses_per_proband, na.rm = TRUE),
        diagnoses_per_proband_ui_low = quantile(diagnoses_per_proband, 0.025, na.rm = TRUE),
        diagnoses_per_proband_ui_high = quantile(diagnoses_per_proband, 0.975, na.rm = TRUE),
        cost_per_diagnosis_cad_mean = mean(cost_per_diagnosis_cad, na.rm = TRUE),
        cost_per_diagnosis_cad_ui_low = quantile(cost_per_diagnosis_cad, 0.025, na.rm = TRUE),
        cost_per_diagnosis_cad_ui_high = quantile(cost_per_diagnosis_cad, 0.975, na.rm = TRUE),
        pr_at_least_one_vus_mean = NA_real_,
        pr_at_least_one_vus_ui_low = NA_real_,
        pr_at_least_one_vus_ui_high = NA_real_,
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
assert_required_columns(summary_stats, summary_required, "cystic_outcomes_by_strategy", exact = TRUE)
assert_no_na(summary_stats, c("strategy_id", "strategy_label"), "cystic_outcomes_by_strategy")
write_csv_validated(summary_stats, SUMMARY_TABLE, "cystic_outcomes_by_strategy")

# ==============================================================================
# 3. Cost Composition (schema aligned with §7.4 for subgroup naming)
# ==============================================================================
comp_wide <- cystic_data %>%
    group_by(strategy_id, strategy_label) %>%
    summarise(
        cost_consultation_mean = mean(cost_consultation_cad, na.rm = TRUE),
        cost_testing_mean = mean(cost_testing_cad, na.rm = TRUE),
        cost_cascade_mean = mean(cost_cascade_cad, na.rm = TRUE),
        cost_vus_followup_mean = mean(cost_vus_followup_cad, na.rm = TRUE),
        total_cost_per_proband_cad_mean = mean(total_cost_per_proband_cad, na.rm = TRUE),
        total_cost_per_proband_cad_ui_low = quantile(total_cost_per_proband_cad, 0.025, na.rm = TRUE),
        total_cost_per_proband_cad_ui_high = quantile(total_cost_per_proband_cad, 0.975, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(total_cost_per_proband_cad_mean)

component_label_map <- c(
    cost_consultation_mean = "Genetics Visits",
    cost_testing_mean = "Laboratory Testing",
    cost_cascade_mean = "Cascade Testing",
    cost_vus_followup_mean = "VUS Follow-up"
)

comp_long <- comp_wide %>%
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

write_csv_validated(comp_long, COMP_TABLE, "cystic_cost_components_by_strategy")

plot_data <- comp_long %>%
    mutate(
        strategy_label = factor(strategy_label, levels = comp_wide$strategy_label),
        component_name = factor(component_name,
            levels = c("Genetics Visits", "Laboratory Testing", "Cascade Testing", "VUS Follow-up")
        )
    )

p_comp <- ggplot(plot_data, aes(x = strategy_label, y = component_cost_per_proband_cad_mean, fill = component_name)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    geom_errorbar(
        data = comp_wide,
        aes(
            x = strategy_label,
            y = total_cost_per_proband_cad_mean,
            ymin = total_cost_per_proband_cad_ui_low,
            ymax = total_cost_per_proband_cad_ui_high
        ),
        inherit.aes = FALSE,
        width = 0.2,
        linewidth = 0.4
    ) +
    scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.9, direction = 1) +
    scale_y_continuous(labels = scales::dollar_format(prefix = "$", suffix = "")) +
    labs(
        title = "Mean Cost Composition (Cystic Subgroup)",
        x = NULL,
        y = "Total Cost per Proband (CAD)",
        fill = "Cost Component"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        plot.background = element_rect(fill = "white", color = NA)
    )

ggsave(paste0(COMP_FIG_BASE, ".png"), p_comp, width = 8, height = 6, dpi = 300)
ggsave(paste0(COMP_FIG_BASE, ".svg"), p_comp, width = 8, height = 6)

# ==============================================================================
# 4. Incremental Analysis & Frontier (§7.5 schema)
# ==============================================================================
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

results <- perform_incremental_analysis(analysis_data)

inc_table <- results %>%
    arrange(cost) %>%
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
        dominance_status
    )

inc_required <- c(
    "strategy_id", "strategy_label",
    "total_cost_per_proband_cad_mean", "total_cost_per_proband_cad_ui_low", "total_cost_per_proband_cad_ui_high",
    "diagnoses_per_proband_mean", "diagnoses_per_proband_ui_low", "diagnoses_per_proband_ui_high",
    "incremental_cost_vs_prev_cad_mean", "incremental_diagnoses_vs_prev_mean",
    "icpd_cad_per_additional_diagnosis_mean", "dominance_status"
)
assert_required_columns(inc_table, inc_required, "incremental_table_cystic", exact = TRUE)
assert_no_na(inc_table, c("strategy_id", "strategy_label", "dominance_status"), "incremental_table_cystic")
assert_frontier_first_incremental_na(inc_table, table_name = "incremental_table_cystic")
write_csv_validated(inc_table, INC_TABLE, "incremental_table_cystic")

plot_df <- results %>%
    mutate(
        is_frontier = dominance_status == "non_dominated",
        alpha_val = ifelse(is_frontier, 1, 0.4),
        size_val = ifelse(is_frontier, 3, 2)
    )

frontier_line <- results %>%
    filter(dominance_status == "non_dominated") %>%
    arrange(effect)

p_frontier <- ggplot(plot_df, aes(x = effect, y = cost)) +
    geom_errorbar(
        aes(ymin = cost_ui_low, ymax = cost_ui_high, alpha = alpha_val),
        width = 0.005, color = "gray60", linewidth = 0.4
    ) +
    geom_errorbarh(
        aes(xmin = effect_ui_low, xmax = effect_ui_high, alpha = alpha_val),
        height = 50, color = "gray60", linewidth = 0.4
    ) +
    geom_line(
        data = frontier_line,
        aes(x = effect, y = cost),
        color = "black", linewidth = 0.8, linetype = "solid"
    ) +
    geom_point(
        aes(color = strategy_label, shape = strategy_label, alpha = alpha_val, size = size_val)
    ) +
    geom_text_repel(
        aes(label = strategy_label, alpha = alpha_val),
        size = 3.5, family = "sans",
        box.padding = 0.7, point.padding = 0.5,
        min.segment.length = 0.1, segment.size = 0.3,
        max.overlaps = 20,
        force = 2,
        force_pull = 0.5
    ) +
    scale_color_viridis_d(option = "D", begin = 0.15, end = 0.85, direction = 1) +
    scale_shape_manual(values = c(15, 16, 17, 18)) +
    scale_alpha_identity() +
    scale_size_identity() +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = 0.15)) +
    scale_y_continuous(labels = scales::dollar_format(prefix = "$", suffix = ""), expand = expansion(mult = c(0.05, 0.1))) +
    labs(
        title = "Efficiency Frontier: Cystic Subgroup",
        subtitle = "Error bars show 95% uncertainty intervals",
        x = "Diagnostic Yield (Proportion of Probands Diagnosed)",
        y = "Total Cost per Proband (CAD)",
        color = NULL,
        shape = NULL
    ) +
    theme_minimal(base_size = 11, base_family = "sans") +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11, face = "bold"),
        plot.title = element_text(size = 13, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 9, color = "gray30", hjust = 0),
        legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
    )

ggsave(paste0(FRONTIER_FIG_BASE, ".png"), p_frontier, width = 7, height = 5.5, dpi = 300, bg = "white")
ggsave(paste0(FRONTIER_FIG_BASE, ".svg"), p_frontier, width = 7, height = 5.5, bg = "white")

cat("Done. Cystic subgroup outputs generated in", OUT_DIR, "\n")
