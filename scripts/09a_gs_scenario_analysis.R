# 09a_gs_scenario_analysis.R
# Purpose: Scenario analysis for Genome Sequencing with enhanced detection (§8.5)

library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
if (!requireNamespace("config", quietly = TRUE)) {
    stop("Package 'config' is required. Install it with install.packages('config').")
}
source("scripts/utils_schema.R")
source("scripts/utils_cea.R")

# ==============================================================================
# 1. Configuration
# ==============================================================================
cfg <- config::get()
GS_YIELD_UPLIFT <- cfg$scenario$gs_yield_uplift

INPUT_BASE_ITERATION <- "outputs/results/base_case/iteration_level/strategy_iteration_outcomes.csv"
INPUT_GS_BASELINE_ITERATION <- "outputs/results/scenario_analysis/gs_uplift/gs_baseline_iteration_outcomes.csv"
OUTPUT_DIR <- "outputs/results/scenario_analysis/gs_uplift"

OUTPUT_SCENARIO_ITER <- file.path(OUTPUT_DIR, "gs_scenario_iteration_outcomes.csv")
OUTPUT_SUMMARY <- file.path(OUTPUT_DIR, "gs_scenario_summary.csv")
OUTPUT_INCREMENTAL <- file.path(OUTPUT_DIR, "gs_scenario_incremental.csv")
OUTPUT_FIG_BASE <- file.path(OUTPUT_DIR, "gs_scenario_comparison")
OUTPUT_FRONTIER_BASE <- file.path(OUTPUT_DIR, "efficiency_frontier_gs_scenario")
OUTPUT_CEAC_TABLE <- file.path(OUTPUT_DIR, "ceac_table_gs_scenario.csv")
OUTPUT_CEAC_FIG <- file.path(OUTPUT_DIR, "ceac_plot_gs_scenario")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ==============================================================================
# 2. Load Inputs
# ==============================================================================
base_iter <- read_csv(INPUT_BASE_ITERATION, show_col_types = FALSE)
gs_baseline <- read_csv(INPUT_GS_BASELINE_ITERATION, show_col_types = FALSE)

base_required <- c(
    "iteration_id", "strategy_id", "strategy_label",
    "total_cost_per_proband_cad", "diagnoses_per_proband", "cost_per_diagnosis_cad", "pr_at_least_one_vus"
)
gs_required <- c(
    "iteration_id", "strategy_id", "strategy_label", "scenario_label",
    "total_cost_per_proband_cad", "diagnoses_per_proband", "cost_per_diagnosis_cad", "pr_at_least_one_vus"
)
assert_required_columns(base_iter, base_required, "base_strategy_iteration_outcomes")
assert_required_columns(gs_baseline, gs_required, "gs_baseline_iteration_outcomes", exact = TRUE)

# ==============================================================================
# 3. Build Scenario Iteration Outcomes
# ==============================================================================
gs_scenario <- gs_baseline %>%
    filter(scenario_label == "baseline") %>%
    mutate(
        scenario_label = sprintf("uplift_%spct", round(GS_YIELD_UPLIFT * 100)),
        diagnoses_per_proband = diagnoses_per_proband * (1 + GS_YIELD_UPLIFT),
        cost_per_diagnosis_cad = ifelse(
            diagnoses_per_proband > 0,
            total_cost_per_proband_cad / diagnoses_per_proband,
            NA_real_
        )
    )

scenario_iterations <- gs_scenario

assert_required_columns(scenario_iterations, gs_required, "gs_scenario_iteration_outcomes", exact = TRUE)
assert_no_na(
    scenario_iterations,
    c("iteration_id", "strategy_id", "strategy_label", "scenario_label", "total_cost_per_proband_cad", "diagnoses_per_proband"),
    "gs_scenario_iteration_outcomes"
)
write_csv_validated(scenario_iterations, OUTPUT_SCENARIO_ITER, "gs_scenario_iteration_outcomes")

# ==============================================================================
# 4. Summaries and Incremental vs Reflex
# ==============================================================================
reflex_iter <- base_iter %>%
    filter(strategy_label == "Panel_Reflex_ES")

gs_scenario_summary <- gs_scenario %>%
    summarise(
        scenario = sprintf("GS Scenario (+%s%% yield)", round(GS_YIELD_UPLIFT * 100)),
        strategy = "GS",
        yield_mean = mean(diagnoses_per_proband),
        cost_mean = mean(total_cost_per_proband_cad),
        cost_per_dx_mean = mean(cost_per_diagnosis_cad, na.rm = TRUE)
    )

reflex_summary <- reflex_iter %>%
    summarise(
        scenario = "Base Case",
        strategy = "Panel_Reflex_ES",
        yield_mean = mean(diagnoses_per_proband),
        cost_mean = mean(total_cost_per_proband_cad),
        cost_per_dx_mean = mean(cost_per_diagnosis_cad, na.rm = TRUE)
    )

summary_table <- bind_rows(reflex_summary, gs_scenario_summary)
write_csv_validated(summary_table, OUTPUT_SUMMARY, "gs_scenario_summary")

incremental <- data.frame(
    comparison = "GS_Scenario_vs_Reflex",
    gs_yield = gs_scenario_summary$yield_mean,
    reflex_yield = reflex_summary$yield_mean,
    delta_yield = gs_scenario_summary$yield_mean - reflex_summary$yield_mean,
    gs_cost = gs_scenario_summary$cost_mean,
    reflex_cost = reflex_summary$cost_mean,
    delta_cost = gs_scenario_summary$cost_mean - reflex_summary$cost_mean,
    stringsAsFactors = FALSE
)
incremental$icer <- ifelse(incremental$delta_yield > 0, incremental$delta_cost / incremental$delta_yield, NA_real_)
write_csv_validated(incremental, OUTPUT_INCREMENTAL, "gs_scenario_incremental")

# ==============================================================================
# 5. Scenario Comparison Figure
# ==============================================================================
panel_summary <- base_iter %>%
    filter(strategy_label == "Panel") %>%
    summarise(strategy = "Panel", yield = mean(diagnoses_per_proband), cost = mean(total_cost_per_proband_cad))

es_summary <- base_iter %>%
    filter(strategy_label == "ES") %>%
    summarise(strategy = "ES", yield = mean(diagnoses_per_proband), cost = mean(total_cost_per_proband_cad))

plot_data <- bind_rows(
    panel_summary,
    es_summary,
    data.frame(strategy = "Reflex", yield = reflex_summary$yield_mean, cost = reflex_summary$cost_mean),
    data.frame(strategy = "GS (+10% yield)", yield = gs_scenario_summary$yield_mean, cost = gs_scenario_summary$cost_mean)
) %>%
    mutate(
        strategy_type = case_when(
            strategy == "GS (+10% yield)" ~ "Scenario",
            TRUE ~ "Base Case"
        )
    )

p <- ggplot(plot_data, aes(x = yield, y = cost, color = strategy_type)) +
    geom_point(aes(shape = strategy_type), size = 4, alpha = 0.9) +
    geom_label_repel(aes(label = strategy), size = 4, fontface = "bold", box.padding = 0.8, point.padding = 0.4) +
    scale_color_manual(values = c("Base Case" = "#2166AC", "Scenario" = "#4DAF4A")) +
    scale_shape_manual(values = c("Base Case" = 16, "Scenario" = 17)) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::dollar_format(prefix = "$")) +
    labs(
        title = "Scenario Analysis: GS with Enhanced Detection",
        subtitle = sprintf("%s%% relative GS yield uplift", round(GS_YIELD_UPLIFT * 100)),
        x = "Diagnostic Yield",
        y = "Total Cost per Proband (CAD)",
        color = NULL,
        shape = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title = element_text(face = "bold", size = 14))

ggsave(paste0(OUTPUT_FIG_BASE, ".png"), p, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(paste0(OUTPUT_FIG_BASE, ".svg"), p, width = 8, height = 6, bg = "white")

# ==============================================================================
# 6. Efficiency Frontier Plot with Scenario GS
# ==============================================================================
vus_panel <- mean(base_iter$pr_at_least_one_vus[base_iter$strategy_label == "Panel"], na.rm = TRUE)
vus_es <- mean(base_iter$pr_at_least_one_vus[base_iter$strategy_label == "ES"], na.rm = TRUE)
vus_reflex <- mean(base_iter$pr_at_least_one_vus[base_iter$strategy_label == "Panel_Reflex_ES"], na.rm = TRUE)
vus_gs <- mean(gs_scenario$pr_at_least_one_vus, na.rm = TRUE)

frontier_raw <- bind_rows(
    data.frame(strategy_label = "Panel", effect = panel_summary$yield, cost = panel_summary$cost, vus_prob = vus_panel),
    data.frame(strategy_label = "ES", effect = es_summary$yield, cost = es_summary$cost, vus_prob = vus_es),
    data.frame(strategy_label = "Panel reflex to ES", effect = reflex_summary$yield_mean, cost = reflex_summary$cost_mean, vus_prob = vus_reflex),
    data.frame(strategy_label = "GS (+10% yield)", effect = gs_scenario_summary$yield_mean, cost = gs_scenario_summary$cost_mean, vus_prob = vus_gs)
)

frontier_dominance <- perform_incremental_analysis(
    frontier_raw %>% select(strategy_label, cost, effect)
) %>%
    select(strategy_label, dominance_status)

frontier_plot_data <- frontier_raw %>%
    left_join(frontier_dominance, by = "strategy_label") %>%
    rename(strategy_display = strategy_label) %>%
    mutate(
        strategy_label_full = paste0(strategy_display, "\n(VUS: ", round(vus_prob * 100, 0), "%)"),
        is_dominated = dominance_status != "non_dominated",
        alpha_val = ifelse(is_dominated, 0.5, 1.0),
        bubble_size = 3 + (vus_prob * 8)
    )

frontier_line_data <- frontier_plot_data %>%
    filter(dominance_status == "non_dominated") %>%
    arrange(effect)

p_frontier <- ggplot(frontier_plot_data, aes(x = effect, y = cost)) +
    geom_line(data = frontier_line_data, aes(x = effect, y = cost), color = "black", linewidth = 1.0) +
    geom_point(aes(color = strategy_display, alpha = alpha_val, size = bubble_size)) +
    geom_label_repel(
        aes(label = strategy_label_full, alpha = alpha_val),
        size = 4.5, family = "sans", fontface = "bold",
        box.padding = 1.0, point.padding = 0.6,
        min.segment.length = 0.1, segment.size = 0.4,
        segment.color = "gray40", max.overlaps = 20,
        force = 3, force_pull = 0.4,
        fill = "white", label.size = 0.2,
        label.padding = unit(0.25, "lines")
    ) +
    scale_color_viridis_d(option = "D", begin = 0.15, end = 0.85, direction = 1) +
    scale_alpha_identity() +
    scale_size_identity() +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = 0.18)) +
    scale_y_continuous(labels = scales::dollar_format(prefix = "$", suffix = ""), expand = expansion(mult = c(0.08, 0.12))) +
    labs(
        title = "Efficiency Frontier: Cost vs Diagnostic Yield",
        subtitle = "Includes GS scenario with +10% diagnostic yield uplift.",
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

ggsave(paste0(OUTPUT_FRONTIER_BASE, ".png"), p_frontier, width = 7, height = 5.5, dpi = 300, bg = "white")
ggsave(paste0(OUTPUT_FRONTIER_BASE, ".svg"), p_frontier, width = 7, height = 5.5, bg = "white")

# ==============================================================================
# 7. CEAC including GS Scenario
# ==============================================================================
wtp_min <- cfg$wtp$min
wtp_max <- cfg$wtp$max
wtp_step <- cfg$wtp$step
wtp_range <- seq(wtp_min, wtp_max, by = wtp_step)

base_for_ceac <- base_iter %>%
    select(iteration_id, strategy_id, strategy_label, total_cost_per_proband_cad, diagnoses_per_proband)

gs_for_ceac <- gs_scenario %>%
    transmute(
        iteration_id,
        strategy_id,
        strategy_label = "GS_scenario",
        total_cost_per_proband_cad,
        diagnoses_per_proband
    )

ceac_input <- bind_rows(base_for_ceac, gs_for_ceac)
ceac_final <- compute_ceac(ceac_input, wtp_range)

ceac_required <- c("wtp_threshold_cad", "strategy_id", "strategy_label", "pr_cost_effective")
assert_required_columns(ceac_final, ceac_required, "ceac_table_gs_scenario", exact = TRUE)
assert_no_na(ceac_final, ceac_required, "ceac_table_gs_scenario")
write_csv_validated(ceac_final, OUTPUT_CEAC_TABLE, "ceac_table_gs_scenario")

ceac_plot <- ceac_final %>%
    mutate(
        strategy_display = case_when(
            strategy_label == "Panel_Reflex_ES" ~ "Panel reflex to ES",
            strategy_label == "GS_scenario" ~ "GS (+10% yield)",
            TRUE ~ strategy_label
        )
    )

p_ceac <- ggplot(ceac_plot, aes(x = wtp_threshold_cad, y = pr_cost_effective, color = strategy_display)) +
    geom_line(linewidth = 1.0) +
    scale_color_viridis_d(option = "D", begin = 0.15, end = 0.85, direction = 1) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0)) +
    scale_x_continuous(labels = scales::dollar_format(prefix = "$", suffix = ""), expand = c(0, 0)) +
    labs(
        title = "Cost-Effectiveness Acceptability Curve (CEAC)",
        subtitle = "Base-case strategies plus GS (+10% yield) scenario",
        x = "Willingness-to-Pay per Additional Diagnosis (CAD)",
        y = "Probability Cost-Effective",
        color = NULL
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
        plot.subtitle = element_text(size = 10, color = "gray30", hjust = 0),
        legend.position = "right",
        legend.text = element_text(size = 10),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
    )

ggsave(paste0(OUTPUT_CEAC_FIG, ".png"), p_ceac, width = 7, height = 5, dpi = 300, bg = "white")
ggsave(paste0(OUTPUT_CEAC_FIG, ".svg"), p_ceac, width = 7, height = 5, bg = "white")

# ==============================================================================
# 8. GS Uplift Threshold Analysis
# ==============================================================================
# Sweep uplift from 0% to 20% in 0.5% steps. At each uplift, compute mean GS
# yield and cost, run dominance analysis, and record whether GS is on the
# efficiency frontier and the ICPD vs Reflex.

cat("\n--- GS Uplift Threshold Analysis ---\n")

OUTPUT_THRESHOLD <- file.path(OUTPUT_DIR, "gs_uplift_threshold.csv")

gs_baseline_iter <- gs_baseline %>% filter(scenario_label == "baseline")
reflex_mean_yield <- mean(reflex_iter$diagnoses_per_proband, na.rm = TRUE)
reflex_mean_cost  <- mean(reflex_iter$total_cost_per_proband_cad, na.rm = TRUE)
gs_base_yield     <- mean(gs_baseline_iter$diagnoses_per_proband, na.rm = TRUE)
gs_mean_cost      <- mean(gs_baseline_iter$total_cost_per_proband_cad, na.rm = TRUE)

uplift_seq <- seq(0, 0.20, by = 0.005)

threshold_results <- lapply(uplift_seq, function(u) {
    gs_yield_u <- gs_base_yield * (1 + u)

    frontier_input <- data.frame(
        strategy_label = c("Panel", "ES", "Reflex", "GS"),
        cost   = c(panel_summary$cost, es_summary$cost, reflex_mean_cost, gs_mean_cost),
        effect = c(panel_summary$yield, es_summary$yield, reflex_mean_yield, gs_yield_u),
        stringsAsFactors = FALSE
    )

    dom <- perform_incremental_analysis(frontier_input)
    gs_row <- dom[dom$strategy_label == "GS", , drop = FALSE]

    data.frame(
        uplift_pct       = u * 100,
        gs_yield         = gs_yield_u,
        reflex_yield     = reflex_mean_yield,
        delta_yield      = gs_yield_u - reflex_mean_yield,
        gs_cost          = gs_mean_cost,
        reflex_cost      = reflex_mean_cost,
        delta_cost       = gs_mean_cost - reflex_mean_cost,
        gs_on_frontier   = gs_row$dominance_status == "non_dominated",
        icpd_vs_reflex   = ifelse(
            gs_row$dominance_status == "non_dominated" & !is.na(gs_row$icer),
            gs_row$icer, NA_real_
        ),
        stringsAsFactors = FALSE
    )
})

threshold_df <- do.call(rbind, threshold_results)
write_csv_validated(threshold_df, OUTPUT_THRESHOLD, "gs_uplift_threshold")

# Identify minimum uplift where GS enters the frontier
min_frontier_row <- threshold_df[threshold_df$gs_on_frontier, , drop = FALSE]
if (nrow(min_frontier_row) > 0) {
    min_row <- min_frontier_row[1, , drop = FALSE]
    cat(sprintf(
        "Minimum uplift for GS on frontier: +%.1f%% (GS yield %.1f%% vs Reflex %.1f%%; ICPD $%s)\n",
        min_row$uplift_pct,
        min_row$gs_yield * 100,
        min_row$reflex_yield * 100,
        formatC(round(min_row$icpd_vs_reflex), format = "f", digits = 0, big.mark = ",")
    ))
} else {
    cat("GS does not enter the frontier at any tested uplift (0-20%).\n")
}

cat("GS scenario analysis complete. Outputs written to", OUTPUT_DIR, "\n")
