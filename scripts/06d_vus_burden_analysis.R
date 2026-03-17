# 06d_vus_burden_analysis.R
# Purpose: Analyze VUS burden and trade-offs (§7.7)
# Author: Renal Genetics CEA Team
# Date: 2025-01-03

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(scales)

# ==============================================================================
# 1. Configuration
# ==============================================================================
INPUT_FILE <- "outputs/results/base_case/iteration_level/strategy_iteration_outcomes.csv"
OUTPUT_DIR <- "outputs/results/base_case/vus_burden"
OUTPUT_TABLE <- file.path(OUTPUT_DIR, "outcome_proportions_by_strategy.csv")
FIG_TRADEOFF <- file.path(OUTPUT_DIR, "diagnosis_vus_tradeoff")
FIG_STACKED <- file.path(OUTPUT_DIR, "outcome_proportions_stacked")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ==============================================================================
# 2. Load Data
# ==============================================================================
cat("Loading iteration outcomes from:", INPUT_FILE, "\n")
iter_data <- read_csv(INPUT_FILE, show_col_types = FALSE)

# Check required columns
req_cols <- c(
    "strategy_id", "strategy_label", "iteration_id",
    "diagnoses_per_proband", "pr_at_least_one_vus", "total_cost_per_proband_cad"
)
stopifnot(all(req_cols %in% names(iter_data)))

# ==============================================================================
# 3. Compute Outcome Proportions (Iteration Level)
# ==============================================================================
# Approximation: pr_vus_only =~ pr_at_least_one_vus * (1 - pr_diagnosed)
# Assumption: VUS occurrence is independent of Diagnosis status (approximation)

iter_data <- iter_data %>%
    mutate(
        pr_diagnosed = diagnoses_per_proband,
        pr_vus_only = pr_at_least_one_vus * (1 - diagnoses_per_proband),
        pr_negative = 1 - pr_diagnosed - pr_vus_only
    )

# Verify summation to 1 (floating point tolerance)
check_sums <- rowSums(iter_data %>% select(pr_diagnosed, pr_vus_only, pr_negative))
if (any(abs(check_sums - 1) > 1e-6)) {
    warning(
        "Outcome proportions do not sum to 1 in some iterations (max error: ",
        max(abs(check_sums - 1)), ")"
    )
}

# ==============================================================================
# 4. Summarize Results by Strategy
# ==============================================================================
get_summary <- function(x) {
    list(
        mean = mean(x),
        ui_low = quantile(x, 0.025),
        ui_high = quantile(x, 0.975)
    )
}

summary_table <- iter_data %>%
    group_by(strategy_id, strategy_label) %>%
    summarise(
        # Cost (used for ordering)
        cost_mean = mean(total_cost_per_proband_cad),

        # Diagnosis (Yield)
        yield_mean = mean(diagnoses_per_proband),
        yield_low = quantile(diagnoses_per_proband, 0.025),
        yield_high = quantile(diagnoses_per_proband, 0.975),

        # VUS Probability (Total)
        vus_prob_mean = mean(pr_at_least_one_vus),
        vus_prob_low = quantile(pr_at_least_one_vus, 0.025),
        vus_prob_high = quantile(pr_at_least_one_vus, 0.975),

        # Outcome Proportions
        pr_diagnosed_mean = mean(pr_diagnosed),
        pr_diagnosed_ui_low = quantile(pr_diagnosed, 0.025),
        pr_diagnosed_ui_high = quantile(pr_diagnosed, 0.975),
        pr_vus_only_mean = mean(pr_vus_only),
        pr_vus_only_ui_low = quantile(pr_vus_only, 0.025),
        pr_vus_only_ui_high = quantile(pr_vus_only, 0.975),
        pr_negative_mean = mean(pr_negative),
        pr_negative_ui_low = quantile(pr_negative, 0.025),
        pr_negative_ui_high = quantile(pr_negative, 0.975),
        .groups = "drop"
    )

# Export Table
write_csv(summary_table, OUTPUT_TABLE)
cat("Summary table exported to:", OUTPUT_TABLE, "\n")

# ==============================================================================
# 5. Visualization 1: Diagnosis-VUS Tradeoff Plot
# ==============================================================================
cat("Generating Tradeoff Plot...\n")

# Clean labels for plot
plot_data_tradeoff <- summary_table %>%
    mutate(
        strategy_display = case_when(
            strategy_label == "Panel_Reflex_ES" ~ "Panel reflex to ES",
            TRUE ~ strategy_label
        )
    )

p1 <- ggplot(plot_data_tradeoff, aes(x = yield_mean, y = vus_prob_mean)) +
    # Error bars
    geom_errorbar(aes(ymin = vus_prob_low, ymax = vus_prob_high),
        width = 0.002, color = "gray70"
    ) +
    geom_errorbarh(aes(xmin = yield_low, xmax = yield_high),
        height = 0.005, color = "gray70"
    ) +
    # Points
    geom_point(aes(color = strategy_display, shape = strategy_display), size = 4) +
    # Labels
    geom_text_repel(
        aes(label = strategy_display),
        size = 4,
        box.padding = 0.5,
        min.segment.length = 0.1,
        force = 2
    ) +
    # Scales
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), name = "Diagnostic Yield (Diagnoses per Proband)") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), name = "Probability of ≥1 VUS") +
    scale_color_viridis_d(option = "D", end = 0.8) +
    scale_shape_manual(values = c(15, 16, 17, 18)) +
    # Theme
    theme_minimal(base_size = 12) +
    theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray30")
    ) +
    labs(
        title = "Diagnosis-VUS Tradeoff",
        subtitle = "Higher yield (right) vs. higher VUS burden (up). Lower-right is optimal."
    )

ggsave(paste0(FIG_TRADEOFF, ".png"), p1, width = 7, height = 5.5, dpi = 300, bg = "white")
ggsave(paste0(FIG_TRADEOFF, ".svg"), p1, width = 7, height = 5.5, bg = "white")


# ==============================================================================
# 6. Visualization 2: Stacked Outcome Proportions (Ordered by Cost)
# ==============================================================================
cat("Generating Stacked Bar Plot...\n")

# Prepare Long Data
# Ordering strategies by Cost
strat_levels <- plot_data_tradeoff %>%
    arrange(cost_mean) %>%
    pull(strategy_display) # Use display names

plot_data_stacked <- plot_data_tradeoff %>%
    select(strategy_display, pr_diagnosed_mean, pr_vus_only_mean, pr_negative_mean) %>%
    pivot_longer(
        cols = c(pr_diagnosed_mean, pr_vus_only_mean, pr_negative_mean),
        names_to = "outcome",
        values_to = "proportion"
    ) %>%
    mutate(
        strategy_display = factor(strategy_display, levels = strat_levels),
        outcome = factor(outcome,
            levels = c("pr_negative_mean", "pr_vus_only_mean", "pr_diagnosed_mean"),
            labels = c("Negative", "VUS only", "Diagnosed (P/LP)")
        )
    )

# Colors
# Diagnosed = Green (Viral/Success) -> roughly #2A9D8F or Viridis equivalent
# VUS = Amber (Warning) -> #E9C46A
# Negative = Gray -> #E0E0E0 or similar

custom_colors <- c(
    "Diagnosed (P/LP)" = "#2ca02c", # Strong Green
    "VUS only" = "#e6ab02", # Mustard/Amber
    "Negative" = "#bababa" # Gray
)

p2 <- ggplot(plot_data_stacked, aes(x = proportion, y = strategy_display, fill = outcome)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.2) +
    scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
    scale_fill_manual(values = custom_colors) +
    theme_minimal(base_size = 12) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        axis.title.y = element_blank()
    ) +
    labs(
        title = "Outcome Proportions by Strategy",
        subtitle = "Strategies ordered by increasing mean cost",
        x = "Proportion of Probands",
        fill = NULL
    )

ggsave(paste0(FIG_STACKED, ".png"), p2, width = 8, height = 5, dpi = 300, bg = "white")
ggsave(paste0(FIG_STACKED, ".svg"), p2, width = 8, height = 5, bg = "white")

cat("Analysis complete. Outputs saved to:", OUTPUT_DIR, "\n")
