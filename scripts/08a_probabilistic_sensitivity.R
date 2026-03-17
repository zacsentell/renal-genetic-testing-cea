# 08a_probabilistic_sensitivity.R
# Purpose: Probabilistic Sensitivity Analysis (§8.1)
# Author: Renal Genetics CEA Team
# Date: 2025-12-30

library(dplyr)
library(readr)
library(ggplot2)

# Source shared utility for incremental analysis logic
source("scripts/utils_cea.R")

# ==============================================================================
# Configuration & IO
# ==============================================================================
INPUT_FILE <- "outputs/results/base_case/iteration_level/strategy_iteration_outcomes.csv"
OUT_DIR <- "outputs/results/uncertainty_sensitivity/probabilistic_sensitivity_analysis"

# Output files
ROBUSTNESS_TABLE <- file.path(OUT_DIR, "psa_decision_robustness.csv")
SCATTER_FIG_BASE <- file.path(OUT_DIR, "psa_cost_diagnosis_scatter")

# Ensure output directory exists
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ==============================================================================
# 1. Load Data
# ==============================================================================
cat("Loading iteration-level outcomes from:", INPUT_FILE, "\n")
iter_data <- read_csv(INPUT_FILE, show_col_types = FALSE)

# Check required columns
required_cols <- c(
    "iteration_id", "strategy_id", "strategy_label",
    "total_cost_per_proband_cad", "diagnoses_per_proband"
)
stopifnot("Missing required columns" = all(required_cols %in% names(iter_data)))

# ==============================================================================
# 2. Decision Robustness (Dominance per Iteration)
# ==============================================================================
cat("\n--- Computing Decision Robustness ---\n")

# Get unique iterations
iteration_ids <- unique(iter_data$iteration_id)
n_iter <- length(iteration_ids)
cat("Analyzing dominance across", n_iter, "iterations...\n")

# Initialize storage for dominance results
dominance_results <- list()

# Loop through iterations (vectorization is hard due to sorting/dominance logic)
# Could use parallel, but 1000 iter x 4 strategies is fast enough.
for (i in iteration_ids) {
    # Subset for this iteration
    df_iter <- iter_data %>%
        filter(iteration_id == i) %>%
        select(
            strategy_id, strategy_label,
            cost = total_cost_per_proband_cad,
            effect = diagnoses_per_proband
        )

    # Calculate dominance using shared utility
    res <- perform_incremental_analysis(df_iter)

    # Store relevant info
    dominance_results[[i]] <- data.frame(
        iteration_id = i,
        strategy_label = res$strategy_label,
        dominance_status = res$dominance_status,
        stringsAsFactors = FALSE
    )

    if (i %% 100 == 0) cat(".")
}
cat("\n")

# Combine results
dominance_df <- do.call(rbind, dominance_results)

# Calculate Probability Non-Dominated
robustness_summary <- dominance_df %>%
    group_by(strategy_label) %>%
    summarise(
        n_non_dominated = sum(dominance_status == "non_dominated"),
        n_total = n(),
        pr_non_dominated = n_non_dominated / n_total,
        .groups = "drop"
    ) %>%
    arrange(desc(pr_non_dominated))

# Determine strategy order for display (match base case cost order if possible, or just alphabetically)
# Let's get strategy IDs to order by ID
strat_map <- iter_data %>%
    select(strategy_id, strategy_label) %>%
    distinct() %>%
    arrange(strategy_id)
robustness_summary <- robustness_summary %>%
    left_join(strat_map, by = "strategy_label") %>%
    select(strategy_id, strategy_label, pr_non_dominated) %>%
    arrange(strategy_id)

# Export Table
write_csv(robustness_summary, ROBUSTNESS_TABLE)
cat("Decision robustness table exported to:", ROBUSTNESS_TABLE, "\n")
print(as.data.frame(robustness_summary))


# ==============================================================================
# 3. PSA Scatter Plot
# ==============================================================================
cat("\n--- Generating PSA Scatter Plot ---\n")

# Calculate mean points for overlay
mean_points <- iter_data %>%
    group_by(strategy_label) %>%
    summarise(
        mean_cost = mean(total_cost_per_proband_cad),
        mean_effect = mean(diagnoses_per_proband),
        .groups = "drop"
    )

# Clean labels for plot
iter_data <- iter_data %>%
    mutate(
        strategy_display = case_when(
            strategy_label == "Panel_Reflex_ES" ~ "Panel reflex to ES",
            TRUE ~ strategy_label
        )
    )

mean_points <- mean_points %>%
    mutate(
        strategy_display = case_when(
            strategy_label == "Panel_Reflex_ES" ~ "Panel reflex to ES",
            TRUE ~ strategy_label
        )
    )

p_scatter <- ggplot(iter_data, aes(x = diagnoses_per_proband, y = total_cost_per_proband_cad)) +
    # Scatter points (semi-transparent)
    geom_point(aes(color = strategy_display), alpha = 0.15, size = 1.5) +
    # Mean centroids (larger, solid)
    geom_point(
        data = mean_points,
        aes(x = mean_effect, y = mean_cost, fill = strategy_display),
        shape = 21, size = 5, color = "black", stroke = 1
    ) +
    # Scales
    scale_color_viridis_d(option = "D", begin = 0.15, end = 0.85) +
    scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.85) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::dollar_format(prefix = "$", suffix = "")) +
    # Labels
    labs(
        title = "Probabilistic Sensitivity Analysis: Cost vs Diagnostic Yield",
        subtitle = "1,000 Monte Carlo iterations. Large points represent mean outcomes.",
        x = "Diagnostic Yield (Proportion of Probands Diagnosed)",
        y = "Total Cost per Proband (CAD)",
        color = "Strategy",
        fill = "Strategy"
    ) +
    # Theme
    theme_minimal(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", color = NA)
    )

ggsave(paste0(SCATTER_FIG_BASE, ".png"), p_scatter, width = 8, height = 6, dpi = 300)
ggsave(paste0(SCATTER_FIG_BASE, ".svg"), p_scatter, width = 8, height = 6)
cat("PSA scatter plot exported.\n")

cat("\nDone. PSA analysis complete.\n")
