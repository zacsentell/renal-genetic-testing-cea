# 08_wtp_analysis.R
# Purpose: Willingness-to-Pay (WTP) Analysis & CEAC (§8.4)

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(scales)
if (!requireNamespace("config", quietly = TRUE)) {
    stop("Package 'config' is required. Install it with install.packages('config').")
}
source("scripts/utils_schema.R")
source("scripts/utils_cea.R")

# ==============================================================================
# 1. Configuration
# ==============================================================================
cfg <- config::get()
INPUT_FILE <- "outputs/results/base_case/iteration_level/strategy_iteration_outcomes.csv"
OUTPUT_DIR <- "outputs/results/uncertainty_sensitivity/willingness_to_pay"
OUTPUT_TABLE <- file.path(OUTPUT_DIR, "ceac_table.csv")
OUTPUT_FIG_BASE <- file.path(OUTPUT_DIR, "ceac_plot")

WTP_MIN <- cfg$wtp$min
WTP_MAX <- cfg$wtp$max
WTP_STEP <- cfg$wtp$step
WTP_RANGE <- seq(WTP_MIN, WTP_MAX, by = WTP_STEP)

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ==============================================================================
# 2. Load Data
# ==============================================================================
outcomes <- read_csv(INPUT_FILE, show_col_types = FALSE)
required_cols <- c(
    "iteration_id", "strategy_id", "strategy_label",
    "total_cost_per_proband_cad", "diagnoses_per_proband"
)
assert_required_columns(outcomes, required_cols, "strategy_iteration_outcomes")
assert_no_na(outcomes, c("iteration_id", "strategy_id", "strategy_label"), "strategy_iteration_outcomes")

cat("Loaded", nrow(outcomes), "rows from", INPUT_FILE, "\n")

# ==============================================================================
# 3. CEAC Calculation
# ==============================================================================
final_ceac <- compute_ceac(outcomes, WTP_RANGE)

ceac_required <- c("wtp_threshold_cad", "strategy_id", "strategy_label", "pr_cost_effective")
assert_required_columns(final_ceac, ceac_required, "ceac_table", exact = TRUE)
assert_no_na(final_ceac, ceac_required, "ceac_table")

write_csv_validated(final_ceac, OUTPUT_TABLE, "ceac_table")

# ==============================================================================
# 4. Figure
# ==============================================================================
plot_data <- final_ceac %>%
    mutate(
        strategy_display = case_when(
            strategy_label == "Panel_Reflex_ES" ~ "Panel reflex to ES",
            TRUE ~ strategy_label
        )
    )

p <- ggplot(plot_data, aes(x = wtp_threshold_cad, y = pr_cost_effective, color = strategy_display)) +
    geom_line(linewidth = 1.0) +
    scale_color_viridis_d(option = "D", begin = 0.15, end = 0.85, direction = 1) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0)) +
    scale_x_continuous(labels = scales::dollar_format(prefix = "$", suffix = ""), expand = c(0, 0)) +
    labs(
        title = "Cost-Effectiveness Acceptability Curve (CEAC)",
        subtitle = "Probability that a strategy is cost-effective at each willingness-to-pay threshold",
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

ggsave(paste0(OUTPUT_FIG_BASE, ".png"), p, width = 7, height = 5, dpi = 300, bg = "white")
ggsave(paste0(OUTPUT_FIG_BASE, ".svg"), p, width = 7, height = 5, bg = "white")

cat("Exported:\n")
cat(" -", OUTPUT_TABLE, "\n")
cat(" -", paste0(OUTPUT_FIG_BASE, ".png"), "\n")
cat(" -", paste0(OUTPUT_FIG_BASE, ".svg"), "\n")
