# scripts/06e_phenotype_figure.R
# Purpose: Phenotype-stratified cost-effectiveness figures
#   Figure 1: Reflex cost vs yield by phenotype (canonical results figure)
#   Figure 2: Trajectory of Value – Reflex strategy across phenotypes
#   Figure 3: All-modes efficiency landscape (small multiples)

library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)

# ==============================================================================
# IO Paths
# ==============================================================================
INPUT_PHENO_CSV <- "outputs/results/supplement/phenotype_stratified/base_case/phenotype_iteration_outcomes.csv"
OUTPUT_DIR      <- "outputs/results/supplement/phenotype_stratified"

OUTPUT_FIG_REFLEX_PNG <- file.path(OUTPUT_DIR, "reflex_cost_yield_by_phenotype.png")
OUTPUT_FIG_REFLEX_SVG <- file.path(OUTPUT_DIR, "reflex_cost_yield_by_phenotype.svg")
OUTPUT_FIG_TRAJ_PNG  <- file.path(OUTPUT_DIR, "trajectory_of_value.png")
OUTPUT_FIG_TRAJ_SVG  <- file.path(OUTPUT_DIR, "trajectory_of_value.svg")
OUTPUT_FIG_ALL_PNG   <- file.path(OUTPUT_DIR, "phenotype_stratified_all_modes.png")
OUTPUT_FIG_ALL_SVG   <- file.path(OUTPUT_DIR, "phenotype_stratified_all_modes.svg")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ==============================================================================
# Load and aggregate
# ==============================================================================
cat("Loading phenotype outcomes from:", INPUT_PHENO_CSV, "\n")
df <- read.csv(INPUT_PHENO_CSV)

pheno_levels <- c("CKDu", "glomerular", "tubulointerstitial", "tubulopathies", "cystic")
pheno_labels <- c("CKDu", "Glomerular", "Tubulointerstitial", "Tubulopathies", "Cystic")

df_summary <- df %>%
    group_by(strategy_label, phenotype_category) %>%
    summarise(
        yield_mean         = mean(diagnoses_per_proband),
        cost_per_diag_mean = mean(cost_per_diagnosis_cad, na.rm = TRUE),
        cost_mean          = mean(total_cost_per_proband_cad),
        .groups = "drop"
    ) %>%
    mutate(
        phenotype_factor = factor(phenotype_category, levels = pheno_levels, labels = pheno_labels)
    )

# ==============================================================================
# Figure 1: Reflex cost vs diagnostic yield by phenotype
# ==============================================================================
df_reflex_scatter <- df_summary %>%
    filter(strategy_label == "Panel_Reflex_ES") %>%
    arrange(yield_mean)

p_reflex <- ggplot(df_reflex_scatter, aes(x = yield_mean, y = cost_mean)) +
    geom_point(aes(fill = phenotype_factor), shape = 21, size = 6, color = "black", stroke = 1.2) +
    geom_label_repel(
        aes(label = paste0(phenotype_factor, "\n", percent(yield_mean, accuracy = 0.1))),
        box.padding   = 0.7,
        point.padding = 0.5,
        size          = 4,
        fontface      = "bold",
        fill          = "white",
        label.size    = 0.2,
        min.segment.length = 0.1
    ) +
    scale_x_continuous(
        labels = percent_format(accuracy = 1),
        name   = "Diagnostic Yield (Proportion Diagnosed)",
        expand = expansion(mult = 0.15)
    ) +
    scale_y_continuous(
        labels = dollar_format(prefix = "$"),
        name   = "Total Cost per Proband (CAD)",
        expand = expansion(mult = c(0.05, 0.12))
    ) +
    scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.85, guide = "none") +
    labs(
        title    = "Reflex Testing: Cost and Yield by Phenotype",
        subtitle = "Each point represents a phenotype category under the Panel -> ES reflex strategy"
    ) +
    theme_minimal(base_size = 13, base_family = "sans") +
    theme(
        panel.grid.minor  = element_blank(),
        panel.grid.major  = element_blank(),
        axis.line         = element_line(color = "black", linewidth = 0.6),
        axis.ticks        = element_line(color = "black", linewidth = 0.5),
        axis.text         = element_text(color = "black", size = 11),
        axis.title        = element_text(color = "black", size = 12, face = "bold"),
        plot.title        = element_text(size = 14, face = "bold", hjust = 0),
        plot.subtitle     = element_text(size = 10, color = "gray30", hjust = 0),
        panel.background  = element_rect(fill = "white", color = NA),
        plot.background   = element_rect(fill = "white", color = NA),
        plot.margin       = margin(15, 15, 10, 10)
    )

cat("Saving Figure 1 (Reflex cost vs yield) to:", OUTPUT_FIG_REFLEX_PNG, "\n")
ggsave(OUTPUT_FIG_REFLEX_PNG, plot = p_reflex, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(OUTPUT_FIG_REFLEX_SVG, plot = p_reflex, width = 8, height = 6, bg = "white")

# ==============================================================================
# Figure 2: Trajectory of Value – Reflex strategy
# ==============================================================================
df_reflex <- df_summary %>%
    filter(strategy_label == "Panel_Reflex_ES") %>%
    arrange(yield_mean)

p_traj <- ggplot(df_reflex, aes(x = yield_mean, y = cost_per_diag_mean)) +
    geom_path(
        arrow     = arrow(length = unit(0.3, "cm"), type = "closed"),
        color     = "grey60",
        linewidth = 1,
        linetype  = "solid"
    ) +
    geom_point(aes(fill = yield_mean, size = yield_mean), shape = 21, color = "black", stroke = 1.2) +
    geom_text_repel(
        aes(label = paste0(phenotype_factor, "\n", percent(yield_mean, accuracy = 0.1))),
        box.padding   = 0.6,
        point.padding = 0.6,
        size          = 4.5,
        fontface      = "bold",
        color         = "black",
        bg.color      = "white",
        bg.r          = 0.15
    ) +
    scale_x_continuous(labels = percent_format(accuracy = 1), name = "Diagnostic Yield") +
    scale_y_continuous(
        labels = dollar_format(prefix = "$", suffix = " CAD"),
        name   = "Cost per Diagnosis",
        limits = c(0, NA)
    ) +
    scale_fill_viridis_c(option = "C", name = "Yield", guide = "none") +
    scale_size_continuous(range = c(6, 12), guide = "none") +
    annotate("text",
        x = max(df_reflex$yield_mean), y = -Inf, vjust = -0.5, hjust = 1.0,
        label = "High Efficiency Zone", size = 5, fontface = "bold", color = "#E65100"
    ) +
    labs(
        title    = "The Trajectory of Value",
        subtitle = "Reflex Strategy: Phenotype enrichment drives higher yield and lower cost per diagnosis"
    ) +
    theme_classic(base_size = 14) +
    theme(
        axis.line.x    = element_line(color = "black"),
        axis.line.y    = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title     = element_text(face = "bold", size = 18),
        plot.subtitle  = element_text(face = "italic", size = 12)
    )

cat("Saving Figure 2 (Trajectory of Value) to:", OUTPUT_FIG_TRAJ_PNG, "\n")
ggsave(OUTPUT_FIG_TRAJ_PNG, plot = p_traj, width = 9, height = 7, dpi = 300, bg = "white")
ggsave(OUTPUT_FIG_TRAJ_SVG, plot = p_traj, width = 9, height = 7, bg = "white")

# ==============================================================================
# Figure 3: All-modes efficiency landscape (small multiples)
# ==============================================================================
strategy_cols <- c(
    "Panel"          = "grey50",
    "Panel_Reflex_ES"= "#E65100",
    "ES"             = "#1976D2",
    "GS"             = "#7B1FA2"
)
strategy_labels_clean <- c(
    "Panel"          = "Panel",
    "Panel_Reflex_ES"= "Reflex",
    "ES"             = "Exome",
    "GS"             = "Genome"
)

arrow_data <- df_summary %>%
    select(phenotype_factor, strategy_label, yield_mean, cost_per_diag_mean) %>%
    filter(strategy_label %in% c("Panel", "Panel_Reflex_ES")) %>%
    tidyr::pivot_wider(names_from = strategy_label, values_from = c(yield_mean, cost_per_diag_mean))

p_all <- ggplot() +
    geom_segment(
        data = arrow_data,
        aes(
            x    = yield_mean_Panel,          y    = cost_per_diag_mean_Panel,
            xend = yield_mean_Panel_Reflex_ES, yend = cost_per_diag_mean_Panel_Reflex_ES
        ),
        arrow     = arrow(length = unit(0.2, "cm")),
        color     = "grey70",
        linewidth = 0.8
    ) +
    geom_point(
        data  = df_summary,
        aes(x = yield_mean, y = cost_per_diag_mean, color = strategy_label),
        size  = 4, alpha = 0.9
    ) +
    facet_wrap(~phenotype_factor, nrow = 1, scales = "fixed") +
    scale_color_manual(values = strategy_cols, labels = strategy_labels_clean, name = "Strategy") +
    scale_x_continuous(labels = percent_format(accuracy = 1), name = "Diagnostic Yield") +
    scale_y_continuous(
        labels = dollar_format(prefix = "$", scale = 1e-3, suffix = "k"),
        name   = "Cost per Diagnosis (CAD)"
    ) +
    theme_bw(base_size = 14) +
    theme(
        strip.background = element_rect(fill = "grey95"),
        strip.text       = element_text(face = "bold"),
        legend.position  = "bottom",
        panel.grid.minor = element_blank()
    ) +
    labs(
        title    = "Efficiency Landscape by Phenotype",
        subtitle = "Comparing strategies across clinical subgroups. Arrow indicates upgrade from Panel to Reflex."
    )

cat("Saving Figure 3 (All Modes) to:", OUTPUT_FIG_ALL_PNG, "\n")
ggsave(OUTPUT_FIG_ALL_PNG, plot = p_all, width = 14, height = 5, dpi = 300, bg = "white")
ggsave(OUTPUT_FIG_ALL_SVG, plot = p_all, width = 14, height = 5, bg = "white")

cat("Done. All figures written to:", OUTPUT_DIR, "\n")
