# 07b_prcc_composite_figure.R
# Rebuild the PRCC composite figure from already-computed CSV results.
# Shows Reflex-to-ES vs Panel only (A = Incremental Cost, B = Incremental Diagnoses).

library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)
library(grid)

OUTPUT_DIR <- "outputs/results/uncertainty_sensitivity/global_sensitivity_analysis_prcc"

# ==============================================================================
# Labels and helpers
# ==============================================================================
param_labels <- c(
    "prev_CKDu"               = "CKDu prevalence",
    "prev_cystic"             = "Cystic prevalence",
    "prev_glomerular"         = "Glomerular prevalence",
    "prev_tubulointerstitial" = "Tubulointerstitial prevalence",
    "prev_tubulopathies"      = "Tubulopathies prevalence",
    "yield_CKDu"              = "CKDu yield",
    "yield_cystic"            = "Cystic yield",
    "yield_glomerular"        = "Glomerular yield",
    "yield_tubulointerstitial"= "Tubulointerstitial yield",
    "yield_tubulopathies"     = "Tubulopathies yield",
    "prob_vus_panel_mean"     = "VUS probability (Panel)",
    "prob_vus_es"             = "VUS probability (ES)",
    "prob_vus_gs"             = "VUS probability (GS)",
    "det_snv_es"              = "SNV detection (ES)",
    "det_cnv_es"              = "CNV detection (ES)",
    "det_snv_gs"              = "SNV detection (GS)",
    "det_cnv_gs"              = "CNV detection (GS)",
    "det_snv_panel_weighted"  = "SNV detection (Panel)",
    "det_cnv_panel_weighted"  = "CNV detection (Panel)",
    "det_pkd1_es"             = "PKD1 detection (ES)",
    "det_pkd1_gs"             = "PKD1 detection (GS)",
    "det_muc1_es"             = "MUC1 detection (ES)",
    "det_muc1_gs"             = "MUC1 detection (GS)",
    "det_pkd1_panel"          = "PKD1 detection (Cystic panel)",
    "det_muc1_panel"          = "MUC1 detection (ADTKD panel)",
    "cost_panel"              = "Panel Cost",
    "cost_es"                 = "Exome Sequencing Cost",
    "cost_gs"                 = "Genome Sequencing Cost",
    "cost_familial"           = "Targeted Familial Test Cost",
    "cost_consult_pre"        = "Pre-test Consultation Cost",
    "cost_consult_post"       = "Post-test Consultation Cost"
)

get_label <- function(param_name) {
    clean_name <- gsub("\\.", "_", param_name)
    if (clean_name %in% names(param_labels)) return(param_labels[[clean_name]])
    tools::toTitleCase(gsub("_", " ", clean_name))
}

assign_category <- function(param_name) {
    if (grepl("^prev_", param_name))              return("Phenotype Prevalence")
    if (grepl("^cost_", param_name))              return("Unit Cost")
    if (grepl("det_pkd1|det_muc1", param_name))  return("Difficult Locus Detection")
    if (grepl("^det_", param_name))               return("Detection Sensitivity")
    if (grepl("^yield_", param_name))             return("Diagnostic Yield")
    if (grepl("^prob_vus", param_name))           return("VUS Probability")
    "Other"
}

# Build a tornado panel; legend is suppressed — call get_legend() separately.
plot_tornado_panel <- function(df, outcome_col, panel_title, panel_subtitle,
                               show_legend = FALSE, top_n = 20) {
    ci_low_col  <- paste0(outcome_col, "_ci_low")
    ci_high_col <- paste0(outcome_col, "_ci_high")

    plot_data <- df %>%
        mutate(
            value          = .data[[outcome_col]],
            ci_low         = .data[[ci_low_col]],
            ci_high        = .data[[ci_high_col]],
            is_significant = sign(ci_low) == sign(ci_high),
            display_label  = sapply(parameter_name, get_label),
            category       = sapply(parameter_name, assign_category)
        ) %>%
        arrange(desc(abs(value))) %>%
        slice_head(n = top_n) %>%
        mutate(display_label = factor(display_label, levels = rev(display_label)))

    max_abs <- max(abs(c(plot_data$ci_low, plot_data$ci_high, plot_data$value)), na.rm = TRUE)
    lim <- max(0.1, ceiling(max_abs * 10) / 10)

    leg_pos <- if (show_legend) "bottom" else "none"

    ggplot(plot_data, aes(x = display_label, y = value, fill = category)) +
        geom_bar(stat = "identity", width = 0.7, aes(alpha = is_significant)) +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        geom_errorbar(
            aes(ymin = ci_low, ymax = ci_high, alpha = is_significant),
            width = 0.2, color = "gray20", linewidth = 0.4
        ) +
        coord_flip() +
        scale_fill_brewer(palette = "Set2", name = "Parameter category") +
        scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1), guide = "none") +
        scale_y_continuous(limits = c(-lim, lim), breaks = scales::pretty_breaks(n = 5)) +
        labs(
            title    = panel_title,
            subtitle = panel_subtitle,
            x        = NULL,
            y        = "Partial Rank Correlation Coefficient (PRCC)"
        ) +
        theme_minimal(base_size = 11) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text        = element_text(color = "black", size = 10),
            axis.title.x     = element_text(color = "black", size = 11, face = "bold"),
            plot.title       = element_text(size = 12, face = "bold", hjust = 0),
            plot.subtitle    = element_text(size = 9, color = "gray30", hjust = 0),
            legend.position  = leg_pos,
            legend.direction = "horizontal",
            legend.text      = element_text(size = 9),
            legend.title     = element_text(size = 9, face = "bold"),
            legend.key.size  = unit(0.45, "cm")
        )
}

# Extract the legend grob from a ggplot
get_legend_grob <- function(p) {
    g <- ggplotGrob(p)
    idx <- which(vapply(g$grobs, function(x) {
        !is.null(x$name) && grepl("guide-box", x$name, fixed = TRUE)
    }, logical(1)))
    if (length(idx) == 0L) stop("Legend grob not found.")
    g$grobs[[idx[1]]]
}

# ==============================================================================
# Load PRCC results (reflex vs panel only)
# ==============================================================================
min_thr <- 0.02

rvp <- read_csv(file.path(OUTPUT_DIR, "prcc_results_reflex_vs_panel.csv"),
                show_col_types = FALSE)

# ==============================================================================
# Build panels (no legend)
# ==============================================================================
p_A <- plot_tornado_panel(
    filter(rvp, abs(prcc_delta_cost_cad) >= min_thr),
    "prcc_delta_cost_cad",
    "Reflex to ES vs Panel",
    "Outcome: Incremental Cost",
    show_legend = FALSE
)

p_B <- plot_tornado_panel(
    filter(rvp, abs(prcc_delta_diagnoses) >= min_thr &
               sapply(parameter_name, assign_category) != "Unit Cost"),
    "prcc_delta_diagnoses",
    "Reflex to ES vs Panel",
    "Outcome: Incremental Diagnoses",
    show_legend = FALSE
)

# Build legend from panel A (with legend shown) — extract grob, then discard plot
p_legend_src <- plot_tornado_panel(
    filter(rvp, abs(prcc_delta_cost_cad) >= min_thr),
    "prcc_delta_cost_cad",
    "", "",
    show_legend = TRUE
)
legend_grob <- get_legend_grob(p_legend_src)

# ==============================================================================
# Composite: panels side by side, legend below
# ==============================================================================
panels <- p_A + p_B +
    plot_layout(ncol = 2) +
    plot_annotation(tag_levels = "A")

composite <- wrap_elements(full = ggplotGrob(
    # dummy — we rebuild as gtable below
    ggplot() + theme_void()
))

# Assemble with patchwork: panels / legend row
composite <- (p_A + p_B +
    plot_layout(ncol = 2) +
    plot_annotation(tag_levels = "A")) /
    wrap_elements(full = legend_grob) +
    plot_layout(heights = c(10, 0.8))

ggsave(
    file.path(OUTPUT_DIR, "prcc_tornado_composite.png"),
    composite,
    width = 13, height = 6.5, dpi = 300, bg = "white"
)
ggsave(
    file.path(OUTPUT_DIR, "prcc_tornado_composite.svg"),
    composite,
    width = 13, height = 6.5, bg = "white"
)

cat("Saved composite PRCC figure (A-B, single legend).\n")
