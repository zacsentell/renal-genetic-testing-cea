# 07b_prcc_composite_figure.R
# Purpose: Build composite PRCC figure (forest plot + heatmap) from CSV outputs.
#
# Design:
#   Panel A — ranked forest plot for the primary comparison (Reflex vs Panel).
#             Shows both incremental cost and yield outcomes as dodged dots with
#             95% CI lines. Parameters are selected as the union of:
#               (i)  top 6 from the primary comparison, and
#               (ii) any parameter with |PRCC| >= 0.2 and CI excluding zero in
#                    any comparison.
#             Parameters significant only in scenario comparisons appear faded
#             at the bottom, communicating they do not drive the primary decision.
#   Panel B — PRCC heatmap for all three comparisons (Reflex vs Panel,
#             ES vs Panel, GS vs Reflex), showing the same parameter rows as
#             Panel A. Fill = signed PRCC on a diverging blue-white-red palette.
#             Numeric values printed in cells where |PRCC| >= 0.1 and CI
#             excludes zero. Grey cells = not material or not significant.
#
# Inputs (produced by 07_global_sensitivity_analysis.R):
#   prcc_results_reflex_vs_panel.csv
#   prcc_results_es_vs_panel.csv
#   prcc_results_genome_vs_reflex.csv
#
# Outputs:
#   prcc_tornado_composite.png  (replaces old 2x2 tornado composite)
#   prcc_tornado_composite.pdf

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(patchwork)

OUTPUT_DIR <- "outputs/results/uncertainty_sensitivity/global_sensitivity_analysis_prcc"

# ==============================================================================
# 1. Parameter labels and categories
# ==============================================================================
param_labels <- c(
    "prev_CKDu"               = "Prevalence: CKDu",
    "prev_cystic"             = "Prevalence: Cystic",
    "prev_glomerular"         = "Prevalence: Glomerular",
    "prev_tubulointerstitial" = "Prevalence: Tubulointerstitial",
    "prev_tubulopathies"      = "Prevalence: Tubulopathies",
    "yield_CKDu"              = "Monogenic probability: CKDu",
    "yield_cystic"            = "Monogenic probability: Cystic",
    "yield_glomerular"        = "Monogenic probability: Glomerular",
    "yield_tubulointerstitial"= "Monogenic probability: Tubulointerstitial",
    "yield_tubulopathies"     = "Monogenic probability: Tubulopathies",
    "prob_vus_panel_mean"     = "VUS probability: Panel",
    "prob_vus_es"             = "VUS probability: ES",
    "prob_vus_gs"             = "VUS probability: GS",
    "det_snv_es"              = "SNV/indel detection: ES",
    "det_cnv_es"              = "CNV/SV detection: ES",
    "det_snv_gs"              = "SNV/indel detection: GS",
    "det_cnv_gs"              = "CNV/SV detection: GS",
    "det_snv_panel_weighted"  = "SNV/indel detection: Panel",
    "det_cnv_panel_weighted"  = "CNV/SV detection: Panel",
    "det_pkd1_es"             = "PKD1 detection: ES",
    "det_pkd1_gs"             = "PKD1 detection: GS",
    "det_muc1_es"             = "MUC1 detection: ES",
    "det_muc1_gs"             = "MUC1 detection: GS",
    "det_pkd1_panel"          = "PKD1 detection: Cystic panel",
    "det_muc1_panel"          = "MUC1 detection: ADTKD panel",
    "cost_panel"              = "Panel unit cost",
    "cost_es"                 = "Exome unit cost",
    "cost_gs"                 = "Genome unit cost",
    "cost_familial"           = "Familial variant test cost",
    "cost_consult_pre"        = "Pre-test consultation cost",
    "cost_consult_post"       = "Post-test consultation cost",
    "uptake_reflex"           = "Reflex ES uptake probability",
    "uptake_cascade"          = "Cascade uptake probability"
)

get_label <- function(param_name) {
    clean <- gsub("\\.", "_", param_name)
    if (clean %in% names(param_labels)) return(param_labels[[clean]])
    tools::toTitleCase(gsub("_", " ", clean))
}

assign_category <- function(param_name) {
    if (grepl("^prev_",            param_name)) return("Phenotype prevalence")
    if (grepl("^cost_",            param_name)) return("Unit cost")
    if (grepl("det_pkd1|det_muc1", param_name)) return("Difficult-locus detection")
    if (grepl("^det_",             param_name)) return("Standard detection")
    if (grepl("^yield_",           param_name)) return("Monogenic probability")
    if (grepl("^prob_vus",         param_name)) return("VUS probability")
    if (grepl("^uptake_",          param_name)) return("Uptake probability")
    "Other"
}

# ==============================================================================
# 2. Load PRCC results
# ==============================================================================
# Comparisons in display order (A = primary, B, C = scenario)
comparison_meta <- list(
    list(label = "reflex_vs_panel", display = "Reflex vs Panel",  file = "prcc_results_reflex_vs_panel.csv"),
    list(label = "es_vs_panel",     display = "ES vs Panel",      file = "prcc_results_es_vs_panel.csv"),
    list(label = "genome_vs_reflex",display = "GS vs Reflex",     file = "prcc_results_genome_vs_reflex.csv")
)

prcc_results <- list()
for (m in comparison_meta) {
    path <- file.path(OUTPUT_DIR, m$file)
    if (!file.exists(path)) {
        message("PRCC result not found — skipping: ", path)
        next
    }
    df <- read_csv(path, show_col_types = FALSE)
    prcc_results[[m$label]] <- list(data = df, display = m$display)
}

if (length(prcc_results) == 0) stop("No PRCC result files found in ", OUTPUT_DIR)

primary_key <- if ("reflex_vs_panel" %in% names(prcc_results)) "reflex_vs_panel" else names(prcc_results)[1]
primary_df  <- prcc_results[[primary_key]]$data

# ==============================================================================
# 3. Determine union parameter set
# ==============================================================================
PRCC_SHOW_THRESH <- 0.10   # minimum |PRCC| to colour a heatmap cell (must also be sig)
PRCC_UNION_THRESH <- 0.20  # minimum |PRCC| for a parameter to enter the union set
N_PRIMARY <- 6             # top-N from primary comparison to always include

is_sig <- function(lo, hi) !is.na(lo) & !is.na(hi) & sign(lo) == sign(hi)

# Top-N from primary comparison
top_primary <- primary_df %>%
    mutate(max_abs = pmax(abs(prcc_delta_cost_cad), abs(prcc_delta_diagnoses), na.rm = TRUE)) %>%
    arrange(desc(max_abs)) %>%
    slice_head(n = N_PRIMARY) %>%
    pull(parameter_name)

# Any parameter with |PRCC| >= threshold and significant in any comparison
top_scenario <- unique(unlist(lapply(names(prcc_results), function(lbl) {
    df <- prcc_results[[lbl]]$data
    df %>%
        filter(
            (abs(prcc_delta_cost_cad)   >= PRCC_UNION_THRESH &
                 is_sig(prcc_delta_cost_cad_ci_low,   prcc_delta_cost_cad_ci_high)) |
            (abs(prcc_delta_diagnoses)  >= PRCC_UNION_THRESH &
                 is_sig(prcc_delta_diagnoses_ci_low, prcc_delta_diagnoses_ci_high))
        ) %>%
        pull(parameter_name)
})))

# Union: primary top-N first, then scenario-specific additions
param_union <- unique(c(top_primary, top_scenario))

# Order scenario additions by their max |PRCC| in any non-primary comparison
scenario_only <- setdiff(top_scenario, top_primary)
if (length(scenario_only) > 0) {
    scenario_ranks <- vapply(scenario_only, function(p) {
        max(sapply(names(prcc_results)[names(prcc_results) != primary_key], function(lbl) {
            df <- prcc_results[[lbl]]$data
            row <- df[df$parameter_name == p, ]
            if (nrow(row) == 0) return(0)
            max(abs(row$prcc_delta_cost_cad[1]), abs(row$prcc_delta_diagnoses[1]), na.rm = TRUE)
        }), na.rm = TRUE)
    }, numeric(1))
    scenario_only <- scenario_only[order(scenario_ranks, decreasing = TRUE)]
}
param_union <- unique(c(top_primary, scenario_only))

# Display labels and categories
param_display  <- setNames(vapply(param_union, get_label,      character(1)), param_union)
param_category <- setNames(vapply(param_union, assign_category, character(1)), param_union)

# y-axis factor levels: top-ranked at top of plot (ggplot reads bottom-up so we reverse)
y_levels <- rev(unname(param_display))

# ==============================================================================
# 4. Forest plot data (primary comparison, both outcomes)
# ==============================================================================
primary_df_reflex <- primary_df %>% filter(parameter_name %in% param_union)

# Rank params by primary comparison max |PRCC| (determines y position)
rank_order <- primary_df_reflex %>%
    mutate(max_abs = pmax(abs(prcc_delta_cost_cad), abs(prcc_delta_diagnoses), na.rm = TRUE),
           display_label = param_display[parameter_name]) %>%
    arrange(desc(max_abs)) %>%
    pull(display_label)
# Append scenario-only params at the bottom (they'll have near-zero primary PRCC)
rank_order_full <- unique(c(rank_order, unname(param_display[scenario_only])))
y_levels <- rev(rank_order_full)  # reversed for ggplot

forest_df <- bind_rows(
    primary_df_reflex %>%
        transmute(parameter_name,
                  outcome  = "Incremental cost",
                  prcc     = prcc_delta_cost_cad,
                  ci_low   = prcc_delta_cost_cad_ci_low,
                  ci_high  = prcc_delta_cost_cad_ci_high),
    primary_df_reflex %>%
        transmute(parameter_name,
                  outcome  = "Incremental yield",
                  prcc     = prcc_delta_diagnoses,
                  ci_low   = prcc_delta_diagnoses_ci_low,
                  ci_high  = prcc_delta_diagnoses_ci_high)
) %>%
    mutate(
        display_label = factor(param_display[parameter_name], levels = y_levels),
        category      = param_category[parameter_name],
        outcome       = factor(outcome, levels = c("Incremental cost", "Incremental yield")),
        significant   = is_sig(ci_low, ci_high)
    )

# ==============================================================================
# 5. Heatmap data (all comparisons)
# ==============================================================================
comparison_order <- vapply(comparison_meta, `[[`, character(1), "display")
comparison_order <- comparison_order[comparison_order %in%
                                         sapply(prcc_results, `[[`, "display")]

heatmap_parts <- list()
for (lbl in names(prcc_results)) {
    info <- prcc_results[[lbl]]
    df   <- info$data
    heatmap_parts[[paste0(lbl, "_cost")]] <- df %>%
        transmute(parameter_name,
                  comparison = info$display,
                  outcome    = "Cost",
                  prcc       = prcc_delta_cost_cad,
                  ci_low     = prcc_delta_cost_cad_ci_low,
                  ci_high    = prcc_delta_cost_cad_ci_high)
    heatmap_parts[[paste0(lbl, "_yield")]] <- df %>%
        transmute(parameter_name,
                  comparison = info$display,
                  outcome    = "Yield",
                  prcc       = prcc_delta_diagnoses,
                  ci_low     = prcc_delta_diagnoses_ci_low,
                  ci_high    = prcc_delta_diagnoses_ci_high)
}

heatmap_df <- bind_rows(heatmap_parts) %>%
    filter(parameter_name %in% param_union) %>%
    complete(parameter_name = param_union,
             comparison     = comparison_order,
             outcome        = c("Cost", "Yield")) %>%
    mutate(
        significant = is_sig(ci_low, ci_high),
        material    = !is.na(prcc) & abs(prcc) >= PRCC_SHOW_THRESH & significant,
        fill_value  = ifelse(material, prcc, NA_real_),
        # PRCC text with explicit sign; white text for high-magnitude cells
        text_label  = ifelse(material,
                             ifelse(prcc >= 0,
                                    sprintf("+%.2f", prcc),
                                    sprintf("%.2f",  prcc)),
                             ""),
        text_color  = ifelse(material & abs(prcc) >= 0.60, "white", "grey15"),
        display_label = factor(
            vapply(parameter_name,
                   function(p) param_display[[gsub("\\.", "_", p)]],
                   character(1)),
            levels = y_levels
        ),
        comparison = factor(comparison, levels = comparison_order),
        outcome    = factor(outcome, levels = c("Cost", "Yield"))
    )

# ==============================================================================
# 6. Forest plot (Panel A)
# ==============================================================================
# Use viridis_d for category colours; Set2 for internal contrast
cats_present  <- unique(param_category)
n_cats        <- length(cats_present)

forest_plot <- ggplot(
    forest_df,
    aes(x = prcc, y = display_label,
        color = category, shape = outcome, group = outcome)
) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey45", linewidth = 0.5) +
    # CI lines (no end-caps for a clean forest look)
    geom_errorbar(
        aes(xmin = ci_low, xmax = ci_high, alpha = significant),
        width       = 0,
        linewidth   = 0.65,
        orientation = "y",
        position    = position_dodge(width = 0.55)
    ) +
    # Point estimates
    geom_point(
        aes(alpha = significant),
        size     = 3.0,
        position = position_dodge(width = 0.55)
    ) +
    scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.22), guide = "none") +
    scale_shape_manual(
        values = c("Incremental cost" = 16, "Incremental yield" = 18),
        name   = "Outcome"
    ) +
    scale_color_viridis_d(
        option = "D",
        end    = 0.88,
        name   = "Parameter category"
    ) +
    scale_x_continuous(
        limits = c(-1, 1),
        breaks = c(-1, -0.5, 0, 0.5, 1),
        labels = c("-1", "-0.5", "0", "+0.5", "+1"),
        expand = expansion(mult = 0.02)
    ) +
    labs(
        title    = "A  Reflex vs Panel (primary comparison)",
        subtitle = paste0(
            "Parameters significant in any comparison, ranked by importance to primary decision; ",
            "faded = not significant here"
        ),
        x = "Partial rank correlation coefficient (PRCC)",
        y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.35),
        axis.text.y        = element_text(size = 9.5, color = "black"),
        axis.text.x        = element_text(size = 9,   color = "black"),
        axis.title.x       = element_text(size = 10,  face = "bold",
                                          margin = margin(t = 6)),
        plot.title         = element_text(size = 12, face = "bold"),
        plot.subtitle      = element_text(size = 8.5, color = "grey40",
                                          margin = margin(b = 4)),
        legend.position    = "bottom",
        legend.box         = "vertical",
        legend.text        = element_text(size = 9),
        legend.title       = element_text(size = 9, face = "bold"),
        legend.key.size    = unit(0.45, "cm"),
        plot.margin        = margin(t = 6, r = 6, b = 6, l = 6)
    ) +
    guides(
        color = guide_legend(order = 1, nrow = 3,
                             override.aes = list(shape = 16, size = 3,
                                                 alpha = 1)),
        shape = guide_legend(order = 2,
                             override.aes = list(color = "grey30",
                                                 size = 3, alpha = 1))
    )

# ==============================================================================
# 7. Heatmap (Panel B)
# ==============================================================================
heatmap_plot <- ggplot(
    heatmap_df,
    aes(x = outcome, y = display_label, fill = fill_value)
) +
    geom_tile(color = "white", linewidth = 0.9) +
    # Numeric PRCC values inside significant cells
    geom_text(
        aes(label = text_label, color = text_color),
        size     = 2.7,
        fontface = "plain",
        show.legend = FALSE
    ) +
    scale_fill_gradient2(
        low      = "#2166AC",   # blue (negative PRCC)
        mid      = "white",
        high     = "#B2182B",   # red  (positive PRCC)
        midpoint = 0,
        limits   = c(-1, 1),
        na.value = "grey94",
        name     = "PRCC",
        breaks   = c(-1, -0.5, 0, 0.5, 1),
        labels   = c("-1", "-0.5", "0", "+0.5", "+1"),
        guide    = guide_colorbar(
            barwidth       = 0.7,
            barheight      = 5.5,
            title.position = "top",
            title.hjust    = 0.5
        )
    ) +
    scale_color_identity(guide = "none") +
    facet_wrap(~ comparison, nrow = 1, strip.position = "top") +
    labs(
        title    = "B  All comparisons",
        subtitle = "Coloured cells: |PRCC| >= 0.10, 95% CI excludes zero. Grey: not material or not significant.",
        x        = NULL,
        y        = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        axis.text.x      = element_text(size = 9, color = "black"),
        strip.text       = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "grey92", color = NA),
        panel.grid       = element_blank(),
        plot.title       = element_text(size = 12, face = "bold"),
        plot.subtitle    = element_text(size = 8.5, color = "grey40",
                                        margin = margin(b = 4)),
        legend.position  = "right",
        legend.title     = element_text(size = 9.5, face = "bold"),
        legend.text      = element_text(size = 9),
        plot.margin      = margin(t = 6, r = 10, b = 6, l = 2)
    )

# ==============================================================================
# 8. Compose and save
# ==============================================================================
composite <- forest_plot + heatmap_plot +
    plot_layout(widths = c(3, 2))

out_png <- file.path(OUTPUT_DIR, "prcc_tornado_composite.png")
out_pdf <- file.path(OUTPUT_DIR, "prcc_tornado_composite.pdf")

ggsave(out_png, composite, width = 16, height = 9.0, dpi = 300, bg = "white")
ggsave(out_pdf, composite, width = 16, height = 9.0, bg = "white")

cat("Saved composite PRCC figure:\n  ", out_png, "\n  ", out_pdf, "\n")
cat("Parameter set (", length(param_union), " parameters):\n",
    paste0("  ", param_union, " -> ", param_display, "\n"), sep = "")
