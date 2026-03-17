# 03_vis_parameters.R
# Purpose: Generate individual parameter and meta-analysis figures (PNG + CSV)
# Components:
# 1. Meta-analysis summary table (CSV)
# 2. Forest plot of diagnostic yields (PNG)
# 3. Beta density distributions (PNG)
# 4. Genetic landscape waterfall (PNG)
# 5. Gene architecture by phenotype (PNG)
# 6. Joint proportions table (CSV)
# Author: Renal Genetics CEA Team
# Date: 2025-12-22

# Set CRAN mirror
local({
    r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org"
    options(repos = r)
})

# Libraries (readxl and meta removed - using pre-computed values from RDS)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tidyr, ggplot2, readr, purrr, stringr, scales, viridis, patchwork, forcats)

# IO Paths
INPUT_PARAMS <- "data/params"
OUTPUT_CSV_DIR <- "outputs/parameters"
OUTPUT_FIG_DIR <- "outputs/figures/parameters"
OUTPUT_META_CSV <- file.path(OUTPUT_CSV_DIR, "meta_analysis_summary.csv")
OUTPUT_JOINT_CSV <- file.path(OUTPUT_CSV_DIR, "joint_proportions_full.csv")
OUTPUT_JOINT_TBL_CSV <- file.path(OUTPUT_CSV_DIR, "empirical_joint_proportions.csv")

# Figure output paths (individual PNGs, 300 dpi)
OUTPUT_FOREST_PNG <- file.path(OUTPUT_FIG_DIR, "forest_diagnostic_yield.png")
OUTPUT_BETA_PNG <- file.path(OUTPUT_FIG_DIR, "beta_density_distributions.png")
OUTPUT_LANDSCAPE_PNG <- file.path(OUTPUT_FIG_DIR, "genetic_landscape_top20.png")
OUTPUT_GENE_ARCH_PNG <- file.path(OUTPUT_FIG_DIR, "gene_architecture_by_phenotype.png")

# Ensure output directories exist
if (!dir.exists(OUTPUT_CSV_DIR)) dir.create(OUTPUT_CSV_DIR, recursive = TRUE)
if (!dir.exists(OUTPUT_FIG_DIR)) dir.create(OUTPUT_FIG_DIR, recursive = TRUE)

cat("--- Starting Parameter & Meta-Analysis Visualization ---\n")

# ==============================================================================
# Helper Functions & Theme
# ==============================================================================
# Helper for Phenotype Labels
format_phenotype <- function(x) {
    x <- gsub("_", " ", x)
    x <- str_to_title(x)
    x <- gsub("Ckdu", "CKDu", x) # Fix capitalization
    return(x)
}

theme_atlas <- theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "gray30", size = 11),
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 15, 10, 10), # Extra right margin for labels
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9)
    )

# ==============================================================================
# 1. Load Data (from validated intermediate RDS)
# ==============================================================================
INPUT_RDS <- "data/intermediate/01_imported_data.rds"

# Check data exists
if (!file.exists(INPUT_RDS)) {
    stop("Validated data not found. Run 01_import_validate.R first: ", INPUT_RDS)
}

data_raw <- readRDS(INPUT_RDS)
cohort <- data_raw$cohort
cat("Loaded cohort data:", nrow(cohort), "rows\n")
cat("Phenotype remapping note:", data_raw$metadata$phenotype_remapping, "\n")

# Aggregate cohort for downstream use
cohort_agg <- cohort %>%
    group_by(study_id, phenotype) %>%
    summarise(
        n_total = sum(n_total, na.rm = TRUE),
        n_monogenic = sum(n_monogenic, na.rm = TRUE),
        .groups = "drop"
    )

# Computed Parameters
yield_betas <- readRDS(file.path(INPUT_PARAMS, "yield_betas.rds"))
gene_variant_joint <- readRDS(file.path(INPUT_PARAMS, "gene_variant_joint.rds"))
difficult_loci <- readRDS(file.path(INPUT_PARAMS, "difficult_loci.rds"))
difficult_by_phenotype <- readRDS(file.path(INPUT_PARAMS, "difficult_by_phenotype.rds"))
phenotypes <- names(yield_betas)

# ==============================================================================
# Component 1: Meta-Analysis Table (Using pre-computed values - NO re-running metaprop)
# ==============================================================================
cat("1. Generating Meta-Analysis Table (from pre-computed yield_betas)...\n")

meta_list <- list()
for (p in phenotypes) {
    yb <- yield_betas[[p]]
    meta_list[[p]] <- data.frame(
        Phenotype = format_phenotype(p),
        k = yb$n_studies,
        N = sum((cohort_agg %>% filter(phenotype == p))$n_total),
        Yield = sprintf("%.1f%%", yb$mean * 100),
        CI_95 = sprintf("[%.1f-%.1f]", yb$ci_lower * 100, yb$ci_upper * 100),
        I2 = sprintf("%.0f%%", yb$i_squared * 100),
        # Raw values for CSV
        Yield_Raw = yb$mean,
        CI_Lower_Raw = yb$ci_lower,
        CI_Upper_Raw = yb$ci_upper,
        Alpha = yb$params["shape1"],
        Beta = yb$params["shape2"]
    )
}
meta_df_full <- bind_rows(meta_list)

# Export Full CSV
write_csv(meta_df_full, OUTPUT_META_CSV)
cat("   -> Exported:", OUTPUT_META_CSV, "\n")


# ==============================================================================
# Component 2: Forest Plot (Robust)
# ==============================================================================
cat("2. Generating Forest Plot...\n")

forest_data <- meta_df_full %>%
    select(Phenotype, Mean = Yield_Raw, Lower = CI_Lower_Raw, Upper = CI_Upper_Raw, N) %>%
    arrange(Mean)

# Add N to label for immediate context
forest_data <- forest_data %>%
    mutate(Label = paste0(Phenotype, " (n=", comma(N), ")"))

forest_data$Label <- factor(forest_data$Label, levels = forest_data$Label)

p_forest <- ggplot(forest_data, aes(y = Label, x = Mean)) +
    # Robust error bar (orientation = y)
    geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 0.2, color = "steelblue", linewidth = 0.8) +
    geom_point(size = 4, color = "darkblue") +
    geom_text(aes(label = sprintf("%.1f%%", Mean * 100)), vjust = -1.2, size = 3.5, fontface = "bold") +
    # Removed "Global Mean" line (misleading)
    scale_x_continuous(labels = percent_format(), limits = c(0, 0.8)) +
    labs(
        title = "Meta-Analysis of Diagnostic Yield",
        x = "Diagnostic Yield (95% CI)", y = NULL
    ) +
    theme_atlas +
    theme(
        panel.grid.major.y = element_line(color = "gray90"),
        plot.margin = margin(10, 20, 20, 10) # Increased bottom margin for x-axis
    )

# ==============================================================================
# Component 3: Beta Density Plots
# ==============================================================================
cat("3. Generating Beta Density Plots...\n")

beta_plot_df <- list()
for (p in phenotypes) {
    params <- yield_betas[[p]]$params
    x_vals <- seq(0, 1, length.out = 300)
    dens <- dbeta(x_vals, params["shape1"], params["shape2"])
    beta_plot_df[[p]] <- data.frame(
        Phenotype = format_phenotype(p),
        x = x_vals,
        y = dens,
        Mean = yield_betas[[p]]$mean
    )
}
beta_df <- bind_rows(beta_plot_df)
# Match order to meta-analysis yield
ordered_phenos <- gsub(" \\(n=.*\\)", "", levels(forest_data$Label))
beta_df$Phenotype <- factor(beta_df$Phenotype, levels = ordered_phenos)

p_betas <- ggplot(beta_df, aes(x, y)) +
    geom_area(fill = "#69b3a2", alpha = 0.4) +
    geom_line(color = "#404080", linewidth = 0.8) +
    geom_vline(aes(xintercept = Mean), linetype = "dashed", color = "black", alpha = 0.6) +
    facet_wrap(~Phenotype, scales = "free_y", ncol = 2) + # 2 columns for Portrait
    scale_x_continuous(labels = percent_format(), breaks = c(0, 0.5, 1)) +
    labs(
        title = "Probabilistic Priors: Fitted Beta Distributions",
        x = "Diagnostic Probability", y = "Density"
    ) +
    theme_atlas +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "gray95", color = NA),
        plot.margin = margin(10, 15, 30, 15) # Extra bottom/right room
    )

# ==============================================================================
# Component 4: Genetic Landscape (3-Panel Waterfall Chart)
# ==============================================================================
cat("4. Generating Genetic Landscape Figure (Waterfall)...\n")

# Step 1: Merge gene_variant_joint with difficult_by_phenotype to include PKD1/PKD2
u_joint <- bind_rows(lapply(names(gene_variant_joint), function(p) {
    df <- gene_variant_joint[[p]]
    df$Phenotype <- format_phenotype(p)
    df$Note <- NA_character_
    df
}))

# Add difficult locus data (PKD1, PKD2)
difficult_formatted <- difficult_by_phenotype %>%
    mutate(Phenotype = format_phenotype(phenotype)) %>%
    select(gene, variant_class, count, proportion, Phenotype) %>%
    mutate(
        modal_inheritance = "AD", # PKD1/PKD2 are AD
        Note = "Difficult Locus (Separate Model)"
    )

# Combine all gene data and recalculate GLOBAL proportions per phenotype
all_gene_data <- bind_rows(u_joint, difficult_formatted) %>%
    group_by(Phenotype) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()

# Export Joint CSV (updated with difficult loci)
write_csv(all_gene_data, OUTPUT_JOINT_CSV)
cat("   -> Exported:", OUTPUT_JOINT_CSV, "\n")

# Step 2: Compute per-gene totals
gene_totals <- all_gene_data %>%
    filter(gene != "Other") %>%
    group_by(gene) %>%
    summarise(total_count = sum(count), .groups = "drop") %>%
    arrange(desc(total_count)) %>%
    slice_head(n = 20)

top_genes <- gene_totals$gene
top_gene_data <- all_gene_data %>% filter(gene %in% top_genes)

# Step 3: Prepare data for each panel
# Gene order (descending by total count)
gene_order <- gene_totals$gene
diff_genes <- unique(difficult_loci$gene)
gene_labels <- ifelse(gene_order %in% diff_genes, paste0(gene_order, "*"), gene_order)

# Variant class distribution per gene (for left panel)
variant_dist <- top_gene_data %>%
    group_by(gene, variant_class) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    group_by(gene) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup() %>%
    mutate(gene_label = factor(
        ifelse(gene %in% diff_genes, paste0(gene, "*"), gene),
        levels = rev(gene_labels)
    ))

# Phenotype distribution per gene (for right panel)
phenotype_dist <- top_gene_data %>%
    group_by(gene, Phenotype) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    group_by(gene) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup() %>%
    mutate(gene_label = factor(
        ifelse(gene %in% diff_genes, paste0(gene, "*"), gene),
        levels = rev(gene_labels)
    ))

# Gene totals for center panel
gene_totals_plot <- gene_totals %>%
    mutate(gene_label = factor(
        ifelse(gene %in% diff_genes, paste0(gene, "*"), gene),
        levels = rev(gene_labels)
    ))

# Step 4: Color palettes
variant_colors <- c(
    "SNV" = "#2C3E50", # Dark Slate Blue
    "INDEL" = "#5D6D7E", # Lighter Slate Blue/Grey
    "CNV" = "#8E44AD", # Wisteria Purple (Unified CNV/SV)
    "Other" = "#95A5A6" # Gray
)

phenotype_colors <- c(
    "Glomerular" = "#3498DB", # Blue
    "Cystic" = "#27AE60", # Green
    "Tubulopathies" = "#F39C12", # Orange
    "Tubulointerstitial" = "#9B59B6", # Purple
    "CKDu" = "#7F8C8D" # Gray
)

# Step 5: Create the three panels
# Left Panel: Variant Class Distribution - Stacked bar
# Use shorter axis labels if needed or rely on legend
p_left <- ggplot(variant_dist, aes(x = prop, y = gene_label, fill = variant_class)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_manual(values = variant_colors, name = "Class") + # Shortened Legend
    scale_x_reverse(labels = percent_format(), limits = c(1, 0)) +
    labs(x = NULL, y = NULL, title = "Variant Class") +
    theme_atlas +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + # Force horizontal rows
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal", # Horizontal box
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 10, hjust = 0.5),
        plot.margin = margin(5, 0, 10, 5)
    )

# Center Panel: Gene Prevalence (Waterfall Bars)
p_center <- ggplot(gene_totals_plot, aes(x = total_count, y = gene_label)) +
    geom_bar(stat = "identity", fill = "#34495E", width = 0.7) +
    geom_text(aes(label = total_count), hjust = -0.2, size = 3, fontface = "bold") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = "n", y = NULL, title = "Diagnoses") +
    theme_atlas +
    theme(
        axis.text.y = element_text(size = 10, face = "bold", color = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "gray90"),
        plot.title = element_text(size = 11, hjust = 0.5),
        plot.margin = margin(5, 5, 5, 5)
    )

# Right Panel: Phenotype Distribution
p_right <- ggplot(phenotype_dist, aes(x = prop, y = gene_label, fill = Phenotype)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_manual(values = phenotype_colors, name = "Phenotype") + # Shortened Legend
    scale_x_continuous(labels = percent_format(), limits = c(0, 1)) +
    labs(x = NULL, y = NULL, title = "Phenotype") +
    theme_atlas +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + # Force horizontal rows
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 10, hjust = 0.5),
        plot.margin = margin(5, 5, 10, 0)
    )

# Step 6: Combine panels with patchwork
p_landscape <- p_left + p_center + p_right +
    plot_layout(widths = c(1.5, 2, 1.5), guides = "collect") +
    plot_annotation(
        title = "Genetic Landscape of Top 20 Causal Genes",
        subtitle = "Genes ranked by total diagnoses. Asterisks (*) denote difficult loci (PKD1, PKD2) modeled separately.",
        theme = theme(
            plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
            plot.subtitle = element_text(color = "gray30", size = 11, hjust = 0.5),
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.margin = margin(0, 0, 0, 0),
            legend.box.margin = margin(-5, 0, 0, 0)
        )
    )

# ==============================================================================
# Component 5: Joint Proportions Table
# ==============================================================================
cat("5. Generating Joint Proportions Table...\n")

joint_tbl_view <- all_gene_data %>%
    arrange(Phenotype, desc(proportion)) %>%
    group_by(Phenotype) %>%
    slice_head(n = 5) %>%
    mutate(
        Proportion = sprintf("%.1f%%", proportion * 100),
        Inheritance = modal_inheritance
    ) %>%
    select(Phenotype, Gene = gene, "Variant Class" = variant_class, Proportion, Inheritance, Note)

# Export Table 3 to CSV with Caption
caption_text <- "# Table 2. Empirical joint (Gene, Variant Class) proportions derived from weighted proband counts for multinomial sampling in the CEA model. Top 5 pairs per phenotype shown."
writeLines(caption_text, OUTPUT_JOINT_TBL_CSV)
write_csv(joint_tbl_view, OUTPUT_JOINT_TBL_CSV, append = TRUE, col_names = TRUE)
cat("   -> Exported:", OUTPUT_JOINT_TBL_CSV, "\n")



# ==============================================================================
# Component 6: Gene and Variant-Class Architecture by Phenotype (standalone PNG)
# ==============================================================================
cat("6. Generating Gene Architecture by Phenotype figure...\n")

# OUTPUT_GENE_ARCH_PNG already defined in IO Paths above

TOP_N_GENES <- 8

# Difficult locus gene names (already loaded as difficult_loci)
diff_gene_names <- unique(difficult_loci$gene)

# Assign clean variant class labels; flag difficult locus genes
arch_raw <- all_gene_data %>%
    filter(!is.na(gene)) %>%
    mutate(vc_clean = case_when(
        gene %in% diff_gene_names                         ~ "Difficult Locus",
        grepl("CNV|SV", variant_class, ignore.case = TRUE) ~ "CNV/SV",
        TRUE                                               ~ variant_class
    ))

# Top N genes per phenotype by summed within-phenotype proportion
top_genes_tbl <- arch_raw %>%
    group_by(Phenotype, gene) %>%
    summarise(gene_prop = sum(proportion), .groups = "drop") %>%
    group_by(Phenotype) %>%
    slice_max(gene_prop, n = TOP_N_GENES, with_ties = FALSE) %>%
    select(Phenotype, gene, gene_prop)

# Collapse non-top genes to "Other"
arch_plot <- arch_raw %>%
    left_join(top_genes_tbl, by = c("Phenotype", "gene")) %>%
    mutate(
        gene_lbl = if_else(!is.na(gene_prop), gene, "Other"),
        vc_lbl   = if_else(!is.na(gene_prop), vc_clean, "Other")
    ) %>%
    group_by(Phenotype, gene_lbl, vc_lbl) %>%
    summarise(prop = sum(proportion), .groups = "drop") %>%
    group_by(Phenotype) %>%
    mutate(prop = prop / sum(prop)) %>%   # re-normalise to 1 within phenotype
    ungroup()

# Order phenotypes by pooled diagnostic yield (highest first)
pheno_yield_vals  <- sapply(names(yield_betas), function(p) yield_betas[[p]]$mean)
pheno_lvls        <- format_phenotype(names(sort(pheno_yield_vals, decreasing = TRUE)))
arch_plot         <- arch_plot %>% mutate(Phenotype = factor(Phenotype, levels = pheno_lvls))

# Per-phenotype gene rank for x-axis ordering
gene_rank_tbl <- arch_plot %>%
    filter(gene_lbl != "Other") %>%
    group_by(Phenotype, gene_lbl) %>%
    summarise(gene_total = sum(prop), .groups = "drop") %>%
    group_by(Phenotype) %>%
    arrange(Phenotype, desc(gene_total)) %>%
    mutate(gene_rank = row_number()) %>%
    ungroup()

arch_plot <- arch_plot %>%
    left_join(gene_rank_tbl %>% select(Phenotype, gene_lbl, gene_rank),
              by = c("Phenotype", "gene_lbl")) %>%
    mutate(
        gene_rank = if_else(gene_lbl == "Other", 999L, as.integer(gene_rank)),
        # Create per-phenotype unique key for free_x facet ordering
        gene_in_pheno = paste(Phenotype, sprintf("%03d", gene_rank), gene_lbl, sep = "|||")
    ) %>%
    arrange(Phenotype, gene_rank) %>%
    mutate(gene_in_pheno = factor(gene_in_pheno, levels = unique(gene_in_pheno)))

# Colour palette
arch_colors <- c(
    "SNV"             = "#2C3E50",
    "INDEL"           = "#5D6D7E",
    "SNV/Indel"       = "#2C3E50",
    "CNV/SV"          = "#8E44AD",
    "Difficult Locus" = "#C0392B",
    "Other"           = "#BDC3C7"
)
# Add any variant classes not covered above
present_vc <- unique(arch_plot$vc_lbl)
missing_vc <- setdiff(present_vc, names(arch_colors))
if (length(missing_vc) > 0)
    arch_colors <- c(arch_colors, setNames(rep("#95A5A6", length(missing_vc)), missing_vc))

p_gene_arch <- ggplot(arch_plot, aes(x = gene_in_pheno, y = prop, fill = vc_lbl)) +
    geom_bar(stat = "identity", position = "stack", width = 0.75) +
    facet_wrap(~Phenotype, scales = "free_x", ncol = 3) +
    scale_x_discrete(labels = function(x) sub(".*\\|\\|\\|", "", x)) +
    scale_fill_manual(values = arch_colors, name = "Variant class") +
    scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        expand = expansion(mult = c(0, 0.05))
    ) +
    labs(x = NULL, y = "Proportion of P/LP diagnoses") +
    theme_atlas +
    theme(
        axis.text.x        = element_text(angle = 45, hjust = 1, size = 9, face = "italic"),
        strip.background   = element_rect(fill = "gray95", color = NA),
        strip.text         = element_text(face = "bold", size = 11),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray92"),
        legend.position    = "bottom"
    )

ggsave(OUTPUT_GENE_ARCH_PNG, p_gene_arch, width = 10, height = 7, dpi = 300, bg = "white")
cat("   -> Exported:", OUTPUT_GENE_ARCH_PNG, "\n")

# ==============================================================================
# Figure Export (Individual PNGs, 300 dpi)
# ==============================================================================
cat("Exporting individual figures...\n")

# Forest plot
ggsave(OUTPUT_FOREST_PNG, p_forest, width = 8, height = 5, dpi = 300, bg = "white")
cat("   -> Exported:", OUTPUT_FOREST_PNG, "\n")

# Beta density distributions
ggsave(OUTPUT_BETA_PNG, p_betas, width = 8, height = 6, dpi = 300, bg = "white")
cat("   -> Exported:", OUTPUT_BETA_PNG, "\n")

# Genetic landscape (wider for 3-panel layout)
# Suppress known patchwork false-positive warning about plot_annotation theme
suppressWarnings(ggsave(OUTPUT_LANDSCAPE_PNG, p_landscape, width = 12, height = 8, dpi = 300, bg = "white"))
cat("   -> Exported:", OUTPUT_LANDSCAPE_PNG, "\n")

cat("SUCCESS: Figures saved to", OUTPUT_FIG_DIR, "\n")
cat("SUCCESS: CSVs saved to", OUTPUT_CSV_DIR, "\n")
