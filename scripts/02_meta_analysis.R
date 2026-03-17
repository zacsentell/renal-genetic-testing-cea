# 02_meta_analysis.R
# Purpose: Generate probabilistic parameters (Beta, Dirichlet) from raw data.
# Author: Zachary Sentell
# Date: 2025-12-10

# Set CRAN mirror for non-interactive sessions
local({
    r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org"
    options(repos = r)
})

# Libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tidyr, meta, readr, purrr, stringr)

# IO Paths
INPUT_RDS <- "data/intermediate/01_imported_data.rds"
OUTPUT_PARAMS <- "data/params"
OUTPUT_AUDIT <- "outputs/parameters/02_parameter_summary.txt"

# Ensure directories exist
if (!dir.exists(OUTPUT_PARAMS)) dir.create(OUTPUT_PARAMS, recursive = TRUE)
if (!dir.exists(dirname(OUTPUT_AUDIT))) dir.create(dirname(OUTPUT_AUDIT), recursive = TRUE)

cat("--- Starting Meta-Analysis & Parameterization ---\n")

# Load Data
data_raw <- readRDS(INPUT_RDS)
cohort <- data_raw$cohort
variant <- data_raw$variant

# ==============================================================================
# 0. Note: Phenotype Remapping
# ==============================================================================
# Phenotype remapping (complement/aHUS → other_mixed) is now centralized in
# 01_import_validate.R. Data is already remapped in the imported RDS.

# ==============================================================================
# 1. Phenotype Prevalence (Dirichlet)
# ==============================================================================
cat("1. Processing Phenotype Prevalence...\n")

# Aggregate before pivot (handles merged phenotypes with multiple rows per study)
phenotype_wide <- cohort %>%
    group_by(study_id, phenotype) %>%
    summarise(n_total = sum(n_total, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = phenotype, values_from = n_total, values_fill = 0)

counts_matrix <- as.matrix(phenotype_wide[, -1])

# Simple Dirichlet estimation (Method of Moments on pooled proportions)
# Note: For CEA, we often use the global sum + effective sample size logic.
# Here we will calculate the global proportions and assign a confidence (N_eff).
# Using total N as proxy for N_eff for now, or we can use the dirichlet package if needed.
# Let's stick to the Plan: Pooled Means + Normalize.

pheno_summary <- cohort %>%
    group_by(phenotype) %>%
    summarise(total_n = sum(n_total, na.rm = TRUE)) %>%
    mutate(prop = total_n / sum(total_n))

# Parameter object: Counts (Alphas) for Dirichlet
phenotype_alphas <- setNames(pheno_summary$total_n, pheno_summary$phenotype)

# ==============================================================================
# 2. Diagnostic Yield (Beta via Meta-Analysis)
# ==============================================================================
cat("2. Meta-Analyzing Diagnostic Yields...\n")

phenotypes <- unique(cohort$phenotype)
yield_params <- list()
yield_audit_rows <- list()

for (p in phenotypes) {
    d_sub <- cohort %>% filter(phenotype == p)

    # Only run meta-analysis if we have >1 study with non-zero counts, else raw pool
    if (nrow(d_sub) > 0) {
        # Run random effects meta-analysis
        # metaprop handles 0 counts gracefully (usually adds 0.5 continuity correction)
        m <- metaprop(
            event = n_monogenic, n = n_total,
            data = d_sub,
            method = "GLMM",
            sm = "PLOGIT"
        ) # Logit transformation for proportions

        # Extract pooled estimate and variance (on logit scale approx, but metaprop gives us CI)
        # We will use the Random Effects Estimate
        mu <- m$TE.random.w # Logit mean
        se <- m$seTE.random # Logit SE

        # But for Beta, we need mu and var on the proportion scale [0,1].
        # metaprop provides 'prop' results in the summary object.
        # Let's simple-fit Beta to the Mean and CI using Method of Moments logic or
        # just use the meta-analyzed proportion and a conservative N_eff.

        pooled_prop <- plogis(m$TE.random) # Back-transform logit to prob

        # Approximate Alpha/Beta from Mean and an assumed information derived from CI width?
        # Or simpler: Alpha = Pooled_Prop * Total_N, Beta = (1-P)*Total_N?
        # That ignores heterogeneity.
        # Better: Use the 'weights' from random effects to get effective N?
        # Let's use: Mean = P_rand, Variance = SE_rand^2 (converted to prop scale)

        # Delta method approximation for variance on prop scale:
        # Var(p) approx (p*(1-p))^2 * Var(logit_p)
        var_prop <- (pooled_prop * (1 - pooled_prop))^2 * (m$seTE.random^2)

        # Method of Moments for Beta:
        # alpha = mu * ( (mu*(1-mu)/var) - 1 )
        # beta = (1 - mu) * ( (mu*(1-mu)/var) - 1 )

        # Guard against zero/near-zero variance BEFORE calculating term (prevents Inf/NaN from division)
        if (is.na(var_prop) || !is.finite(var_prop) || var_prop < 1e-10) {
            warning(sprintf("Phenotype '%s': var_prop invalid (%.2e), using fallback", p, var_prop))
            var_prop <- pooled_prop * (1 - pooled_prop) / 10 # Conservative fallback
        }

        term <- (pooled_prop * (1 - pooled_prop) / var_prop) - 1

        # Guard against negative or invalid term (if var is huge)
        if (is.na(term) || !is.finite(term) || term <= 0) term <- 1 # Fallback

        alpha_est <- pooled_prop * term
        beta_est <- (1 - pooled_prop) * term

        # Store I² for downstream use (avoids re-running meta-analysis)
        i_squared <- m$I2

        yield_params[[p]] <- list(
            dist = "beta",
            params = c(shape1 = alpha_est, shape2 = beta_est),
            mean = pooled_prop,
            ci_lower = plogis(m$lower.random),
            ci_upper = plogis(m$upper.random),
            i_squared = i_squared,
            n_studies = nrow(d_sub)
        )

        yield_audit_rows[[p]] <- data.frame(
            Phenotype = p,
            Pooled_Yield = sprintf("%.1f%%", pooled_prop * 100),
            CI_95_Lower = sprintf("%.1f%%", plogis(m$lower.random) * 100),
            CI_95_Upper = sprintf("%.1f%%", plogis(m$upper.random) * 100),
            N_Studies = nrow(d_sub),
            Alpha = round(alpha_est, 1),
            Beta = round(beta_est, 1)
        )
    }
}

yield_audit_df <- bind_rows(yield_audit_rows)

# ==============================================================================
# 3. Joint (Gene, Variant_Class) Proportions for Multinomial Sampling
# ==============================================================================
cat("3. Computing Joint (Gene, Variant_Class) Proportions...\n")

# Prepare weighted counts (default to 1 if missing)
variant <- variant %>%
    mutate(weight = replace_na(n_probands_with_variant, 1))

# Clean inheritance for metadata attachment
variant <- variant %>%
    mutate(inheritance_clean = ifelse(is.na(inheritance_mode), "Unknown", inheritance_mode))

# Separate difficult locus variants (handled distinctly in Section 4)
# Update: Capture BOTH 1 (PKD1) and 2 (MUC1) as difficult
variant_standard <- variant %>%
    filter(is.na(difficult_locus_flag) | !difficult_locus_flag %in% c(1, 2))

variant_difficult <- variant %>%
    filter(!is.na(difficult_locus_flag) & difficult_locus_flag %in% c(1, 2))

cat("   -> Standard variants:", nrow(variant_standard), "\n")
cat("   -> Difficult locus variants:", nrow(variant_difficult), "\n")

# Audit Exports
AUDIT_DIR <- "data/audit"
if (!dir.exists(AUDIT_DIR)) dir.create(AUDIT_DIR, recursive = TRUE)
write_csv(variant_standard, file.path(AUDIT_DIR, "02_variant_standard.csv"))
write_csv(variant_difficult, file.path(AUDIT_DIR, "02_variant_difficult.csv"))
cat("   -> Audit CSVs exported to", AUDIT_DIR, "\n")

# Compute modal inheritance per gene (for metadata attachment)
gene_inheritance_mode <- variant_standard %>%
    group_by(gene, inheritance_clean) %>%
    summarise(inh_count = sum(weight), .groups = "drop") %>%
    group_by(gene) %>%
    slice_max(inh_count, n = 1, with_ties = FALSE) %>%
    select(gene, modal_inheritance = inheritance_clean)

# Joint distribution: (Gene, Variant_Class) per Phenotype
phenos_var <- unique(variant_standard$phenotype)
gene_variant_joint <- list()

for (p in phenos_var) {
    v_sub <- variant_standard %>% filter(phenotype == p)

    # Step 1: Aggregate counts by (gene, variant_class)
    joint_raw <- v_sub %>%
        group_by(gene, variant_class) %>%
        summarise(count = sum(weight), .groups = "drop")

    # Step 2: Compute total and proportions
    total_count <- sum(joint_raw$count)
    joint_raw <- joint_raw %>%
        mutate(proportion = count / total_count) %>%
        arrange(desc(proportion))

    # Step 3: Attach modal inheritance as metadata
    joint_final <- joint_raw %>%
        left_join(gene_inheritance_mode, by = "gene") %>%
        mutate(modal_inheritance = ifelse(is.na(modal_inheritance), "Mixed", modal_inheritance))

    gene_variant_joint[[p]] <- joint_final

    cat(sprintf("   [%s] %d gene-variant pairs\n", p, nrow(joint_final)))
}

# ==============================================================================
# 4. Difficult Loci (Separate Treatment)
# ==============================================================================
cat("4. Modeling Difficult Loci (Separate)...\n")

difficult_loci <- variant_difficult %>%
    group_by(gene) %>%
    summarise(
        n_variants = n(),
        total_weight = sum(weight),
        primary_variant_class = names(which.max(table(variant_class))),
        modal_inheritance = names(which.max(table(inheritance_clean))),
        .groups = "drop"
    )

# Also compute difficult loci proportions by phenotype for conditional sampling
difficult_by_phenotype <- variant_difficult %>%
    group_by(phenotype, gene, variant_class) %>%
    summarise(count = sum(weight), .groups = "drop") %>%
    group_by(phenotype) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()

# ==============================================================================
# 5. Save Outputs
# ==============================================================================
cat("5. Saving Parameters...\n")

# Object Construction (new structure)
final_params <- list(
    phenotype_alphas = phenotype_alphas,
    yield_betas = yield_params,
    gene_variant_joint = gene_variant_joint,
    gene_inheritance_mode = gene_inheritance_mode,
    difficult_loci = difficult_loci,
    difficult_by_phenotype = difficult_by_phenotype
)

# Save each parameter as a separate object
cat("Saving individual parameters to", OUTPUT_PARAMS, "...\n")
for (param_name in names(final_params)) {
    fpath <- file.path(OUTPUT_PARAMS, paste0(param_name, ".rds"))
    saveRDS(final_params[[param_name]], fpath)
    cat("  - Saved:", fpath, "\n")
}

# Timestamp for the run
saveRDS(Sys.time(), file.path(OUTPUT_PARAMS, "timestamp.rds"))

# Audit File
sink(OUTPUT_AUDIT)
cat("=== RENAL GENETICS CEA: PARAMETER AUDIT ===\n")
cat("Generated:", as.character(Sys.time()), "\n\n")

cat("--- 1. DIAGNOSTIC YIELDS (Pooled Random Effects) ---\n")
print(yield_audit_df, row.names = FALSE)

cat("\n--- 2. PHENOTYPE PREVALENCE (Dirichlet Alphas) ---\n")
print(as.data.frame(phenotype_alphas))

cat("\n--- 3. JOINT (GENE, VARIANT_CLASS) PROPORTIONS BY PHENOTYPE ---\n")
cat("Proportions should sum to 1.0.\n\n")
for (p in names(gene_variant_joint)) {
    cat(paste0(
        "[", toupper(p), "] (Sum = ",
        round(sum(gene_variant_joint[[p]]$proportion), 4), ")\n"
    ))
    # Show top 10 pairs
    top_pairs <- head(gene_variant_joint[[p]], 10)
    print(as.data.frame(top_pairs), row.names = FALSE)
    cat("\n")
}

cat("\n--- 4. DIFFICULT LOCI (Separate Treatment) ---\n")
if (nrow(difficult_loci) > 0) {
    print(as.data.frame(difficult_loci), row.names = FALSE)
    cat("\nDifficult Loci by Phenotype:\n")
    if (nrow(difficult_by_phenotype) > 0) {
        print(as.data.frame(difficult_by_phenotype), row.names = FALSE)
    }
} else {
    cat("No difficult loci flagged in input dataset.\n")
}

sink()

cat("SUCCESS: Parameters saved to", OUTPUT_PARAMS, "\n")
cat("AUDIT: Written to", OUTPUT_AUDIT, "\n")
