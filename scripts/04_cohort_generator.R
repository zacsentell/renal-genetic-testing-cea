# 04_cohort_generator.R
# Purpose: Generate a synthetic cohort of CKD probands with phenotypes,
#          monogenic status, joining causal variant architecture.
# Author: Renal Genetics CEA Team
# Date: 2025-12-24

# Set CRAN mirror for non-interactive sessions
local({
    r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org"
    options(repos = r)
})

if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("readr")) install.packages("readr")
if (!require("purrr")) install.packages("purrr")

library(dplyr)
library(tidyr)
library(readr)
library(purrr)

# ==============================================================================
# IO Paths
# ==============================================================================
PARAMS_DIR <- "data/params"
OUTPUT_RDS <- "data/intermediate/04_synthetic_cohort.rds"
OUTPUT_CSV <- "outputs/cohort/synthetic_cohort.csv"
OUTPUT_SUMMARY <- "outputs/cohort/cohort_summary.txt"

# Ensure output directories exist
if (!dir.exists(dirname(OUTPUT_RDS))) dir.create(dirname(OUTPUT_RDS), recursive = TRUE)
if (!dir.exists(dirname(OUTPUT_CSV))) dir.create(dirname(OUTPUT_CSV), recursive = TRUE)

cat("--- Starting Synthetic Cohort Generation ---\n")

# ==============================================================================
# 1. Load Parameters
# ==============================================================================
cat("Loading parameters from", PARAMS_DIR, "...\n")

read_param <- function(name) {
    path <- file.path(PARAMS_DIR, paste0(name, ".rds"))
    if (!file.exists(path)) stop("Parameter file not found: ", path)
    readRDS(path)
}

# Load all required parameters
params <- list(
    phenotype_alphas = read_param("phenotype_alphas"),
    yield_betas = read_param("yield_betas"),
    gene_variant_joint = read_param("gene_variant_joint"),
    difficult_loci = read_param("difficult_loci"),
    difficult_by_phenotype = read_param("difficult_by_phenotype")
    # timestamp = read_param("timestamp") # Optional, might not exist if 02 wasn't fully run with it
)

cat("Parameters loaded successfully.\n")

# ==============================================================================
# 2. Key Sampling Functions
# ==============================================================================

#' Generate random Dirichlet samples using Gamma distribution
#'
#' @param n Number of samples (NOTE: usually 1 for PSA outer loop, but here we
#'          might just want the probability vector)
#' @param alpha Vector of concentration parameters
#' @return A matrix where each row is a sample summing to 1
rdirichlet_custom <- function(n, alpha) {
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- rowSums(x)
    x / as.vector(sm)
}

#' Sample Phenotypes for N probands
#'
#' @param n_probands Integer, size of cohort
#' @param alphas Named numeric vector of Dirichlet alphas
#' @param seed Integer seed
#' @return List with phenotypes (character vector) and sampled_proportions (named numeric vector)
sample_phenotype <- function(n_probands, alphas, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    # Step 1: Sample ONE set of probabilities from Dirichlet (representing this cohort's true prevalence)
    # In a full PSA, this would be updated per outer iteration.
    # Here we sample the probabilities once for this specific cohort realization.
    probs <- rdirichlet_custom(1, alphas)[1, ]
    names(probs) <- names(alphas)

    # Step 2: Sample n_probands from these probability weights
    # We use sample() with replacement
    phenotypes <- sample(names(probs), n_probands, replace = TRUE, prob = probs)

    # Return both the phenotypes and the sampled proportions for parameter tracing
    list(
        phenotypes = phenotypes,
        sampled_proportions = probs
    )
}

#' Sample Monogenic Status
#'
#' @param phenotypes Character vector of assigned phenotypes
#' @param yield_betas List of beta parameters (shape1, shape2) by phenotype
#' @param sampled_yields Named numeric vector of pre-sampled yields by phenotype
#' @param seed Integer seed
#' @return Logical vector (TRUE = Monogenic)
sample_monogenic_status <- function(phenotypes, yield_betas, sampled_yields = NULL, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    # Get unique phenotypes to optimize
    unique_phenos <- unique(phenotypes)

    if (is.null(sampled_yields)) {
        # Sample the "True Yield" for this cohort/iteration for each phenotype.
        pheno_yields <- sapply(unique_phenos, function(p) {
            bp <- yield_betas[[p]]$params
            if (is.null(bp)) {
                return(0)
            }
            rbeta(1, shape1 = bp["shape1"], shape2 = bp["shape2"])
        })
    } else {
        # Use externally sampled yields to keep one stochastic truth per iteration.
        pheno_yields <- sampled_yields[unique_phenos]
        pheno_yields[is.na(pheno_yields)] <- 0
    }

    # Map back to full vector
    prob_vector <- pheno_yields[phenotypes]

    # Bernoulli trial for each patient
    runif(length(phenotypes)) < prob_vector
}

#' Harmonize Variant Class
#' Maps raw variant classes to detection categories
harmonize_variant_class <- function(raw_class, is_difficult) {
    if (is.null(raw_class) || is.na(raw_class)) {
        return(NA_character_)
    }

    # If it's a difficult locus gene, it is ALWAYS "Difficult_Locus" for detection purposes
    if (is_difficult) {
        return("Difficult_Locus")
    }

    # Map standard classes
    if (raw_class %in% c("SNV", "INDEL", "SNV_INDEL", "SNV_Indel")) {
        return("SNV_Indel")
    }
    if (raw_class %in% c("CNV", "SV", "CNV_SV")) {
        return("CNV_SV")
    }

    return(raw_class) # Fallback
}

#' Sample Causal Architecture
#'
#' @param phenotypes Vector of phenotypes
#' @param is_monogenic Vector of monogenic status
#' @param gene_variant_joint List of joint proportions
#' @param difficult_by_phenotype Tibble of difficult locus proportions
#' @param seed Integer seed
#' @return Data frame of architecture details
sample_causal_architecture <- function(phenotypes, is_monogenic,
                                       gene_variant_joint, difficult_by_phenotype,
                                       difficult_loci,
                                       seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    n <- length(phenotypes)

    # Initialize output vectors
    out_gene <- rep(NA_character_, n)
    out_variant_raw <- rep(NA_character_, n)
    out_inheritance <- rep(NA_character_, n)
    out_difficult_bool <- rep(FALSE, n)
    out_difficult_gene <- rep(NA_character_, n)

    # Process only monogenic cases
    indices <- which(is_monogenic)

    if (length(indices) == 0) {
        return(data.frame(
            causal_gene = out_gene,
            variant_class_raw = out_variant_raw,
            inheritance_mode = out_inheritance,
            is_difficult_locus = out_difficult_bool,
            difficult_locus_gene = out_difficult_gene,
            variant_class = rep(NA_character_, n),
            stringsAsFactors = FALSE
        ))
    }

    # Process by phenotype to vectorize sampling where possible
    # Or simply loop for clarity since N=1000 is small

    # Pre-calculate difficult vs standard probabilities per phenotype
    # Diff proportion = Sum of difficult locus proportions in difficult_by_phenotype
    phenos_with_diff <- unique(difficult_by_phenotype$phenotype)

    for (i in indices) {
        p <- phenotypes[i]

        # Check if this phenotype has difficult loci
        prob_diff <- 0

        # Get standard variant counts
        joint_sub <- gene_variant_joint[[p]]
        count_std <- if (is.null(joint_sub)) 0 else sum(joint_sub$count)

        # Get difficult variant counts
        count_diff <- 0
        diff_sub <- NULL

        if (p %in% phenos_with_diff) {
            diff_sub <- difficult_by_phenotype %>% filter(phenotype == p)
            count_diff <- sum(diff_sub$count)
        }

        # Total count for this phenotype
        total_variants <- count_std + count_diff

        if (total_variants > 0) {
            prob_diff <- count_diff / total_variants
        }

        # Decide: Difficult or Standard?
        is_diff <- runif(1) < prob_diff

        sampled_row <- NULL

        if (is_diff) {
            # Sample from Difficult table
            diff_sub <- difficult_by_phenotype %>% filter(phenotype == p)
            # Normalize proportions to sum to 1 within difficult set
            probs <- diff_sub$proportion / sum(diff_sub$proportion)
            chosen_idx <- sample(seq_len(nrow(diff_sub)), 1, prob = probs)
            sampled_row <- diff_sub[chosen_idx, ]

            out_gene[i] <- sampled_row$gene
            out_variant_raw[i] <- sampled_row$variant_class
            # Difficult table doesn't have inheritance, pull from difficult_loci summary or default to AD (PKD1/MUC1 mostly AD)
            # Looking at params$difficult_loci for modal inheritance
            mod_inh <- difficult_loci %>%
                filter(gene == sampled_row$gene) %>%
                pull(modal_inheritance)
            if (length(mod_inh) == 0) mod_inh <- "AD" # Fallback
            out_inheritance[i] <- mod_inh

            out_difficult_bool[i] <- TRUE
            out_difficult_gene[i] <- sampled_row$gene
        } else {
            # Sample from Standard Joint table
            joint_sub <- gene_variant_joint[[p]]
            if (is.null(joint_sub) || nrow(joint_sub) == 0) {
                # Fallback if no data (should not happen if parameterized correctly)
                out_gene[i] <- "Unknown_Gene"
                out_variant_raw[i] <- "SNV"
                out_inheritance[i] <- "AD"
            } else {
                # Sample
                probs <- joint_sub$proportion
                # normalize just in case (though they should sum to 1 approx)
                if (sum(probs) == 0) probs <- rep(1, length(probs))
                chosen_idx <- sample(seq_len(nrow(joint_sub)), 1, prob = probs)
                sampled_row <- joint_sub[chosen_idx, ]

                out_gene[i] <- sampled_row$gene
                out_variant_raw[i] <- sampled_row$variant_class
                out_inheritance[i] <- sampled_row$modal_inheritance
            }
        }
    }

    # Harmonize variant class
    out_variant_class <- mapply(harmonize_variant_class, out_variant_raw, out_difficult_bool)

    data.frame(
        causal_gene = out_gene,
        variant_class_raw = out_variant_raw,
        inheritance_mode = out_inheritance,
        is_difficult_locus = out_difficult_bool,
        difficult_locus_gene = out_difficult_gene,
        variant_class = out_variant_class,
        stringsAsFactors = FALSE
    )
}

#' Generate Full Cohort
#'
#' @param n_probands Integer
#' @param params List of parameter objects
#' @param sampled_yields Optional named numeric vector of per-phenotype yields
#' @param seed Master seed
generate_cohort <- function(n_probands, params, sampled_yields = NULL, seed = 12345) {
    set.seed(seed)

    # 1. Phenotypes (returns list with phenotypes and sampled_proportions)
    # Generate a seed for this step to ensure independence validity
    seed_pheno <- sample.int(1e7, 1)
    pheno_result <- sample_phenotype(n_probands, params$phenotype_alphas, seed = seed_pheno)
    phenotypes <- pheno_result$phenotypes
    sampled_proportions <- pheno_result$sampled_proportions

    # 2. Monogenic Status
    if (is.null(sampled_yields)) {
        sampled_yields <- sapply(names(params$phenotype_alphas), function(p) {
            bp <- params$yield_betas[[p]]$params
            if (is.null(bp)) {
                return(0)
            }
            rbeta(1, shape1 = bp["shape1"], shape2 = bp["shape2"])
        })
    }

    seed_mono <- sample.int(1e7, 1)
    is_monogenic <- sample_monogenic_status(
        phenotypes = phenotypes,
        yield_betas = params$yield_betas,
        sampled_yields = sampled_yields,
        seed = seed_mono
    )

    # 3. Causal Architecture
    seed_arch <- sample.int(1e7, 1)
    architecture <- sample_causal_architecture(
        phenotypes, is_monogenic,
        params$gene_variant_joint,
        params$difficult_by_phenotype,
        params$difficult_loci,
        seed = seed_arch
    )

    # Combine
    cohort <- data.frame(
        proband_id = 1:n_probands,
        phenotype = phenotypes,
        is_monogenic = is_monogenic,
        stringsAsFactors = FALSE
    ) %>%
        bind_cols(architecture)

    list(
        cohort = cohort,
        sampled_proportions = sampled_proportions, # For PRCC parameter tracing
        sampled_yields = sampled_yields, # Keep cohort and PRCC yield draws synchronized
        metadata = list(
            n = n_probands,
            seed = seed,
            timestamp = Sys.time()
        )
    )
}

# ==============================================================================
# 3. Main Execution
# ==============================================================================

if (sys.nframe() == 0) {
    # User configuration (could be moved to config file later)
    N_PROBANDS <- 1000
    MASTER_SEED <- 999

    cat(sprintf("Generating cohort of %d probands (Seed: %d)...\n", N_PROBANDS, MASTER_SEED))

    results <- generate_cohort(N_PROBANDS, params, seed = MASTER_SEED)
    cohort <- results$cohort

    # ==============================================================================
    # 4. Validation & Export
    # ==============================================================================

    cat("Cohort generation complete. Running validation checks...\n")

    sink(OUTPUT_SUMMARY)
    cat("=== SYNTHETIC COHORT SUMMARY ===\n")
    cat("Generated:", as.character(results$metadata$timestamp), "\n")
    cat("N Probands:", N_PROBANDS, "\n")
    cat("Seed:", MASTER_SEED, "\n\n")

    # 1. Phenotype Distribution
    cat("--- Phenotype Distribution ---\n")
    pheno_stats <- cohort %>%
        count(phenotype) %>%
        mutate(percent = n / sum(n) * 100)
    print(pheno_stats)
    cat("\nCompare to Dirichlet Alphas (Approx):\n")
    print(params$phenotype_alphas / sum(params$phenotype_alphas) * 100)

    # 2. Monogenic Yield by Phenotype
    cat("\n--- Monogenic Yield by Phenotype ---\n")
    yield_stats <- cohort %>%
        group_by(phenotype) %>%
        summarise(
            n_total = n(),
            n_monogenic = sum(is_monogenic),
            yield_percent = mean(is_monogenic) * 100
        )
    print(yield_stats)

    cat("\n--- Variant Class Rules (Monogenic Only) ---\n")
    # Check how many difficult loci were flagged per detection category
    class_stats <- cohort %>%
        filter(is_monogenic) %>%
        group_by(phenotype, variant_class) %>%
        summarise(n = n(), .groups = "drop") %>%
        pivot_wider(names_from = variant_class, values_from = n, values_fill = 0)
    print(class_stats)

    # 3. Difficult Locus Check
    cat("\n--- Difficult Locus Assignments ---\n")
    diff_check <- cohort %>%
        filter(is_difficult_locus) %>%
        count(phenotype, difficult_locus_gene, variant_class)
    print(diff_check)


    # 4. Top 20 Diagnostic Genes
    cat("\n--- Top 20 Diagnostic Genes ---\n")
    top_genes <- cohort %>%
        filter(is_monogenic) %>%
        group_by(causal_gene, phenotype, inheritance_mode, variant_class) %>%
        summarise(n = n(), .groups = "drop") %>%
        mutate(percent = n / sum(n) * 100) %>%
        arrange(desc(n)) %>%
        head(20)

    # Print as a nice data frame
    print(as.data.frame(top_genes))

    sink()

    # Save objects
    saveRDS(results, OUTPUT_RDS)
    write_csv(cohort, OUTPUT_CSV)

    cat("\nSUCCESS: Cohort generated.\n")
    cat("- RDS:", OUTPUT_RDS, "\n")
    cat("- CSV:", OUTPUT_CSV, "\n")
    cat("- Summary:", OUTPUT_SUMMARY, "\n")
}
