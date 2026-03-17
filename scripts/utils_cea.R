# utils_cea.R
# Purpose: Shared utility functions for Cost-Effectiveness Analysis
# Author: Renal Genetics CEA Team
# Date: 2025-12-30

library(dplyr)

#' Perform incremental cost-effectiveness analysis with dominance checking
#'
#' @param data Data frame with columns: strategy_label, cost, effect
#' @return Data frame with incremental analysis and dominance status
perform_incremental_analysis <- function(data) {
    # Ensure cost is numeric to avoid arranging issues
    data$cost <- as.numeric(data$cost)
    data$effect <- as.numeric(data$effect)

    # Sort by cost (ascending)
    data <- data %>% arrange(cost)

    n <- nrow(data)
    data$dominance_status <- "non_dominated"
    data$incremental_cost <- NA_real_
    data$incremental_effect <- NA_real_
    data$icer <- NA_real_

    # Step 1: Identify strictly dominated strategies
    # A strategy is strictly dominated if it costs more than another and has equal or lower effect
    for (i in 1:n) {
        for (j in 1:n) {
            if (i != j) {
                if (data$cost[i] >= data$cost[j] && data$effect[i] <= data$effect[j] &&
                    !(data$cost[i] == data$cost[j] && data$effect[i] == data$effect[j])) {
                    data$dominance_status[i] <- "strictly_dominated"
                    break
                }
            }
        }
    }

    # Step 2: Extended dominance among non-strictly-dominated strategies
    # Iteratively remove strategies with decreasing or equal ICER
    repeat {
        # Get non-dominated strategies
        frontier <- data %>%
            filter(dominance_status == "non_dominated") %>%
            arrange(cost, effect)

        if (nrow(frontier) <= 2) break

        # Calculate ICERs for frontier
        frontier$temp_incr_cost <- c(NA, diff(frontier$cost))
        frontier$temp_incr_effect <- c(NA, diff(frontier$effect))
        frontier$temp_icer <- frontier$temp_incr_cost / frontier$temp_incr_effect

        # Check for extended dominance (ICER decreasing)
        extended_dom_found <- FALSE
        for (i in 2:(nrow(frontier) - 1)) {
            if (!is.na(frontier$temp_icer[i]) && !is.na(frontier$temp_icer[i + 1])) {
                if (frontier$temp_icer[i + 1] < frontier$temp_icer[i]) {
                    # Strategy i is extendedly dominated
                    strat_label <- frontier$strategy_label[i]
                    data$dominance_status[data$strategy_label == strat_label] <- "extendedly_dominated"
                    extended_dom_found <- TRUE
                    break
                }
            }
        }

        if (!extended_dom_found) break
    }

    # Step 3: Calculate final incremental values for frontier
    frontier_final <- data %>%
        filter(dominance_status == "non_dominated") %>%
        arrange(cost)

    for (i in 1:nrow(frontier_final)) {
        idx <- which(data$strategy_label == frontier_final$strategy_label[i])

        if (i == 1) {
            # First frontier strategy has no previous comparator by spec.
            data$incremental_cost[idx] <- NA_real_
            data$incremental_effect[idx] <- NA_real_
        } else {
            # Incremental vs. previous non-dominated strategy
            data$incremental_cost[idx] <- frontier_final$cost[i] - frontier_final$cost[i - 1]
            data$incremental_effect[idx] <- frontier_final$effect[i] - frontier_final$effect[i - 1]
        }

        # Calculate ICER
        if (!is.na(data$incremental_effect[idx]) && data$incremental_effect[idx] > 0) {
            data$icer[idx] <- data$incremental_cost[idx] / data$incremental_effect[idx]
        }
    }

    return(data)
}

#' Compute Cost-Effectiveness Acceptability Curve (CEAC)
#'
#' @param ceac_input Data frame with iteration_id, strategy_id, strategy_label,
#'   total_cost_per_proband_cad, diagnoses_per_proband
#' @param wtp_range Numeric vector of WTP thresholds
#' @return Data frame with wtp_threshold_cad, strategy_id, strategy_label, pr_cost_effective
compute_ceac <- function(ceac_input, wtp_range) {
    strategies <- ceac_input %>%
        distinct(strategy_id, strategy_label) %>%
        arrange(strategy_id, strategy_label)

    n_iter <- dplyr::n_distinct(ceac_input$iteration_id)

    ceac_rows <- vector("list", length(wtp_range))
    for (k in seq_along(wtp_range)) {
        wtp <- wtp_range[k]
        best <- ceac_input %>%
            mutate(nmb = (wtp * diagnoses_per_proband) - total_cost_per_proband_cad) %>%
            group_by(iteration_id) %>%
            slice_max(nmb, n = 1, with_ties = FALSE) %>%
            ungroup()
        ceac_rows[[k]] <- best %>%
            count(strategy_id, strategy_label, name = "n_optimal") %>%
            mutate(
                wtp_threshold_cad = wtp,
                pr_cost_effective = n_optimal / n_iter
            ) %>%
            select(wtp_threshold_cad, strategy_id, strategy_label, pr_cost_effective)
    }

    ceac_raw <- bind_rows(ceac_rows)

    ceac_final <- tidyr::expand_grid(
        wtp_threshold_cad = wtp_range,
        strategy_id = strategies$strategy_id
    ) %>%
        left_join(strategies, by = "strategy_id") %>%
        left_join(ceac_raw, by = c("wtp_threshold_cad", "strategy_id", "strategy_label")) %>%
        mutate(pr_cost_effective = dplyr::coalesce(pr_cost_effective, 0)) %>%
        arrange(wtp_threshold_cad, strategy_id, strategy_label)

    # Integrity check: probabilities must sum to 1 at each WTP threshold
    prob_sums <- ceac_final %>%
        group_by(wtp_threshold_cad) %>%
        summarise(total = sum(pr_cost_effective), .groups = "drop")
    if (any(abs(prob_sums$total - 1) > 1e-6)) {
        stop("CEAC probabilities do not sum to 1 for one or more WTP thresholds")
    }

    ceac_final
}
