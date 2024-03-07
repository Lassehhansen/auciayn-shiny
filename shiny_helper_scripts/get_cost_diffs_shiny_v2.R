### Filtering AUROC and AUPRC values with relevance
library(tidyverse)
results_df = read_csv("simulation_data/results_simulated_data_v8.csv")
# Filter out models that don't meet the criteria

calculate_overall_prevalence <- function(prevalence_a, prevalence_b, attribute_ratio) {
  R <- attribute_ratio
  overall_prevalence <- (prevalence_a * R / (R + 1)) + (prevalence_b * 1 / (R + 1))
  return(overall_prevalence)
}


filtered_results_df <- results_df %>%
  mutate(
    overall_prevalence = calculate_overall_prevalence(p_pos_space, p_pos_space, attribute_ratio)
  ) %>% 
  filter(auroc_total < 0.90, auprc_total < 0.90) %>% 
  filter(auroc_total > 0.5, auprc_total > overall_prevalence)

### write/read csv


library(readr)



### Creating auroc auprc ratio


filtered_results_df$auroc_prc_ratio_a = filtered_results_df$auroc_a/filtered_results_df$auprc_a
filtered_results_df$auroc_prc_ratio_b = filtered_results_df$auprc_a/filtered_results_df$auprc_b
filtered_results_df$auroc_prc_ratio_total = filtered_results_df$auroc_total/filtered_results_df$auprc_total


## Plot 0

tpr <- function(tau, pdf_pos) {
  result <- 1 - pbeta(tau, pdf_pos$alpha, pdf_pos$beta)
  if (any(is.nan(result)) || any(is.infinite(result))) {
    cat("Invalid TPR: alpha =", pdf_pos$alpha, "beta =", pdf_pos$beta, "tau =", tau, "\n")
  }
  return(result)
}

fpr <- function(tau, pdf_neg) {
  result <- 1 - pbeta(tau, pdf_neg$alpha, pdf_neg$beta)
  if (any(is.nan(result)) || any(is.infinite(result))) {
    cat("Invalid FPR: alpha =", pdf_neg$alpha, "beta =", pdf_neg$beta, "tau =", tau, "\n")
  }
  return(result)
}


# Positive Predictive Value
ppv <- function(tau, scenario) {
  if (tau == 0) {
    return(0)  # PPV is not defined for tau = 0, set to 0 or another appropriate value
  } else {
    tpr_val <- tpr(tau, scenario$pdf_pos)
    fpr_val <- fpr(tau, scenario$pdf_neg)
    numerator <- scenario$p_pos * tpr_val
    denominator <- scenario$p_pos * tpr_val + scenario$p_neg * fpr_val
    return(ifelse(denominator == 0, 1, numerator / denominator))  # Avoid division by zero
  }
}


## Plot 1:

### Defining different types of expected cost functions:

min_cost_thresh <- function(pdf_pos, pdf_neg, prevalence, cost_fp, cost_fn, n) {
  tau <- seq(0, 1, length.out = 100)  # Thresholds
  N <- n  # Total number of samples (adjust as needed)
  
  # Calculate expected FPs and FNs for each threshold
  expected_fps <- sapply(tau, function(t) N * (1 - prevalence) * fpr(t, pdf_neg))
  expected_fns <- sapply(tau, function(t) N * prevalence * (1 - tpr(t, pdf_pos)))
  
  # Choose a threshold that minimizes the expected cost
  min_cost_threshold <- tau[which.min((cost_fp * expected_fps + cost_fn * expected_fns) / N)]
  
  return(min_cost_threshold)
}

library(ggplot2)
library(dplyr)
library(pracma) 


set.seed(1997)
cost_threshold_data <- function(alpha_pos, beta_pos, alpha_neg, beta_neg, p_pos_space, cost_fp, cost_fn, n) {
  tau <- round(seq(0, 1, length.out = 100), digits = 2)
  N <- n  # Total number of samples
  
  # Ensure all parameters are positive
  if(alpha_pos <= 0 || beta_pos <= 0 || alpha_neg <= 0 || beta_neg <= 0) {
    stop("Alpha and beta parameters must be positive")
  }
  
  # Calculate TPR and FPR using the CDF (pbeta) instead of integrating the PDF
  tpr_values <- sapply(tau, function(t) 1 - pbeta(t, alpha_pos, beta_pos))
  fpr_values <- sapply(tau, function(t) 1 - pbeta(t, alpha_neg, beta_neg))
  
  # Calculate expected false positives and false negatives for each threshold
  expected_fps <- N * (1 - p_pos_space) * fpr_values
  expected_fns <- N * p_pos_space * (1 - tpr_values)
  
  # Calculate the expected cost
  expected_costs <- (cost_fp * expected_fps + cost_fn * expected_fns) / N
  
  # Create a data frame for plotting or further analysis
  plot_data <- data.frame(threshold = tau, expected_cost = expected_costs)
  return(plot_data)
}

### Choose only the top 1

filtered_results_df$ev_pos_a = round(filtered_results_df$ev_pos_a, digits = 2)
filtered_results_df$ev_pos_b = round(filtered_results_df$ev_pos_b, digits = 2)

best_auprc_model <- filtered_results_df %>%
  group_by(p_pos_space, attribute_ratio, ev_pos_a, ev_pos_b) %>% 
  filter(auprc_total == max(auprc_total)) %>%
  slice(1) %>% 
  mutate(model_choice = "AUPRC")

best_auroc_model <- filtered_results_df %>%
  group_by(p_pos_space, attribute_ratio, ev_pos_a, ev_pos_b,) %>% 
  filter(auroc_total == max(auroc_total)) %>%
  slice(1) %>% 
  mutate(model_choice = "AUROC")

best_auroc_auprc_models = rbind(best_auroc_model, best_auprc_model)

best_auroc_auprc_models = best_auroc_auprc_models %>% filter(ev_pos_a > 0.56,
                                                             ev_pos_b > 0.56)

# Create the cost grid
mesh_density <- 5
cost_fn_space <- seq(from = 0, to = 8, by = 2)
cost_fp_space <- seq(from = 0, to = 8, by = 2)
model_choice = c("AUROC", "AUPRC")
subgroup = c("A", "B")
# Create an empty data frame to store the cost grid

cost_grid = expand.grid(cost_fn = cost_fn_space, 
                        cost_fp = cost_fp_space, 
                        model_choice = model_choice,
                        ev_pos_a = unique(best_auroc_auprc_models$ev_pos_a),
                        ev_pos_b = unique(best_auroc_auprc_models$ev_pos_b),
                        attribute_ratio = unique(best_auroc_auprc_models$attribute_ratio),
                        p_pos_space = unique(best_auroc_auprc_models$p_pos_space),
                        Subgroup = subgroup)



# Left join the best model data onto the cost grid
cost_grid_with_models <- cost_grid %>%
  left_join(best_auroc_auprc_models, by = c("ev_pos_a",
                                            "ev_pos_b",
                                            "model_choice",
                                            "p_pos_space", 
                                            "attribute_ratio"))


# Function to apply cost_threshold_data and return a data frame with additional model information
get_expected_costs <- function(row, n) {
  # Select the right subgroup based on the 'Subgroup' column
  subgroup_suffix <- ifelse(row$Subgroup == "A", "_a", "_b")
  alpha_pos <- row[[paste0("alpha_pos", subgroup_suffix)]]
  beta_pos <- row[[paste0("beta_pos", subgroup_suffix)]]
  alpha_neg <- row[[paste0("alpha_neg", subgroup_suffix)]]
  beta_neg <- row[[paste0("beta_neg", subgroup_suffix)]]
  p_pos_space <- row[["p_pos_space"]]
  
  # Calculate the expected costs
  expected_costs_df <- cost_threshold_data(alpha_pos, beta_pos, alpha_neg, beta_neg, p_pos_space, row$cost_fp, row$cost_fn, n)
  
  # Add additional columns
  expected_costs_df$model_choice <- row$model_choice
  expected_costs_df$cost_fn <- round(row$cost_fn, digits = 2)
  expected_costs_df$cost_fp <- round(row$cost_fp, digits = 2)
  expected_costs_df$subgroup <- row$Subgroup
  expected_costs_df$p_pos <- row$p_pos_space
  expected_costs_df$attribute_ratio <- row$attribute_ratio
  expected_costs_df$alpha_pos_a <- row$alpha_pos_a
  expected_costs_df$alpha_pos_b <- row$alpha_pos_b
  expected_costs_df$alpha_neg_a <- row$alpha_neg_a
  expected_costs_df$alpha_neg_b <- row$alpha_neg_b
  
  expected_costs_df$beta_pos_a <- row$beta_pos_a
  expected_costs_df$beta_pos_b <- row$beta_pos_b
  expected_costs_df$beta_neg_a <- row$beta_neg_a
  expected_costs_df$beta_neg_b <- row$beta_neg_b
  expected_costs_df$ev_pos_a <-  row$ev_pos_a
  expected_costs_df$ev_pos_b <-  row$ev_pos_b
  return(expected_costs_df)
}

# Initialize an empty data frame to store all the results
all_expected_costs <- data.frame()

# Loop over each row of the cost grid and calculate expected costs
for (i in 1:nrow(cost_grid_with_models)) {
  row_costs <- get_expected_costs(cost_grid_with_models[i, ], n = 100)
  all_expected_costs <- rbind(all_expected_costs, row_costs)
}

write_csv(all_expected_costs, "simulation_data/thresholds_exp_costs_shiny_v2.csv")


