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

set.seed(1997)
cost_threshold_data <- function(alpha_pos, beta_pos, alpha_neg, beta_neg, p_pos_space, cost_fp, cost_fn, n, tau_max) {
  tau <- seq(0, tau_max, length.out = 100)  # Thresholds
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


### Calculate auec

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


calculate_auec <- function(alpha_pos, beta_pos, alpha_neg, beta_neg, p_pos_space, cost_fp, cost_fn, n, tau_max) {
  # Ensure the input is a data frame with the expected structure
  cost_data <- cost_threshold_data(alpha_pos, beta_pos, alpha_neg, beta_neg, p_pos_space, cost_fp, cost_fn, n, tau_max)
  
  # Sort the data by threshold to ensure the integration is done correctly
  cost_data <- cost_data[order(cost_data$threshold), ]
  
  # Calculate the area under the curve using the trapezoidal rule
  auc <- pracma::trapz(cost_data$threshold, cost_data$expected_cost)
  #auc2 <- DescTools::AUC(cost_data$threshold, cost_data$expected_cost)
  
  return(auc)
}


n = 1000
# Calculate expected cost for the best AUROC and AUPRC models across the cost grid

### Something like this:

# Prepare a grid for cost_fn and cost_fp combinations
mesh_density <- 5
cost_fn_space <- seq(from = 0, to = 8, by = 2)
cost_fp_space <- seq(from = 0, to = 8, by = 2)
model_choice = c("AUROC", "AUPRC")
# Create an empty data frame to store the cost grid

#testingcost_grrid = expand.grid(cost_fn = cost_fn_space, cost_fp = cost_fp_space)
#testingcost_grrid$ratio = testingcost_grrid$cost_fn/testingcost_grrid$cost_fp
cost_grid <- expand.grid(cost_fn = cost_fn_space, 
                         cost_fp = cost_fp_space,
                         ev_pos_a = unique(best_auroc_auprc_models$ev_pos_a),
                         ev_pos_b = unique(best_auroc_auprc_models$ev_pos_b),
                         p_pos_space = unique(filtered_results_df$p_pos_space),
                         attribute_ratio = unique(filtered_results_df$attribute_ratio),
                         model_choice = model_choice)



cost_grid_merge = left_join(cost_grid, best_auroc_auprc_models, by = c(
  "ev_pos_a",
  "ev_pos_b",
  "p_pos_space",
  "attribute_ratio",
  "model_choice"))

total_samples = 1000
cost_grid_merge <- cost_grid_merge %>%
  mutate(
    auec_a_10 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_a, 
                       beta_pos = beta_pos_a, 
                       alpha_neg = alpha_neg_a, 
                       beta_neg = beta_neg_a, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.1),
    auec_b_10 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_b, 
                       beta_pos = beta_pos_b, 
                       alpha_neg = alpha_neg_b, 
                       beta_neg = beta_neg_b, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.1),
    
    auec_a_20 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_a, 
                       beta_pos = beta_pos_a, 
                       alpha_neg = alpha_neg_a, 
                       beta_neg = beta_neg_a, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.2),
    auec_b_20 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_b, 
                       beta_pos = beta_pos_b, 
                       alpha_neg = alpha_neg_b, 
                       beta_neg = beta_neg_b, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.2),
    
    auec_a_30 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_a, 
                       beta_pos = beta_pos_a, 
                       alpha_neg = alpha_neg_a, 
                       beta_neg = beta_neg_a, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.3),
    auec_b_30 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_b, 
                       beta_pos = beta_pos_b, 
                       alpha_neg = alpha_neg_b, 
                       beta_neg = beta_neg_b, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.3),
    
    auec_a_40 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_a, 
                       beta_pos = beta_pos_a, 
                       alpha_neg = alpha_neg_a, 
                       beta_neg = beta_neg_a, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.4),
    auec_b_40 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_b, 
                       beta_pos = beta_pos_b, 
                       alpha_neg = alpha_neg_b, 
                       beta_neg = beta_neg_b, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.4),
    auec_a_50 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_a, 
                       beta_pos = beta_pos_a, 
                       alpha_neg = alpha_neg_a, 
                       beta_neg = beta_neg_a, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.5),
    auec_b_50 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_b, 
                       beta_pos = beta_pos_b, 
                       alpha_neg = alpha_neg_b, 
                       beta_neg = beta_neg_b, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.5),
    
    auec_a_60 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_a, 
                       beta_pos = beta_pos_a, 
                       alpha_neg = alpha_neg_a, 
                       beta_neg = beta_neg_a, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.6),
    auec_b_60 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_b, 
                       beta_pos = beta_pos_b, 
                       alpha_neg = alpha_neg_b, 
                       beta_neg = beta_neg_b, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.6),
    
    auec_a_70 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_a, 
                       beta_pos = beta_pos_a, 
                       alpha_neg = alpha_neg_a, 
                       beta_neg = beta_neg_a, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.7),
    auec_b_70 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_b, 
                       beta_pos = beta_pos_b, 
                       alpha_neg = alpha_neg_b, 
                       beta_neg = beta_neg_b, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.7),
    
    auec_a_80 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_a, 
                       beta_pos = beta_pos_a, 
                       alpha_neg = alpha_neg_a, 
                       beta_neg = beta_neg_a, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.8),
    auec_b_80 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_b, 
                       beta_pos = beta_pos_b, 
                       alpha_neg = alpha_neg_b, 
                       beta_neg = beta_neg_b, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.8),
    auec_a_90 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_a, 
                       beta_pos = beta_pos_a, 
                       alpha_neg = alpha_neg_a, 
                       beta_neg = beta_neg_a, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.9),
    auec_b_90 = mapply(calculate_auec, 
                       alpha_pos = alpha_pos_b, 
                       beta_pos = beta_pos_b, 
                       alpha_neg = alpha_neg_b, 
                       beta_neg = beta_neg_b, 
                       p_pos_space = p_pos_space, 
                       cost_fp = cost_fp, 
                       cost_fn = cost_fn,
                       n = 1000,
                       tau_max = 0.9),
    auec_a_100 = mapply(calculate_auec, 
                        alpha_pos = alpha_pos_a, 
                        beta_pos = beta_pos_a, 
                        alpha_neg = alpha_neg_a, 
                        beta_neg = beta_neg_a, 
                        p_pos_space = p_pos_space, 
                        cost_fp = cost_fp, 
                        cost_fn = cost_fn,
                        n = 1000,
                        tau_max = 1),
    auec_b_100 = mapply(calculate_auec, 
                        alpha_pos = alpha_pos_b, 
                        beta_pos = beta_pos_b, 
                        alpha_neg = alpha_neg_b, 
                        beta_neg = beta_neg_b, 
                        p_pos_space = p_pos_space, 
                        cost_fp = cost_fp, 
                        cost_fn = cost_fn,
                        n = 1000,
                        tau_max = 1)
    
  )

cost_grid_merge$cost_ratio_space = cost_grid_merge$cost_fn/cost_grid_merge$cost_fp
library(tidyr)
# Pivot the data longer for each tau_max setting
library(tidyverse)
cost_grid_merge <- cost_grid_merge %>%
  mutate(
    across(starts_with("auec_"), ~round(., 3))  # Assuming you want to round all auec values
  )

plot_new_wide <- cost_grid_merge %>%
  select(cost_fn, cost_fp, attribute_ratio, p_pos_space, model_choice, ev_pos_a, ev_pos_b, 37:57) %>% 
  pivot_wider(
    id_cols = c(cost_fn, cost_fp, attribute_ratio, p_pos_space, ev_pos_a, ev_pos_b),
    names_from = model_choice,
    values_from = starts_with("auec")
  )

#write_csv(plot_new_wide, "simulation_data/auec_values_shiny.csv")


# Now, calculate the differences after grouping
plot_new_wide <- plot_new_wide %>%
  mutate(
    auec_mean_diff_10_a = auec_a_10_AUROC - auec_a_10_AUPRC,
    auec_mean_diff_20_a = auec_a_20_AUROC - auec_a_20_AUPRC,
    auec_mean_diff_30_a = auec_a_30_AUROC - auec_a_30_AUPRC,
    auec_mean_diff_40_a = auec_a_40_AUROC - auec_a_40_AUPRC,
    auec_mean_diff_50_a = auec_a_50_AUROC - auec_a_50_AUPRC,
    
    auec_mean_diff_60_a = auec_a_60_AUROC - auec_a_60_AUPRC,
    auec_mean_diff_70_a = auec_a_70_AUROC - auec_a_70_AUPRC,
    auec_mean_diff_80_a = auec_a_80_AUROC - auec_a_80_AUPRC,
    auec_mean_diff_90_a = auec_a_90_AUROC - auec_a_90_AUPRC,
    auec_mean_diff_100_a = auec_a_100_AUROC - auec_a_100_AUPRC,
    
    
    auec_mean_diff_10_b = auec_b_10_AUROC - auec_b_10_AUPRC,
    auec_mean_diff_20_b = auec_b_20_AUROC - auec_b_20_AUPRC,
    auec_mean_diff_30_b = auec_b_30_AUROC - auec_b_30_AUPRC,
    auec_mean_diff_40_b = auec_b_40_AUROC - auec_b_40_AUPRC,
    auec_mean_diff_50_b = auec_b_50_AUROC - auec_b_50_AUPRC,
    auec_mean_diff_60_b = auec_b_60_AUROC - auec_b_60_AUPRC,
    auec_mean_diff_70_b = auec_b_70_AUROC - auec_b_70_AUPRC,
    auec_mean_diff_80_b = auec_b_80_AUROC - auec_b_80_AUPRC,
    auec_mean_diff_90_b = auec_b_90_AUROC - auec_b_90_AUPRC,
    
    auec_mean_diff_100_b = auec_b_100_AUROC - auec_b_100_AUPRC
    
  )

write_csv(plot_new_wide, "simulation_data/auec_data_shiny_v3.csv")

