library(purrr)
library(dplyr)
library(yardstick)

set.seed(1997)

# Define the parameter space for expected value and confidence
mesh_density <- 4 # Adjust mesh density as needed
ev_space_pos <- seq(0.55, 0.80, length.out = mesh_density)  # Example range for expected values
ev_space_neg <- seq(0.20, 0.45, length.out = mesh_density)  # Example range for expected values
confidence_space <- seq(0.1, 4, length.out = 3)  # Example range for confidence
p_pos_space <- seq(0.05, 0.50, by = 0.15000)
attribute_ratio <- seq(0.05, 0.25, by = 0.10000)
total_samples <- 10000



# Create parameter grid with constraints on expected values
parameter_grid <- expand.grid(
  ev_pos_a = ev_space_pos,
  ev_pos_b = ev_space_pos,
  confidence_pos_a = confidence_space,
  confidence_pos_b = confidence_space,
  confidence_neg_a = confidence_space,
  confidence_neg_b = confidence_space,
  ev_neg_a = ev_space_neg,
  ev_neg_b = ev_space_neg,
  p_pos_space = p_pos_space,
  attribute_ratio = attribute_ratio,
  total_samples = total_samples
)

# Apply constraint: ev_pos >= ev_neg
parameter_grid <- parameter_grid[parameter_grid$ev_pos_a >= parameter_grid$ev_neg_a &
                                   parameter_grid$ev_pos_b >= parameter_grid$ev_neg_b, ]


convert_ev_confidence_to_alpha_beta <- function(ev, confidence) {
  alpha <- ev * confidence
  beta <- (1 - ev) * confidence
  return(list(alpha = alpha, beta = beta))
}

# Function to simulate a subgroup
simulate_subgroup <- function(total_samples, ev_pos, ev_neg, confidence_pos, confidence_neg, prevalence) {
  
  pdf_pos_params <- convert_ev_confidence_to_alpha_beta(ev_pos, confidence_pos)
  pdf_neg_params <- convert_ev_confidence_to_alpha_beta(ev_neg, confidence_neg)
  
  alpha_pos = pdf_pos_params$alpha
  beta_pos = pdf_pos_params$beta
  
  alpha_neg = pdf_neg_params$alpha
  beta_neg = pdf_neg_params$beta
  
  n_pos <- round(total_samples * prevalence)
  n_neg <- total_samples - n_pos
  
  # Simulate positives
  positives <- rbeta(n_pos, alpha_pos, beta_pos)
  
  # Simulate negatives
  negatives <- rbeta(n_neg, alpha_neg, beta_neg)
  
  probability <- c(positives, negatives)
  labels <- c(rep(1, n_pos), rep(0, n_neg))
  
  return(list(probability = probability, labels = labels))
}

calculate_metrics <- function(dataset) {
  # Ensure labels are a factor and the levels are correctly ordered to reflect the event of interest
  dataset$labels <- factor(dataset$labels, levels = c(1, 0))
  
  # Calculate AUROC
  auroc_value <- yardstick::roc_auc(dataset, labels, probability)$.estimate
  
  # Calculate AUPRC
  auprc_value <- yardstick::pr_auc(dataset, labels, probability)$.estimate
  
  brier_value <- yardstick::brier_class(dataset, labels, probability)$.estimate
  
  return(list(auroc = auroc_value, auprc = auprc_value, brier = brier_value))
}

sample_population <- function(total_samples, attribute_ratio, 
                              ev_pos_a, ev_neg_a, ev_pos_b, ev_neg_b, 
                              confidence_pos_a, confidence_neg_a, confidence_pos_b, confidence_neg_b, 
                              p_pos_space) {
  
  # Total samples for simulation
  total_samples <- total_samples # Adjust as needed
  
  # Calculate samples for each subgroup based on the specified ratio
  total_samples_a <- round(total_samples * attribute_ratio)
  total_samples_b <- total_samples - total_samples_a
  
  # Simulate for ethnic Danes
  sim_a <- simulate_subgroup(total_samples_a, ev_pos_a, ev_neg_a, confidence_pos_a, confidence_neg_a, p_pos_space)
  
  # Simulate for non-ethnic Danes
  sim_b <- simulate_subgroup(total_samples_b, ev_pos_b, ev_neg_b, confidence_pos_b, confidence_neg_b, p_pos_space)
  
  # Merge the simulated data
  data_combined <- c(sim_a$probability, sim_b$probability)
  labels_combined <- c(sim_a$labels, sim_b$labels)
  subgroup_combined <- c(rep("a", length(sim_a$probability)), rep("b", length(sim_b$probability)))
  
  dataset_total <- data.frame(
    probability = data_combined,
    labels = as.factor(labels_combined),
    subgroup = subgroup_combined
  )
  
  return(dataset_total)
}


# Adjusted get_metrics function to include parameters
get_metrics <- function(ev_pos_a, ev_pos_b, confidence_pos_a, confidence_pos_b, 
                        confidence_neg_a, confidence_neg_b, ev_neg_a, ev_neg_b, 
                        p_pos_space, attribute_ratio, total_samples) {
  # Convert expected values and confidence to alpha and beta
  pdf_pos_a_params <- convert_ev_confidence_to_alpha_beta(ev_pos_a, confidence_pos_a)
  pdf_neg_a_params <- convert_ev_confidence_to_alpha_beta(ev_neg_a, confidence_neg_a)
  
  pdf_pos_b_params <- convert_ev_confidence_to_alpha_beta(ev_pos_b, confidence_pos_b)
  pdf_neg_b_params <- convert_ev_confidence_to_alpha_beta(ev_neg_b, confidence_neg_b)
  
  # Sample population
  dataset_total <- sample_population(total_samples, attribute_ratio, 
                                     ev_pos_a, ev_neg_a, ev_pos_b, ev_neg_b, 
                                     confidence_pos_a, confidence_neg_a, confidence_pos_b, confidence_neg_b, 
                                     p_pos_space)
  
  # Calculate metrics
  dataset_a <- filter(dataset_total, subgroup == "a")
  dataset_b <- filter(dataset_total, subgroup == "b")
  
  metrics_a <- calculate_metrics(dataset_a)
  metrics_b <- calculate_metrics(dataset_b)
  metrics_total <- calculate_metrics(dataset_total)
  
  # Create a dataframe for results including parameters
  data <- data.frame(
    auroc_a = metrics_a$auroc,
    auprc_a = metrics_a$auprc,
    auroc_b = metrics_b$auroc,
    auprc_b = metrics_b$auprc,
    auroc_total = metrics_total$auroc,
    auprc_total = metrics_total$auprc,
    brier_a = metrics_a$brier,
    brier_b = metrics_b$brier,
    brier_total = metrics_total$brier,
    
    alpha_pos_a = pdf_pos_a_params$alpha,
    beta_pos_a = pdf_pos_a_params$beta,
    alpha_neg_a = pdf_neg_a_params$alpha,
    beta_neg_a = pdf_neg_a_params$beta,
    alpha_pos_b = pdf_pos_b_params$alpha,
    beta_pos_b = pdf_pos_b_params$beta,
    alpha_neg_b = pdf_neg_b_params$alpha,
    beta_neg_b = pdf_neg_b_params$beta,
    p_pos_space = p_pos_space,
    attribute_ratio = attribute_ratio,
    total_samples, total_samples,
    ev_pos_a = ev_pos_a, 
    ev_pos_b = ev_pos_b, 
    confidence_pos_a = confidence_pos_a, 
    confidence_pos_b = confidence_pos_b,
    ev_neg_a = ev_neg_a, 
    ev_neg_b = ev_neg_b, 
    confidence_neg_a = confidence_neg_a, 
    confidence_neg_b = confidence_neg_b
  )
  
  return(data)
}

library(furrr)
plan(multisession)

# Apply 'get_metrics' across the parameter grid using furrr::future_pmap
results_list <- future_pmap(parameter_grid, get_metrics,
                            .progress = T)

# Combine results into a single dataframe
results_df <- bind_rows(results_list)
library(tidyverse)
write_csv(results_df, "simulation_data/results_simulated_data_v8.csv")

