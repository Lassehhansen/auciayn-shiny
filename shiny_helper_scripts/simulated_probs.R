library(purrr)
library(dplyr)

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


set.seed(1997)


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

model_choice = c("AUROC", "AUPRC")
cost_grid <- expand.grid(
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



convert_ev_confidence_to_alpha_beta <- function(ev, confidence) {
  alpha <- ev * confidence
  beta <- (1 - ev) * confidence
  return(list(alpha = alpha, beta = beta))
}

# Function to simulate a subgroup
simulate_subgroup <- function(total_samples, alpha_pos, beta_pos, alpha_neg, beta_neg, prevalence) {
  
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
  
  attribute_ratio = rep(attribute_ratio, total_samples) 
  ev_pos_a = rep(ev_pos_a, total_samples)
  ev_neg_a = rep(ev_neg_a, total_samples)
  ev_pos_b = rep(ev_pos_b, total_samples)
  ev_neg_b = rep(ev_neg_b, total_samples)
  confidence_pos_a = rep(confidence_pos_a, total_samples)
  confidence_neg_a = rep(confidence_neg_a, total_samples)
  confidence_pos_b = rep(confidence_pos_b, total_samples) 
  confidence_neg_b = rep(confidence_neg_b, total_samples)
  p_pos_space = rep(p_pos_space, total_samples)
  
  dataset_total <- data.frame(
    probability = data_combined,
    subgroup = subgroup_combined,
    labels_combined = labels_combined,
    attribute_ratio = attribute_ratio, 
    ev_pos_a = ev_pos_a, 
    ev_neg_a = ev_neg_a, 
    ev_pos_b = ev_pos_b, 
    ev_neg_b = ev_neg_b,
    confidence_pos_a = confidence_pos_a, 
    confidence_neg_a = confidence_neg_a, 
    confidence_pos_b = confidence_pos_b, 
    confidence_neg_b = confidence_neg_b,
    p_pos_space = p_pos_space
  )
  
  return(dataset_total)
}

# Use apply to iterate over the rows of the dataframe
total_samples <- 5000

cost_grid_merge_auroc = cost_grid_merge %>% filter(model_choice == "AUROC")
cost_grid_merge_auprc = cost_grid_merge %>% filter(model_choice == "AUPRC")

# Use apply to iterate over the rows of the dataframe
simulation_results_auroc <- apply(cost_grid_merge_auroc, 1, function(row) {
  sample_population(
    total_samples = total_samples,
    attribute_ratio = as.numeric(row['attribute_ratio']),
    ev_pos_a = as.numeric(row['ev_pos_a']),
    ev_neg_a = as.numeric(row['ev_neg_a']),
    ev_pos_b = as.numeric(row['ev_pos_b']),
    ev_neg_b = as.numeric(row['ev_neg_b']),
    confidence_pos_a = as.numeric(row['confidence_pos_a']),
    confidence_neg_a = as.numeric(row['confidence_neg_a']),
    confidence_pos_b = as.numeric(row['confidence_pos_b']),
    confidence_neg_b = as.numeric(row['confidence_neg_b']),
    p_pos_space = as.numeric(row['p_pos_space'])
  )
})

simulation_results_auprc <- apply(cost_grid_merge_auprc, 1, function(row) {
  sample_population(
    total_samples = total_samples,
    attribute_ratio = as.numeric(row['attribute_ratio']),
    ev_pos_a = as.numeric(row['ev_pos_a']),
    ev_neg_a = as.numeric(row['ev_neg_a']),
    ev_pos_b = as.numeric(row['ev_pos_b']),
    ev_neg_b = as.numeric(row['ev_neg_b']),
    confidence_pos_a = as.numeric(row['confidence_pos_a']),
    confidence_neg_a = as.numeric(row['confidence_neg_a']),
    confidence_pos_b = as.numeric(row['confidence_pos_b']),
    confidence_neg_b = as.numeric(row['confidence_neg_b']),
    p_pos_space = as.numeric(row['p_pos_space'])
  )
})

# Combine the list of data frames into a single data frame
simulation_results_df_auroc <- do.call(rbind, simulation_results_auroc)
simulation_results_df_auprc <- do.call(rbind, simulation_results_auprc)

simulation_results_df_auroc$model_choice = "AUROC"
simulation_results_df_auprc$model_choice = "AUPRC"

simulation_results_tot = rbind(simulation_results_df_auroc, simulation_results_df_auprc)

write_csv(simulation_results_tot, "simulation_data/simulated_probs_shiny.csv")

# Combine the list of data frames into a single data frame
