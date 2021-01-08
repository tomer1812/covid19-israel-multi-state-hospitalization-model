
compute_statistics = function(monte_carlo_runs) {
  # Death probability
  p_death = probability_of_death(monte_carlo_runs)
  
  # The probability of being SEVERE. 
  p_severe = probability_of_severe(monte_carlo_runs)
  
  # Quantiles of length of hospitalization - 10%, 25%, 50%, 75%, 90%. 
  q_time_hospitalized = quantiles_of_time_hospitalized(monte_carlo_runs)
  
  # Quantiles of days at severe state among the MC visits at severe - 10%, 25%, 50%, 75%, 90%
  # if time at severe was 0, ie sample was not in severe - sample was excluded
  q_time_at_severe = quantiles_of_time_at_severe(monte_carlo_runs)
  
  next_row = data.frame('probability_of_death' = c(p_death), 
                        'probability_of_visit_to_severe_state' = c(p_severe))
  next_row[, c(paste('time_hospitalized_remaining_quantile', QUANTILES, sep='_'))] = c(q_time_hospitalized)
  next_row[, c(paste('time_remaining_in_severe_state_quantile', QUANTILES, sep='_'))] = c(q_time_at_severe)
  
  return(next_row)
}

compute_sd_using_bootstrap = function(weighted_bootstrap_models, 
                                      covariates, 
                                      current_state, 
                                      total_days_since_hospitalization, 
                                      m_monte_carlo_runs, 
                                      max_transitions) {
  bootstrap_df = data.frame()
  for (model in weighted_bootstrap_models) {
    
    monte_carlo_runs = model$run_monte_carlo_simulation(covariates, 
                                                        to_model_state(current_state),
                                                        current_time = total_days_since_hospitalization,
                                                        n_random_samples = m_monte_carlo_runs, 
                                                        max_transitions = max_transitions)
    
    bootstrap_df = rbind(bootstrap_df, compute_statistics(monte_carlo_runs))
  }
  
  return(lapply(bootstrap_df, function(col) sd(col)))
}


plot_cdf_of_remaining_time_in_hospital = function(monte_carlo_runs) {
  df = data.frame(time_hospitalized_per_run = sapply(monte_carlo_runs, 
                                                     function(run) sum(run$time_at_each_state)))
  ggplot(df, 
         aes(time_hospitalized_per_run)) + 
    stat_ecdf(geom = "step") +
    ggtitle("CDF of remaining time in hospital") + 
    xlab("t (days)") +
    ylab("F(t)")
}

