NO_STATE_CHANGE = 0
RECOVERED_OR_OOHQ = 16
MILD_OR_MODERATE = 23
SEVERE = 4
DECEASED = 5

TERMINAL_STATES = c(DECEASED)

MAX_PATH_LENGTH = 9


MILD_OR_MODERATE = 23

to_model_state = function(state) {
  if (state == MILD || state == MODERATE) {
    return(MILD_OR_MODERATE)
  } else {
    return(state)
  }
}


construct_multi_state_model_dataset = function() {
  df = read.csv('./malka_dataframe_right_before_analysis.csv')
  
  df = df[!is.na(df$age),]
  df = df[!is.na(df$sex),]
  
  # these 0 values are only used for initialization 
  # values for each state are updated according to the transition function when fitting the model.
  df$was_severe = 0
  df$cumulative_time = 0
  initial_covariates_df = as.data.frame(model.matrix(~age*sex+age*as.factor(type0)+age*cumulative_time+age*was_severe, df))
  
  covariate_names = tail(names(initial_covariates_df), n=-1) # exclude intercept generated from model.matrix
  non_covariate_column_names = c('source_state', 'target_state', 'tstop', 'tstart', 'id', 'first_date', 'hospital')
  
  df = cbind(initial_covariates_df[,covariate_names], df[,non_covariate_column_names])
  
  head(df, n=5)
  
  
  construct_one_object = function(id) {
    
    patient_data = df[df$id == id,]
    
    # construct the series of states visited:
    states = c(patient_data$source_state)
    last_state = tail(patient_data$target_state, n=1)
    if (last_state != NO_STATE_CHANGE) states = c(states, last_state)
    
    time_at_each_state = patient_data$tstop - patient_data$tstart
    
    # covariates at time of hospitalization:
    covariates = as.numeric(patient_data[1, covariate_names])
    
    object = list(
      covariates = covariates,
      states = states,
      time_at_each_state = time_at_each_state,
      id = id,
      
      date_of_hospitalization = as.Date(patient_data$first_date[1]),
      hospital = as.character(patient_data$hospital[1])
    )
    
    return(object)
  }
  
  
  dataset = lapply(unique(df$id), function(id) construct_one_object(id))
}


get_corona_model_covariate_names = function() {
  # df = read.csv('./malka_dataframe_right_before_analysis.csv')
  # 
  # df = df[!is.na(df$age),]
  # df = df[!is.na(df$sex),]
  # 
  # # these 0 values are only used for initialization
  # # values for each state are updated according to the transition function when fitting the model.
  # df$was_severe = 0
  # df$cumulative_time = 0
  # initial_covariates_df = as.data.frame(model.matrix(~age*sex+age*as.factor(type0)+age*cumulative_time+age*was_severe, df))
  # 
  # covariate_names = tail(names(initial_covariates_df), n=-1) # exclude intercept generated from model.matrix
  # 
  # return(covariate_names)
  return(readRDS('israeli_corona_model_covariate_names.RDS'))
}


AGE_COVARIATE_INDEX = 1
CUMULATIVE_TIME_COVARIATE_INDEX = 5
CUMULATIVE_TIME_INTERACTION_WITH_AGE_COVARIATE_INDEX = 10
PATIENT_WAS_SEVERE_COVARIATE_INDEX = 6
PATIENT_WAS_SEVERE_INTERACTION_WITH_AGE_COVARIATE_INDEX = 11


update_covariates = function(sample_covariates, origin_state, target_state, time_at_origin_state, absolute_time_of_entry_to_target_state=NULL) {
  new_covariates = sample_covariates
  
  # update 'patient was _severe'
  if (origin_state == SEVERE) {
    new_covariates[PATIENT_WAS_SEVERE_COVARIATE_INDEX] = 1
    new_covariates[PATIENT_WAS_SEVERE_INTERACTION_WITH_AGE_COVARIATE_INDEX] = sample_covariates[AGE_COVARIATE_INDEX]
  }
  
  # update cumulative time
  new_covariates[CUMULATIVE_TIME_COVARIATE_INDEX] = new_covariates[CUMULATIVE_TIME_COVARIATE_INDEX] + time_at_origin_state
  
  if (!is.null(absolute_time_of_entry_to_target_state) && 
      absolute_time_of_entry_to_target_state > new_covariates[CUMULATIVE_TIME_COVARIATE_INDEX]) {
    new_covariates[CUMULATIVE_TIME_COVARIATE_INDEX] = absolute_time_of_entry_to_target_state
  }
  
  
  new_covariates[CUMULATIVE_TIME_INTERACTION_WITH_AGE_COVARIATE_INDEX] = new_covariates[CUMULATIVE_TIME_COVARIATE_INDEX]*new_covariates[AGE_COVARIATE_INDEX]
  
  return(new_covariates)
}



probability_of_death = function(all_runs) {
  return(mean(sapply(all_runs, function(run) tail(run$states, n=1) == DECEASED)))
}

probability_of_severe = function(all_runs) {
  return(mean(sapply(all_runs, function(run) SEVERE %in% run$states)))
}


time_at_hospital = function(monte_carlo_run) {
  states = monte_carlo_run$states
  time_at_each_state = monte_carlo_run$time_at_each_state
  
  if (length(states) > length(time_at_each_state)){
    states = head(states, -1)
  }
  
  return(sum(time_at_each_state[states != RECOVERED_OR_OOHQ]))
}


truncated_mean_time_hospitalized = function(all_runs) {
  t = sapply(all_runs, function(run) time_at_hospital(run) )
  return(mean(t, trim=0.1))
}


QUANTILES = c(0.10,
              0.25,
              0.5,
              0.75,
              0.90)

quantiles_of_time_hospitalized = function(all_runs) {
  t = sapply(all_runs, function(run) time_at_hospital(run) )
  return(quantile(t, probs=QUANTILES))
}


quantile_of_time_hospitalized = function(all_runs, q) {
  t = sapply(all_runs, function(run) time_at_hospital(run) )
  return(quantile(t, probs=q))
}


time_at_severe = function(monte_carlo_run) {
  states = monte_carlo_run$states
  time_at_each_state = monte_carlo_run$time_at_each_state
  
  if (length(states) > length(time_at_each_state)){
    states = head(states, -1)
  }
  
  return(sum(time_at_each_state[states == SEVERE]))
}

truncated_mean_time_at_severe = function(all_runs) {
  t = sapply(all_runs, function(run) time_at_severe(run) )
  return(mean(t, trim=0.1))
}


quantile_of_time_at_severe = function(all_runs, q) {
  t = sapply(all_runs, function(run) time_at_severe(run) )
  if (length(t[t > 0]) == 0) return(0)
  return(quantile(t[t > 0], probs=q, na.rm = TRUE))
}

# if time at severe was 0, ie sample was not in severe - sample was excluded
quantiles_of_time_at_severe = function(all_runs) {
  t = sapply(all_runs, function(run) time_at_severe(run) )
  return(quantile(t[t > 0], probs=QUANTILES, na.rm = TRUE))
}



#---------------------------------------------------
#
# Construct a single patient's covariates
#
# EXAMPLE USE:
# male, age=55, type0=severe, cumulative_time=7
# cumulative_time = 7
# days_spent_in_mild_state = 3
# 
# covariates = construct_covariates(sex = MALE, # 1
#                                   age = 55,
#                                   entry_state = SEVERE, # 4
#                                   cumulative_time = cumulative_time,
#                                   was_severe = 1)
# 
# # current state=mild, spent 3 days in mild
# model$run_monte_carlo_simulation(covariates,
#                                  origin_state = MILD_OR_MODERATE, # 23
#                                  current_time = cumulative_time + days_spent_in_mild_state,
#                                  n_random_samples = 3,
#                                  max_transitions = MAX_TRANSITIONS)
#---------------------------------------------------

MALE = 1
FEMALE = 0

MILD = 2
MODERATE = 3
SEVERE = 4

construct_covariates = function(sex, entry_state, age, cumulative_time = 0, was_severe = 0) {
  covariates = rep(0, 11)
  
  covariates[1] = age
  covariates[2] = (sex == MALE)
  covariates[3] = if (entry_state == MODERATE) 1 else 0
  covariates[4] = if (entry_state == SEVERE) 1 else 0
  covariates[5] = cumulative_time # cumulative time
  covariates[6] = was_severe # was severe
  covariates[7] = age * sex
  covariates[8] = age * covariates[3] # age * started moderate
  covariates[9] = age * covariates[4] # age * started severe
  covariates[10] = age * cumulative_time
  covariates[11] = age * was_severe
  
  return(covariates)
}