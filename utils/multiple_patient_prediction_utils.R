source('./utils/utils.R')

load_dataset = function(global_start, global_censor) {
  
  global_start = "2018-03-01" # arbitrarily set, to include as many patients as possible
  global_censor = as.Date("2020-05-03") - 1
  
  dataset = construct_multi_state_model_dataset()
  
  #Filter subject that don't have time spent in non terminal state
  fix_last_transition <- function(obj) {
    end_point = as.Date(obj$date_of_hospitalization) + sum(obj$time_at_each_state)
    if (end_point < global_censor) {
      if ((as.Date(obj$date_of_hospitalization) < as.Date(global_start)) || (length(obj$time_at_each_state) < length(obj$states) || ((as.Date(obj$date_of_hospitalization) + sum(obj$time_at_each_state) < global_censor))) && !tail(obj$states, n=1) %in% c(TERMINAL_STATES,RECOVERED_OR_OOHQ)) {
        return(NULL)
      }
    }
    n_states = length(obj$states)
    if (n_states > 1) {
      for (i in 1:(n_states - 1)) {
        ins = obj$states[i]
        outs = obj$states[i+1]
        for (bts in BAD_TRANSITIONS) {
          if ((ins == bts[1]) && (outs == bts[2])) {
            return(NULL)
          }
        }
      }
    }
    return(obj)
  }
  dataset = sapply(dataset, function(obj) fix_last_transition(obj))
  dataset = Filter(dataset, f = function(x) !is.null(x))
  return(dataset)
}


# Code used for fitting the model:
#
# covariate_names = get_corona_model_covariate_names()
# dataset = load_dataset(global_start, global_censor)
# model = IsraeliMultiStateModel()
# model$fit(dataset = dataset,
#           terminal_states = TERMINAL_STATES, 
#           update_covariates_function = update_covariates, 
#           covariate_names = covariate_names)
# 
# remove(dataset)
# 
# save.image('multiple_patient_prediction_model.Rdata')


# including
dates_between = function(first_date, last_date) {
  sapply(1:(as.numeric(as.Date(last_date) - as.Date(first_date)) + 1), 
         function(i) as.character(as.Date(first_date) + (i-1)))
}


# obj (input to the function) is some object with attributes:
#   $states 
#   $time_at_each_state
#   $t_start (the process described by states/time_at_each_state begins at t_start)
get_state_per_day = function(obj, first_date, last_date,Ndays=NULL) {
  if (is.null(Ndays)) {
    all_dates = dates_between(first_date, last_date)
    Ndays = length(all_dates)
  }
  day_diff = as.integer(first_date - as.Date(obj$t_start))
  state_per_day = c()
  
  # time_at_each_state can be null if start at recovered/oohq and stay there forever
  if (!is.null(obj$time_at_each_state)) {
    for (i in 1:length(obj$time_at_each_state)) {
      state_per_day = c(state_per_day, rep(obj$states[[i]], round(obj$time_at_each_state[[i]])))
    }
  }
  
  # As opposed to early stopping paths, which could end in a non-terminal state
  # full paths end in a terminal state and dont have a time_at_state value for this state
  # in such cases fill the rest of the time with the final state
  if (length(obj$states) > length(obj$time_at_each_state)) {
    
    last_state = tail(obj$states, n=1)
    #stopifnot(last_state %in% c(DECEASED, RECOVERED_OR_OOHQ))
    
    # it could be the case that we have already filled values up to last_date
    days_left_in_final_state_until_last_date = max(as.numeric((as.Date(last_date) - as.Date(obj$t_start))) - length(state_per_day) + 1, 0)
    state_per_day = c(state_per_day, 
                      rep(last_state, days_left_in_final_state_until_last_date))
  }
  ret = rep(NA,Ndays)
  all_dates = seq(1:Ndays)
  filled_dates = seq(1 - day_diff, length(state_per_day) - day_diff)
  mask = filled_dates[filled_dates %in% all_dates]
  #print(day_diff)
  #print(mask)
  ret[mask] = state_per_day[mask + day_diff]
  
  return(ret)
}



# returns a dataframe where:
# each row contains the patients state for the corresponding day (NULLs if no state)
#
# each object must have attributes as described in get_state_per_day()
construct_state_per_day_matrix = function(obj_list, first_date, last_date) {
  dates = dates_between(first_date, last_date)
  Ndays = length(dates)
  
  tmp = t(sapply(obj_list,function(obj) {return(c(get_state_per_day(obj, first_date, last_date, Ndays)))}))
  tmp = as.data.frame(tmp)
  colnames(tmp) = dates
  return(tmp)
}



# returns a dataframe such that:
# each row coresponds with a single run. 
# each row holds the number of hospitalized per day, for a certain run. 
#
# the obj_list is the main data_object list constructed in the beggining.
construct_n_hospitalized_per_run_matrix = function(obj_list, first_date, last_date, t_start = NULL,count_states=c(SEVERE, MILD_OR_MODERATE)) {
  
  dates = dates_between(first_date, last_date)
  n_hospitalized_per_run_matrix = data.frame()
  for (i in 1:N_MONTE_CARLO_RUNS) {
    
    # take the ith run sampled per each patient and construct a "patient style object"
    ith_run_obj_list = lapply(obj_list, function(obj) {
      return(list(
        states = obj$all_runs[[i]]$states,
        time_at_each_state = obj$all_runs[[i]]$time_at_each_state,
        t_start = if (is.null(t_start)) as.Date(obj$t_start) else max(as.Date(t_start),as.Date(obj$t_start))
      ))
    })
    
    # use this run to compute the number of hospitalized patients each day
    ith_run_state_per_day_matrix = construct_state_per_day_matrix(ith_run_obj_list, first_date, last_date)
    n_hospitalized_per_day = as.numeric(sapply(ith_run_state_per_day_matrix, function(col) sum(col %in% count_states)))
    
    
    next_row = as.data.frame(matrix(rep(NA, length(dates)), nrow=1))
    colnames(next_row) = dates
    next_row[1, dates] = c(n_hospitalized_per_day)
    
    n_hospitalized_per_run_matrix = rbind(n_hospitalized_per_run_matrix, next_row)
  }
  return(n_hospitalized_per_run_matrix)
}

# constract_cov <-function(age,sex,type) {
#   return(c(age,sex,type==3,type==4,0,0,age*sex,age*(type==3),age*(type==4),0,0))
# }

construct_patient_object = function(sex, 
                                    age, 
                                    state_at_hospitalization,
                                    date_of_hospitalization,
                                    current_state, # the state at time T 
                                    days_spent_in_current_state, # number of days already in this state
                                    was_severe, # at any time during hospitalization
                                    T
) {
  return(
    list(
      covariates_at_T = construct_covariates(sex = sex, 
                                             entry_state = state_at_hospitalization, 
                                             age = age, 
                                             # days since hospitalization at entry to last state:
                                             cumulative_time = as.numeric(T - date_of_hospitalization - days_spent_in_current_state),
                                             was_severe = if(was_severe) 1 else 0),
      state_at_T = to_model_state(current_state),
      days_since_hospitalization_at_T = as.numeric(T - date_of_hospitalization),
      date_of_hospitalization = date_of_hospitalization
    )
  )
}


construct_future_arrival_patient_object = function(sex, 
                                                   age, 
                                                   state_at_hospitalization,
                                                   date_of_hospitalization) {
  return(
    list(
      covariates_at_T = construct_covariates(sex = sex, 
                                             entry_state = state_at_hospitalization, 
                                             age = age, 
                                             # days since hospitalization at entry to last state:
                                             cumulative_time = 0,
                                             was_severe = 0),
      state_at_T = to_model_state(state_at_hospitalization),
      days_since_hospitalization_at_T = 0,
      date_of_hospitalization = as.Date(date_of_hospitalization)
    )
  )
}

run_monte_carlo_per_patient = function(patients, N_MONTE_CARLO_RUNS, MAX_PATH_LENGTH) {
  patients = lapply(patients, function(obj) {
    obj$all_runs = model$run_monte_carlo_simulation(obj$covariates_at_T,
                                                    origin_state = obj$state_at_T,
                                                    current_time = obj$days_since_hospitalization_at_T,
                                                    n_random_samples = N_MONTE_CARLO_RUNS,
                                                    max_transitions = MAX_PATH_LENGTH)
    obj$t_start = obj$date_of_hospitalization
    return(obj)
  })
  
  return(patients)
  
}


plot_expected_number_of_patients_each_day = function(patients, # make sure to run: patients = run_monte_carlo_per_patient(patients)
                                                     first_date,
                                                     last_date,
                                                     states_to_include = c(SEVERE, MILD_OR_MODERATE) # default is all "hospitalized" states
) {
  
  n_hospitalized_per_run_matrix = construct_n_hospitalized_per_run_matrix(patients,
                                                                          first_date = first_date,
                                                                          last_date = last_date, 
                                                                          t_start=first_date,
                                                                          count_states=states_to_include)
  
  estimated_mean_n_hospitalized = sapply(n_hospitalized_per_run_matrix, mean)
  
  estimated_q75_n_hospitalized =  sapply(n_hospitalized_per_run_matrix, function(col) quantile(col, probs=0.75))
  estimated_q25_n_hospitalized =  sapply(n_hospitalized_per_run_matrix, function(col) quantile(col, probs=0.25))
  
  estimated_q90_n_hospitalized =  sapply(n_hospitalized_per_run_matrix, function(col) quantile(col, probs=0.9))
  estimated_q10_n_hospitalized =  sapply(n_hospitalized_per_run_matrix, function(col) quantile(col, probs=0.1))
  
  
  dates = dates_between(first_date, last_date)
  
  df = data.frame(mean = estimated_mean_n_hospitalized,
                  estimated_q75_n_hospitalized = estimated_q75_n_hospitalized,
                  estimated_q25_n_hospitalized = estimated_q25_n_hospitalized,
                  
                  estimated_q90_n_hospitalized = estimated_q90_n_hospitalized,
                  estimated_q10_n_hospitalized = estimated_q10_n_hospitalized,
                  
                  date=dates,
                  metric=rep('n patients', length(dates)))
  
  ggplot(df, aes(date, mean)) +
    geom_crossbar(aes(ymin=estimated_q10_n_hospitalized, ymax=estimated_q90_n_hospitalized), show.legend = FALSE, fill="#eaac8b") +
    geom_crossbar(aes(ymin=estimated_q25_n_hospitalized, ymax=estimated_q75_n_hospitalized), show.legend = FALSE, fill="#e56b6f") +
    # ggtitle('n Hospitalized Patients Per Day') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ylab("n patients")  
}


generate_random_hospitalized_patients_df = function(n_hospitalized_patients, T, PREDICT_N_DAYS_AHEAD) {
  # random sex: male or female
  sex = sample(c(MALE, FEMALE), 
               size = n_hospitalized_patients, 
               replace = TRUE)
  
  # random age: uniformally chosen between 10 and 90
  age = 10 + round(runif(n = n_hospitalized_patients)*80)
  
  # some list of between 1 to 5 states:
  states = lapply(1:n_hospitalized_patients, function(i) {
    sample(c(MILD, MODERATE, SEVERE), size = sample(1:4), replace = TRUE)
  })
  
  # corresponding number of days spent at each state
  time_at_each_state = lapply(states, function(patient_states) {
    1 + round(runif(length(patient_states))*4)
  })
  
  # the first state:
  state_at_hospitalization = sapply(states, function(patient_states) {
    patient_states[[1]]
  })
  
  # T minus the total time passed in the path:
  date_of_hospitalization = as.Date(sapply(time_at_each_state, function(patient_times) {
    as.character(T - sum(patient_times))
  }))
  
  # last state:
  current_state = sapply(states, function(patient_states) {
    as.vector(tail(patient_states, n=1))
  })
  
  # last time:
  days_spent_in_current_state = sapply(time_at_each_state, function(patient_times) {
    as.vector(tail(patient_times, n=1))
  })
  
  # was severe before the current state:
  was_severe = sapply(states, function(patient_states) {
    SEVERE %in% patient_states[-length(patient_states)]
  })
  
  
  
  hospitalized_patients_df = data.frame(
    sex = sex,
    age = age,
    state_at_hospitalization = state_at_hospitalization,
    date_of_hospitalization = date_of_hospitalization,
    current_state = current_state,
    days_spent_in_current_state = days_spent_in_current_state,
    was_severe = was_severe
  )
  
  return(hospitalized_patients_df)
}


generate_future_arrival_patients_df = function(n_future_arrivals, T, PREDICT_N_DAYS_AHEAD) {
  # random sex: male or female
  sex = sample(c(MALE, FEMALE), size = n_future_arrivals, replace = TRUE)
  
  # random age: uniformally chosen between 10 and 90
  age = 10 + round(runif(n = n_future_arrivals)*80)
  
  # random state: 
  state_at_hospitalization = sample(c(MILD, MODERATE, SEVERE), size = n_future_arrivals, replace = TRUE)
  
  # random arrival date: chosen from between T and T+PREDICT_N_DAYS_AHEAD
  date_of_hospitalization = T + round(runif(n = n_future_arrivals)*PREDICT_N_DAYS_AHEAD)
  
  future_arrival_patients_df = data.frame(
    sex = sex,
    age = age,
    state_at_hospitalization = state_at_hospitalization,
    date_of_hospitalization = date_of_hospitalization
  )
  
  
}


construct_expected_deaths_table = function(patients, T, PREDICT_N_DAYS_AHEAD){
  n_dead_per_run = construct_n_hospitalized_per_run_matrix(patients,
                                                           first_date = T,
                                                           last_date = T+PREDICT_N_DAYS_AHEAD, 
                                                           t_start=T,
                                                           count_states=c(DECEASED))
  
  every_fifth_column_idx = (1:(length(n_dead_per_run)/5))*5
  n_dead_on_fifth_day = n_dead_per_run[,every_fifth_column_idx]
  expected_cumulative_number_of_deaths = colMeans(n_dead_on_fifth_day)
  
  return(data.frame(expected_cumulative_number_of_deaths = expected_cumulative_number_of_deaths,
                    expected_number_of_deaths_in_last_5_days = diff(c(0,as.numeric(expected_cumulative_number_of_deaths)))))
}




# Below you can see the expected proportion of deaths by days since hospitalization. 
# That is, we are counting deaths among all patients 5 days since hospitalization, 10 days hospitalization etc.; averaging over all sampled monte carlo paths.
expected_deaths_by_days_since_hospitalization = function(patients, T, PREDICT_N_DAYS_AHEAD) {
  
  expected_deaths_df = data.frame(days_since_hospitalization = seq(5, PREDICT_N_DAYS_AHEAD, by=5))
  expected_deaths = c()
  
  for (n_days in expected_deaths_df$days_since_hospitalization) {
    
    number_of_deaths_before_n_days_since_hospitalization_per_run = sapply(1:N_MONTE_CARLO_RUNS, function(run_i) {
      
      # count the number of patients dead by n_days in run_i
      number_of_deaths_before_n_days_since_hospitalization = sum(sapply(patients, function(patient) {
        
        patient_run_i = patient$all_runs[[run_i]]
        time_of_run = sum(patient_run_i$time_at_each_state)
        
        patient_dies_before_n_days = DECEASED %in% patient_run_i$states && # patient dies
          time_of_run < n_days && # dies within n_days
          patient$date_of_hospitalization + time_of_run < T + PREDICT_N_DAYS_AHEAD # dies within prediction period
        
        return(patient_dies_before_n_days)
      }))
      
      return(number_of_deaths_before_n_days_since_hospitalization)
    })
    
    expected_deaths_by_n_days = sum(number_of_deaths_before_n_days_since_hospitalization_per_run) / N_MONTE_CARLO_RUNS
    expected_deaths = c(expected_deaths, expected_deaths_by_n_days)
  }
  
  expected_deaths_df$n_expected_deaths = as.vector(expected_deaths)
  expected_deaths_df$expected_proportion_of_deaths = as.vector(expected_deaths) / length(patients)
  return(expected_deaths_df)
}