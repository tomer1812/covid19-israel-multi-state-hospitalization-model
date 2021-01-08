source('./utils/wave2_utils.R')
library(parallel)
# ----------------------------------------------------------------------------------------

# Construct Group 1 and 2 Data: Positive With/without History

# ----------------------------------------------------------------------------------------


generate_random_positive_with_history_df = function(n, T) {
  
  is_male = sample(c(1,0), size = n, replace = TRUE)
  age = 10 + round(runif(n = n)*80) # random int from 10 to 90
  sector = sample(c("unkown", "rest", "mixed", "haredim"), size = n, replace = TRUE)
  first_positive_test_date = sample(T - c(1:20), size=n, replace = TRUE)  
  
  return(
    data.frame(
      is_male=is_male,
      age=age,
      sector=sector,
      first_positive_test_date=first_positive_test_date
    )
  )
}


generate_random_positive_with_NO_history_df = function(n, T, predict_n_days_ahead) {
  
  is_male = sample(c(1,0), size = n, replace = TRUE)
  age = 10 + round(runif(n = n)*80) # random int from 10 to 90
  sector = sample(c("unkown", "rest", "mixed", "haredim"), size = n, replace = TRUE)
  first_positive_test_date = sample(T + c(1:predict_n_days_ahead-1), size=n, replace = TRUE)  
  
  return(
    data.frame(
      is_male=is_male,
      age=age,
      sector=sector,
      first_positive_test_date=first_positive_test_date
    )
  )
}


# ----------------------------------------------------------------------------------------

# Construct Groups 3 and 4 Data: Hospitalized with/without History

# ----------------------------------------------------------------------------------------


# not full / correct etc. Simply a stub
generate_hospitalized_with_history_df = function(n, T) {
  # random sex: male or female
  is_male = sample(c(1,0), size=n, replace = TRUE)
  
  # random age: uniformally chosen between 10 and 90
  age = 10 + round(runif(n=n)*80)
  
  # some list of between 1 to 5 states:
  states = lapply(1:n, function(i) {
    sample(c(MILD_OR_MODERATE, SEVERE, CRITICAL), size = sample(1:4), replace = TRUE)
  })
  
  # corresponding number of days spent at each state
  time_at_each_state = lapply(states, function(patient_states) {
    1 + round(runif(length(patient_states))*4)
  })
  
  # the first state:
  state_at_hospitalization = sapply(states, function(patient_states) {
    first_state = patient_states[[1]]
    
    if (first_state == MILD_OR_MODERATE) {
      return(sample(c(MILD, MODERATE), size=1))
    } else {
      return(first_state)
    }
  })
  
  # T minus the total time passed in the path:
  date_of_hospitalization = as.Date(sapply(time_at_each_state, function(patient_times) {
    as.character(T - sum(patient_times))
  }))
  
  # last state:
  state_at_t_start = sapply(states, function(patient_states) {
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
  
  # was severe before the current state:
  was_critical = sapply(states, function(patient_states) {
    CRITICAL %in% patient_states[-length(patient_states)]
  })
  
  
  hospitalized_patients_df = data.frame(
    is_male = is_male,
    age = age,
    state_at_hospitalization = state_at_hospitalization,
    date_of_hospitalization = date_of_hospitalization,
    state_at_t_start = state_at_t_start,
    days_spent_in_current_state = days_spent_in_current_state,
    was_severe = was_severe,
    was_critical = was_critical
  )
  
  return(hospitalized_patients_df)
}

generate_hospitalized_with_NO_history_df = function(n=30, T=T, predict_n_days_ahead=PREDICT_N_DAYS_AHEAD) {
  # random sex: male or female
  is_male = sample(c(1,0), size=n, replace = TRUE)
  
  # random age: uniformally chosen between 10 and 90
  age = 10 + round(runif(n=n)*80)
  
  # random state: 
  state_at_hospitalization = sample(c(MILD, MODERATE, SEVERE, CRITICAL), size=n, replace = TRUE)
  
  # random arrival date: chosen from between T and T+PREDICT_N_DAYS_AHEAD
  date_of_hospitalization = T + round(runif(n=n)*predict_n_days_ahead)
  
  future_arrival_patients_df = data.frame(
    is_male=is_male,
    age = age,
    state_at_hospitalization = state_at_hospitalization,
    date_of_hospitalization = date_of_hospitalization
  ) 
}
# ----------------------------------------------------------------------------------------

# Perform the Monte-Carlo Sampling:

# ----------------------------------------------------------------------------------------

run_monte_carlo_per_patient = function(patients, n_monte_carlo_runs, max_path_length, n_cores=1) {
  patients = lapply(patients, function(obj) {
    obj$all_runs = model$run_monte_carlo_simulation(obj$covariates_at_t_start,
                                                    origin_state = obj$state_at_t_start,
                                                    current_time = obj$cumulative_time_at_t_start,
                                                    n_random_samples = n_monte_carlo_runs,
                                                    max_transitions = max_path_length, n_cores=n_cores)
    return(obj)
  })
  
  return(patients)
}

# ----------------------------------------------------------------------------------------

# Analysis of Monte-Carlo Runs:

# ----------------------------------------------------------------------------------------

# obj (input to the function) is some object with attributes:
#   $states 
#   $time_at_each_state
#   $t_start (the process described by states/time_at_each_state begins at t_start)
# 
# last_date included
get_state_per_day = function(obj, first_date, last_date) {
  
  stopifnot(last_date > first_date)
  stopifnot(obj$t_start >= first_date)
  stopifnot(obj$t_start <= last_date)
  
  n_days = as.integer(last_date - first_date) + 1
  
  # 1. fill days until t_start with NA
  t_start_offset_from_first_date = as.integer(obj$t_start - first_date)
  states_per_day = c(rep(NA, t_start_offset_from_first_date))
  
  # 2. add non-terminal states (have a given duration)
  n_times = length(obj$time_at_each_state)
  if (n_times > 0) {
    states_per_day_from_t_start = unlist(sapply(c(1:n_times), function(i) {
      rep(obj$states[[i]], round(obj$time_at_each_state[[i]]) )
    }))
    states_per_day = c(states_per_day, states_per_day_from_t_start)    
  }
  
  # can stop if already filled n_days:
  if (length(states_per_day) >= n_days) return(states_per_day[1:n_days])
  
  # 3. If path ends with a "forever" state, fill to the end with that state:
  if (length(obj$states) > length(obj$time_at_each_state)) {
    
    last_state = tail(obj$states, n=1)  
    stopifnot(last_state %in% c(TERMINAL_STATES, DISCHARGED))
    
    states_per_day = c(states_per_day, rep(last_state, n_days - length(states_per_day)))
  }
  
  return(states_per_day[1:n_days])
}



# --- TESTS: ---

test_get_state_per_day = function(verbose=TRUE) {
  
  
  first_date = as.Date("2020-10-1")
  last_date = as.Date("2020-10-10")
  n_days = as.integer(last_date - first_date)
  n_days
  
  
  test_cases = list(
    list(
      obj = list( # starts on first day, last state is a "forever state" before end
        states = c(MILD_OR_MODERATE, SEVERE, DISCHARGED),
        time_at_each_state = c(3,2),
        t_start = first_date
      ),
      expected_output = c(rep(MILD_OR_MODERATE, 3), rep(SEVERE,2), rep(DISCHARGED, 5))
    ),
    list(
      obj = list( # starts on first day, last state is a "forever state" after end
        states = c(MILD_OR_MODERATE, SEVERE, DISCHARGED),
        time_at_each_state = c(5,6),
        t_start = first_date
      ),
      expected_output = c(rep(MILD_OR_MODERATE, 5), rep(SEVERE,5))
    ),
    list(
      obj = list( # starts on 3rd day, last state is a "forever state" before end
        states = c(MILD_OR_MODERATE, SEVERE, DISCHARGED),
        time_at_each_state = c(2,3),
        t_start = first_date + 2
      ),
      expected_output = c(rep(NA,2), rep(MILD_OR_MODERATE, 2), rep(SEVERE,3), rep(DISCHARGED, 3))
    ),
    list(
      obj = list( # starts on 3rd day, last state is a "forever state" after end
        states = c(MILD_OR_MODERATE, SEVERE, DISCHARGED),
        time_at_each_state = c(5,6),
        t_start = first_date + 2
      ),
      expected_output = c(rep(NA,2), rep(MILD_OR_MODERATE, 5), rep(SEVERE,3))
    ),
    list(
      obj = list( # no last terminal state (might occur in an "early stop" path), starts no first day, stops before end
        states = c(MILD_OR_MODERATE, SEVERE),
        time_at_each_state = c(3,2),
        t_start = first_date
      ),
      expected_output = c(rep(MILD_OR_MODERATE, 3), rep(SEVERE,2), rep(NA, 5))
    ),
    list(
      obj = list( # no last terminal state (might occur in an "early stop" path), starts on 3rd day, stops before end
        states = c(MILD_OR_MODERATE, SEVERE),
        time_at_each_state = c(3,2),
        t_start = first_date + 2
      ),
      expected_output = c(rep(NA, 2), rep(MILD_OR_MODERATE, 3), rep(SEVERE,2), rep(NA, 3))
    ),
    list(
      obj = list( # starts on day 1 at DISCHARGED and remains DISCHARGED forever
        states = c(DISCHARGED),
        time_at_each_state = NULL,
        t_start = first_date
      ),
      expected_output = c(rep(DISCHARGED,10))
    ),
    list(
      obj = list( # starts on day 3 at DISCHARGED and remains DISCHARGED forever
        states = c(DISCHARGED),
        time_at_each_state = NULL,
        t_start = first_date + 2
      ),
      expected_output = c(rep(NA,2), rep(DISCHARGED,8))
    )
  )
  
  for (test_case in test_cases) {
    
    output = get_state_per_day(test_case$obj, first_date, last_date)
    if (verbose) {
      print(test_case)
      print(output)
      print(
        (output == test_case$expected_output) |
          (is.na(output) & is.na(test_case$expected_output))
      )
      old_output = old_get_state_per_day(test_case$obj, first_date, last_date-1)
      print(old_output)
      print(output)
    }
    
    stopifnot( # either equal, or both na; in all locations
      (output == test_case$expected_output) |
        (is.na(output) == is.na(test_case$expected_output))
    )
  }
  
}

test_get_state_per_day(verbose = FALSE)

# including last date:
dates_between = function(first_date, last_date) {
  sapply(1:(as.numeric(as.Date(last_date) - as.Date(first_date)) + 1), 
         function(i) as.character(as.Date(first_date) + (i-1)))
}

# returns a dataframe where:
# each row contains the patients state for the corresponding day (NULLs if no state)
#
# each object must have attributes as described in get_state_per_day()
#
# includes last date
         
# construct_state_per_day_matrix = function(obj_list, first_date, last_date) {
  
#   dates = dates_between(first_date, last_date)
#   state_per_day_matrix = data.frame(matrix(nrow=length(obj_list), ncol=length(dates)))
#   colnames(state_per_day_matrix) = dates
  
  
#   for (i in 1:length(obj_list)) {
#     state_per_day_matrix[i,] = c(get_state_per_day(obj_list[[i]], first_date, last_date))
#   }
  
#   return(state_per_day_matrix)
# }

construct_state_per_day_matrix = function(obj_list, first_date, last_date) {
  
  mat = t(sapply(obj_list, function(obj){
    return(
      c(get_state_per_day(obj, first_date, last_date))
      )
    }))
  
  state_per_day_matrix = as.data.frame(mat)
  
  return(state_per_day_matrix)
}
         

# returns a dataframe such that:
# each row coresponds with a single run. 
# each row holds the number of hospitalized per day, for a certain run. 
construct_n_hospitalized_per_run_matrix = function(patients_path, output_dir, first_date, last_date, count_states, n_monte_carlo_runs, n_cores=1) {
  
  dates = dates_between(first_date, last_date)
  #n_hospitalized_per_run_matrix = data.frame(matrix(nrow=n_monte_carlo_runs, ncol=length(dates)))
  #colnames(n_hospitalized_per_run_matrix) = dates
  
  #for (i in 1:n_monte_carlo_runs) {
  n_hospitalized_per_run_matrix = mcmapply(function(i) {
    patient_objects = readRDS(sprintf("%s/%d.R",patients_path,i))
    # take the ith run sampled per each patient and construct a "patient style object"
    ith_run_obj_list = lapply(patient_objects, function(obj) {
      return(list(
        states = obj$all_runs[[1]]$states,
        time_at_each_state = obj$all_runs[[1]]$time_at_each_state,
        t_start = as.Date(obj$t_start)
      ))
    })
    
    # use this run to compute the number of hospitalized patients each day
    ith_run_state_per_day_matrix = construct_state_per_day_matrix(ith_run_obj_list, first_date, last_date)
    write.csv(ith_run_state_per_day_matrix, sprintf("%s/%d_MC_matrix.csv",output_dir,i), row.names = FALSE)
    n_hospitalized_per_day = as.numeric(sapply(ith_run_state_per_day_matrix, function(col) sum(col %in% count_states)))

    stopifnot(length(n_hospitalized_per_day) == length(dates))    
    
    return(n_hospitalized_per_day)
  },1:n_monte_carlo_runs,mc.cores=n_cores)
  n_hospitalized_per_run_matrix = data.frame(t(n_hospitalized_per_run_matrix))
  colnames(n_hospitalized_per_run_matrix) = dates
  return(n_hospitalized_per_run_matrix)
}


plot_expected_number_of_patients_each_day = function(patients_path, # make sure to run: patients = run_monte_carlo_per_patient(patients)
                                                     first_date,
                                                     last_date,
                                                     states_to_include, # default is all "hospitalized" states
                                                     n_monte_carlo_runs,
                                                     true_counts, n_cores=1, prefix="", 
				  	output_dir="~/") {
  
  n_hospitalized_per_run_matrix = construct_n_hospitalized_per_run_matrix(patients_path, output_dir,
                                                                          first_date = first_date,
                                                                          last_date = last_date, 
                                                                          count_states=states_to_include,
                                                                          n_monte_carlo_runs=n_monte_carlo_runs, n_cores=n_cores)
  
  
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
                  
                  true_counts = true_counts,
                  
                  date=dates,
                  metric=rep('n patients', length(dates)))
  
  top_limit = max(c(estimated_q90_n_hospitalized, true_counts)) + 20
  bottom_limit = min(c(estimated_q10_n_hospitalized, true_counts)) - 20
  write.csv(df, file = paste0(output_dir, prefix, paste(states_to_include, collapse = "-"), ".csv"), row.names = FALSE)
}


# NOTE!
# Function assumes previously hospitalized patients.
#
# Below you can see the expected proportion of deaths by days since hospitalization.
# That is, we are counting deaths among all patients 5 days since hospitalization, 10 days hospitalization etc.; averaging over all sampled monte carlo paths.
#

expected_deaths_by_days_since_hospitalization = function(run_cache_dir, T, PREDICT_N_DAYS_AHEAD, output_dir) {


  expected_deaths_df = data.frame(days_since_hospitalization = seq(5, PREDICT_N_DAYS_AHEAD, by=5))
  expected_deaths = c()

  for (n_days in expected_deaths_df$days_since_hospitalization) {

    number_of_deaths_before_n_days_since_hospitalization_per_run = sapply(1:N_MONTE_CARLO_RUNS, function(run_i) {

      patients = readRDS(sprintf("%s/%d.R",run_cache_dir, run_i))
      stopifnot(all(sapply(patients, function(patient) !is.na(patient$date_of_hospitalization) ))) # date is available for hospitalized patients only (with or without history)

      # count the number of patients dead by n_days in run_i
      number_of_deaths_before_n_days_since_hospitalization = sum(sapply(patients, function(patient) {

        patient_run_i = patient$all_runs[[1]]
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

  patients = readRDS(sprintf("%s/1.R",run_cache_dir))

  expected_deaths_df$n_expected_deaths = as.vector(expected_deaths)
  expected_deaths_df$expected_proportion_of_deaths = as.vector(expected_deaths) / length(patients)


  write.csv(expected_deaths_df, file = paste0(output_dir, "expected_deaths.csv"), row.names = FALSE)
  # return(expected_deaths_df)
}
