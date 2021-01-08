source('./model/competing_risks_model.R')

RIGHT_CENSORING = 0
ID = 'id'
WEIGHT = 'weight'

MultiStateModel = setRefClass('MultiStateModel',
                              fields = list(
                                state_specific_models = "list",
                                update_covariates = "function", 
                                covariate_names = "character",
                                time_is_discrete = "logical",
                                terminal_states = "numeric"
                              )
)

MultiStateModel$methods(
  
  # fit()
  # 
  # Description:
  # ------------------------------------------------------------------------------------------------------------
  # This method fits a competing risks model per state, that is, it treats all state transitions as competing risks. 
  # See the CompetingRisksModel class.
  #
  # Arguments:
  # ------------------------------------------------------------------------------------------------------------
  # dataset: list of objects 
  #   Each object should have the following attributes:
  #     a. obj$covariates: the sample covariates at the initial state.
  #     b. obj$states : a vector of the states visitied (encoded as positive integers, 0 is saved for censoring), 
  #       in the order visited. 
  #     c. obj$time_at_each_state : a vector with the duration spent at each state.
  # 
  # optionally - obj$id : an identification of this sample
  # 
  # Note: if the last state is a terminal state, then the vector of times should be shorter than the vector 
  # of states by 1. 
  # Conversely, if the last state is not a terminal state, then the length of vector times should be the same 
  # as that of the states.
  #
  # update_covariates_function: function
  #   A state-transition function, which updates the time dependent variables. This function is used in fitting 
  #   the model so that the user doesn't have to manually compute the feautres at each state, it is also used
  #   in the monte carlo method of predicting the survival curves per sample.
  #
  #    Arguments:
  #       covariates_upon_entering_origin_state
  #
  #       optional named arguments:
  #         a. origin_state
  #         b. target_state
  #         c. time_at_origin_state
  #     
  #    Returns:
  #     covariates_upon_entering_target_state
  #     
  #    Note: if the function is not specified, the default shall be the identity function, 
  #          corresponding with no time dependent variables.
  #
  # terminal_states: numeric vector
  #   the states which a sample does not leave
  #
  # optional: covariate_names 
  #   list of covariate names to be used in prints
  #
  fit = function(dataset, 
                 terminal_states, 
                 update_covariates_function = default_update_covariates_function,
                 covariate_names=NULL,
                 verbose=1) {
    
    assert_valid_fit_input(dataset, 
                           terminal_states,
                           update_covariates_function,
                           covariate_names)
    
    terminal_states <<- terminal_states
    covariate_names <<- if(!is.null(covariate_names)) covariate_names else get_covariate_names(dataset)
    update_covariates <<- update_covariates_function
    competing_risks_dataset = prepare_dataset_for_competing_risks_fit(dataset, 
                                                                      terminal_states,
                                                                      update_covariates)
    
    time_is_discrete <<- check_if_time_is_discrete(competing_risks_dataset)
      
    for (state in unique(competing_risks_dataset$origin_state)) {
      if (verbose >= 1) print(paste("Fitting Model at State: ", state))
      
      model = fit_state_specific_model(state, competing_risks_dataset, verbose=verbose)
      state_specific_models[[state]] <<- model
    }
  }
)


MultiStateModel$methods(
  
  # a simple identity function:
  # Typical update functions should include the following named arguments:
  # 1. origin_state
  # 2. target_state 
  # 3. time_at_origin_state
  #
  default_update_covariates_function = function(covariates_entering_origin_state, ...) {
    return(covariates_entering_origin_state)
  },
  
  assert_valid_fit_input = function(dataset, terminal_states, update_covariates_function, covariate_names) {
    assert_valid_dataset(dataset, terminal_states)
    assert_valid_update_covariates_function(update_covariates_function)
    assert_valid_covariate_names(dataset, covariate_names)
  },
  
  assert_valid_dataset = function(dataset, terminal_states) {
    
    # number of states should equal number of times, or exceed it by one
    for (obj in dataset) {
      n_states = length(obj$states)
      n_times = length(obj$time_at_each_state)
      stopifnot(n_states == n_times | n_states == n_times + 1)
      
      if ( n_states == 1 & obj$states[[1]] %in% terminal_states) {
        print("Error: encountered a sample with a single state that is a terminal state.\n The object:")
        print(obj)
        stop()
      }
    }
    
    # either all objects have an id, or none have
    has_id = sapply(dataset, function(x) !is.null(x$id))
    stopifnot((!any(has_id)) | all(has_id))
    
    # all covariates are of the same length:
    l = length(dataset[[1]]$covariates)
    have_same_length_as_first_sample = sapply(dataset, function(obj) length(obj$covariates) == l)
    stopifnot(all(have_same_length_as_first_sample))
  },
  
  assert_valid_update_covariates_function = function(update_covariates_function) {
    return()
  },
  
  assert_valid_covariate_names = function(dataset, covariate_names) {
    if (is.null(covariate_names)) return()
    stopifnot(length(covariate_names) == length(dataset[[1]]$covariates))
  },
  
  # construct a single dataframe where each row is composed of the following, concatenated:
  # origin_state, target_state / censoring, time_at_origin_state, covariates_at_entry_to_origin_state
  slow_prepare_dataset_for_competing_risks_fit = function(dataset, terminal_states, update_covariates_function) {

    competing_risks_dataset = data.frame()
    
    for (obj in dataset) {
      
      n_states = length(obj$states)
      origin_state = obj$states[[1]]
      covariates_entering_origin_state = obj$covariates
      t_entry_to_origin_state = 0
      
      for (i in 1:n_states) {

        time_at_origin_state = obj$time_at_each_state[[i]]
        t_transition_to_target_state = t_entry_to_origin_state + time_at_origin_state
        
        if (i+1 <= n_states) {
          target_state = obj$states[[i+1]]
        } else {
          target_state = RIGHT_CENSORING
        }
        
        # append row corresponding to this transition:
        next_row = data.frame(origin_state = c(origin_state), 
                              target_state = c(target_state), 
                              t_entry_to_origin_state = c(t_entry_to_origin_state),
                              t_transition_to_target_state = c(t_transition_to_target_state))
        next_row[, covariate_names] = c(covariates_entering_origin_state)
        if (!is.null(obj$id)) next_row[, ID] = c(obj$id)
        if (!is.null(obj$weight)) next_row[, WEIGHT] = c(obj$weight)
        
        competing_risks_dataset = rbind(competing_risks_dataset, next_row)
        
        # if no more transitions:
        if ( target_state == RIGHT_CENSORING | target_state %in% terminal_states ) {
          break
          
        # else, setup next iteration:
        } else {
          covariates_entering_origin_state = update_covariates_function(covariates_entering_origin_state,
                                                                        origin_state = origin_state,
                                                                        target_state = target_state,
                                                                        time_at_origin_state = time_at_origin_state)
          origin_state = target_state
          t_entry_to_origin_state = t_transition_to_target_state
        }
      }
    }
    
    return(competing_risks_dataset)
  },
   
  # construct a single dataframe where each row is composed of the following, concatenated:
  # origin_state, target_state / censoring, time_at_origin_state, covariates_at_entry_to_origin_state
  prepare_dataset_for_competing_risks_fit = function(dataset, terminal_states, update_covariates_function) {
    
    nrows = sum(sapply(dataset, function(obj) length(obj$time_at_each_state))) # number of transitions
    ncols = 4 + length(covariate_names) + (!is.null(dataset[[1]]$id)) + (!is.null(dataset[[1]]$weight)) # see construction of row below:

    competing_risks_dataset = matrix(nrow=nrows, ncol=ncols)
    row_idx = 1
    
    for (obj in dataset) {
      
      n_states = length(obj$states)
      origin_state = obj$states[[1]]
      covariates_entering_origin_state = obj$covariates
      t_entry_to_origin_state = 0
      
      for (i in 1:n_states) {

        time_at_origin_state = obj$time_at_each_state[[i]]
        t_transition_to_target_state = t_entry_to_origin_state + time_at_origin_state
        
        if (i+1 <= n_states) {
          target_state = obj$states[[i+1]]
        } else {
          target_state = RIGHT_CENSORING
        }
        
        # append row corresponding to this transition:
        next_row = data.frame(origin_state = c(origin_state), 
                              target_state = c(target_state), 
                              t_entry_to_origin_state = c(t_entry_to_origin_state),
                              t_transition_to_target_state = c(t_transition_to_target_state))
        next_row[, covariate_names] = c(covariates_entering_origin_state)
        if (!is.null(obj$id)) next_row[, ID] = c(obj$id)
        if (!is.null(obj$weight)) next_row[, WEIGHT] = c(obj$weight)
        
        competing_risks_dataset[row_idx,] = as.numeric(next_row)
        row_idx = row_idx + 1
        
        
        # if no more transitions:
        if ( target_state == RIGHT_CENSORING | target_state %in% terminal_states ) {
          break
          
        # else, setup next iteration:
        } else {
          covariates_entering_origin_state = update_covariates_function(covariates_entering_origin_state,
                                                                        origin_state = origin_state,
                                                                        target_state = target_state,
                                                                        time_at_origin_state = time_at_origin_state)
          origin_state = target_state
          t_entry_to_origin_state = t_transition_to_target_state
        }
      }
    }
    
    competing_risks_dataset = data.frame(competing_risks_dataset)
    colnames(competing_risks_dataset) = colnames(next_row)
    
    return(competing_risks_dataset)
  },
  
  check_if_time_is_discrete = function(competing_risks_dataset) {
    times = c(competing_risks_dataset$t_entry_to_origin_state, competing_risks_dataset$t_transition_to_target_state)
    if(all(round(times) == times)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  },
  
  get_covariate_names = function(dataset) {
    if (identical(covariate_names, character(0))) {
      return(c(paste('x', 1:length(dataset[[1]]$covariates), sep='')))  
    } else {
      return(covariate_names)
    }
  },
  
  fit_state_specific_model = function(state, competing_risks_dataset, covariate_names_ = covariate_names, verbose=1) {
    
    # prepare arguments for model:
    mask = (competing_risks_dataset$origin_state == state)
    
    t_entry =       competing_risks_dataset$t_entry_to_origin_state[mask]
    t_transition =  competing_risks_dataset$t_transition_to_target_state[mask]
    target_state =  competing_risks_dataset$target_state[mask]
    covariates_X =  competing_risks_dataset[mask, covariate_names_]
    id =            if (ID %in% names(competing_risks_dataset)) competing_risks_dataset$id[mask] else NULL
    weights =       if (WEIGHT %in% names(competing_risks_dataset)) competing_risks_dataset$weight[mask] else NULL
    
    # fit model:
    model = CompetingRisksModel()
    model$fit(t_transition, target_state, covariates_X, sample_ids = id, t_start=t_entry, sample_weights = weights, verbose=verbose)
    
    return(model)
  },
  
  run_monte_carlo_simulation = function(sample_covariates,
                                        origin_state,
                                        current_time = 0,
                                        n_random_samples = 100,
                                        max_transitions = 10, n_cores = 1) {
    
    lapply(1:n_random_samples, function(run_i) one_monte_carlo_run(sample_covariates,
                                                                   origin_state,
                                                                   max_transitions,
                                                                   current_time = current_time))
  },
  
  one_monte_carlo_run = function(sample_covariates,
                                 origin_state,
                                 max_transitions,
                                 current_time = 0) {
    
    run = list(states = c(), 
               time_at_each_state = c())
    
    current_state = origin_state
    for (i in 1:max_transitions) {

      next_state = sample_next_state(current_state,
                                       sample_covariates,
                                       current_time)

      # stop if no transition from current state was seen for this entry time:
      if (is.null(next_state)) {
        run$stopped_early = TRUE
        return(run)
      }
      
      if (next_state == HOME) {
        return(list(states = c(POSITIVE, HOME), 
               time_at_each_state = c(1)))
      }
      
      time_to_next_state = sample_time_to_next_state(current_state,
                                                     next_state,
                                                     sample_covariates,
                                                     current_time)
      run$states = c(run$states, current_state)
      run$time_at_each_state = c(run$time_at_each_state, time_to_next_state)
      
      if (next_state %in% terminal_states) {
        run$states = c(run$states, next_state)
        break
      } else {
        sample_covariates = update_covariates(sample_covariates,
                                              origin_state = current_state, 
                                              target_state = next_state, 
                                              time_at_origin_state = time_to_next_state,
                                              absolute_time_of_entry_to_target_state = current_time + time_to_next_state)
        current_state = next_state
        current_time = current_time + time_to_next_state
      }
    }
    
    return(run)
  },
  
  
  
  probability_for_next_state = function(next_state, competing_risks_model, sample_covariates, t_entry_to_current_state=0) {
    
    unique_event_times = competing_risks_model$unique_event_times(next_state)
    
    # ensure discrete variables are sampled from the next time unit
    if (time_is_discrete) {
      mask = (unique_event_times > floor(t_entry_to_current_state + 1))  
    } else {
      mask = (unique_event_times > t_entry_to_current_state)  
    }
    
    # hazard for the failure type corresponding to 'state': 
    hazard = competing_risks_model$hazard_at_unique_event_times(sample_covariates, next_state)[mask]
    
    # overall survival function evaluated at time of failures corresponding to 'state'
    survival = competing_risks_model$survival_function(unique_event_times[mask], sample_covariates)
    
    probability_for_state = sum(hazard*survival)
    
    return(probability_for_state)
  },
  
  
  sample_next_state = function(current_state, sample_covariates, t_entry_to_current_state) {
    
    competing_risks_model = state_specific_models[[current_state]] 
    possible_next_states = competing_risks_model$failure_types
    
    # compute probabilities for multinomial distribution:
    probabilities = lapply(possible_next_states, function(state) probability_for_next_state(state, 
                                                                                            competing_risks_model, 
                                                                                            sample_covariates,
                                                                                            t_entry_to_current_state))
    # when no transition after t_entry_to_current_state was seen:
    if (all(probabilities == 0)) return(NULL)
    
    mult = rmultinom(1, 1, prob = probabilities)
    next_state = possible_next_states[which(mult==max(mult))]
    return(next_state)
  },

  sample_time_to_next_state = function(current_state, next_state, sample_covariates, t_entry_to_current_state) {
    
    competing_risks_model = state_specific_models[[current_state]]
    
    unique_event_times = competing_risks_model$unique_event_times(next_state)
    
    # ensure discrete variables are sampled from the next time unit
    if (time_is_discrete) {
      mask = (unique_event_times > floor(t_entry_to_current_state + 1))
    } else {
      mask = (unique_event_times > t_entry_to_current_state)  
    }
    
    unique_event_times = unique_event_times[mask]
    
    # hazard for the failure type corresponding to 'state': 
    hazard = competing_risks_model$hazard_at_unique_event_times(sample_covariates, next_state)[mask]
    
    # overall survival function evaluated at time of failures corresponding to 'state'
    survival = competing_risks_model$survival_function(unique_event_times, sample_covariates)
    
    probability_for_each_t = cumsum(hazard*survival)
    probability_for_each_t_given_next_state = probability_for_each_t / max(probability_for_each_t)
    
    eps = runif(1)
    # take the first event time whose probability is less than or equal to eps
    # if we drew a very small eps, use the minimum observed time
    time_to_next_state = max(c(unique_event_times[probability_for_each_t_given_next_state <= eps], unique_event_times[1]))
    time_to_next_state = time_to_next_state - t_entry_to_current_state
    
    return(time_to_next_state)
  }
)
