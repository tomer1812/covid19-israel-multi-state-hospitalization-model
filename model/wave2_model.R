source('./model/multi_state_competing_risks_model.R')
source('./utils/wave2_utils.R')


Wave2MultiStateModel = setRefClass('Wave2MultiStateModel', contains="MultiStateModel")

# specific patches needed for israeli data:
Wave2MultiStateModel$methods(
  one_monte_carlo_run = function(sample_covariates,
                                 origin_state,
                                 max_transitions,
                                 current_time = 0) {
    
    run = list(states = c(), 
               time_at_each_state = c())
    
    current_state = origin_state
    for (i in 1:max_transitions) {
      
      if (current_state == DISCHARGED && 
          stays_discharged_forever(sample_covariates, current_time)) {
        
          run$states = c(run$states, current_state)
          return(run)
          
      } else {
        next_state = sample_next_state(current_state, 
                                       sample_covariates, 
                                       current_time)  
      }
      
      # stop if no transition from current state was seen for this entry time:
      if (is.null(next_state)) {
        run$stopped_early = TRUE
        return(run)
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
        # Note: From POSITIVE state, next_state might be "mild" or "moderate" (which aren't "states")
        sample_covariates = update_covariates(sample_covariates,
                                              origin_state = current_state, 
                                              target_state = next_state, 
                                              time_at_origin_state = time_to_next_state,
                                              absolute_time_of_entry_to_target_state = current_time + time_to_next_state)
        
        if (current_state == POSITIVE) {
          current_state = ifelse(next_state %in% c(MILD, MODERATE), MILD_OR_MODERATE, next_state)
          current_time = 0 # reset T=0 for day of hospitalization
        } else {
          current_state = next_state
          current_time = current_time + time_to_next_state          
        }
      }
    }
    
    return(run)
  }, 
  
  
  probability_stays_discharged_forever = function(sample_covariates, t_entry_to_state) {
    
    competing_risks_model = state_specific_models[[DISCHARGED]]
    possible_next_states = competing_risks_model$failure_types
    
    probabilities = sapply(possible_next_states, function(state){
                            probability_for_next_state(state, 
                                                       competing_risks_model, 
                                                       sample_covariates, 
                                                       t_entry_to_state)
                          })
    
    return(1-sum(probabilities))
  },
  
  
  stays_discharged_forever = function(sample_covariates, t_entry_to_state) {
    return(rbinom(1,1,probability_stays_discharged_forever(sample_covariates, t_entry_to_state)) == 1)
  },
  
  # construct a single dataframe where each row is composed of the following, concatenated:
  # origin_state, target_state / censoring, time_at_origin_state, covariates_at_entry_to_origin_state
  prepare_dataset_for_competing_risks_fit = function(dataset, terminal_states, update_covariates_function) {
    
    covariate_column_names = get_covariate_names(dataset)

    nrows = sum(sapply(dataset, function(obj) length(obj$time_at_each_state))) # number of transitions
    ncols = 4 + length(covariate_column_names) + (!is.null(dataset[[1]]$id)) + (!is.null(dataset[[1]]$weight)) # see construction of row below:
    
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
        next_row[, covariate_column_names] = c(covariates_entering_origin_state)
        if (!is.null(obj$id)) next_row[, ID] = c(obj$id)
        if (!is.null(obj$weight)) next_row[, WEIGHT] = c(obj$weight)
        
        # PATCH for our data exclude rows with bad transitions
        #--------------------------------------------------------------------------------------------
        if (!any(sapply(EXCLUDED_TRANSITIONS, function(transition) transition[1] == origin_state & 
                                                                   transition[2] == target_state))) {
          
          competing_risks_dataset[row_idx,] = as.numeric(next_row)
          row_idx = row_idx + 1
          
        } else {
          break
        }
        #--------------------------------------------------------------------------------------------
        
        
        # PATCH for our data, non-terminal states with no follow-up time. 
        # there should be no non-terminal states without time at that state
        #--------------------------------------------------------------------------------------------
        
        if (length(obj$time_at_each_state) == i) break
        
        # -------------------------------------------------------------------------------------------
        
        # if no more transitions:
        if ( target_state == RIGHT_CENSORING | target_state %in% terminal_states) {
          break
          
          # else, setup next iteration:
        } else {
          covariates_entering_origin_state = update_covariates_function(covariates_entering_origin_state,
                                                                        origin_state = origin_state,
                                                                        target_state = target_state,
                                                                        time_at_origin_state = time_at_origin_state,
                                                                        absolute_time_of_entry_to_target_state = t_transition_to_target_state)
          origin_state = target_state
          t_entry_to_origin_state = t_transition_to_target_state
        }
      }
    }
    
    competing_risks_dataset = data.frame(competing_risks_dataset)
    colnames(competing_risks_dataset) = colnames(next_row)
    
    competing_risks_dataset = competing_risks_dataset[!is.na(competing_risks_dataset$origin_state), ]
    
    return(competing_risks_dataset)
  }
)
