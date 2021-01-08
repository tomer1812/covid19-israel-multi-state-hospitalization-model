source('./model/multi_state_competing_risks_model.R')


RECOVERED_OR_OOHQ = 16
MILD_OR_MODERATE = 23
SEVERE = 4
DECEASED = 5

BAD_TRANSITIONS = list(c(RECOVERED_OR_OOHQ, DECEASED), 
                       c(RECOVERED_OR_OOHQ, SEVERE),
                       c(SEVERE, RECOVERED_OR_OOHQ))


IsraeliMultiStateModel = setRefClass('IsraeliMultiStateModel', contains="MultiStateModel")


# function for overwriting the mild or moderate model, in order to deal with non-convergence.
get_modified_mild_or_moderate_to_deceased_cox_model = function(model, dataset, terminal_states, update_covariates){

  competing_risks_dataset = model$prepare_dataset_for_competing_risks_fit(dataset = dataset,
                                                                          terminal_states = terminal_states,
                                                                          update_covariates_function = update_covariates)
  
  type0_3_or_4 = competing_risks_dataset[,"as.factor(type0)3"] | competing_risks_dataset[,"as.factor(type0)4"]
  modified_covariates = data.frame(
    'age' = competing_risks_dataset$age,
    'sexMale' = competing_risks_dataset$sexMale,
    'type0-3or4' = type0_3_or_4,
    'cumulative_time' = competing_risks_dataset$cumulative_time,
    'age:sexMale' = competing_risks_dataset$`age:sexMale`,
    'age:type0-3or4' = competing_risks_dataset$age * type0_3_or_4,
    'age:cumulative-time' = competing_risks_dataset$age * competing_risks_dataset$cumulative_time
  )
  
  
  modified_covariate_names = colnames(modified_covariates)
  
  non_covariate_column_names = c('origin_state', 'target_state', 't_entry_to_origin_state', 't_transition_to_target_state', ID, WEIGHT)
  non_covariate_column_names = non_covariate_column_names[non_covariate_column_names %in% colnames(competing_risks_dataset)] # might not have weight
  modified_competing_risks_dataset = cbind(competing_risks_dataset[, non_covariate_column_names], modified_covariates)
  
  modified_mild_or_moderate_cr_model = model$fit_state_specific_model(state = MILD_OR_MODERATE,
                                                                      competing_risks_dataset = modified_competing_risks_dataset,
                                                                      covariate_names_ = modified_covariate_names)
  
  print("The fitted Cox Model from State Mild or Moderate to DECEASED")
  print(modified_mild_or_moderate_cr_model$event_specific_models[[DECEASED]])
  
  coefs = modified_mild_or_moderate_cr_model$event_specific_models[[DECEASED]]$coefficients
  # The reorganized regression coefficient vector= (b1 , b2 , b3 , b3 , b4 , 0 , b5 , b6 , b6 , b7 , 0):
  new_coefs = c(coefs[c(1,2,3,3,4)], 0, coefs[c(5,6,6,7)],0)
  modified_mild_or_moderate_cr_model$event_specific_models[[DECEASED]]$coefficients = new_coefs
  
  return(modified_mild_or_moderate_cr_model$event_specific_models[[DECEASED]])
}


# function for overwriting the recovered_or_oohq to mild_or_moderate model, in order to deal with non-convergence.
get_modified_recovered_or_oohq_to_mild_or_moderate_cox_model = function(model, dataset, terminal_states, update_covariates){
  
  competing_risks_dataset = model$prepare_dataset_for_competing_risks_fit(dataset = dataset,
                                                                          terminal_states = terminal_states,
                                                                          update_covariates_function = update_covariates)
  
  type0_3_or_4 = competing_risks_dataset[,"as.factor(type0)3"] | competing_risks_dataset[,"as.factor(type0)4"]
  modified_covariates = data.frame(
    'age' = competing_risks_dataset$age,
    'sexMale' = competing_risks_dataset$sexMale,
    'type0-3or4' = type0_3_or_4,
    'cumulative_time' = competing_risks_dataset$cumulative_time,
    'age:sexMale' = competing_risks_dataset$`age:sexMale`,
    'age:type0-3or4' = competing_risks_dataset$age * type0_3_or_4,
    'age:cumulative-time' = competing_risks_dataset$age * competing_risks_dataset$cumulative_time
  )
  
  
  modified_covariate_names = colnames(modified_covariates)
  
  non_covariate_column_names = c('origin_state', 'target_state', 't_entry_to_origin_state', 't_transition_to_target_state', ID, WEIGHT)
  non_covariate_column_names = non_covariate_column_names[non_covariate_column_names %in% colnames(competing_risks_dataset)] # might not have weight
  modified_competing_risks_dataset = cbind(competing_risks_dataset[, non_covariate_column_names], modified_covariates)
  
  modified_recovered_or_oohq_cr_model = model$fit_state_specific_model(state = RECOVERED_OR_OOHQ,
                                                                       competing_risks_dataset = modified_competing_risks_dataset,
                                                                       covariate_names_ = modified_covariate_names)
  
  print("The fitted Cox Model from State RECOVERED_OR_OOHQ to MILD_OR_MODERATE")
  print(modified_recovered_or_oohq_cr_model$event_specific_models[[MILD_OR_MODERATE]])
  
  coefs = modified_recovered_or_oohq_cr_model$event_specific_models[[MILD_OR_MODERATE]]$coefficients
  # The reorganized regression coefficient vector= (b1 , b2 , b3 , b3 , b4 , 0 , b5 , b6 , b6 , b7 , 0):
  new_coefs = c(coefs[c(1,2,3,3,4)], 0, coefs[c(5,6,6,7)],0)
  modified_recovered_or_oohq_cr_model$event_specific_models[[MILD_OR_MODERATE]]$coefficients = new_coefs
  
  return(modified_recovered_or_oohq_cr_model$event_specific_models[[MILD_OR_MODERATE]])
}



IsraeliMultiStateModel$methods(
  fit = function(dataset,
                 terminal_states,
                 update_covariates_function = default_update_covariates_function,
                 covariate_names=NULL) {
    
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
      print(paste("Fitting Model at State: ", state))
      model = fit_state_specific_model(state, competing_risks_dataset)
      state_specific_models[[state]] <<- model
    }
    
    # PATCH for our data, the following transitions had too little samples
    # For model convergence, specific models were fit to these transitions
    #--------------------------------------------------------------------------------------------
    
    print('---------------------------------------------------------')
    .self$state_specific_models[[MILD_OR_MODERATE]]$event_specific_models[[DECEASED]] =
      get_modified_mild_or_moderate_to_deceased_cox_model(.self,
                                                          dataset,
                                                          terminal_states,
                                                          update_covariates)
    
    print('---------------------------------------------------------')
    .self$state_specific_models[[RECOVERED_OR_OOHQ]]$event_specific_models[[MILD_OR_MODERATE]] =
      get_modified_recovered_or_oohq_to_mild_or_moderate_cox_model(.self,
                                                                   dataset,
                                                                   terminal_states,
                                                                   update_covariates)
    
    #--------------------------------------------------------------------------------------------
  }
)


# specific patches needed for israeli data:
IsraeliMultiStateModel$methods(
  one_monte_carlo_run = function(sample_covariates,
                                 origin_state,
                                 max_transitions,
                                 current_time = 0) {
    
    run = list(states = c(), 
               time_at_each_state = c())
    
    current_state = origin_state
    for (i in 1:max_transitions) {
      
      if (current_state == RECOVERED_OR_OOHQ) {
        if (stays_recovered_forever(sample_covariates, current_time)) {
          run$states = c(run$states, current_state)
          return(run)
        } else {
          next_state = MILD_OR_MODERATE
        }
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
  
  
  probability_stays_recovered_forever = function(sample_covariates, t_entry_to_state) {
    competing_risks_model = state_specific_models[[RECOVERED_OR_OOHQ]]
    p_stays_recovered_forever = 1-probability_for_next_state(MILD_OR_MODERATE, competing_risks_model, sample_covariates, t_entry_to_state)
    
    return(p_stays_recovered_forever)
  },
  
  
  stays_recovered_forever = function(sample_covariates, t_entry_to_state) {
    # print(paste("probability for MILD OR MODERATE:", 1-p_stays_recovered_forever))
    return(rbinom(1,1,probability_stays_recovered_forever(sample_covariates, t_entry_to_state)) == 1)
  },
  
  # construct a single dataframe where each row is composed of the following, concatenated:
  # origin_state, target_state / censoring, time_at_origin_state, covariates_at_entry_to_origin_state
  prepare_dataset_for_competing_risks_fit = function(dataset, terminal_states, update_covariates_function) {
    
    covariate_column_names = get_covariate_names(dataset)
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
        next_row[, covariate_column_names] = c(covariates_entering_origin_state)
        if (!is.null(obj$id)) next_row[, ID] = c(obj$id)
        if (!is.null(obj$weight)) next_row[, WEIGHT] = c(obj$weight)
        
        # PATCH for our data exclude rows with bad transitions
        #--------------------------------------------------------------------------------------------
        if (!any(sapply(BAD_TRANSITIONS, function(transition) transition[1] == origin_state & transition[2] == target_state))) {
          competing_risks_dataset = rbind(competing_risks_dataset, next_row)
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
    
    return(competing_risks_dataset)
  }
)




