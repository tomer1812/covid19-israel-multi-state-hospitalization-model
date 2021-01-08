library(plyr)
library(survival)
library(methods)


CompetingRisksModel = setRefClass('CompetingRisksModel',
                                  fields = list(
                                      failure_types = "numeric",
                                      
                                      # Each element of "event_specific_models"is a list 
                                      # with the following attributes:
                                      # 1. coefficients
                                      # 2. unique_event_times
                                      # 3. baseline_hazard
                                      # 4. cumulative_baseline_hazard_function
                                      event_specific_models = "list"
                                    )
                                  )
# MAIN API:
CompetingRisksModel$methods(
  # fit()
  # 
  # Description:
  # ------------------------------------------------------------------------------------------------------------
  # This method fits a cox proportional hazards model for each failure type, treating others as censoring events. 
  # Tied event times are dealt with by adding an epsilon to tied event times.
  #
  # Arguments:
  # ------------------------------------------------------------------------------------------------------------
  # t : numeric vector 
  #   A length n vector of positive times of events
  #
  # failure_types: numeric vector 
  #   The event type corresponding to the time in vector t.
  #   Failure types are encoded as integers from 1 to m. 
  #   Right-censoring events (the only kind of censoring supported) are encoded as 0. 
  #   Thus, the failure type argument holds integers from 0 to m, where m is the number of distinct failure types
  #
  # covariates_X: numeric dataframe
  #   an n by #(covariates) numerical matrix
  #   All columns are used in the estimate.
  #
  # OPTIONAL:
  #
  # sample_ids:
  #   used inside the coxph model in order to identify subjects with repeating entries.
  #
  # t_start: 
  #   A length n vector of positive start times, used in case of interval data.
  #   In that case: left=t_start, and right=t
  #
  # epsilon_min/max: 
  #   epsilon is added to events with identical times to break ties. 
  #   epsilon is sampled from a uniform distribution in the range (epsilon_min, epsilon_max)
  #   these values should be chosen so that they do not change the order of the events. 
  fit = function(t, 
                 failure_types, 
                 covariates_X, 
                 sample_ids=NULL,
                 t_start=NULL, 
                 break_ties=TRUE, 
                 sample_weights=NULL,
                 epsilon_min=0.0, 
                 epsilon_max=0.0001,
                 verbose=1) {

    assert_valid_dataset(t, failure_types, covariates_X)
    
    if (break_ties) t = break_ties_by_adding_epsilon(t, epsilon_min, epsilon_max)
    
    failure_types <<- unique(failure_types[failure_types > 0])
    for (type in .self$failure_types) {
      cox_model = fit_event_specific_model(t, failure_types, covariates_X, type, sample_ids, t_start, sample_weights, verbose)
      event_specific_models[[type]] <<- extract_necessary_attributes(cox_model)
    }
  },
  
  
  # predict_CIF()
  #
  # Description:
  # ------------------------------------------------------------------------------------------------------------
  # This method computes the failure-type-specific cumulative incidence function, given that 'time_passed' time
  # has passed (default is 0)
  # 
  # Arguments:
  # ------------------------------------------------------------------------------------------------------------
  # predict_at_t: numeric vector
  #   times at which the cif will be computed 
  #
  # sample_covariates: numeric vector
  #   a numerical vector of same length as the covariate matrix the model was fit to.
  #
  # failure_type: integer
  #   integer corresponing to the failure type, as given when fitting the model   
  # 
  # time_passed: numeric
  #   compute the cif conditioned on the fact that this amount of time has already passed. 
  # 
  # Returns:
  # ------------------------------------------------------------------------------------------------------------
  # the predicted cumulative incidence values for the given sample_covariates at times predict_at_t.   
  # 
  #
  predict_CIF = function(predict_at_t, sample_covariates, failure_type, time_passed = 0) {
    cif_function = compute_cif_function(sample_covariates, failure_type)
    
    predictions = cif_function(predict_at_t)
    
    # re-normalize the probability to account for the time passed
    if (time_passed  > 0) {
      predictions = (predictions - cif_function(time_passed)) / survival_function(time_passed, sample_covariates)
    }
    
    return(predictions)
  }
  
)



# These are inner, helper functions. "outside of the API"
CompetingRisksModel$methods(

  assert_valid_dataset = function(t, failure_types, covariates_X) { 
    
      # t should be non-negative
      stopifnot(t >= 0)

      # failure types should be integers from 0 to m, not necessarily consecutive
      stopifnot(failure_types %% 1 == 0) # integers
      stopifnot(min(failure_types) >= 0)

      # covariates should all be numerical
      stopifnot(sapply(covariates_X, is.numeric))
      
      # all 3 arguments should have the same length of n
      stopifnot(length(t) == length(failure_types))
      stopifnot(nrow(covariates_X) == length(t))
      
  },
  
  
  break_ties_by_adding_epsilon = function(t, epsilon_min = 0.0, epsilon_max = 0.0001) {
    set.seed(42)
    counts = count(t)
    non_unique_times = counts$x[counts$freq > 1]
    eps = runif(length(t), epsilon_min, epsilon_max)
    
    # note: leave zero times as-is, to avoid affecting model fit
    t + ((t != 0) *(t %in% non_unique_times) * eps)
  }
)


# These are the model-specific methods needed to be overriden when implementing RSF
CompetingRisksModel$methods(
  
  # Treat all 'failure_types' except 'type' as censoring events
  fit_event_specific_model = function(t, failure_types, covariates_X, type, sample_ids=NULL, t_start=NULL, sample_weights=NULL, verbose=TRUE) {
    is_event = (failure_types == type)
    
    if (verbose >= 1) print(paste(">>> Fitting Transition to State: ", type, ", n events: ", sum(is_event)))
    
    surv_object = if (is.null(t_start)) Surv(t, is_event) else Surv(t_start, t, is_event)
    cox_model = coxph(surv_object ~ . + cluster(sample_ids), 
                      weights = sample_weights, 
                      data = covariates_X)
    
    if (verbose >= 2) print(cox_model)
    
    return(cox_model)
  },
  
  
  # constructs a cif step function
  compute_cif_function = function(sample_covariates, failure_type) {
    cif_x = unique_event_times(failure_type)
    cif_y = cumsum(hazard_at_unique_event_times(sample_covariates, failure_type)*survival_function(cif_x, sample_covariates))
    return(stepfun(cif_x, c(0, cif_y)))
  },
  
  # the hazard is given by multiplying the baseline hazard (which has value per unique event time) by the partial hazard 
  hazard_at_unique_event_times = function(sample_covariates, failure_type) {
    hazard = baseline_hazard(failure_type) * c(partial_hazard(failure_type, sample_covariates))
    
    stopifnot(length(hazard) == length(unique_event_times(failure_type)))
    return(hazard)
  },
  

  # the cumulative baseline hazard is given as a non-paramateric function, whose values are given at the times of observed events
  # the cumulative baseline hazard is the sum of hazards at observed event times
  cumulative_baseline_hazard = function(cox_model) {
    
    cumulative_baseline_hazard = basehaz(cox_model, centered = FALSE)
    
    # step > 0 corresponds exactly to unique time events
    mask = (diff(c(0, cumulative_baseline_hazard$hazard)) > 0)
    
    return(cumulative_baseline_hazard$hazard[mask])
  },
  
  cumulative_baseline_hazard_step_function = function(cox_model) {
    return(stepfun(coxph.detail(cox_model)$time, 
                   c(0,cumulative_baseline_hazard(cox_model))))
  },
  # the baseline hazard is given as a non-paramateric function, whose values are given at the times of observed events
  # the cumulative hazard is the sum of hazards at times of events, the hazards are then the diffs 
  baseline_hazard = function(failure_type) {
    return(event_specific_models[[failure_type]]$baseline_hazard)
  },
  
  # cache all relevant model attributes in one method:
  # ALWAYS CACHE AN UNALTERED COX MODEL! THE COMPUTATION OF THE BASELINE HAZARD IS BASED ON THE COVARIATES!
  extract_necessary_attributes = function(cox_model) {
      
    if (0 %in% cox_model$coefficients) print(" ------ WARNING ------: Are you caching a modified cox model?")
    
     model = list(
        coefficients = cox_model$coefficients,
        unique_event_times = coxph.detail(cox_model)$time,
        baseline_hazard = diff(c(0, cumulative_baseline_hazard(cox_model))),
        cumulative_baseline_hazard_function = cumulative_baseline_hazard_step_function(cox_model)
      )
     
     return(model)
  },
  
  # simply e^x_dot_beta for the chosen failure type's coefficients  
  partial_hazard = function(failure_type, sample_covariates) {
    coefs = event_specific_models[[failure_type]]$coefficients
    x_dot_beta = as.numeric(sample_covariates) %*% coefs
    return(exp(x_dot_beta)) 
  },
  
  # uses a coxph function which returns unique times, regardless of the original fit which might have tied times. 
  unique_event_times = function(failure_type) {
    return(event_specific_models[[failure_type]]$unique_event_times) 
  },
  
  # simply: exp( sum of cumulative hazards of all types )
  survival_function = function(t, sample_covariates) {
    exponent = rep(0, length(t))
    for (type in failure_types) {
      exponent = exponent - ( event_specific_models[[type]]$cumulative_baseline_hazard_function(t) * c(partial_hazard(type, sample_covariates)) )
    }
    
    survival_function_at_t = exp(exponent)
    
    stopifnot(length(survival_function_at_t) == length(t))
    return(survival_function_at_t)
  }


)
