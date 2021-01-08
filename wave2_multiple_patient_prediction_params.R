#path = "data/input_data/mean_input_2020-10-10"
args <- commandArgs(TRUE)
path = args[1]

if (!file.exists(paste0(path, "/model_file.Rds"))){
  rmarkdown::render("wave2-fit.Rmd", params = list(positive_to_hospitalized_csv=paste0(path, "/msm_fit/final_pos_for_model.csv"), 
                    hospitalized_onward_csv=paste0(path, "/msm_fit/final_for_model.csv"), 
                    model_cache_file=paste0(path, "/model_file.Rds")), 
                    output_file=paste0(path, "/model_output.html"))
}

args <- commandArgs(TRUE)
path = args[1]
analysis_date = strsplit(path, "input_")[[1]][2]  # expecting a defualt format of PATH/msm_input_ANALYSIS-DATE
cache_dir = paste0(path, "/run_cache")
system(paste("mkdir -p ", paste0(path, "/out" )))
system(paste("rm -rf", cache_dir))
system(paste("mkdir", cache_dir))

params <-
list(n_cores = 32L, n_monte_carlo_runs = 1000L, max_path_length = 14L,
    test_run = FALSE, model_cache_file = paste0(path, "/model_file.Rds"),
    analysis_date = analysis_date, predict_n_days_ahead = 8,
    seed = strtoi(42), 
    output_dir = paste0(path, "/out/"),
    positive_with_history_input_csv = paste0(path, "/input_group1.csv"), 
    positive_no_history_input_csv = paste0(path, "/input_group2.csv"),
    hospitalized_with_history_input_csv = paste0(path, "/input_group3.csv"), 
    hospitalized_no_history_input_csv = paste0(path, "/input_group4.csv"), 
    true_counts_csv = paste0(path, "/true_counts_all_groups.csv"),
    prefix = 'all_groups_')

## ------------------------------------------------------------------------------------------------------------
params

OUTPUT_DIR = params$output_dir
N_CORES = params$n_cores
N_MONTE_CARLO_RUNS = params$n_monte_carlo_runs
MAX_PATH_LENGTH = params$max_path_length
TEST_RUN = params$test_run
SEED = params$seed
PREFIX = params$prefix

MODEL_CACHE_FILE = params$model_cache_file

POSITIVE_WITH_HISTORY_INPUT_CSV = params$positive_with_history_input_csv
POSITIVE_NO_HISTORY_INPUT_CSV = params$positive_no_history_input_csv

HOSPITALIZED_WITH_HISTORY_INPUT_CSV = params$hospitalized_with_history_input_csv
HOSPITALIZED_NO_HISTORY_INPUT_CSV = params$hospitalized_no_history_input_csv

TRUE_COUNTS_CSV = params$true_counts_csv

T = as.Date(params$analysis_date)
PREDICT_N_DAYS_AHEAD = params$predict_n_days_ahead


## ----include=FALSE-------------------------------------------------------------------------------------------

source('./utils/wave2_utils.R')
source('./utils/wave2_multiple_patient_prediction_utils.R')
source('./model/wave2_model.R')


## ------------------------------------------------------------------------------------------------------------
model = readRDS(MODEL_CACHE_FILE)


## ------------------------------------------------------------------------------------------------------------
all_patients = c()


## ------------------------------------------------------------------------------------------------------------
if (!is.null(POSITIVE_WITH_HISTORY_INPUT_CSV)) {
  
  positive_with_history_df = read.csv(POSITIVE_WITH_HISTORY_INPUT_CSV,stringsAsFactors=FALSE)
  positive_with_history_df$patient_id <- NULL
  positive_with_history_df$pred_id <- NULL
  stopifnot(positive_with_history_df$first_positive_test_date <= T)
  stopifnot(sum(complete.cases(positive_with_history_df)) == nrow(positive_with_history_df))
  
  print(paste('Group 1: Positive With a History, ', nrow(positive_with_history_df), ' patients', sep = ''))

  positive_with_history_patient_objects = lapply(1:nrow(positive_with_history_df), function(row_idx) {
    
    row = positive_with_history_df[row_idx,]
    
    list(
      t_start = T,
      covariates_at_t_start = construct_positive_covariates(row$age, row$is_male, row$sector),
      covariates_at_t_start = construct_positive_covariates(row$age, row$is_male, row$sector),
      cumulative_time_at_t_start = as.integer(T - as.Date(row$first_positive_test_date)),
      state_at_t_start = POSITIVE
    )
  })
  
  positive_with_history_patient_objects[[1]]
  
  all_patients = c(all_patients, positive_with_history_patient_objects)
}


## ------------------------------------------------------------------------------------------------------------
if (!is.null(POSITIVE_NO_HISTORY_INPUT_CSV)) {
  
  positive_with_NO_history_df = read.csv(POSITIVE_NO_HISTORY_INPUT_CSV,stringsAsFactors=FALSE)
  positive_with_NO_history_df$patient_id <- NULL
  positive_with_NO_history_df$pred_id <- NULL
  stopifnot(positive_with_NO_history_df$first_positive_test_date > T)
  stopifnot(sum(complete.cases(positive_with_NO_history_df)) == nrow(positive_with_NO_history_df))
  print(paste('Group 2: Positive With NO History, ', nrow(positive_with_NO_history_df), ' patients', sep = ''))

  positive_with_NO_history_patient_objects = lapply(1:nrow(positive_with_NO_history_df), function(row_idx) {
    
    row = positive_with_NO_history_df[row_idx,]
    
    stopifnot(T <= row$first_positive_test_date)
    
    list(
      t_start = row$first_positive_test_date,
      covariates_at_t_start = construct_positive_covariates(row$age, row$is_male, row$sector),
      cumulative_time_at_t_start = 0,
      state_at_t_start = POSITIVE
    )
  })
  
  positive_with_NO_history_patient_objects[[1]]
  
  all_patients = c(all_patients, positive_with_NO_history_patient_objects)
}


## ------------------------------------------------------------------------------------------------------------
if (!is.null(HOSPITALIZED_WITH_HISTORY_INPUT_CSV)) {
  
  hospitalized_with_history_df = read.csv(HOSPITALIZED_WITH_HISTORY_INPUT_CSV,stringsAsFactors=FALSE)  
  hospitalized_with_history_df$patient_id <- NULL
  hospitalized_with_history_df$pred_id <- NULL
  stopifnot(hospitalized_with_history_df$date_of_hospitalization <= T)
  stopifnot(sum(complete.cases(hospitalized_with_history_df)) == nrow(hospitalized_with_history_df))
  stopifnot(sum(hospitalized_with_history_df$state_at_t_start %in% c(MILD, MODERATE) ) == 0)
  stopifnot(sum(hospitalized_with_history_df$state_at_hospitalization %in% c(MILD, MODERATE) ) != 0)

  print(paste('Group 3: Hospitalized with History, ', nrow(hospitalized_with_history_df), ' patients', sep = ''))

  hospitalized_with_history_patient_objects = lapply(1:nrow(hospitalized_with_history_df), function(row_idx) {
    
    row = hospitalized_with_history_df[row_idx,]
    
    covariates_at_entry_to_state_at_t_start = construct_hospitalized_patient_covariates(
      age = row$age, 
      is_male = row$is_male, 
      state_at_hospitalization = row$state_at_hospitalization, 
      
      # the cumulative time when entered this state (that was the moment covariates were updated):
      cumulative_time_since_hospitalization = as.integer( (T-as.Date(row$date_of_hospitalization)) - row$days_spent_in_current_state),
      was_severe = row$was_severe,
      was_critical = row$was_critical
    )
    
    list(
      t_start = T,
      covariates_at_t_start = covariates_at_entry_to_state_at_t_start,
      
      # the time since hospitalization including the time in current state (predictions are made based on this time forward):
      cumulative_time_at_t_start = as.integer(T-as.Date(row$date_of_hospitalization)), 
      state_at_t_start = row$state_at_t_start,
      date_of_hospitalization = as.Date(row$date_of_hospitalization)
    )
  })
  
  hospitalized_with_history_patient_objects[[2]]
 
  all_patients = c(all_patients, hospitalized_with_history_patient_objects)
 
}



## ------------------------------------------------------------------------------------------------------------
if (!is.null(HOSPITALIZED_NO_HISTORY_INPUT_CSV)) {
  hospitalized_with_NO_history_df = read.csv(HOSPITALIZED_NO_HISTORY_INPUT_CSV,stringsAsFactors=FALSE)
  hospitalized_with_NO_history_df$patient_id <- NULL
  hospitalized_with_NO_history_df$pred_id <- NULL
  stopifnot(hospitalized_with_NO_history_df$date_of_hospitalization > T)
  stopifnot(sum(complete.cases(hospitalized_with_NO_history_df)) == nrow(hospitalized_with_NO_history_df))
  stopifnot(sum(hospitalized_with_NO_history_df$state_at_hospitalization %in% c(MILD, MODERATE) ) != 0)
  
  print(paste('Group 4: Hospitalized NO History, ', nrow(hospitalized_with_NO_history_df), ' patients', sep = ''))

  
  hospitalized_with_NO_history_patient_objects = lapply(1:nrow(hospitalized_with_NO_history_df), function(row_idx) {
    
    row = hospitalized_with_NO_history_df[row_idx,]
    
    covariates_at_entry_to_state_at_t_start = construct_hospitalized_patient_covariates(
      age = row$age, 
      is_male = row$is_male, 
      state_at_hospitalization = row$state_at_hospitalization, 
      cumulative_time_since_hospitalization = 0,
      was_severe = 0,
      was_critical = 0
    )
    
    list(
      t_start = row$date_of_hospitalization,
      covariates_at_t_start = covariates_at_entry_to_state_at_t_start,
      cumulative_time_at_t_start = 0,
      state_at_t_start = ifelse(row$state_at_hospitalization %in% c(MILD, MODERATE), 
                                MILD_OR_MODERATE, 
                                row$state_at_hospitalization),
      date_of_hospitalization = as.Date(row$date_of_hospitalization)
    )
    
  })
  
  all_patients = c(all_patients, hospitalized_with_NO_history_patient_objects)
}



## ------------------------------------------------------------------------------------------------------------
length(all_patients)

start_time <- Sys.time()
if (TEST_RUN) all_patients = all_patients[1:1000]
tmp = mclapply(1:N_MONTE_CARLO_RUNS,function(i) {
    set.seed(i + SEED)
    all_patients = run_monte_carlo_per_patient(all_patients,
                                           1,
                                           MAX_PATH_LENGTH,
                                           n_cores = 1)
    saveRDS(all_patients, paste0(path, sprintf("/run_cache/%d.R",i)))
    return(1)
},mc.cores=N_CORES)
end_time <- Sys.time()
print("run monte carlo per patient")
print(end_time - start_time)


## ------------------------------------------------------------------------------------------------------------
true_counts_df = read.csv(TRUE_COUNTS_CSV,stringsAsFactors=FALSE)
true_counts_df = true_counts_df[(true_counts_df$date >= T) & (true_counts_df$date <= T + PREDICT_N_DAYS_AHEAD), ]
true_counts_df$severe_or_critical <- true_counts_df$severe + true_counts_df$critical


## ----echo=FALSE----------------------------------------------------------------------------------------------
start_time <- Sys.time()

plot_expected_number_of_patients_each_day(cache_dir,
                                          first_date = T,
                                          last_date = T + PREDICT_N_DAYS_AHEAD,
                                          states_to_include = c(SEVERE, CRITICAL),
                                          n_monte_carlo_runs = N_MONTE_CARLO_RUNS,
                                          true_counts = true_counts_df$severe_or_critical,
					  n_cores = N_CORES, prefix = PREFIX, output_dir=OUTPUT_DIR)


end_time <- Sys.time()
print("plot severe+critical")
print(end_time - start_time)

## ----echo=FALSE----------------------------------------------------------------------------------------------
# start_time <- Sys.time()

# plot_expected_number_of_patients_each_day(cache_dir,
#                                          first_date = T,
#                                          last_date = T + PREDICT_N_DAYS_AHEAD,
#                                          states_to_include = c(DECEASED),
#                                          n_monte_carlo_runs = N_MONTE_CARLO_RUNS,
#                                          true_counts = true_counts_df$cum_deaths_from_analysis_date,
# 					  n_cores = N_CORES, prefix = PREFIX, output_dir=OUTPUT_DIR)


# end_time <- Sys.time()
# print("plot deceased")
# print(end_time - start_time)

## ----echo=FALSE----------------------------------------------------------------------------------------------
# start_time <- Sys.time()

# plot_expected_number_of_patients_each_day(cache_dir,
#                                          first_date = T,
#                                          last_date = T + PREDICT_N_DAYS_AHEAD,
#                                          states_to_include = c(CRITICAL),
#                                          n_monte_carlo_runs = N_MONTE_CARLO_RUNS,
#                                          true_counts = true_counts_df$critical,
#					  n_cores = N_CORES, prefix = PREFIX, output_dir=OUTPUT_DIR)

# end_time <- Sys.time()
# print("plot critical")
# print(end_time - start_time)


## ------------------------------------------------------------------------------------------------------------
# start_time <- Sys.time()

# plot_expected_number_of_patients_each_day(cache_dir,
#                                          first_date = T,
#                                          last_date = T + PREDICT_N_DAYS_AHEAD,
#                                          states_to_include = c(SEVERE),
#                                          n_monte_carlo_runs = N_MONTE_CARLO_RUNS,
#                                          true_counts = true_counts_df$severe,
#					  n_cores = N_CORES, prefix = PREFIX, output_dir=OUTPUT_DIR)

# end_time <- Sys.time()
# print("plot severe")
# print(end_time - start_time)

## ------------------------------------------------------------------------------------------------------------
# start_time <- Sys.time()

# plot_expected_number_of_patients_each_day(cache_dir,
#                                          first_date = T,
#                                          last_date = T + PREDICT_N_DAYS_AHEAD,
#                                          states_to_include = c(MILD_OR_MODERATE),
#                                          n_monte_carlo_runs = N_MONTE_CARLO_RUNS,
#                                          true_counts = true_counts_df$moderate + true_counts_df$mild,
# 					  n_cores = N_CORES, prefix = PREFIX, output_dir=OUTPUT_DIR)
# end_time <- Sys.time()
# print("plot mild moderate")
# print(end_time - start_time)
