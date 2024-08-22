library(data.table)
library(plyr)
library(dplyr)
library(tidyr)

usage <- function() {
  cat("
description: summarize prediction accuracy results after k-fold cross validation of simulated traits.

usage: Rscript summarize_prediction_accuracies_cor-curve.R [folder_base] [traits] [outfile_name] [...]

positional arguments:
  folder_base                 path to folder with results of cross validation
  traits                      comma-separated list of trait names to summarize
  sim_n_markers               comma-separated list of number of markers used to simulate each trait
  sim_effects                 comma-separated list of effect size of markers used to simulate each trait
  outfile_name                output filename

optional argument:
  --help                      show this helpful message
  --models=[LIST]             comma-separated list of models
  --var-structures=[LIST]     comma-separated list of variance structures
  --pred-iters=[LIST]         number of predictors iterations used

"
  )
}

getArgValue <- function(arg) {
  
  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(arg[2])
  
}

getAccuracy <- function(folder_accuracy, trait, model, var_str, sim_n_marker, sim_effect, pred_iter) {
  
  # create empty df to store results of CVs
  df_results <- data.frame(stringsAsFactors = FALSE)
  
  for (cv in c("CV1", "CV2")) {
    
    # try to access file with results
    pred_accuracy <- try(fread(paste0(folder_accuracy, "/prediction_accuracy.", cv, ".txt"),
                               header = TRUE, data.table = FALSE))
    
    if (class(pred_accuracy) != "try-error") {
      
      # get avg accuracy across all environments
      pred_accuracy_means <- data.frame(t(rowMeans(pred_accuracy[, -1])))
      colnames(pred_accuracy_means) <- c("mean_accuracy_envs", "mean_se", "mean_lowerCI", "mean_upperCI")
      
      # keep accuracy about each environment as well -- just need to reformat
      pred_accuracy_envs <- pred_accuracy %>%
        pivot_longer(-stat, names_to = "env", values_to = "vals") %>%
        arrange(env) %>%
        unite(env:stat, col = "env_info", sep = "_") %>%
        t() %>% as.data.frame(stringsAsFactors = FALSE)
      colnames(pred_accuracy_envs) <- apply(pred_accuracy_envs[1, ], MARGIN = 1, function(x) as.character(gsub("_CI", "CI", x)))
      pred_accuracy_envs <- data.frame(apply(pred_accuracy_envs[-1, ], MARGIN = c(1, 2), function(x) as.numeric(x)))
      rownames(pred_accuracy_envs) <- NULL
      
      # add metadata
      df_results <- rbind(df_results,
                          data.frame(trait = trait, model = model, var_str = var_str,
                                     sim_n_marker = sim_n_marker, sim_effect = sim_effect,
                                     pred_iter = pred_iter, cv = cv,
                                     pred_accuracy_means, pred_accuracy_envs))
      
    }
  }
  
  # only return values if both CV1 and CV2 results exist
  if (NROW(df_results) == 2) {
    return(df_results)
  } else {
    # if model didn't converge in at least one of the CVs, add NAs
    df_results <- data.frame(trait = trait, model = model, var_str = var_str,
                             sim_n_marker = sim_n_marker, sim_effect = sim_effect,
                             pred_iter = pred_iter, cv = c("CV1", "CV2"))
    return(df_results)
  }
  
}



#### command line options ----

# set default
models <- 'A,D,AD'
var_str <- 'Gstr-idv_Rstr-idv,Gstr-diag_Rstr-us'
pred_iters <- 1

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 5) stop(usage(), "missing positional argument(s)")

if (length(args) > 5) {
  
  opt_args <- args[-1:-5]
  opt_args_allowed <- c("--models", "--var-structures", "--pred-iters")
  opt_args_requested <- as.character(sapply(opt_args, function(x) unlist(strsplit(x, split = "="))[1]))
  if (any(!opt_args_requested %in% opt_args_allowed)) stop(usage(), "wrong optional argument(s)")
  
  # change default based on the argument provided
  for (argument in opt_args_allowed) {
    if (any(grepl(argument, opt_args_requested))) {
      arg_name <- gsub("-", "_", gsub("--", "", argument))
      arg_value <- getArgValue(opt_args[grep(argument, opt_args)])
      assign(arg_name, arg_value)
    }
  }
  
}

# adjust format from optional arguments
models <- unlist(strsplit(models, split = ","))
var_structures <- unlist(strsplit(var_structures, split = ","))
pred_iters <- unlist(strsplit(pred_iters, split = ","))

# get positional arguments
folder_base <- args[1]
traits <- unlist(strsplit(args[2], split = ","))
sim_n_markers <- as.numeric(unlist(strsplit(args[3], split = ",")))
sim_effects <- as.numeric(unlist(strsplit(args[4], split = ",")))
outfile_name <- args[5]

#### summarize prediction results ----

# get results simulated traits
for (trait in traits) {

  # create empty df to store results
  prediction_results <- data.frame(stringsAsFactors = FALSE)
  
  for (model in models) {
    for (var_str in var_structures) {
      for (sim_n_marker in sim_n_markers) {
        for (sim_effect in sim_effects) {
          for (pred_iter in 1:pred_iters) {
            
            # get folder with prediction accuracy results
            folder_accuracy <- paste0(folder_base, "/sim_traits/", trait, "/traits/michael_method/n_markers_",
                                      sim_n_marker, "/effects_", sim_effect, "/prediction_results/iter",
                                      pred_iter, "/model-", model, "/", var_str)
            
            # set accuracy results
            prediction_results <- rbind.fill(prediction_results,
                                             getAccuracy(folder_accuracy, trait = trait, model = model,
                                                         var_str = var_str, sim_n_marker = sim_n_marker,
                                                         sim_effect = sim_effect, pred_iter = pred_iter))
            
          }
        }
      }
    }
  }
  
  # summarize results
  prediction_summary <- prediction_results %>%
    group_by(trait, model, var_str, sim_n_marker, sim_effect, cv) %>%
    summarize(iters = n(),
              envs = length(c_across(contains("_mean"))),
              envs_not_na = sum(!is.na(c_across(contains("_mean")))),
              across(c("mean_accuracy_envs", "mean_se", "mean_lowerCI", "mean_upperCI", contains("env")), ~ mean(.x, na.rm = TRUE)),
              accuracy_se = sd(c_across(contains("_mean")), na.rm = TRUE) / sqrt(iters * envs_not_na),
              accuracy_lowerCI = mean_accuracy_envs - (qt(p = 0.05/2, df = (iters * envs) - 1, lower.tail = FALSE) * accuracy_se),
              accuracy_upperCI = mean_accuracy_envs + (qt(p = 0.05/2, df = (iters * envs) - 1, lower.tail = FALSE) * accuracy_se))
  
  # round numbers
  prediction_results[, 8:NCOL(prediction_results)] <- apply(prediction_results[, 8:NCOL(prediction_results)],
                                                            MARGIN = 2, function(x) round(x, digits = 4))
  prediction_summary[, 10:NCOL(prediction_summary)] <- apply(prediction_summary[, 10:NCOL(prediction_summary)],
                                                             MARGIN = 2, function(x) round(x, digits = 4))
  # transform NaN to NA
  prediction_summary[sapply(prediction_summary, is.nan)] <- NA
  
  # write results
  fwrite(prediction_results, file = paste0(outfile_name, '.sim.iter1-', pred_iters, '.no-AD.', trait, '.full.txt'),
         quote = FALSE, sep = "\t", row.names = FALSE, na = NA)
  fwrite(prediction_summary, file = paste0(outfile_name, '.sim.iter1-', pred_iters, '.no-AD.', trait, '.summary.txt'),
         quote = FALSE, sep = "\t", row.names = FALSE, na = NA)
  
}


# get results real traits
for (trait in traits) {
  
  # create empty df to store results
  prediction_results <- data.frame(stringsAsFactors = FALSE)
  
  for (model in models) {
    for (var_str in var_structures) {
      for (pred_iter in 1:pred_iters) {
        
        # get folder with prediction accuracy results
        folder_accuracy <- paste0(folder_base, "/empirical_traits/", trait, "/prediction_results/iter",
                                  pred_iter, "/model-", model, "/", var_str)
        
        # set accuracy results
        prediction_results <- rbind.fill(prediction_results,
                                         getAccuracy(folder_accuracy, trait = trait, model = model,
                                                     var_str = var_str, sim_n_marker = NA,
                                                     sim_effect = NA, pred_iter = pred_iter))
        
      }
    }
  }
  
  # summarize results
  prediction_summary <- prediction_results %>%
    group_by(trait, model, var_str, sim_n_marker, sim_effect, cv) %>%
    summarize(iters = n(),
              envs = length(c_across(contains("_mean"))),
              envs_not_na = sum(!is.na(c_across(contains("_mean")))),
              across(c("mean_accuracy_envs", "mean_se", "mean_lowerCI", "mean_upperCI", contains("env")), ~ mean(.x, na.rm = TRUE)),
              accuracy_se = sd(c_across(contains("_mean")), na.rm = TRUE) / sqrt(iters * envs_not_na),
              accuracy_lowerCI = mean_accuracy_envs - (qt(p = 0.05/2, df = (iters * envs) - 1, lower.tail = FALSE) * accuracy_se),
              accuracy_upperCI = mean_accuracy_envs + (qt(p = 0.05/2, df = (iters * envs) - 1, lower.tail = FALSE) * accuracy_se))
  
  # round numbers
  prediction_results[, 8:NCOL(prediction_results)] <- apply(prediction_results[, 8:NCOL(prediction_results)],
                                                            MARGIN = 2, function(x) round(x, digits = 4))
  prediction_summary[, 10:NCOL(prediction_summary)] <- apply(prediction_summary[, 10:NCOL(prediction_summary)],
                                                             MARGIN = 2, function(x) round(x, digits = 4))
  # transform NaN to NA
  prediction_summary[sapply(prediction_summary, is.nan)] <- NA
  
  # write results
  fwrite(prediction_results, file = paste0(outfile_name, '.real.iter1-', pred_iters, '.no-AD.', trait, '.full.txt'),
         quote = FALSE, sep = "\t", row.names = FALSE, na = NA)
  fwrite(prediction_summary, file = paste0(outfile_name, '.real.iter1-', pred_iters, '.no-AD.', trait, '.summary.txt'),
         quote = FALSE, sep = "\t", row.names = FALSE, na = NA)
  
}



#### debug ----

# folder_base <- "analysis"
# # folder_base <- "analysis/real_traits/cor_curve"
# outfile_name <- "prediction_results"
# # outfile_name <- "analysis/real_traits/cor_curve/test_sim_prediction_results.txt"
# pred_iters <- 3
# # models <- c("A", "D", "AD")
# models <- c("A", "D")
# traits <- c("YLD", "EHT", "PHT", "Moisture")
# sim_n_markers <- c(5,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350)
# sim_effects <- c(1, 0.1)
# # var_structures <- c("Gstr-idv_Rstr-idv", "Gstr-diag_Rstr-us")
# var_structures <- "Gstr-diag_Rstr-diag"
