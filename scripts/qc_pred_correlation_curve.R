library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

folder_results <- "analysis/sim_traits_gwas/cor_curve"

traits <- c("YLD", "EHT", "PHT", "Moisture")
# avg_ranks <- c("top", "all")
avg_ranks <- c("all")
n_markers <- c(1, 3, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
               200, 300, 400, 500, 600, 700, 800, 900, 1000)
#n_markers <- c(1, 3, 5, 10, 30, 50, 100, 300, 500, 1000)
effects <- c(0.1, 1)

pred_iters <- 3
models <- c("A", "D")
var_structures <- "Gstr-diag_Rstr-us"

env_names <- list("EHT" = list("env1" = "COR19", "env2" = "MIN19",
                               "env3" = "MIN20","env4" = "URB19"),
                  "PHT" = list("env1" = "COR19", "env2" = "COR20", "env3" = "MIN19", 
                               "env4" = "MIN20", "env5" = "URB19"),
                  "YLD" = list("env4" = "COR19", "env5" = "COR20", "env6" = "MIN19",
                               "env7" = "MIN20", "env8" = "SYN19", "env9" = "SYN20"),
                  "Moisture" = list("env5" = "COR19", "env6" = "COR20", "env7" = "MIN19",
                                    "env8" = "MIN20", "env9" = "SYN19", "env10" = "SYN20"))


# create empty df to store results
cor_results <- data.frame(stringsAsFactors = FALSE)

for (trait in traits) {
  
  print(trait)
  
  for (avg_rank in avg_ranks) {
    for (n_marker in n_markers) {
      for (effect in effects) {
        
        # get folder with sim data
        folder_sim <- paste0("analysis/sim_traits_gwas/cor_curve/", trait,
                         "/traits/avg_rank_", avg_rank, "/n_markers_", n_marker,
                         "/effects_", effect)
        # get filename
        sim_real_data_file <- paste0(folder_sim, "/sim_vs_real_data.txt")
        
        for (model in models) {
          for (var_str in var_structures) {
            for (pred_iter in 1:pred_iters) {
              for (CV in c("CV2", "CV1")) {
                
                # get folder with prediction accuracy results
                folder_accuracy <- paste0(folder_sim, "/prediction_results/iter",
                                          pred_iter, "/model-", model, "/", var_str)
                # get filenames
                gebv_file <- list.files(folder_accuracy, pattern = paste0("GEBVs.", CV, ".txt"), full.names = TRUE)
                
                if (file.exists(sim_real_data_file) & length(gebv_file) > 0) {
                  
                  # load sim data
                  sim_real_data <- fread(sim_real_data_file, header = TRUE, data.table = FALSE)
                  # load prediction results
                  gebv_data <- fread(gebv_file, header = TRUE, data.table = FALSE)
                  # transform to long format and take average per iteration in each env
                  gebv_data <- gebv_data %>% 
                    pivot_longer(-genotype, names_to = "env_iter", values_to = "pred_pheno") %>% 
                    separate(env_iter, into = c("env", "iter"), sep = "_") %>% 
                    group_by(genotype, env) %>% 
                    summarize(pred_pheno_avg = mean(pred_pheno, na.rm = TRUE))
                  # rename envs
                  gebv_data$env <- sapply(gebv_data$env, function(x) env_names[[trait]][[x]])
                  
                  # filter sim data to have the same envs as pred data
                  sim_real_data_filtered <- subset(sim_real_data, env %in% as.character(unlist(env_names[[trait]])))
                  colnames(sim_real_data_filtered)[1] <- "genotype"
                  
                  # merge data
                  sim_real_pred_data <- merge(x = sim_real_data_filtered, y = gebv_data,
                                              by = c("genotype", "env"), all = TRUE)
                  
                  # get correlations across and within envs
                  cor_summary <- sim_real_pred_data %>% 
                    mutate(sim_vs_real__across = cor(sim_value_avg, real_pheno, method = "spearman", use = "complete.obs"),
                           sim_vs_pred__across = cor(sim_value_avg, pred_pheno_avg, method = "spearman", use = "complete.obs"),
                           real_vs_pred__across = cor(real_pheno, pred_pheno_avg, method = "spearman", use = "complete.obs")) %>%
                    group_by(env) %>%
                    summarize(sim_vs_real__across = unique(sim_vs_real__across),
                              sim_vs_pred__across = unique(sim_vs_pred__across),
                              real_vs_pred__across = unique(real_vs_pred__across),
                              sim_vs_real__within = cor(sim_value_avg, real_pheno, method = "spearman", use = "complete.obs"),
                              sim_vs_pred__within = cor(sim_value_avg, pred_pheno_avg, method = "spearman", use = "complete.obs"),
                              real_vs_pred__within = cor(real_pheno, pred_pheno_avg, method = "spearman", use = "complete.obs")) %>% 
                    ungroup() %>% 
                    pivot_longer(-env, names_to = c("cor_type", "cor_env_avg"), names_sep = "__", values_to = "cor_value") 
                  
                  # append to main df
                  cor_results <- rbind(cor_results,
                                       data.frame(trait = trait, avg_rank = avg_rank,
                                                  n_marker = n_marker, effect = effect,
                                                  model = model, var_str = var_str, 
                                                  pred_iter = pred_iter, CV = CV,
                                                  cor_summary))
                  
                } else {
                  
                  # append empty df
                  cor_results <- rbind(cor_results,
                                       data.frame(trait = trait, avg_rank = avg_rank,
                                                  n_marker = n_marker, effect = effect,
                                                  model = model, var_str = var_str, 
                                                  pred_iter = pred_iter, CV = CV,
                                                  env = NA, cor_type = NA,
                                                  cor_env_avg = NA, cor_value = NA))
                  
                }
                
              }
            }
          }
        }
      }
    }
  }
}

# write files
fwrite(cor_results, file = paste0(folder_results, "/summary_cor_sim-vs-real-vs-pred.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
