library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

usage <- function() {
  cat("
description: calculate correlations between simulated and empirical traits, and
             summarize percent variance explained from ANOVA of each simulated trait.

usage: Rscript qc_sim_correlation_curve.R [folder_base] [traits] [outfile_name] [...]

positional arguments:
  folder_base                   path to folder with results of cross validation
  traits                        comma-separated list of trait names to summarize
  n_markers                     comma-separated list of number of causative variants for each trait
  effect_sizes                  comma-separated list of effect sizes of causative variants for each trait

optional argument:
  --help                        show this helpful message
  --avg-rank-markers=[VALUE]    marker selection was based on average rank among 'all' markers or
                                'top' markers in each environment? (default: 'all')

"
  )
}

getArgValue <- function(arg) {
  
  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(arg[2])
  
}



#### command line options ----

# set default
avg_rank_markers <- "all"

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 4) stop(usage(), "missing positional argument(s)")

if (length(args) > 4) {
  
  opt_args <- args[-1:-4]
  opt_args_allowed <- c("--avg-rank-markers")
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

if (!avg_rank_markers %in% c("top", "all")) {
  stop("Optional argument '--avg-rank-markers' should be 'all' or 'top'")
}

# get positional arguments
folder_base <- args[1]
traits <- unlist(strsplit(args[2], split = ","))
n_markers <- unlist(strsplit(args[3], split = ","))
effect_sizes <- unlist(strsplit(args[4], split = ","))



#### summarize results ----

# create empty df to store results
cor_results <- data.frame(stringsAsFactors = FALSE)
across_envs_results <- data.frame(stringsAsFactors = FALSE)
within_envs_results <- data.frame(stringsAsFactors = FALSE)
pve_results <- data.frame(stringsAsFactors = FALSE)

for (trait in traits) {
  
  print(trait)
  
  for (avg_rank in avg_rank_markers) {
    for (n_marker in n_markers) {
      for (effect in effect_sizes) {
        
        # get file
        folder <- paste0(folder_base, "/", trait, "/traits/michael_method",
                         "/n_markers_", n_marker,
                         "/effects_", effect)
        sim_real_data <- paste0(folder, "/sim_vs_real_data.txt")
        pve_data <- paste0(folder, "/simulated_traits_pve.txt")
        
        if (file.exists(sim_real_data)) {
          
          # load data
          sim_real_data <- fread(sim_real_data, header = TRUE, data.table = FALSE)
          pve_data <- fread(pve_data, header = TRUE, data.table = FALSE)
          
          # get correlations across and within envs
          cor_summary <- sim_real_data %>% 
            mutate(cor_across_envs = cor(sim_value_avg, real_pheno, method = "spearman", use = "complete.obs")) %>% 
            group_by(env) %>% 
            summarize(cor_across_envs = unique(cor_across_envs),
                      cor_within_envs = cor(sim_value_avg, real_pheno, method = "spearman", use = "complete.obs"))
          # append to main df
          cor_results <- rbind(cor_results,
                               data.frame(trait = trait, avg_rank = avg_rank,
                                          n_marker = n_marker, effect = effect,
                                          cor_summary))
        
          # get summary across envs
          summary_across_envs <- sim_real_data %>% 
            pivot_longer(-c(hybrid, env), names_to = "pheno", values_to = "value") %>% 
            mutate(pheno = if_else(pheno == "sim_value_avg", true = "sim", false = "real")) %>% 
            group_by(pheno) %>% 
            summarize(min = as.numeric(summary(value)[1]),
                      q1 = as.numeric(summary(value)[2]),
                      median = as.numeric(summary(value)[3]),
                      mean = as.numeric(summary(value)[4]),
                      q3 = as.numeric(summary(value)[5]),
                      max = as.numeric(summary(value)[6]),
                      n = n(),
                      sd = sd(value, na.rm = TRUE),
                      se = sd / sqrt(n)) %>% 
            pivot_longer(-pheno, names_to = "stat", values_to = "value")
          # append to main df
          across_envs_results <- rbind(across_envs_results,
                                       data.frame(trait = trait, avg_rank = avg_rank,
                                                  n_marker = n_marker, effect = effect,
                                                  summary_across_envs))
          
          # get summary within envs
          summary_within_envs <- sim_real_data %>% 
            pivot_longer(-c(hybrid, env), names_to = "pheno", values_to = "value") %>% 
            mutate(pheno = if_else(pheno == "sim_value_avg", true = "sim", false = "real")) %>% 
            group_by(pheno, env) %>% 
            summarize(min = as.numeric(summary(value)[1]),
                      q1 = as.numeric(summary(value)[2]),
                      median = as.numeric(summary(value)[3]),
                      mean = as.numeric(summary(value)[4]),
                      q3 = as.numeric(summary(value)[5]),
                      max = as.numeric(summary(value)[6]),
                      n = n(),
                      sd = sd(value, na.rm = TRUE),
                      se = sd / sqrt(n)) %>% 
            pivot_longer(-c(pheno, env), names_to = "stat", values_to = "value")
          # append to main df
          within_envs_results <- rbind(within_envs_results,
                                       data.frame(trait = trait, avg_rank = avg_rank,
                                                  n_marker = n_marker, effect = effect,
                                                  summary_within_envs))
          
          # get pve
          pve_results <- rbind(pve_results,
                               data.frame(trait = trait, avg_rank = avg_rank,
                                          n_marker = n_marker, effect = effect,
                                          pve_data))
          
          
        } else {
          
          # append empty df
          cor_results <- rbind(cor_results,
                               data.frame(trait = trait, avg_rank = avg_rank,
                                          n_marker = n_marker, effect = effect,
                                          env = NA, cor_across_envs = NA, cor_within_envs = NA))
          across_envs_results <- rbind(across_envs_results,
                                       data.frame(trait = trait, avg_rank = avg_rank,
                                                  n_marker = n_marker, effect = effect,
                                                  pheno = NA, stat = NA, value = NA))
          within_envs_results <- rbind(within_envs_results,
                                       data.frame(trait = trait, avg_rank = avg_rank,
                                                  n_marker = n_marker, effect = effect,
                                                  pheno = NA, env = NA, stat = NA, value = NA))
          pve_results <- rbind(pve_results,
                               data.frame(trait = trait, avg_rank = avg_rank,
                                          n_marker = n_marker, effect = effect,
                                          source = NA, pve = NA, pval = NA, signif = NA))
          
        }
        
      }
    }
  }
}
rm(sim_real_data, cor_summary, summary_across_envs, summary_within_envs)

# write files
fwrite(cor_results, file = paste0(folder_base, "/summary_cor_across-within_envs.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
fwrite(across_envs_results, file = paste0(folder_base, "/summary_stats_across_envs.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
fwrite(within_envs_results, file = paste0(folder_base, "/summary_stats_within_envs.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
fwrite(pve_results, file = paste0(folder_base, "/summary_pve.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)


#### debug ----

# folder_base <- "analysis/sim_traits"
# traits <- c("EHT", "PHT", "Moisture", "YLD")
# n_markers <- c(1, 3, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
#                200, 300, 400, 500, 600, 700, 800, 900, 1000)
# effect_sizes <- c(0.1, 1)
# # avg_rank_markers <- c("top", "all")
# avg_rank_markers <- "all"
