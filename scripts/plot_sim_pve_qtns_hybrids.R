library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gtools)

usage <- function() {
  cat("
description: plot percent variance explained per QTN for simulated traits.

usage: Rscript plot_sim_pve_qtns_hybrids.R [dir_sim_traits] [sv_list] [...]

positional arguments:
  dir_sim_traits              path to folder containing simulated traits from simplePHENOTYPES
  sv_list                     single-column file containing only SV IDs

optional argument:
  --help                      show this helpful message

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

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# get positional arguments
sim_folder <- args[1]
sv_list <- args[2]
# sim_folder <- "analysis/sim_traits_gwas/EHT"
# sv_list <- "data/SVs_IDs_poly.txt"



#### load info about QTNs ----

# load file with simulations
qtn_info <-  data.frame(stringsAsFactors = FALSE)
for (effect_type in c("Additive", "Dominance")) {
  
  # find file in folder
  file_qtn_info <- list.files(sim_folder, pattern = effect_type, full.names = TRUE, ignore.case = TRUE)
  
  if (length(file_qtn_info) > 0) {
    
    # load file
    qtn_info_effect <- fread(file_qtn_info, header = TRUE, data.table = FALSE)
    # format columns
    colnames(qtn_info_effect)[3] <- "env"
    colnames(qtn_info_effect)[4] <- "effect_size"
    qtn_info_effect$env <- gsub("trait_", "env", qtn_info_effect$env)
    qtn_info_effect <- cbind(qtn_info_effect[, 1:3],
                             effect_type = tolower(effect_type),
                             qtn_info_effect[, 4:ncol(qtn_info_effect)],
                             stringsAsFactors = FALSE)
    # merge add and dom files
    qtn_info <- rbind(qtn_info, qtn_info_effect)
    rm(qtn_info_effect)
    
  }
}

# load SV information
SVs <- fread(sv_list, header = FALSE, data.table = FALSE)
SVs <- SVs[, 1]

# identify which QTNs are SVs
qtn_info$var <- case_when(qtn_info$snp %in% SVs ~ "SV", !qtn_info$snp %in% SVs ~ "SNP")
qtn_info <- relocate(var, .before = snp, .data = qtn_info)



#### variance explained by qtn ----

# find files with percent variation explained
files_var_qtn <- list.files(sim_folder, pattern = "PVE", full.names = TRUE)
# create empty df to store results
pve_qtns <- data.frame()
# count number of environments
n_envs <- length(unique(qtn_info$env))

for (environment in 1:n_envs) {
  
  # get inital qtn effects for this environment
  qtn_info_env <- subset(qtn_info, env == paste0("env", environment))
  
  # create empty variables
  var_qtn_env <- data.frame(stringsAsFactors = FALSE)
  n_qtns_env <- 0
  for (qtn_effect in c("ADD", "DOM")) {
    
    # load pve of qtns for that specific effect type
    var_qtn_env_effect <- files_var_qtn[grep(paste0("PVE_of_", qtn_effect, "_QTNs_Trait_", environment, ".txt"), files_var_qtn)]
    
    if (length(var_qtn_env_effect) > 0) {
      
      # read file
      var_qtn_env_effect <- fread(var_qtn_env_effect, header = FALSE, data.table = FALSE)
      # count how many qtns
      n_qtns_env <- n_qtns_env + (ncol(var_qtn_env_effect) - 1)
      # merge tables
      if (nrow(var_qtn_env) == 0) {
        var_qtn_env <- rbind(var_qtn_env, var_qtn_env_effect)
      } else {
        var_qtn_env <- cbind(var_qtn_env, var_qtn_env_effect[, -1])
        var_qtn_env[1, -1] <- paste0("QTN_", 1:(n_qtns_env))
      }
      rm(var_qtn_env_effect)
      
    } 
  }
  
  # add qtn names and type
  var_qtn_env <- rbind(c("names", qtn_info_env$snp), c("type", qtn_info_env$effect_type), c("var", qtn_info_env$var),
                       c("maf", qtn_info_env$maf), c("effect", qtn_info_env$effect_size), var_qtn_env)
  # format data frame
  var_qtn_env <- data.frame(env = environment, t(var_qtn_env[, -1]), stringsAsFactors = FALSE, row.names = NULL)
  colnames(var_qtn_env)[1:7] <- c("env", "qtn_name", "qtn_type", "var", "qtn_maf", "qtn_effect", "qtn_number")
  colnames(var_qtn_env)[8:NCOL(var_qtn_env)] <- paste0("rep", 1:length(8:NCOL(var_qtn_env)))
  var_qtn_env <- var_qtn_env %>% 
    pivot_longer(contains("rep"), names_to = "rep", values_to = "pve") %>%
    mutate(env = paste0("env", env),
           qtn_number = as.numeric(gsub("QTN_", "", qtn_number)),
           qtn_effect = as.numeric(qtn_effect),
           qtn_maf = as.numeric(qtn_maf),
           pve = as.numeric(pve)) %>% 
    relocate(rep, .after = env)
  # append to main df
  pve_qtns <- rbind(pve_qtns, var_qtn_env)
  
}
rm(var_qtn_env, qtn_info_env)

# calculate correlation of initial qtn effects and final variance explained
pearson <- round(cor(pve_qtns$qtn_effect, pve_qtns$pve, method = "pearson", use = "complete.obs"), digits = 2)
spearman <- round(cor(pve_qtns$qtn_effect, pve_qtns$pve, method = "spearman", use = "complete.obs"), digits = 2)

# reorder envs to plot
pve_qtns$env <- factor(pve_qtns$env, levels = mixedsort(unique(pve_qtns$env)))
pve_qtns$var <- factor(pve_qtns$var, levels = c("SV", "SNP"))
pve_qtns$qtn_type <- factor(pve_qtns$qtn_type, levels = c("additive", "dominance"))

# plot distribution of initial effects -- check if add and dom effects have same distribution
plot_qtn_eff_dist <- pve_qtns %>% 
  group_by(env, qtn_name, qtn_type) %>% 
  summarize(qtn_effect = mean(qtn_effect)) %>%
  ungroup() %>% 
  ggplot(aes(x = qtn_effect, fill = qtn_type)) +
  geom_density() +
  facet_grid(qtn_type ~ env, scales = "free") +
  labs(x = "QTN effect", y = "density")

ggsave(filename = paste0(sim_folder, "/dist_initial_qtn_effects.pdf"),
       plot = plot_qtn_eff_dist, device = "pdf", width = 15)

# boxplot of initial effects
plot_qtn_eff_boxplot <- pve_qtns %>% 
  group_by(env, qtn_name, qtn_type) %>% 
  summarize(qtn_effect = mean(qtn_effect)) %>%
  ungroup() %>% 
  ggplot() +
  geom_boxplot(aes(x = env, y = qtn_effect, color = qtn_type), show.legend = FALSE) +
  facet_grid(~ qtn_type) +
  labs(x = "QTN effect", y = "density")

ggsave(filename = paste0(sim_folder, "/boxplot_initial_qtn_effects.pdf"),
       plot = plot_qtn_eff_boxplot, device = "pdf", width = 15)

# plot percent variance explained by qtns in each environment
plot_pve_qtns <- ggplot(pve_qtns, aes(x = as.factor(qtn_number), y = pve, color = interaction(var, qtn_type))) +
  geom_violin() +
  geom_text(aes(y = 0, label = qtn_type), vjust = 1.5, size = 2) +
  facet_wrap(~env) +
  labs(x = "QTNs", y = "variance explained",
       caption = paste0("Correlation between QTN effects and variance explained: ",
                        pearson, " (Pearson) / ", spearman, " (Spearman)")) +
  # theme(axis.text.x = element_text(size = 5))
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(filename = paste0(sim_folder, "/var_explained_by_each_qtn.pdf"),
       plot = plot_pve_qtns, device = "pdf", width = 25)

# plot percent variance explained from all qtns in each environment
# (separated by qtn type)
total_pve <- pve_qtns %>%
  group_by(env, rep, qtn_type, var) %>%
  summarize(var_all_qtns = sum(pve))

plot_total_pve <- ggplot(total_pve, aes(x = env, y = var_all_qtns, color = interaction(var, qtn_type))) +
  geom_boxplot() +
  labs(x = "environment", y = "variance explained (all QTNs)")

ggsave(filename = paste0(sim_folder, "/var_explained_all_qtns.pdf"),
       plot = plot_total_pve, device = "pdf", width = 12)

# # plot percent variance explained by MAF
# plot_pve_maf <- ggplot(pve_qtns, aes(x = as.numeric(qtn_maf), y = as.numeric(pve), color = rep)) +
#   geom_point() +
#   facet_wrap(~env, nrow = length(unique(pve_qtns$env)), scales = "free_y") +
#   labs(x = "MAF", y = "PVE")
# ggsave(filename = paste0(sim_folder, "/var_explained_by_qtn_maf_v2.pdf"),
#        plot = plot_pve_maf, device = "pdf", width = 10)

# get average variance across reps
pve_qtns_avg <- pve_qtns %>% 
  group_by(env, qtn_name, qtn_type, var, qtn_maf, qtn_effect, qtn_number) %>% 
  summarize(var_exp_mean = mean(pve), var_exp_se = sd(pve)/sqrt(n())) %>% 
  ungroup() %>% 
  arrange(env, qtn_number)

# write summary of PVE per QTNs (averaged across reps)
fwrite(x = pve_qtns_avg, file = paste0(sim_folder, "/summary_var_explained_per_qtn.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
