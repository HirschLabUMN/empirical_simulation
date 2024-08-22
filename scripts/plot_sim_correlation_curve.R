library(data.table)
library(ggplot2)
library(dplyr)
library(ggh4x)

usage <- function() {
  cat("
description: plot summary of correlation curve between simulated and empirical traits.

usage: Rscript plot_sim_correlation_curve.R [gwas_folder] [chr_info] [centromere_info] [...]

positional arguments:
  folder_results          folder with correlation curve results
  cor_results             file with sim vs empirical correlations
  across_envs_results     file with summary distribution stats across envs
  within_envs_results     file with summary distribution stats within envs
  pve_results             file with summary PVE from ANOVA

optional argument:
  --help              show this helpful message

"
  )
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 5) stop(usage(), "missing positional argument(s)")

# assign arguments to variables
folder_results <- args[1]
cor_results <- args[2]
across_envs_results <- args[3]
within_envs_results <- args[4]
pve_results <- args[5]
# folder_results <- "analysis/sim_traits"
# cor_results <- paste0(folder_results, "/summary_cor_across-within_envs.txt")
# across_envs_results <- paste0(folder_results, "/summary_stats_across_envs.txt")
# within_envs_results <- paste0(folder_results, "/summary_stats_within_envs.txt")
# pve_results <- paste0(folder_results, "/summary_pve.txt")

# load data
cor_results <- fread(cor_results, header = TRUE, data.table = FALSE)
across_envs_results <- fread(across_envs_results, header = TRUE, data.table = FALSE)
within_envs_results <- fread(within_envs_results, header = TRUE, data.table = FALSE)
pve_results <- fread(pve_results, header = TRUE, data.table = FALSE)


#### cor across envs ----

cor_results %>%
  filter(avg_rank == "all") %>% 
  group_by(trait, n_marker, effect) %>% 
  summarize(cor_across_envs = mean(cor_across_envs, na.rm = TRUE)) %>% 
  replace(is.na(.), NA) %>% 
  ggplot(aes(x = n_marker, y = cor_across_envs, color = factor(effect))) +
  facet_wrap(~ trait) +
  geom_point() +
  geom_smooth() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "cor across envs", subtitle = "all")
ggsave(filename = paste0(folder_results, "/cor_across_envs.all.pdf"),
       device = "pdf", width = 8, height = 8)



#### cor within envs ----

cor_results %>%
  filter(avg_rank == "all") %>% 
  group_by(trait, n_marker, effect) %>% 
  summarize(mean_cor_within_envs = mean(cor_within_envs, na.rm = TRUE)) %>% 
  replace(is.na(.), NA) %>% 
  ggplot(aes(x = n_marker, y = mean_cor_within_envs, color = factor(effect))) +
  facet_wrap(~ trait) +
  geom_point() +
  geom_smooth() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "mean cor within envs", subtitle = "all")
ggsave(filename = paste0(folder_results, "/cor_mean_within_envs.all.pdf"),
       device = "pdf", width = 8, height = 8)


#### cor per envs ----

cor_results %>%
  filter(avg_rank == "all", !is.na(env)) %>% 
  ggplot(aes(x = n_marker, y = cor_within_envs, color = factor(effect))) +
  facet_grid(env ~ trait) +
  geom_point() +
  geom_smooth() +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  labs(title = "cor per envs", subtitle = "all")
ggsave(filename = paste0(folder_results, "/cor_per_env.all.pdf"),
       device = "pdf", width = 8, height = 10)


#### summary across envs ----

across_envs_results %>% 
  filter(avg_rank == "all", !is.na(pheno),
         stat %in% c("min", "q1", "mean", "median", "q3", "max")) %>% 
  mutate(stat = factor(stat, levels = c("min", "q1", "mean", "median", "q3", "max"))) %>% 
  ggplot(aes(x = n_marker, y = value, color = stat)) +
  facet_grid(trait ~ effect + pheno, scales = "free_y") +
  # geom_point() +
  geom_smooth() +
  labs(title = "distribution summary across envs", subtitle = "all") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(filename = paste0(folder_results, "/dist_summary_across_envs.all.pdf"),
       device = "pdf", width = 8, height = 8)


#### summary within envs ----

for (trait_name in c("EHT", "PHT", "Moisture", "YLD")) {

  within_envs_results %>% 
    filter(avg_rank == "all", !is.na(pheno), trait == trait_name,
           stat %in% c("min", "q1", "mean", "median", "q3", "max")) %>% 
    mutate(stat = factor(stat, levels = c("min", "q1", "mean", "median", "q3", "max"))) %>% 
    ggplot(aes(x = n_marker, y = value, color = stat)) +
    # facet_grid(trait ~ effect + stat, scales = "free_y") +
    facet_nested(env ~ effect + pheno, scales = "free_y",
                 nest_line = element_line(linetype = 1, color = "gray80")) +
    # geom_point() +
    geom_smooth() +
    labs(title = paste0("distribution summary within ", trait_name, " envs"), subtitle = "all") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  ggsave(filename = paste0(folder_results, "/dist_summary_within_envs.", trait_name, ".all.pdf"),
         device = "pdf", width = 8, height = 10)
  
}



#### summary pve ----

# create empty df to store results
pve_real <- data.frame(stringsAsFactors = FALSE)

for (trait in c("EHT", "PHT", "Moisture", "YLD")) {
  # get filename
  pve_file <- paste0("data/NIFA_CompleteDataset.", trait, ".pve.txt")
  # load pve for trait
  pve_file <- fread(pve_file, header = TRUE, data.table = FALSE)
  # append to main df
  pve_real <- rbind(pve_real, data.frame(trait = trait, pve_file, stringsAsFactors = FALSE))
}
rm(pve_file)

pve_real$source <- factor(pve_real$source,
                          levels = c("genotype", "environment", "rep:environment", "genotype:environment", "residual"),
                          labels = c("Genotype", "Env", "Env/Rep", "Genotype x Env", "Residual"))

pve_results %>% 
  filter(avg_rank == "all", !is.na(source)) %>% 
  mutate(source = factor(source, levels = c("genotype", "environment", "rep:environment", "genotype:environment", "residual"),
                         labels = c("Genotype", "Env", "Env/Rep", "Genotype x Env", "Residual"))) %>% 
  ggplot(aes(x = n_marker, y = pve, color = factor(effect))) +
  facet_grid(source ~ trait) +
  # geom_point() +
  geom_smooth() +
  coord_cartesian(ylim = c(0, 1)) +
  geom_hline(data = pve_real, aes(yintercept = pve)) +
  labs(title = "anova pve", subtitle = "all")
ggsave(filename = paste0(folder_results, "/pve_summary.all.pdf"),
       device = "pdf", width = 8, height = 8)
