library(data.table)
library(ggplot2)
library(dplyr)
library(ggh4x)

folder_results <- "analysis/sim_traits_gwas/cor_curve"

# load data
cor_results <- paste0(folder_results, "/summary_cor_sim-vs-real-vs-pred.txt")
cor_results <- fread(cor_results, header = TRUE, data.table = FALSE)


#### cor across envs ----

cor_results %>%
  filter(avg_rank == "all", cor_env_avg == "across", model == "A") %>% 
  group_by(trait, n_marker, effect, cor_type) %>% 
  summarize(cor_across_envs = mean(cor_value, na.rm = TRUE),
            se_across_envs = sd(cor_value) / sqrt(n())) %>% 
  replace(is.na(.), NA) %>% 
  ggplot(aes(x = n_marker, y = cor_across_envs, color = factor(effect))) +
  facet_grid(cor_type ~ trait) +
  geom_point() +
  geom_smooth(se = FALSE) +
  geom_errorbar(aes(ymin = cor_across_envs - se_across_envs, ymax = cor_across_envs + se_across_envs),
                width = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(title = "cor across envs", subtitle = "add")
# ggsave(filename = paste0(folder_results, "/cor_sim-real-pred_across_envs.add.pdf"),
#        device = "pdf", width = 8, height = 8)

cor_results %>%
  filter(avg_rank == "all", cor_env_avg == "across", model == "D") %>% 
  group_by(trait, n_marker, effect, cor_type) %>% 
  summarize(cor_across_envs = mean(cor_value, na.rm = TRUE),
            se_across_envs = sd(cor_value) / sqrt(n())) %>% 
  replace(is.na(.), NA) %>% 
  ggplot(aes(x = n_marker, y = cor_across_envs, color = factor(effect))) +
  facet_grid(cor_type ~ trait) +
  geom_point() +
  geom_smooth(se = FALSE) +
  geom_errorbar(aes(ymin = cor_across_envs - se_across_envs, ymax = cor_across_envs + se_across_envs),
                width = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(title = "cor across envs", subtitle = "dom")
# ggsave(filename = paste0(folder_results, "/cor_sim-real-pred_across_envs.dom.pdf"),
#        device = "pdf", width = 8, height = 8)


cor_results %>%
  filter(avg_rank == "all", cor_env_avg == "across", model == "A", effect == "0.1") %>% 
  group_by(trait, n_marker, cor_type) %>% 
  summarize(cor_across_envs = mean(cor_value, na.rm = TRUE),
            se_across_envs = sd(cor_value) / sqrt(n())) %>% 
  replace(is.na(.), NA) %>% 
  ggplot(aes(x = n_marker, y = cor_across_envs, color = cor_type)) +
  facet_wrap( ~ trait) +
  geom_point() +
  geom_smooth(se = FALSE) +
  geom_errorbar(aes(ymin = cor_across_envs - se_across_envs, ymax = cor_across_envs + se_across_envs),
                width = 0.2) +
  coord_cartesian(ylim = c(-0.5, 1)) +
  labs(title = "cor across envs", subtitle = "add")

cor_results %>%
  filter(avg_rank == "all", cor_env_avg == "across", model == "A", effect == "1") %>% 
  group_by(trait, n_marker, cor_type) %>% 
  summarize(cor_across_envs = mean(cor_value, na.rm = TRUE),
            se_across_envs = sd(cor_value) / sqrt(n())) %>% 
  replace(is.na(.), NA) %>% 
  ggplot(aes(x = n_marker, y = cor_across_envs, color = cor_type)) +
  facet_wrap( ~ trait) +
  geom_point() +
  geom_smooth(se = FALSE) +
  geom_errorbar(aes(ymin = cor_across_envs - se_across_envs, ymax = cor_across_envs + se_across_envs),
                width = 0.2) +
  coord_cartesian(ylim = c(-0.5, 1)) +
  labs(title = "cor across envs", subtitle = "add")


#### cor within envs ----

cor_results %>%
  filter(avg_rank == "all", cor_env_avg == "within", model == "A") %>%
  group_by(trait, n_marker, effect, cor_type) %>% 
  summarize(mean_cor_within_envs = mean(cor_value, na.rm = TRUE),
            se_within_envs = sd(cor_value) / sqrt(n())) %>% 
  replace(is.na(.), NA) %>% 
  ggplot(aes(x = n_marker, y = mean_cor_within_envs, color = factor(effect))) +
  facet_grid(cor_type ~ trait) +
  geom_point() +
  geom_smooth(se = FALSE) +
  geom_errorbar(aes(ymin = mean_cor_within_envs - se_within_envs, ymax = mean_cor_within_envs + se_within_envs),
                width = 0.2) +
  coord_cartesian(ylim = c(-0.5, 1)) +
  labs(title = "cor within envs", subtitle = "add")
# ggsave(filename = paste0(folder_results, "/cor_sim-real-pred_within_envs.add.pdf"),
#        device = "pdf", width = 8, height = 8)

cor_results %>%
  filter(avg_rank == "all", cor_env_avg == "within", model == "D") %>%
  group_by(trait, n_marker, effect, cor_type) %>% 
  summarize(mean_cor_within_envs = mean(cor_value, na.rm = TRUE)) %>% 
  replace(is.na(.), NA) %>% 
  ggplot(aes(x = n_marker, y = mean_cor_within_envs, color = factor(effect))) +
  facet_grid(cor_type ~ trait) +
  geom_point() +
  geom_smooth(se = FALSE) +
  coord_cartesian(ylim = c(-0.5, 1)) +
  labs(title = "cor within envs", subtitle = "dom")
# ggsave(filename = paste0(folder_results, "/cor_sim-real-pred_within_envs.dom.pdf"),
#        device = "pdf", width = 8, height = 8)



#### cor within envs (avg effect results) ----

cor_results %>%
  filter(avg_rank == "all", cor_env_avg == "within", model == "A") %>%
  group_by(trait, n_marker, cor_type) %>% 
  summarize(mean_cor_within_envs = mean(cor_value, na.rm = TRUE),
            se_within_envs = sd(cor_value) / sqrt(n())) %>% 
  replace(is.na(.), NA) %>% 
  ggplot(aes(x = n_marker, y = mean_cor_within_envs, color = cor_type)) +
  facet_wrap(~ trait) +
  geom_point() +
  geom_smooth(se = FALSE) +
  geom_errorbar(aes(ymin = mean_cor_within_envs - se_within_envs, ymax = mean_cor_within_envs + se_within_envs),
                width = 0.2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "cor within envs", subtitle = "add")

cor_results %>%
  filter(avg_rank == "all", cor_env_avg == "within", model == "D") %>%
  group_by(trait, n_marker, cor_type) %>% 
  summarize(mean_cor_within_envs = mean(cor_value, na.rm = TRUE),
            se_within_envs = sd(cor_value) / sqrt(n())) %>% 
  replace(is.na(.), NA) %>% 
  ggplot(aes(x = n_marker, y = mean_cor_within_envs, color = cor_type)) +
  facet_wrap(~ trait) +
  geom_point() +
  geom_smooth(se = FALSE) +
  geom_errorbar(aes(ymin = mean_cor_within_envs - se_within_envs, ymax = mean_cor_within_envs + se_within_envs),
                width = 0.2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "cor within envs", subtitle = "dom")



#### cor per envs ----

cor_results %>%
  filter(avg_rank == "top", !is.na(env)) %>% 
  ggplot(aes(x = n_marker, y = cor_within_envs, color = factor(effect))) +
  facet_grid(env ~ trait) +
  geom_point() +
  geom_smooth() +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  labs(title = "cor per envs", subtitle = "top")
ggsave(filename = paste0(folder_results, "/cor_per_env.top.pdf"),
       device = "pdf", width = 8, height = 10)

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


