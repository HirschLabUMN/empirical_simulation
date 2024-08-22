library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

usage <- function() {
  cat("
description: plot prediction accuracy results after k-fold cross validation of simulated traits.

usage: Rscript plot_prediction_accuracy_cor-curve.R [prediction_summary_file] [...]

positional arguments:
  prediction_summary_file       path to file with prediction summary

optional argument:
  --help                        show this helpful message
  --error-bars=[VALUE]          error bars can take one of the following values:
                                'CI_accuracy' (confidence intervals of the mean prediction accuracy; default),
                                'CI_error' (confidence intervals of the mean prediction error),
                                'SE_accuracy' (standard error of the mean prediction accuracy), or
                                'SE_error' (standard error of the mean prediction error)

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
error_bars <- "CI_accuracy"

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 1) stop(usage(), "missing positional argument(s)")

if (length(args) > 1) {
  
  opt_args <- args[-1:-2]
  opt_args_allowed <- c("--error-bars")
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

# get positional arguments
prediction_summary_file <- args[1]
outfolder <- args[2]


#### load prediction summary ----

# load data
prediction_summary <- fread(prediction_summary_file, header = TRUE, data.table = FALSE)

# reorder levels of marker type for better visualization
prediction_summary$var_str <- factor(prediction_summary$var_str,
                                     levels = unique(prediction_summary$var_str),
                                     labels = gsub("_", " / ", unique(prediction_summary$var_str)))

# create output folder if it doesn't exist
if (!dir.exists(outfolder)) dir.create(outfolder)



#### plot results ----

# plot
results_plot <- ggplot(data = prediction_summary,
                       aes(x = as.factor(sim_n_marker), y = mean_accuracy_envs, color = as.factor(sim_effect))) +
  geom_point() +
  # facet_grid(var_str ~ model) +
  facet_grid(cv ~ model) +
  coord_cartesian(ylim = c(-0.5, 1)) +
  scale_x_discrete(drop = FALSE) +
  labs(title = unique(prediction_summary$trait),
       x = "Number of markers used for simulations",
       y = "Prediction accuracy") +
  # scale_color_manual(values = c("#DE8282FF", "#AD0000FF")) +
  guides(color = guide_legend("Marker effect size")) +
  theme_bw() +
  theme(title = element_text(size = 12),
        text = element_text(size = 12),
        panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(vjust = -1.5),
        axis.title.y = element_text(vjust = 1.5),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# ggplot(data = prediction_summary,
#        aes(x = sim_n_marker, y = mean_accuracy_envs, color = as.factor(sim_effect))) +
#   geom_point() +
#   geom_smooth() +
#   facet_grid(cv ~ model) +
#   coord_cartesian(ylim = c(0, 1))

# add error bars
if (error_bars == "SE_accuracy") {
  results_plot <- results_plot +
    geom_errorbar(aes(ymin = mean_accuracy_envs - accuracy_se, ymax = mean_accuracy_envs + accuracy_se), width = 0.2)
}

if (error_bars == "CI_accuracy") {
  results_plot <- results_plot +
    geom_errorbar(aes(ymin = accuracy_lowerCI, ymax = accuracy_upperCI), width = 0.2)
}

if (error_bars == "SE_error") {
  results_plot <- results_plot +
    geom_errorbar(aes(ymin = mean_accuracy_envs - mean_SE, ymax = mean_accuracy_envs + mean_SE), width = 0.2)
}

if (error_bars == "CI_error") {
  results_plot <- results_plot +
    geom_errorbar(aes(ymin = mean_lowerCI, ymax = mean_upperCI), width = 0.2)
}

# print(results_plot)
# Sys.sleep(3)




###############

# NEED TO RUN PREDICTION IWITH 3 ITERS OF 500 MARKERS ON REAL TRAITS
# ONCE THAT'S DONE, ADD LINE WITH REAL PREDICTION VALUES HERE

# add prediction accuracy values from real data
prediction_summary_real <- fread(prediction_summary_real, header = TRUE, data.table = FALSE)

results_plot <- results_plot +
  geom_hline(data = prediction_summary_real, aes(yintercept = mean_accuracy_envs), linetype = "dashed")


# save plot
plot_name <- paste0(outfolder, "/pred_accuracy.", unique(prediction_summary$trait), ".pdf")
ggsave(results_plot, filename = plot_name, device = "pdf", width = 8, height = 6)
###############












#### choosing colors ----

# library(paletteer)
# library(tinter)
# paletteer_c(`"viridis::inferno"`, n = 20)
# for (col in c("#360961FF", "#8C2369FF", "#D84D3EFF", "#F1731DFF", "#FCA108FF")) {
#   print(prismatic::color(rev(tinter(col, direction = "tints"))[c(1,3)]))
# }



#### debug ----

# # prediction_summary_file <- "analysis/sim_traits_gwas/test_prediction_results.EHT.summary.txt"
# # prediction_summary_file <- "analysis/sim_traits_gwas/test_prediction_results.Moisture.summary.txt"
# # prediction_summary_file <- "analysis/sim_traits_gwas/test_prediction_results.PHT.summary.txt"
# # prediction_summary_file <- "analysis/sim_traits_gwas/test_prediction_results.YLD.summary.txt"
# # prediction_summary_file <- "analysis/sim_traits_gwas/all_qtls_gwas_reduced/results_prediction.iter1-3.no-AD.YLD.summary.txt"
# prediction_summary_file <- "analysis/sim_traits_gwas/cor_curve/test_sim_prediction_results.YLD.summary.txt"
# prediction_summary_real <- "analysis/real_traits/cor_curve/test_sim_prediction_results.YLD.summary.txt"
# prediction_summary_file <- "analysis/sim_traits_gwas/cor_curve/test_sim_prediction_results.EHT.summary.txt"
# prediction_summary_real <- "analysis/real_traits/cor_curve/test_sim_prediction_results.EHT.summary.txt"
# prediction_summary_file <- "analysis/sim_traits_gwas/cor_curve/test_sim_prediction_results.PHT.summary.txt"
# prediction_summary_real <- "analysis/real_traits/cor_curve/test_sim_prediction_results.PHT.summary.txt"
prediction_summary_file <- "analysis/sim_traits_gwas/cor_curve/test_sim_prediction_results.Moisture.summary.txt"
prediction_summary_real <- "analysis/real_traits/cor_curve/test_sim_prediction_results.Moisture.summary.txt"
outfolder <- "analysis/sim_traits_gwas/cor_curve/test_prediction_plots"
error_bars <- "CI_accuracy"
