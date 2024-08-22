library(data.table)
library(readxl)
library(gtools)
library(dplyr)
library(tidyr)
library(asreml)

usage <- function() {
  cat("
description: format empirical phenotypic from an excel sheet to multiple tab-delimited
             files (one for each trait and each environment).

usage: Rscript get_BLUEs_from_empirical-data.R [pheno_file] [...]

positional arguments:
  pheno_file          excel file with phenotypic data
  trait               name of trait to extract phenotypic values
  outfolder           name of output file
  folder_residuals    folder to save plots of residuals

optional argument:
  --help              show this helpful message
  --envs=VALUE        comma-separated list of environments to keep (default is to
                      keep all environments)
  --checks=VALUE      comma-separated list of checks to remove from data

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

# assign arguments to variables
pheno_file <- args[1]
trait <- args[2]
outfolder <- args[3]
folder_residuals <- args[4]
if (!dir.exists(folder_residuals)) dir.create(folder_residuals)

# set default
envs <- NULL
checks <- NULL

# assert to have the correct optional arguments
pos_args <- 4
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {

  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--envs", "--checks")
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

# make sure optional arguments are valid
if (!is.null(envs)) envs <- unlist(strsplit(envs, split = ","))
if (!is.null(checks)) checks <- unlist(strsplit(checks, split = ","))



#### format table ----

# read excel sheet with all trait information
data <- read_excel(pheno_file, sheet = "Combined", guess_max = 15000)

# format columns e remove unuseful info
data$Experiment <- gsub("NIFA-20", "", data$Experiment)
data$environment <- paste0(data$Location, data$Experiment)
data <- data[ c("Hybrid", "environment", "Replication", trait)]
colnames(data)[c(1, 3)] <- c("genotype", "rep")

# remove checks, if requested
if (!is.null(checks)) data <- data[!data$genotype %in% checks, ]

# make sure all genotypes have a rep number
data <- data %>%
  group_by(genotype, environment) %>%
  mutate(rep = 1:n())

# # transform to wide format
# data <- pivot_wider(data, names_from = "environment", values_from = all_of(trait))
#
# # keep only user-defined environments and make sure rep column is the last
# if (!is.null(envs)) {
#   data <- data[, c("genotype", envs, "rep")]
# } else {
#   # if no envs defined, use all of them
#   data <- data.frame(data[, c(1, 3:ncol(data))], rep = data$rep)
# }
#
# # reorder data
# data <- data[mixedorder(data$genotype), ]
#
# # transform df into long format again, but now NAs will be added for hybrids without info in an env
# data <- pivot_longer(data, !c(genotype, rep), names_to = "environment", values_to = "trait_value") %>%
#   mutate(rep = paste0("rep", rep)) %>%
#   relocate(rep, .after = environment) %>%
#   slice(mixedorder(genotype)) %>%
#   mutate(genotype = factor(genotype),
#          environment = factor(environment, levels = mixedsort(unique(environment))),
#          rep = factor(rep, levels = mixedsort(unique(rep))))



#### get BLUEs for each environment ----

# create empty data frame to store h2
h2_traits <- data.frame(stringsAsFactors = FALSE)

# get blues for each environment
for (env in envs) {

  print(env)

  # filter data by environment
  data_env <- subset(data, environment == env)
  # adjust data
  data_env$genotype <- factor(data_env$genotype)
  data_env$rep <- factor(data_env$rep)
  colnames(data_env)[4] <- "trait_value"
  data_env$trait_value <- as.numeric(data_env$trait_value)

  # run mixed models
  sim_trait_env <- asreml(fixed = trait_value ~ genotype,
                          random = ~ rep,
                          maxit = 100,
                          data = data_env,
                          trace = FALSE)
  # update model until converge and % change in parameters is lower than 0.1
  while (!sim_trait_env$converge | any(summary(sim_trait_env)$varcomp$`%ch` > 0.1, na.rm = TRUE)) {
    cat("updating model...\n")
    sim_trait_env <- update.asreml(sim_trait_env, trace = FALSE)
    cat("...done\n")
  }

  # inspect resitduals
  pdf(paste0(folder_residuals, "/residuals_", trait, ".", env, ".pdf"))
  plot(sim_trait_env)
  dev.off()

  # get predicted values
  predicted_values_env <- predict(sim_trait_env, classify = "genotype", sed = TRUE, maxit = 100)
  predicted_values_env <- data.frame(predicted_values_env$pvals[, c(1, 2)])

  # reorder table
  predicted_values_env$genotype <- as.character(predicted_values_env$genotype)
  predicted_values_env <- predicted_values_env[mixedorder(predicted_values_env$genotype), c("genotype", "predicted.value")]
  colnames(predicted_values_env) <- c("hybrid", "BLUE")

  # write file
  outfile <- paste0(outfolder, "/1stStage_BLUEs.", trait, "_", env, ".txt")
  fwrite(predicted_values_env, file = outfile, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

  # heritability
  h2_trait_env <- asreml(fixed = trait_value ~ 1,
                         random = ~ genotype + rep,
                         maxit = 100,
                         data = data_env,
                         trace = FALSE)
  # update model until converge and % change in parameters is lower than 0.1
  while (!h2_trait_env$converge | any(summary(h2_trait_env)$varcomp$`%ch` > 0.1, na.rm = TRUE)) {
    cat("updating model...\n")
    h2_trait_env <- update.asreml(h2_trait_env, trace = FALSE)
    cat("...done\n")
  }

  # calculate heritability
  h2 <- vpredict(h2_trait_env, h2 ~ V2 / (V1 + V2 + V3))
  h2 <- data.frame(trait = trait, env = env, h2 = round(as.numeric(h2$Estimate), digits = 3))

  # append h2 to final df
  h2_traits <- rbind(h2_traits, h2)

}

# write h2 file
outfile <- paste0(outfolder, "/1stStage_heritabilities.", trait, ".txt")
fwrite(h2_traits, file = outfile, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

#### debug ----

# pheno_file <- "data/NIFA_CompleteDataset.xlsx"
# # trait <- "YLD"
# outfolder <- "data"
# folder_residuals <- "analysis/pheno_BLUES_qc"
# # envs <- c("BEC-BL19", "BEC-BL20", "COR19", "COR20", "MIN19", "MIN20", "SYN19", "SYN20", "URB19")
# checks <- c("UIUC-CHECK1", "UIUC-CHECK2", "UIUC-CHECK3", "UIUC-CHECK4", "6049V2P")
# trait <- "EHT"
# envs <- c("COR19", "MIN19", "MIN20", "URB19")