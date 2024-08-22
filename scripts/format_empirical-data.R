library(data.table)
library(readxl)
library(gtools)
library(dplyr)
library(tidyr)

usage <- function() {
  cat("
description: format empirical phenotypic from an excel sheet to a tab-delimited file
             in the wide format (first column are genotypes, remaining columns are
             environments).

usage: Rscript format_empirical-data_for_predictions.R [excel_file] [...]

positional arguments:
  excel_file          excel file with phenotypic data, one trait per sheet
  trait               name of trait to extract phenotypic values
  
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
excel_file <- args[1]
trait <- args[2]

# set default
envs <- NULL
checks <- NULL

# assert to have the correct optional arguments
pos_args <- 2
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
data <- read_excel(excel_file, sheet = "Combined", guess_max = 15000)

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

# transform to wide format
data <- pivot_wider(data, names_from = "environment", values_from = all_of(trait))

# keep only user-defined environments and make sure rep column is the last
if (!is.null(envs)) {
  data <- data[, c("genotype", envs, "rep")]
} else {
  # if no envs defined, use all of them
  data <- data.frame(data[, c(1, 3:ncol(data))], rep = data$rep)
}

# reorder data
data <- data[mixedorder(data$genotype), ]

# write file
outfile <- gsub(".xlsx", paste0(".", trait, ".txt"), excel_file)
fwrite(data, file = outfile, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)



#### debug ----

# excel_file <- "data/NIFA_CompleteDataset.xlsx"
# trait <- "YLD"
# envs <- c("BEC-BL19", "BEC-BL20", "BEC-EP20", "COR19", "COR20", "MIN19", "MIN20", "SYN19", "SYN20", "URB19")
# checks <- c("UIUC-CHECK1", "UIUC-CHECK2", "UIUC-CHECK3", "UIUC-CHECK4", "6049V2P")
