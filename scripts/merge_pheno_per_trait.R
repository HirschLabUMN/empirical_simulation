library(data.table)
library(gtools)
library(tidyr)

usage <- function() {
  cat("
description: merge phenotypic data for each environment of a single trait.

usage: Rscript merge_pheno_per_trait.R [trait] [envs] real_pheno_folder] [real_pheno_prefix] [...]

positional arguments:
  trait                 name of trait to merge data from each environment
  envs                  comma-separated list of all environments
  real_pheno_folder     folder with real phenotypic data
  real_pheno_prefix     prefix of filename in common for all real phenotypic data
  
optional argument:
  --help              show this helpful message

"
  )
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 4) stop(usage(), "missing positional argument(s)")

# assign arguments to variables
trait <- args[1]
envs <- unlist(strsplit(args[2], split = ","))
real_pheno_folder <- args[3]
real_pheno_prefix <- args[4]
# trait <- "EHT"
# envs <- c("BAY19", "BEC-BL19", "BEC-BL20", "BEC-EP20", "COR19", "COR20", "MIN19", "MIN20", "SYN19", "SYN20", "URB19")
# real_pheno_folder <- "data"
# real_pheno_prefix <- "1stStage_BLUEs"


#### merge phenotypic data ----

# get filenames of real phenotypic values
real_pheno_files <- list.files(path = real_pheno_folder,
                               pattern = paste0(real_pheno_prefix, ".", trait, "_"),
                               full.names = TRUE)

# create empty df to store real phenotypic data
real_pheno <- data.frame(stringsAsFactors = FALSE)

for (env in envs) {
  
  # look for trait-env file
  real_pheno_env_file <- real_pheno_files[grep(env, real_pheno_files)]
  
  if (length(real_pheno_env_file) > 0) {
    
    # load real phenotypic values
    real_pheno_env <- fread(real_pheno_env_file, header = TRUE, data.table = FALSE)
    
    if (nrow(real_pheno) == 0) {
      # add both hybrid and pheno columns
      real_pheno <- rbind(real_pheno, real_pheno_env)
    } else {
      # merge tables based on hybrid id
      real_pheno <- merge(x = real_pheno, y = real_pheno_env, by = "hybrid", all = TRUE)
    }
    
    # change column name to env name
    colnames(real_pheno)[ncol(real_pheno)] <- env
    # order rows
    real_pheno <- real_pheno[mixedorder(real_pheno$hybrid), ]
    rm(real_pheno_env)
    
  }
}

# transform to long format
real_pheno <- pivot_longer(real_pheno, cols = -hybrid, names_to = "env", values_to = "real_pheno")
real_pheno <- real_pheno[mixedorder(real_pheno$env), ]

# write merged table
fwrite(real_pheno, file = paste0(real_pheno_folder, "/", real_pheno_prefix, ".", trait, "-per-env.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
