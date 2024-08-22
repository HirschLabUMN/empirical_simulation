library(data.table)

usage <- function() {
  cat("
description: filter number of environments from the 1st stage of genomic prediction.

usage: Rscript filter_envs_1st_stage.R [blues_file] [envs_to_keep] [...]

positional arguments:
  blues_file          file with BLUEs in long format
  envs_to_keep        comma-separated list of environment names to keep
  
optional argument:
  --help              show this helpful message

"
  )
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 2) stop(usage(), "missing positional argument(s)")

# get positional arguments
blues_file <- args[1]
envs_to_keep <- unlist(strsplit(args[2], split = ","))



#### filter blues file ----

# load data
blues <- fread(blues_file, header = TRUE, data.table = FALSE)

# filter data
blues_filtered <- subset(blues, environment %in% envs_to_keep)

# overwrite file
fwrite(blues_filtered, file = blues_file, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)



#### debug ---- 

# blues_file <- "analysis/sim_traits_gwas/YLD/blues_1st_stage.txt"
# envs_to_keep <- "env4,env5,env6,env7,env8,env9"
# envs_to_keep <- unlist(strsplit(envs_to_keep, split = ","))
