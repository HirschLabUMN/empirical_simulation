library(data.table)
#####################
# Loading Libraries #
#####################
library("multtest")
library("snpStats")
library("gplots")
library("LDheatmap")
library("genetics")
library("ape")
library("EMMREML")
library("compiler")
library("scatterplot3d")
library("bigmemory")
library("biganalytics")

##############################
# Sourcing GAPIT and FarmCPU #
##############################
source("/home/hirschc1/burns756/empirical_sim/GAPIT_Source_Code/GAPIT_Source_Code.txt")

usage <- function() {
  cat("
description: transform a marker dataset from hapmap to numeric format (.num.txt),
             and generates an additional file with marker information (.map.txt).
             Files are save in same folder as input file.

usage: Rscript hmp2numeric.R [hmp_file] [...]

positional arguments:
  hmp_file            hapmap file with genotypic data
  
optional argument:
  --help              show this helpful message
  --model=VALUE       choose to numericalize using either 'additive' (default)
                      or 'dominant' model. The former transforms homozygous
                      alleles into 0 or 2, and heterozygous alleles to 1. The
                      latter transforms homozygous alleles to 0 and hets to 1.

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
model <- "additive"

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 1) stop(usage(), "missing positional argument(s)")

if (length(args) > 1) {

  opt_args <- args[-1]
  opt_args_allowed <- c("--model")
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
if (!model %in% c("additive", "dominant")) {
  stop("Optional argument '--model' should be either 'additive' or 'dominant'")
}

# get positional arguments
hmp_file <- args[1]



#### transform hapmap to numeric ----

# load data
hmp <- fread(hmp_file, header = FALSE, data.table = FALSE)

# adjust numericalization parameters according to model
if (model == "additive") {
  marker_effect = "Add"
  marker_impute = "Major"
}

if (model == "dominant") {
  marker_effect = "Dom"
  marker_impute = "Middle"
}

# transform to numeric
hmp_num <- GAPIT.HapMap(G = hmp, SNP.effect = marker_effect, SNP.impute = marker_impute)
rm(hmp)
# format numeric marker data to be compatible to GAPIT
markers <- data.frame(hmp_num$GD)
colnames(markers) <- hmp_num$GI$SNP
markers <- cbind(taxa = hmp_num$GT, markers)
# generate genetic map with marker names, chr and position
map <- hmp_num$GI
rm(hmp_num)

# name output files
outfile_num <- gsub(".hmp.txt", paste0(".", model, ".num.txt"), hmp_file)
outfile_map <- gsub(".hmp.txt", paste0(".", model, ".map.txt"), hmp_file)
# write files
fwrite(x = markers, file = outfile_num, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
fwrite(x = map, file = outfile_map, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)



#### debug ----

# hmp_file <- "tests/usda_hybrids.all_markers.adjusted-n-markers.hmp.txt"
# model <- "additive"
# # model <- "dominant"
