library(data.table)
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
source("/home/hirschc1/burns756/empirical_sim/GAPIT_Source_Code/GAPIT_Source_Code.txt")

usage <- function() {
  cat("
description: run PCA on genotypic data and generate diagnostic plots to help
             decide how many principal components should be used in GWAS.

usage: Rscript pca_prior_gwas.R [geno_file] [output_folder] [...]

positional arguments:
  geno_file           file with numeric genotypic data
  output_folder       folder to save files
  
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
geno_file <- args[1]
output_folder <- args[2]



#### run PCA ----

# load data
geno <- fread(geno_file, head = TRUE, data.table = FALSE)

# working directory
init_folder <- getwd()
# change directories for plotting
setwd(paste0(init_folder, "/", output_folder))
# plot pca files
pca <- GAPIT.PCA(geno[, -1], taxa = geno[, 1], PCA.total = 10)
# return to initial directory
setwd(init_folder)



#### debug ----

# geno_file <- "tests/usda_hybrids.all_markers.adjusted-n-markers.additive.num.txt"
# output_folder <- "tests"
