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
description: calculate the kinship matrix among markers using the VanRaden method.

usage: Rscript kinship_VanRaden.R [geno_file] [output_folder] [...]

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



#### calculate kinship ----

# load data
geno_data <- fread(geno_file, header = TRUE, data.table = FALSE)

# format genotypic information
rownames(geno_data) <- geno_data[, 1]
geno_data <- geno_data[, -1]
# calculate kinship with the VanRaden method
kinship_matrix <- GAPIT.kinship.VanRaden(geno_data)
# format table
kinship_matrix <- cbind(taxa = rownames(kinship_matrix), kinship_matrix)
rownames(kinship_matrix) <- 1:NROW(kinship_matrix)

# write kinship matrix
fwrite(kinship_matrix, file = paste0(output_folder, "/GAPIT.Kin.VanRaden.csv"),
       sep = ",", quote = FALSE, na = NA, row.names = FALSE, col.names = FALSE)



#### debug ----

# geno_file <- "tests/usda_hybrids.all_markers.adjusted-n-markers.additive.num.txt"
# output_folder <- "tests"
