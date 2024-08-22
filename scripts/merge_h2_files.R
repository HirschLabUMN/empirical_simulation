library(data.table)
library(tidyr)

usage <- function() {
  cat("
description: merge heritability files from different traits.

usage: Rscript merge_h2_files.R [h2_folder] [h2_prefix] [...]

positional arguments:
  h2_folder          folder with heritability files
  h2_prefix          prefix of h2 filename
  
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
h2_folder <- args[1]
h2_prefix <- args[2]
# h2_folder <- "data"
# h2_prefix <- "1stStage_heritabilities"



#### merge h2 files ----

# get list of h2 files
h2_files <- list.files(path = h2_folder, pattern = h2_prefix, full.names = TRUE)

# merge files
h2 <- data.frame(stringsAsFactors = FALSE)
for (h2_file in h2_files) {
  h2 <- rbind(h2, fread(h2_file, header = TRUE, data.table = FALSE))
}

# transform to wide format
h2 <- pivot_wider(h2, names_from = "trait", values_from = "h2")

# write file
outfile <- paste0(h2_folder, "/", h2_prefix, ".csv")
fwrite(h2, file = outfile, quote = FALSE, sep = ",", na = NA, row.names = FALSE)
