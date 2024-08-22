library(data.table)
library(ggplot2)

usage <- function() {
  cat("
description: plot LD distibution of significant GWAS hits.

usage: Rscript plot_ld_sig-markers.R [gwas_folder] [...]

positional arguments:
  gwas_folder         folder to with gwas results with the internal structure as
                      'trait_folder/model_folder/'

optional argument:
  --help              show this helpful message

"
  )
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 1) stop(usage(), "missing positional argument(s)")

# assign arguments to variables
gwas_folder <- args[1]
# gwas_folder <- "analysis/gwas"


#### summarize gwas ----

# create empty dataframe to store results
ld_summary <- data.frame(stringsAsFactors = FALSE)

# get traits used in gwas
traits <- list.dirs(gwas_folder, full.names = FALSE, recursive = FALSE)

for (trait in traits) {
  
  # get trait folder
  trait_folder <- paste0(gwas_folder, "/", trait)
  
  for (model in c("additive", "dominant")) {
    
    for (pval_method in c("FDR", "permutation")) {
      
      # get model folder
      model_folder <- list.dirs(trait_folder, recursive = TRUE)
      model_folder <- model_folder[grep(model, model_folder)]
      model_folder <- model_folder[grep(pval_method, model_folder)]
      
      # open file with significant hits
      ld_sig_hits <- list.files(model_folder, full.names = TRUE)
      ld_sig_hits <- ld_sig_hits[grep("ld_sig-markers.ld", ld_sig_hits)]
      
      # get significant hit info, if any
      if (length(ld_sig_hits) > 0) {
        
        ld_sig_hits <- fread(ld_sig_hits, header = TRUE, data.table = FALSE)
        ld_sig_hits <- cbind(trait = trait, model = model, pval_method = pval_method, ld_sig_hits)
        
      }
    }
    
    # append to final summary table
    ld_summary <- rbind(ld_summary, ld_sig_hits)
    
  }
}

# write summary
outfile_summary <- paste0(gwas_folder, "/summary_ld_sig-hits.txt")
fwrite(ld_summary, file = outfile_summary, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# check which ld values are from different chromosomes
extra_ld_info <- apply(ld_summary, MARGIN = 1, function(marker) {
  if (marker["CHR_A"] == marker["CHR_B"]) {
    inter_chr <- "no"
    marker_dist <- abs(as.numeric(marker["BP_A"]) - as.numeric(marker["BP_B"]))
  } else {
    inter_chr <- "yes"
    marker_dist <- NA
  }
  return(c(inter_chr = inter_chr, marker_dist = marker_dist))
})
extra_ld_info <- data.frame(t(extra_ld_info))

# append extra ld info to main df
ld_summary <- cbind(ld_summary, extra_ld_info)
rm(extra_ld_info)

# plot distribution
ld_summary$inter_chr <- factor(ld_summary$inter_chr, levels = c("no", "yes"),
                               labels = c("same chromosomes", "different chromosomes"))
ld_inter_chr <- ggplot(ld_summary) +
  facet_wrap(~ inter_chr) +
  geom_histogram(aes(x = R2))
ggsave(ld_inter_chr, filename = paste0(gwas_folder, "/dist_ld_sig-markers.pdf"),
       height = 6, width = 14)

ld_same_chr <- ggplot(subset(ld_summary, inter_chr == "same chromosomes")) +
  geom_point(aes(x = as.numeric(marker_dist), y = R2)) +
  labs(x = "Marker distance (bp)", y = "R2")
ggsave(ld_same_chr, filename = paste0(gwas_folder, "/dist_ld-bp_sig-markers.pdf"),
       height = 6, width = 14)



