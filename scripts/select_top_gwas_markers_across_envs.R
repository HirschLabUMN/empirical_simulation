library(data.table)
library(dplyr)
library(tidyr)

usage <- function() {
  cat("
description: filter gwas peaks according to the number of markers desired for correlation curve test.

usage: Rscript select_top_gwas_markers_across_envs.R [gwas_file] [ld_file] [n_markers_to_keep] [outfolder] [...]

positional arguments:
  gwas_file               file with non-redundant top GWAS hits summarized by environment (output of
                          script 'merge_gwas_per_trait.R')
  ld_file                 ld file among top markers in all environments
  n_markers_to_keep       number of top non-redudant markers across envs to keep
  outfolder               folder to save results

optional argument:
  --help                  show this helpful message
  --ld-threshold          LD threshold (default: 0.9)
  --avg-rank-markers      calculate average rank among 'all' markers or 'top' markers in each environment?
                          (default: 'all')


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

# get positional arguments
gwas_file <- args[1]
ld_file <- args[2]
n_markers_to_keep <- as.numeric(args[3])
outfolder <- args[4]

# set default for optional arguments
ld_threshold <- "0.9"
avg_rank_markers <- "all"

# assert to have the correct optional arguments
pos_args <- 4
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {

  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--ld-threshold", "--avg-rank-markers")
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
if (suppressWarnings(!is.na(as.numeric(ld_threshold)))) {
  ld_threshold <- as.numeric(ld_threshold)
} else {
  stop("Optional argument '--ld-threshold' should be a number")
}

if (!avg_rank_markers %in% c("top", "all")) {
  stop("Optional argument '--avg-rank-markers' should be 'all' or 'top'")
}



#### filter gwas top peak file ----

# load files
gwas_hits <- fread(gwas_file, header = TRUE, data.table = FALSE)
ld_markers <- fread(ld_file, header = TRUE, data.table = FALSE)

# create df to save results
gwas_hits_ld <- data.frame(stringsAsFactors = FALSE)

# need to rank additive and dom models separately
for (model in c("additive", "dominant")) {

  # subset by add or dom model
  gwas_hits_model <- subset(gwas_hits, type == model)

  # create df to save results
  gwas_peaks_model <- data.frame(stringsAsFactors = FALSE)
  # create list to keep non-redudant markers -- this list will start with all
  # markers before the loop, and after each iteration, the list will be updated
  # and markers belonging to the same GWAS peak will be removed -- at the end of
  # the loop, only non-redudant markers should remain
  non_redundant_markers_list <- gwas_hits_model$marker

  for (marker in gwas_hits_model$marker) {

    # if (marker == gwas_hits_model$marker[1]) stop("debug mode: run line by line")
    # marker <- "snp.1.256386607"

    if (marker %in% non_redundant_markers_list) {

      # make sure only markers on non-redudant list stay here
      ld_marker <- subset(ld_markers, SNP_A %in% non_redundant_markers_list & SNP_B %in% non_redundant_markers_list)
      # filter ld table to have that sig hit
      ld_marker <- subset(ld_marker, SNP_A == marker | SNP_B == marker)
      # keep markers that belong to the same GWAS peak
      ld_marker <- subset(ld_marker, R2 >= ld_threshold)

      # if marker is not in ld with anything else...
      if (nrow(ld_marker) == 0) {

        # leave it by itself
        gwas_hits_model[gwas_hits_model$marker == marker, "marker"] <- marker

        # remove redundant markers from list
        non_redundant_markers_list <- non_redundant_markers_list[!non_redundant_markers_list == marker]

      } else {

        # otherwise...
        # gather info for each marker in the peak
        markers_in_peak <- unique(c(ld_marker$SNP_A, ld_marker$SNP_B))
        markers_in_peak <- gwas_hits_model[match(markers_in_peak, gwas_hits_model$marker), ]

        # get marker with highest effect - use abs() function to get rid of effect's sign
        marker_to_keep <- markers_in_peak[abs(markers_in_peak$effect) == max(abs(markers_in_peak$effect)), ]
        # if more than one marker with same effect...
        if (nrow(marker_to_keep) > 1) {
          # keep the with lowest pvalue
          marker_to_keep <- marker_to_keep[marker_to_keep$pval == min(marker_to_keep$pval), ]
          # if still more than one marker...
          if (nrow(marker_to_keep) > 1) {
            # keep the with highest maf
            marker_to_keep <- marker_to_keep[marker_to_keep$maf == max(marker_to_keep$maf), ]
            # if still more than one marker...
            if (nrow(marker_to_keep) > 1) {
              # get first marker on the list
              marker_to_keep <- marker_to_keep[1, ]
            }
          }
        }

        # add marker info to results table
        for (redundant_marker in markers_in_peak$marker) {
          gwas_hits_model[gwas_hits_model$marker == redundant_marker, "marker"] <- marker
        }

        # remove redundant markers from list
        non_redundant_markers_list <- non_redundant_markers_list[!non_redundant_markers_list %in% markers_in_peak$marker]

      }
    }
  }

  # append to main df
  gwas_hits_ld <- rbind(gwas_hits_ld, gwas_hits_model)

}
# clean environment
rm(gwas_hits, gwas_hits_model, gwas_peaks_model, ld_marker, ld_markers, marker_to_keep, markers_in_peak)

# if markers within same environment are in ld, keep only the one with lowest pval
gwas_hits_ld <- gwas_hits_ld %>%
  group_by(env, marker, type) %>%
  arrange(pval, .by_group = TRUE) %>%
  filter(row_number() == 1)

# create empty df to save results
rank_markers <- data.frame(stringsAsFactors = FALSE)
gwas_non_redudant <- data.frame(stringsAsFactors = FALSE)

# select the top markers based on highest average rank across envs
for (model in c("additive", "dominant")) {

  if (avg_rank_markers == "top") {

    # get average rank of top markers markers across environments
    gwas_hits_ld_rank <- gwas_hits_ld %>%
      # get marker rank
      filter(type == model) %>%
      group_by(env) %>%
      arrange(pval, .by_group = TRUE) %>%
      filter(row_number() <= n_markers_to_keep) %>%
      mutate(rank = rank(pval)) %>%
      # if there's linked markers in the same env, get avg rank
      group_by(env, marker) %>%
      summarize(type = unique(type), rank = mean(rank)) %>%
      ungroup() %>%
      # transform to wide format so each marker has a value in each env
      pivot_wider(names_from = "env", values_from = "rank") %>%
      # replace NAs with lowest rank possible
      replace(is.na(.), values = (n_markers_to_keep + 1)) %>%
      # get average rank across envs
      rowwise() %>%
      mutate(avg_rank = mean(c(COR19, MIN19, MIN20, URB19), na.rm = TRUE)) %>%
      ungroup() %>%
      # keep only top markers
      arrange(avg_rank) %>%
      filter(row_number() <= n_markers_to_keep)

  }

  if (avg_rank_markers == "all") {

    # get avg rank of all markers across all envs
    gwas_hits_ld_rank <- gwas_hits_ld %>%
      # get marker rank
      filter(type == model) %>%
      group_by(env) %>%
      arrange(pval, .by_group = TRUE) %>%
      # filter(row_number() <= n_markers_to_keep) %>%
      mutate(rank = rank(pval)) %>%
      # if there's linked markers in the same env, get avg rank
      group_by(env, marker) %>%
      summarize(type = unique(type), rank = mean(rank)) %>%
      ungroup() %>%
      # transform to wide format so each marker has a value in each env
      pivot_wider(names_from = "env", values_from = "rank")
    # replace NAs with lowest rank possible
    gwas_hits_ld_rank <- replace(gwas_hits_ld_rank, is.na(gwas_hits_ld_rank), values = nrow(gwas_hits_ld_rank) + 1)
    # get average rank across envs
    gwas_hits_ld_rank$avg_rank <- rowMeans(gwas_hits_ld_rank[, 3:ncol(gwas_hits_ld_rank)], na.rm = TRUE)
    # keep only top markers
    gwas_hits_ld_rank <- gwas_hits_ld_rank[order(gwas_hits_ld_rank$avg_rank), ]
    gwas_hits_ld_rank <- head(gwas_hits_ld_rank, n = n_markers_to_keep)

  }

  # keep gwas file only with top markers
  gwas_hits_ld_filtered <- subset(gwas_hits_ld, type == model & marker %in% gwas_hits_ld_rank$marker)
  gwas_hits_ld_filtered <- gwas_hits_ld_filtered[order(gwas_hits_ld_filtered$signif, gwas_hits_ld_filtered$chr, gwas_hits_ld_filtered$pos,
                                                       decreasing = c(TRUE, FALSE, FALSE)), ]

  # append to main df
  rank_markers <- rbind(rank_markers, as.data.frame(gwas_hits_ld_rank))
  gwas_non_redudant <- rbind(gwas_non_redudant, as.data.frame(gwas_hits_ld_filtered))

}
rm(gwas_hits_ld, gwas_hits_ld_rank, gwas_hits_ld_filtered)

# create folders
outfolder_ranks <- paste0(outfolder, "/rank_markers")
if (!dir.exists(outfolder_ranks)) dir.create(outfolder_ranks, recursive = TRUE)
# save marker ranks
outfile_ranks <- paste0(outfolder_ranks, "/top", n_markers_to_keep, "_add-dom_markers.avg-rank-", avg_rank_markers, ".txt")
fwrite(rank_markers, file = outfile_ranks, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# create folders
outfolder_gwas <- paste0(outfolder, "/gwas_markers")
if (!dir.exists(outfolder_gwas)) dir.create(outfolder_gwas, recursive = TRUE)
# save marker ranks
outfile_gwas <- paste0(outfolder_gwas, "/gwas_top", n_markers_to_keep, "-peaks.avg-rank-", avg_rank_markers, ".txt")
fwrite(gwas_non_redudant, file = outfile_gwas, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)



#### debug ----

# gwas_file <- "analysis/gwas/gwas_top200-peaks_EHT-per-env.txt"
# ld_file <- "analysis/sim_traits_gwas/cor_curve/EHT/ld/ld_top-markers.ld.gz"
# n_markers_to_keep <- 10
# outfolder <- "analysis/sim_traits_gwas/cor_curve/EHT"
# ld_threshold <- 0.9
# avg_rank_markers <- "top"
# # avg_rank_markers <- "all"
