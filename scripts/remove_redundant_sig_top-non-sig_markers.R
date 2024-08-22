library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

usage <- function() {
  cat("
description: remove redudant significant and top non-significant GWAS hits due to LD.

usage: Rscript remove_redundant_sig_top-non-sig_markers.R [gwas_all_file] [ld_top_hits_file] [output_folder] [...]

positional arguments:
  gwas_all_file                   full GWAS results table from GAPIT
  ld_top_hits_file                plink LD file from significant and top non-significant hits
  output_folder                   folder to save files

optional argument:
  --help                          show this helpful message
  --n-top-non-sig=VALUE           number of top non-significant hits to keep (default: 1000)
  --pval-threshold=VALUE          p-value threshold used for significance in GWAS (default: 0.05)
  --ld-threshold=VALUE            R2 value threshold to determine which markers are in LD (default: 0.9)


credits:
  plots inspired by Yan Holtz's tutorial at https://www.r-graph-gallery.com/101_Manhattan_plot.html
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
gwas_all_file <- args[1]
ld_top_hits_file <- args[2]
output_folder <- args[3]

# set default for optional arguments
n_top_non_sig <- "1000"
pval_threshold <- "0.05"
ld_threshold <- "0.9"

# assert to have the correct optional arguments
pos_args <- 3
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {

  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--n-top-non-sig", "--pval-threshold", "--ld-threshold")
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
if (suppressWarnings(!is.na(as.numeric(n_top_non_sig)))) {
  n_top_non_sig <- as.numeric(n_top_non_sig)
} else {
  stop("Optional argument '--n-top-non-sig' should be a number")
}

if (suppressWarnings(!is.na(as.numeric(pval_threshold)))) {
  pval_threshold <- as.numeric(pval_threshold)
} else {
  stop("Optional argument '--pval-threshold' should be a number")
}

if (suppressWarnings(!is.na(as.numeric(ld_threshold)))) {
  ld_threshold <- as.numeric(ld_threshold)
} else {
  stop("Optional argument '--ld-threshold' should be a number")
}



#### remove redundant markers ----

# load data
gwas_all <- fread(gwas_all_file, header = TRUE, data.table = FALSE)
# get significant hits
sig_hits <- gwas_all[gwas_all$P.value < pval_threshold, ]
# adjust df depending on number of significant markers
if (nrow(sig_hits) > 0) {
  sig_hits <- sig_hits[order(sig_hits$P.value), ]
  sig_hits$sig <- TRUE
} else {
  sig_hits <- data.frame(matrix(nrow = 0, ncol = 11))
  colnames(sig_hits) <- c(colnames(gwas_all), "sig")
}
# get top non-significant hits
non_sig_hits <- gwas_all[order(gwas_all$P.value)[(nrow(sig_hits) + 1):(n_top_non_sig + nrow(sig_hits))], ]
non_sig_hits$sig <- FALSE
# merge sig and non sig hits
gwas_top <- rbind(sig_hits, non_sig_hits)
rm(gwas_all, sig_hits, non_sig_hits)

# load ld results
ld_top_hits <- fread(ld_top_hits_file, header = TRUE, data.table = FALSE)

# create df to save results
gwas_peaks <- data.frame(stringsAsFactors = FALSE)
# create list to keep non-redudant markers -- this list will start with all
# markers before the loop, and after each iteration, the list will be updated
# and markers belonging to the same GWAS peak will be removed -- at the end of
# the loop, only non-redudant markers should remain
non_redundant_markers_list <- gwas_top$SNP

for (marker in gwas_top$SNP) {

  # if (marker == "snp.1.268984334") stop("debug mode: run line by line")

  if (marker %in% non_redundant_markers_list) {

    # make sure only markers on non-redudant list stay here
    ld_marker <- subset(ld_top_hits, SNP_A %in% non_redundant_markers_list & SNP_B %in% non_redundant_markers_list)
    # filter ld table to have that sig hit
    ld_marker <- subset(ld_marker, SNP_A == marker | SNP_B == marker)
    # keep markers that belong to the same GWAS peak
    ld_marker <- subset(ld_marker, R2 >= ld_threshold)

    # if marker is not in ld with anything else...
    if (nrow(ld_marker) == 0) {

      # leave it by itself
      gwas_peaks <- rbind(gwas_peaks, data.frame(marker = gwas_top[gwas_top$SNP == marker, "SNP"],
                                                 chr = gwas_top[gwas_top$SNP == marker, "Chromosome"],
                                                 pos = gwas_top[gwas_top$SNP == marker, "Position"],
                                                 maf = gwas_top[gwas_top$SNP == marker, "maf"],
                                                 effect = gwas_top[gwas_top$SNP == marker, "effect"],
                                                 pval = gwas_top[gwas_top$SNP == marker, "P.value"],
                                                 signif = gwas_top[gwas_top$SNP == marker, "sig"]))

      # remove redundant markers from list
      non_redundant_markers_list <- non_redundant_markers_list[!non_redundant_markers_list == marker]

    } else {

      # otherwise...
      # gather info for each marker in the peak
      markers_in_peak <- unique(c(ld_marker$SNP_A, ld_marker$SNP_B))
      markers_in_peak <- gwas_top[match(markers_in_peak, gwas_top$SNP), ]

      # get marker with highest effect - use abs() function to get rid of effect's sign
      marker_to_keep <- markers_in_peak[abs(markers_in_peak$effect) == max(abs(markers_in_peak$effect)), ]
      # if more than one marker with same effect...
      if (nrow(marker_to_keep) > 1) {
        # keep the with lowest pvalue
        marker_to_keep <- marker_to_keep[marker_to_keep$P.value == min(marker_to_keep$P.value), ]
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
      gwas_peaks <- rbind(gwas_peaks, data.frame(marker = marker_to_keep$SNP,
                                                 chr = marker_to_keep$Chromosome,
                                                 pos = marker_to_keep$Position,
                                                 maf = marker_to_keep$maf,
                                                 effect = marker_to_keep$effect,
                                                 pval = marker_to_keep$P.value,
                                                 signif = marker_to_keep$sig))

      # remove redundant markers from list
      non_redundant_markers_list <- non_redundant_markers_list[!non_redundant_markers_list %in% markers_in_peak$SNP]

    }
  }
}

# remove NAs if generated
gwas_peaks <- gwas_peaks[!is.na(gwas_peaks$pval), ]
# sort table by chr and position
gwas_peaks <- gwas_peaks[order(gwas_peaks$signif, gwas_peaks$chr, gwas_peaks$pos, decreasing = c(TRUE, FALSE, FALSE)), ]
# add peak numbers
if (sum(gwas_peaks$signif == TRUE) > 0) {
  gwas_peaks[gwas_peaks$signif == TRUE, "peak"] <- paste0(1:sum(gwas_peaks$signif == TRUE), "_sig")
}
gwas_peaks[gwas_peaks$signif == FALSE, "peak"] <- paste0(1:sum(gwas_peaks$signif == FALSE), "_nonsig")

# write file
outfile_peaks <- paste0(output_folder, "/top-gwas-peaks_non-redundant_markers.txt")
fwrite(gwas_peaks, file = outfile_peaks, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)




#### create manhattan plot ----

# read gwas results again
gwas_all <- fread(gwas_all_file, header = TRUE, data.table = FALSE)

# add cumulative chr info
gwas_all <- gwas_all %>%
  # compute chromosome size
  group_by(Chromosome) %>%
  summarize(chr_len = max(Position)) %>%
  # calculate cumulative position of each chromosome
  mutate(chr_cum = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  # add this info to the initial dataset
  left_join(gwas_all, ., by = c("Chromosome" = "Chromosome")) %>%
  # add a cumulative position of each SNP,
  arrange(Chromosome, Position) %>%
  mutate(pos_cum = Position + chr_cum)

# prepare x axis to not want to display the cumulative position of SNP in bp
# (just show the chromosome name instead)
axis_plot <- gwas_all %>%
  group_by(Chromosome) %>%
  summarize(center = (as.numeric(max(pos_cum)) + as.numeric(min(pos_cum))) / 2)

# keep only columns that matter
gwas_all <- gwas_all[, c("SNP", "chr_cum", "pos_cum")]
# merge gwas info with peaks
gwas_plot_info <- merge(x = gwas_peaks, y = gwas_all, by.x = "marker", by.y = "SNP")
# add log10 pval column
gwas_plot_info$log_pval <- -log10(gwas_plot_info$pval)

# basic manhattan plot
gwas_plot <- gwas_plot_info %>%
  # # uncomment this if need to speed up plotting (skip markers with very high pvalue)
  # filter(-log10(P.value) > 1) %>%
  ggplot(aes(x = pos_cum, y = log_pval)) +
  # show all points
  geom_point(aes(color = as.factor(chr)), alpha = 0.8, size = 3) +
  # scale_color_manual(values = rep(c("grey", "skyblue"), length(unique(gwas_plot_info$chr)))) +
  # add pvalue threshold
  geom_hline(yintercept = -log10(pval_threshold), color = "firebrick") +
  # # add highlighted points
  # geom_point(data = subset(gwas_plot_info, is_sv == "yes"), color = "orange", size = 2) +
  # add label using ggrepel to avoid overlapping
  geom_label_repel(data = subset(gwas_plot_info, signif == TRUE), aes(label = marker), size = 2) +
  # custom X axis
  scale_x_continuous(label = axis_plot$Chromosome, breaks = axis_plot$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ceiling(max(gwas_plot_info$log_pval, na.rm = TRUE)))) +
  # add labels
  labs(subtitle = "Non-redundant markers from top GWAS peaks",
       x = "Chromosome",
       y = "-log10(p-value)") +
  # custom theme
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        text = element_text(size = 15))

# save plot
outfile_plot <- paste0(output_folder, "/top-gwas-peaks_non-redundant_markers.pdf")
# write file
ggsave(filename = outfile_plot, plot = gwas_plot, device = "pdf", width = 16, height = 8)



#### debug ----

# gwas_all_file <- "analysis/gwas/YLD_MIN20/additive_model/permutation/GAPIT.MLM.BLUE.GWAS.Results.csv"
# ld_top_hits_file <- "analysis/gwas/YLD_MIN20/additive_model/permutation/ld_top5000_non-sig-markers.ld.gz"
# output_folder <- "analysis/gwas/YLD_MIN20/additive_model/permutation"
# n_top_non_sig <- 5000
# pval_threshold <- 5.4810743568268e-06
# ld_threshold <- 0.9
