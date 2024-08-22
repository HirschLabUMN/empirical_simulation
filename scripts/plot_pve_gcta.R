library(data.table)
library(ggplot2)

usage <- function() {
  cat("
description: plot percent phenotypic variance explained by all markers and
             significant markers from GWAS.

usage: Rscript plot_pve_gcta.R [gcta_folder] [...]

positional arguments:
  gcta_folder         folder with GCTA results where each trait is a subfolder


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
gcta_folder <- args[1]
# gcta_folder <- "analysis/gcta"


# create empty dataframe to store results
summary_gcta <- data.frame(stringsAsFactors = FALSE)

# get traits used in gcta
traits <- list.dirs(gcta_folder, full.names = FALSE, recursive = FALSE)

for (trait in traits) {

  # get trait folder
  trait_folder <- paste0(gcta_folder, "/", trait)

  # get gcta files for all markers and for significant gwas hits
  gcta_files <- list.files(trait_folder, pattern = ".hsq", full.names = TRUE)
  all_markers_files <- gcta_files[grep("pve_all-markers", gcta_files)]
  top_markers_files <- gcta_files[grep("pve_top_non-redudant_markers", gcta_files)]

  for (model in c("add", "dom")) {

    # get pve results for all markers
    all_markers <- all_markers_files[grep(model, all_markers_files)]
    all_markers <- fread(all_markers, header = TRUE, fill = TRUE, data.table = FALSE)
    all_markers <- data.frame(trait = trait, model = model, markers = "all",
                              pve = all_markers[all_markers$Source == "V(G)/Vp", "Variance"],
                              se = all_markers[all_markers$Source == "V(G)/Vp", "SE"],
                              pval = all_markers[all_markers$Source == "Pval", "Variance"])

    # get pve results for significant markers
    top_markers <- top_markers_files[grep(model, top_markers_files)]
    if (length(top_markers) > 0) {

      top_markers <- fread(top_markers, header = TRUE, fill = TRUE, data.table = FALSE)
      top_markers <- data.frame(trait = trait, model = model, markers = "top",
                             pve = top_markers[top_markers$Source == "V(G)/Vp", "Variance"],
                             se = top_markers[top_markers$Source == "V(G)/Vp", "SE"],
                             pval = top_markers[top_markers$Source == "Pval", "Variance"])

    } else {

      top_markers <- data.frame(trait = trait, model = model, markers = "top",
                             pve = NA, se = NA, pval = NA)

    }

    # append to final summary table
    summary_gcta <- rbind(summary_gcta, rbind(all_markers, top_markers))

  }
}

# write summary
outfile_summary <- paste0(gcta_folder, "/summary_pve_gcta.txt")
fwrite(summary_gcta, file = outfile_summary, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# transform columns to factor for plotting
summary_gcta$model <- factor(summary_gcta$model, labels = c("additive", "dominant"))
summary_gcta$markers <- factor(summary_gcta$markers, levels = c("all", "top"))

# plot summary
summary_plot <- ggplot(summary_gcta) +
  geom_errorbar(aes(x = trait, ymin = pve - 0.05, ymax = pve + se, color = markers),
                width = 0.4, show.legend = FALSE, position = position_dodge(0.9)) +
  geom_bar(aes(x = trait, y = pve, fill = markers), stat = "identity", position = "dodge") +
  facet_wrap(~ model, nrow = 2) +
  coord_cartesian(y = c(0, 1)) +
  scale_y_continuous(labels = function(x) x * 100, expand = c(0,0)) +
  labs(x = "Trait/env", y = "Percent variance explained (%)", fill = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "top")
ggsave(summary_plot, filename = paste0(gcta_folder, "/pve_all-vs-top_markers.pdf"),
       height = 6, width = 14)
