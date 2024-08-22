library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(shades)
library(gtools)

print('Libraries Loaded')

usage <- function() {
  cat("
description: summarize GWAS results.

usage: Rscript summarize_gwas_results.R [gwas_folder] [chr_info] [centromere_info] [...]

positional arguments:
  gwas_folder         folder to with gwas results with the internal structure as
                      'trait_folder/model_folder/'
  chr_info            tab-delimited file with chromosome coordinates (format: chr, length)
  centromere_info     tab-delimited file with centromere coordinates (format: chr, start_pos, end_pos)

optional argument:
  --help              show this helpful message

"
  )
}

print('Reading in CLI options')

#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 3) stop(usage(), "missing positional argument(s)")

# assign arguments to variables
gwas_folder <- args[1]
chr_info <- args[2]
centromere_info <- args[3]
# gwas_folder <- "analysis/gwas"
# chr_info <- "data/B73_RefGen_V4_chrm_info.txt"
# centromere_info <- "data/centromeres_Schneider-2016-pnas_v4.bed"

print('Summarizing gwas')

#### summarize gwas ----

# create empty dataframe to store results
summary_gwas <- data.frame(stringsAsFactors = FALSE)

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
      sig_hits <- list.files(model_folder, full.names = TRUE)
      sig_hits <- sig_hits[grep("sig-hits.csv", sig_hits)]

      # get significant hit info, if any
      if (length(sig_hits) > 0) {

        sig_hits <- fread(sig_hits, header = TRUE, data.table = FALSE)
        sig_hits <- cbind(trait = trait, model = model, sig_hits)

      } else {

        sig_hits <- data.frame(trait = trait, model = model, SNP = NA, Chromosome = NA,
                               Position = NA, maf = NA, effect = NA, P.value = NA,
                               P.value.threshold = NA, P.value.method = NA)

      }
    }

    # append to final summary table
    summary_gwas <- rbind(summary_gwas, sig_hits)

  }
}

print('Writing summary of data')

# write summary
outfile_summary <- paste0(gwas_folder, "/summary_gwas.txt")
fwrite(summary_gwas, file = outfile_summary, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

print('Plotting summary')

# plot summary
summary_plot <- summary_gwas %>%
  filter(P.value.method == "permutation") %>%
  group_by(trait, model) %>%
  summarize(SNPs = length(grep("^del|^dup|^ins|^inv|^tra", SNP, perl = TRUE, invert = TRUE)) - sum(is.na(SNP)),
            SVs = length(grep("^del|^dup|^ins|^inv|^tra", SNP, perl = TRUE))) %>%
  pivot_longer(cols = c(SNPs, SVs), names_to = "sig_markers", values_to = "count") %>%
  ungroup() %>%
  ggplot() +
  geom_col(aes(x = trait, y = count, fill = sig_markers), position = "dodge") +
  facet_wrap(~ model) +
  labs(x = "Trait/Env",
       y = "Number of significant markers") +
  scale_y_continuous(n.breaks = 10) +
  scale_fill_discrete(name = "Marker type") +
  theme_bw() +
  theme(text = element_text(size = 15),
        legend.position = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
ggsave(summary_plot, filename = paste0(gwas_folder, "/gwas_sig_hits.pdf"),
       height = 6, width = 14)


# plot distribution of effect sizes
effects_plot <- summary_gwas %>%
  filter(P.value.method == "permutation") %>%
  separate(col = trait, into = c("trait", "env"), sep = "_") %>%
  filter(!is.na(effect)) %>%
  group_by(trait) %>%
  mutate(to_keep = if_else(length(effect) > 1, "yes", "no")) %>%
  ungroup() %>%
  filter(to_keep == "yes") %>%
  ggplot() +
  geom_density(aes(x = effect, color = model)) +
  facet_wrap(~ trait, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(effects_plot, filename = paste0(gwas_folder, "/gwas_effects_sig-hits.per-trait.pdf"),
       height = 6, width = 14)




print('Making karyoplot')

#### base for karyoplot ----

print('Defining chromosomes')
# chromosomes
chrms <- fread(chr_info, header = TRUE, data.table = FALSE)

print('Converting to data frame')
chrms <- data.frame(chr = chrms$chr, start_pos = 0, end_pos = chrms$length)

print('Defining centromeres')
# centromeres
centros <- read.delim(centromere_info, sep = "\t", header = TRUE)
centros <- data.frame(chr = centros$chr, start_pos = centros$start_pos, end_pos = centros$end_pos)

print('Formatting summary for plotting')
# format summary for plotting
gwas_hits <- subset(summary_gwas, P.value.method == "permutation")
gwas_hits <- gwas_hits[, c("SNP", "Chromosome", "Position", "trait", "model")]
gwas_hits <- separate(gwas_hits, col = "trait", sep = "_", into = c("trait", "env"))
colnames(gwas_hits) <- c("id", "chr", "pos", "trait", "env", "model")
gwas_hits <- gwas_hits[!is.na(gwas_hits$id), ]

print('Converting attributes to factors')
# make sure attributes are factors
gwas_hits$trait <- factor(gwas_hits$trait)
# get levels of each attribute
trait_levels <- levels(gwas_hits$trait)

print('Plotting blank karyotypes')
# plot blank karyotypes
karyoplot <- ggplot() +
  geom_segment(data = chrms,
               aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
               lineend = "round", color = "Gainsboro", size = 10) +
  geom_segment(data = centros,
               aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
               lineend = "round", color = "DimGray", size = 10) +
  facet_grid(~chr, switch = "y") +
  scale_y_reverse(breaks = seq(0, 3.5e8, 0.50e8), labels = c(1, seq(50, 350, 50))) +
  scale_x_discrete(position = "top") +
  coord_cartesian(xlim = c(-3, 3)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.caption = element_text(size = rel(1.1), color = "DimGray"),
        axis.text = element_text(size = rel(2)),
        axis.title = element_text(size = rel(2)),
        strip.text.x = element_text(size = rel(2)),
        legend.position = "bottom") +
  labs(x = "Chromosomes",
       y = "Genomic positions (Mb)")


print('Karyoplot for all hits')
#### karyoplot for all hits ----

# get colors and positions for each attribute
trait_colors <- brewer.pal(n = 8, name = 'Dark2')[2:6]
trait_positions <- seq(-1, 1, length.out = length(trait_levels))

karyoplot_all_hits <- karyoplot
for (trait in 1:length(trait_levels)) {

  # filter data to have only specific attribute
  gwas_hits_filtered <- gwas_hits[which(gwas_hits[, "trait"] == trait_levels[trait]), ]
  trait_color <- trait_colors[trait]
  gwas_hits_filtered$point_pos <- trait_positions[trait]

  # add attributes one at a time
  if (NROW(gwas_hits_filtered) > 0) {

    karyoplot_all_hits <- karyoplot_all_hits +
      geom_point(data = gwas_hits_filtered, aes(x = point_pos, y = pos, shape = model),
                 color = trait_color, size = 2, show.legend = FALSE) +
      scale_shape_manual(values = c(16, 1))

  }
}

# open graphic device
pdf(paste0(gwas_folder, "/karyoplot_gwas_hits.pdf"), width = 10, height = 8)

grid.newpage()
grid.draw(arrangeGrob(karyoplot_all_hits))

for (trait in 1:length(trait_levels)) {
  grid.rect(x = 0.9, y = 0.35 - (trait / 20), width = 0.03, height = 0.01,
            gp = gpar(col = NA, fill = trait_colors[trait]))
  grid.text(x = 0.93, y = 0.35 - (trait / 20), label = trait_levels[trait], hjust = 0,
            gp = gpar(fontsize = 12))
}

grid.circle(x = 0.75, y = 0.3, r = 0.01, gp = gpar(col = NA, fill = "gray50"))
grid.text(x = 0.78, y = 0.3, label = "additive", hjust = 0, gp = gpar(fontsize = 12))
grid.circle(x = 0.75, y = 0.25, r = 0.01, gp = gpar(col = "gray50", fill = NA))
grid.text(x = 0.78, y = 0.25, label = "dominant", hjust = 0, gp = gpar(fontsize = 12))

# close graphic device for saving plot
dev.off()



print('Karyoplot for trait hits')
#### karyoplot for trait hits ----

# get colors for each trait
trait_colors <- brewer.pal(n = 8, name = 'Dark2')[2:6]

for (trait in 1:length(trait_levels)) {

  # create blank karyoplot
  karyoplot_trait_hits <- karyoplot

  # filter data to have only specific trait
  gwas_hits_filtered <- gwas_hits[which(gwas_hits[, "trait"] == trait_levels[trait]), ]
  trait_color <- trait_colors[trait]

  # define colors and positions for envs
  gwas_hits_filtered$env <- factor(gwas_hits_filtered$env)
  env_levels <- mixedsort(levels(gwas_hits_filtered$env))
  env_colors <- lightness(trait_color, scalefac(seq(0.5, 1.5, length.out = length(env_levels))))
  env_positions <- seq(-1, 1, length.out = length(env_levels))

  for (env in 1:length(env_levels)) {

    # filter data to have only specific env
    gwas_hits_filtered_env <- gwas_hits_filtered[which(gwas_hits_filtered[, "env"] == env_levels[env]), ]
    env_color <- env_colors[env]
    gwas_hits_filtered_env$point_pos <- env_positions[env]

    # add attributes one at a time
    if (NROW(gwas_hits_filtered_env) > 0) {

      karyoplot_trait_hits <- karyoplot_trait_hits +
        geom_point(data = gwas_hits_filtered_env, aes(x = point_pos, y = pos, shape = model),
                   color = env_color, size = 2, show.legend = FALSE) +
        scale_shape_manual(values = c(16, 1))

    }
  }

  # open graphic device
  pdf(paste0(gwas_folder, "/karyoplot_gwas_hits.", trait_levels[trait], ".pdf"), width = 10, height = 8)

  grid.newpage()
  grid.draw(arrangeGrob(karyoplot_trait_hits))

  for (env in 1:length(env_levels)) {

    if (length(env_levels) <= length(seq(0.35, 0.05, by = -0.05))) {
      # if not many envs, it's possible to put legend equally spaced
      grid.rect(x = 0.85, y = 0.35 - (env / 20), width = 0.03, height = 0.01,
                gp = gpar(col = NA, fill = env_colors[env]))
      grid.text(x = 0.88, y = 0.35 - (env / 20), label = env_levels[env], hjust = 0,
                gp = gpar(fontsize = 12))
    } else {
      y_pos <- seq(0.35, 0.05, length.out = length(env_levels))
      grid.rect(x = 0.85, y = y_pos[env], width = 0.03, height = 0.01,
                gp = gpar(col = NA, fill = env_colors[env]))
      grid.text(x = 0.88, y = y_pos[env], label = env_levels[env], hjust = 0,
                gp = gpar(fontsize = 12))
    }

  }

  grid.circle(x = 0.70, y = 0.30, r = 0.01, gp = gpar(col = NA, fill = "gray50"))
  grid.text(x = 0.73, y = 0.30, label = "additive", hjust = 0, gp = gpar(fontsize = 12))
  grid.circle(x = 0.70, y = 0.25, r = 0.01, gp = gpar(col = "gray50", fill = NA))
  grid.text(x = 0.73, y = 0.25, label = "dominant", hjust = 0, gp = gpar(fontsize = 12))

  # close graphic device for saving plot
  dev.off()


}
print('Finished')
