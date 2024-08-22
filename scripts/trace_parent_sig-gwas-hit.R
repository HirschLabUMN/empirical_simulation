library(data.table)
library(tibble)
library(readxl)
library(gtools)
library(tidyr)
library(ggplot2)

usage <- function() {
  cat("
description: trace the parental inbred line that contributed to a favorable allele
             from a marker in hybrids

usage: trace_parent_sig-gwas-hit.R [ld_file] [sig_marker] [outfolder] [...]

positional arguments:
  geno_parents_file         genotypes of parental lines (row: genotype calls for
                            that marker, col: parental lines)
  geno_hybrids_file         genotypes of hybrids (row: genotype calls for
                            that marker, col: hybrids)
  geno_hybrids_gwas_file    two-column table with taxa name in the first column,
                            and numeric genotype calls for the marker in the
                            second column (should be the same numeric genotypes
                            from GWAS)
  hybrids_info_file         excel file containing information about parents of
                            each hybrid
  pheno_file                phenotypic values for a trait across multiple environments
                            (first col: taxa, remaining cols: environments)
  outfolder                 output folder to save plots

optional argument:
  --help                    show this helpful message

"
  )
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 6) stop(usage(), "missing positional argument(s)")

# get arguments
geno_parents_file <- args[1]
geno_hybrids_file <- args[2]
geno_hybrids_gwas_file <- args[3]
hybrids_info_file <- args[4]
pheno_file <- args[5]
outfolder <- args[6]
if (!dir.exists(outfolder)) dir.create(outfolder)
# geno_parents_file <- "analysis/sig_marker_yield/geno_parents.markers-ld-sig-hit.hmp.txt"
# geno_hybrids_file <- "analysis/sig_marker_yield/geno_hybrids.markers-ld-sig-hit.hmp.txt"
# geno_hybrids_gwas_file <- "analysis/sig_marker_yield/geno_hybrids.markers-ld-sig-hit.add-major-imputed.num.txt"
# hybrids_info_file <- "data/NIFA2020_CombinedData.xlsx"
# pheno_file <- "data/1stStage_BLUEs.YLD-per-env.txt"
# outfolder <- "analysis/sig_marker_yield"



#### load and format data ----

# load genotypic data
geno_parents <- fread(geno_parents_file, header = TRUE, data.table = FALSE)
geno_hybrids <- fread(geno_hybrids_file, header = TRUE, data.table = FALSE)
# reformat genotypic data
geno_parents <- rownames_to_column(data.frame(geno = t(geno_parents)), var = "hybrid")
geno_hybrids <- rownames_to_column(data.frame(geno = t(geno_hybrids)), var = "hybrid")

# load genotypic data used in gwas
geno_hybrids_gwas <- fread(geno_hybrids_gwas_file, header = TRUE, data.table = FALSE)
colnames(geno_hybrids_gwas) <- c("hybrid", "geno_gwas")
# merge hybrid geno data
geno_hybrids <- merge(x = geno_hybrids, y = geno_hybrids_gwas, by = "hybrid", all = TRUE)
geno_hybrids <- geno_hybrids[mixedorder(geno_hybrids$hybrid), ]
rm(geno_hybrids_gwas)

# # identify major/minor allele
# allele_count <- as.character(geno_hybrids[geno_hybrids$geno != "NN", "geno"])
# allele_count <- paste0(allele_count, collapse = "")
# allele_count <- data.frame(table(strsplit(allele_count, split = "")))

# load phenotypic data
pheno <- fread(pheno_file, header = TRUE, data.table = FALSE)

# load hybrid information
hybrid_info <- read_excel(hybrids_info_file)
# keep only pedigree information
hybrid_info <- hybrid_info[, c("Hybrid", "ParentA", "ParentB")]
hybrid_info <- hybrid_info[!duplicated(hybrid_info$Hybrid), ]
# remove checks
hybrid_info <- hybrid_info[grep("check", hybrid_info$Hybrid, invert = TRUE, ignore.case = TRUE), ]
# adjust format
hybrid_info <- as.data.frame(hybrid_info[mixedorder(hybrid_info$Hybrid), ])
colnames(hybrid_info) <- c("hybrid", "p1", "p2")



#### merge information and plot ----

# merge all hybrid info
hybrid_info <- merge(x = hybrid_info, y = geno_hybrids, by = "hybrid", all = TRUE)
hybrid_info <- merge(x = hybrid_info, y = pheno, by = "hybrid", all = TRUE)
hybrid_info <- hybrid_info[mixedorder(hybrid_info$hybrid), ]
rm(geno_hybrids, pheno)

# transform to long format
hybrid_info <- pivot_longer(hybrid_info, cols = !hybrid:geno_gwas, names_to = "env", values_to = "pheno")

# get trait name
pheno_name <- rev(unlist(strsplit(pheno_file, split = "/")))[1]
pheno_name <- gsub(".txt", "", pheno_name)

# match genotype to respective numeric genotype
names(hybrid_info$geno_gwas) <- hybrid_info$geno
labels_geno_gwas <- hybrid_info$geno_gwas[!duplicated(names(hybrid_info$geno_gwas))]
labels_geno_gwas <- labels_geno_gwas[names(labels_geno_gwas) != "NN"]
# reorder factors
hybrid_info$geno <- factor(hybrid_info$geno,
                           levels = sort(unique(as.character(hybrid_info$geno))))
hybrid_info$geno_gwas <- factor(hybrid_info$geno_gwas,
                                levels = as.character(sort(unique(hybrid_info$geno_gwas))),
                                labels = names(sort(labels_geno_gwas)))
hybrid_info$geno_gwas <- factor(hybrid_info$geno_gwas,
                                levels = sort(levels(hybrid_info$geno_gwas)))

# function to count number of data points per group
n_per_bin <- function(x) {
  return(data.frame(y = Inf, label = paste('n =', length(x))))
}

# plot results
plot_non_imp <- ggplot(hybrid_info, aes(x = geno, y = pheno)) +
  geom_boxplot() +
  facet_wrap(~ env, nrow = 2) +
  labs(x = "Genotypes before imputation",
       y = pheno_name) +
  stat_summary(fun.data = n_per_bin, geom = "text", size = 2, vjust = 1)
ggsave(filename = paste0(outfolder, "/sig-hit_non-imp-geno_vs_", pheno_name, ".pdf"),
       plot = plot_non_imp, device = "pdf", width = 13, height = 9)

plot_imp <- ggplot(hybrid_info, aes(x = geno_gwas, y = pheno)) +
  geom_boxplot() +
  facet_wrap(~ env, nrow = 2) +
  labs(x = "Genotypes after imputation",
       y = pheno_name) +
  stat_summary(fun.data = n_per_bin, geom = "text", size = 2, vjust = 1)
ggsave(filename = paste0(outfolder, "/sig-hit_imp-geno_vs_", pheno_name, ".pdf"),
       plot = plot_imp, device = "pdf", width = 13, height = 9)
