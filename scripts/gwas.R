library(data.table)
library(bigmemory)
library(biganalytics)
library(compiler)
library(doParallel)
library(foreach)
library(doRNG)
library(multtest)
library(snpStats)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler)
library(scatterplot3d)
library(bigmemory)
library(biganalytics)
source("/home/hirschc1/burns756/empirical_sim/GAPIT_Source_Code/GAPIT_Source_Code.txt")
source("/home/hirschc1/burns756/empirical_sim/GAPIT_Source_Code/FarmCPU_Source_Code.txt")
library(ggplot2)
library(dplyr)
library(ggrepel)

usage <- function() {
  cat("
description: run GWAS with GAPIT.

usage: Rscript gwas.R [pheno_file] [geno_file] [geno_info] [pca_file] [output_folder]

positional arguments:
  pheno_file                file with phenotypic data (2 columns: taxa, value)
  geno_file                 file with numeric genotypic data (first column should be taxa)
  geno_info                 file with information about markers (3 columns: name, chr, pos)
  kinship                   file with kinship among individuals (first column should be taxa; no column names)
  pca_file                  file with PCA results from genotypic data (first column should be taxa)
  output_folder             folder to save files

optional argument:
  --help                    show this helpful message
  --PCs=VALUE               number of principal components to be used in GWAS (default: 5)
  --gwas-model=VALUE        type of GWAS model available in GAPIT. Available models:
                            'GLM', 'MLM' (default), 'CMLM', 'SUPER', 'MLMM', 'FarmCPU', 'Blink'
  --pval-method=VALUE       define significance threshold by 'FDR' (default) or 'permutation'


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

# set default
PCs <- "5"
gwas_model <- "MLM"
pval_method <- "FDR"

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 6) stop(usage(), "missing positional argument(s)")

if (length(args) > 6) {

  opt_args <- args[-1:-6]
  opt_args_allowed <- c("--PCs", "--gwas-model", "--pval-method")
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
if (suppressWarnings(!is.na(as.integer(PCs)))) {
  PCs <- as.integer(PCs)
} else {
  stop("Optional argument '--PCs' should be an integer")
}

if (!gwas_model %in% c("GLM", "MLM", "CMLM", "SUPER", "MLMM", "FarmCPU", "Blink")) {
  stop("Optional argument '--gwas-model' should be one of these models: 'GLM', 'MLM', 'CMLM', 'SUPER', 'MLMM', 'FarmCPU', 'Blink'")
}

if (!pval_method %in% c("FDR", "permutation")) {
  stop("Optional argument '--gwas-model' should be one of these methods: 'FDR' or 'permutation'")
}

# get positional arguments
pheno_file <- args[1]
geno_file <- args[2]
geno_info <- args[3]
kinship_file <- args[4]
pca_file <- args[5]
output_folder <- args[6]



#### load data ----

# read files
pheno <- fread(pheno_file, header = TRUE, data.table = FALSE)
geno_map <- fread(geno_info, head = TRUE, data.table = FALSE)
geno_data <- fread(geno_file, header = TRUE, data.table = FALSE)
kinship <- fread(kinship_file, header = FALSE, data.table = FALSE)
pca <- fread(pca_file, head = TRUE, data.table = FALSE)
# save initial path
init_folder <- getwd()

# adjust pca to number of desired PCs and exclude taxa column
pca <- pca[, 1:(PCs + 1)]



#### p-value threshold ----

pvalThreshold <- function(GD = NULL, GM = NULL, Y = NULL, CV = NULL, KI = NULL,
                          model = "MLM", trait = "", theRep = 100) {

  # I adapted the 'FarmCPU.P.Threshold()' function created by Xiaolei Liu
  # to do permutations on other GAPIT models as well

  if (is.null(GD)) return(NULL)
  if (is.null(GM)) return(NULL)
  if (is.null(Y)) return(NULL)

  # run permutations in parallel
  num.cores <- detectCores()
  registerDoParallel(cores = num.cores)
  registerDoRNG(12345)
  pvalue.final <- foreach(i = 1:theRep, .combine = c) %dopar% {

    # shuffle pheno data
    index = 1:nrow(Y)
    index.shuffle = sample(index, length(index), replace = FALSE)
    Y.shuffle = Y
    Y.shuffle[, 2] = Y.shuffle[index.shuffle, 2]
    # GWAS with GAPIT
    myGAPIT <- GAPIT(Y = Y.shuffle[, c(1, 2)],  # phenotype
                     GD = GD,                   # genotype
                     GM = GM,                   # map information
                     KI = KI,                   # kinship
                     CV = CV,                   # covariates
                     model = model,             # gwas model to use
                     file.output = FALSE)
    # get lowest pvalue
    pvalue <- min(myGAPIT$GWAS[, 4], na.rm = TRUE)
    # return pvalue
    pvalue

  }
  stopImplicitCluster()
  pvalue.final <- as.numeric(pvalue.final)

  # determine pvalue threshold
  pval.to.use <- sort(pvalue.final)[ceiling(theRep * 0.05)]
  print("The p.threshold of this data set should be:")
  print(pval.to.use)

  return(pval.to.use)

}

if (pval_method == "permutation") {
  # From FarmCPU manual: 'In this function, the phenotypes are permuted
  # to break the relationship with the genotypes. The experiment is
  # replicated for a number of times. A vector of minimum p value of
  # each experiment is outputted and the 95% quantile value of the vector
  # is recommended for p.threshold in FarmCPU model.'

  # here I adapted original FarmCPU.P.Threshold() function to do
  # permutations on other GAPIT models as well

  # calcuate p-value threshold based on trait of interest
  # with the number of permutation times set to 30

  tmp_folder <- paste0(output_folder, "/tmp")
  if (!dir.exists(tmp_folder)) dir.create(tmp_folder)
  setwd(tmp_folder)
  pval <- pvalThreshold(Y = pheno[, c(1, 2)],
                        GD = geno_data,
                        GM = geno_map,
                        KI = kinship,
                        CV = pca,
                        model = gwas_model,
                        trait = colnames(pheno)[2],
                        theRep = 30)
  setwd(init_folder)
  system(paste0("rm ", tmp_folder, "/*"))
  system(paste0("rmdir ", tmp_folder))
}

if (pval_method == "FDR") {
  # FDR-adjusted threshold of 0.05
  pval <- 0.05
}



#### gwas ----

# run gwas
setwd(output_folder)
gwas <- GAPIT(Y = pheno,
              GD = geno_data,
              GM = geno_map,
              KI = kinship,
              CV = pca,
              model = gwas_model)
setwd(init_folder)

# get gwas table
gwas_results <- gwas$GWAS
colnames(gwas_results) <- gsub(" ", "", colnames(gwas_results))

# adjust pvalue if FDR
if (pval_method == "FDR") {
  gwas_results$P.value <- p.adjust(gwas_results$P.value, method = "fdr")
}

# select signficant hits
sig_hits <- subset(gwas_results, P.value <= pval)

if (NROW(sig_hits) > 0) {
  sig_hits <- cbind(sig_hits[, c("SNP", "Chromosome", "Position", "maf", "effect", "P.value")],
                    P.value.threshold = pval,
                    P.value.method = pval_method)
  # create output name
  outfile_sig_hits <- paste0(output_folder, "/", gwas_model, ".", colnames(pheno)[2], ".GWAS.Results.sig-hits.csv")
  # write file
  fwrite(sig_hits, file = outfile_sig_hits, quote = FALSE, sep = ",", na = NA, row.names = FALSE)
}



#### plot ----

gwas_results <- gwas_results %>%
  # compute chromosome size
  group_by(Chromosome) %>%
  summarize(chr_len = max(Position)) %>%
  # calculate cumulative position of each chromosome
  mutate(chr_cum = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  # add this info to the initial dataset
  left_join(gwas_results, ., by = c("Chromosome" = "Chromosome")) %>%
  # add a cumulative position of each SNP,
  arrange(Chromosome, Position) %>%
  mutate(pos_cum = Position + chr_cum) %>%
  # also add which markers are SVs and which ones should receive annotation
  mutate(log_pval = -log10(P.value),
         is_sv = ifelse(grepl("^del|^dup|^ins|^inv|^tra", SNP, perl = TRUE), "yes", "no"),
         is_annotate = ifelse(log_pval > -log10(pval), "yes", "no"))

# prepare x axis to not want to display the cumulative position of SNP in bp
# (just show the chromosome name instead)
axis_plot <- gwas_results %>%
  group_by(Chromosome) %>%
  summarize(center = (as.numeric(max(pos_cum)) + as.numeric(min(pos_cum))) / 2)

# basic manhattan plot
gwas_plot <- gwas_results %>%
  # # uncomment this if need to speed up plotting (skip markers with very high pvalue)
  # filter(-log10(P.value) > 1) %>%
  ggplot(aes(x = pos_cum, y = log_pval)) +
  # show all points
  geom_point(aes(color = as.factor(Chromosome)), alpha = 0.8, size = 3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), length(unique(gwas_results$Chromosome)))) +
  # add pvalue threshold
  geom_hline(yintercept = -log10(pval), color = "firebrick") +
  # # add highlighted points
  # geom_point(data = subset(gwas_results, is_sv == "yes"), color = "orange", size = 2) +
  # add label using ggrepel to avoid overlapping
  geom_label_repel(data = subset(gwas_results, is_annotate == "yes"), aes(label = SNP), size = 2) +
  # custom X axis
  scale_x_continuous(label = axis_plot$Chromosome, breaks = axis_plot$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ceiling(max(gwas_results$log_pval, na.rm = TRUE)))) +
  # add labels
  labs(title = colnames(pheno)[2],
       x = "Chromosome",
       y = "-log10(p-value)") +
  # custom theme
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        text = element_text(size = 20))

# save plot
outfile_plot <- paste0(output_folder, "/", gwas_model, ".", colnames(pheno)[2], ".Manhattan.Plot.customized.pdf")
# write file
ggsave(filename = outfile_plot, plot = gwas_plot, device = "pdf", width = 16, height = 8)



#### debug ----

# pheno_file <- "tests/YLD-MIN20.txt"
# geno_file <- "data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.prunned-1000kb.dominant.num.txt"
# geno_info <- "data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.prunned-1000kb.dominant.map.txt"
# pca_file <- "analysis/pca_gwas/pruning_1000kb/dominant_model/GAPIT.PCA.csv"
# kinship_file <- "analysis/kinship/pruning_1000kb/dominant_model/GAPIT.Kin.VanRaden.csv"
# output_folder <- "tests/gwas_models"
# PCs <- 5
# gwas_model <- "MLM"
# pval_method <- "permutation"
