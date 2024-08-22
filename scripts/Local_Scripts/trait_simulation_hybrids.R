library(data.table)
library(tidyr)
library(ggplot2)
library(gtools)
library(dplyr)
suppressWarnings(suppressMessages(library(simplePHENOTYPES)))
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
description: simulate trait in hybrids based on a GWAS genetic architecture.

usage: Rscript trait_simulation_hybrid.R [trait] [gwas_file] [geno_file] [sv_list]
                                         [res_cor_matrix_file] [h2_file] [means_file]
                                         [outfolder] [real_pheno_file] [...]
positional arguments:
  trait                           name of trait to be simulated (need to match the name on GWAS)
  gwas_file                       file with non-redundant top GWAS hits summarized by environment (output of
                                  script 'merge_gwas_per_trait.R')
  geno_file                       hapmap file containing SNPs and SVs
  sv_list                         single-column file containing only SV IDs
  res_cor_matrix_file             residual correlation matrix among environments to be simulated (output
                                  of script 'get_usda_env-types.R')
  h2_file                         file with heritabilities for each trait (column) in each environment (row)
  means_file                      file with trait means in each environment (three columns: trait, env, mean)
  outfolder                       path to folder to simulated traits and plots
  real_pheno_file                 file with real phenotypic values summarized by environment (output of
                                  script 'merge_gwas-and-pheno_per_trait.R')

optional argument:
  --help                          show this helpful message
  --effect-size-reduction=VALUE   decrease the effect size of extra QTLs by 0 < VALUE <= 1 (default: 0.1)
  --reps=VALUE                    number of reps (default: 3)
  --impute-type=VALUE             the marker type ('Major', 'Middle' or 'Minor') when imputing missing data
                                  at the hapmap numericalization step (default: 'Middle')
  --architecture=VALUE            genetic architecture to be simulated ('pleiotropic' for traits being
                                  controlled by the same QTNs, 'partially' for traits being controlled by
                                  pleiotropic and trait-specific QTNs, 'LD' for traits being exclusively
                                  controlled by different QTNs)
  --seed=VALUE                    value for set.seed (default: NULL; random number is selected)
  --QTN-variance                  add this option to write files with percent variance explained per QTN
  --minor-qtls                    add this option to add non-significant markers to simulations
  --reduce-effect-all             add this option to reduce the effects of both major and minor QTLs


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
trait <- args[1]
gwas_file <- args[2]
geno_file <- args[3]
sv_list <- args[4]
res_cor_matrix_file <- args[5]
h2_file <- args[6]
means_file <- args[7]
outfolder <- args[8]
real_pheno_file <- args[9]

# set default
effect_size_reduction <- "0.1"
reps <- "3"
impute_type <- "Middle"
architecture <- "pleiotropic"
seed <- NULL
QTN_variance <- FALSE
minor_qtls <- FALSE
reduce_effect_all <- FALSE

# assert to have the correct optional arguments
pos_args <- 9
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {

  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--effect-size-reduction", "--reps", "--impute-type",
                        "--architecture", "--seed", "--QTN-variance",
                        "--minor-qtls", "--reduce-effect-all")
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
if (suppressWarnings(!is.na(as.numeric(effect_size_reduction)))) {
  effect_size_reduction <- as.numeric(effect_size_reduction)
  if (effect_size_reduction <= 0 | effect_size_reduction > 1) {
    stop("Optional argument '--effect-size-reduction' should be 0 < VALUE <= 1")
  }
} else {
  stop("Optional argument '--effect-size-reduction' should be 0 < VALUE <= 1")
}

if (suppressWarnings(!is.na(as.integer(reps)))) {
  reps <- as.integer(reps)
} else {
  stop("Optional argument '--reps' should be an integer")
}

if (!impute_type %in% c("Major", "Middle", "Minor")) {
  stop("Optional argument '--impute-type' should be 'Major', 'Middle' or 'Minor'")
}

if (!architecture %in% c("pleiotropic", "partially", "LD")) {
  stop("Optional argument '--architecture' should be 'pleiotropic', 'partially' or 'LD'")
}

if (is.null(seed)) {
  seed <- ceiling(runif(1, 0, 1000000))
} else {
  if (suppressWarnings(!any(is.na(as.numeric(seed))))) {
    seed <- as.numeric(seed)
  } else {
    stop("Optional argument '--seed' should be a number")
  }
}



#### define genetic architecture ----

cat("loading GWAS data\n")

# load gwas data
gwas_hits <- fread(gwas_file, header = TRUE, data.table = FALSE)
# add whether marker is SNP or SV
gwas_hits$marker_type <- sapply(gwas_hits$marker, function(marker) {
  type <- ifelse(grepl("^del|^dup|^inv|^ins|^tra", marker, perl = TRUE),
                 yes = "SV", no = "SNP")
  return(type)
})

# get all envs with significant gwas hits
envs_names <- unique(gwas_hits$env)
envs <- length(envs_names)


cat("defining genetic architecture\n")

# get vector with type of marker effects from significant hits
sig_add_dom_hits <- gwas_hits$type[gwas_hits$signif]
# define proportion of signficant gwas hits that are add vs dom
# this option will be used to determine the type of model to simulate traits
if (length(sig_add_dom_hits) > 0) {
  add_dom_ratio <- data.frame(model = c("additive", "dominant"))
  add_dom_ratio <- merge(x = add_dom_ratio, y = data.frame(table(sig_add_dom_hits)),
                         by.x = "model", by.y = "sig_add_dom_hits", all = TRUE)
  add_dom_ratio <- as.numeric(add_dom_ratio[add_dom_ratio$model == "additive", "Freq"]) / sum(add_dom_ratio$Freq, na.rm = TRUE)
  add_dom_ratio <- ifelse(is.na(add_dom_ratio), yes = 0, no = round(add_dom_ratio, digits = 4))
} else {
  # if there are no significant hits for this trait, then use a 50:50 add-dom ratio
  add_dom_ratio <- 0.5
}
rm(sig_add_dom_hits)

# define causal variants with large effect size (gwas hits)
causal_vars_add_major <- unique(gwas_hits[gwas_hits$type == "additive" & gwas_hits$signif == TRUE, "marker"])
causal_vars_dom_major <- unique(gwas_hits[gwas_hits$type == "dominant" & gwas_hits$signif == TRUE, "marker"])

# define major marker effects
add_effect_value_major <- list()
dom_effect_value_major <- list()
for (env in 1:envs) {
  # add additive effects
  gwas_hits_env_add <- gwas_hits[gwas_hits$env == envs_names[env] & gwas_hits$type == "additive", ]
  gwas_hits_env_add <- gwas_hits_env_add[match(causal_vars_add_major, gwas_hits_env_add$marker), ]
  add_effect_value_major[[env]] <- gwas_hits_env_add[, "effect"]
  add_effect_value_major[[env]] <- sapply(add_effect_value_major[[env]], function(x) if (is.na(x)) return(0) else return(x))
  # add dominance effects
  gwas_hits_env_dom <- gwas_hits[gwas_hits$env == envs_names[env] & gwas_hits$type == "dominant", ]
  gwas_hits_env_dom <- gwas_hits_env_dom[match(causal_vars_dom_major, gwas_hits_env_dom$marker), ]
  dom_effect_value_major[[env]] <- gwas_hits_env_dom[, "effect"]
  dom_effect_value_major[[env]] <- sapply(dom_effect_value_major[[env]], function(x) if (is.na(x)) return(0) else return(x))
}
rm(gwas_hits_env_add, gwas_hits_env_dom)

if (reduce_effect_all) {
  # decrease effect size of major qtls, if requested
  add_effect_value_major <- lapply(add_effect_value_major, function(x) if(length(x) > 0) return(x * effect_size_reduction) else return(list()))
  dom_effect_value_major <- lapply(dom_effect_value_major, function(x) if(length(x) > 0) return(x * effect_size_reduction) else return(list()))
}


if (minor_qtls) {

  # define causal variants with smaller effect size (top non-significant gwas hits)
  causal_vars_add_minor <- unique(gwas_hits[gwas_hits$type == "additive" & gwas_hits$signif == FALSE, "marker"])
  causal_vars_dom_minor <- unique(gwas_hits[gwas_hits$type == "dominant" & gwas_hits$signif == FALSE, "marker"])

  # define minor marker effects
  add_effect_value_minor <- list()
  dom_effect_value_minor <- list()
  for (env in 1:envs) {
    # add additive effects
    gwas_hits_env_add_minor <- gwas_hits[gwas_hits$env == envs_names[env] & gwas_hits$type == "additive", ]
    gwas_hits_env_add_minor <- gwas_hits_env_add_minor[match(causal_vars_add_minor, gwas_hits_env_add_minor$marker), ]
    add_effect_value_minor[[env]] <- gwas_hits_env_add_minor[, "effect"]
    add_effect_value_minor[[env]] <- sapply(add_effect_value_minor[[env]], function(x) if (is.na(x)) return(0) else return(x))
    # add dominance effects
    gwas_hits_env_dom_minor <- gwas_hits[gwas_hits$env == envs_names[env] & gwas_hits$type == "dominant", ]
    gwas_hits_env_dom_minor <- gwas_hits_env_dom_minor[match(causal_vars_dom_minor, gwas_hits_env_dom_minor$marker), ]
    dom_effect_value_minor[[env]] <- gwas_hits_env_dom_minor[, "effect"]
    dom_effect_value_minor[[env]] <- sapply(dom_effect_value_minor[[env]], function(x) if (is.na(x)) return(0) else return(x))
  }
  rm(gwas_hits_env_add_minor, gwas_hits_env_dom_minor)

  # decrease effect size of extra qtls
  add_effect_value_minor <- lapply(add_effect_value_minor, function(x) if(length(x) > 0) return(x * effect_size_reduction) else return(list()))
  dom_effect_value_minor <- lapply(dom_effect_value_minor, function(x) if(length(x) > 0) return(x * effect_size_reduction) else return(list()))

  # adjust effects of markers that were significant in one environment, but not others
  add_markers_to_adjust <- intersect(causal_vars_add_major, causal_vars_add_minor)
  if (length(add_markers_to_adjust) > 0) {
    for (marker in add_markers_to_adjust) {

      # get index of causal variant
      causal_var_major_idx <- which(causal_vars_add_major == marker)
      # get names of enviroments to adjust
      envs_marker_to_adjust <- gwas_hits[gwas_hits$marker == marker & gwas_hits$type == "additive" & gwas_hits$signif == FALSE, "env"]
      # adjust effect of marker in that particular environment where it is not significant
      for (env in which(envs_names %in% envs_marker_to_adjust)) {
        add_effect_value_major[[env]][causal_var_major_idx] <- add_effect_value_major[[env]][causal_var_major_idx] * effect_size_reduction
      }

      # remove marker from list of minor causal vars
      causal_var_minor_idx <- which(causal_vars_add_minor == marker)
      causal_vars_add_minor <- causal_vars_add_minor[-causal_var_minor_idx]
      # remove marker from list of minor effects
      for (env in 1:length(add_effect_value_minor)) {
        add_effect_value_minor[[env]] <- add_effect_value_minor[[env]][-causal_var_minor_idx]
      }

    }
    rm(add_markers_to_adjust, envs_marker_to_adjust, marker, env, causal_var_major_idx, causal_var_minor_idx)
  }

  # do the same for dominant markers
  dom_markers_to_adjust <- intersect(causal_vars_dom_major, causal_vars_dom_minor)
  if (length(dom_markers_to_adjust) > 0) {
    for (marker in dom_markers_to_adjust) {

      # get index of causal variant
      causal_var_major_idx <- which(causal_vars_dom_major == marker)
      # get names of enviroments to adjust
      envs_marker_to_adjust <- gwas_hits[gwas_hits$marker == marker & gwas_hits$type == "dominant" & gwas_hits$signif == FALSE, "env"]
      # adjust effect of marker in that particular environment where it is not significant
      for (env in which(envs_names %in% envs_marker_to_adjust)) {
        dom_effect_value_major[[env]][causal_var_major_idx] <- dom_effect_value_major[[env]][causal_var_major_idx] * effect_size_reduction
      }

      # remove marker from list of minor causal vars
      causal_var_minor_idx <- which(causal_vars_dom_minor == marker)
      causal_vars_dom_minor <- causal_vars_dom_minor[-causal_var_minor_idx]
      # remove marker from list of minor effects
      for (env in 1:length(dom_effect_value_minor)) {
        dom_effect_value_minor[[env]] <- dom_effect_value_minor[[env]][-causal_var_minor_idx]
      }

    }
    rm(dom_markers_to_adjust, envs_marker_to_adjust, marker, env, causal_var_major_idx, causal_var_minor_idx)
  }

}

# add effects to list of gwas qnts
add_effect_value <- list()
dom_effect_value <- list()
for (env in 1:envs) {

  if (minor_qtls) {
    add_effect_value[[env]] <- c(add_effect_value_major[[env]], add_effect_value_minor[[env]])
    dom_effect_value[[env]] <- c(dom_effect_value_major[[env]], dom_effect_value_minor[[env]])
  } else {
    add_effect_value[[env]] <- add_effect_value_major[[env]]
    dom_effect_value[[env]] <- dom_effect_value_major[[env]]
  }

}
rm(gwas_hits, add_effect_value_major, dom_effect_value_major, add_effect_value_minor, dom_effect_value_minor)

# transform to list for compatibility with simplePHENOTYPES
if (minor_qtls) {
  QTN_list <- list(add = list(c(causal_vars_add_major, causal_vars_add_minor)),
                   dom = list(c(causal_vars_dom_major, causal_vars_dom_minor)))
} else {
  QTN_list <- list(add = list(causal_vars_add_major),
                   dom = list(causal_vars_dom_major))
}




#### load marker data ----

cat("loading marker data\n")

# load SV information
SVs <- fread(sv_list, header = FALSE, data.table = FALSE)
SVs <- SVs[, 1]

# load genotypic data
geno_data <- fread(geno_file, header = TRUE, data.table = FALSE)
# filter to have only markers from selected qtl list
geno_data <- subset(geno_data, geno_data[, 1] %in% as.character(unlist(QTN_list)))



#### load environmental data ----

cat("loading environmental data\n")

# load env data
res_cor_matrix <- fread(res_cor_matrix_file, header = TRUE, data.table = FALSE, stringsAsFactors = FALSE)
# keep only envs in gwas
rownames(res_cor_matrix) <- colnames(res_cor_matrix)
res_cor_matrix <- res_cor_matrix[envs_names, envs_names]
# transform to matrix
res_cor_matrix <- as.matrix(res_cor_matrix)

# load heritabilities
h2 <- read.csv(h2_file, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
# keep only envs in gwas for trait of interest
h2 <- h2[envs_names, trait]

# load means
means <- fread(means_file, header = TRUE, data.table = FALSE, stringsAsFactors = FALSE)
# keep only envs in gwas for trait of interest
means <- means[means$trait == trait & means$env %in% envs_names, ]
# match env order and get only the means
means <- means[match(envs_names, means$env), "mean"]



#### simulate traits ----

cat("simulating traits\n")

# adjust models and imputation method to be used
if (add_dom_ratio == 1) {
  model <- "A"
  impute_effect <- "Add"
} else if (add_dom_ratio == 0) {
  model <- "D"
  impute_effect <- "Dom"
} else {
  model <- "AD"
  impute_effect <- "Add"
}

# get total number of QTNs
if (minor_qtls) {
  add_QTN_num <- length(causal_vars_add_major) + length(causal_vars_add_minor)
  dom_QTN_num <- length(causal_vars_dom_major) + length(causal_vars_dom_minor)
} else {
  add_QTN_num <- length(causal_vars_add_major)
  dom_QTN_num <- length(causal_vars_dom_major)
}

# simulate traits
sim_pheno <- create_phenotypes(geno_obj = geno_data,
                               QTN_list = QTN_list,
                               rep = reps,
                               # architecture
                               ntraits = envs,
                               h2 = h2,
                               model = model,
                               mean = means,
                               add_QTN_num = add_QTN_num,
                               add_effect = add_effect_value,
                               dom_QTN_num = dom_QTN_num,
                               dom_effect = dom_effect_value,
                               architecture = architecture,
                               sim_method = "custom",
                               vary_QTN = FALSE,
                               cor = NULL,
                               cor_res = res_cor_matrix,
                               seed = seed,
                               # ouput
                               export_gt = TRUE,
                               home_dir = getwd(),
                               output_dir = 'sim',
                               to_r = TRUE,
                               out_geno = NULL,
                               QTN_variance = QTN_variance,
                               quiet = TRUE,
                               # numericalization
                               SNP_effect = impute_effect,
                               SNP_impute = impute_type)



#### plot sim vs real traits ----

cat("plotting results\n")

# create folder to save plots
plots_folder <- paste0(outfolder, "/plots")
if (!dir.exists(plots_folder)) dir.create(plots_folder)

# format sim data
colnames(sim_pheno)[1] <- "hybrid"
colnames(sim_pheno)[2:(length(envs_names) + 1)] <- envs_names
sim_pheno <- pivot_longer(sim_pheno, cols = all_of(envs_names), names_to = "env", values_to = "sim_value")

# plot sim data per rep
plot_sim_per_rep <- ggplot(sim_pheno) +
  geom_density(aes(x = sim_value, color = env, fill = env), alpha = 0.1) +
  facet_grid(~ Rep) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
ggsave(filename = paste0(plots_folder, "/sim_data_per_rep.pdf"),
       plot = plot_sim_per_rep, device = "pdf", width = 10, height = 4)

# load real phenotypic values
real_pheno <- fread(real_pheno_file, header = TRUE, data.table = FALSE)

# get average phenotypic values for simulated data
sim_pheno_avg <- sim_pheno %>%
  group_by(hybrid, env) %>%
  summarize(sim_value_avg = mean(sim_value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(hybrid = as.character(hybrid))
sim_pheno_avg <- sim_pheno_avg[mixedorder(sim_pheno_avg$hybrid), ]

# merge simulated and real data
sim_real_phenos <- merge(x = sim_pheno_avg, y = real_pheno, by = c("hybrid", "env"), all = TRUE)
sim_real_phenos <- sim_real_phenos[mixedorder(sim_real_phenos$hybrid), ]
rm(sim_pheno_avg, real_pheno)
# write file
fwrite(sim_real_phenos, file = paste0(outfolder, "/sim_vs_real_data.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)


# correlation all envs
cor_all_pearson <- round(cor(sim_real_phenos$sim_value_avg, sim_real_phenos$real_pheno,
                             method = "pearson", use = "complete.obs"), digits = 2)
cor_all_spearman <- round(cor(sim_real_phenos$sim_value_avg, sim_real_phenos$real_pheno,
                              method = "spearman", use = "complete.obs"), digits = 2)

# plot correlations
plot_cor <- ggplot(sim_real_phenos) +
  geom_point(aes(x = sim_value_avg, y = real_pheno, color = env)) +
  labs(title = trait,
       subtitle = paste0("overall correlation: ", cor_all_spearman),
       x = "Simulated phenotypic value",
       y = "Real phenotypic value") +
  scale_color_discrete(name = "Environments") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
ggsave(filename = paste0(plots_folder, "/correlation_sim_vs_real_data.pdf"),
       plot = plot_cor, device = "pdf", width = 10, height = 6)

# correlation per env
cor_per_env <- data.frame(stringsAsFactors = FALSE)
for (env in 1:length(envs_names)) {
  sim_real_phenos_env <- sim_real_phenos[sim_real_phenos$env == envs_names[env], ]
  cor_per_env <- rbind(cor_per_env,
                       data.frame(env = envs_names[env],
                                  cor_pearson = cor(sim_real_phenos_env$sim_value_avg, sim_real_phenos_env$real_pheno,
                                                    method = "pearson", use = "complete.obs"),
                                  cor_spearman = cor(sim_real_phenos_env$sim_value_avg, sim_real_phenos_env$real_pheno,
                                                     method = "spearman", use = "complete.obs")))
}
cor_per_env$cor_pearson <- round(cor_per_env$cor_pearson, digits = 2)
cor_per_env$cor_spearman <- round(cor_per_env$cor_spearman, digits = 2)

# get thrshold for top 10 percent lines in real data
top10_perc_threshold <- ceiling(length(unique(sim_real_phenos$hybrid)) * 0.1)

# plot sim vs real data
plot_sim_vs_real <- sim_real_phenos %>%
  group_by(env) %>%
  mutate(rank_env_real = rank(-real_pheno, na.last = "keep"),
         top10_env_real = if_else(condition = (rank_env_real <= top10_perc_threshold & !is.na(rank_env_real)),
                                  true = TRUE, false = FALSE)) %>%
  ungroup() %>%
  pivot_longer(cols = c("sim_value_avg", "real_pheno"), names_to = "data_type", values_to = "pheno_value") %>%
  mutate(data_type = factor(data_type, levels = c("sim_value_avg", "real_pheno"), labels = c("Simulated", "Real"))) %>%
  ggplot() +
  geom_density(aes(x = pheno_value, color = data_type)) +
  geom_text(data = cor_per_env,
            mapping = aes(x = -Inf, y = Inf, label = paste0("env correlation: ", cor_spearman)),
            hjust = -0.05, vjust = 1.5, size = 3) +
  facet_wrap(~ env, nrow = 2) +
  geom_rug(data = function(x) subset(x, top10_env_real == TRUE & data_type == "Simulated"),
           aes(x = pheno_value, color = data_type), alpha = 0.3) +
  labs(title = trait,
       subtitle = paste0("overall correlation: ", cor_all_spearman),
       x = "Phenotypic value",
       y = "Density") +
  scale_color_discrete(name = "Data type") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
ggsave(filename = paste0(plots_folder, "/sim_vs_real_data.pdf"),
       plot = plot_sim_vs_real, device = "pdf", width = 10, height = 6)

cat("simulated vs real data correlation:", cor_all_pearson, "(Pearson) /", cor_all_spearman, "(Spearman)\n\n")
print(cor_per_env)



#### debug ----

# trait <- "PHT"
# gwas_file <- "data/gwas_top-peaks_PHT-per-env.txt"
# geno_file <- "data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.hmp.txt"
# sv_list <- "data/SVs_IDs_poly.txt"
# res_cor_matrix_file <- "data/usda_envs_cor_matrix.txt"
# h2_file <- "data/1stStage_heritabilities.csv"
# means_file <- "data/1st_stage_means.txt"
# outfolder <- paste0("tests/sim_traits_gwas/", trait, "_per_env")
# real_pheno_file <- "data/1stStage_BLUEs.PHT-per-env.txt"
# effect_size_reduction <- 0.1
# seed <- 823162
# reps <- 3
# impute_type <- "Middle"
# architecture <- "pleiotropic"
# QTN_variance <- TRUE
# 
# minor_qtls <- FALSE
# reduce_effect_all <- FALSE
# outfolder <- paste0("tests/sim_traits_gwas/no_minor_qtls/", trait, "_per_env")
