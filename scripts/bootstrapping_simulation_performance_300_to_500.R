# Libraries
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
library(tidyverse)

# Set working directory
setwd('~/empirical_sim')

# Load environment
#load('trait_simulation_n400_issue.RData')

# Load files and parameters
geno_file <- "data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.hmp.txt"
sv_list <- "data/SVs_IDs_poly.txt"
res_cor_matrix_file <- "data/usda_envs_cor_matrix.txt"
h2_file <- "data/1stStage_heritabilities.csv"
means_file <- "data/1st_stage_means.txt"
outfolder <- 'test'
real_pheno_file <- "data/1stStage_BLUEs.YLD-per-env.txt"
effect_size_reduction <- 1
seed <- 823162
reps <- 3
impute_type <- "Middle"
architecture <- "pleiotropic"
QTN_variance <- TRUE

minor_qtls <- TRUE
reduce_effect_all <- FALSE

# load genotypic data
geno_data <- fread(geno_file, header = TRUE, data.table = FALSE)

# Create storage table
data_storage = tibble()

# Loop through simulation code
for(trait in c('EHT', 'PHT', 'YLD', 'Moisture')){
  for(n_marker in c(3,5,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000)){
  # Show iteration progress
    print(paste0('--- ', trait, ' ', n_marker, ' ---'))
  gwas_file <- paste0("analysis/sim_traits/", trait, "/gwas_markers/gwas_top-", n_marker, "_markers.txt")
  #gwas_hits_500 = fread(gwas_file, header = TRUE, data.table = FALSE)
  gwas_hits = fread(gwas_file, header = TRUE, data.table = FALSE)
  #for(effect_size in c(0.1, 1, 10)){
    
    #for(i in sort(c(seq(397,403,1), seq(350, 450, 10)))){
    #  print(paste0('--- ', trait, ' ', effect_size, ' ' , i, ' ---'))
    #  gwas_hits = head(gwas_hits_500, i*length(unique(gwas_hits_500$env))*2)
    #  gwas_hits$marker_type <- sapply(gwas_hits$marker, function(marker) {
    #    type <- ifelse(grepl("^del|^dup|^inv|^ins|^tra", marker, perl = TRUE),
    #                   yes = "SV", no = "SNP")
    #    return(type)
    #  })
      #gwas_hits[,"effect"] = effect_size
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
          add_effect_value[[env]] <- c(unlist(add_effect_value_major[[env]]), add_effect_value_minor[[env]])
          dom_effect_value[[env]] <- c(unlist(dom_effect_value_major[[env]]), dom_effect_value_minor[[env]])
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
      
      # filter to have only markers from selected qtl list
      geno_data_sub <- subset(geno_data, geno_data[, 1] %in% as.character(unlist(QTN_list)))
      
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
      sim_pheno <- create_phenotypes(geno_obj = geno_data_sub,
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
                                     output_dir = outfolder,
                                     to_r = TRUE,
                                     out_geno = NULL,
                                     QTN_variance = QTN_variance,
                                     quiet = TRUE,
                                     # numericalization
                                     SNP_effect = impute_effect,
                                     SNP_impute = impute_type)
      
      
      
      #### plot sim vs real traits ----
    
      # format sim data
      colnames(sim_pheno)[1] <- "hybrid"
      colnames(sim_pheno)[2:(length(envs_names) + 1)] <- envs_names
      sim_pheno <- pivot_longer(sim_pheno, cols = all_of(envs_names), names_to = "env", values_to = "sim_value")
      
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
      
      data_storage = data_storage %>%
        bind_rows(tibble(n_markers = n_marker,
                         correlation = cor(sim_real_phenos$sim_value_avg, sim_real_phenos$real_pheno, use = 'complete.obs'),
                         #effect_size = effect_size,
                         trait = trait))
   # }
  }
}

#data_storage

#data_storage %>%
#  write_csv('performance_drop_across_traits_fake_effects.csv')

plot = data_storage %>%
  ggplot(aes(x = n_markers, y = correlation))+
  geom_point()+
  geom_line()+
  #labs(color = 'Effect Size\nof all\nMarkers')+
  facet_wrap(~trait)+
  theme_classic()

#ggsave('performance_drop_across_traits_fake_effects.png', plot, device = 'png', width = 3.75, height = 3.5, units = 'in')

#write_delim(geno_data, 'geno_data.hmp.txt', delim = '\t')
#write_delim(gwas_hits_400, 'gwas_hits_400.txt', delim = '\t')
#save(res_cor_matrix, file = 'environmental_similarity_matrix.Rdata')
