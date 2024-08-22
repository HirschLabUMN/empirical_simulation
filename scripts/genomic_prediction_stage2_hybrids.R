library(data.table)
library(asreml)
library(cvTools)
library(plyr)
library(rrBLUP)
library(MASS)
library(reshape)
library(tidyr)
library(gtools)
#source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("~/michael/empirical_sim/scripts/GAPIT_3.4_Source_Code.txt")
suppressWarnings(suppressMessages(library(doParallel)))
library(AGHmatrix)


usage <- function() {
  cat("
description: simulate trait based on user-defined genetic architecture.

usage: Rscript genomic_prediction_stage2_hybrids.R [markers_file] [blues_file] [output_folder] [...]

positional arguments:
  markers_file            hapmap file with genotypic data
  blues_file              file with BLUEs in long format
  output_folder           output folder name

optional argument:
  --help                    show this helpful message
  --cv-type=VALUE           choose 'CV2 or 'CV1' cross-validation scheme (default: CV2)
  --n-folds=VALUE           number of folds to be used in cross validation
  --cv-iter=VALUE           number of cross validataion iterations
  --total-envs=VALUE        choose how many environments will be used (default: all)
  --impute-effect=VALUE     the marker effect ('Add' or 'Dom') when imputing missing data at the hapmap
                            numericalization step (default: 'Dom')
  --impute-type=VALUE       the marker type ('Major', 'Middle' or 'Minor') when imputing missing data at
                            the hapmap numericalization step (default: 'Middle')
  --seed=VALUE              value for set.seed (default: NULL; random number is selected)
  --model-effect=VALUE      choose if prediction models will use only additive ('A'), only dominance ('D'),
                            or both additive and dominance ('AD') effects (default: 'A')
  --model-Gstr=VALUE        type of variance-covariance structure to use with random effects. Valid options are
                            'idv' (homogenous, default), 'diag' (heterogenous), 'fa' (factor analytic) or 'us'
                            (unstructured)
  --model-Rstr=VALUE        type of variance-covariance structure to use with residual effects. Valid options are
                            'idv' (homogenous, default), 'diag' (heterogenous), or 'us' (unstructured)
  --model-Gstr-init=VALUE   type of random variance-covariance structure to be obtain initial variance values. Take
                            same values as '--model-Gstr'
  --model-Rstr-init=VALUE   type of residual variance-covariance structure to be obtain initial variance values. Take
                            same values as '--model-Rstr'
  --envs-weight             incorporate weights (inverse of the variance of the predictions
                            from first-stage analysis) into the multi-environment model

credits:
  part of this script was modified from Fernandes et al. 2018
  (TAG, doi:10.1007/s00122-017-3033-y), available at
  https://github.com/samuelbfernandes/Trait-assisted-GS
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
markers_file <- args[1]
blues_file <- args[2]
output_folder <- args[3]

# set default
cv_type <- "CV2"
n_folds <- "5"
cv_iter <- "3"
total_envs <- NULL
impute_effect <- "Dom"
impute_type <- "Middle"
model_effect <- "A"
model_Gstr <- "idv"
model_Rstr <- "idv"
model_Gstr_init <- "idv"
model_Rstr_init <- "idv"
seed <- NULL
envs_weight <- FALSE

# assert to have the correct optional arguments
if (length(args) < 3) stop(usage(), "missing positional argument(s)")

if (length(args) > 3) {

  opt_args <- args[-1:-3]
  opt_args_allowed <- c("--cv-type", "--n-folds", "--cv-iter", "--total-envs","--impute-effect", "--impute-type",
                        "--model-effect", "--model-Gstr", "--model-Rstr", "--model-Gstr-init", "--model-Rstr-init",
                        "--seed", "--envs-weight")
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
if (!cv_type %in% c("CV2", "CV1")) {
  stop("Optional argument '--cv-type' should be either 'CV2' or 'CV1'")
}

if (suppressWarnings(!is.na(as.numeric(n_folds)))) {
  n_folds <- as.numeric(n_folds)
} else {
  stop("Optional argument '--n-folds' should be a number")
}

if (suppressWarnings(!is.na(as.numeric(cv_iter)))) {
  cv_iter <- as.numeric(cv_iter)
} else {
  stop("Optional argument '--cv-iter' should be a number")
}

if (!is.null(total_envs)) {
  if (suppressWarnings(!is.na(as.numeric(total_envs)))) {
    total_envs <- as.numeric(total_envs)
  } else {
    stop("Optional argument '--total-envs' should be a number")
  }
}

if (!impute_effect %in% c("Add", "Dom")) {
  stop("Optional argument '--impute-effect' should be either 'Add' or 'Dom'")
}

if (!impute_type %in% c("Major", "Middle", "Minor")) {
  stop("Optional argument '--impute-type' should be either 'Major', 'Middle' or 'Minor'")
}

if (!model_effect %in% c("A", "D", "AD")) {
  stop("Optional argument '--model-effect' should be either 'A', 'D' or 'AD'")
}

if (!model_Gstr %in% c("idv", "diag", "fa", "us")) {
  stop("Optional argument '--model-Gstr' should be either 'idv', 'diag', 'fa' or 'us'")
}

if (!model_Rstr %in% c("idv", "diag", "us")) {
  stop("Optional argument '--model-Rstr' should be either 'idv', 'diag' or 'us'")
}

if (!model_Gstr_init %in% c("idv", "diag", "fa", "us")) {
  stop("Optional argument '--model-Gstr-init' should be either 'idv', 'diag', 'fa' or 'us'")
}

if (!model_Rstr_init %in% c("idv", "diag", "us")) {
  stop("Optional argument '--model-Rstr-init' should be either 'idv', 'diag' or 'us'")
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



#### load data ----

# "header = FALSE" for compatibility with GAPIT
hapmap <- fread(markers_file, header = FALSE, data.table = FALSE)
# numericalize hapmap
hapmap <- GAPIT.HapMap(G = hapmap, SNP.effect = impute_effect, SNP.impute = impute_type)
#   hm$GT: vector with line names
#   hm$GD: matrix where each row is a line and each column is a marker (the numbers of the cells are the numeric genotypes)
#   hm$GI: data frame with marker information (name, chromosome and position)

# length(hapmap$GI$SNP) == NCOL(markers)
# length(hapmap$GT) == NROW(markers)
markers <- data.frame(hapmap$GD - 1)
colnames(markers) <- hapmap$GI$SNP
rownames(markers) <- hapmap$GT
rm(hapmap)

# remove duplicated markers to avoid singularities on relationship matrix
markers <- markers[!duplicated(markers), ]

# load and format phenotypic data
means <- fread(blues_file, header = TRUE, data.table = FALSE)

if (envs_weight) {
  # create a separate data frame with environmental weights
  weights <- aggregate(weight ~ environment, FUN = mean, data = means)
  weights <- weights[mixedorder(weights$environment), ]
  # remove weight column for compatibility
  means <- means[, colnames(means) != "weight"]
}

# transform data to wide format
means <- pivot_wider(means, names_from = environment, values_from = predicted_trait_value)
means <- data.frame(means)

# filter total number of environments
if (is.null(total_envs)) total_envs <- NCOL(means) - 1
means <- means[, c(1:(total_envs + 1))]
if (envs_weight) weights <- weights[1:total_envs, ]

# using only phenotypes with snp information
geno_to_keep <- intersect(means$genotype, rownames(markers))
markers <- markers[geno_to_keep, ]
means <- means[means$genotype %in% geno_to_keep, ]
means <- means[match(rownames(markers), means$genotype), ]
if (!all(rownames(markers) == means$genotype)) stop("genotype names don't match between marker and pheno data")
rm(geno_to_keep)

# changing names for numbers, but save ids for later as well
geno_names <- means$genotype
geno_ids <- 1:length(geno_names)
means$genotype <- geno_ids
means <- transform(means, genotype = factor(genotype))



#### create relationship matrix ----

# additive relationship matrix
kin_add <- Gmatrix(SNPmatrix = as.matrix(markers) + 1, method = "VanRaden", verify.posdef = TRUE)  # NOT positive definite
kin_add <- kin_add + diag(0.001, ncol(kin_add)) # add noise to the diag to reduce saverage information matrix signularities
ginv_A <- ginv(kin_add)
ginv_A <- data.frame(formatmatrix(ginv_A, round.by = 4, exclude.0 = FALSE, return = TRUE, save = FALSE))
colnames(ginv_A)<- c("Row", "Column", "GINV")
attr(ginv_A, "rowNames") <- geno_ids

# dominance relationship matrix
kin_dom <- Gmatrix(SNPmatrix = as.matrix(markers) + 1, method = "Vitezica", verify.posdef = TRUE)  # positive definite
kin_dom <- kin_dom + diag(0.001, ncol(kin_dom)) # add noise to the diag to reduce saverage information matrix signularities
ginv_D <- ginv(kin_dom)
ginv_D <- data.frame(formatmatrix(ginv_D, round.by = 4, exclude.0 = FALSE, return = TRUE, save = FALSE))
colnames(ginv_D)<- c("Row", "Column", "GINV")
attr(ginv_D, "rowNames") <- geno_ids
rm(markers, kin_add, kin_dom)



#### create groups for 5-fold CV ----

# get names of envs
envs_names <- colnames(means)[-1]

# i'll do a few extra CV iterations than what user supplied in the command line
# to account for possible problems during model convergence by asreml

k_folds <- list()
for(i in 1:(cv_iter * 3)){

  # set seed for reproducibility
  set.seed(seed + i)

  # divide genotype into folds according to CV scheme
  if (cv_type == "CV1") {
    folds <- cvFolds(nlevels(means$genotype), type = "random", K = n_folds, R = 1)
    sample <- cbind(folds$which, folds$subsets)
    cv <- split(sample[, 2], f = sample[, 1])
  }

  if (cv_type == "CV2") {
    folds <- cvFolds(nlevels(means$genotype), type = "random", K = n_folds, R = total_envs)
    sample <- cbind(folds$which, folds$subsets)
    # different environments in the same fold should have different genotypes selected
    cv <- list()
    for (env in 1:total_envs) {
      cv_env <- split(sample[, 1 + env], f = sample[, 1])
      names(cv_env) <- paste0("fold", names(cv_env), "_", envs_names[env])
      cv <- append(cv, cv_env)
    }
    rm(cv_env, env)
  }

  # add folds to iteration
  k_folds[[i]] <- cv

}
rm(i, folds, cv, sample)



#### genomic prediction ----

# format inverted A and D matrix for asreml compatibility
attr(ginv_A, "INVERSE") <- TRUE
attr(ginv_D, "INVERSE") <- TRUE
attr(ginv_A, "rowNames") <- as.character(attr(ginv_A, "rowNames"))
attr(ginv_D, "rowNames") <- as.character(attr(ginv_D, "rowNames"))

# start list to store results
accuracy <- list()
gebv <- geno_ids

# run k-fold cross-validation if CV1 or CV2
for(j in 1:cv_iter) {

  cat(paste0("\n---- ", n_folds, "-fold ", cv_type," iteration ", j, " ---\n"))

  extra_cv_iter <- 0
  models_converged <- FALSE
  while (models_converged == FALSE) {

    registerDoParallel(cores = n_folds)
    results_folds <- try(foreach(i = 1:n_folds, .combine = c) %dopar% {

      cat("  fold ", i, ":\n", sep = "")

      cat("    preparing data\n")

      # get pheno data
      means_test <- means
      # mask genotypes in test population according to CV scheme
      if (cv_type == "CV1") {
        means_test[which(means[, "genotype"] %in% k_folds[[j+extra_cv_iter]][[i]]), 2:NCOL(means)] <- NA
      }
      if (cv_type == "CV2") {
        for (env in 1:total_envs) {
          means_test[which(means[, "genotype"] %in% k_folds[[j+extra_cv_iter]][[paste0("fold", i, "_", envs_names[env])]]), 1 + env] <- NA
        }
      }
      # transform data to long format
      means_test <- pivot_longer(means_test, -genotype, names_to = "environment", values_to = "sim_trait")
      means_test <- means_test[order(means_test$genotype, means_test$environment), ]
      means_test$environment <- as.factor(means_test$environment)
      if (envs_weight) means_test$weight <- weights[match(means_test$environment, weights$environment), "weight"]

      cat("    obtaining initial variance values\n")

      # define asreml random variance-covariance structures
      Gstr <- paste0(model_Gstr_init, "(environment)")
      if (model_Gstr_init == "fa") Gstr <- paste0(model_Gstr_init, "(environment, 1)")
      # define asreml residual variance-covariance structures
      Rstr <- paste0(model_Rstr_init, "(environment)")
      # define asreml random terms
      if (model_effect == "A") random_terms <- paste0("~ ", Gstr, ":vm(genotype, source = ginv_A)")
      if (model_effect == "D") random_terms <- paste0("~ ", Gstr, ":vm(genotype, source = ginv_D)")
      if (model_effect == "AD") random_terms <- paste0("~ ", Gstr, ":vm(genotype, source = ginv_A) + ", Gstr, ":vm(genotype, source = ginv_D)")
      # define asreml residual terms
      residual_terms <- paste0("~ id(units):", Rstr)
      # add or not env weights into asreml
      if (envs_weight) {
        envs_weight_terms <- "weights = weight, family = asr_gaussian(dispersion = 1),"
      } else {
        envs_weight_terms <- ""
      }
      # set up model
      model_init_cmd <- paste0('asreml(fixed = sim_trait ~ environment,
                                       random = ', random_terms, ',
                                       residual = ', residual_terms, ',
                                       asmv = environment, ',
                                       envs_weight_terms,
                                       'workspace = "1000mb",
                                       maxit = 100,
                                       trace = FALSE,
                                       data = means_test)')
      # run initial model
      model_init <- eval(parse(text = model_init_cmd))
      # get initial values
      if (model_effect == "A") {
        G_init <- data.frame(Component = c(names(model_init$G.param$`environment:vm(genotype, source = ginv_A)`$environment$initial)),
                             Value = c(as.numeric(model_init$G.param$`environment:vm(genotype, source = ginv_A)`$environment$initial)),
                             Constraint = c(model_init$G.param$`environment:vm(genotype, source = ginv_A)`$environment$con))
        R_init <- data.frame(Component = names(model_init$R.param$`units:environment`$environment$initial),
                             Value = as.numeric(model_init$R.param$`units:environment`$environment$initial),
                             Constraint = model_init$R.param$`units:environment`$environment$con)
      }
      if (model_effect == "D") {
        G_init <- data.frame(Component = c(names(model_init$G.param$`environment:vm(genotype, source = ginv_D)`$environment$initial)),
                             Value = c(as.numeric(model_init$G.param$`environment:vm(genotype, source = ginv_D)`$environment$initial)),
                             Constraint = c(model_init$G.param$`environment:vm(genotype, source = ginv_D)`$environment$con))
        R_init <- data.frame(Component = names(model_init$R.param$`units:environment`$environment$initial),
                             Value = as.numeric(model_init$R.param$`units:environment`$environment$initial),
                             Constraint = model_init$R.param$`units:environment`$environment$con)
      }
      if (model_effect == "AD") {
        G_init <- data.frame(Component = c(names(model_init$G.param$`environment:vm(genotype, source = ginv_A)`$environment$initial),
                                           names(model_init$G.param$`environment:vm(genotype, source = ginv_D)`$environment$initial)),
                             Value = c(as.numeric(model_init$G.param$`environment:vm(genotype, source = ginv_A)`$environment$initial),
                                       as.numeric(model_init$G.param$`environment:vm(genotype, source = ginv_D)`$environment$initial)),
                             Constraint = c(model_init$G.param$`environment:vm(genotype, source = ginv_A)`$environment$con,
                                            model_init$G.param$`environment:vm(genotype, source = ginv_D)`$environment$con))
        R_init <- data.frame(Component = names(model_init$R.param$`units:environment`$environment$initial),
                             Value = as.numeric(model_init$R.param$`units:environment`$environment$initial),
                             Constraint = model_init$R.param$`units:environment`$environment$con)
      }
      rm(model_init_cmd, model_init)

      cat("    fitting model\n")

      # define asreml random variance-covariance structures
      Gstr <- paste0(model_Gstr, "(environment)")
      if (model_Gstr == "fa") Gstr <- paste0(model_Gstr, "(environment, 1)")
      # define asreml residual variance-covariance structures
      Rstr <- paste0(model_Rstr, "(environment)")
      # define asreml random terms
      if (model_effect == "A") random_terms <- paste0("~ ", Gstr, ":vm(genotype, source = ginv_A)")
      if (model_effect == "D") random_terms <- paste0("~ ", Gstr, ":vm(genotype, source = ginv_D)")
      if (model_effect == "AD") random_terms <- paste0("~ ", Gstr, ":vm(genotype, source = ginv_A) + ", Gstr, ":vm(genotype, source = ginv_D)")
      # define asreml residual terms
      residual_terms <- paste0("~ id(units):", Rstr)
      # add or not env weights into asreml
      if (envs_weight) {
        envs_weight_terms <- "weights = weight, family = asr_gaussian(dispersion = 1),"
      } else {
        envs_weight_terms <- ""
      }
      # set up model
      model_cmd <- paste0('asreml(fixed = sim_trait ~ environment,
                                  random = ', random_terms, ',
                                  residual = ', residual_terms, ',
                                  asmv = environment, ',
                                  envs_weight_terms,
                                  'workspace = "1000mb",
                                  maxit = 100,
                                  trace = FALSE,
                                  data = means_test,
                                  G.param = G_init,
                                  R.param = R_init)')
      # run more complex model with initial values
      model <- eval(parse(text = model_cmd))

      # update model until converge -- do it max 5 times
      try <- 1
      while (!model$converge & try <= 5) {
        cat("    updating model until converge...\n")
        model <- update.asreml(model)
        try <- try + 1
      }
      if (!model$converge & try == 6) stop("Model didn't converge")


      cat("    running predictions...\n")

      # predict phenotypes
      pred <- predict(model, classify = "environment:genotype", pworkspace = "1000mb", trace = FALSE)
      pred <- data.frame(pred$pvals)
      pred <- pred[, c("genotype", "environment", "predicted.value")]
      pred <- pivot_wider(pred, names_from = environment, values_from = predicted.value)

      # return predicted phenotypes for that fold according to CV scheme
      if (cv_type == "CV1") {
        if (all(k_folds[[j+extra_cv_iter]][[i]] == pred[k_folds[[j+extra_cv_iter]][[i]], "genotype"])) {
          # cv1 table will be in wide format
          pred_fold <- pred[k_folds[[j+extra_cv_iter]][[i]], ]
        } else {
          stop("genotypes don't match")
        }
      }
      if (cv_type == "CV2") {
        pred_fold <- data.frame()
        for (env in 1:total_envs) {
          fold_env <- names(k_folds[[j+extra_cv_iter]])[grep(paste0("fold", i, "_", envs_names[env], "$"), names(k_folds[[j+extra_cv_iter]]), perl = TRUE)]
          fold_values <- pred[k_folds[[j+extra_cv_iter]][[fold_env]], c("genotype", envs_names[env])]
          colnames(fold_values) <- c("genotype", "pred_values")
          # cv2 table will be in long format to keep track of genotypes masked in each env
          pred_fold <- rbind(pred_fold, cbind(fold = paste0("fold", i), environment = envs_names[env], fold_values))
        }
      }

      cat("    done!\n\n")

      return(list(pred_fold))

    })
    stopImplicitCluster()

    # if any of the folds had singularity/convergence issues, try again with
    # an another set of folds (i.e. another CV iteration)
    if (class(results_folds) == "try-error") {

      models_converged <- FALSE
      extra_cv_iter <- extra_cv_iter + cv_iter
      if (extra_cv_iter >= cv_iter * 3) stop("ASREML couldn't converge models after 3 extra CV iterations")

    } else {

      # if there was no issue, proceed to calculate accuracy across folds
      models_converged <- TRUE

    }
  }

  # combine predicted phenotypes from all folds
  pred_pheno <- do.call(rbind, results_folds)
  if (cv_type == "CV1") {
    # sort by genotype
    pred_pheno <- pred_pheno[order(pred_pheno$genotype), ]
  }
  if (cv_type == "CV2") {
    pred_pheno_merged <- data.frame()
    # combine one environment at a time
    for (env in 1:total_envs) {
      pred_pheno_env <- subset(pred_pheno, environment == envs_names[env])
      # sort by genotype
      pred_pheno_env <- pred_pheno_env[order(pred_pheno_env$genotype), c("genotype", "pred_values")]
      colnames(pred_pheno_env) <- c("genotype", envs_names[env])
      if (env == 1) {
        # keep genotype column
        pred_pheno_merged <- pred_pheno_env
      } else {
        # skip genotype column
        pred_pheno_merged <- cbind(pred_pheno_merged, pred_pheno_env[, 2, drop = FALSE])
      }
    }
    # overwrite original table for compatibility with the rest of the script
    pred_pheno <- pred_pheno_merged
    rm(pred_pheno_merged, pred_pheno_env)
  }

  # calculate prediction accuracy
  accuracy[[j]] <- list()
  for (env in colnames(means)[2:NCOL(means)]) {
    accuracy[[j]] <- append(accuracy[[j]], cor(pred_pheno[, env], means[, env]))
  }
  accuracy[[j]] <- unlist(accuracy[[j]])

  # also get the prediction values in each iteration
  colnames(pred_pheno)[-1] <- paste0(colnames(pred_pheno)[-1], "_iter", j)
  if (all(geno_ids == pred_pheno[, "genotype"])) {
    gebv <- cbind(gebv, pred_pheno[, 2:NCOL(pred_pheno)])
  } else {
    stop("genotypes don't match")
  }

}

# transform into a data frame
accuracy <- do.call(rbind, accuracy)

# calculate summary stats after many iterations of CV
summary <- data.frame(apply(accuracy, MARGIN = 2, function(x) {
  # mean
  mean <- mean(x, na.rm = TRUE)
  # std error
  se <- sd(x, na.rm = TRUE)/sqrt(length(x))
  # 95% confidence interval
  t_value <- qt(p = 0.05/2, df = NROW(accuracy) - 1, lower.tail = FALSE)
  lower_CI <- mean - (t_value * se)
  upper_CI <- mean + (t_value * se)
  return(c(mean = mean, se = se, lower_CI = lower_CI, upper_CI = upper_CI))
}))
colnames(summary) <- colnames(means)[2:NCOL(means)]




#### write results ----

# make sure folder exists
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# prediction accuracy
summary <- cbind(stat = rownames(summary), summary)
rownames(summary) <- NULL
fwrite(summary, file = paste0(output_folder, "/prediction_accuracy.", cv_type, ".txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# GEBVs
colnames(gebv)[1] <- "genotype"
if (any(geno_ids != gebv$genotype)) gebv[intersect(geno_ids, gebv$genotype), ]
gebv$genotype <- geno_names
fwrite(gebv, file = paste0(output_folder, "/GEBVs.", cv_type, ".txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# get mean predicted performance of each genotype at each environment
gebv_top <- pivot_longer(data = gebv, cols = -genotype, names_to = c("env", "iter"), names_sep = "_", values_to = "pred_pheno")
gebv_top <- aggregate(pred_pheno ~ genotype + env, FUN = mean, data = gebv_top)

gebv_ranks <- gebv$genotype
for (i in 1:total_envs) {

  # sort by best to worst performance
  gebv_top_env <- subset(gebv_top, env == paste0("env", i))
  gebv_top_env <- gebv_top_env[order(gebv_top_env$pred_pheno, decreasing = TRUE), c("genotype", "pred_pheno")]
  rownames(gebv_top_env) <- NULL
  fwrite(gebv_top_env, file = paste0(output_folder, "/GEBVs_best-env", i, ".", cv_type, ".txt"),
         quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

  # also save ranks of genotypes by environments
  gebv_ranks <- cbind(gebv_ranks, match(gebv$genotype, gebv_top_env$genotype))
  colnames(gebv_ranks)[NCOL(gebv_ranks)] <- paste0("env", i)

}
colnames(gebv_ranks)[1] <- "genotype"
fwrite(gebv_ranks, file = paste0(output_folder, "/GEBVs_ranks.", cv_type, ".txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)



# #### debug ----
# setwd('~/Desktop/Grad_School/Research/empirical_sim/')
# # markers_file <- "analysis/test_prediction/linux/usda_hybrids.all_markers.adjusted-n-markers.iter1.hmp.txt"
# markers_file <- "analysis/sim_traits/predictors/iter1/hybrids_geno.iter1.hmp.txt"
# # blues_file <- "analysis/test_prediction/linux/blups_1st_stage_hybrids.10-add-snps_0.3-h2.pop1.txt"
# blues_file <- "analysis/sim_traits/YLD/traits/michael_method/n_markers_3/effects_1/blues_1st_stage.txt"
# # output_folder <- "analysis/test_prediction/linux/prediction_all_markers_hybrids"
# output_folder <- "analysis/sim_traits/YLD/traits/michael_method/n_markers_3/effects_1/prediction_results/iter1/model-A/Gstr-diag_Rstr-us"
# n_folds <- 5
# cv_iter <- 3
# total_envs <- NULL
# seed <- 2021
# cv_type <- "CV1"
# envs_weight <- FALSE
# # envs_weight <- TRUE
# impute_effect <- "Add"
# impute_type <- "Middle"
# model_effect <- "A"
# model_Gstr <- "fa" # Was idv
# model_Rstr <- "diag"
# model_Gstr_init <- "diag" # Was idv
# model_Rstr_init <- "diag" # Was idv
