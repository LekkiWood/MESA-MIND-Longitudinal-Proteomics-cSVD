#rm(list=ls())
#library(tidyverse)
#library(targets)
#library(glmnet)
#cleaned_proteins = tar_read(Proteins_long_clean)
#traits_db <- tar_read(traits_db)
#fa_factor_covariates = c("gender", "race", "site", "edu", "htnmeds", "smoking", "E4", "AFprevalent", "diabetes", "MIprevalent", "CHFprevalent")
#fa_numeric_covariates = c("age", "egfr", "BMI", "sbp", "ldl")
#factor_covariates = fa_factor_covariates
#numeric_covariates = fa_numeric_covariates
#outcome = "fa"
#chosen_exam = 1 
#MWASres = tar_read(FA_E1_PWAS)
#threshold = 0.05
#library(fastDummies)



cross_sectional_LASSO_function <- function(cleaned_proteins,
                                  traits_db, numeric_covariates, factor_covariates, outcome, chosen_exam, MWASres, threshold){
  
  
  
  #Make db
  
  #Select proteins
  
  results_of_interest <- paste(paste(outcome, chosen_exam, sep="_"), "fdr", sep="_")
  
  included_prots <- MWASres |>
    dplyr::filter(get(results_of_interest) < threshold) |>
    pull(OlinkID)
  
  rm(results_of_interest)
  
  vars_to_keep <- c(outcome, numeric_covariates, factor_covariates, "Exam", included_prots)
  rm(included_prots)
  
  common_keys <- intersect(names(cleaned_proteins), names(traits_db))
  
  
  dat <- traits_db |>
    dplyr::left_join(cleaned_proteins, by = common_keys) |>
    dplyr::select(dplyr::all_of(vars_to_keep)) |>
    dplyr::filter(Exam == chosen_exam) |> 
    droplevels() |>
    dplyr::select(-Exam) |>
    na.omit() |> #Remove missing data for dummy coding so no NA category 
    fastDummies::dummy_cols(remove_first_dummy = TRUE) |>
    dplyr::select(!all_of(factor_covariates))
  
  rm(vars_to_keep)

  
  
  # confirm required columns exist - doesn't work as dummy coded vars now... 
  # stopifnot(all(c(outcome, factor_covariates, numeric_covariates) %in% names(dat)))
  

  
  
  ########################################
  #Good to check when first writing LASSO
  ########################################
  
  #class <- as.data.frame(sapply(dat, class))
  #colnames(class) <- "class"
  #table(class$class) #one is an integer
  #rm(class)

  ########################################
  

  # Split into training and testing datasets 
  
  set.seed(11042012)
  #x <- model.matrix(outcome ~. , dat)[,-1]
  #Divide rows by 2 for 50/50 training:validation split
  sample <- sample(1:nrow(dat), nrow(dat)/2)
  
  train <- dat[sample, ]
  test <- dat[-sample, ]
  
  #Split into X and Y  from training and validation
  
  mdlY <- as.matrix(train[outcome])
  mdlX <- as.matrix(train[setdiff(colnames(dat), outcome)])
  newY <- as.matrix(test[outcome])
  newX <- as.matrix(test[setdiff(colnames(dat), outcome)])
  
  
  #Regularized model using elastic net
  
  registerDoParallel(cores=4)
  
  a <- seq(0.1, 0.9, 0.05)
  search <- foreach::foreach(i = a, .combine = rbind) %dopar% {
    cv <- cv.glmnet(mdlX, mdlY, family = "gaussian", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
    data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
  }
  cv3 <- search[search$cvm == min(search$cvm), ]
  md3 <- glmnet(mdlX, mdlY, family = "gaussian", lambda = cv3$lambda.1se, alpha = cv3$alpha)
  
  #Create list of vars selected by EN
  
  all_vars <- as.data.frame(as.matrix(coef(md3)))
  var_names <- rownames(all_vars)
  all_vars$feature <- var_names
  rownames(all_vars) <-c()
  colnames(all_vars) <- c("Value", "Feature")
  important_vars <- subset(all_vars, Value !=0) 
  #remove intercept
  n <- nrow(important_vars)
  important_vars <- important_vars[2:n,]
  
  list(
    LASSO_output_all_vars = all_vars,
    LASSO_output_important_vars = important_vars
    )

  
}