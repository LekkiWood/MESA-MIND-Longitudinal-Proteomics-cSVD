longitudinal_LASSO_function_BLUP <- function(
    cleaned_proteins,
    traits_db,
    numeric_covariates,
    factor_covariates,
    outcome,
    MWASres,
    threshold,
    exam_mri,
    min_exams_per_subject,
    nfolds,
    alpha_grid,
    cores
)
{
  
  # Packages assumed loaded by caller:
  # dplyr, tidyr, lme4, glmnet, fastDummies, foreach, doParallel
  
  #---------------------------
  # Select proteins
  #---------------------------
  results_of_interest <- paste0("P_time_", outcome, "_int_fdr")
  
  included_prots <- MWASres |>
    dplyr::filter(.data[[results_of_interest]] < threshold) |>
    dplyr::pull(protein) |>
    unique()
  
  if (length(included_prots) == 0) {
    stop("No proteins passed the MWAS threshold. Check 'threshold' and MWASres column: ", results_of_interest)
  }
  
  #---------------------------
  # Build long dataset for slope extraction
  #---------------------------
  dat <- cleaned_proteins |>
    dplyr::select(sidno, Exam, dplyr::all_of(included_prots)) |>
    dplyr::left_join(traits_db, dplyr::join_by(sidno, Exam)) |>
    dplyr::select(
      sidno, Exam, time,
      BL_age, egfr, icv,
      gender, race, edu, smoking, E4,
      dplyr::all_of(outcome),
      dplyr::all_of(included_prots)
    ) |>
    dplyr::group_by(sidno) |>
    dplyr::filter(dplyr::n_distinct(Exam) >= min_exams_per_subject) |>
    dplyr::ungroup()
  
  # Covariates to adjust in slope model:
  # (You can remove icv here if you prefer it only at MRI stage.)
  covs_for_slope <- c("BL_age", "egfr", "icv", "gender", "race", "edu", "smoking", "E4")
  
  #---------------------------
  # Helper: BLUP slope extraction for one protein
  #---------------------------
  get_blup_slopes_one_protein <- function(df, protein, covs) {
    
    # Scale time once per protein-model dataset
    df2 <- df |>
      dplyr::mutate(
        time_sc   = as.numeric(scale(time)),
        BL_age_sc = as.numeric(scale(BL_age)),
        egfr_sc   = as.numeric(scale(egfr)),
        icv_sc    = as.numeric(scale(icv))
      ) |>
      dplyr::select(
        sidno, time_sc,
        BL_age_sc, egfr_sc, icv_sc,
        gender, race, edu, smoking, E4,
        dplyr::all_of(protein)
      )
    
    model_vars <- c("sidno", "time_sc", "BL_age_sc", "egfr_sc", "icv_sc",
                    "gender", "race", "edu", "smoking", "E4", protein)
    
    df_m <- df2 |>
      dplyr::filter(stats::complete.cases(dplyr::across(dplyr::all_of(model_vars)))) |>
      droplevels()
    
    # Need at least 2 observations per person to estimate a slope;
    # with 3 timepoints expected, this mostly guards against heavy missingness.
    if (nrow(df_m) < 2 || dplyr::n_distinct(df_m$sidno) < 2) {
      out <- data.frame(sidno = unique(df$sidno), slope = NA_real_)
      names(out)[2] <- protein
      return(out)
    }
    
    # Mixed model with random intercept + random slope (no correlation term)
    fml <- stats::as.formula(
      paste0(
        protein,
        " ~ time_sc + BL_age_sc + egfr_sc + icv_sc + gender + race + edu + smoking + E4",
        " + (1 | sidno) + (0 + time_sc | sidno)"
      )
    )
    
    fit <- lme4::lmer(fml, data = df_m, REML = TRUE)
    
    # BLUP subject-specific slope: fixed slope + random slope deviation
    re <- lme4::ranef(fit)$sidno
    if (!("time_sc" %in% colnames(re))) {
      # Shouldn't happen, but keep safe
      out <- data.frame(sidno = unique(df$sidno), slope = NA_real_)
      names(out)[2] <- protein
      return(out)
    }
    
    subj_slope <- lme4::fixef(fit)[["time_sc"]] + re[, "time_sc"]
    
    out <- data.frame(sidno = rownames(re), slope = as.numeric(subj_slope))
    names(out)[2] <- protein
    out
  }
  
  #---------------------------
  # Extract BLUP slopes for all proteins
  #---------------------------
  slopes_list <- lapply(included_prots, function(p) get_blup_slopes_one_protein(dat, p, covs_for_slope))
  slopes_df <- Reduce(function(x, y) merge(x, y, by = "sidno", all = TRUE), slopes_list)
  
  #---------------------------
  # Build MRI-stage dataset (one row per subject at Exam 6)
  #---------------------------
  slopes_df_LASSO <- traits_db |>
    dplyr::filter(Exam == exam_mri) |>
    dplyr::select(sidno, dplyr::all_of(outcome), dplyr::all_of(numeric_covariates), dplyr::all_of(factor_covariates)) |>
    dplyr::left_join(slopes_df, dplyr::join_by(sidno))
  
  # Drop time at MRI if it sneaks in via numeric_covariates and you don't want it:
  if ("time" %in% names(slopes_df_LASSO)) {
    slopes_df_LASSO <- slopes_df_LASSO |> dplyr::select(-time)
  }
  
  # Keep complete cases for now (as you noted, this can shrink N a lot)
  slopes_df_LASSO_cc <- slopes_df_LASSO |>
    dplyr::select(-sidno) |>
    tidyr::drop_na()
  
  # Diagnostics
  diag <- list(
    N_after_complete_case = nrow(slopes_df_LASSO_cc),
    P_total_columns = ncol(slopes_df_LASSO_cc),
    outcome_sd = stats::sd(slopes_df_LASSO_cc[[outcome]]),
    any_na = anyNA(slopes_df_LASSO_cc)
  )
  
  if (diag$N_after_complete_case < 20) {
    warning("Very small N after complete-case filtering: ", diag$N_after_complete_case,
            ". LASSO/EN may select nothing. Consider imputing slopes or filtering proteins by missingness.")
  }
  
  #---------------------------
  # Encode factors + build X/Y
  #---------------------------
  set.seed(11042012)
  #x <- model.matrix(outcome ~. , dat)[,-1]
  #Divide rows by 2 for 50/50 training:validation split
  sample <- sample(1:nrow(slopes_df_LASSO_cc), nrow(slopes_df_LASSO_cc)/2)
  
  LASSO_N = nrow(slopes_df_LASSO_cc)
  Outcome_SD = sd(slopes_df_LASSO_cc[[outcome]])
  N_NA =  anyNA(slopes_df_LASSO_cc)
  
  
  train <- slopes_df_LASSO_cc[sample, ]
  test <- slopes_df_LASSO_cc[-sample, ]
  
  #Split into X and Y  from training and validation
  
  mdlY <- as.matrix(train[outcome])
  mdlX <- as.matrix(train[setdiff(colnames(slopes_df_LASSO_cc), outcome)])
  newY <- as.matrix(test[outcome])
  newX <- as.matrix(test[setdiff(colnames(slopes_df_LASSO_cc), outcome)])
  
  
  #Regularized model using elastic net
  
  registerDoParallel(cores=4)
  
  a <- seq(0.1, 0.9, 0.05)
  search <- foreach::foreach(i = a, .combine = rbind) %dopar% {
    cv <- cv.glmnet(mdlX, mdlY, family = "gaussian", nfold = 10, type.measure = "deviance", parallel = TRUE, alpha = i)
    #data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
    # }
    data.frame(
      cvm_min = min(cv$cvm),
      lambda.min = cv$lambda.min,
      lambda.1se = cv$lambda.1se,
      alpha = i
    )
  }
  
  cv3 <- search[which.min(search$cvm_min), ]
  #md3 <- glmnet(mdlX, mdlY, family="gaussian", lambda=cv3$lambda.min, alpha=cv3$alpha)
  cv_final <- cv.glmnet(mdlX, mdlY, alpha=cv3$alpha)
  coef_min <- coef(cv_final, s="lambda.min")
  coef_1se <- coef(cv_final, s="lambda.1se")
  
  
  # cv3 <- search[search$cvm == min(search$cvm), ]
  #  md3 <- glmnet(mdlX, mdlY, family = "gaussian", lambda = cv3$lambda.1se, alpha = cv3$alpha)
  
  #Create list of vars selected by EN
  
  #####Stringent lambda
  all_vars_coef_1se <- as.data.frame(as.matrix(coef_1se))
  var_names <- rownames(all_vars_coef_1se)
  all_vars_coef_1se$feature <- var_names
  rownames(all_vars_coef_1se) <-c()
  colnames(all_vars_coef_1se) <- c("Value", "Feature")
  important_vars_coef_1se <- subset(all_vars_coef_1se, Value !=0) 
  #remove intercept
  important_vars_coef_1se <- subset(important_vars_coef_1se, Feature !="(Intercept)")
  
  #####Relaxed lambda
  all_vars_coef_min <- as.data.frame(as.matrix(coef_min))
  var_names <- rownames(all_vars_coef_min)
  all_vars_coef_min$feature <- var_names
  rownames(all_vars_coef_min) <-c()
  colnames(all_vars_coef_min) <- c("Value", "Feature")
  important_vars_coef_min <- subset(all_vars_coef_min, Value !=0) 
  #remove intercept
  important_vars_coef_min <- subset(important_vars_coef_min, Feature !="(Intercept)")
  
  #---------------------------
  # Return
  #---------------------------
  list(
    diagnostics = list(LASSO_N = LASSO_N,
                       Outcome_SD = Outcome_SD,
                       N_NA = N_NA),
    data = list(
      LASSO_output_all_vars_coef_min = all_vars_coef_min,
      LASSO_output_important_vars_coef_min = important_vars_coef_min,
      LASSO_output_all_vars_coef_1se = all_vars_coef_1se,
      LASSO_output_important_vars_coef_1se = important_vars_coef_1se)
  )
  
  
}