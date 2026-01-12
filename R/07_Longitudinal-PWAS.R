# Proteins are outcomes (LHS), predictor (e.g., "fa") is RHS.
# Scales BOTH the protein outcome and the predictor.
longitudinal_PWAS_function <- function(cleaned_proteins, protein_mapping,
                                              traits_db, numeric_covariates,
                                              factor_covariates, predictor) {
  
  common_keys <- intersect(names(cleaned_proteins), names(traits_db))
  
  PWAS_file <- traits_db |>
    dplyr::left_join(cleaned_proteins, by = common_keys)
  
  proteins <- protein_mapping |>
    dplyr::filter(OlinkID %in% names(PWAS_file)) |>
    dplyr::pull(OlinkID)
  
  stopifnot(all(c(predictor, factor_covariates, numeric_covariates) %in% names(PWAS_file)))
  
  # keep subjects with all exams (you were using >= 3)
  dat <- PWAS_file |>
    dplyr::group_by(idno) |>
    dplyr::filter(n_distinct(Exam) >= 3) |>
    dplyr::ungroup()
  
  out <- data.frame(
    Protein = proteins,
    Nobs = NA_integer_,
    Model = NA_character_,
    Beta_fa = NA_real_,
    SE_fa = NA_real_,
    P_fa = NA_real_,
    Beta_time = NA_real_,
    SE_time = NA_real_,
    P_time = NA_real_,
    Beta_int = NA_real_,
    SE_int = NA_real_,
    P_int = NA_real_
  )
  
  for (k in seq_along(proteins)) {
    
    prot <- proteins[k]
    
    # RHS: predictor * time + other covariates
    fixed_terms <- c(
      paste0("scale(", predictor, ") * scale(time)"),
      paste0("scale(", setdiff(numeric_covariates, "time"), ")"),
      factor_covariates
    )
    
    rhs_fixed <- paste(fixed_terms, collapse = " + ")
    
    # LHS is now the protein outcome (scaled)
    fml_fixed <- as.formula(
      paste0("scale(", prot, ") ~ ", rhs_fixed)
    )
    
    fml_mixed <- as.formula(
      paste0("scale(", prot, ") ~ ", rhs_fixed, " + (1 | idno)")
    )
    
    res <- fit_mixed_or_fixed(fml_mixed, fml_fixed, dat)
    fit <- res$fit
    
    out$Nobs[k]  <- stats::nobs(fit)
    out$Model[k] <- res$model_type
    
    # term names to extract
    fa_term   <- paste0("scale(", predictor, ")")
    time_term <- "scale(time)"
    int_term  <- paste0(fa_term, ":scale(time)")
    
    fa_est   <- extract_term(fit, fa_term)
    time_est <- extract_term(fit, time_term)
    int_est  <- extract_term(fit, int_term)
    
    out[k, c("Beta_fa", "SE_fa", "P_fa")] <- fa_est
    out[k, c("Beta_time", "SE_time", "P_time")] <- time_est
    out[k, c("Beta_int", "SE_int", "P_int")] <- int_est
    
    
  }
  

  out <- out |>
    dplyr::mutate(pred_fdr = p.adjust(P_fa, method="fdr"),
                  time_fdr = p.adjust(P_time, method="fdr"),
                  int_fdr = p.adjust(P_int, method="fdr")) |>
    dplyr::rename(protein = Protein) |>
    dplyr::rename_with(~paste("Nobs", predictor, sep="_"), Nobs) |>
    dplyr::rename_with(~paste("model", predictor, sep="_"), Model) |>
    dplyr::rename_with(~paste("beta", predictor, sep="_"), Beta_fa) |>
    dplyr::rename_with(~paste("SE", predictor, sep="_"), SE_fa) |>
    dplyr::rename_with(~paste("P", predictor, sep="_"), P_fa) |>
    dplyr::rename_with(~paste(paste("P", predictor, sep="_"), "fdr", sep="_"), pred_fdr) |>
    dplyr::rename_with(~paste("beta_time", predictor, sep="_"), Beta_time) |>
    dplyr::rename_with(~paste("SE_time", predictor, sep="_"), SE_time) |>
    dplyr::rename_with(~paste("P_time", predictor, sep="_"), P_time) |>
    dplyr::rename_with(~paste(paste("P_time", predictor, sep="_"), "fdr", sep="_"), time_fdr) |>
    dplyr::rename_with(~paste("beta_int_time", predictor, sep="_"), Beta_int) |>
    dplyr::rename_with(~paste("SE_int_time", predictor, sep="_"), SE_int) |>
    dplyr::rename_with(~paste("P_int_time", predictor, sep="_"), P_int) |>
    dplyr::rename_with(~paste(paste(paste("P_time", predictor, sep="_"), "int", sep="_"), "fdr", sep="_"), int_fdr)
  
  

  
  out
}
