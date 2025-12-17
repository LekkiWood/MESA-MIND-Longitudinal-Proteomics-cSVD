cross_sectional_PWAS_function <- function(cleaned_proteins, protein_mapping,
                                          traits_db, numeric_covariates, factor_covariates, outcome, chosen_exam) {
  
  
  common_keys <- intersect(names(cleaned_proteins), names(traits_db))
  
  PWAS_file <- traits_db |>
    dplyr::left_join(cleaned_proteins, by = common_keys)
  
  proteins <- protein_mapping |>
    dplyr::filter(OlinkID %in% names(PWAS_file)) |>
    dplyr::pull(OlinkID)
  
  
  
  # confirm required columns exist
  stopifnot(all(c(outcome, factor_covariates, numeric_covariates, "Exam") %in% names(PWAS_file)))
  
  # subset once (faster + avoids repeated subset() inside loop)
  dat <- PWAS_file |> 
    dplyr::filter(Exam == chosen_exam) |> 
    droplevels()
  
  # storage
  nvar <- length(proteins)
  out <- data.frame(
    protein = proteins,
    nobs    = NA_real_,
    beta    = NA_real_,
    se      = NA_real_,
    pvalue  = NA_real_
  )
  
  for (k in seq_along(proteins)) {
    predictor <- proteins[k]
    
    rhs_terms <- c(
      paste0("scale(", predictor, ")"),
      paste0("scale(", numeric_covariates, ")"),
      factor_covariates
    )
    
    fml <- stats::as.formula(
      paste0("scale(", outcome, ") ~ ", paste(rhs_terms, collapse = " + "))
    )
    
    fit <- stats::lm(fml, data = dat)
    coefs <- summary(fit)$coefficients
    
    # row 2 is the (scaled) predictor term; safe here because we put it first
    out$nobs[k]   <- nobs(fit)
    out$beta[k]   <- unname(coefs[2, 1])
    out$se[k]     <- unname(coefs[2, 2])
    out$pvalue[k] <- unname(coefs[2, 4])
  }
  
  out$fdr <- p.adjust(out$pvalue, method = "fdr")
  names(out) <- paste(paste(outcome, chosen_exam, sep="_"), names(out), sep="_")
  
  results <- as.data.frame(out) |>
    dplyr::rename(OlinkID = paste(paste(outcome, chosen_exam, sep="_"), "protein", sep="_")) |>
    dplyr::left_join(protein_mapping, dplyr::join_by(OlinkID))
  
  ######Output
  results
}