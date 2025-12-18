cross_sectional_logistic_PWAS_function <- function(cleaned_proteins, protein_mapping,
                                                   traits_db, numeric_covariates, factor_covariates, outcome, chosen_exam) {
  
  common_keys <- intersect(names(cleaned_proteins), names(traits_db))
  
  PWAS_file <- traits_db |>
    dplyr::left_join(cleaned_proteins, by = common_keys)
  
  proteins <- protein_mapping |>
    dplyr::filter(OlinkID %in% names(PWAS_file)) |>
    dplyr::pull(OlinkID)
  
  factor_covs <- factor_covariates
  
  numeric_covs <- numeric_covariates
  
  # confirm required columns exist
  stopifnot(all(c(outcome, factor_covs, numeric_covs, "Exam") %in% names(PWAS_file)))
  
  # subset once (faster + avoids repeated subset() inside loop)
  dat <- PWAS_file |> 
    dplyr::filter(Exam == chosen_exam) |> 
    droplevels()
  
  
  # --- ensure outcome is 0/1 numeric for logistic ---
  y <- dat[[outcome]]
  if (is.factor(y)) y <- droplevels(y)
  
  if (is.logical(y)) {
    dat[[outcome]] <- as.integer(y)
  } else if (is.factor(y)) {
    # assumes 2-level factor; uses second level as "1"
    if (nlevels(y) != 2) stop("Outcome factor does not have 2 levels: ", outcome)
    dat[[outcome]] <- as.integer(y == levels(y)[2])
  } else {
    dat[[outcome]] <- as.integer(y)
  }
  
  # basic validation: must be only 0/1 (allow NA)
  vals <- unique(dat[[outcome]][!is.na(dat[[outcome]])])
  if (!all(vals %in% c(0L, 1L))) {
    stop("Outcome must be binary (0/1) after coercion. Found values: ",
         paste(vals, collapse = ", "))
  }
  
  # storage
  nvar <- length(proteins)
  out <- data.frame(
    protein  = proteins,
    nobs     = NA_real_,
    OR       = NA_real_,
    CI_lower = NA_real_,
    CI_upper = NA_real_,
    pvalue   = NA_real_
  )
  
  for (k in seq_along(proteins)) {
    predictor <- proteins[k]
    
    rhs_terms <- c(
      paste0("scale(", predictor, ")"),
      paste0("scale(", numeric_covs, ")"),
      factor_covs
    )
    
    # DO NOT scale outcome for logistic regression
    fml <- stats::as.formula(
      paste0(outcome, " ~ ", paste(rhs_terms, collapse = " + "))
    )
    
    fit <- stats::glm(fml, family = stats::binomial(link = "logit"), data = dat)
    
    # coefficient for the (scaled) protein term is the first RHS term => row 2
    b  <- stats::coef(fit)[2]
    se <- sqrt(diag(sandwich::vcovHC(fit, type = "HC0")))[2]
    
    OR       <- exp(b)
    CI_lower <- exp(b - 1.96 * se)
    CI_upper <- exp(b + 1.96 * se)
    z_value  <- b / se
    p_value  <- 2 * stats::pnorm(abs(z_value), lower.tail = FALSE)
    
    out$nobs[k]     <- stats::nobs(fit)
    out$OR[k]       <- OR
    out$CI_lower[k] <- CI_lower
    out$CI_upper[k] <- CI_upper
    out$pvalue[k]   <- p_value
  }
  
  out$fdr <- p.adjust(out$pvalue, method = "fdr")
  names(out) <- paste(paste(outcome, chosen_exam, sep="_"), names(out), sep="_")
  
  results <- as.data.frame(out) |>
    dplyr::rename(OlinkID = paste(paste(outcome, chosen_exam, sep="_"), "protein", sep="_")) |>
    dplyr::left_join(protein_mapping, dplyr::join_by(OlinkID))
  
  ######Output
  results
}
