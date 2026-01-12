extract_term <- function(fit, term) {
  coefs <- summary(fit)$coefficients
  
  if (!term %in% rownames(coefs)) {
    return(c(beta = NA_real_, se = NA_real_, p = NA_real_))
  }
  
  beta <- unname(coefs[term, "Estimate"])
  se   <- unname(coefs[term, "Std. Error"])
  
  # lm: "Pr(>|t|)"
  # lmerTest: usually "Pr(>|t|)" (or sometimes "Pr(>|z|)" depending on method)
  pcol <- intersect(colnames(coefs), c("Pr(>|t|)", "Pr(>|z|)"))
  
  if (length(pcol) == 1) {
    p <- unname(coefs[term, pcol])
  } else {
    # No p-values available (plain lme4). Return NA (or see option B below).
    p <- NA_real_
  }
  
  c(beta = beta, se = se, p = p)
}
