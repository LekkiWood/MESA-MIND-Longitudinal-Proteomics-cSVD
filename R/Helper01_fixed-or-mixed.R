fit_mixed_or_fixed <- function(formula_mixed, formula_fixed, data) {
  
  ctrl <- lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5),
    check.conv.singular = "ignore",
    check.conv.grad = "ignore",
    check.conv.hess = "ignore"
  )
  
  # Prefer lmerTest if installed (adds p-values to summary)
  lmer_fun <- if (requireNamespace("lmerTest", quietly = TRUE)) {
    lmerTest::lmer
  } else {
    lme4::lmer
  }
  
  fit <- try(
    lmer_fun(formula_mixed, data = data, control = ctrl),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error") || lme4::isSingular(fit, tol = 1e-4)) {
    fit <- stats::lm(formula_fixed, data = data)
    model_type <- "lm"
  } else {
    model_type <- if (inherits(fit, "lm")) "lm" else "lmer"
  }
  
  list(fit = fit, model_type = model_type)
}

