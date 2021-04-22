predict.ziP0 <- function(object, newdata, se.fit = FALSE, ...) {

  ## linear predictors
  pred_gamma <- predict.gam(object, newdata = newdata, type = "link",
                            se.fit = FALSE, ...)
  ## binary part
  theta <- object$family$getTheta(TRUE)
  eta <- theta[1] + theta[2] * pred_gamma
  et <- exp(eta)
  pred_p <- 1 - exp(-et)

  ## based on predict.gam
  pred_resp <- predict.gam(object, newdata = newdata, se.fit = se.fit,
                           type = "response", ...)
  pred_y <- if (se.fit) pred_resp$fit else pred_resp

  ## Poisson part
  pred_lambda <- pred_y / pred_p

  fv <- list(fit = pred_resp, fit.p = pred_p, fit.lambda = pred_lambda)
  if (se.fit) fv[["se.fit"]] <- pred_resp$se.fit

  fv
}
