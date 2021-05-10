predict.ziP0 <- function(object, newdata, se.fit = FALSE, ...) {

  if (substr(object$family$family, 1, 22) != "Joint binomial-Poisson")
    stop("predict.ziP0 only works for objects with family = ziP0")

  ## linear predictors
  pred_gamma <- predict.gam(object, newdata = newdata, type = "link",
                            se.fit = FALSE, ...)
  ## binary part
  omega <- object$family$getTheta(TRUE)
  eta <- omega[1] + omega[2] * pred_gamma
  et <- exp(eta)
  pred_p <- 1 - exp(-et)

  ## based on predict.gam
  pred_resp <- predict.gam(object, newdata = newdata, se.fit = se.fit,
                           type = "response", ...)
  pred_y <- if (se.fit) pred_resp$fit else pred_resp

  ## Poisson part
  pred_lambda <- pred_y / pred_p

  fv <- list(fit = pred_y, fit.p = pred_p, fit.lambda = pred_lambda)

  if (se.fit) {
    fv[["se.fit"]] <- pred_resp$se.fit

    # Need to extract se(gamma) to calculate se(p) and se(lambda)
    dp.dg <- exp(-et) * et * omega[2]
    mu <- dmu.dg <- exp(pred_gamma)
    se.fact <- abs(dp.dg * mu + dmu.dg * pred_p)
    se.gamma <- pred_resp$se.fit / se.fact

    fv[["se.p"]] <- abs(dp.dg) * se.gamma
    fv[["se.lambda"]] <- mu * se.gamma
  }

  fv
}
