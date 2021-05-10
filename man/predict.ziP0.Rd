\name{predict.ziP0}
\alias{predict.ziP0}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Prediction from fitted GAM model with joint binomial--Poisson (ziP0) family}

\description{ Takes a fitted \code{gam} object produced by \code{gam()} with
\code{family = ziP0} and produces predictions (of the outcome, the Poisson part and the
binomial part) given a new set of values for the model covariates or the original values
used for the model fit. Predictions can be accompanied by standard errors, based on the
posterior distribution of the model coefficients.}

\usage{
predict.ziP0(object, newdata, ,se.fit=FALSE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{

 \item{object}{ a fitted \code{gam} object as produced by \code{gam()} with \code{family = ziP0}.}
 \item{newdata}{ A data frame or list containing the values of the model covariates at which predictions
  are required. If this is not provided then predictions corresponding to the
                original data are returned. If \code{newdata} is provided then
                it should contain all the variables needed for prediction: a
                warning is generated if not. See details for use with \code{link{linear.functional.terms}}. }

\item{se.fit}{ when this is TRUE (not default) standard error estimates are returned for each prediction.}

\item{...}{ other arguments, passed to \code{predict.gam}.}

}

\value{A list with items \code{fit}, \code{fit.p} and \code{fit.lambda}, which are each
vectors of the expected value of the outcome, the probability of a non-fixed zero, and
the expected value of the Poisson part, respectively.

If \code{se.fit} is \code{TRUE}, the list also contains vectors \code{se.fit}, \code{se.p} and
\code{se.lambda}, their associated estimated standard errors.
}

\details{\code{fit} and \code{se.fit} are the same as what would be obtained from
\code{predict.gam(object, type = "response")}; further details can be found
at \code{\link[mgcv]{predict.gam}}.

The other components are calculated by transforming these.
}

\author{ Mark W. Donoghoe \email{markdonoghoe@gmail.com}


}

\seealso{  \code{\link[mgcv]{predict.gam}}, \code{\link{ziP0}}}

\examples{
library(ziP0)
set.seed(44)
n <- 200
dat <- gamSim(1, n = n)
dat$y <- rziP0(n = n, lambda = dat$f0, omega = c(0.6, 0.2))

b <- gam(y ~ s(x0), data = dat, family = ziP0)

b.pred <- predict.ziP0(b, se.fit = TRUE)

pp <- par()$mfrow
par(mfrow = c(2, 2))
plot(b, main = "plot.gam")

x0o <- order(dat$x0)
plot(dat$x0[x0o], b.pred$fit[x0o], type = "l", ylim = c(0, 2.9),
     xlab = "x0", ylab = "Exp outcome", main = "fit")
lines(dat$x0[x0o], b.pred$fit[x0o] + qnorm(0.975) * b.pred$se.fit[x0o], lty = 2)
lines(dat$x0[x0o], b.pred$fit[x0o] - qnorm(0.975) * b.pred$se.fit[x0o], lty = 2)

plot(dat$x0[x0o], b.pred$fit.p[x0o], type = "l", ylim = c(0.65, 0.9),
     xlab = "x0", ylab = "Prob of non-fixed zero", main = "fit.p")
lines(dat$x0[x0o], b.pred$fit.p[x0o] + qnorm(0.975) * b.pred$se.p[x0o], lty = 2)
lines(dat$x0[x0o], b.pred$fit.p[x0o] - qnorm(0.975) * b.pred$se.p[x0o], lty = 2)

plot(dat$x0[x0o], b.pred$fit.lambda[x0o], type = "l", ylim = c(0, 2.9),
     xlab = "x0", ylab = "Exp latent count", main = "fit.lambda")
lines(dat$x0[x0o], b.pred$fit.lambda[x0o] + qnorm(0.975) * b.pred$se.lambda[x0o], lty = 2)
lines(dat$x0[x0o], b.pred$fit.lambda[x0o] - qnorm(0.975) * b.pred$se.lambda[x0o], lty = 2)

# Reset the graph window
par(mfrow = pp)
}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..

