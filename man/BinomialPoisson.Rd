\name{BinomialPoisson}
\alias{BinomialPoisson}
\alias{dziP0}
\alias{pziP0}
\alias{qziP0}
\alias{rziP0}
\title{A Joint Binomial--Poisson Distribution}
\description{Density, distribution function, quantile function and random generation
  for a joint binomial--Poisson distribution, with parameters \code{lambda} (for the Poisson part)
  and \code{omega} (which relate the binomial probability to the Poisson mean).}
\usage{
  dziP0(x, lambda, omega, log = FALSE)
  pziP0(q, lambda, omega, lower.tail = TRUE, log.p = FALSE)
  qziP0(p, lambda, omega, lower.tail = TRUE, log.p = FALSE)
  rziP0(n, lambda, omega)
}
\arguments{
  \item{x}{vector of (integer, \code{>= -1}) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of random values to return.}
  \item{lambda}{vector of (non-negative) means for the Poisson count.}
  \item{omega}{vector of length 2 controlling the relationship between binomial probability and Poisson mean. See details.}
  \item{theta}{the 2 parameters controlling the slope and intercept of the
    linear transform of the mean controlling the zero inflation rate. If supplied
    then treated as fixed parameters (\eqn{\theta_1}{theta_1} and \eqn{\theta_2}{theta_2}), otherwise estimated.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{Pr[X \le x]}, otherwise, \eqn{Pr[X > x]}.}
}
\value{
  \code{dziP0} gives the (log) density,
  \code{pziP0} gives the (log) distribution function,
  \code{qziP0} gives the quantile function, and
  \code{rziP0} generates random deviates.

  Invalid \code{lambda} will result in return value \code{NaN}, with a warning.

  The length of the result is determined by \code{n} for
  \code{rziP0}, and is the maximum of the lengths of the
  numerical arguments (except \code{omega}) for the other functions.

  The numerical arguments other than \code{n} and \code{omega} are recycled to the
  length of the result. Only the first elements of the logical
  arguments are used.

  \code{rziP0} returns a vector of type \link{integer} unless generated
  values exceed the maximum representable integer when \code{\link{double}}
  values are returned since R version 4.0.0.
}
\details{The outcome takes value \code{-1} (representing a fixed zero) with
probability \eqn{1-\pi}{1-\pi}, and with probability \eqn{\pi}{\pi} follows a Poisson
distribution with mean \code{lambda} = \eqn{\lambda}{\lambda}.

So its density is given by \eqn{Pr(X = -1) = 1-\pi} and

\deqn{Pr(X = x) = \pi \frac{\lambda^x e^{-\lambda}}{x!}}{Pr(X = x) = \pi \lambda^x exp(-\lambda)/x!}
  for \eqn{x = 0, 1, 2, \ldots}.

The Poisson mean \eqn{\lambda} is related to the non-fixed-zero probability \eqn{\pi}
via the parameter \code{omega} = (\eqn{\omega_1}, \eqn{\omega_2}) such that
\deqn{\log(-\log(1 - \pi)) = \omega_1 + \omega_2 \log(\lambda)}{log(-log(1 - \pi)) = \omega_1 + \omega_2 log(\lambda)}.

We would typically expect a positive association between \eqn{\lambda} and \eqn{\pi},
that is, \eqn{\omega_2 > 0}. A warning will be produced if this is violated.

Note that \eqn{\lambda = 0} is really a limit case resulting in a point mass
at \eqn{-1}.

If an element of \code{x} is not integer, the result of \code{dziP0} is zero,
with a warning.

The quantile is right continuous: \code{qziP0(p, lambda)} is the smallest
integer \eqn{x} such that \eqn{Pr(X \le x) \ge p}.
}
\author{
  Mark W. Donoghoe \email{markdonoghoe@gmail.com}.
}

\seealso{
  \code{\link{Binomial}}, \code{\link{Poisson}}.

  \code{\link{ziP0}} for the \code{extended.family} to be used in \code{\link[mgcv]{gam}}.
}
\examples{
require(graphics)

# See what happens as the relationship between binomial probability and Poisson mean changes
pp <- par()$mfrow
par(mfrow = c(2, 2))
x <- seq(-1, 10)
plot(x, dziP0(x, 3, c(0.5, 0)), type = "h", ylim = c(0, 0.5),
     ylab = "Pr(X = x)", main = "lambda = 3, omega = (0.5, 0)")
plot(x, dziP0(x, 3, c(0.5, 0.3)), type = "h", ylim = c(0, 0.5),
     ylab = "Pr(X = x)", main = "lambda = 3, omega = (0.5, 0.3)")
plot(x, dziP0(x, 3, c(-0.35, 0)), type = "h", ylim = c(0, 0.5),
     ylab = "Pr(X = x)", main = "lambda = 3, omega = (-0.35, 0)")
plot(x, dziP0(x, 3, c(-0.35, 0.3)), type = "h", ylim = c(0, 0.5),
     ylab = "Pr(X = x)", main = "lambda = 3, omega = (-0.35, 0.3)")
# Reset the graph window
par(mfrow = pp)
}

\keyword{models} \keyword{regression}
