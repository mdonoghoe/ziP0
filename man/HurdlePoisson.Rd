\name{HurdlePoisson}
\alias{HurdlePoisson}
\alias{dziP}
\alias{pziP}
\alias{qziP}
\alias{rziP}
\title{A Zero-Inflated (Hurdle) Poisson Distribution}
\description{Density, distribution function, quantile function and random generation
  for a zero-inflated (hurdle) Poisson distribution, with parameters \code{lambda} (for the Poisson part)
  and \code{omega} (which relate the hurdle probability to the Poisson mean).}
\usage{
  dziP(x, lambda, omega, log = FALSE)
  pziP(q, lambda, omega, lower.tail = TRUE, log.p = FALSE)
  qziP(p, lambda, omega, lower.tail = TRUE, log.p = FALSE)
  rziP(n, lambda, omega)
}
\arguments{
  \item{x}{vector of (non-negative integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of random values to return.}
  \item{lambda}{vector of (non-negative) means for the Poisson count.}
  \item{omega}{vector of length 2 controlling the relationship between hurdle probability and Poisson mean. See details.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{Pr[X \le x]}, otherwise, \eqn{Pr[X > x]}.}
}
\value{
  \code{dziP} gives the (log) density,
  \code{pziP} gives the (log) distribution function,
  \code{qziP} gives the quantile function, and
  \code{rziP} generates random deviates.

  Invalid \code{lambda} will result in return value \code{NaN}, with a warning.

  The length of the result is determined by \code{n} for
  \code{rziP}, and is the maximum of the lengths of the
  numerical arguments (except \code{omega}) for the other functions.

  The numerical arguments other than \code{n} and \code{omega} are recycled to the
  length of the result. Only the first elements of the logical
  arguments are used.

  \code{rziP} returns a vector of type \link{integer} unless generated
  values exceed the maximum representable integer when \code{\link{double}}
  values are returned since R version 4.0.0.
}
\details{The probability of a zero count is \eqn{1-\pi}{1-\pi}, and with probability \eqn{\pi}{\pi} follows a
truncated Poisson distribution with parameter \code{lambda} = \eqn{\lambda}{\lambda}.

So its density is given by \eqn{Pr(X = 0) = 1-\pi} and

\deqn{Pr(X = x) = \pi \frac{\lambda^x}{(exp(\lambda) - 1) x!}}{Pr(X = x) = \pi \lambda^x/((exp(\lambda) - 1) x!)}
  for \eqn{x = 1, 2, \ldots}.

The Poisson parameter \eqn{\lambda} is related to the non-fixed-zero probability \eqn{\pi}
via the parameter \code{omega} = (\eqn{\omega_1}, \eqn{\omega_2}) such that
\deqn{\log(-\log(1 - \pi)) = \omega_1 + \omega_2 \log(\lambda)}{log(-log(1 - \pi)) = \omega_1 + \omega_2 log(\lambda)}.

We would typically expect a positive association between \eqn{\lambda} and \eqn{\pi},
that is, \eqn{\omega_2 > 0}. A warning will be produced if this is violated.

Note that \eqn{\lambda = 0} is really a limit case resulting in a point mass
at \eqn{0}.

If an element of \code{x} is not integer, the result of \code{dziP} is zero,
with a warning.

The quantile is right continuous: \code{qziP(p, lambda)} is the smallest
integer \eqn{x} such that \eqn{Pr(X \le x) \ge p}.
}
\author{
  Mark W. Donoghoe \email{markdonoghoe@gmail.com}.
}

\seealso{
  \code{\link{Binomial}}, \code{\link{Poisson}}.

  \code{\link[mgcv]{ziP}} for the \code{extended.family} to be used in \code{\link[mgcv]{gam}}.
}
\examples{
require(graphics)

# See what happens as the relationship between binomial probability and Poisson mean changes
pp <- par()$mfrow
par(mfrow = c(2, 2))
x <- seq(0, 10)
plot(x, dziP(x, 3, c(0.5, 0)), type = "h", ylim = c(0, 0.5),
     ylab = "Pr(X = x)", main = "lambda = 3, omega = (0.5, 0)")
plot(x, dziP(x, 3, c(0.5, 0.3)), type = "h", ylim = c(0, 0.5),
     ylab = "Pr(X = x)", main = "lambda = 3, omega = (0.5, 0.3)")
plot(x, dziP(x, 3, c(-0.35, 0)), type = "h", ylim = c(0, 0.5),
     ylab = "Pr(X = x)", main = "lambda = 3, omega = (-0.35, 0)")
plot(x, dziP(x, 3, c(-0.35, 0.3)), type = "h", ylim = c(0, 0.5),
     ylab = "Pr(X = x)", main = "lambda = 3, omega = (-0.35, 0.3)")
# Reset the graph window
par(mfrow = pp)
}

\keyword{models} \keyword{regression}
