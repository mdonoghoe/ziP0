\name{ziP0}
\alias{ziP0}
\alias{zip0ll}
\title{GAM joint binomial-Poisson regression family}
\description{Family for use with \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}, implementing regression for joint binomial-Poisson data
when the complimentary log log of the binomial (fixed) zero probability is linearly dependent on the log of the Poisson parameter. This is distinct from the usual zero-inflated Poisson model (\code{family=\link[mgcv]{ziP}}) because we must be able to identify which zeros are 'fixed' zeroes from the binomial part, and which are 'true' zeroes from the Poisson part. Fixed zeroes should appear as \code{-1} in the data; true zeroes as \code{0}.

\code{zip0ll} evaluates the log-likelihood and its derivatives. It is a helper function that should not be called directly.
}
\usage{
ziP0(theta = NULL, link = "identity", b = 0)
zip0ll(y, g, eta, deriv = 0)
}
\arguments{
\item{theta}{the 2 parameters controlling the slope and intercept of the
linear transform of the mean controlling the zero inflation rate. If supplied
then treated as fixed parameters (\eqn{\theta_1}{theta_1} and \eqn{\theta_2}{theta_2}), otherwise estimated.}
\item{link}{The link function: only the \code{"identity"} is currently supported.}
\item{b}{a non-negative constant, specifying the minimum dependence of the zero inflation rate on the linear predictor.}
\item{y}{Observed values. Fixed zeroes should appear as \code{-1}.}
\item{g}{Value of \eqn{\gamma = \log(\mu)}{gamma = log(mu)} parameter.}
\item{eta}{Value of \eqn{\eta}{eta} parameter.}
\item{deriv}{Which derivatives are desired. 0 = evaluate log-likelihood only; 1 = also grad and Hessian; 2 = also third derivatives; 4 = also fourth derivatives.}
}
\details{The outcome is a fixed zero (\code{-1}) with probability \eqn{1-\pi}{1-\pi},
and with probability \eqn{\pi}{\pi} follows a Poisson distribution with mean \eqn{\mu}{mu}.

The linear predictor gives \eqn{\gamma = \log \mu}{gamma = log(mu)}, while
\eqn{\eta = \log(-\log(1-\pi)) }{eta=log(-log(1-\pi))} and
\eqn{\eta = \theta_1 + \{b+\exp(\theta_2)\} \log \mu }{eta = theta_1 + (b+exp(theta_2)) gamma}.

The \code{theta} parameters are estimated alongside the smoothing parameters.
Increasing the \code{b} parameter from zero can greatly reduce identifiability
problems, particularly when there are very few non-fixed-zero data.

The fitted values for this model are the log of the Poisson parameter. Use the
\code{predict} function with \code{type="response"} to get the predicted expected
response (where fixed zeroes are treated as zeroes), or \code{\link{predict.ziP0}}
to get the binomial and Poisson expectations separately.

Note that the theta parameters reported in model summaries (or with
\code{getTheta(trans = TRUE)}) are
\eqn{\omega_1 = \theta_1}{omega_1 = theta_1} and \eqn{\omega_2 = b + \exp(\theta_2)}{omega_2 = b + exp(theta_2)}.
}
\value{
An object of class \code{extended.family}.
}
\references{
Wood, S.N., N. Pya and B. Saefken (2016), Smoothing parameter and
model selection for general smooth models.
Journal of the American Statistical Association 111, 1548-1575
\doi{10.1080/01621459.2016.1180986}
}
\author{
Mark W. Donoghoe \email{markdonoghoe@gmail.com}; based heavily on \code{\link[mgcv]{ziP}} by Simon N. Wood \email{simon.wood@r-project.org}.
}

\seealso{
\code{\link[mgcv]{ziP}}, \code{\link[mgcv]{ziplss}}
}
\examples{
set.seed(99)
n <- 400
dat <- gamSim(1, n = n)
dat$y <- rziP0(n = n, lambda = exp(dat$f/4 - 1), omega = c(-2, 1.3))

b <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), family = ziP0(), data = dat)

b
plot(b,pages=1)
}

\keyword{models} \keyword{regression}
