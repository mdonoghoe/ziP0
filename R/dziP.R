dziP <- function(x, lambda, omega, log = FALSE) {
  omega <- check_omega(omega)

  nn <- max(length(x), length(lambda))
  x <- rep_len(x, nn)
  lambda <- rep_len(lambda, nn)

  gamma <- log(lambda)
  prob <- 1 - exp(-exp(omega[1] + omega[2]*gamma))

  if (log) {
    ifelse(x == 0, dbinom(0, size = 1, prob = prob, log = TRUE),
           dbinom(1, size = 1, prob = prob, log = TRUE) +
             dpois(x, lambda = lambda, log = TRUE) -
             ppois(0, lambda = lambda, lower.tail = FALSE, log.p = TRUE))
  } else {
    ifelse(x == 0, dbinom(0, size = 1, prob = prob, log = FALSE),
           dbinom(1, size = 1, prob = prob, log = FALSE) *
             dpois(x, lambda = lambda, log = FALSE) /
             ppois(0, lambda = lambda, lower.tail = FALSE, log.p = FALSE))
  }
}

pziP <- function(q, lambda, omega, lower.tail = TRUE, log.p = FALSE) {
  omega <- check_omega(omega)

  nn <- max(length(q), length(lambda))
  q <- rep_len(q, nn)
  lambda <- rep_len(lambda, nn)

  gamma <- log(lambda)
  prob <- 1 - exp(-exp(omega[1] + omega[2]*gamma))

  lt <- ifelse(q >= 0, dbinom(0, size = 1, prob = prob, log = FALSE), 0) +
    ifelse(q >= 1, dbinom(1, size = 1, prob = prob, log = FALSE) *
             (1 - ppois(q, lambda = lambda, lower.tail = FALSE, log.p = FALSE) /
                ppois(0, lambda = lambda, lower.tail = FALSE, log.p = FALSE)),
           0)

  pz <- if(!lower.tail) 1 - lt else lt
  if (log.p) log(pz) else pz
}

qziP <- function(p, lambda, omega, lower.tail = TRUE, log.p = FALSE) {
  omega <- check_omega(omega)

  nn <- max(length(p), length(lambda))
  p <- rep_len(p, nn)
  lambda <- rep_len(lambda, nn)

  gamma <- log(lambda)
  prob <- 1 - exp(-exp(omega[1] + omega[2]*gamma))

  if (log.p) p <- exp(p)

  gt0 <- ppois(0, lambda = lambda, lower.tail = FALSE, log.p = FALSE)

  if (lower.tail) {
    ifelse(1 - prob >= p, ifelse(p < 0, NaN, 0),
           qpois(1 - gt0 * (1 - p) / prob, lambda = lambda, lower.tail = TRUE,
                 log.p = FALSE))
  } else {
    ifelse(prob <= p, ifelse(p > 1, NaN, 0),
           qpois(p / prob * gt0, lambda = lambda, lower.tail = FALSE,
                 log.p = FALSE))
  }
}

rziP <- function(n, lambda, omega) {
  omega <- check_omega(omega)

  gamma <- log(lambda)
  prob <- 1 - exp(-exp(omega[1] + omega[2]*gamma))

  # Binary part
  z <- rbinom(n, size = 1, prob = prob)
  # Truncated Poisson part
  x <- qpois(runif(n, dpois(0, lambda), 1), lambda)

  # Combined outcome
  x[z == 0L] <- 0L
  x
}
