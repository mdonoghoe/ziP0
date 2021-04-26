# Density, distribution, quantile and random generation for ziP0
check_omega <- function(omega) {
  if (length(omega) != 2L) stop("omega should have length 2")
  if (omega[2] < 0) warning("omega[2] < 0 implies probability of non-zero is inversely related to lambda")
  omega
}

dziP0 <- function(x, lambda, omega, log = FALSE) {
  omega <- check_omega(omega)

  nn <- max(length(x), length(lambda))
  x <- rep_len(x, nn)
  lambda <- rep_len(lambda, nn)

  gamma <- log(lambda)
  prob <- 1 - exp(-exp(omega[1] + omega[2]*gamma))

  if (log) {
    ifelse(x == -1, dbinom(0, size = 1, prob = prob, log = TRUE),
           dbinom(1, size = 1, prob = prob, log = TRUE) +
             dpois(x, lambda = lambda, log = TRUE))
  } else {
    ifelse(x == -1, dbinom(0, size = 1, prob = prob, log = FALSE),
           dbinom(1, size = 1, prob = prob, log = FALSE) *
             dpois(x, lambda = lambda, log = FALSE))
  }
}

pziP0 <- function(q, lambda, omega, lower.tail = TRUE, log.p = FALSE) {
  omega <- check_omega(omega)

  nn <- max(length(q), length(lambda))
  q <- rep_len(q, nn)
  lambda <- rep_len(lambda, nn)

  gamma <- log(lambda)
  prob <- 1 - exp(-exp(omega[1] + omega[2]*gamma))

  lt <- ifelse(q >= -1, dbinom(0, size = 1, prob = prob, log = FALSE), 0) +
    dbinom(1, size = 1, prob = prob, log = FALSE) *
    ppois(q, lambda = lambda, lower.tail = TRUE, log.p = FALSE)

  pz <- if(!lower.tail) 1 - lt else lt
  if (log.p) log(pz) else pz
}

qziP0 <- function(p, lambda, omega, lower.tail = TRUE, log.p = FALSE) {
  omega <- check_omega(omega)

  nn <- max(length(p), length(lambda))
  p <- rep_len(p, nn)
  lambda <- rep_len(lambda, nn)

  gamma <- log(lambda)
  prob <- 1 - exp(-exp(omega[1] + omega[2]*gamma))

  if (log.p) p <- exp(p)

  if (lower.tail) {
    ifelse(1 - prob >= p, ifelse(p < 0, NaN, -1),
           qpois(1 - (1 - p) / prob, lambda = lambda, lower.tail = TRUE, log.p = FALSE))
  } else {
    ifelse(prob <= p, ifelse(p > 1, NaN, -1),
           qpois(p / prob, lambda = lambda, lower.tail = FALSE, log.p = FALSE))
  }
}

rziP0 <- function(n, lambda, omega) {
  omega <- check_omega(omega)

  gamma <- log(lambda)
  prob <- 1 - exp(-exp(omega[1] + omega[2]*gamma))

  # Binary part
  z <- rbinom(n, size = 1, prob = prob)
  # Poisson part
  x <- rpois(n, lambda)

  # Combined outcome
  x[z == 0L] <- -1L
  x
}
