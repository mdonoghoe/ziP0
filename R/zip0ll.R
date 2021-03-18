zip0ll <- function(y,g,eta,deriv=0) {
  ## function to evaluate joint binomial-Poisson log likelihood
  ## and its derivatives w.r.t. g/gamma and eta where
  ## 1-p = exp(-exp(eta)) and lambda = exp(gamma), for each datum in vector y.
  ## p is probability of non-fixed-zero. lambda is Poisson mean
  ## given non-fixed-zero.
  ## deriv: 0 - eval
  ##        1 - grad (l,p) and Hess (ll,lp,pp)
  ##        2 - third derivs lll,llp,lpp,ppp
  ##        4 - 4th derivs. llll,lllp,llpp,lppp,pppp
  l1 <- El2 <- l2 <- l3 <- l4 <- NULL
  zind <- y == -1 ## index of fixed zeros
  yp <- y[!zind]
  l <- et <- exp(eta)
  eg <- exp(g)
  l[zind] <- -et[zind]  # -exp(eta[ind])
  l[!zind] <- mgcv:::l1ee(eta[!zind]) + yp*g[!zind] - eg[!zind] - lgamma(yp+1)
  p <- 1-exp(-et) ## probability of non-fixed-zero

  if (deriv > 0) { ## get first and second derivs
    n <- length(y)
    l1 <- matrix(0, n, 2)
    le <- mgcv:::lde(eta,deriv) ## derivs of ll wrt eta
    l1[!zind,1] <- yp - eg[!zind]      ## l_gamma, y >= 0
    l1[zind,2] <- l[zind]              ## l_eta, y == -1
    l1[!zind,2] <- le$l1[!zind]        ## l_eta, y >= 0

    El2 <- l2 <- matrix(0, n, 3)
    ## order gg, ge, ee
    l2[!zind,1] <- -eg[!zind]          ## l_gg, y >= 0
    l2[zind,3] <- l[zind]              ## l_ee, y == -1
    l2[!zind,3] <- le$l2[!zind]        ## l_ee, y >= 0
    El2[,1] <- -p*eg                   ## E(l_gg)
    El2[,3] <- -(1-p)*et + p*le$l2     ## E(l_ee)
  }
  if (deriv > 1) {
    ## the third derivatives
    ## order ggg,gge,gee,eee
    l3 <- matrix(0, n, 4)
    l3[!zind,1] <- -eg[!zind]          ## l_ggg, y >= 0
    l3[zind,4] <- l[zind]              ## l_eee, y == -1
    l3[!zind,4] <- le$l3[!zind]        ## l_eee, y >= 0
  }
  if (deriv > 3) {
    ## the fourth derivatives
    ## order gggg,ggge,ggee,geee,eeee
    l4 <- matrix(0, n, 5)
    l4[!zind,1] <- -eg[!zind]          ## l_gggg, y >= 0
    l4[zind,5] <- l[zind]              ## l_eeee, y == -1
    l4[!zind,5] <- le$l4[!zind]        ## l_eeee, y >= 0
  }
  list(l=l,l1=l1,l2=l2,l3=l3,l4=l4,El2=El2)
}
