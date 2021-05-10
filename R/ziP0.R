## extended families for mgcv, standard components.
## family - name of family character string
## link - name of link character string
## linkfun - the link function
## linkinv - the inverse link function
## mu.eta - d mu/d eta function (derivative of inverse link wrt eta)
## note: for standard links this information is supplemented using
##       function fix.family.link.extended.family with functions
##       gkg where k is 2,3 or 4 giving the kth derivative of the
##       link over the first derivative of the link to the power k.
##       for non standard links these functions muct be supplied.
## dev.resids - function computing deviance residuals.
## Dd - function returning derivatives of deviance residuals w.r.t. mu and theta.
## aic - function computing twice - log likelihood for 2df to be added to.
## initialize - expression to be evaluated in gam.fit4 and initial.spg
##              to initialize mu or eta.
## preinitialize - optional expression evaluated in estimate.gam to
##                 e.g. initialize theta parameters (see e.g. ocat)
## postproc - expression to evaluate in estimate.gam after fitting (see e.g. betar)
## ls - function to evaluated log saturated likelihood and derivatives w.r.t.
##      phi and theta for use in RE/ML optimization. If deviance used is just -2 log
##      lik. can njust return zeroes.
## validmu, valideta - functions used to test whether mu/eta are valid.
## n.theta - number of theta parameters.
## no.r.sq - optional TRUE/FALSE indicating whether r^2 can be computed for family
## ini.theta - function for initializing theta.
## putTheta, getTheta - functions for storing and retriving theta values in function
##                      environment.
## rd - optional function for simulating response data from fitted model.
## residuals - optional function for computing residuals.
## predict - optional function for predicting from model, called by predict.gam.
## family$data - optional list storing any family specific data for use, e.g. in predict
##               function.

ziP0 <- function (theta = NULL, link = "identity", b = 0) {
  ## joint binomial-Poisson parameterized in terms of the log Poisson parameter, gamma.
  ## eta = theta[1] + exp(theta[2])*gamma, and 1-p = exp(-exp(eta)) where p is
  ## probability of non-fixed-zero

  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("identity")) {
    stats <- make.link(linktemp)
  } else  stop(linktemp, " link not available for zero inflated; available link for `lambda' is only  \"loga\"")
  ## Theta <-  NULL;
  n.theta <- 2
  if (!is.null(theta)) {
    ## fixed theta supplied
    iniTheta <-  c(theta[1],theta[2])
    n.theta <- 0 ## no thetas to estimate
  } else iniTheta <- c(0,0) ## inital theta value - start at 1 - p = Pr(y > 0)

  env <- new.env(parent = environment(ziP0))# new.env(parent = .GlobalEnv)

  if (b<0) b <- 0; assign(".b", b, envir = env)
  assign(".Theta", iniTheta, envir = env)
  getTheta <- function(trans=FALSE) {
    ## trans transforms to the original scale...
    th <- get(".Theta")
    if (trans) {
      th[2] <- get(".b") + exp(th[2])
    }
    th
  }

  putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))

  validmu <- function(mu) all(is.finite(mu))

  dev.resids <- function(y, mu, wt,theta=NULL) {
    ## this version ignores saturated likelihood
    if (is.null(theta)) theta <- get(".Theta")
    b <- get(".b")
    p <- theta[1] + (b + exp(theta[2])) * mu ## l.p. for prob present
    -2*zip0ll(y,mu,p,deriv=0)$l
  }

  Dd <- function(y, mu, theta, wt=NULL, level=0) {
    ## here mu is lin pred for Poisson mean so E(y) = exp(mu)
    ## Deviance for log lik of joint binomial-Poisson.
    ## code here is far more general than is needed - could deal
    ## with any 2 parameter mapping of lp of mean to lp of prob presence.
    if (is.null(theta)) theta <- get(".Theta")
    deriv <- 1; if (level==1) deriv <- 2 else if (level>1) deriv <- 4
    b <- get(".b")
    g <- mgcv:::lind(mu,theta,level,b) ## the derviatives of the transform mapping mu to p
    z <- zip0ll(y,mu,g$p,deriv)
    oo <- list();n <- length(y)
    if (is.null(wt)) wt <- rep(1,n)
    oo$Dmu <- -2*wt*(z$l1[,1] + z$l1[,2]*g$p.l)
    oo$Dmu2 <- -2*wt*(z$l2[,1] + 2*z$l2[,2]*g$p.l + z$l2[,3]*g$p.l^2 + z$l1[,2]*g$p.ll)
    ## WARNING: following requires z$El1 term to be added if transform modified so
    ##          that g$p.ll != 0....
    oo$EDmu2 <- -2*wt*(z$El2[,1] + 2*z$El2[,2]*g$p.l + z$El2[,3]*g$p.l^2)

    if (level>0) { ## l,p - ll,lp,pp -  lll,llp,lpp,ppp - llll,lllp,llpp,lppp,pppp
      oo$Dth <- -2*wt*z$l1[,2]*g$p.th ## l_p p_th
      oo$Dmuth <- -2*wt*(z$l2[,2]*g$p.th + z$l2[,3]*g$p.l*g$p.th + z$l1[,2]*g$p.lth)
      oo$Dmu2th <- -2*wt*(z$l3[,2]*g$p.th + 2*z$l3[,3]*g$p.l*g$p.th + 2* z$l2[,2]*g$p.lth +
                            z$l3[,4]*g$p.l^2*g$p.th + z$l2[,3]*(2*g$p.l*g$p.lth + g$p.th*g$p.ll) + z$l1[,2]*g$p.llth)
      oo$Dmu3 <- -2*wt*(z$l3[,1] + 3*z$l3[,2]*g$p.l + 3*z$l3[,3]*g$p.l^2 + 3*z$l2[,2]*g$p.ll +
                          z$l3[,4]*g$p.l^3 +3*z$l2[,3]*g$p.l*g$p.ll + z$l1[,2]*g$p.lll)
    }
    if (level>1) {
      p.thth <- matrix(0,n,3);p.thth[,1] <- g$p.th[,1]^2
      p.thth[,2] <- g$p.th[,1]*g$p.th[,2];p.thth[,3] <- g$p.th[,2]^2
      oo$Dth2 <- -2*wt*(z$l2[,3]*p.thth + z$l1[,2]*g$p.th2)
      p.lthth <- matrix(0,n,3);p.lthth[,1] <- g$p.th[,1]*g$p.lth[,1]*2
      p.lthth[,2] <- g$p.th[,1]*g$p.lth[,2] + g$p.th[,2]*g$p.lth[,1];
      p.lthth[,3] <- g$p.th[,2]*g$p.lth[,2]*2
      oo$Dmuth2 <- -2*wt*( z$l3[,3]*p.thth + z$l2[,2]*g$p.th2 + z$l3[,4]*g$p.l*p.thth +
                             z$l2[,3]*(g$p.th2*g$p.l + p.lthth) + z$l1[,2]*g$p.lth2)
      p.lthlth <- matrix(0,n,3);p.lthlth[,1] <- g$p.lth[,1]*g$p.lth[,1]*2
      p.lthlth[,2] <- g$p.lth[,1]*g$p.lth[,2] + g$p.lth[,2]*g$p.lth[,1];
      p.lthlth[,3] <- g$p.lth[,2]*g$p.lth[,2]*2
      p.llthth <- matrix(0,n,3);p.llthth[,1] <- g$p.th[,1]*g$p.llth[,1]*2
      p.llthth[,2] <- g$p.th[,1]*g$p.llth[,2] + g$p.th[,2]*g$p.llth[,1];
      p.llthth[,3] <- g$p.th[,2]*g$p.llth[,2]*2

      oo$Dmu2th2 <- -2*wt*(z$l4[,3]*p.thth + z$l3[,2]*g$p.th2 + 2*z$l4[,4] * p.thth *g$p.l + 2*z$l3[,3]*(g$p.th2*g$p.l + p.lthth) +
                             2*z$l2[,2]*g$p.lth2 + z$l4[,5]*p.thth*g$p.l^2 + z$l3[,4]*(g$p.th2*g$p.l^2 + 2*p.lthth*g$p.l + p.thth*g$p.ll) +
                             z$l2[,3]*(p.lthlth + 2*g$p.l*g$p.lth2 + p.llthth + g$p.th2*g$p.ll) + z$l1[,2]*g$p.llth2)

      oo$Dmu3th <- -2*wt*(z$l4[,2]*g$p.th + 3*z$l4[,3]*g$p.th*g$p.l + 3*z$l3[,2]*g$p.lth + 2*z$l4[,4]*g$p.th*g$p.l^2 +
                            z$l3[,3]*(6*g$p.lth*g$p.l + 3*g$p.th*g$p.ll) + 3*z$l2[,2]*g$p.llth + z$l4[,4]*g$p.th*g$p.l^2 +
                            z$l4[,5]*g$p.th*g$p.l^3 + 3*z$l3[,4]*(g$p.l^2*g$p.lth + g$p.th*g$p.l*g$p.ll) +
                            z$l2[,3]*(3*g$p.lth*g$p.ll + 3*g$p.l*g$p.llth + g$p.th*g$p.lll) + z$l1[,2]*g$p.lllth)
      oo$Dmu4 <- -2*wt*(z$l4[,1] + 4*z$l4[,2]*g$p.l + 6*z$l4[,3]*g$p.l^2 + 6*z$l3[,2]*g$p.ll +
                          4*z$l4[,4]*g$p.l^3 + 12*z$l3[,3]*g$p.l*g$p.ll + 4*z$l2[,2]*g$p.lll + z$l4[,5] * g$p.l^4 +
                          6*z$l3[,4]*g$p.l^2*g$p.ll + z$l2[,3] *(4*g$p.l*g$p.lll + 3*g$p.ll^2) + z$l1[,2]*g$p.llll)

    }
    oo
  } ## end Dd for ziP0 (all the same as ziP)

  aic <- function(y, mu, theta=NULL, wt, dev) {
    if (is.null(theta)) theta <- get(".Theta")
    b <- get(".b")
    p <- theta[1] + (b + exp(theta[2])) * mu ## l.p. for prob complete cycle
    sum(-2*wt*zip0ll(y,mu,p,0)$l)
  }

  ls <- function(y,w,theta,scale) {
    ## the log saturated likelihood function.
    ## ls is defined as zero for REML/ML expression as deviance is defined as -2*log.lik
    #vec <- !is.null(attr(theta, "vec.grad"))
    #lsth1 <- if (vec) matrix(0,length(y),2) else c(0,0)
    list(ls=0,## saturated log likelihood
         lsth1=c(0,0),  ## first deriv vector w.r.t theta - last element relates to scale
         lsth2=matrix(0,2,2)) ##Hessian w.r.t. theta
  }

  initialize <- expression({
    if (any(y < -1)) stop("values for the joint binomial-Poisson family must be >= -1")
    if (all.equal(y,round(y))!=TRUE) {
      stop("Non-integer response variables are not allowed with ziP ")
    }
    if (min(y) == 0) stop("No unobserved counts: won't be able to estimate a model")
    ##n <- rep(1, nobs)
    mustart <- log(pmax(0, y) + (y<=0)/5)
  })

  postproc <- function(family,y,prior.weights,fitted,linear.predictors,offset,intercept) {
    posr <- list()
    posr$family <-
      paste("Joint binomial-Poisson(",paste(round(family$getTheta(TRUE),3),collapse=","),")",sep="")
    ## need to fix deviance here!!
    ## wts <- object$prior.weights
    lf <- family$saturated.ll(y,family,prior.weights)
    ## storing the saturated loglik for each datum...
    ##object$family$data <- list(ls = lf)
    l2 <- family$dev.resids(y,linear.predictors,prior.weights)
    posr$deviance <- sum(l2-lf)
    fnull <- function(gamma,family,y,wt) {
      ## evaluate deviance for single parameter model
      sum(family$dev.resids(y, rep(gamma, length(y)), wt))
    }
    meany <- mean(y[y >= 0])
    posr$null.deviance <-
      optimize(fnull,interval=c(meany/5,meany*3),family=family,y=y,wt = prior.weights)$objective - sum(lf)

    ## object$weights <- pmax(0,object$working.weights) ## Fisher can be too extreme
    ## E(y) = p * E(y) - but really can't mess with fitted.values if e.g. rd is to work.
    posr
  } ## postproc

  rd <- function(mu,wt,scale) {
    ## simulate data given fitted latent variable in mu
    rzip0 <- function(gamma,theta) { ## generate ziP0 deviates according to model and lp gamma
      y <- gamma; n <- length(y)
      lambda <- exp(gamma)
      mlam <- max(c(lambda[is.finite(lambda)],.Machine$double.eps^.2))
      lambda[!is.finite(lambda)] <- mlam
      b <- get(".b")
      eta <- theta[1] + (b+exp(theta[2]))*gamma
      p <- 1 - exp(-exp(eta))
      ## generate observed count indicator
      ind <- p > runif(n)
      y[!ind] <- -1
      np <- sum(ind)
      ## simulate from Poisson, given count was observed...
      y[ind] <- rpois(np, lambda[ind])
      y
    }
    rzip0(mu,get(".Theta"))
  }

  saturated.ll <- function(y,family,wt=rep(1,length(y))) {
    ## function to get saturated ll for ziP0 -
    ## actually computes -2 sat ll.
    pind <- y>=0 ## only these are interesting
    wt <- wt[pind]
    y <- y[pind];
    mu <- log(y+0.5)
    keep.on <- TRUE
    theta <- family$getTheta()
    r <- family$Dd(y,mu,theta,wt)
    l <- family$dev.resids(y,mu,wt,theta)
    lmax <- max(abs(l))
    ucov <- abs(r$Dmu) > lmax*1e-7
    k <- 0
    while (keep.on) {
      step <- -r$Dmu/r$Dmu2
      step[!ucov] <- 0
      mu1 <- mu + step
      l1 <- family$dev.resids(y,mu1,wt,theta)
      ind <- l1>l & ucov
      kk <- 0
      while (sum(ind)>0&&kk<50) {
        step[ind] <- step[ind]/2
        mu1 <- mu + step
        l1 <- family$dev.resids(y,mu1,wt,theta)
        ind <- l1>l & ucov
        kk <- kk + 1
      }
      mu <- mu1;l <- l1
      r <- family$Dd(y,mu,theta,wt)
      ucov <- abs(r$Dmu) > lmax*1e-7
      k <- k + 1
      if (all(!ucov)||k==100) keep.on <- FALSE
    }
    l1 <- rep(0,length(pind));l1[pind] <- l
    l1
  } ## saturated.ll

  residuals <- function(object,type=c("deviance","working","response")) {
    if (type == "working") {
      res <- object$residuals
    } else if (type == "response") {
      # For the purpose of response residuals, unobserved counts are 0 (not -1)
      res <- pmax(object$y, 0) - predict.gam(object,type="response")
    } else if (type == "deviance") {
      y <- object$y
      mu <- object$linear.predictors
      wts <- object$prior.weights
      res <- object$family$dev.resids(y,mu,wts)
      ## next line is correct as function returns -2*saturated.log.lik
      res <- res - object$family$saturated.ll(y,object$family,wts)
      fv <- predict.gam(object,type="response")
      s <- attr(res,"sign")
      if (is.null(s)) s <- sign(y-fv)
      res <- as.numeric(sqrt(pmax(res,0)) * s)
    }
    res
  } ## residuals (ziP0)

  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                      beta=NULL,off=NULL,Vb=NULL) {
    ## optional function to give predicted values - idea is that
    ## predict.gam(...,type="response") will use this, and that
    ## either eta will be provided, or {X, beta, off, Vb}. family$data
    ## contains any family specific extra information.

    theta <- family$getTheta()

    if (is.null(eta)) { ## return probabilities
      discrete <- is.list(X)
      ## linear predictor for poisson parameter...
      gamma <- off + if (discrete) mgcv:::Xbd(X$Xd,beta,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop) else drop(X%*%beta)
      if (se) {
        se <- if (discrete) sqrt(pmax(0,mgcv:::diagXVXd(X$Xd,Vb,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,nthreads=1))) else
          sqrt(pmax(0,rowSums((X%*%Vb)*X))) ## se of lin pred
      } else se <- NULL
      #gamma <- drop(X%*%beta + off) ## linear predictor for poisson parameter
      #se <- if (se) drop(sqrt(pmax(0,rowSums((X%*%Vb)*X)))) else NULL ## se of lin pred
    } else { se <- NULL; gamma <- eta}
    ## now compute linear predictor for probability of presence...
    b <- get(".b")
    eta <- theta[1] + (b+exp(theta[2]))*gamma
    et <- exp(eta)
    p <- 1 - exp(-et)
    fv <- lambda <- exp(gamma)
    ind <- gamma < log(.Machine$double.eps)/2
    mu <- lambda
    #mu[!ind] <- lambda[!ind]
    #mu[ind] <- 1
    fv <- list(p*mu)    ## E(y)
    if (is.null(se)) return(fv) else {
      dp.dg <- p
      ind <- eta < log(.Machine$double.xmax)/2
      dp.dg[!ind] <- 0
      dp.dg <- exp(-et)*et*exp(b+theta[2])
      dmu.dg <- mu
      fv[[2]] <- abs(dp.dg*mu+dmu.dg*p)*se
      names(fv) <- c("fit","se.fit")
      return(fv)
    }
  } ## predict

  environment(saturated.ll) <- environment(dev.resids) <- environment(Dd) <-
    environment(aic) <- environment(getTheta) <- environment(rd) <- environment(predict) <-
    environment(putTheta) <- env

  structure(list(family = "joint binomial-Poisson", link = linktemp, linkfun = stats$linkfun,
                 linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd, rd=rd,residuals=residuals,
                 aic = aic, mu.eta = stats$mu.eta, g2g = stats$g2g,g3g=stats$g3g, g4g=stats$g4g,
                 #preinitialize=preinitialize,
                 initialize = initialize,postproc=postproc,ls=ls,no.r.sq=TRUE,
                 validmu = validmu, valideta = stats$valideta,n.theta=n.theta,predict=predict,
                 ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta,saturated.ll = saturated.ll),
            class = c("extended.family","family"))

} ## ziP0
