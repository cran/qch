#' @import utils
## quiets concerns of R CMD check re: the .'s that appear in pipelines
#if(getRversion() >= "2.15.1")  utils::globalVariables(c(".",">"))
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".",":=",".x",">"))

###################################################################
#' FastKerFdr unsigned
#'
#' @param X a vector of probit-transformed p-values (corresponding to a p-value serie)
#' @param p0 a priori proportion of H0 hypotheses
#' @param plotting boolean, should some diagnostic graphs be plotted. Default is FALSE.
#' @param NbKnot The (maximum) number of knot for the kde procedure. Default is 1e5
#' @param tol a tolerance value for convergence. Default is 1e-5
#' @param max_iter the maximum number of iterations allowed for the algorithm to converge or complete its process.(Default is 1e4.)
#' @import ks graphics stats
#'
#' @return A list of 3 objects.
#' Object 'p0' is an estimate of the proportion of H0 hypotheses,
#' Object 'tau' is the vector of H1 posteriors,
#' Object 'f1' is a numeric vector, each coordinate i corresponding to the evaluation of the H1 density at point xi, where xi is the ith item in X.
#' Object 'F1' is a numeric vector, each coordinate i corresponding to the evaluation of the H1 ;cdf at point xi, where xi is the ith item in X.

FastKerFdr_unsigned <- function(X,p0=NULL,plotting=FALSE,NbKnot=1e5,tol = 1e-5, max_iter=1e4){
  
  n <- length(X)
  
  ## Get a p0 estimate
  if(is.null(p0)){
    p0 <- min(2*sum(X < 0)/n,1-1/n)
  }
  p1 <- 1 - p0
  
  ## Knots, counts and initialization
  if(length(X)>NbKnot){
    Hist <- hist(X, breaks=NbKnot, plot=FALSE)
    Knots <- Hist$mids
    ActualNbKnot <- length(Knots)
    Counts <- Hist$counts
  } else {
    Knots <- X
    ActualNbKnot <- length(X)
    Counts <- rep(1,ActualNbKnot)
  }
  if (plotting){
    Order <- order(Knots)
    Knots <- Knots[Order]
    Counts <- Counts[Order]
  }
  
  h <- ks::hpi(X)
  
  ## Initialize the taus using GM
  phi <- dnorm(Knots)
  MaxKnot <- max(abs(range(Knots))) 
  tau <- 0.8*(purrr::map_dbl(Knots,~max(0,.x)))/MaxKnot+0.1 
  
  ## Get the weighted kernel density estimate
  diff <- 2*tol
  iter <- 0
  while(diff > tol & iter<= 1000){
    iter <- iter + 1
    weights <- tau*Counts
    weights <- ActualNbKnot * weights / sum(weights)
    
    f1 <- ks::kde(x=Knots, w=weights, eval.points=Knots, h = h)$estimate
    
    f1[f1<0] <- 0
    tauNew <- p1*f1/(p0*phi + p1*f1)
    
    ## Dirty job 1: get rid of the f1 mass on the left
    tauNew[Knots< -3] <- 0
    diff <- max(abs(tau - tauNew))
    tau <- tauNew
  }
  if (iter > max_iter & diff > tol) {
    message(paste0("Warning: The algorithm did not converge within max_iter=", max_iter, "."))
  }
  if(plotting){
    Hist.fig <- hist(X, freq=TRUE, breaks=sqrt(n), main='', border=8,
                     xlab="Q-transformed pvalues", ylab="Densities")
    bin.width <- mean(diff(Hist.fig$breaks))
    lines(Knots, n*bin.width*p0*phi, type='l', col=4, lwd=2);
    lines(Knots, n*bin.width*p1*f1, col=2,lwd=2);
    lines(Knots, n*bin.width*(p0*phi+p1*f1), lwd=2)
    legend("topright", legend=c("H0 dist", "H1 dist","Mixture Dist"),
           col=c("blue","red", "black"), lty=c(1,1,2), cex=0.8)
  }
  
  ## Now get the f1 estimate
  KDE <- ks::kde(x=Knots, w=weights)
  f1 <- ks::dkde(x = X,fhat = KDE)
  F1 <- pkde_adapted(fhat = KDE,q = X)
  
  ##test plot 
  plot(X, F1)
  
  ## Dirty job 2: get rid of numeric problems
  f1[f1<0] <- 1e-30
  F1[F1 < 0] <- 1e-30
  
  ## Get better estimate of p0
  f0 <- dnorm(X)
  tau0 <- p0*f0 / (p1*f1 + p0*f0)
  diff <- abs(p0-sum(tau0)/n); iter <- 0;
  while(diff !=0 & iter <1000){
    iter <- iter + 1
    p0 <- sum(tau0)/n
    tau0 <- p0*f0 / ((1-p0)*f1 + p0*f0)
    diff <- abs(p0-sum(tau0)/n)
  }
  
  p1 <- 1 - p0
  tau1 <- p1*f1 / (p1*f1 + p0*f0)
  
  return(list(p0=p0,tau=tau1,f1=f1,F1=F1))
}


###################################################################
#' FastKerFdr signed
#'
#' @param X a vector of probit-transformed p-values (corresponding to a p-value serie)
#' @param p0 a priori proportion of H0 hypotheses
#' @param plotting boolean, should some diagnostic graphs be plotted. Default is FALSE.
#' @param NbKnot The (maximum) number of knot for the kde procedure. Default is 1e5
#' @param tol a tolerance value for convergence. Default is 1e-5
#' @param max_iter the maximum number of iterations allowed for the algorithm to converge or complete its process.(Default is 1e4.)
#' @import ks graphics stats
#' 
#' @return A list of 3 objects.
#' Object 'p0' is an estimate of the proportion of H0 hypotheses,
#' Object 'tau' is the vector of H1 posteriors,
#' Object 'f1' is a numeric vector, each coordinate i corresponding to the evaluation of the H1 density at point xi, where xi is the ith item in X.
#' Object 'F1' is a numeric vector, each coordinate i corresponding to the evaluation of the H1 ;cdf at point xi, where xi is the ith item in X.

FastKerFdr_signed <- function(X,p0=NULL,plotting=FALSE,NbKnot=1e5,tol = 1e-5, max_iter = 1e4){
  
  X <- as.numeric(as.matrix(X))
  n <-  length(X)
  
  ## Get a p0 estimate
  if(is.null(p0)){
    p0 = min(2*sum(X < qnorm(0.5))/n,1-1/n);
  }
  p1 = 1 - p0
  
  ## Knots and counts 
  if(length(X)>NbKnot){
    Hist = hist(X, breaks=NbKnot, plot=FALSE) 
    Knots = Hist$mids
    ActualNbKnot = length(Knots)
    Counts = Hist$counts
  } else {
    Knots = X
    ActualNbKnot = length(X)
    Counts = rep(1,ActualNbKnot) 
  }
  if (plotting){
    Order <- order(Knots)
    Knots <- Knots[Order]
    Counts <- Counts[Order]
  }
  
  h <- ks::hpi(X)
  
  ## Initialize the taus 
  phi <- dnorm(Knots) 
  MaxKnot <- max(abs(range(Knots)))
  tau <- sign(Knots)*0.8*Knots/MaxKnot+0.1 
  
  ## Get the weighted kernel density estimate
  diff <- 2*tol; iter <- 0
  while(diff > tol & iter <= max_iter){
    iter <- iter + 1
    weights <- tau*Counts
    weights <- ActualNbKnot * weights / sum(weights) 
    f1 <- ks::kde(x=Knots, w=weights, eval.points=Knots, h=h)$estimate
    f1[f1<0] <- 0
    tauNew <- p1*f1/(p0*phi + p1*f1)
    diff <- max(abs(tau - tauNew))
    tau <- tauNew
  }
  if (iter > max_iter & diff > tol) {
    message(paste0("Warning: The algorithm did not converge within max_iter=", max_iter, "."))
  }
  
  if(plotting){
    Hist.fig <- hist(X, freq=TRUE, breaks=sqrt(n), main='', border=8, 
                     xlab="Q-transformed pvalues", ylab="Densities")
    bin.width <- mean(diff(Hist.fig$breaks))
    lines(Knots, n*bin.width*p0*phi, type='l', col=4, lwd=2); 
    lines(Knots, n*bin.width*p1*f1, col=2,lwd=2); 
    lines(Knots, n*bin.width*(p0*phi+p1*f1), lwd=2)
    legend("topright", legend=c("H0 dist", "H1 dist","Mixture Dist"),
           col=c("blue","red", "black"), lty=c(1,1,2), cex=0.8)  
  }
  
  ## Now get the f1 estimate
  KDE <- ks::kde(x=Knots, w=weights)
  f1 <- ks::dkde(x = X,fhat = KDE)
  F1 <- pkde_adapted(fhat = KDE,q = X)
  
  
  ## Dirty job 2: get rid of numeric problems
  f1[f1<0] <- 1e-30
  F1[F1<0] <- 1e-30
  
  ## Get better estimate of p0
  f0 <- dnorm(X)
  tau0 <- p0*f0 / (p1*f1 + p0*f0)
  diff <- abs(p0-sum(tau0)/n); iter <- 0;
  while(diff !=0 & iter <1000){
    iter <- iter + 1
    p0 <- sum(tau0)/n
    tau0 <- p0*f0 / ((1-p0)*f1 + p0*f0)
    diff <- abs(p0-sum(tau0)/n)
  }
  
  p1 <- 1 - p0
  tau1 <- p1*f1 / (p1*f1 + p0*f0)
  
  return(list(p0=p0,tau=tau1,f1=f1,F1=F1))
}  

###################################################################
#' Signed case function: Separate f1 into f+ and f-
#'
#' @param XMat a matrix of probit-transformed p-values, each column corresponding to a p-value serie.
#' @param f0Mat a matrix containing the evaluation of the marginal density functions under H0 at each items, each column corresponding to a p-value serie.
#' @param f1Mat a matrix containing the evaluation of the marginal density functions under H1 at each items, each column corresponding to a p-value serie.
#' @param p0 the proportions of H0 items for each series.  
#' @param plotting boolean, should some diagnostic graphs be plotted. Default is FALSE. 
#' 
#' @importFrom dplyr %>%
#'
#' @return A list of 4 objects 'f1plusMat', 'f1minusMat', 'p1plus', 'p1minus'.
#' Object 'f1plusMat' is a matrix containing the evaluation of the marginal density functions under H1plus at each items, each column corresponding to a p-value serie.
#' Object 'f1minusMat' is a matrix containing the evaluation of the marginal density functions under H1minus at each items, each column corresponding to a p-value serie.
#' Object 'p1plus' is an estimate of the proportions of H1plus items for each series.
#' Object 'p1minus' is an estimate of the proportions of H1minus items for each series.

f1_separation_signed <- function(XMat, f0Mat, f1Mat, p0, plotting=FALSE){
  
  n <- nrow(XMat)
  Q <- ncol(XMat)
  
  ### Compute of the integral for normalization
  XMat_temp <- apply(XMat,2, sort)
  XMat_tempm1 <- matrix(1,n,Q)
  f1Mat_temp <- matrix(1,n,Q)
  for(q in 1:Q){
    f1Mat_temp[,q] <- f1Mat[order(XMat[,q]),q]
    XMat_tempm1[,q] <- c(0,XMat_temp[2:n-1,q])
  }
  Xplustemp <- XMat_temp > 0
  lambda_plus <- colSums((XMat_temp[-1,] - XMat_tempm1[-1,]) * (f1Mat_temp[-1,] * Xplustemp[-1,]))
  lambda_minus <- 1 - lambda_plus
  
  ### Compute the distributions f+ and f-
  Xplus <- XMat > 0
  fplusMat <- (f1Mat * Xplus) %>% sweep( 2 , lambda_plus, `/`)
  fminusMat <- (f1Mat * abs(Xplus - 1)) %>% sweep( 2 , lambda_minus, `/`)
  
  ### Compute the proportions p+ and p-
  pplus <- (1 - p0)*lambda_plus
  pminus <- 1 - p0 - pplus
  
  if(plotting){
    for(q in 1:Q){
      Knots <- XMat[,q]
      ActualNbKnot <- length(XMat[,q])
      Counts <- rep(1,ActualNbKnot)
      
      Order <- order(Knots)
      Knots <- Knots[Order]
      Counts <- Counts[Order]
      
      f0 <- f0Mat[Order,q] 
      fplus <- fplusMat[Order,q]
      fminus <- fminusMat[Order,q]
      f1 <- f1Mat[Order,q]
      
      
      Hist.fig <- hist(XMat[,q], freq=TRUE, breaks=sqrt(n), main='Marginal distributions inferred', border=8, 
                       xlab="Q-transformed pvalues", ylab="Densities")
      bin.width <- mean(diff(Hist.fig$breaks))
      lines(Knots, n*bin.width*p0[q]*f0, type='l', col=4, lwd=2); 
      lines(Knots, n*bin.width*pplus[q]*fplus, col=2,lwd=2);
      lines(Knots, n*bin.width*pminus[q]*fminus, col=3,lwd=2); 
      lines(Knots, n*bin.width*(p0[q]*f0+pplus[q]*fplus+pminus[q]*fminus), lwd=2,lty=2)
      legend("topright", legend=c("f0 dist", "f+ dist","f- dist","Mixture dist"),
             col=c(4,2,3,1), lty=c(1,1,1,2), cex=0.8)  
      
    }
  }
  return(list(f1plusMat=fplusMat,f1minusMat=fminusMat,p1plus=pplus,p1minus=pminus))
  
}


#############################################################################
## Cumulative integral for KDE
#############################################################################
integral.kde_adapted <- function(q, fhat, density, tol = 1e-10)
{
  gridsize <- length(fhat$eval.points)
  
  ## Use Simpson's rule to compute numerical integration
  simp.rule <- rep(0, gridsize-1)
  for (i in 1:(gridsize-1))
  {
    del <- fhat$eval.points[i+1] - fhat$eval.points[i]
    simp.rule[i] <- min(fhat$estimate[i], fhat$estimate[i+1])*del + 1/2*abs(fhat$estimate[i+1] - fhat$estimate[i])*del 
  }
  
  ## add last incomplete trapezoid
  q.ind <- findInterval(x=q, vec=fhat$eval.points)
  q.prob <- rep(0, length(q))
  i <- 0
  
  for (qi in q.ind)
  {
    i <- i+1
    
    if (qi==0)
      q.prob[i] <- 0
    else if (qi < gridsize)
    {
      ## linearly interpolate kde 
      fhat.estqi <- (fhat$est[qi+1] - fhat$est[qi])/(fhat$eval[qi+1] - fhat$eval[qi]) * (q[i] - fhat$eval[qi]) + fhat$est[qi]
      delqi <- q[i] - fhat$eval[qi] 
      
      simp.ruleqi <- min(fhat.estqi, fhat$est[qi])*delqi + 1/2*abs(fhat.estqi - fhat$est[qi])*delqi
      q.prob[i] <- sum(simp.rule[1:qi]) + simp.ruleqi
    }
    else
    {
      if (density) q.prob[i] <- 1
      else q.prob[i] <- sum(simp.rule) 
    }
  }
  
  if (density) q.prob[q.prob>=1] <- 1
  
  ## remove possible decreasing values in q.prob
  ##---------- current version ---------------
  # dec.ind <- which(diff(q.prob)<0)
  # while (length(dec.ind)>0) 
  # {
  #     dec.ind <- dec.ind+1 
  #     for (i in dec.ind) q.prob[i] <- q.prob[i-1] 
  #     dec.ind <- which(diff(q.prob)<0)
  # }
  ##------------------------------------------
  
  ##------ proposition of new version --------
  order.q <- order(q)
  dec.ind <- which(diff(q.prob[order.q])<(-tol))
  while (length(dec.ind)>0) 
  {
    dec.ind <- dec.ind+1 
    for (i in dec.ind) q.prob[order.q][i] <- q.prob[order.q][i-1] 
    dec.ind <- which(diff(q.prob[order.q])<(-tol))
  }
  ##------------------------------------------
  
  return(q.prob)
}

## cumulative probability P(fhat <= q)
pkde_adapted <- function(q, fhat)
{
  return(integral.kde_adapted(q=q, fhat=fhat, density=TRUE, tol= 1e-10))
}
