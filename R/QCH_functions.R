###################################################################
#' Generate H0/H1 configurations and specify the ones corresponding to the composed H1
#'
#' @param Q number of test series to be combined
#' @param AtLeast How many H1 hypotheses at least for the item to be of interest ?
#' @param Consecutive Should the significant test series be consecutive ? Default=FALSE
#'
#' @export
#' @return A list of two objects 'Hconfig' and 'Hconfig.H1'.
#' Hconfig is the list of all possible combination of H0 and H1 hypotheses among Q hypotheses tested.
#' Hconfig.H1 is the vector of components of Hconfig that correspond to the 'AtLeast' specification.
#'
#' @examples
#' GetHinfo(4,2)
#'
#' @seealso [GetHinfoEqual()]
GetHinfo <- function(Q,AtLeast,Consecutive=FALSE){

  ## Build H configurations
  Hconfig <- as.matrix(expand.grid(lapply(1:Q, function(q) 0:1)))
  Hconfig <- split(Hconfig, seq(2^Q))

  ## Find the ones that match H1
  if (!Consecutive){
    MatchingH1 <- sapply(Hconfig, function(h){sum(h)>=AtLeast})
    Hconfig.H1 <- which(MatchingH1)
    names(Hconfig.H1) <- NULL
  } else {
    Consec <- paste(rep(1,AtLeast),collapse='')
    Hconcat <- sapply(Hconfig, function(hh){paste(hh,collapse='')})
    Hconfig.H1 <- grep(pattern = Consec,x = Hconcat)
  }

  ## Collect results
  return(list(Hconfig=Hconfig,Hconfig.H1=Hconfig.H1))

}


###################################################################
#' Generate H0/H1 configurations and specify the ones corresponding to the composed H1
#'
#' @param Q number of test series to be combined
#' @param Equal How many H1 hypotheses exactly for the item to be of interest ?
#' @param Consecutive Should the significant test series be consecutive ? Default=FALSE
#'
#' @export
#' @return A list of two objects 'Hconfig' and 'Hconfig.H1'.
#' Hconfig is the list of all possible combination of H0 and H1 hypotheses among Q hypotheses tested.
#' Hconfig.H1 is the vector of components of Hconfig that correspond to the 'Equal' specification.
#'
#' @examples
#' GetHinfoEqual(4,2)
#' @seealso [GetHinfo()]
GetHinfoEqual <- function(Q,Equal,Consecutive=FALSE){

  ## Build H configurations
  Hconfig <- as.matrix(expand.grid(lapply(1:Q, function(q) 0:1)))
  Hconfig <- split(Hconfig, seq(2^Q))

  ## Find the ones that match H1
  if (!Consecutive){
    MatchingH1 <- sapply(Hconfig, function(h){sum(h)==Equal})
    Hconfig.H1 <- which(MatchingH1)
    names(Hconfig.H1) <- NULL
  } else {
    Consec <- paste(rep(1,Equal),collapse='')
    ConsecP1 <- paste(rep(1,Equal+1),collapse='')
    Hconcat <- sapply(Hconfig, function(hh){paste(hh,collapse='')})
    Hconfig.H1 <- intersect(grep(pattern = Consec,x = Hconcat),which(sapply(Hconfig, function(h){sum(h)==Equal})))
  }

  ## Collect results
  return(list(Hconfig=Hconfig,Hconfig.H1=Hconfig.H1))

}



###################################################################
#' FastKerFdr
#'
#' @param Pval a vector of p-values (corresponding to a p-value serie)
#' @param p0 a priori proportion of H0 hypotheses
#' @param plotting boolean, should some diagnostic graphs be plotted. Default is FALSE.
#' @param NbKnot The (maximum) number of knot for the kde procedure. Default is 1e5
#' @param tol a tolerance value for convergence. Default is 1e-5
#' @import ks mclust graphics stats
#'
#' @return A list of 3 objects. Object p0 is an estimate of the proportion of H0 hypotheses.,
#' tau is the vector of H1 posteriors.
#' f1 is a numeric vector, each coordinate i corresponding to the evaluation of the H1 density at point pi, where pi is the ith p-value in Pval.
FastKerFdr <- function(Pval,p0=NULL,plotting=FALSE,NbKnot=1e5,tol = 1e-5){

  ## Transform pvalues into N(0,1) quantiles
  n = length(Pval)
  X = -qnorm(Pval)

  ## Get a p0 estimate
  if(is.null(p0)){
    p0 = min(2*sum(X < 0)/n,1-1/n);
  }
  p1 = 1 - p0

  ## Knots, counts and initialization (using Mclust)
  if(length(X)>NbKnot){
    Hist = hist(X, breaks=NbKnot, plot=FALSE)
    Knots = Hist$mids; ActualNbKnot = length(Knots); Counts = Hist$counts;
    Xsample = sample(X, NbKnot)
    GM = mclust::Mclust(Xsample, G=3, modelNames='E');
    mu = max(GM$parameters$mean)
  } else {
    Knots = X; ActualNbKnot = length(X); Counts = rep(1,ActualNbKnot);
    GM = mclust::Mclust(X, G=3, modelNames='E');
    mu = max(GM$parameters$mean)
  }
  if (plotting){
    Order <- order(Knots)
    Knots <- Knots[Order]
    Counts <- Counts[Order]
  }

  ## Initialize the taus using GM
  phi = dnorm(Knots); f1 = dnorm(Knots, mean=mu, sd=1)
  tau = p1*f1/(p0*phi + p1*f1)

  ## Get the weighted kernel density estimate
  diff = 2*tol; iter = 0
  while(diff > tol){
    iter = iter + 1
    weights = tau*Counts; weights = ActualNbKnot * weights / sum(weights)
    f1 = ks::kde(x=Knots, w=weights, eval.points=Knots)$estimate
    tauNew = p1*f1/(p0*phi + p1*f1)
    ## Dirty job 1: get rid of the f1 mass on the left
    tauNew[Knots< -3] <- 0
    diff = max(abs(tau - tauNew))
    tau = tauNew
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
  KDE = ks::kde(x=Knots, w=weights, eval.points=X)
  f1 = KDE$estimate

  ## Dirty job 2: get rid of numeric problems
  f1[f1<0] <- 1e-30
  tau = p1*f1 / (p1*f1 + p0*dnorm(X))

  return(list(p0=p0,tau=tau,f1=f1))
}


###################################################################
#' Infer Hconfig posteriors
#'
#' @param pValMat a matrix of p-values, each column corresponding to a p-value serie.
#' @param Hconfig an Hconfig list as generated by the [GetHinfo()] function.
#' @param plotting a boolean. Should some diagnostic graphs be plotted ? Default is FALSE.
#' @import stats
#'
#' @return A list of 2 objects 'prior' and 'posterior'.
#' Object 'prior' is a vector of estimated prior probabilities for each of the H-configurations.
#' Object 'posterior' is a matrix providing for each item (in row) its posterior probability to belong to each of the H-configurations (in columns).
#' @export
#'
#' @examples
#' data(PvalSets)
#' PvalMat <- as.matrix(PvalSets[,-3])
#' ## Build the Hconfig objects
#' Q <- 2
#' AtLeast <- 2
#' Hconfig <- GetHinfo(Q,AtLeast)$Hconfig
#'
#' ## Run the function
#' res.fit <- qch.fit(PvalMat,Hconfig)
#'
#' ## Display the prior of each class of items
#' res.fit$prior
#'
#' ## Display the first posteriors
#' head(res.fit$posterior)
qch.fit <- function(pValMat,Hconfig, plotting=FALSE){

  n <- nrow(pValMat)
  Q <- ncol(pValMat)

  #### Step 1: Marginal density estimation

  ## Get p0 estimates
  p0 <- rep(0, Q)
  for (q in 1:Q){
    p0[q] = min(2*sum(pValMat[,q] > 0.5)/n,1-1/n);
  }
  SomeH1 <-  which(p0<1)
  NoH1 <- which(p0==1-1/n);
  if(length(NoH1)==1){
    message(paste("Pvalue serie",NoH1, "may have very few H1 (or a weird distribution)"))
  }
  if(length(NoH1)>1){
    message(paste("Pvalue series",paste(NoH1,collapse=' '), "may have very few H1 (or a weird distribution)"))
  }

  ## Fit a 2-component mixture to each test serie using kerFdr
  f1Mat <- matrix(1, n, Q);
  for(q in SomeH1){
    ker <- FastKerFdr(pValMat[, q], p0=p0[q], plotting=FALSE)
    f1Mat[,q] <- ker$f1
  }
  f0Mat <- matrix(dnorm(-qnorm(pValMat)),ncol=Q)

  #### Step 2: transform marginal densities into config densities

  Logf0Mat <- log(f0Mat);
  Logf1Mat <- log(f1Mat);
  f.Hconfig <- sapply(Hconfig, function(h){
    f <- rep(0,nrow(Logf0Mat))
    if (length(which(h==1)) > 0){f <- f + rowSums(Logf1Mat[, which(h==1), drop=FALSE])}
    if (length(which(h==0)) > 0){f <- f + rowSums(Logf0Mat[, which(h==0), drop=FALSE])}
    return(exp(f))
  })

  #### Step 3: Infer prior estimation using an EM procedure

  ## Initialization: simple product of marginal priors estimator
  NewPrior <- sapply(1:length(Hconfig), function(c){
    prod(p0[which(Hconfig[[c]]==0)]) * prod(1-p0[which(Hconfig[[c]]==1)])
  })
  PriorsAt0 <- which(NewPrior==0)

  ## EM calibration
  NotOK <- TRUE
  Precision <- 1e-6
  NoLowerThan <- 1e-7
  while(NotOK){

    ## E step
    Tau <- f.Hconfig*(tcrossprod(rep(1:n),NewPrior))
    Tau <- Tau/rowSums(Tau)

    ## M step
    OldPrior <- NewPrior
    NewPrior <- colMeans(Tau)
    if(length(PriorsAt0)==0){
      NewPrior[NewPrior<NoLowerThan] <- NoLowerThan
    } else {
      NoLowerCoord <- setdiff(which(NewPrior<NoLowerThan),PriorsAt0)
      if(length(NoLowerCoord)>0){
        NewPrior[NoLowerCoord] <- NoLowerThan
      }
    }
    NewPrior <- NewPrior/sum(NewPrior)
    NotOK <- max((OldPrior-NewPrior)^2) > Precision

  }
  priorHconfigEM <- NewPrior

  #### Step 4: Posterior computation
  posterior <- f.Hconfig*(tcrossprod(rep(1:n),priorHconfigEM))
  posterior <- posterior/rowSums(posterior)

  #### Last but not least: output results
  Res <- list(prior=priorHconfigEM, posterior=posterior)
  return(Res)
}


###################################################################
#' Perform composed hypothesis testing with FDR control
#'
#' @param posterior a matrix of posterior probabilities for each item to belong the different H-configurations, as provided by the [qch.fit()] function.
#' @param Hconfig.H1 a list of H1 config, as created by the [GetHinfo()] function.
#' @param Alpha the nominal Type I error rate for FDR control.
#'
#' @return A list of 2 objects 'Rejection' and 'lFDR'.
#' Object 'Rejection' is a vector providing for each item the result of the composed hypothesis test, after multiple testing correction.
#' Object 'lFDR' is a vector providing for each item its local FDR estimate.
#' @export
#'
#' @examples
#' data(PvalSets)
#' PvalMat <- as.matrix(PvalSets[,-3])
#' Truth <- PvalSets[,3]
#'
#' ## Build the Hconfig objects
#' Q <- 2
#' AtLeast <- 2
#' Hconfig <- GetHinfo(Q,AtLeast)$Hconfig
#' Hconfig.H1 <- GetHinfo(Q,AtLeast)$Hconfig.H1
#'
#' ## Infer the posteriors
#' res.fit <- qch.fit(PvalMat,Hconfig)
#'
#' ## Run the test procedure with FDR control
#' res.test <- qch.test(res.fit$posterior,Hconfig.H1)
#' table(res.test$Rejection,Truth==4)
qch.test <- function(posterior,Hconfig.H1,Alpha=0.05){

  n <- nrow(posterior)
  localFDR <- 1-rowSums(posterior[,Hconfig.H1,drop=FALSE])
  Order <- order(localFDR)
  FDR <- cumsum(localFDR[Order])/(1:n)
  NbReject <- max(which(FDR<=Alpha))
  Rejection <- rep(0,n)
  if (NbReject>0){
    Rejection[Order[1:NbReject]] <- 1
  }
  return(list(Rejection=Rejection,lFDR=localFDR))

}

