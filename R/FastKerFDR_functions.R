#' @import utils
## quiets concerns of R CMD check re: the .'s that appear in pipelines
# if(getRversion() >= "2.15.1")  utils::globalVariables(c(".",">"))
if (getRversion() >= "2.15.1") utils::globalVariables(c(".", ":=", ".x", ">"))

###################################################################
#' FastKerFdr unsigned
#'
#' Kernel estimation of the density in a two-components mixture model
#' where one component are a standard Gaussian density.
#' Here we suppose that the density to estimate lives in \eqn{R^+}.
#'
#' @param X a vector of probit-transformed p-values (corresponding to a p-value serie)
#' @param p0 a priori proportion of \eqn{H_0} hypotheses
#' @param plotting boolean, should some diagnostic graphs be plotted. (Default is FALSE.)
#' @param NbKnot The (maximum) number of knot for the \code{kde} procedure. (Default is 1e5.)
#' @param tol a tolerance value for convergence. (Default is 1e-5.)
#' @param max_iter the maximum number of iterations allowed for the algorithm to converge or complete its process.(Default is 1e4.)
#' @import ks graphics stats
#'
#' @return A list with the following elements:
#' \tabular{ll}{
#' \code{p0} \tab vector of  the estimated proportions of \eqn{H_0} hypotheses
#' for each of p-value serie. \cr
#' \code{tau} \tab the vector of \eqn{H_1} posteriors. \cr
#' \code{f1} \tab a numeric vector, each coordinate \eqn{i}
#' corresponding to the evaluation of the \eqn{H_1} density at point \eqn{x_i},
#' where \eqn{x_i} is the \eqn{i}th item in \code{X}. \cr
#' \code{F1} \tab  a numeric vector, each coordinate \eqn{i}
#' corresponding to the evaluation of the \eqn{H_1} cdf at point \eqn{x_i},
#' where \eqn{x_i} is the \eqn{i}th item in \code{X}.
#' }

FastKerFdr_unsigned <- function(X, p0 = NULL, plotting = FALSE, NbKnot = 1e5, tol = 1e-5, max_iter = 1e4) {
  n <- length(X)

  ## Get a p0 estimate
  if (is.null(p0)) {
    p0 <- min(2 * sum(X < 0) / n, 1 - 1 / n)
  }
  p1 <- 1 - p0

  ## Knots, counts and initialization
  if (length(X) > NbKnot) {
    Hist <- hist(X, breaks = NbKnot, plot = FALSE)
    Knots <- Hist$mids
    ActualNbKnot <- length(Knots)
    Counts <- Hist$counts
  } else {
    Knots <- X
    ActualNbKnot <- length(X)
    Counts <- rep(1, ActualNbKnot)
  }

  if (plotting) {
    Order <- order(Knots)
    Knots <- Knots[Order]
    Counts <- Counts[Order]
  }

  h <- ks::hpi(X)

  ## Initialize the taus using GM
  phi <- dnorm(Knots)
  MaxKnot <- max(abs(range(Knots)))
  tau <- 0.8 * (purrr::map_dbl(Knots, ~ max(0, .x))) / MaxKnot + 0.1

  ## Get the weighted kernel density estimate
  diff <- 2 * tol
  iter <- 0
  while (diff > tol & iter <= max_iter) {
    iter <- iter + 1
    weights <- tau * Counts
    weights <- ActualNbKnot * weights / sum(weights)

    f1 <- ks::kde(x = Knots, w = weights, eval.points = Knots, h = h, density = TRUE)$estimate
    f1[f1 < 0] <- 0

    tauNew <- p1 * f1 / (p0 * phi + p1 * f1)

    ## Dirty job 1: get rid of the f1 mass on the left
    tauNew[Knots < -3] <- 0
    diff <- max(abs(tau - tauNew))
    tau <- tauNew
  }
  if (iter > max_iter & diff > tol) {
    message(paste0("Warning: The algorithm did not converge within max_iter=", max_iter, "."))
  }
  if (plotting) {
    Hist.fig <- hist(X,
      freq = TRUE, breaks = sqrt(n), main = "", border = 8,
      xlab = "Q-transformed pvalues", ylab = "Densities"
    )
    bin.width <- mean(diff(Hist.fig$breaks))
    lines(Knots, n * bin.width * p0 * phi, type = "l", col = 4, lwd = 2)
    lines(Knots, n * bin.width * p1 * f1, col = 2, lwd = 2)
    lines(Knots, n * bin.width * (p0 * phi + p1 * f1), lwd = 2)
    legend("topright",
      legend = c("H0 dist", "H1 dist", "Mixture Dist"),
      col = c("blue", "red", "black"), lty = c(1, 1, 2), cex = 0.8
    )
  }

  ## Now get the f1 estimate
  KDE <- ks::kde(x = Knots, w = weights, density = TRUE)
  f1 <- ks::dkde(x = X, fhat = KDE)
  Fhat <- ks::kcde(x = Knots, w = weights, h = KDE$h)
  F1 <- predict(Fhat, x = X)


  ## Dirty job 2: get rid of numeric problems
  f1[f1 < 0] <- 1e-30

  ## Get better estimate of p0
  f0 <- dnorm(X)
  tau0 <- p0 * f0 / (p1 * f1 + p0 * f0)
  diff <- abs(p0 - sum(tau0) / n)
  iter <- 0
  while (diff != 0 & iter < 1000) {
    iter <- iter + 1
    p0 <- sum(tau0) / n
    tau0 <- p0 * f0 / ((1 - p0) * f1 + p0 * f0)
    diff <- abs(p0 - sum(tau0) / n)
  }

  p1 <- 1 - p0
  tau1 <- p1 * f1 / (p1 * f1 + p0 * f0)

  return(list(p0 = p0, tau = tau1, f1 = f1, F1 = F1))
}


###################################################################
#' FastKerFdr signed
#'
#' Kernel estimation of the density in a two-components mixture model
#' where one component are a standard Gaussian density.
#'
#' @param X a vector of probit-transformed p-values (corresponding to a p-value serie).
#' @param p0 a priori proportion of \eqn{H_0} hypotheses.
#' @param plotting boolean, should some diagnostic graphs be plotted. (Default is FALSE.)
#' @param NbKnot The (maximum) number of knot for the \code{kde} procedure. (Default is 1e5.)
#' @param tol a tolerance value for convergence. (Default is 1e-5.)
#' @param max_iter the maximum number of iterations allowed for the algorithm to converge or complete its process.(Default is 1e4.)
#' @import ks graphics stats
#'
#' @return A list with the following elements:
#' \tabular{ll}{
#' \code{p0} \tab vector of  the estimated proportions of \eqn{H_0} hypotheses
#' for each of p-value serie. \cr
#' \code{tau} \tab the vector of \eqn{H_1} posteriors. \cr
#' \code{f1} \tab a numeric vector, each coordinate \eqn{i}
#' corresponding to the evaluation of the \eqn{H_1} density at point \eqn{x_i},
#' where \eqn{x_i} is the \eqn{i}th item in \code{X}. \cr
#' \code{F1} \tab  a numeric vector, each coordinate \eqn{i}
#' corresponding to the evaluation of the \eqn{H_1} cdf at point \eqn{x_i},
#' where \eqn{x_i} is the \eqn{i}th item in \code{X}.
#' }

FastKerFdr_signed <- function(X, p0 = NULL, plotting = FALSE, NbKnot = 1e5, tol = 1e-5, max_iter = 1e4) {
  X <- as.numeric(as.matrix(X))
  n <- length(X)

  ## Get a p0 estimate
  if (is.null(p0)) {
    p0 <- min(2 * sum(X < qnorm(0.5)) / n, 1 - 1 / n)
  }
  p1 <- 1 - p0

  ## Knots and counts
  if (length(X) > NbKnot) {
    Hist <- hist(X, breaks = NbKnot, plot = FALSE)
    Knots <- Hist$mids
    ActualNbKnot <- length(Knots)
    Counts <- Hist$counts
  } else {
    Knots <- X
    ActualNbKnot <- length(X)
    Counts <- rep(1, ActualNbKnot)
  }
  if (plotting) {
    Order <- order(Knots)
    Knots <- Knots[Order]
    Counts <- Counts[Order]
  }
  h <- ks::hpi(X)

  ## Initialize the taus
  phi <- dnorm(Knots)
  MaxKnot <- max(abs(range(Knots)))
  tau <- sign(Knots) * 0.8 * Knots / MaxKnot + 0.1

  ## Get the weighted kernel density estimate
  diff <- 2 * tol
  iter <- 0
  while (diff > tol & iter <= max_iter) {
    iter <- iter + 1
    weights <- tau * Counts
    weights <- ActualNbKnot * weights / sum(weights)
    f1 <- ks::kde(x = Knots, w = weights, eval.points = Knots, h = h, density = TRUE)$estimate
    f1[f1 < 0] <- 0
    tauNew <- p1 * f1 / (p0 * phi + p1 * f1)
    diff <- max(abs(tau - tauNew))
    tau <- tauNew
  }
  if (iter > max_iter & diff > tol) {
    message(paste0("Warning: The algorithm did not converge within max_iter=", max_iter, "."))
  }
  if (plotting) {
    Hist.fig <- hist(X,
      freq = TRUE, breaks = sqrt(n), main = "", border = 8,
      xlab = "Q-transformed pvalues", ylab = "Densities"
    )
    bin.width <- mean(diff(Hist.fig$breaks))
    lines(Knots, n * bin.width * p0 * phi, type = "l", col = 4, lwd = 2)
    lines(Knots, n * bin.width * p1 * f1, col = 2, lwd = 2)
    lines(Knots, n * bin.width * (p0 * phi + p1 * f1), lwd = 2)
    legend("topright",
      legend = c("H0 dist", "H1 dist", "Mixture Dist"),
      col = c("blue", "red", "black"), lty = c(1, 1, 2), cex = 0.8
    )
  }

  ## Now get the f1 estimate
  KDE <- ks::kde(x = Knots, w = weights, density = TRUE)
  f1 <- ks::dkde(x = X, fhat = KDE)
  Fhat <- ks::kcde(x = Knots, w = weights, h = KDE$h)
  F1 <- predict(Fhat, x = X)

  ## Dirty job 2: get rid of numeric problems
  f1[f1 < 0] <- 1e-30

  ## Get better estimate of p0
  f0 <- dnorm(X)
  tau0 <- p0 * f0 / (p1 * f1 + p0 * f0)
  diff <- abs(p0 - sum(tau0) / n)
  iter <- 0
  while (diff != 0 & iter < 1000) {
    iter <- iter + 1
    p0 <- sum(tau0) / n
    tau0 <- p0 * f0 / ((1 - p0) * f1 + p0 * f0)
    diff <- abs(p0 - sum(tau0) / n)
  }

  p1 <- 1 - p0
  tau1 <- p1 * f1 / (p1 * f1 + p0 * f0)

  return(list(p0 = p0, tau = tau1, f1 = f1, F1 = F1))
}

###################################################################
#' Signed case function: Separate f1 into f+ and f-
#'
#' @param XMat a matrix of probit-transformed p-values, each column corresponding to a p-value serie.
#' @param f0Mat a matrix containing the evaluation of the marginal density functions under \eqn{H_0} at each items,
#' each column corresponding to a p-value serie.
#' @param f1Mat a matrix containing the evaluation of the marginal density functions under \eqn{H_1} at each items,
#' each column corresponding to a p-value serie.
#' @param p0 the proportions of \eqn{H_0} items for each series.
#' @param plotting boolean, should some diagnostic graphs be plotted. (Default is FALSE.)
#'
#' @importFrom dplyr %>%
#'
#' @return A list with the following elements:
#' \tabular{ll}{
#' \code{f1plusMat} \tab a matrix containing the evaluation of the marginal density functions under \eqn{H_1^+}
#' at each items, each column corresponding to a p-value serie. \cr
#' \code{f1minusMat} \tab a matrix containing the evaluation of the marginal density functions under \eqn{H_1^-}
#' at each items, each column corresponding to a p-value serie. \cr
#' \code{p1plus} \tab an estimate of the proportions of \eqn{H_1^+} items for each series. \cr
#' \code{p1minus} \tab an estimate of the proportions of \eqn{H_1^-} items for each series.
#' }

f1_separation_signed <- function(XMat, f0Mat, f1Mat, p0, plotting = FALSE) {
  n <- nrow(XMat)
  Q <- ncol(XMat)

  ### Compute of the integral for normalization
  XMat_temp <- apply(XMat, 2, sort)
  XMat_tempm1 <- matrix(1, n, Q)
  f1Mat_temp <- matrix(1, n, Q)
  for (q in 1:Q) {
    f1Mat_temp[, q] <- f1Mat[order(XMat[, q]), q]
    XMat_tempm1[, q] <- c(0, XMat_temp[2:n - 1, q])
  }
  Xplustemp <- XMat_temp > 0
  lambda_plus <- colSums((XMat_temp[-1, ] - XMat_tempm1[-1, ]) * (f1Mat_temp[-1, ] * Xplustemp[-1, ]))
  lambda_minus <- 1 - lambda_plus

  ### Compute the distributions f+ and f-
  Xplus <- XMat > 0
  fplusMat <- (f1Mat * Xplus) %>% sweep(2, lambda_plus, `/`)
  fminusMat <- (f1Mat * abs(Xplus - 1)) %>% sweep(2, lambda_minus, `/`)

  ### Compute the proportions p+ and p-
  pplus <- (1 - p0) * lambda_plus
  pminus <- 1 - p0 - pplus

  if (plotting) {
    for (q in 1:Q) {
      Knots <- XMat[, q]
      ActualNbKnot <- length(XMat[, q])
      Counts <- rep(1, ActualNbKnot)

      Order <- order(Knots)
      Knots <- Knots[Order]
      Counts <- Counts[Order]

      f0 <- f0Mat[Order, q]
      fplus <- fplusMat[Order, q]
      fminus <- fminusMat[Order, q]
      f1 <- f1Mat[Order, q]


      Hist.fig <- hist(XMat[, q],
        freq = TRUE, breaks = sqrt(n), main = "Marginal distributions inferred", border = 8,
        xlab = "Q-transformed pvalues", ylab = "Densities"
      )
      bin.width <- mean(diff(Hist.fig$breaks))
      lines(Knots, n * bin.width * p0[q] * f0, type = "l", col = 4, lwd = 2)
      lines(Knots, n * bin.width * pplus[q] * fplus, col = 2, lwd = 2)
      lines(Knots, n * bin.width * pminus[q] * fminus, col = 3, lwd = 2)
      lines(Knots, n * bin.width * (p0[q] * f0 + pplus[q] * fplus + pminus[q] * fminus), lwd = 2, lty = 2)
      legend("topright",
        legend = c("f0 dist", "f+ dist", "f- dist", "Mixture dist"),
        col = c(4, 2, 3, 1), lty = c(1, 1, 1, 2), cex = 0.8
      )
    }
  }
  return(list(f1plusMat = fplusMat, f1minusMat = fminusMat, p1plus = pplus, p1minus = pminus))
}
