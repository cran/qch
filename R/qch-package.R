#' @keywords internal
#' 
#' @details The main functions of the package \code{\link{GetHconfig}}, \code{\link{GetH1AtLeast}},
#' \code{\link{GetH1Equal}},
#' \code{\link{qch.fit}} and \code{\link{qch.test}} correspond to the 
#' 4 steps for querying a composite hypothesis:
#' \itemize{
#' \item Building all possible combination of simple hypotheses \eqn{H_0}/\eqn{H_1}
#' \item Composite alternative hypothesis formulation
#' \item Inferring the null distribution
#' \item Testing the composite null hypothesis
#' }
#' 
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib qch, .registration = TRUE
## usethis namespace: end
NULL
