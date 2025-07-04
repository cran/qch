#' Synthetic example to illustrate the main qch functions
#'
#' PvalSets is a data.frame with 10,000 rows and 3 columns. Each row corresponds to an item,
#' columns 'Pval1' and 'Pval2'  each correspond to a test serie over the items, and column 'Class'
#' provides the truth, i.e. if item \eqn{i} belongs to class 1 then the H0 hypothesis is true for the 2 tests,
#' if item \eqn{i} belongs to class 2 (resp. 3) then the H0 hypothesis is true for the first (resp. second)
#' test only, and if item \eqn{i} belongs to class 4 then both H0 hypotheses are false (for the first
#' and the second test).
#'
#' @format A data.frame
#'
"PvalSets"



#' Synthetic example to illustrate the main qch functions using Gaussian copula
#'
#' PvalSets_cor is a data.frame with 10,000 rows and 3 columns. Each row corresponds to an item,
#' columns \code{Pval1} and \code{Pval2}  each correspond to a test serie over the items, and column 'Class'
#' provides the truth, i.e. if item \eqn{i} belongs to class 1 then the \eqn{H_0} hypothesis is true for the 2 tests,
#' if item \eqn{i} belongs to class 2 (resp. 3) then the \eqn{H_0} hypothesis is true for the first (resp. second)
#' test only, and if item \eqn{i} belongs to class 4 then both H0 hypotheses are false (for the first
#' and the second test). The correlation between the two pvalues series within each class is 0.3.
#'
#' @format A data.frame
#'
"PvalSets_cor"
