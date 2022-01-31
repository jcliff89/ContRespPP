#' Example Continuous Response ANOVA Dataset.
#'
#' A dataset containing an example test set result to run a predictive probability model
#' with a continuous response variable.
#' This data set is a result of a 2^4 full factorial with five replicates that included five
#' main effects, with two-way interactions (excluding two-way interactions for the third factor).
#' This design has a power of 80%, with an 80% confidence level, to detect a difference of 50
#' units with a standard deviation of 100 units.
#'
#' Please refer to the vignette for an example analysis of this data using ContRespPP.
#'
#' @format A matrix with 80 rows and 14 columns:
#' \describe{
#'   \item{Column 1}{Continuous Response Variable}
#'   \item{Columns 2-14}{Design Matrix}
#' }
#' @source Sieck VRC, Christensen FGW. A framework for improving the efficiency of operational testing through Bayesian adaptive design. Quality and Reliability Engineering International. 2021; 3018-3033.
"exData"

