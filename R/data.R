################################################################################
## Small toy data for examples
################################################################################


#' Fake data at the individual level
#' 
#' @format A data frame with 2500 rows and 7 variables:
#' \describe{
#'    \item{X1, X2, X3, X4}{Categorical variables}
#'    \item{response}{Whether the individual responded to the survey}
#'    \item{intarget}{Whether the individual is in the target population (all 1)}
#'    \item{y}{Binary answer to survey question}
#' }
"data_individual"


#' Fake data at the cell level 
#' 
#' @format A data frame with 91 rows and 6 variables:
#' \describe{
#'    \item{X1, X2, X3, X4}{Categorical variables}
#'    \item{sample_count}{Number of respondents in the cell}
#'    \item{target_count}{Total population size in the cell}
#'    \item{y}{Proportion that answered "1" to a binary survey question in the cell}
#' }
"data_cell"