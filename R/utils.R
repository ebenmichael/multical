################################################################################
## Utilities to examine the weights
################################################################################

#' Compute the imbalance for a given design matrix
#' @param D Design matrix, output from \code{\link{create_design_matrix}}
#' @param weights Vector of weights
#' @param sample_counts Vector of sample counts for each cell
#' @param target_counts Vector of target counts for each cell
#' 
#' @keywords internal
compute_balance <- function(D, weights, sample_counts, target_counts) {
  as.numeric(Matrix::t(D) %*% (weights * sample_counts - target_counts))
}

#' Get the difference between the reweighted sample and the target
#'
#' For each interaction term, computes the difference between the reweighted 
#' sample count and the target count
#'
#' @param output Output of \code{\link{multical}}
#' @param order Integer. Order of interactions to include in the balance measure
#' 
#' @return A data frame with the difference for each term and value of the 
#' hyper-parameter
get_balance <- function(output, order) {

  # get distinct cells
  cells <- output %>% select(-sample_count, -target_count, -lambda, -weight) %>%
    distinct()
  
  D <- Matrix::sparse.model.matrix(~ . - 1, cells)
  if(order > 1) {
    D <- cbind(D, create_design_matrix(cells, order))
  }

  output %>%
    nest(data = colnames(output)[colnames(output) != "lambda"]) %>%
    mutate(imbalance = lapply(data,
      function(df) {
        data.frame(
          term = colnames(D),
          difference = compute_balance(D, df$weight, df$sample_count,
                                       df$target_count)
        )
      })) %>%
    unnest(imbalance) %>%
    select(-data)
}


#' Compute the overall level of balance versus the effective sample size
#' 
#' For each value of the hyper-parameter, computes the effective sample size 
#' and the overall root sum of squared differences in all interactions,
#' normalized by the size of the target population
#' 
#' @inheritParams get_balance
#' 
#' @return A data frame with the imbalance and effective sample size for each
#' value of the hyper-parameter
get_balance_v_sample_size <- function(output, order) {

  target_pop <- sum(output$target_count)
  # get imbalances
  imbals <- get_balance(output, order) %>%
    group_by(lambda) %>%
    summarise(imbalance = sqrt(sum(difference ^ 2) / target_pop)) %>%
    ungroup()

  # get effective sample sizes
  neff <- output %>%
    group_by(lambda) %>%
    summarise(n_eff = sum(weight * sample_count) ^ 2 /
                      sum(weight ^ 2 * sample_count)) %>%
    ungroup()
  
  return(inner_join(imbals, neff, by = "lambda"))
}

