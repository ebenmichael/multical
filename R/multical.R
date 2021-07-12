################################################################################
## Multilevel calibration
################################################################################

#' multical
#' 
#' @description A package for multilevel calibration weighting
#' @docType package
#' @name multical-package
#' @importFrom magrittr "%>%"
#' @import dplyr
#' @import tidyr
#' @importFrom stats terms
#' @importFrom stats formula
NULL



#' Calibrate sample to target via multilevel calibration
#'
#' Finds weights that exactly calibrate first order margins between respondents and the target population. Requires individual-level of cell-level data.
#'
#' @param formula Formula of the form \code{sample_count ~ covariates}, where 
#' \code{sample_count} is whether or not an individual responded 
#' (with individual-level data) or the number of respondents in the cell 
#' (with cell-level data) and \code{covariates} are the covariates to calibrate 
#' on
#' @param target_count Name of column with indicators for whether an individual 
#' is in the target population (with individual-level data) or the target counts for each cell (with cell-level data)
#' @param data Dataframe with covariate information, sample and target counts
#' @param order Integer. What order interactions to balance
#' @param lambda Numeric. Regularization hyperparamter, by default fits weights 
#' for a range of values
#' @param lambda_max Numeric. Maximum hyperparameter to fit weights with, 
#' default is the root sum of squared differences between the (unweighted) sample and the target
#' @param n_lambda Integer. Number of hyper-parameters to fit weights for, from 
#' lambda_max to lambda_max * 1e-5, equally spaced on the log scale. Default, 20
#' @param lowlim Lower bound on weights, default 0
#' @param uplim Upper bound on weights, default Inf
#' @param verbose Boolean. Show optimization information, default False
#' @param ... Additional parameters for osqp
#' 
#' @return data frame with the weight for each distinct cell, for each value of 
#' the hyperparameter \code{lambda}
#' @export
multical <- function(formula, target_count, data,
                      order = NULL, lambda = NULL,
                      lambda_max = NULL, n_lambda = 20,
                      lowlim = 0, uplim = Inf,
                      verbose = FALSE, ...) {
  
  # create distinct cells for all interactions
  if(verbose) message("Creating table of cell counts")
  cells <- create_cells(formula, enquo(target_count), data)

  # get weights
  weights <- calibrate_(cells %>% select(-sample_count, -target_count),
                        cells$sample_count, cells$target_count,
                        order, lambda, lambda_max, n_lambda,
                        lowlim, uplim, verbose, ...)
  # combine back in and return
  cells %>% filter(sample_count != 0) %>%
    bind_cols(weights) %>%
    select(-sample_count, -target_count) %>%
    right_join(cells,
               by = cells %>% select(-sample_count, -target_count) %>%
                    names()) %>%
    pivot_longer(!names(cells), names_to = "lambda", values_to = "weight") %>%
    mutate(weight = replace_na(weight, 0),
           lambda = as.numeric(lambda)) %>%
    return()
}


#' Create data frame with all distinct cells, and their sample and target counts
#' @inheritParams multical
#'
#' @keywords internal
create_cells <- function(formula, target_count, data) {
  covs <- all.vars(terms(Formula::Formula(formula), rhs = 1)[[3]])
  sample_count <- terms(Formula::Formula(formula), rhs = 1)[[2]]
  cells <- data %>%
    mutate(across(covs, as.factor)) %>%
    group_by(across(covs)) %>%
    summarise(sample_count = sum(!!sample_count),
              target_count = sum(!!target_count)) %>%
    ungroup()

  return(cells)
}


#' Internal function to perform multilevel calibration
#' @param cells Dataframe of distinct cells
#' @param sample_counts Vector of sample counts for each cell
#' @param target_counts Vector of target counts for each cell
#' @inheritParams multical
#'
#' @keywords internal
calibrate_ <- function(cells, sample_counts, target_counts, order = NULL,
                      lambda = 1, lambda_max = NULL, n_lambda = 100,
                      lowlim = 0, uplim = Inf, verbose = FALSE,
                      ...) {

  if(verbose) message("Creating design matrix")
  # get design matrix for the number of interactions
  D <- create_design_matrix(cells, order)

  # create constraints for raking
  constraints <- create_rake_constraints(cells, D, sample_counts,
                                         target_counts, lowlim, uplim,
                                         verbose)


  # P matrix and q vector
  if(verbose) message("Creating quadratic term matrix")
  if(is.null(lambda)) {
    if(is.null(lambda_max)) {
      unif_imbal <- Matrix::t(D) %*% (sample_counts - target_counts)
      lambda_max <- sqrt(sum(unif_imbal ^ 2))
    }
    lam_seq <- lambda_max * 10 ^ seq(0, -5, length.out = n_lambda)
    P <- create_rake_Pmat(D, sample_counts, 0)
  } else {
    P <- create_rake_Pmat(D, sample_counts, lambda)
  }
  
  
  if(verbose) message("Creating linear term vector")
  qvec <- create_rake_qvec(D, sample_counts, target_counts, lambda)


  settings <- do.call(osqp::osqpSettings,
                        c(list(...), list(verbose = verbose)))
  if(verbose) message("Setting up optimization")
  solver <- osqp::osqp(P, qvec, constraints$A, constraints$l, constraints$u,
                       pars = settings)
  if(verbose) message("Optimizing")

  if(is.null(lambda)) {
      
    x <- NULL
    y <- NULL
    solution <- matrix(, nrow = nrow(P), ncol = n_lambda)
    # solve for values of lambda
    for(i in 1:length(lam_seq)) {
      if(verbose) message(paste0("Solving with lambda = ", lam_seq[i]))

      P_new <- update_rake_Pmat(P, sample_counts, lam_seq[i])
      solver$Update(Px = Matrix::triu(P_new)@x)

      # use warm start with previous solution to speed up optimization
      solver$WarmStart(x = x, y = y)
      sol_lam <- solver$Solve()
      x <- sol_lam$x
      y <- sol_lam$y
      solution[, i] <- x
    }
    solution <- as.data.frame(solution)
    names(solution) <- lam_seq
  } else {
    solution <- data.frame(solver$Solve()$x)
    names(solution) <- lambda
  }
  return(solution)

}

#' Creates the underlying design matrix to solve the multilevel calibration
#' optimization problem
#' 
#' @param cells Dataframe of distinct cells
#' @param order Integer. What order interactions to balance
#'
#' @keywords internal
create_design_matrix <- function(cells, order) {

  if(is.null(order)) {
    order <- ncol(cells)
  }
  if(order == 1) {
    return(Matrix::sparseMatrix(dims = c(nrow(cells), 1), i = {}, j = {}))
  }

  contrast_cells <- lapply(cells, contrasts, contrasts = F)

  form <- as.formula(paste("~ . ^ ", order, " - . + 0"))

  D <- Matrix::sparse.model.matrix(form, data = cells)
  return(D)
}


#' Create underlying constraints for optimization problem
#' @param cells Dataframe of distinct cells
#' @param D Design matrix, output from \code{\link{create_design_matrix}}
#' @param sample_counts Vector of sample counts for each cell
#' @param target_counts Vector of target counts for each cell
#' @inheritParams multical
#'
#' @keywords internal
create_rake_constraints <- function(cells, D, sample_counts, target_counts,
                                    lowlim, uplim, verbose) {

  if(verbose) message("Creating constraint matrix")

  # number of non-empty cells
  n_nempty <- sum(sample_counts != 0)
  sample_counts_nempty <- sample_counts[sample_counts != 0]

  # sum to N constraint
  if(verbose) message("\tx Sum to N constraint")

  A_sumN <- Matrix::Matrix(sample_counts_nempty, nrow = 1, ncol = n_nempty)
  l_sumN <- sum(target_counts)
  u_sumN <- l_sumN


  # non-negative constraint
  if(verbose) message("\tx Non-negativity constraint")
  A_nn <- Matrix::Diagonal(n_nempty)
  l_nn <- rep(lowlim, n_nempty)
  u_nn <- rep(uplim, n_nempty)

  # marginal constraints
  if(verbose) message("\tx Exact marginal balance constraints")
  design_mat <- Matrix::sparse.model.matrix(~ . - 1, cells)

  A_marg <- Matrix::t(design_mat[sample_counts != 0, , drop = F] * 
                      sample_counts_nempty)
  l_marg <- as.numeric(Matrix::t(design_mat) %*% target_counts)
  u_marg <- l_marg

  # combine
  A <- rbind(A_sumN, A_nn, A_marg)
  l <- c(l_sumN, l_nn, l_marg)
  u <- c(u_sumN, u_nn, u_marg)

  return(list(A = A, l = l, u = u))

}

#' Create matrix for quadratic term in QP
#' @param D Design matrix, output from \code{\link{create_design_matrix}}
#' @param sample_counts Vector of sample counts for each cell
#' @param lambda Numeric. Hyperparameter
#'
#' @keywords internal
create_rake_Pmat <- function(D, sample_counts, lambda) {


  nnz_idxs <- which(sample_counts != 0)
  sqrtP <- D[nnz_idxs,, drop = F] * sample_counts[nnz_idxs]
  P <- sqrtP %*% Matrix::t(sqrtP) +
    lambda * sample_counts[nnz_idxs] * Matrix::Diagonal(length(nnz_idxs))
  return(P)

}

#' Update the quadratic term matrix with a new hyperparameter value
#' @param P Original quadratic term matrix
#' @inheritParams create_rake_Pmat
#'
#' @keywords internal
update_rake_Pmat <- function(P, sample_counts, lambda) {

  nnz_idxs <- which(sample_counts != 0)
  P <- P +
    lambda * sample_counts[nnz_idxs] * Matrix::Diagonal(length(nnz_idxs))
  return(P)
}



#' Create linear term vector in QP
#' @inheritParams create_rake_constraints
#'
#' @keywords internal
create_rake_qvec <- function(D, sample_counts, target_counts, lambda) {

  nnz_idxs <- which(sample_counts != 0)
  rhs <- Matrix::t(D) %*% target_counts
  lhs <- D[nnz_idxs,, drop = F] * sample_counts[nnz_idxs]
  return(-c(as.numeric(lhs %*% rhs)))
}