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
#' @importFrom stats as.formula
#' @importFrom stats weights
#' @importFrom stats lm.fit
#' @importFrom stats model.matrix
#' @importFrom stats predict
#' @importFrom glmnet cv.glmnet
#' @importFrom rlang .data
#' @importFrom rlang eval_tidy
#' @importFrom rlang enquo
#' @importFrom xgboost xgboost
#' @importFrom xgboost xgb.cv
NULL



#' Calibrate sample to target via multilevel calibration
#'
#' Finds weights that exactly calibrate first order margins between respondents
#' and the target population, operating at the individual level.
#'
#' @param formula Right-hand-side formula of the form \code{~ covariates}
#' specifying the covariates to calibrate on
#' @param sample_data Individual-level data frame of respondents. One row per
#' respondent; each row receives sample weight 1 in the optimization.
#' @param pop_data Data frame describing the target population. Can be
#' individual-level (one row per population member) or pre-aggregated to
#' the cell level. Cells present in \code{pop_data} but absent from
#' \code{sample_data} are handled automatically.
#' @param target_count Optional. Name of column in \code{pop_data} giving the
#' count or weight for each row. Defaults to 1 per row (row-count semantics).
#' @param order Integer. What order interactions to balance. Default is all orders
#' @param lambda Numeric. Regularization hyperparamter, by default fits weights 
#' for a range of values
#' @param lambda_max Numeric. Maximum hyperparameter to fit weights with, 
#' default is the root sum of squared differences between the (unweighted) sample and the target
#' @param n_lambda Integer. Number of hyper-parameters to fit weights for, from 
#' lambda_max to lambda_max * lambda_min_ratio, equally spaced on the log scale. Default, 20
#' @param lambda_min_ratio Numeric. Ratio of min to max lambda to consider.
#' @param lowlim Lower bound on weights, default 0
#' @param uplim Upper bound on weights, default Inf
#' @param base_weights Optional. Name of column in \code{sample_data} giving
#' individual-level base weights to regularize towards. When supplied, the
#' regularization penalty becomes
#' \eqn{\lambda \sum_i s_i (w_i - b_i)^2 / 2} instead of
#' \eqn{\lambda \sum_i s_i w_i^2 / 2}. Defaults to 1 for all respondents.
#' @param exact_order Integer. Order of interactions to calibrate exactly
#'   (i.e. constrain the weighted sample to match the target). Default 1
#'   (first-order margins only). Setting to 2 also exactly matches all
#'   two-way interaction margins, etc. Must be between 1 and \code{order}.
#' @param scale_by_order Logical. If \code{TRUE} (default), each column of the
#'   interaction design matrix is scaled by \eqn{1/\sqrt{n_K}}, where
#'   \eqn{n_K} is the total number of design-matrix columns corresponding to
#'   order-\eqn{K} interactions. This makes the objective equivalent to summing
#'   over orders and then averaging across moment conditions within each order,
#'   rather than summing all moment conditions equally. Set to \code{FALSE} to
#'   remove scaling
#' @param verbose Boolean. Show optimization information, default FALSE
#' @param ... Additional parameters for osqp
#' 
#' @return A \code{multical} object with the following fields:
#' \itemize{
#'   \item \code{weights}: numeric matrix (n_respondents x n_lambdas) of
#'     calibration weights, column names are lambda values
#'   \item \code{lambda}: numeric vector of lambda values
#'   \item \code{cells}: data frame with covariate columns plus
#'     \code{sample_count}, \code{target_count}, and \code{base_weight};
#'     includes pop-only rows (\code{sample_count = 0})
#'   \item \code{formula}: the formula passed to \code{multical}
#'   \item \code{order}: resolved order of interactions used
#'   \item \code{exact_order}: order of interactions constrained exactly
#'   \item \code{n_respondents}: number of respondents
#'   \item \code{default_lambda_idx}: index of the auto-selected best lambda
#'   \item \code{balance_threshold}: threshold used for lambda selection
#'   \item \code{prop_uncovered}: data frame with columns \code{order} (integer,
#'     1 through the calibration order) and \code{prop_uncovered} (numeric in
#'     \eqn{[0, 1]}), giving the proportion of total target count in cells that have
#'     positive target counts but no sample observations, for each interaction
#'     order
#' }
#' Use \code{\link{weights}} to extract respondent weights at the default
#' lambda, and \code{\link{get_balance}} to assess calibration quality.
#'
#' @export
multical <- function(formula, sample_data, pop_data, target_count = NULL,
                      order = NULL, exact_order = 1, lambda = NULL,
                      lambda_max = NULL, n_lambda = 20,
                      lambda_min_ratio = 1e-10,
                      lowlim = 0, uplim = Inf,
                      base_weights = NULL,
                      scale_by_order = TRUE,
                      verbose = FALSE, ...) {

  # ungroup both inputs if grouped
  if(is_grouped_df(sample_data)) sample_data <- ungroup(sample_data)
  if(is_grouped_df(pop_data))    pop_data    <- ungroup(pop_data)

  # extract base weights from sample_data (must be done before create_units_sep
  # to preserve 1:1 row correspondence with sample_data)
  bw_quo <- enquo(base_weights)
  tc_quo <- enquo(target_count)
  if(rlang::quo_is_null(bw_quo)) {
    bw_vec <- rep(1, nrow(sample_data))
  } else {
    bw_vec <- sample_data %>% pull(!!bw_quo)
  }

  # rescale bw_vec so that it sums to the total target count
  if(rlang::quo_is_null(tc_quo)) {
    total_target <- nrow(pop_data)
  } else {
    total_target <- pop_data %>%
      summarise(total = sum(!!tc_quo), .groups = "drop") %>%
      pull("total")
  }
  # bw_vec <- bw_vec / sum(bw_vec) * total_target

  # build unit-level data frame: respondent rows + pop-only rows
  if(verbose) message("Creating table of unit counts")
  units <- create_units_sep(formula, sample_data, pop_data, tc_quo)
  units <- units %>% mutate(.row_id = row_number())

  # pad bw_vec to the full length of `units` (respondents + pop-only rows).
  # Pop-only rows have sample_count = 0, so their base weight value is never
  # used in the objective, but the vectors must be the same length.
  n_pop_only <- sum(units$sample_count == 0)
  bw_vec_full <- c(bw_vec, rep(1, n_pop_only))

  # attach base_weight as a column so it appears in the output (NA for pop-only)
  units <- units %>%
    mutate(base_weight = c(bw_vec, rep(NA_real_, n_pop_only)))

  # resolve order now so it can be stored in the multical object
  if (is.null(order)) order <- length(all.vars(formula))

  # get weights (one per respondent, i.e. per unit with sample_count != 0)
  weights <- calibrate_(units %>% select(-c("sample_count", "target_count",
                                             "base_weight", ".row_id")),
                        units$sample_count,
                        units$target_count,
                        order, exact_order, lambda, lambda_max, n_lambda,
                        lambda_min_ratio,
                        lowlim, uplim,
                        bw_vec_full, scale_by_order, verbose, ...)

  # extract the weights matrix (respondents only, rows in sample_data order)
  weights_matrix <- as.matrix(weights)
  lam_vec <- as.numeric(colnames(weights_matrix))

  # cells = units without the internal row-id helper column;
  # pop-only rows (sample_count = 0) are retained so that get_balance() can
  # compute correct marginal targets.
  cells <- units %>% select(-".row_id")

  new_multical(weights_matrix, lam_vec, cells, formula, order, exact_order,
               scale_by_order)
}


#' Construct a multical object
#'
#' @param weights_matrix Numeric matrix (n_respondents x n_lambdas) of
#'   calibration weights. Column names are the lambda values.
#' @param lambda Numeric vector of lambda values.
#' @param cells Data frame with covariate columns, \code{sample_count},
#'   \code{target_count}, and \code{base_weight}. One row per respondent
#'   followed by one row per pop-only cell.
#' @param formula The formula passed to \code{\link{multical}}.
#' @param order Integer. Resolved order of interactions used in the calibration.
#' @param exact_order Integer. Order of interactions constrained exactly in the
#'   optimization. Default 1.
#' @param scale_by_order Logical. Whether the interaction design matrix was
#'   scaled by \eqn{1/\sqrt{n_K}} within each order. Must match the value used
#'   during optimization so that balance is measured on the same scale.
#'   Default \code{TRUE}.
#' @param balance_threshold Numeric in (0, 1). Fraction of the total possible
#'   balance gain required to qualify a lambda for selection. Default 0.95.
#'
#' @keywords internal
new_multical <- function(weights_matrix, lambda, cells, formula, order,
                         exact_order = 1,
                         scale_by_order = TRUE,
                         balance_threshold = 0.95) {
  n_respondents <- sum(cells$sample_count != 0)
  default_lambda_idx <- select_default_lambda(
    weights_matrix, cells, lambda, order, balance_threshold, scale_by_order
  )
  structure(
    list(
      weights            = weights_matrix,
      lambda             = lambda,
      cells              = cells,
      formula            = formula,
      order              = order,
      exact_order        = exact_order,
      scale_by_order     = scale_by_order,
      n_respondents      = n_respondents,
      default_lambda_idx = default_lambda_idx,
      balance_threshold  = balance_threshold,
      prop_uncovered     = compute_prop_uncovered(cells, order)
    ),
    class = "multical"
  )
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


#' Create data frame with one row per unit (no aggregation), with sample and target indicators
#' @inheritParams multical
#'
#' @keywords internal
create_units <- function(formula, target_count, data) {
  covs <- all.vars(terms(Formula::Formula(formula), rhs = 1)[[3]])
  sample_count <- terms(Formula::Formula(formula), rhs = 1)[[2]]
  units <- data %>%
    mutate(across(all_of(covs), as.factor),
           sample_count = !!sample_count,
           target_count = !!target_count) %>%
    select(all_of(covs), sample_count, target_count)

  return(units)
}


#' Build the combined unit-level data frame from separate sample and population data
#'
#' Returns one row per respondent (\code{sample_count = 1},
#' \code{target_count = pop_count_cell / n_respondents_cell}) followed by one
#' row per population-only cell (\code{sample_count = 0},
#' \code{target_count = pop_count_cell}). Population-only rows have no
#' optimization variables but contribute to the RHS of the marginal constraints.
#'
#' @param formula Right-hand-side formula specifying covariates
#' @param sample_data Individual-level respondent data frame
#' @param pop_data Population data frame (individual or cell level)
#' @param target_count Quosure of the target count column in \code{pop_data},
#' or a null quosure to use row counts
#'
#' @keywords internal
create_units_sep <- function(formula, sample_data, pop_data, target_count) {

  covs <- all.vars(formula)

  # cast covariates to factors in both datasets
  pop_data    <- pop_data    %>% mutate(across(all_of(covs), as.factor))
  sample_data <- sample_data %>% mutate(across(all_of(covs), as.factor))

  # aggregate pop_data to one count per cell
  if(rlang::quo_is_null(target_count)) {
    pop_agg <- pop_data %>%
      group_by(across(all_of(covs))) %>%
      summarise(.pop_count = n(), .groups = "drop")
  } else {
    pop_agg <- pop_data %>%
      group_by(across(all_of(covs))) %>%
      summarise(.pop_count = sum(!!target_count), .groups = "drop")
  }

  # count respondents per cell
  sample_agg <- sample_data %>%
    group_by(across(all_of(covs))) %>%
    summarise(.n_respondents = n(), .groups = "drop")

  # respondent rows: sample_count = 1, target_count = pop_count / n_respondents
  # (cells absent from pop_data get target_count = 0)
  respondent_units <- sample_data %>%
    left_join(pop_agg,    by = covs) %>%
    left_join(sample_agg, by = covs) %>%
    mutate(sample_count = 1L,
           target_count = replace_na(.data$.pop_count, 0) / .data$.n_respondents) %>%
    select(all_of(covs), "sample_count", "target_count")

  # population-only rows: cells present in pop but absent from sample
  # sample_count = 0; they contribute to constraint RHS but have no weight variable
  pop_only <- pop_agg %>%
    anti_join(sample_agg, by = covs) %>%
    mutate(sample_count = 0L,
           target_count = .data$.pop_count) %>%
    select(all_of(covs), "sample_count", "target_count")

  bind_rows(respondent_units, pop_only)
}


#' Internal function to perform multilevel calibration
#' @param cells Dataframe of distinct cells
#' @param sample_counts Vector of sample counts for each cell
#' @param target_counts Vector of target counts for each cell
#' @inheritParams multical
#'
#' @keywords internal
calibrate_ <- function(cells, sample_counts, target_counts, order = NULL,
                      exact_order = 1,
                      lambda = 1, lambda_max = NULL, n_lambda = 100,
                      lambda_min_ratio = 1e-10,
                      lowlim = 0, uplim = Inf, base_weights = NULL,
                      scale_by_order = TRUE,
                      verbose = FALSE,
                      ...) {

  # default base weights to 1 (recovers standard regularization)
  if(is.null(base_weights)) base_weights <- rep(1, length(sample_counts))
  # rescale target_counts and base weights to sum to the sample size
  target_counts <- target_counts / sum(target_counts) * sum(sample_counts)
  base_weights <- base_weights / sum(base_weights) * sum(sample_counts)

  if(is.null(order)) {
    order <- ncol(cells)
  }
  if(verbose) message("Creating design matrix")
  # get design matrix for the number of interactions
  D <- create_design_matrix(cells, order)
  if(scale_by_order) D <- scale_design_matrix(D)

  # create constraints for raking
  constraints <- create_rake_constraints(cells, D, sample_counts,
                                         target_counts, lowlim, uplim,
                                         exact_order, verbose)


  # P matrix and q vector
  if(verbose) message("Creating quadratic term matrix")
  if(is.null(lambda) & order > 1) {
    if(is.null(lambda_max)) {
      unif_imbal <- Matrix::t(D) %*% (base_weights - target_counts)
      lambda_max <- sqrt(sum(unif_imbal ^ 2))
    }
    lam_seq <- lambda_max * 10 ^ seq(0, log10(lambda_min_ratio),
                                     length.out = n_lambda)
    P <- create_rake_Pmat(D, sample_counts, 0)
  } else if(is.null(lambda) & order == 1) {
    lambda <- 0
    P <- create_rake_Pmat(D, sample_counts, lambda)
  } else {
    P <- create_rake_Pmat(D, sample_counts, lambda)
  }
  
  
  if(verbose) message("Creating linear term vector")
  qvec <- create_rake_qvec(D, sample_counts, target_counts,
                           if(is.null(lambda)) 0 else lambda,
                           base_weights)

  
  P <- P
  qvec <- qvec
  settings <- do.call(osqp::osqpSettings,
                        c(list(...),
                          list(verbose = verbose)))
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
      q_new <- create_rake_qvec(D, sample_counts, target_counts, lam_seq[i], base_weights)
      solver$Update(Px = Matrix::triu(P_new)@x, q = q_new)

      # use warm start with previous solution to speed up optimization
      solver$WarmStart(x = x, y = y)
      sol_lam <- solver$Solve()
      check_osqp_status_(sol_lam, context = paste0("lambda = ", lam_seq[i]))
      x <- sol_lam$x
      y <- sol_lam$y
      solution[, i] <- pmin(pmax(lowlim, x), uplim)
    }
    solution <- as.data.frame(solution)
    names(solution) <- lam_seq
  } else {
    sol <- solver$Solve()
    check_osqp_status_(sol)
    solution <- data.frame(pmin(pmax(lowlim, sol$x), uplim))
    names(solution) <- lambda
  }
  return(solution)

}

#' Check OSQP solver status and stop with an informative error if the problem
#' was not solved.
#' @param sol Solution object returned by \code{osqp::osqp(...)$Solve()}
#' @param context Optional string appended to the error message (e.g. the
#'   lambda value being solved).
#' @keywords internal
check_osqp_status_ <- function(sol, context = NULL) {
  status <- sol$info$status
  ok_statuses <- c("solved", "solved_inaccurate")
  if (!status %in% ok_statuses) {
    ctx <- if (!is.null(context)) paste0(" [", context, "]") else ""
    stop(
      "Calibration optimization failed", ctx, ": ", status, ".\n",
      "This usually means the calibration constraints are infeasible. ",
      "Consider relaxing the constraints (e.g. collapsing covariate ",
      "categories) or checking that the ",
      "target population is compatible with the sample)."
    )
  }
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


  form <- as.formula(paste("~ . ^ ", order, " - . + 0"))

  D <- Matrix::sparse.model.matrix(form, data = cells)
  return(D)
}


#' Scale design matrix columns inversely by the number of columns at each interaction order
#'
#' For each column \eqn{j} of \code{D} whose interaction order is \eqn{K}
#' (determined by counting \code{:} in the column name), the column is
#' multiplied by \eqn{1 / \sqrt{n_K}}, where \eqn{n_K} is the total number of
#' columns of order \eqn{K}. This makes the objective equivalent to summing
#' over orders and averaging moment conditions within each order.
#'
#' @param D Sparse design matrix from \code{\link{create_design_matrix}}
#' @return A sparse matrix of the same dimensions as \code{D}, with columns
#'   scaled as described.
#'
#' @keywords internal
scale_design_matrix <- function(D) {
  if(ncol(D) == 0 || is.null(colnames(D))) return(D)
  # interaction order = number of ':' in column name (e.g. "A:B" -> order 2)
  col_orders <- nchar(colnames(D)) - nchar(gsub(":", "", colnames(D), fixed = TRUE))
  n_k <- tabulate(col_orders + 1L)[col_orders + 1L]
  D <- D %*% Matrix::Diagonal(x = 1 / sqrt(n_k))
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
                                    lowlim, uplim, exact_order = 1, verbose) {

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
  if (exact_order == 1) {
    exact_form <- formula(~ . - 1)
  } else {
    exact_form <- as.formula(paste("~ . ^", exact_order, "- 1"))
  }
  design_mat <- Matrix::sparse.model.matrix(exact_form, cells)

  A_marg <- Matrix::t(design_mat[sample_counts != 0, , drop = F] *
                      sample_counts_nempty) / sum(sample_counts)
  l_marg <- as.numeric(Matrix::t(design_mat) %*% target_counts) / sum(sample_counts)
  u_marg <- l_marg

  # Check for infeasible constraints (all-zero rows with non-zero RHS)
  all_zero_rows <- Matrix::rowSums(A_marg != 0) == 0
  infeasible_idx <- which(all_zero_rows & l_marg != 0)
  if(length(infeasible_idx) > 0) {
    row_names <- rownames(A_marg)[infeasible_idx]
    warning("Infeasible raking constraint(s). The following values never appear ",
         "in the sample but do appear in the target population:\n\n",
         paste(row_names, collapse = "\n"),
         "\n\nIt is not possible to calibrate the margin for these values. ",
         "These constraints have been dropped from the calibration procedure.\n\n",
         "Consider collapsing categories or changing the target population.")
    A_marg <- A_marg[-infeasible_idx, , drop = F]
    l_marg <- l_marg[-infeasible_idx]
    u_marg <- u_marg[-infeasible_idx]
  }

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
#' @param base_weights Vector of base weights for each unit (length = all units
#' including pop-only rows). Only the respondent entries are used.
#'
#' @keywords internal
create_rake_qvec <- function(D, sample_counts, target_counts, lambda,
                             base_weights = NULL) {

  if(is.null(base_weights)) base_weights <- rep(1, length(sample_counts))

  nnz_idxs <- which(sample_counts != 0)
  rhs <- Matrix::t(D) %*% target_counts
  lhs <- D[nnz_idxs,, drop = F] * sample_counts[nnz_idxs]

  # data-fit term: -D_nz^T * (s_nz * (D_nz * D^T * t))
  q_fit <- -c(as.numeric(lhs %*% rhs))

  # base-weight regularization term: -lambda * s_i * b_i for each respondent i
  q_base <- -lambda * sample_counts[nnz_idxs] * base_weights[nnz_idxs]

  return(q_fit + q_base)
}