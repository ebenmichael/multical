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
  as.numeric(Matrix::t(D) %*% (weights * sample_counts / sum(weights * sample_counts)  - target_counts / sum(target_counts)))
}

#' Get the difference between the reweighted sample and the target
#'
#' For each interaction term, computes the difference between the reweighted
#' sample count and the target count.
#'
#' @param x A \code{\link{multical}} object
#' @param order Integer. Order of interactions to include in the balance measure
#' @param ... Additional arguments (unused)
#'
#' @return A data frame with the difference for each term and value of the
#' hyper-parameter
#'
#' @export
get_balance <- function(x, order, ...) UseMethod("get_balance")

#' @method get_balance multical
#' @export
get_balance.multical <- function(x, order, ...) {

  cov_cols <- setdiff(colnames(x$cells),
                      c("sample_count", "target_count", "base_weight"))
  unit_covs <- x$cells %>% select(all_of(cov_cols))

  D <- Matrix::sparse.model.matrix(~ . - 1, unit_covs)
  if (order > 1) {
    D_int <- create_design_matrix(unit_covs, order)
    if (isTRUE(x$scale_by_order)) D_int <- scale_design_matrix(D_int)
    D <- cbind(D, D_int)
  }

  n_units <- nrow(x$cells)

  result_list <- lapply(seq_along(x$lambda), function(i) {
    # reconstruct full-length weight vector (respondents + pop-only at 0)
    w_full <- rep(0, n_units)
    w_full[x$cells$sample_count != 0] <- x$weights[, i]
    data.frame(
      lambda     = x$lambda[i],
      term       = colnames(D),
      difference = compute_balance(D, w_full, x$cells$sample_count,
                                   x$cells$target_count)
    )
  })

  do.call(rbind, result_list)
}


#' Compute the overall level of balance versus the effective sample size
#'
#' For each value of the hyper-parameter, computes the effective sample size
#' and the overall root sum of squared differences in all interactions,
#' normalized by the size of the target population.
#'
#' @param x A \code{\link{multical}} object
#' @param order Integer. Order of interactions to include in the balance measure
#' @param ... Additional arguments (unused)
#'
#' @return A data frame with the imbalance and effective sample size for each
#' value of the hyper-parameter
#'
#' @export
get_balance_v_sample_size <- function(x, order, ...) UseMethod("get_balance_v_sample_size")

#' @method get_balance_v_sample_size multical
#' @export
get_balance_v_sample_size.multical <- function(x, order, ...) {

  target_pop <- sum(x$cells$target_count)

  imbals <- get_balance(x, order) %>%
    group_by(.data$lambda) %>%
    summarise(imbalance = sqrt(sum(.data$difference ^ 2))) %>%
    ungroup()

  neff <- do.call(rbind, lapply(seq_along(x$lambda), function(i) {
    w <- x$weights[, i]
    data.frame(
      lambda = x$lambda[i],
      n_eff  = sum(w) ^ 2 / sum(w ^ 2)
    )
  }))

  inner_join(imbals, neff, by = "lambda")
}


#' Select the default lambda index for a multical object
#'
#' Picks the lambda with the largest effective sample size that achieves at
#' least \code{balance_threshold} of the total possible balance gain (measured
#' using the same \code{order} as the original \code{multical} call).
#'
#' @param weights_matrix Numeric matrix (n_respondents x n_lambdas)
#' @param cells Data frame with covariate columns plus \code{sample_count},
#'   \code{target_count}, \code{base_weight}
#' @param lambda Numeric vector of lambda values
#' @param order Integer. Order of interactions for computing balance
#' @param balance_threshold Numeric in (0, 1). Default 0.95
#'
#' @return Integer index into \code{lambda}
#' @keywords internal
select_default_lambda <- function(weights_matrix, cells, lambda, order,
                                  balance_threshold = 0.95,
                                  scale_by_order = TRUE) {
  # with a single lambda there is nothing to select
  if (length(lambda) == 1L) return(1L)

  cov_cols <- setdiff(colnames(cells),
                      c("sample_count", "target_count", "base_weight"))
  unit_covs <- cells %>% select(all_of(cov_cols))

  D <- Matrix::sparse.model.matrix(~ . - 1, unit_covs)
  if (order > 1) {
    D_int <- create_design_matrix(unit_covs, order)
    if (scale_by_order) D_int <- scale_design_matrix(D_int)
    D <- cbind(D, D_int)
  }

  target_pop <- sum(cells$target_count)
  n_units    <- nrow(cells)

  imbalances <- numeric(length(lambda))
  n_effs     <- numeric(length(lambda))

  for (i in seq_along(lambda)) {
    w_full <- rep(0, n_units)
    w_full[cells$sample_count != 0] <- weights_matrix[, i]
    imb          <- compute_balance(D, w_full, cells$sample_count,
                                    cells$target_count)
    imbalances[i] <- sqrt(sum(imb ^ 2) / target_pop)
    w            <- weights_matrix[, i]
    n_effs[i]    <- sum(w) ^ 2 / sum(w ^ 2)
  }

  # lambda[1] is the largest (least adjustment); lambda[end] is smallest
  base_imbal <- imbalances[1]
  best_imbal <- imbalances[length(lambda)]
  target_imbal <- base_imbal - balance_threshold * (base_imbal - best_imbal)

  meets <- which(imbalances <= target_imbal)
  if (length(meets) == 0L) return(length(lambda))
  meets[which.max(n_effs[meets])]
}


#' Compute the proportion of target population in uncovered interactions
#'
#' For each highest-order interaction term (i.e. the columns of the design
#' matrix at order \code{order}), aggregates \code{target_count} and
#' \code{sample_count} across all cells that contribute to that term via
#' \eqn{D_K^\top \mathbf{t}} and \eqn{D_K^\top \mathbf{s}}.  A term is
#' \emph{uncovered} if its aggregated target count is positive but its
#' aggregated sample count is zero.  The statistic is the share of the total
#' (positive) target mass that falls in uncovered terms.
#'
#' @param cells Data frame with covariate columns plus \code{sample_count},
#'   \code{target_count}, and \code{base_weight}, as stored in a
#'   \code{\link{multical}} object.
#' @param order Integer.  Order of interactions used in the calibration
#'   (i.e. \code{x$order} for a \code{multical} object \code{x}).
#'
#' @return A scalar in \eqn{[0, 1]}.
#' @keywords internal
compute_prop_uncovered <- function(cells, order) {
  cov_cols <- setdiff(colnames(cells),
                      c("sample_count", "target_count", "base_weight"))
  unit_covs <- cells[, cov_cols, drop = FALSE]

  if (order == 1) {
    # Highest-order terms are the main effects; use full one-hot encoding so
    # every factor level gets its own column (mirrors get_balance.multical).
    D_K <- Matrix::sparse.model.matrix(~ . - 1, unit_covs)
  } else {
    D_int <- create_design_matrix(unit_covs, order)
    # Column interaction order = number of ':'  (e.g. "X1a:X2b" -> 1 colon -> order 2)
    n_colons <- nchar(colnames(D_int)) -
      nchar(gsub(":", "", colnames(D_int), fixed = TRUE))
    D_K <- D_int[, n_colons == order - 1L, drop = FALSE]
  }

  if (ncol(D_K) == 0L) return(0)

  target_by_term <- as.numeric(Matrix::t(D_K) %*% cells$target_count)
  sample_by_term <- as.numeric(Matrix::t(D_K) %*% cells$sample_count)

  pos_target <- target_by_term > 0
  if (!any(pos_target)) return(0)

  total_target     <- sum(target_by_term[pos_target])
  uncovered_target <- sum(target_by_term[pos_target & sample_by_term == 0])
  uncovered_target / total_target
}


#' Print a multical object
#'
#' @param x A \code{multical} object
#' @param ... Ignored
#' @method print multical
#' @export
print.multical <- function(x, ...) {
  cat("A multical object\n")
  cat("Formula : "); print(x$formula)
  cat("Order   :", x$order, "\n")
  cat("Respondents :", x$n_respondents, "\n")
  cat(sprintf("Prop. target pop. in uncovered cells: %.1f%%\n",
              100 * x$prop_uncovered))
  cat("Lambda values:", length(x$lambda), "\n")
  if (length(x$lambda) > 1) {
    cat(sprintf(
      "Default lambda: %.4g (index %d, balance threshold = %.0f%%)\n",
      x$lambda[x$default_lambda_idx],
      x$default_lambda_idx,
      100 * x$balance_threshold
    ))
  } else {
    cat("Lambda :", x$lambda, "\n")
  }
  cat("Use weights() to extract respondent weights,",
      "get_balance() to assess calibration.\n")
  invisible(x)
}


#' Plot balance vs. effective sample size for a multical object
#'
#' Calls \code{\link{get_balance_v_sample_size}} and plots the trade-off curve,
#' highlighting the auto-selected default lambda in red.
#'
#' @param x A \code{multical} object
#' @param order Integer. Order of interactions for the balance measure;
#'   defaults to \code{x$order}
#' @param ... Passed to \code{\link[graphics]{plot}} when \pkg{ggplot2} is
#'   unavailable
#' @method plot multical
#' @export
plot.multical <- function(x, order = NULL, ...) {
  if (is.null(order)) order <- x$order
  bvn         <- get_balance_v_sample_size(x, order)
  default_row <- bvn[bvn$lambda == x$lambda[x$default_lambda_idx], , drop = FALSE]

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(bvn,
           ggplot2::aes(x = .data$n_eff, y = .data$imbalance)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::geom_point(data = default_row, colour = "red", size = 3) +
      ggplot2::labs(x = "Effective sample size",
                    y = "Normalized imbalance",
                    title = "Balance vs. effective sample size") +
      ggplot2::theme_bw()
    print(p)
  } else {
    graphics::plot(bvn$n_eff, bvn$imbalance,
         xlab = "Effective sample size", ylab = "Normalized imbalance",
         type = "b", ...)
    graphics::points(default_row$n_eff, default_row$imbalance,
                     col = "red", pch = 19, cex = 1.5)
  }
  invisible(x)
}


#' Extract calibration weights from a multical object
#'
#' Returns a numeric vector of calibration weights for respondents at the
#' selected lambda. Weights are in the same row order as the \code{sample_data}
#' passed to \code{\link{multical}}.
#'
#' @param object A \code{multical} object
#' @param lambda_idx Integer index into \code{object$lambda} selecting which
#'   lambda's weights to return. Defaults to \code{object$default_lambda_idx}.
#' @param balance_threshold Numeric in (0, 1). If supplied, chooses
#'   regularization hyper-parameter lambda by selecting the one with the
#'   largest effective sample size that achieves at least
#'   \code{balance_threshold} of the total possible balance gain (measured
#'   using the same \code{order} as the original \code{multical} call). If
#'   supplied, overrides any explicit value of \code{lambda_idx}.
#'   selection with this threshold and ignores \code{lambda_idx}.
#' @param ... Ignored
#' @method weights multical
#' @export
weights.multical <- function(object, lambda_idx = NULL,
                             balance_threshold = NULL, ...) {
  if (!is.null(balance_threshold)) {
    lambda_idx <- select_default_lambda(
      object$weights, object$cells, object$lambda, object$order,
      balance_threshold
    )
  }
  if (is.null(lambda_idx)) lambda_idx <- object$default_lambda_idx
  object$weights[, lambda_idx]
}


#' Estimate the population mean of an outcome from a multical object
#'
#' Uses either the weighted mean or a model-assisted estimator
#' and returns a point estimate and standard error.
#'
#' @param object A \code{\link{multical}} object
#' @param ... Arguments passed to the method (see \code{\link{estimate.multical}})
#'
#' @export
estimate <- function(object, ...) UseMethod("estimate")

#' @describeIn estimate Method for \code{multical} objects
#'
#' @param y Outcome variable. Can be:
#'   \itemize{
#'     \item A bare column name evaluated in \code{data}
#'       (e.g. \code{estimate(cal, y, data = sample_df)})
#'     \item Any expression that evaluates to a numeric vector of length
#'       \code{object$n_respondents} (e.g. \code{estimate(cal, sample_df$y)})
#'     \item A \strong{factor} or \strong{character} vector of length
#'       \code{object$n_respondents}. Each level is treated as a binary
#'       indicator and estimated separately. Character vectors are coerced
#'       to factor (levels in alphabetical order) before processing.
#'   }
#' @param data Optional data frame in which to evaluate \code{y}.
#' @param method Character string selecting the estimator:
#'   \describe{
#'     \item{\code{"linearized"} (default)}{Calibration-weighted mean with a
#'       regression-assisted SE. OLS residuals from a first-order design matrix
#'       are used in the sandwich SE formula.}
#'     \item{\code{"hajek"}}{Plain calibration-weighted mean (Horvitz-Thompson
#'       style) with Taylor linearization SE.}
#'     \item{\code{"greg"}}{GREG (generalized regression) estimator with
#'       cross-fitted ridge regression. The point estimate is the weighted mean
#'       of ridge-predicted values in the target population plus the
#'       calibration-weighted mean of residuals. The interaction order matches
#'       the order used in the \code{multical} call.}
#'     \item{\code{"drp"}}{Double Regression with Post-Stratification estimator. Uses
#'       cross-fitted gradient-boosted trees (xgboost) in place of OLS.
#'       Requires the \pkg{xgboost} package.}
#'   }
#' @param lambda_idx Integer index into \code{object$lambda} selecting which
#'   lambda's weights to use. Defaults to \code{object$default_lambda_idx}.
#' @param balance_threshold Numeric in (0, 1). If supplied, re-runs lambda
#'   selection with this threshold (overrides \code{lambda_idx}).
#' @param nfolds Integer. Number of cross-fitting folds used when
#'   \code{method = "greg"} or \code{method = "drp"}. Default \code{3}.
#' @param subset Optional vector for restricting estimation to a subset of
#'   respondents. Can be:
#'   \itemize{
#'     \item \code{NULL} (default) --- all respondents are used.
#'     \item A \strong{logical} vector of length \code{object$n_respondents}.
#'       \code{TRUE} retains the corresponding respondent.
#'     \item An \strong{integer} index vector. Each element is a 1-based
#'       position of a respondent to retain.
#'   }
#'   Any respondent whose outcome \code{y} is \code{NA} and who would otherwise
#'   be retained by \code{subset} is automatically excluded, and a warning is
#'   issued stating how many observations were dropped.
#' @param ... Additional arguments to pass to xgboost when \code{method = "drp"}. Ignored for other methods.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{level}{(Only present when \code{y} is a factor or character.) The
#'     factor level. Rows are in \code{levels(y)} order.}
#'   \item{estimate}{Point estimate of the population mean (or proportion for
#'     a factor/character outcome).}
#'   \item{se}{Standard error.}
#'   \item{lambda}{The lambda value used.}
#'   \item{method}{The estimator name.}
#' }
#' For numeric \code{y} the data frame has one row. For factor or character
#' \code{y} it has one row per level.
#' @method estimate multical
#' @export
estimate.multical <- function(object, y, data = NULL, method = "linearized",
                              lambda_idx = NULL, balance_threshold = NULL,
                              nfolds = 3, subset = NULL, ...) {
  y_quo <- enquo(y)
  y_vec <- eval_tidy(y_quo, data = data)

  # coerce character to factor
  if (is.character(y_vec)) y_vec <- factor(y_vec)

  # single length check (applies before both factor and numeric paths)
  if (length(y_vec) != object$n_respondents) {
    stop("`y` must have one value per respondent (",
         object$n_respondents, " expected, ",
         length(y_vec), " supplied).")
  }

  # --- build keep mask from subset ---
  n_resp <- object$n_respondents
  if (is.null(subset)) {
    keep <- rep(TRUE, n_resp)
  } else if (is.logical(subset)) {
    if (length(subset) != n_resp) {
      stop("`subset` must have length equal to the number of respondents (",
           n_resp, " expected, ", length(subset), " supplied).")
    }
    keep <- subset
  } else if (is.numeric(subset) || is.integer(subset)) {
    idx <- as.integer(subset)
    if (any(idx < 1L | idx > n_resp)) {
      stop("`subset` contains out-of-range indices (must be between 1 and ",
           n_resp, ").")
    }
    keep <- seq_len(n_resp) %in% idx
  } else {
    stop("`subset` must be NULL, a logical vector, or an integer index vector.")
  }

  # --- NA detection: warn once (before factor loop), extend keep ---
  na_in_keep <- is.na(y_vec) & keep
  if (any(na_in_keep)) {
    n_na <- sum(na_in_keep)
    warning(n_na, " observation(s) with missing outcome (NA) excluded from estimation.")
    keep[na_in_keep] <- FALSE
  }

  # fan out over factor levels, passing the already-finalized keep mask
  if (is.factor(y_vec)) {
    lvls <- levels(y_vec)
    results <- lapply(lvls, function(lvl) {
      y_bin <- as.numeric(y_vec == lvl)
      out   <- estimate.multical(object, y_bin, data = NULL,
                                 method = method,
                                 lambda_idx = lambda_idx,
                                 balance_threshold = balance_threshold,
                                 nfolds = nfolds,
                                 subset = keep, ...)
      cbind(level = lvl, out, stringsAsFactors = FALSE)
    })
    return(do.call(rbind, results))
  }

  if (!is.numeric(y_vec)) stop("`y` must evaluate to a numeric, factor, or character vector.")

  # apply subset mask to outcome and weights
  y_vec <- y_vec[keep]
  w_vec <- weights(object, lambda_idx = lambda_idx,
                   balance_threshold = balance_threshold)[keep]

  lam_idx   <- if (!is.null(balance_threshold)) {
    select_default_lambda(object$weights, object$cells, object$lambda,
                          object$order, balance_threshold)
  } else if (!is.null(lambda_idx)) {
    lambda_idx
  } else {
    object$default_lambda_idx
  }
  lam_val <- object$lambda[lam_idx]

  if (method == "hajek") {
    out <- hajek_estimate_(y_vec, w_vec)
  } else if (method == "linearized") {
    cov_cols   <- setdiff(colnames(object$cells),
                          c("sample_count", "target_count", "base_weight"))
    cells_resp <- object$cells[object$cells$sample_count != 0, cov_cols,
                               drop = FALSE][keep, , drop = FALSE]
    out <- linearized_estimate_(y_vec, w_vec, cells_resp, order = 1,
                                use_ridge = FALSE)
  } else if (method == "greg") {
    cov_cols   <- setdiff(colnames(object$cells),
                          c("sample_count", "target_count", "base_weight"))
    cells_resp <- object$cells[object$cells$sample_count != 0, cov_cols,
                               drop = FALSE][keep, , drop = FALSE]
    cells_target <- object$cells[object$cells$target_count != 0, cov_cols,
                               drop = FALSE]
    target_counts <- object$cells$target_count[object$cells$target_count != 0]
    out <- greg_estimate_(y_vec, w_vec, cells_resp, cells_target, target_counts,
                          object$order, use_ridge = TRUE, nfolds = nfolds)
  } else if(method == "drp") {
    cov_cols   <- setdiff(colnames(object$cells),
                          c("sample_count", "target_count", "base_weight"))
    cells_resp <- object$cells[object$cells$sample_count != 0, cov_cols,
                               drop = FALSE][keep, , drop = FALSE]
    cells_target <- object$cells[object$cells$target_count != 0, cov_cols,
                                 drop = FALSE]
    target_counts <- object$cells$target_count[object$cells$target_count != 0]
    out <- drp_estimate_(y_vec, w_vec, cells_resp, cells_target, target_counts,
                         nfolds = nfolds, ...)
  } else {
    stop('Unknown method "', method, '". Supported: "hajek", "linearized", "greg", "drp".')
  }
  out$lambda <- lam_val
  return(out)
}


#' Internal Hajek (weighted mean) estimator
#'
#' @param y_vec Numeric outcome vector (length = n_respondents)
#' @param w_vec Numeric weight vector (length = n_respondents)
#'
#' @return A one-row data frame with columns estimate, se, method
#' @keywords internal
hajek_estimate_ <- function(y_vec, w_vec) {
  w_sum  <- sum(w_vec)
  mu_hat <- sum(w_vec * y_vec) / w_sum
  # Taylor linearization variance: sum(w_i^2 (y_i - mu_hat)^2) / (sum w_i)^2
  se     <- sqrt(sum(w_vec ^ 2 * (y_vec - mu_hat) ^ 2)) / w_sum
  data.frame(estimate = mu_hat, se = se, method = "hajek")
}


#' Internal estimator with linearized standard error
#'
#' Point estimate is the calibration-weighted mean of \code{y}. The SE uses
#' residuals from a regression of \code{y} on the covariate design matrix built
#' from \code{cells_resp} at interaction order \code{order}. When
#' \code{use_ridge = TRUE}, ridge regression is fitted via \code{cv.glmnet}
#' (\code{alpha = 0}, penalty at \code{lambda.min}); otherwise plain OLS via
#' \code{lm.fit} is used.
#'
#' @param y_vec Numeric outcome vector (length = n_respondents)
#' @param w_vec Numeric weight vector (length = n_respondents)
#' @param cells_resp Data frame of covariate columns for respondents only
#' @param order Integer. Order of interactions for the regression model
#' @param use_ridge Logical. Use \code{cv.glmnet} ridge regression to fit
#'   the outcome model. Default \code{FALSE}.
#'
#' @return A one-row data frame with columns estimate, se, method
#' @keywords internal
linearized_estimate_ <- function(y_vec, w_vec, cells_resp, order,
                                 use_ridge = FALSE) {
  form_str <- if (order == 1) "~ ." else paste0("~ .^", order)
  X        <- model.matrix(as.formula(form_str), data = cells_resp)
  if (use_ridge & ncol(X) > 2) {
    X_noInt <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    cv_fit <- cv.glmnet(X_noInt, y_vec, alpha = 0)
    fitted <- as.numeric(predict(cv_fit, newx = X_noInt, s = "lambda.min"))
    resids <- y_vec - fitted
  } else {
    resids <- lm.fit(X, y_vec)$residuals
  }

  w_sum  <- sum(w_vec)
  mu_hat <- sum(w_vec * y_vec) / w_sum
  se     <- sqrt(sum(w_vec ^ 2 * resids ^ 2)) / w_sum
  data.frame(estimate = mu_hat, se = se, method = "linearized")
}



#' Internal GREG (generalized regression) estimator
#'
#' The point estimate is the population-weighted mean of regression predictions
#' on the target cells plus the calibration-weighted mean of residuals:
#' \eqn{\hat{\mu} = \bar{y}_{\text{target}} + \sum_i w_i e_i / \sum_i w_i},
#' where \eqn{e_i = y_i - \hat{y}_i}. The SE uses a sandwich estimator based
#' on the same residuals. The outcome model is OLS (via \code{lm.fit}) by
#' default; linearly-dependent columns (NA coefficients) are dropped
#' symmetrically from both the respondent and target design matrices. When
#' \code{use_ridge = TRUE}, cross-fitted ridge regression is used instead
#' (via \code{\link{crossfit_ridge}}, which calls \code{cv.glmnet} with
#' \code{alpha = 0} on each training fold).
#'
#' @param y_vec Numeric outcome vector (length = n_respondents)
#' @param w_vec Numeric weight vector (length = n_respondents)
#' @param cells_resp Data frame of covariate columns for respondents only
#' @param cells_target Data frame of covariate columns for target population only
#' @param target_counts Vector of target counts for each cell in \code{cells_target}
#' @param order Integer. Order of interactions for the regression model
#' @param use_ridge Logical. Use cross-fitted ridge regression (via
#'   \code{cv.glmnet} with \code{alpha = 0}) to fit the outcome model.
#'   Default \code{TRUE}.
#' @param nfolds Integer. Number of cross-fitting folds when
#'   \code{use_ridge = TRUE}. Default \code{3}. Ignored when
#'   \code{use_ridge = FALSE}.
#'
#' @return A one-row data frame with columns estimate, se, method
#' @keywords internal
greg_estimate_ <- function(y_vec, w_vec, cells_resp, cells_target,
                                 target_counts, order,
                                 use_ridge = TRUE, nfolds = 3) {
  form_str <- if (order == 1) "~ ." else paste0("~ .^", order)
  X        <- model.matrix(as.formula(form_str), data = cells_resp)
  X_noInt  <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  X_target <- model.matrix(as.formula(form_str), data = cells_target)
  X_target_noInt <- X_target[, colnames(X_target) != "(Intercept)", drop = FALSE]
  if (use_ridge) {
    cf_out      <- crossfit_ridge(X_noInt, y_vec, X_target_noInt, nfolds = nfolds)
    resids      <- y_vec - cf_out$pred
    target_mean <- sum(cf_out$target_pred * target_counts) / sum(target_counts)
  } else {
    coefs         <- lm.fit(X, y_vec)$coefficients
    keep          <- !is.na(coefs)
    fitted        <- as.numeric(X[, keep, drop = FALSE]        %*% coefs[keep])
    target_pred   <- as.numeric(X_target[, keep, drop = FALSE] %*% coefs[keep])
    resids        <- y_vec - fitted
    target_mean   <- sum(target_pred * target_counts) / sum(target_counts)
  }

  w_sum  <- sum(w_vec)
  mu_hat <- target_mean + sum(w_vec * resids) / w_sum
  se     <- sqrt(sum(w_vec ^ 2 * resids ^ 2)) / w_sum
  data.frame(estimate = mu_hat, se = se, method = "greg")
}




#' Estimate outcome model with cross-fitting using ridge regression
#'
#' Splits the respondent data into \code{nfolds} folds. For each fold, fits a
#' ridge regression model (via \code{cv.glmnet} with \code{alpha = 0}) on the
#' remaining folds and predicts on (a) the held-out fold and (b) all
#' target-population cells. The target-population predictions are averaged
#' across the \code{nfolds} models.
#'
#' @param X Numeric matrix of respondent covariates (no intercept)
#' @param y Numeric vector of observed outcomes (length = nrow(X))
#' @param target_X Numeric matrix of target-population covariates (no intercept)
#' @param nfolds Integer. Number of cross-fitting folds (default: 3)
#' @param ... Additional arguments forwarded to \code{glmnet::cv.glmnet()}
#'
#' @return A named list with:
#'   \describe{
#'     \item{pred}{Numeric vector of cross-fitted predictions for respondents}
#'     \item{target_pred}{Numeric vector of averaged predictions for target cells}
#'   }
#'
#' @keywords internal
crossfit_ridge <- function(X, y, target_X, nfolds = 3, ...) {

  n <- nrow(X)
  folds <- split(sample(n), cut(1:n, nfolds))

  pred <- numeric(n)
  target_pred <- matrix(0, nrow = nrow(target_X), ncol = length(folds))
  k <- 1
  for (fold in folds) {
    tr_idx <- setdiff(1:n, fold)
    cv_fit <- cv.glmnet(X[tr_idx, , drop = FALSE], y[tr_idx], alpha = 0, ...)
    pred[fold] <- as.numeric(predict(cv_fit, newx = X[fold, , drop = FALSE],
                                     s = "lambda.min"))
    target_pred[, k] <- as.numeric(predict(cv_fit, newx = target_X,
                                           s = "lambda.min"))
    k <- k + 1
  }
  return(list(pred = pred, target_pred = rowMeans(target_pred)))
}



#' Estimate outcome model with cross-fitting using xgboost
#'
#' Splits the respondent data into \code{nfolds} folds. For each fold, fits an
#' xgboost model on the remaining folds and predicts on (a) the held-out fold
#' and (b) all target-population cells. The target-population predictions are
#' averaged across the \code{nfolds} models.
#'
#' @param X Numeric matrix of respondent covariates (no intercept)
#' @param y Numeric vector of observed outcomes (length = nrow(X))
#' @param target_X Numeric matrix of target-population covariates (no intercept)
#' @param nfolds Integer. Number of cross-fitting folds (default: 3)
#' @param nrounds Number of boosting rounds for xgboost (default: 1000)
#' @param verbose Verbosity level passed to xgboost (default: 0, silent)
#' @param early_stopping_rounds Number of consecutive non-improving rounds
#'   before stopping (default: 20)
#' @param ... Additional arguments forwarded to \code{xgboost::xgboost()}
#'
#' @return A named list with:
#'   \describe{
#'     \item{pred}{Numeric vector of cross-fitted predictions for respondents}
#'     \item{target_pred}{Numeric vector of averaged predictions for target cells}
#'   }
#'
#' @keywords internal
crossfit_xgboost <- function(X, y, target_X, nfolds = 3, nrounds = 1000,
                             verbose = 0, early_stopping_rounds = 20, ...) {

  n <- nrow(X)
  folds <- split(sample(n), cut(1:n, nfolds))

  pred <- numeric(n)
  target_pred <- matrix(0, nrow = nrow(target_X), ncol = length(folds))
  k <- 1
  for (fold in folds) {

    tr_idx <- setdiff(1:n, fold)
    Xtr <- X[tr_idx, , drop = FALSE]
    ytr <- y[tr_idx]
    gb <- xgboost(data = Xtr, label = ytr, objective = "reg:squarederror",
                  nrounds = nrounds, verbose = verbose, early_stopping_rounds = early_stopping_rounds, ...)
    pred[fold] <- predict(gb, X[fold, , drop = FALSE])
    target_pred[, k] <- predict(gb, target_X)
    k <- k + 1
  }
  return(list(pred = pred, target_pred = rowMeans(target_pred)))
}



#' Internal double regression and post-stratification (DRP) estimator
#'
#' Uses cross-fitted xgboost predictions as the outcome model. The
#' point estimate is the target-population mean of the cross-fitted predictions
#' plus the calibration-weighted mean of residuals (same GREG-style formula as
#' \code{\link{greg_estimate_}}). The SE uses a sandwich estimator based on
#' the cross-fitted residuals. Additional arguments in \code{...} are forwarded
#' to \code{\link{crossfit_xgboost}} and then to \code{xgboost::xgboost()}.
#'
#' @param y_vec Numeric outcome vector (length = n_respondents)
#' @param w_vec Numeric weight vector (length = n_respondents)
#' @param cells_resp Data frame of covariate columns for respondents only
#' @param cells_target Data frame of covariate columns for target population only
#' @param target_counts Vector of target counts for each cell in \code{cells_target}
#' @param nfolds Integer. Number of cross-fitting folds (default: 3)
#' @param ... Additional arguments to pass to xgboost
#'
#' @return A one-row data frame with columns estimate, se, method
#' @keywords internal
drp_estimate_ <- function(y_vec, w_vec, cells_resp, cells_target,
                          target_counts, nfolds = 3, ...) {
  X        <- model.matrix(~ . - 1, data = cells_resp)
  X_target <- model.matrix(~ . - 1, data = cells_target)

  xboost_out <- crossfit_xgboost(X, y_vec, X_target, nfolds = nfolds, ...)
  target_mean <- sum(xboost_out$target_pred * target_counts) / sum(target_counts)
  resids <- y_vec - xboost_out$pred
  w_sum  <- sum(w_vec)
  mu_hat <- target_mean + sum(w_vec * resids) / w_sum
  se     <- sqrt(sum(w_vec ^ 2 * resids ^ 2)) / w_sum
  data.frame(estimate = mu_hat, se = se, method = "drp")
}
