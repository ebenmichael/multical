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
    D <- cbind(D, create_design_matrix(unit_covs, order))
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
    summarise(imbalance = sqrt(sum(.data$difference ^ 2) / target_pop)) %>%
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
#' @param balance_threshold Numeric in (0, 1). Default 0.9
#'
#' @return Integer index into \code{lambda}
#' @keywords internal
select_default_lambda <- function(weights_matrix, cells, lambda, order,
                                  balance_threshold = 0.9) {
  # with a single lambda there is nothing to select
  if (length(lambda) == 1L) return(1L)

  cov_cols <- setdiff(colnames(cells),
                      c("sample_count", "target_count", "base_weight"))
  unit_covs <- cells %>% select(all_of(cov_cols))

  D <- Matrix::sparse.model.matrix(~ . - 1, unit_covs)
  if (order > 1) {
    D <- cbind(D, create_design_matrix(unit_covs, order))
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