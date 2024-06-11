#' @export
print.bivinereg <- function(x, ...) {
  cat("Y-vine regression model: ")
  n_predictors <- length(x$order)
  if (n_predictors <= 10) {
    predictors <- paste(x$order, collapse = ", ")
  } else {
    predictors <- paste(x$order[1:10], collapse = ", ")
    predictors <- paste0(predictors, ", ... (", n_predictors - 10, " more)")
  }
  names <- colnames(model.extract(x$model_frame, "response"))
  cat(paste(names[1], names[2], sep = ", "), "|", predictors, "\n")
  stats <- unlist(x$stats[1:5])
  stats <- paste(names(stats), round(stats, 2), sep = " = ")
  cat(paste(stats, collapse = ", "), "\n")
  invisible(x)
}

#' @export
summary.bivinereg <- function(object, ...) {
  names <- colnames(model.extract(object$model_frame, "response"))
  data.frame(
    var = c(paste(names[1], names[2], sep = ", "), object$order),
    edf = object$stats$var_edf,
    cll = object$stats$var_cll,
    caic = object$stats$var_caic,
    cbic = object$stats$var_cbic,
    p_value = object$stats$var_p_value
  )
}

#' Summary for Y-vines
#'
#' A summary for Y-vines tailored to Y-vine regression.
#'
#' @param object An object of class "vinecop". The fitted vine from `bivinereg()`.
#' @param digits Digits in summary.
#' @param trees The trees for which to print the summary. Either a vector of
#'    numbers or "ALL" (default) to print all trees.
#' @param rows Specific rows of the summary to be printed. Either a vector of
#'    numbers or "ALL" (default) for all rows.
#' @param names_cols A vector of names of columns to be printed in the summary.
#' @param ... Additional parameters for \code{print.data.frame}.
#'
#' @export
summary_yvine <- function(object, digits = 2, trees = "ALL", rows = "ALL",
                         names_cols = c("tree", "edge", "cond'd",
                                        "conditioning", "family", "rota",
                                        "parameters", "df", "tau", "loglik"),
                         ...) {
  summary_df_ <- summary(object)
  # Change name of columns conditioned to cond'd and rotation to rota,
  # then the summary should fit into a page without line breaks
  names(summary_df_)[c(3, 7)] <- c("cond'd", "rota")

  names <- object$names
  d <- object$structure$d - 1

  resp_alloc <- paste0("Responses:  1 = ", names[1], ", 2 = ", names[2], "\n")
  cov_alloc <- "Covariates: "
  for (i in 3:length(names)) {
    cov_alloc <- paste0(cov_alloc, i, " = ", names[i])
    if (i < length(names)) {
      cov_alloc <- paste0(cov_alloc, ", ")
    } else {
      cov_alloc <- paste0(cov_alloc, "\n")
    }
  }

  if (any(rows == "ALL")) {
    if (any(trees == "ALL")) {
      rows <- 1:dim(summary_df_)[1]
    } else if (all(trees %in% 1:d)) {
      rows <- NULL
      k <- 0
      for (tree in 1:d) {
        if (tree %in% trees) {
          rows <- c(rows, (k+1):(k+d-tree+1))
        }
        k <- k + d - tree + 1
      }
    } else {
      stop("Invalid tree(s) selected.")
    }
  } else {
    if (!all(rows %in% 1:dim(summary_df_)[1])) {
      stop("Invalid row(s) selected.")
    }
  }

  cat("A data.frame:", dim(summary_df_)[1], "x", dim(summary_df_)[2], "\n")
  cat(resp_alloc)
  cat(cov_alloc)
  print.data.frame(summary_df_[rows, names_cols], digits = digits,
                   row.names = FALSE, ...)
}
