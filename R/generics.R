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
