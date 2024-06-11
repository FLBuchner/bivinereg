#' brings newdata in a format appropriate for applying rvinecopulib functions.
#' @noRd
prepare_newdata <- function(newdata, object, use_response = FALSE) {
  newdata <- as.data.frame(newdata)
  if (!use_response) {
    object$model_frame <- object$model_frame[-1]
    names_ <- names(object$model_frame)
  } else {
    names_ <- object$vine$names
  }
  check_var_availability(newdata, names_)
  newdata <- remove_unused(newdata, object, use_response)
}

#' checks if all *selected* covariates are in newdata.
#' @noRd
check_var_availability <- function(newdata, vars) {
  vars_avail <- match(vars, colnames(newdata))
  if (any(is.na(vars_avail))) {
    vars_missing <- paste(vars[is.na(vars_avail)], collapse = ", ")
    stop("'newdata' is missing variables ", vars_missing)
  }
}

#' removes unused variables and returns newdata in the order used for fitting.
#' @noRd
remove_unused <- function(newdata, object, use_response = FALSE) {
  # x must be sorted in the order of the data used for fitting
  fit_order <- object$order[order(object$selected_vars)]
  if (use_response) {
    fit_order <- c(colnames(model.extract(object$model_frame, "response")), fit_order)
  }
  newdata[, fit_order, drop = FALSE]
}

#' @importFrom stats model.matrix
#' @noRd
expand_factors <- function(data) {
  if (is.data.frame(data)) {
    data <- lapply(data, function(x) {
      if (is.numeric(x) | is.ordered(x)) {
        return(x)
      }
      lvs <- levels(x)
      x <- model.matrix(~x)[, -1, drop = FALSE]
      x <- as.data.frame(x)
      x <- lapply(x, function(y) ordered(y, levels = 0:1))
      names(x) <- lvs[-1]
      x
    })
  }
  as.data.frame(data)
}
