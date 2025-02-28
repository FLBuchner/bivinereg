#' Sample conditionally from a Y-vine regression model
#'
#' @param object an object of class \code{bivinereg}.
#' @param newdata data.frame of covariate values for which to sample.
#' @param n integer; the number of samples to generate.
#' @param cores integer; the number of cores to use for computations.
#'
#' @return A data.frame of samples with two columns corresponding to the
#' responses.
#'
#' @examples
#' # load sample data
#' data(data)
#'
#' # fit vine regression model
#' (fit <- bivinereg(cbind(U1, U4) ~ U2 + U3 + U5 + U6,
#'                   data = data,
#'                   family_set = "parametric",
#'                   selcrit = "bic",
#'                   uscale = TRUE))
#'
#' # Sample conditionally
#' bicsample(fit, 10, data[1,])
#'
#' @export
#'
#' @importFrom stats runif
#'
bicsample <- function(object, n, newdata, cores = 1) {
  newdata <- prepare_newdata(newdata, object)
  newdata <- to_uscale(newdata, object$margins[-(1:2)], add_response = TRUE)
  v <- matrix(runif(2 * n, 0, 1), nrow = n, ncol = 2)
  preds <- as.data.frame(cond_sample_cpp(v, as.matrix(newdata), object$vine, cores))
  names(preds) <- colnames(model.extract(object$model_frame, "response"))
  preds <- to_yscale(preds, object)

  preds
}
