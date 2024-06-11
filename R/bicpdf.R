#' Bivariate conditional PDF
#'
#' Calculates the bivariate conditional density of the responses given the
#' covariates.
#'
#' @param object an object of class \code{bivinereg}.
#' @param newdata matrix of responses and covariate values for which to compute
#'   the bivariate conditional density.
#' @param cores integer; the number of cores to use for computations.
#'
#' @examples
#' # load sample data
#' data(data)
#'
#' # fit vine regression model
#' (fit <- bivinereg(cbind(U1,U4) ~ U2 + U3 + U5 + U6,
#'                   data = data,
#'                   family_set = "parametric",
#'                   selcrit = "bic"))
#'
#' # calculate bivariate conditional PDF
#' bicpdf(fit, data[1:5,])
#'
#' @export
bicpdf <- function(object, newdata, cores = 1) {
  newdata <- prepare_newdata(newdata, object, use_response = TRUE)
  cond_bi_dens_cpp(as.matrix(newdata), object$vine, cores)
}

#' Marginal conditional PDF
#'
#' Calculates the marginal conditional density of a response given the
#' covariates.
#'
#' @param object an object of class \code{bivinereg}.
#' @param newdata matrix of response and covariate values for which to compute
#'   the bivariate conditional density.
#' @param margin integer; the margin for which to calculate the conditional PDF.
#' @param cores integer; the number of cores to use for computations.
#'
#' @examples
#' # load sample data
#' data(data)
#'
#' # fit vine regression model
#' (fit <- bivinereg(cbind(U1,U4) ~ U2 + U3 + U5 + U6,
#'                   data = data,
#'                   family_set = "parametric",
#'                   selcrit = "bic"))
#'
#' # calculate marginal conditional PDF of U1 given the covariates
#' mcpdf(fit, data[1:5,], margin = 1)
#'
#' @export
mcpdf <- function(object, newdata, margin, cores = 1) {
  newdata <- prepare_newdata(newdata, object, use_response = TRUE)
  cond_m_dens_cpp(as.matrix(newdata), object$vine, margin - 1, cores)
}

#' Marginal conditional CDF
#'
#' Calculates the marginal conditional distribution of a response given the
#' covariates and other response.
#'
#' @param object an object of class \code{bivinereg}.
#' @param newdata matrix of response and covariate values for which to compute
#'   the marginal conditional distribution.
#' @param margin integer; the margin for which to calculate the conditional CDF.
#' @param cores integer; the number of cores to use for computations.
#'
#' @examples
#' # load sample data
#' data(data)
#'
#' # fit vine regression model
#' (fit <- bivinereg(cbind(U1,U4) ~ U2 + U3 + U5 + U6,
#'                   data = data,
#'                   family_set = "parametric",
#'                   selcrit = "bic"))
#'
#' # calculate the marginal conditional CDF of U1 given U4 and the covariates
#' mccdf(fit, data[1:5,], margin = 1)
#'
#' @export
mccdf <- function(object, newdata, margin, cores = 1) {
  newdata <- prepare_newdata(newdata, object, use_response = TRUE)
  cond_m_dist_cpp(as.matrix(newdata), object$vine, margin, cores)
}

#' Bivariate conditional CDF
#'
#' Calculates the bivariate conditional distribution of the responses given the
#' covariates.
#'
#' @param object an object of class \code{bivinereg}.
#' @param newdata matrix of response and covariate values for which to compute
#'   the bivariate conditional distribution.
#' @param cores integer; the number of cores to use for computations.
#'
#' @importFrom stats integrate
#'
#' @examples
#' # load sample data
#' data(data)
#'
#' # fit vine regression model
#' (fit <- bivinereg(cbind(U1,U4) ~ U2 + U3 + U5 + U6,
#'                   data = data,
#'                   family_set = "parametric",
#'                   selcrit = "bic"))
#'
#' # calculate the bivariate conditional CDF
#' biccdf(fit, data[1:5,])
#'
#' @export
biccdf <- function(object, newdata, cores = 1) {
  newdata <- prepare_newdata(newdata, object, use_response = TRUE)
  out <- NA
  for (i in 1:nrow(newdata)) {
    out[i] <- integrate(biccdf_integrant, 0, newdata[i, 1],
                        newdata = newdata[i,], object = object,
                        cores = cores)$value
  }
  return(out)
}

#' @noRd
biccdf_integrant <- function(v, newdata, object, cores) {
  newdata <- as.matrix(cbind.data.frame(v, newdata[1, -1], row.names = NULL))
  return(cond_m_dens_cpp(newdata, object$vine, 0, cores) *
           cond_m_dist_cpp(newdata, object$vine, 2, cores))
}

#' Bivariate conditional probability
#'
#' Calculates the bivariate conditional probability of the responses given the
#' covariates.
#'
#' @param object an object of class \code{bivinereg}.
#' @param newdata matrix of covariate values for which to compute
#'   the bivariate conditional probability.
#' @param v1_min lower limit of first response.
#' @param v1_max upper limit of first response.
#' @param v2_min lower limit of second response.
#' @param v2_max upper limit of second response.
#' @param cores integer; the number of cores to use for computations.
#'
#' @details
#' This function double integrates the bivariate conditional density using the
#' function \code{integral2} from the \code{pracma} package.
#'
#' @seealso \code{\link{bicprob2}}
#' @importFrom pracma integral2
#'
#' @examples
#' # load sample data
#' data(data)
#'
#' # fit vine regression model
#' (fit <- bivinereg(cbind(U1,U4) ~ U2 + U3 + U5 + U6,
#'                   data = data,
#'                   family_set = "parametric",
#'                   selcrit = "bic"))
#'
#' # calculate the conditional probability that 0.3 <= U1 <= 0.6,
#' # 0.3 <= U4 <= 0.6, given the covariates
#' bicprob(fit, data[1:5,], rep(0.3, 5), rep(0.6, 5), rep(0.3, 5), rep(0.6, 5))
#'
#' @export
bicprob <- function(object, newdata, v1_min, v1_max, v2_min, v2_max, cores = 1)
{
  newdata <- prepare_newdata(newdata, object)
  out <- NA
  for (i in 1:nrow(newdata)) {
    out[i] <- integral2(bicprob_integrant,
                        v1_min[i], v1_max[i], v2_min[i], v2_max[i],
                        newdata = newdata[i,], object = object, cores = cores)$Q
  }
  return(out)
}

#' @noRd
bicprob_integrant <- function(v1, v2, newdata, object, cores) {
  out <- matrix(NA, dim(v1)[1], dim(v1)[2])
  for (i in 1:dim(v1)[2]) {
    newdata_ <- as.matrix(cbind.data.frame(v1[, i], v2[, i],
                                           newdata[rep(1, dim(v1)[1]),],
                                           row.names = NULL))
    out[, i] <- cond_bi_dens_cpp(newdata_, object$vine, cores)
  }
  return(out)
}

#' Bivariate conditional probability
#'
#' Calculates the bivariate conditional probability of the responses given the
#' covariates.
#'
#' @param object an object of class \code{bivinereg}.
#' @param newdata matrix of covariate values for which to compute
#'   the bivariate conditional probability.
#' @param v1_min lower limit of first response.
#' @param v1_max upper limit of first response.
#' @param v2_min lower limit of second response.
#' @param v2_max upper limit of second response.
#' @param cores integer; the number of cores to use for computations.
#'
#' @details
#' This function double integrates the bivariate conditional density using the
#' function \code{hcubature} from the \code{cubature} package.
#'
#' @seealso \code{\link{bicprob}}
#'
#'
#' @importFrom cubature hcubature
#'
#' @examples
#' # load sample data
#' data(data)
#'
#' # fit vine regression model
#' (fit <- bivinereg(cbind(U1,U4) ~ U2 + U3 + U5 + U6,
#'                   data = data,
#'                   family_set = "parametric",
#'                   selcrit = "bic"))
#'
#' # calculate the conditional probability that 0.3 <= U1 <= 0.6,
#' # 0.3 <= U4 <= 0.6, given the covariates
#' bicprob2(fit, data[1:5,], rep(0.3, 5), rep(0.6, 5), rep(0.3, 5), rep(0.6, 5))
#'
#' @export
bicprob2 <- function(object, newdata, v1_min, v1_max, v2_min, v2_max, cores = 1)
{
  newdata <- prepare_newdata(newdata, object)
  out <- NA
  for (i in 1:nrow(newdata)) {
    out[i] <- hcubature(bicprob_integrant2,
                        c(v1_min[i], v2_min[i]), c(v1_max[i], v2_max[i]),
                        newdata = newdata[i,], object = object,
                        cores = cores)$integral
  }
  return(out)
}

#' @noRd
bicprob_integrant2 <- function(v, newdata, object, cores) {
  newdata_ <- as.matrix(cbind.data.frame(v[1], v[2], newdata[1,],
                                         row.names = NULL))
  return(cond_bi_dens_cpp(newdata_, object$vine, cores))
}
