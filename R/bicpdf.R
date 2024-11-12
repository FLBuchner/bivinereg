#' Bivariate conditional PDF
#'
#' Calculates the bivariate conditional density \eqn{c_{V_1 V_2 | \mathbf{U}}}
#' of the responses \eqn{(V_1, V_2)} given covariates \eqn{\mathbf{U}}.
#'
#' @param object an object of class \code{bivinereg}.
#' @param newdata data.frame of new response and covariate values for which to
#'   compute the bivariate conditional density.
#' @param cores integer; the number of cores to use for computations.
#'
#' @examples
#' # load sample data
#' data(data)
#'
#' # fit vine regression model
#' (fit <- bivinereg(cbind(U1, U4) ~ U2 + U3 + U5 + U6,
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
  dens_marg <- if (inherits(object$margins[[1]], "kde1d")) {
    kde1d::dkde1d(newdata[, 1], object$margins[[1]]) *
      kde1d::dkde1d(newdata[, 2], object$margins[[2]])
  } else {
    1
  }
  newdata <- to_uscale(newdata, object$margins)
  cond_bi_dens_cpp(as.matrix(newdata), object$vine, cores) * dens_marg
}

#' Marginal conditional PDF
#'
#' Calculates the marginal conditional density \eqn{c_{V_i | \mathbf{U}}} of a
#' response \eqn{V_i, i = 1, 2,} given covariates \eqn{\mathbf{U}}.
#'
#' @param object an object of class \code{bivinereg}.
#' @param newdata data.frame of new response and covariate values for which to
#'   compute the marginal conditional density.
#' @param margin integer; the margin (1 or 2) for which to calculate the
#'   conditional PDF.
#' @param cores integer; the number of cores to use for computations.
#'
#' @examples
#' # load sample data
#' data(data)
#'
#' # fit vine regression model
#' (fit <- bivinereg(cbind(U1, U4) ~ U2 + U3 + U5 + U6,
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
  dens_marg <- if (inherits(object$margins[[margin]], "kde1d")) {
    kde1d::dkde1d(newdata[, margin], object$margins[[margin]])
  } else {
    1
  }
  newdata <- to_uscale(newdata, object$margins)
  cond_m_dens_cpp(as.matrix(newdata), object$vine, margin - 1, cores) * dens_marg
}

#' Marginal conditional CDF
#'
#' Either calculates the marginal conditional distribution
#' \eqn{C_{V_i | V_j, \mathbf{U}}} of a response \eqn{V_i, i = 1, 2,} given
#' covariates \eqn{\mathbf{U}} and the other response
#' \eqn{V_j, j \in \{1, 2\} \setminus \{i\}}, or the marginal conditional
#' distribution \eqn{C_{V_i | \mathbf{U}}} of a response \eqn{V_i, i = 1, 2,}
#' given covariates \eqn{\mathbf{U}}.
#'
#' @param object an object of class \code{bivinereg}.
#' @param newdata data.frame of new response and covariate values for which to
#'   compute the marginal conditional distribution.
#' @param margin integer; the margin (1 or 2) for which to calculate the
#'   conditional CDF.
#' @param inc_resp boolean; if set to \code{TRUE} includes the other response,
#'   which is not \code{margin} in the conditioning set to calculate
#'   \eqn{C_{V_i | V_j, \mathbf{U}}}.
#' @param cores integer; the number of cores to use for computations.
#'
#' @examples
#' # load sample data
#' data(data)
#'
#' # fit vine regression model
#' (fit <- bivinereg(cbind(U1, U4) ~ U2 + U3 + U5 + U6,
#'                   data = data,
#'                   family_set = "parametric",
#'                   selcrit = "bic"))
#'
#' # calculate the marginal conditional CDF of U1 given U4 and the covariates
#' mccdf(fit, data[1:5,], margin = 1, inc_resp = FALSE)
#'
#' # calculate the marginal conditional CDF of U1 given the covariates
#' mccdf(fit, data[1:5,], margin = 1, inc_resp = TRUE)
#'
#' @export
mccdf <- function(object, newdata, margin, inc_resp, cores = 1) {
  newdata <- prepare_newdata(newdata, object, use_response = TRUE)
  newdata <- to_uscale(newdata, object$margins)
  if (inc_resp) {
    return(cond_m_dist_cpp(as.matrix(newdata), object$vine, margin, cores))
  } else {
    return(cond_m2_dist_cpp(as.matrix(newdata), object$vine, margin - 1, cores))
  }
}

#' Bivariate conditional CDF
#'
#' Calculates the bivariate conditional distribution \eqn{C_{V_1 V_2 | \mathbf{U}}}
#' of the responses \eqn{(V_1, V_2)} given covariates \eqn{\mathbf{U}}.
#'
#' @param object an object of class \code{bivinereg}.
#' @param newdata data.frame of new response and covariate values for which to
#'   compute the bivariate conditional distribution.
#' @param cores integer; the number of cores to use for computations.
#'
#' @importFrom stats integrate
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
#' # calculate the bivariate conditional CDF
#' biccdf(fit, data[1:5,])
#'
#' @export
biccdf <- function(object, newdata, cores = 1) {
  newdata <- prepare_newdata(newdata, object, use_response = TRUE)
  newdata <- as.data.frame(to_uscale(newdata, object$margins))
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

#' #' Bivariate conditional probability
#' #'
#' #' Calculates the bivariate conditional probability
#' #' \deqn{P(v_{1, min} \leq V_1 \leq v_{1, max},
#' #' v_{2, min} \leq V_2 \leq v_{2, max} | \mathbf{U} = \mathbf{u})}
#' #' of the responses \eqn{(V_1, V_2)} given covariates with values
#' #' \eqn{\mathbf{U} = \mathbf{u}} in the rectangle
#' #' \eqn{[v_{1, min}, v_{1, max}] \times [v_{2, min}, v_{2, max}]}.
#' #'
#' #' @param object an object of class \code{bivinereg}.
#' #' @param newdata data.frame of new covariate values for which to compute
#' #'   the bivariate conditional probability.
#' #' @param v1_min vector of lower limits for first response.
#' #' @param v1_max vector of upper limit for first response.
#' #' @param v2_min vector of lower limit for second response.
#' #' @param v2_max vector of upper limit for second response.
#' #' @param cores integer; the number of cores to use for computations.
#' #'
#' #' @details
#' #' This function double integrates the bivariate conditional density using the
#' #' function \code{integral2} from the \code{pracma} package.
#' #'
#' #' @seealso \code{\link{bicprob2}}
#' #' @importFrom pracma integral2
#' #'
#' #' @examples
#' #' # load sample data
#' #' data(data)
#' #'
#' #' # fit vine regression model
#' #' (fit <- bivinereg(cbind(U1, U4) ~ U2 + U3 + U5 + U6,
#' #'                   data = data,
#' #'                   family_set = "parametric",
#' #'                   selcrit = "bic",
#' #'                   uscale = TRUE))
#' #'
#' #' # calculate the conditional probability that 0.3 <= U1 <= 0.6,
#' #' # 0.3 <= U4 <= 0.6, given the covariates
#' #' bicprob(fit, data[1:5,], rep(0.3, 5), rep(0.6, 5), rep(0.3, 5), rep(0.6, 5))
#' #'
#' #' @export
#' bicprob <- function(object, newdata, v1_min, v1_max, v2_min, v2_max, cores = 1)
#' {
#'   newdata <- prepare_newdata(newdata, object)
#'   newdata <- as.data.frame(to_uscale(newdata, object$margins))
#'   out <- NA
#'   for (i in 1:nrow(newdata)) {
#'     out[i] <- integral2(bicprob_integrant,
#'                         v1_min[i], v1_max[i], v2_min[i], v2_max[i],
#'                         newdata = newdata[i,], object = object, cores = cores)$Q
#'   }
#'   return(out)
#' }
#'
#' #' @noRd
#' bicprob_integrant <- function(v1, v2, newdata, object, cores) {
#'   out <- matrix(NA, dim(v1)[1], dim(v1)[2])
#'   for (i in 1:dim(v1)[2]) {
#'     newdata_ <- as.matrix(cbind.data.frame(v1[, i], v2[, i],
#'                                            newdata[rep(1, dim(v1)[1]),],
#'                                            row.names = NULL))
#'     out[, i] <- cond_bi_dens_cpp(newdata_, object$vine, cores)
#'   }
#'   return(out)
#' }
#'
#' #' Bivariate conditional probability
#' #'
#' #' Calculates the bivariate conditional probability
#' #' \deqn{P(v_{1, min} \leq V_1 \leq v_{1, max},
#' #' v_{2, min} \leq V_2 \leq v_{2, max} | \mathbf{U} = \mathbf{u})}
#' #' of the responses \eqn{(V_1, V_2)} given covariates with values
#' #' \eqn{\mathbf{U} = \mathbf{u}} in the rectangle
#' #' \eqn{[v_{1, min}, v_{1, max}] \times [v_{2, min}, v_{2, max}]}.
#' #'
#' #' @param object an object of class \code{bivinereg}.
#' #' @param newdata data.frame of new covariate values for which to compute
#' #'   the bivariate conditional probability.
#' #' @param v1_min vector of lower limits for first response.
#' #' @param v1_max vector of upper limit for first response.
#' #' @param v2_min vector of lower limit for second response.
#' #' @param v2_max vector of upper limit for second response.
#' #' @param cores integer; the number of cores to use for computations.
#' #'
#' #' @details
#' #' This function double integrates the bivariate conditional density using the
#' #' function \code{hcubature} from the \code{cubature} package.
#' #'
#' #' @seealso \code{\link{bicprob}}
#' #'
#' #'
#' #' @importFrom cubature hcubature
#' #'
#' #' @examples
#' #' # load sample data
#' #' data(data)
#' #'
#' #' # fit vine regression model
#' #' (fit <- bivinereg(cbind(U1, U4) ~ U2 + U3 + U5 + U6,
#' #'                   data = data,
#' #'                   family_set = "parametric",
#' #'                   selcrit = "bic",
#' #'                   uscale = TRUE))
#' #'
#' #' # calculate the conditional probability that 0.3 <= U1 <= 0.6,
#' #' # 0.3 <= U4 <= 0.6, given the covariates
#' #' bicprob2(fit, data[1:5,], rep(0.3, 5), rep(0.6, 5), rep(0.3, 5), rep(0.6, 5))
#' #'
#' #' @export
#' bicprob2 <- function(object, newdata, v1_min, v1_max, v2_min, v2_max, cores = 1)
#' {
#'   newdata <- prepare_newdata(newdata, object)
#'   newdata <- as.data.frame(to_uscale(newdata, object$margins))
#'   out <- NA
#'   for (i in 1:nrow(newdata)) {
#'     out[i] <- hcubature(bicprob_integrant2,
#'                         c(v1_min[i], v2_min[i]), c(v1_max[i], v2_max[i]),
#'                         newdata = newdata[i,], object = object,
#'                         cores = cores)$integral
#'   }
#'   return(out)
#' }
#'
#' #' @noRd
#' bicprob_integrant2 <- function(v, newdata, object, cores) {
#'   newdata_ <- as.matrix(cbind.data.frame(v[1], v[2], newdata[1,],
#'                                          row.names = NULL))
#'   return(cond_bi_dens_cpp(newdata_, object$vine, cores))
#' }
