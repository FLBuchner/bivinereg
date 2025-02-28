#' Bivariate conditional highest density region
#'
#' Plot the highest density region (HDR) of the bivariate conditional density for a set of covariate values.
#'
#' @param object an object of class \code{bivinereg}.
#' @param newdata data.frame with one row of covariate values.
#' @param scale One of \code{"xscale", "uscale", "zscale"}. The scale on which the HDR should be plotted. \code{scale = "xscale"} only works if \code{object} was fitted with using \code{uscale = FALSE} in [bivinereg()].
#' @param xlim Vector of two integers. Interval on the x-axis in which the HDR will be calculated. Automatically set to \code{c(0, 1)} if \code{scale = "uscale"}.
#' @param ylim Vector of two integers. Interval on the y-axis in which the HDR will be calculated. Automatically set to \code{c(0, 1)} if \code{scale = "uscale"}.
#' @param ... Further arguments to [geom_hdr_fun()].
#'
#' @returns A [ggplot2] object.
#'
#' @seealso [bivinereg()] [ggplot2] [geom_hdr_fun()]
#'
#'@examples
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
#' # Plot HDR
#' bichdr(fit, data[1,], scale = "zscale")
#'
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggdensity geom_hdr_fun
#'
bichdr <- function(object, newdata, scale = "xscale", xlim = c(-5, 5), ylim = c(-5, 5), ...) {

  if (!(scale %in% c("xscale", "uscale", "zscale"))) {
    stop("scale has to be one of 'xscale', 'uscale', 'zscale'!")
  }

  if (scale == "xscale" & any(sapply(object$margins, length) == 2)) {
    stop("scale = 'xscale' can only be set when uscale = FALSE was used in bivinereg!")
  }

  if (nrow(newdata) > 1) {
    warning(paste0("newdata has ", nrow(newdata), " rows. Only the first row will be used."))
    newdata <- newdata[1,]
  }

  newdata <- prepare_newdata(newdata, object, use_response = FALSE)

  if (scale == "xscale") {
    fun <- bicpdf_hdr_xscale
  } else if (scale == "uscale") {
    fun <- bicpdf_hdr_uscale
    xlim <- ylim <- c(0, 1)
  } else if (scale == "zscale") {
    fun <- bicpdf_hdr_zscale
  }

  args <- list(
    newdata = newdata,
    object = object
  )

  ggplot() +
    geom_hdr_fun(fun = fun, args = args, xlim = xlim, ylim = ylim, ...)
}

#' bivariate density for hdr on xscale
#' @noRd
bicpdf_hdr_xscale <- function(x, y, newdata, object, cores = 1) {
  newdata <- cbind.data.frame(x, y, newdata, row.names = NULL)
  colnames(newdata) <- object$vine$names

  dens_marg <- kde1d::dkde1d(newdata[, 1], object$margins[[1]]) *
    kde1d::dkde1d(newdata[, 2], object$margins[[2]])

  newdata <- to_uscale(newdata, object$margins)

  cond_bi_dens_cpp(as.matrix(newdata), object$vine, cores) * dens_marg
}

#' bivariate density for hdr on uscale
#' @noRd
bicpdf_hdr_uscale <- function(x, y, newdata, object, cores = 1) {
  newdata <- cbind.data.frame(0.5, 0.5, newdata)
  newdata <- data.frame(to_uscale(newdata, object$margins))

  newdata <- cbind.data.frame(x, y, newdata[, -c(1, 2)], row.names = NULL)
  colnames(newdata) <- object$vine$names

  cond_bi_dens_cpp(as.matrix(newdata), object$vine, cores)
}

#' bivariate density for hdr on zscale
#' @noRd
bicpdf_hdr_zscale <- function(x, y, newdata, object, cores = 1) {
  newdata <- cbind.data.frame(0.5, 0.5, newdata)
  newdata <- data.frame(to_uscale(newdata, object$margins))

  dens_marg <- stats::dnorm(x) * stats::dnorm(y)

  newdata <- cbind.data.frame(stats::pnorm(x), stats::pnorm(y), newdata[, -c(1, 2)], row.names = NULL)
  colnames(newdata) <- object$vine$names

  cond_bi_dens_cpp(as.matrix(newdata), object$vine, cores) * dens_marg
}
