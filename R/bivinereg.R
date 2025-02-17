#' Y-vine regression models
#'
#' Sequential estimation of a regression Y-vine as described in Tepegjozova and
#' Czado (2022).
#'
#' @param formula an object of class "formula"; same as [lm()].
#' @param data data frame (or object coercible by [as.data.frame()]) containing
#'   the variables in the model.
#' @param family_set see `family_set` argument of [rvinecopulib::bicop()].
#' @param selcrit selection criterion based on conditional log-likelihood.
#'   \code{"loglik"} (default) imposes no correction; other choices are
#'   \code{"aic"} and \code{"bic"}.
#' @param par_1d list of options passed to [kde1d::kde1d()], must be one value
#'   for each margin, e.g. `list(xmin = c(0, 0, NaN))` if the responses have
#'   non-negative support.
#' @param weights optional vector of weights for each observation.
#' @param cores integer; the number of cores to use for computations.
#' @param ... further arguments passed to [rvinecopulib::bicop()].
#' @param uscale if TRUE, bivinereg assumes that marginal distributions have been
#'   taken care of in a preliminary step.
#'
#' @return An object of class bivinereg. It is a list containing the elements
#'   \describe{ \item{formula}{the formula used for the fit.}
#'   \item{selcrit}{criterion used for variable selection.}
#'   \item{model_frame}{the data used to fit the regression model.}
#'   \item{margins}{list of marginal models.}
#'   \item{vine}{an [rvinecopulib::vinecop_dist()] object containing the fitted
#'   Y-vine.} \item{stats}{fit statistics such as conditional
#'   log-likelihood/AIC/BIC and p-values for each variable's contribution.}
#'   \item{order}{order of the covariates chosen by the variable selection
#'   algorithm.} \item{selected_vars}{indices of selected variables.} }
#'   `summary.bivinereg()` shows the contribution of each selected variable with
#'   the associated p-value derived from a likelihood ratio test.
#'   `summary_yvine()` gives a more detailed and flexible summary of the fitted
#'   Y-vine.
#'
#' @references Tepegjozova and Czado (2022), Bivariate vine copula based
#' regression, bivariate level and quantile curves,
#' https://arxiv.org/abs/2205.02557
#'
#' @examples
#' # load sample data
#' data(data)
#'
#' # fit vine regression model
#' (fit <- bivinereg(cbind(U1,U4) ~ U2 + U3 + U5 + U6,
#'                   data = data,
#'                   family_set = "parametric",
#'                   selcrit = "bic",
#'                   uscale = TRUE))
#'
#' # inspect model
#' summary(fit)
#'
#' @export
#'
#' @importFrom kde1d kde1d pkde1d
#' @importFrom stats model.frame logLik model.extract
#' @importFrom utils modifyList
#' @importFrom rvinecopulib bicop vinecop
#' @importFrom Rcpp sourceCpp
#' @useDynLib bivinereg, .registration = TRUE
bivinereg <- function(formula, data, family_set = "parametric", selcrit = "aic",
                      par_1d = list(), weights = numeric(),
                      cores = 1, ..., uscale = FALSE) {
  # remove unused variables
  if (!missing(data)) {
    mf <- model.frame(formula, data)
  } else {
    mf <- model.frame(formula, parent.frame())
  }
  if (!is.numeric(mf[[1]][,1]) | !is.numeric(mf[[1]][,2]))
    stop("responses must be numeric")
  if (any(sapply(mf, is.factor)) && uscale)
    stop("factors are not allowed with uscale = TRUE")

  mfx <- expand_factors(mf)
  colnames(mfx)[1:2] <- colnames(model.extract(mf, "response"))
  d <- ncol(mfx)
  var_types <- rep("c", d)
  var_types[sapply(mfx, is.ordered)] <- "d"

  ## prepare fit controls (little hack calling bicop() for all checks)
  arg <- list(
    data = t(c(0.5, 0.5)),
    family_set = family_set,
    selcrit = selcrit,
    cores = cores,
    par_method = "mle",
    nonpar_method = "quadratic",
    mult = 1,
    psi0 = 0.9,
    presel = TRUE,
    keep_data = FALSE
  )
  ctrl <- do.call(
    bicop,
    modifyList(arg, list(...))
  )$controls
  ctrl$weights <- numeric()
  # We assume no order given
  # Can add given order later
  if (!uscale) {
    par_1d <- process_par_1d(mfx, par_1d)
    margins <- fit_margins_cpp(prep_for_kde1d(mfx),
                               xmin = par_1d$xmin,
                               xmax = par_1d$xmax,
                               type = par_1d$type,
                               mult = par_1d$mult,
                               bw = par_1d$bw,
                               deg = par_1d$deg,
                               weights = weights,
                               num_threads = cores)
    margins <- finalize_margins(margins, mfx)
    u <- to_uscale(mfx, margins)
  } else {
    margins <- lapply(1:d, function(x) list(edf = NA, loglik = NA))
    u <- as.matrix(mfx)
  }

  args <- modifyList(
    ctrl,
    list(data = u, var_types = var_types, cores = cores, weights = weights)
  )

  fit <- do.call(select_yvine_cpp, args)

  if (!uscale)
    margins <- margins[c(1:2, sort(fit$selected_vars))] # other margins useless

  finalize_bivinereg_object(
    formula = formula,
    selcrit = selcrit,
    model_frame = mf,
    margins = margins,
    vine = fit$vine,
    selected_vars = fit$selected_vars,
    var_nms = colnames(mfx)
  )
}

#' @noRd
#' @importFrom stats pchisq
#' @importFrom rvinecopulib as_rvine_structure
finalize_bivinereg_object <- function(formula, selcrit, model_frame, margins, vine,
                                    selected_vars, var_nms) {
  vine$names <- c(var_nms[1], var_nms[2], var_nms[sort(selected_vars)])
  nobs <- nrow(model_frame)
  vine$nobs <- nobs
  var_edf <- c(
    margins[[1]]$edf + margins[[2]]$edf,
    vapply(vine$pair_copulas[1:length(selected_vars)], function(pcs) {pcs[[1]]$npars + pcs[[2]]$npars}, numeric(1))
  )
  var_cll <- c(
    margins[[1]]$loglik + margins[[2]]$loglik,
    vapply(vine$pair_copulas[1:length(selected_vars)], function(pcs) {pcs[[1]]$loglik + pcs[[2]]$loglik}, numeric(1))
  )

  var_caic <- -2 * var_cll + 2 * var_edf
  var_cbic <- -2 * var_cll + log(nobs) * var_edf
  var_p_value <- suppressWarnings(
    pchisq(2 * var_cll, var_edf, lower.tail = FALSE)
  )
  var_p_value[1] <- NA
  cll <- sum(var_cll, na.rm = TRUE)
  edf <- sum(var_edf, na.rm = TRUE)
  caic <- sum(var_caic, na.rm = TRUE)
  cbic <- sum(var_cbic, na.rm = TRUE)

  stats <- list(
    nobs = nobs,
    edf = edf,
    acll = cll,
    acaic = caic,
    acbic = cbic,
    var_edf = var_edf,
    var_cll = var_cll,
    var_caic = var_caic,
    var_cbic = var_cbic,
    var_p_value = var_p_value
  )

  out <- list(
    formula = formula,
    selcrit = selcrit,
    model_frame = model_frame,
    margins = margins,
    vine = vine,
    stats = stats,
    order = var_nms[selected_vars],
    selected_vars = selected_vars
  )
  class(out) <- "bivinereg"
  out
}

check_order <- function(order, var_nms) {
  # stopifnot(length(order) > 0)
  # if (!all(order %in% var_nms)) {
  #   stop(
  #     "unknown variable name in 'order'; ",
  #     "allowed values: '", paste(var_nms[-1], collapse = "', '"), "'."
  #   )
  # }
  # if (any(order == var_nms[1])) {
  #   stop(
  #     "response variable '", var_nms[1],
  #     "' must not appear in 'order'."
  #   )
  # }
}
