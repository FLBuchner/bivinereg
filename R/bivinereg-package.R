## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib bivinereg, .registration = TRUE
## usethis namespace: end
NULL

#' Something
#' @param i something
#' @export
test <- function(i) {
  do_it(2*i)
}
