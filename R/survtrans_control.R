#' Ancillary arguments for controlling survtrans fitting
#' @param abstol the absolute tolerance for ADMM primal/dual residuals.
#' Default is 1e-4.
#' @param reltol the relative tolerance for ADMM primal/dual residuals.
#' Default is 1e-3.
#' @param nthreads the number of OpenMP threads used by a single model fit.
#' Default is 1.
#' @param maxit the maximum number of iterations for the proposed algorithm.
#'  Default is 300.
#' @param verbose a logical value indicating whether to print messages
#'  during the fitting process. Default is \code{FALSE}.
#' @return A list with components \code{abstol}, \code{reltol},
#'  \code{nthreads}, \code{maxit}, and \code{verbose}.
#' @keywords internal
survtrans_control <- function(
  abstol = 1e-4, reltol = 1e-3, nthreads = 1L, maxit = 300,
  verbose = FALSE
) {
  if (abstol <= 0) stop("Invalid absolute tolerance")
  if (reltol <= 0) stop("Invalid relative tolerance")
  if (!is.numeric(nthreads) || length(nthreads) != 1L || nthreads < 1) {
    stop("Invalid value for nthreads")
  }
  if (maxit < 0) stop("Invalid value for iterations")
  if (!is.logical(verbose)) stop("Invalid value for verbose")
  list(
    abstol = abstol, reltol = reltol,
    nthreads = as.integer(nthreads),
    maxit = as.integer(maxit), verbose = verbose
  )
}
