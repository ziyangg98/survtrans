#' Ancillary arguments for controlling survtrans fitting
#' @param abstol the absolute tolerance for ADMM primal/dual residuals.
#' Default is 1e-4.
#' @param reltol the relative tolerance for ADMM primal/dual residuals.
#' Default is 1e-3.
#' @param fdev the minimum fractional change of the augmented Lagrangian for
#' convergence. The algorithm stops when
#' \eqn{|L^{k} - L^{k-1}| / (|L^{k-1}| + 1) < fdev}, where \eqn{L} is the
#' augmented Lagrangian. This provides a fallback when primal-dual residuals
#' oscillate under non-convex penalties. Default is 1e-5.
#' @param maxit the maximum number of iterations for the proposed algorithm.
#'  Default is 300.
#' @param verbose a logical value indicating whether to print messages
#'  during the fitting process. Default is \code{FALSE}.
#' @return A list with components \code{abstol}, \code{reltol}, \code{fdev},
#'  \code{maxit}, and \code{verbose}.
#' @keywords internal
survtrans_control <- function(
  abstol = 1e-4, reltol = 1e-3, fdev = 1e-5, maxit = 300,
  verbose = FALSE
) {
  if (abstol <= 0) stop("Invalid absolute tolerance")
  if (reltol <= 0) stop("Invalid relative tolerance")
  if (fdev <= 0) stop("Invalid value for fdev")
  if (maxit < 0) stop("Invalid value for iterations")
  if (!is.logical(verbose)) stop("Invalid value for verbose")
  list(
    abstol = abstol, reltol = reltol, fdev = fdev,
    maxit = as.integer(maxit), verbose = verbose
  )
}
