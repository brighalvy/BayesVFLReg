#' Generates Sketching matrix Phi
#'
#' @description Takes in n (number of samples) and m (desired sketch size) to generate Phi.
#'
#' @param n number of samples
#' @param m desired sketch size
#'
#' @return Phi: A m x n matrix.
#'
#'
#' @export
gnerate_Phi <- function(n, m){
  Phi <- matrix(rnorm(n*m, sd = sqrt(1/n)), nrow = m, ncol = n)
  return(Phi)
}
