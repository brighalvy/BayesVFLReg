#' Post processing
#'
#' @description Takes posterior draws to
#'
#' @param samps The output object of BayesVFLReg() of this package (default NULL).
#' @param B Posterior draws of the B matrix (nsamps x P x Q) (default NULL).
#' @param A Posterior draws of the A matrix (nsamps x Q x Q) (default NULL).
#' @param B_est Estimate of the B matrix (P x Q)  (default NULL).
#' @param A_est Estimate of the A matrix (Q x Q) (default NULL).
#' @param C_est Estimate of the C (coefficient matrix) (P x Q)  (default NULL).
#'
#' @return A vector of length P of TRUE indicating a significant coefficient.
#'
#'
#' @export
post_process <- function(samps = NULL, B = NULL,
                         A = NULL, B_est = NULL,
                         A_est = NULL, C_est = NULL){
  if(!is.null(samps)){
    B <- samps$B
    A <- samps$A
  }
  if(!is.null(A) & !is.null(B)){
    if(is.null(dim(A)) | length(dim(A)) != 3){
      stop("A must be a 3d Array")
    }
    if(is.null(dim(B)) | length(dim(B)) != 3){
      stop("B must be a 3d Array")
    }
    q <- ncol(A[1,,])
    p <- ncol(dat$X)
    iter <- dim(B)[1]
    C_draws <- array(NA, dim = c(iter, p, q))
    for(i in 1:iter){
      C_draws[i,,] <- B[i, ,] %*% t(A[i,,])
    }
    C_est <- colMeans(C_draws, dims = 1)
  }
  if(!is.null(A_est) & !is.null(B_est)){
    C_est <- B_est %*% t(A_est)
  }
  if(is.null(C_est)){
    stop("Invalid input")
  }
  clust_res <- mclust::Mclust((C_est), 2)
  min_group <- which.min(colSums(abs(clust_res$parameters$mean)))
  prob <- clust_res$z[,min_group]
  # find significant coefficients (non_zero)
  non_zero_index <- which(prob < .05)
  return(non_zero_index)
}
