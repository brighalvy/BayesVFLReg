#' Runs MCMC for BayesVFLReg
#'
#' @description Takes in the globally sketched data and performs the MCMC for the High-Dimensional Sparse Reduced Rank Regression.
#'
#' @param Phi_Y An m x Q matrix containing the globally sketched response matrix. (elements must be numeric)
#' @param Phi_X A m x P matrix containing the globally sketched covariate matrix. (elements must be numeric)
#' @param n_samps Number of desired post burn-in samples.
#' @param burnin Number of burn-in iterations, default 500
#'
#' @return A named list containing the following components:
#' * **B**: An array of the posterior samples for the B matrix (n_samps x P x Q).
#' * **A**: An array of the posterior samples for the A matrix (n_samps x Q x Q).
#' * **sigma**: A matrix of posterior samples of sigma^2 (n_samps x Q).
#' * **lambda**: The posterior local shrinkage parameters.
#' * **tau**: The posterior global shrinkage parameters.
#'
#' @export
BayesVFLReg <- function(Phi_Y, Phi_X, nsamps,
                  burnin = 500){
  ## Add checks:
  if(is.null(dim(Phi_Y)) | length(dim(Phi_Y)) != 2){
    stop("The 'Phi_Y' argument must be a matrix of 2-D array.")
  }
  if(!is.numeric(Phi_Y)){
    stop("The 'Phi_Y' argument must contain only numeric variables.")
  }
  if(is.null(dim(Phi_X)) | length(dim(Phi_X)) != 2){
    stop("The 'Phi_X' argument must be a matrix of 2-D array.")
  }
  if(!is.numeric(Phi_X)){
    stop("The 'Phi_X' argument must contain only numeric variables.")
  }


  # Get m, p, q
  m <- nrow(Phi_X)
  p <- ncol(Phi_X)
  q <- ncol(Phi_Y)
  # Initialize y and Phi_X
  y <- c(t(Phi_Y))
  # Initialize parameter storage:
  A <- array(NA, dim = c(nsamps + burnin, q, q))
  lambda <- matrix(NA, nrow = nsamps + burnin, ncol = p * q)
  tau <- sigma <- matrix(NA, nrow = nsamps + burnin, ncol = q)
  B <- array(NA, dim = c(nsamps + burnin, p, q))
  # Set starting values for parameters with random draws from respective priors:
  A[1, , ] <- matrix(rnorm(q^2), nrow = q, ncol = q)
  lambda[1, ] <- extraDistr::rhcauchy(p * q, 1)
  tau[1, ] <- extraDistr::rhcauchy(q, 1)
  B[1, , ] <- matrix(rnorm(p * q, 0, sd = (lambda[1, ] * rep(tau[1, ], p))), nrow = p, ncol = q)
  sigma[1, ] <- invgamma::rinvgamma(q, 2, 2)
  # Run Gibbs Sampler:
  for(b in 2:(nsamps + burnin)){
    # Create Lambda & Sigma_tilde
    Lambda <- diag(lambda[b - 1, ]^2 * rep(tau[b - 1, ]^2, p))
    Sigma_tilde <- diag(rep(sigma[b - 1, ], m))
    # Draw B:
    B[b, , ] <- matrix(sample_B(Lambda, Sigma_tilde, Phi_X, A[b - 1, , ], y, p,
                                q, m), nrow = p, ncol = q, byrow = TRUE)
    # Draw A:
    A[b, , ] <- matrix(sample_A(Sigma_tilde, Phi_X, B[b, , ],
                                y, q), nrow = q, ncol = q, byrow = FALSE)
    # Draw sigma
    sigma[b, ] <- draw_sigma(m, q, Phi_Y, Phi_X, B[b, , ],
                             A[b, , ])
    # Draw Lambda
    lambda[b, ] <- draw_lambda(lambda[b - 1, ], tau[b - 1, ],
                               B[b, , ], p, q)

    # Draw tau
    tau[b, ] <- draw_tau(lambda[b, ], tau[b - 1, ],
                         B[b, , ], p, q)

  }
  # Return posterior samples:
  return(list(B = B[-c(1:burnin), , ],
              A = A[-c(1:burnin), , ],
              sigma = sigma[-c(1:burnin), ],
              lambda = lambda[-c(1:burnin), ],
              tau = tau[-c(1:burnin), ]))
}
