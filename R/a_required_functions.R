## Required functions:

## Draw Sigma function:
draw_sigma <- function(m, q, Y, Phi_X, B, A){
  diff_mat <- Y - Phi_X %*% B %*% t(A)
  y_stars <- colSums(diff_mat^2)
  sigma <- invgamma::rinvgamma(q, m / 2, rate = y_stars / 2)
  return(sigma)
}

## Draw Lambda Function:
draw_lambda <- function(lambda, tau, B, p, q){
  eta <- 1/(lambda^2)
  beta <- c(t(B))
  mu <- beta^2 / rep(tau^2, p)
  ## Sample u:
  u <- runif(p * q, 0, 1/(1 + eta))
  ## Sample new eta:
  params <- cbind(mu, u)
  eta_new <- apply(params, 1, \(x) {TruncExpFam::rtruncexp(1, rate = 2/x[1], a = 0, b = (1 - x[2]) / x[2],
                                              faster = TRUE)})
  return(sqrt(1 / eta_new))
}

