
// [[Rcpp::depends(RcppArmadillo)]]
#include "Sample_A_B.h"
// [[Rcpp::export]]
arma::vec sample_B(const arma::mat& Lambda, 
                   const arma::mat& Sigma_tilde, 
                   const arma::mat& Phi_X, 
                   const arma::mat& A, 
                   const arma::vec& y,
                   int p, int q, int m) {
  
  int pq = p * q;
  int mq = m * q;
  
  // Step 1: Draw u from N(0, Lambda), delta from N(0, I_mq) 
  //arma::mat Lambda_chol = arma::chol(Lambda, "lower");
  arma::mat Lambda_sqrt = arma::diagmat(arma::sqrt(Lambda.diag()));
  arma::vec u = Lambda_sqrt * arma::randn<arma::vec>(pq);
  arma::vec delta = arma::randn<arma::vec>(mq);
  
  // Step 2: Set v = X_tilde * u + delta
  arma::vec sigma_inv_sqrt_vec = 1.0 / arma::sqrt(Sigma_tilde.diag());
  arma::mat X_tilde = arma::kron(Phi_X, A);
  X_tilde.each_col() %= sigma_inv_sqrt_vec;
  arma::vec v = X_tilde * u + delta;
  
  // Step 3: Solve for w
  arma::vec y_tilde = sigma_inv_sqrt_vec % y;
  arma::vec Lambda_diag = Lambda.diag();
  arma::mat Lambda_X = X_tilde;
  Lambda_X.each_row() %= Lambda_diag.t();
  arma::mat XtLXt = Lambda_X * X_tilde.t();
  XtLXt.diag() += 1.0;
  //arma::vec w = arma::pinv(XtLXt) * (y_tilde - v);
  arma::vec w = arma::solve(XtLXt, y_tilde - v, arma::solve_opts::fast + arma::solve_opts::likely_sympd);
  
  // Step 4: Get beta
  arma::vec beta = u + Lambda * X_tilde.t() * w;
  
  // Step 5: Return B matrix
  // arma::mat B(beta.memptr(), q, p, false); // no copy, column-major
  return beta; // Transpose to get byrow=TRUE effect
}

// [[Rcpp::export]]
arma::vec sample_A(const arma::mat& Sigma_tilde,
                   const arma::mat& Phi_X,
                   const arma::mat& B,
                   const arma::vec& y,
                   int q) {
  
  // int m = Phi_X.n_rows;
  int q2 = q * q;
  
  // Step 1: Compute y_tilde
  // arma::mat Sigma_inv_sqrt = arma::diagmat(1.0 / arma::sqrt(Sigma_tilde.diag()));
  arma::vec sigma_inv_sqrt_vec = 1.0 / arma::sqrt(Sigma_tilde.diag());
  arma::vec y_tilde = sigma_inv_sqrt_vec % y;

  // Step 2: Compute X_star
  arma::mat X_star = arma::kron(Phi_X * B, arma::eye(q, q));
  X_star.each_col() %= sigma_inv_sqrt_vec;

  // Step 3: Compute Omega_A and its inverse
  arma::mat Omega_A = X_star.t() * X_star ;
  Omega_A.diag() += 1.0;
  arma::mat Omega_A_inv = arma::inv_sympd(Omega_A);  // Or use pinv if potentially ill-conditioned

  // Step 4: Draw a using spectral decomposition of Omega_A_inv
  // Eigen decomposition: Sigma = Q * D * Q^T
  arma::vec eigval;
  arma::mat eigvec;
  
  arma::eig_sym(eigval, eigvec, Omega_A_inv);
  
  // Ensure non-negative eigenvalues (numerical safety)
  eigval = arma::clamp(eigval, 0, arma::datum::inf);
  
  // Square root of eigenvalues matrix
  arma::mat D_sqrt = arma::diagmat(arma::sqrt(eigval));
  
  // Construct the transformation matrix: Q * sqrt(D)
  arma::mat transform = eigvec * D_sqrt;
  arma::vec a = Omega_A_inv * X_star.t() * y_tilde + transform * arma::randn<arma::vec>(q2);

  // Step 5: Reshape and return A
  //arma::mat A(a.memptr(), q, q, false);  // No copy, column-major
  return a;  // Transpose for byrow = TRUE
}
