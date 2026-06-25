
# BayesVFLReg

<!-- badges: start -->
<!-- badges: end -->

[![R-CMD-check](https://github.com/brighalvy/BayesVFLReg/workflows/R-CMD-check/badge.svg)](https://github.com/brighalvy/BayesVFLReg/actions)
`BayesVFLReg` is an R package for implementing Bayesian Vector Functional Linear Regression. The package utilizes a highly optimized Gibbs sampler written in C++ via `Rcpp` and `RcppArmadillo` to efficiently sample from high-dimensional posterior distributions with shrinkage priors.

## Features

* **Fast MCMC Sampling:** Core sampling steps (`sample_A` and `sample_B`) are compiled in C++ for maximum performance.
* **Robust Priors:** Built-in support for global-local shrinkage structures using Half-Cauchy and Truncated Exponential distributions.
* **High-Dimensional Efficiency:** Optimized matrix operations specifically designed to handle large Kronecker products without unnecessary memory allocation.

## Installation

You can install the development version of BayesVFLReg from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("brighalvy/BayesVFLReg")
```

## Key Functions

| Function | Description |
| :--- | :--- |
| `generate_Phi()`  | Generates the random sketching matrix Phi. |
| `BayesVFLReg()`  | Takes in the Globally sketched data as arguments and runs the MCMC algorithm. |
| `post_process()`  | Takes posterior samples or estimates to run the post processing procedure and returns an indicator of non-zero covariates. |


## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(BayesVFLReg)
## Generate Small Sample data:
## Set q, p, n
  q <- 5
  p <- 10
  n <- 1,000
## Set number of non-zero covariates:
s <- 3
## Simulate X with 0.5 correlation among covariates:
  cor = 0.5
  Sigma_X <- matrix(cor, nrow = p, ncol = p)
  diag(Sigma_X) <- 1
  Sigma_X_c <- chol(Sigma_X)
  draw_x <- function(x){Sigma_X_c %*% rnorm(p)}
  X <- t(sapply(1:n, draw_x))
## Simulate Sigma:
  sigma <- rinvgamma(q, 2, 2)
  
## Simulate A:
  A <- matrix(rnorm(q^2), nrow = q, ncol = q)
  
## Simulate B:
  values <- rnorm(s*q)
  row_index <- sample(1:p, s)
  B <- matrix(0,
              nrow = p, ncol = q)
  B[row_index,] <- values
  
## Get C:
  C <- B %*% t(A)
  
## Get Y:
  Sigma_E <- diag(sigma)
  Sigma_E_full <- bdiag(replicate(n, Sigma_E, simplify = FALSE))
  Sigma_E_c <- chol(Sigma_E_full)
  y <- c(t(X %*% C)) + Sigma_E_c %*% rnorm(n * q)
  Y <- matrix(y, nrow = n, ncol = q, byrow = TRUE)
  
## Generate Phi (m = 200):
Phi <- gnerate_Phi(n , m)

## Sketch the Data:
Phi_X <- Phi %*% X
Phi_Y <- Phi %*% Y

## Run analysis:
mcmc_results <- BayesVFLReg(Phi_X = Phi_X, Phi_Y = Phi_Y,
                            nsampe = 1000, burnin = 10)
                            
## Get which values are non-zero:
non_zero <- post_process(samps = mcmc_results)
```

