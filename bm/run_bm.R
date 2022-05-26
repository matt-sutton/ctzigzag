###### Config ####
iters <- 50000 ## Number of samples / number of events (simulations used 1e6)
warmup <- 10000  ## Number samples/events to burn (simulations used 1e3)

alpha_1 <- 0.2
library(rstan)

## Load the data (original data located in the .nzp files)
# also see https://github.com/matt-graham/continuously-tempered-hmc/tree/master/data/gaussian-bmr)

load("bm/bm_processed.Rdata")

## Run ct hmc
stan_fit_eval <- stan("bm/model_eval.stan", data = data,
                      iter = 1, chains = 1)
target <- function(x){
  stan_ev <- grad_log_prob(stan_fit_eval, c(x,1))
  d_log_q <- as.numeric(stan_ev)[1:p]
  log_q <- attr(stan_ev, "log_prob")
  
  return(list(log_q = log_q, d_log_q = d_log_q))
}
temper <- function(x){
  stan_ev <- grad_log_prob(stan_fit_eval, c(x,0))
  d_log_q <- as.numeric(stan_ev)[1:p]
  log_q <- attr(stan_ev, "log_prob")
  
  return(list(log_q = log_q, d_log_q = d_log_q))
}

source("code/temper_zigzag_hess.R")
## Evaluate H^+
H_p <- Reduce("+",lapply(1:dim(data$q)[1], function(j) {
  qq <- tcrossprod(data$q[j,])
  qq[qq<0] <-0
  return(qq)
}))
## Evaluate H^-
H_n <- Reduce("+",lapply(1:dim(data$q)[1], function(j) {
  qq <- tcrossprod(data$q[j,])
  qq[qq>0] <-0
  qq
}))

H <- pmax(H_p, diag(1,p,p) + abs(H_n))## H upper bounding the hessian

theta_0 <- rep(1,p)
hess_q1 <- matrix(H%*%theta_0[1:p],  nrow = 1)
hess_q0 <- matrix(abs(solve(tcrossprod(data$chol_sigma)))%*%theta_0[1:p], nrow = 1, byrow = T)
(hess <- abs(matrix(c(hess_q0,
                      hess_q1), nrow = 2,ncol = p, byrow = T)))
rownames(hess) <- c("q0", "q1")

x_init <- data$mu
theta_init <- rep(1,p)

## Standard Zig-Zag
set.seed(1);zigzag_a1 <- zigzag_temp(max_events = iters, 
                                      x0 = c(x_init, 1),
                                      theta0 = c(theta_init, 0),
                                      alpha = 1, tau_max = 1,
                                      poly_order = length(poly_coef), echo = F,
                                      poly_coef = poly_coef)

## Tempered Zig-Zag
set.seed(1);zigzag_a.2 <- zigzag_temp(max_events = iters, 
                                      x0 = c(x_init, 0.1),
                                      theta0 = c(theta_init, 1),
                                      alpha = .2, tau_max = 1,
                                      poly_order = length(poly_coef), echo = F,
                                      poly_coef = poly_coef)

## Importance sampling tempered Zig-Zag
## Adjust the log_q definition - equivalent to kappa(t) = zeta^(1-t)
## log_zeta given in data file (from a variational approximation)
temper <- function(x){
  stan_ev <- grad_log_prob(stan_fit_eval, c(x,0))
  d_log_q <- as.numeric(stan_ev)[1:p]
  log_q <- attr(stan_ev, "log_prob")
  
  return(list(log_q = log_q + data_cts$log_zeta, d_log_q = d_log_q))
}

set.seed(1);zigzag_is <- zigzag_temp(max_events = iters, 
                                     x0 = c(x_init,0.1),
                                     theta0 = c(theta_init,1),
                                     alpha = 0, tau_max = 1,
                                     poly_order = 3, echo = F,
                                     poly_coef = rep(0,2))

## Plotting for the 26 and 27 bivariate marginals also the inverse temperature
plot_pdmp_multiple(list(is = zigzag_is, 
                        point_mass = zigzag_a.2, 
                        standard = zigzag_a1), 
                   nsamples = 5e3, inds = 1:1, pch = 20,
                   coords = c(p,p-1, p+1))

