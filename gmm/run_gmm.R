###### Config ####
iters <- 5000 ## Number of samples / number of events (simulations used 30000)
warmup <- 1000  ## Number samples/events to burn (simulations used 20000)

set.seed(1)
p <- 2
x_init <- rnorm(p, mean = 5, sd = 5)
theta_init <- rep(1, p)

## Stan is used for gradient calculations and log density evaluations not sampling

library(rstan)
load("gmm/gmm_data.Rdata")
model_eval <- stan("gmm/model_eval.stan", data = dat, iter = 1, chains = 1)

## Target and Temper return the target and base evaluations of the gradient 
## and density evaluations

target <- function(x){
  stan_ev <- grad_log_prob(model_eval, c(x,1))
  d_log_q <- as.numeric(stan_ev)[1:p]
  log_q <- attr(stan_ev, "log_prob")

  return(list(log_q = log_q, d_log_q = d_log_q))
}
temper <- function(x){
  stan_ev <- grad_log_prob(model_eval, c(x,0))
  d_log_q <- as.numeric(stan_ev)[1:p]
  log_q <- attr(stan_ev, "log_prob")

  return(list(log_q = log_q, d_log_q = d_log_q))
}

source("code/temper_zigzag_hess.R")
mean_diffs <- (apply(mu_mat, 1, max) - apply(mu_mat, 1, min))^2/4

C_1 <- 1/dat$sigma[1]
hess_q1 <- abs(C_1*(1 + C_1*mean_diffs))
hess_q0 <- 1/dat$sigma0
(hess <- abs(matrix(c(hess_q0,
                      hess_q1), 2,2, byrow = T)))
rownames(hess) <- c("q0", "q1")


## Warmup (burnin) iterations for the Zig-Zag process alpha = 0
set.seed(1);zigzag_warm <- zigzag_temp(max_events = warmup, 
                                       x0 = c(x_init,0.1),
                                       theta0 = c(theta_init,1),
                                       alpha = 0, tau_max = 1,
                                       poly_order = 3, echo = F,
                                       poly_coef = rep(0,2))
x_init <- tail(t(zigzag_warm$positions),1)
theta_init <- tail(t(zigzag_warm$thetas),1)

## Build approximation for kappa using an 8th order poly
zigzag_samples <- gen_samples(zigzag_warm$positions, zigzag_warm$times,
                              zigzag_warm$thetas, nsample = warmup, burn = 10)

phi_eval <- apply(zigzag_samples$samples,MARGIN = 2, phi)
pq<-path_quad(zigzag_samples$samples[2+1,], phi_eval)

t_vals_fit <- pq$t;  t_vals_log_z <- pq$log_z
t_vals_poly <- sapply(0:(8-1), function(s) t_vals_fit^s) 

poly_df <- data.frame(t_vals_poly, logz = t_vals_log_z)
poly_coef <- rev(coef(lm(logz ~ 0 + ., data = poly_df)))

# plot eval log_z and poly fit
plot(pq$t, pq$log_z, type = 'l')
lines(pq$t, pracma::polyval(poly_coef, pq$t), col ='red')


## Tempered Zig-Zag
set.seed(1);zigzag_a.2 <- zigzag_temp(max_events = iters, 
                                     x0 = x_init,
                                     theta0 = theta_init,
                                     alpha = .2, tau_max = 1,
                                     poly_order = length(poly_coef), echo = F,
                                     poly_coef = poly_coef)

## Importance sampling tempered Zig-Zag
set.seed(1);zigzag_is <- zigzag_temp(max_events = iters, 
                                     x0 = x_init,
                                     theta0 = theta_init,
                                     alpha = 0, tau_max = 1,
                                     poly_order = 3, echo = F,
                                     poly_coef = rep(0,2))

## Adjust the initial warmup for a=1 
## (this was done separately in the simulations)
a1_x_init <- c(x_init[-3], 1)
a1_theta_init <- c(theta_init[-3], 0)

## Standard Zig-Zag
set.seed(1);zigzag_alpha_a1 <- zigzag_temp(max_events = iters, 
                                           x0 = a1_x_init,
                                           theta0 = a1_theta_init,
                                           alpha = 1, tau_max = 1,
                                           poly_order = 3, echo = F,
                                           poly_coef = rep(0,2))

np <- 2e2
x_v <- seq(from = -1, to = 12, length.out = np)
y_v <- seq(from = -1, to = 12, length.out = np)
xy <- expand.grid(x=x_v, y=y_v)
f <- function(x){exp(target(x)$log_q)}
z <- matrix(apply(as.matrix(xy), 1, f), length(x_v), length(y_v))
par(mfrow = c(1,3))
image(x_v, y_v, z, las=1)

lines(t(zigzag_alpha_a1$positions)[,1:2],
      col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
legend('topleft', legend = c('Zig Zag'), col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5), lwd = 1)

image(x_v, y_v, z, las=1)
lines(t(zigzag_is$positions)[,1:2],
      col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
legend('topleft', legend = c('ct Zig Zag (IS)'), col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5), lwd = 1)

image(x_v, y_v, z, las=1)
lines(t(zigzag_a.2$positions)[,1:2],
      col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
legend('topleft', legend = c('ct Zig Zag'), col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5), lwd = 1)



