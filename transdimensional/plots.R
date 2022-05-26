### Experiment 2 - Spike-slab separated modes
p <- 2;

## Target Likelihood (No data so no addition of log_q or d_log_q)
## Automatically adds a normal slab prior
target <- function(x){
  p <- length(x)
  d_log_q <- rep(0,p)
  log_q <- 0
  
  return(list(log_q = log_q, d_log_q = d_log_q))
}

source("code/temper_rjzigzag_hess.R")
## The code for zigzag_hess will use the hess matrix to construct bounds for thinning
# hess must contain 1 rows for likelihood with temp = 1

mu_prior <- rep(4, p)
sigmas <- rep(.5, p)
true_prob_inc <- rep(.5, p)
(hess <- matrix(rep(0,p),  nrow = 1))

set.seed(1)
theta_0 <- sample(c(1,-1), size = p, replace = T)
x_0 <- rnorm(p, mean = mu_prior, sd = sqrt(sigmas[1]))

print(x_0)
print(theta_0)
true_marginal_means <- mu_prior*true_prob_inc

## Standard method ---
set.seed(1);z_standard <- zigzag_temp(max_events = 1e4, 
                                      x0 = c(x_0,1), 
                                      theta0 = c(theta_0,0),
                                      prior_prob_inc = true_prob_inc,
                                      mus = mu_prior, sigmas = sigmas,
                                      alpha = 1, tau_max = 1, poly_order = 3,
                                      poly_coef = c(0,0,0))


s_standard <- gen_samples(z_standard$positions, z_standard$times,
                          z_standard$thetas, nsample = 1e4)

ppi_est <- rowMeans(abs(s_standard$samples)<1e-14)[1:p]
mean_est <- path_marginal_mean(z_standard)

## Tempered method 
set.seed(1);z_temp <- zigzag_temp(max_events = 1e4, 
                                  x0 = c(x_0,1), 
                                  theta0 = c(theta_0,0),
                                  prior_prob_inc = rep(.5, p),
                                  mus = mu_prior, sigmas = sigmas,
                                  alpha = .2,
                                  tau_max = 1, poly_order = 3,
                                  poly_coef = c(0,0,0))

plot_pdmp(z_temp, coords = c(1,2,p+1), pch = '.', inds = 1:1e4)
s_temp <- gen_samples(z_temp$positions, z_temp$times, z_temp$thetas, nsample = 1e4)
inds_1 <- which(s_temp$samples[p+1,] == 1)
ppi_est <- rowMeans(s_temp$samples[, inds_1] == 0)[1:p]
mean_est <- path_marginal_mean(z_temp)

library(ggplot2)
library(patchwork)
library(ccpdmp)
library(latex2exp)
z <- z_temp
df <- data.frame(x = z$positions[1,], y = z$positions[2,], temperature =  z$positions[p+1,])
pl <- ggplot(df, aes(x,y)) + geom_path(aes(colour=temperature), size = 1.2)+ theme_bw()
paths <- pl + scale_colour_gradient(low = "red", high = "black", na.value = NA)+ plot_annotation(
  title = "ct Zig-Zag"
)+ guides(colour=guide_legend(title=TeX(r"($\beta$)"), keyheight = 2))

paths

exact_samples <- matrix(rnorm(2*1e3, mean = mu_prior, sd = sqrt(sigmas[1])), nrow = 2)
exact_samples <- exact_samples*matrix(sample(c(0,1), size = 2*1e3, replace = T), nrow = 2)

df <- data.frame(x = exact_samples[1,], y = exact_samples[2,], temperature =  rep(1,1e3))
pl <- ggplot(df, aes(x,y)) + geom_point( size = 1.2)+ theme_bw() + plot_annotation(
  title = "Target"
)
samples <- pl
samples


z <- z_standard
df <- data.frame(x = z$positions[1,], y = z$positions[2,], temperature =  z$positions[p+1,])
pl <- ggplot(df, aes(x,y)) + geom_path(size = 1.2)+ theme_bw()
paths_stand <- pl+ plot_annotation(
  title = "Zig-Zag"
)
paths_stand

pl <- paths | ( samples/ paths_stand)
pl+ plot_annotation(tag_levels = "I") + plot_layout(guides = "collect")

