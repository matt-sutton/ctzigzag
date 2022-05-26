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

##
M <- 5; Reps <- 20
m_vals <- seq(from = 0, to = 4, length.out = M)
(hess <- matrix(rep(0,p),  nrow = 1))

## ZZ Mean, CT-ZZ Mean, ZZ PPI, CT-ZZ PPI
rmse <- array(0, dim = c(M,4,Reps))
dimnames(rmse) <- list(paste("m",round(m_vals,2)),
                       c("ZZ Mean", "ZZ PPI", "CT-ZZ Mean", "CT-ZZ PPI"), 1:Reps)

abs_err <- array(0, dim = c(M,4,Reps))
dimnames(abs_err) <- list(paste("m",round(m_vals,2)),
                       c("ZZ Mean", "ZZ PPI", "CT-ZZ Mean", "CT-ZZ PPI"), 1:Reps)
rmse1 <- array(0, dim = c(M,4,Reps))
dimnames(rmse1) <- list(paste("m",round(m_vals,2)),
                          c("ZZ Mean", "ZZ PPI", "CT-ZZ Mean", "CT-ZZ PPI"), 1:Reps)


est <- array(0, dim = c(M,4,2,Reps))
dimnames(est) <- list(paste("m",round(m_vals,2)),
                       c("ZZ Mean", "ZZ PPI", "CT-ZZ Mean", "CT-ZZ PPI"), 1:2, 1:Reps)

for(m in 1:length(m_vals)){
  for(i in 1:Reps){

    mu_prior <- rep(m_vals[m], p)
    sigmas <- rep(.5, p)
    true_prob_inc <- rep(.5, p)

    set.seed(i)
    theta_0 <- sample(c(1,-1), size = p, replace = T)
    x_0 <- rnorm(p, mean = m_vals[m], sd = sqrt(sigmas[1]))

    print(x_0)
    print(theta_0)
    true_marginal_means <- mu_prior*true_prob_inc

    ## Standard method ---
    set.seed(i);z_standard <- zigzag_temp(max_events = 1e4, 
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

    est[m, 1, , i] <- mean_est
    est[m, 2, , i] <- ppi_est

    rmse[m, 1, i] <-  mean((mean_est - true_marginal_means)^2)
    rmse[m, 2, i] <-  mean((ppi_est - true_prob_inc)^2)

    rmse1[m, 1, i] <-  sqrt(mean((mean_est[1] - true_marginal_means[1])^2))
    rmse1[m, 2, i] <-  sqrt(mean((ppi_est[1] - true_prob_inc[1])^2))

    abs_err[m,1,i]<- mean(abs(mean_est - true_marginal_means))
    abs_err[m,2,i]<- mean(abs(ppi_est - true_prob_inc))

    ## Tempered method ---
    set.seed(i);z_temp <- zigzag_temp(max_events = 1e4, x0 = c(x_0,1), theta0 = c(theta_0,0),
                                      prior_prob_inc = rep(.5, p),
                                      mus = mu_prior, sigmas = sigmas,
                                      alpha = .2,
                                      tau_max = 1, poly_order = 3,
                                      poly_coef = c(0,0,0))

    # plot_pdmp(z_temp, coords = c(1,2,p+1), pch = '.', inds = 1:1e4)
    s_temp <- gen_samples(z_temp$positions, z_temp$times, z_temp$thetas, nsample = 1e4)
    inds_1 <- which(s_temp$samples[p+1,] == 1)
    ppi_est <- rowMeans(s_temp$samples[, inds_1] == 0)[1:p]
    mean_est <- path_marginal_mean(z_temp)

    rmse[m, 3, i] <-  mean((true_marginal_means - mean_est)^2)
    rmse[m, 4, i] <-  mean((ppi_est - true_prob_inc)^2)

    rmse1[m, 3, i] <-  sqrt(mean((mean_est[1] - true_marginal_means[1])^2))
    rmse1[m, 4, i] <-  sqrt(mean((ppi_est[1] - true_prob_inc[1])^2))

    abs_err[m,1,i]<- mean(abs(mean_est - true_marginal_means))
    abs_err[m,2,i]<- mean(abs(ppi_est - true_prob_inc))

    est[m, 3, , i] <- mean_est
    est[m, 4, , i] <- ppi_est
  }
  gc()
}

avg_rmse <- apply(rmse1, c(1,2), mean)
stargazer::stargazer(t(avg_rmse))

avg_rmse <- apply(rmse, c(1,2), mean)
plot(m_vals, y = avg_rmse[,1], type = 'l')
lines(m_vals, avg_rmse[,3], col = 2)

plot(m_vals, y = avg_rmse[,2], type = 'l')
lines(m_vals, avg_rmse[,4], col = 2)

stargazer::stargazer(t(avg_rmse), digits = 1, digits.extra = 1)


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

sza <- gen_samples(z$positions,z$times, nsample = 1e4 , burn = 10)
ind1 <- which(sza$samples[p+1,] == 1)
df <- data.frame(x = sza$samples[1,ind1], y = sza$samples[2,ind1], temperature =  sza$samples[p+1,ind1])
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

