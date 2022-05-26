load("GMM/gmm_data.Rdata")
library(rstan)
## ct HMC
stan_fit_eval <- stan("GMM/ct_hmc_eval.stan", data = dat, iter = 1, chains = 1)

## Setup ct-ZZ
target <- function(x){
  stan_ev <- grad_log_prob(stan_fit_eval, c(x,1))
  d_log_q <- as.numeric(stan_ev)[1:2]
  log_q <- attr(stan_ev, "log_prob")

  return(list(log_q = log_q, d_log_q = d_log_q))
}
temper <- function(x){
  stan_ev <- grad_log_prob(stan_fit_eval, c(x,0))
  d_log_q <- as.numeric(stan_ev)[1:2]
  log_q <- attr(stan_ev, "log_prob")

  return(list(log_q = log_q, d_log_q = d_log_q))
}

## only take first 20,000 iters
np <- 3e2
x_v <- seq(from = -1, to = 12, length.out = np)
y_v <- seq(from = -1, to = 12, length.out = np)
xy <- expand.grid(x=x_v, y=y_v)
f <- function(x){exp(target(x)$log_q)}
z <- matrix(apply(as.matrix(xy), 1, f), length(x_v), length(y_v))
par(mfrow = c(1,2), mar = c(2,2,1,2))
image(x_v, y_v, z, las=1, xlab = "x1", ylab = "x2")
load("~/Gits/ctzigzag/GMM/fix_alpha/gmm_alpha_1_iter_1.Rdata")
lines(t(zigzag_fit_ada$positions)[1:30000,1:2],
      col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
#legend('topleft', legend = c('Zig Zag'), col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5), lwd = 1)

colnames(z) = y_v
rownames(z) = x_v

library(tidyverse)
dat = as.data.frame(z) %>%
  rownames_to_column(var="x_v") %>%
  gather(y_v, value, -x_v) %>%
  mutate(y_v=as.numeric(y_v),
         x_v=as.numeric(x_v),
         value_range = cut(value, 8))

plot_1 <-
  ggplot(dat, aes(x_v, y_v, fill=value_range)) +
  geom_raster(show.legend = FALSE) +
  scale_fill_manual(values=colorRampPalette(c("white","black"))(10)) +
  theme_bw()+ theme(axis.title.x=element_blank(), axis.title.y=element_blank())

load("~/Gits/ctzigzag/GMM/fix_alpha/gmm_alpha_1_iter_1.Rdata")
df_z <- data.frame(x = t(zigzag_fit_ada$positions)[1:30000,1],
                   y = t(zigzag_fit_ada$positions)[1:30000,2],
                   value_range = rep(dat$value_range[1],30000) )

plot_1z <- plot_1 + geom_path(data=df_z,aes(x=x,y=y, fill=value_range, colour = 1), show.legend = F)+
  scale_colour_gradient(low = "darkred", high = "darkred", na.value = NA)

plot_1 + geom_path(data=df_z,aes(x=x,y=y, fill=value_range, alpha=0.2))

load("~/Gits/ctzigzag/GMM/fix_alpha/gmm_alpha_0.7_iter_1.Rdata")
df_z <- data.frame(x = t(zigzag_fit_ada$positions)[1:30000,1],
                   y = t(zigzag_fit_ada$positions)[1:30000,2],
                   value_range = rep(dat$value_range[1],30000),
                   col = zigzag_fit_ada$positions[3,])

plot_2z <- plot_1 + geom_path(data=df_z,aes(x=x,y=y, fill=value_range, colour=col)) +
  scale_colour_gradient(low = "pink", high = "darkred", na.value = NA)
plot_2z

library(patchwork)
library(latex2exp)
(plot_1z | plot_2z)+ guides(colour=guide_legend(title=TeX(r"($\beta$)"), keyheight = 2))

df <- data.frame(x=x_v, y=y_v, inv_temp = z)
ggplot(df, aes(x=x,y=y))
library(ggmap)
ggimage(z)

image(x_v, y_v, z, las=1, xlab = "x1", ylab = "x2")
load("~/Gits/ctzigzag/GMM/fix_alpha/gmm_alpha_0.7_iter_1.Rdata")
lines(t(zigzag_fit_ada$positions)[1:30000,1:2],
      col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
# legend('topleft', legend = c('Tempered Zig Zag'), col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5), lwd = 1)

par(mfrow = c(2,1))
hist(zigzag_samples$samples[3,zigzag_samples$samples[3,]<1-1e-8], probability = T,
     main = "Marginal inverse teperature")
plot(zigzag_fit_ada$times[5e3:(5e3+200)],
     zigzag_fit_ada$positions[3,5e3:(5e3+200)],ylim = c(0,1), type = 'l')


c_mm <- function(pdmp, burnin = 1){
  ## Calc Marginal
  d <- nrow(pdmp$positions)-1
  in_model <- (pdmp$positions[d+1, ] == 1 & pdmp$thetas[d+1, ] == 0)
  times <- pdmp$times
  positions <- pdmp$positions
  thetas <- pdmp$thetas

  maxIter <- length(times)
  marg_mean <- rep(0, d)

  for( mi in 1:d){
    beta_mean <- 0
    total_time <- 0
    for(i in (burnin):(maxIter-1)){
      if(in_model[i]){
        tauv <- (times[i+1] - times[i])
        total_time <- total_time + tauv
        beta_mean = beta_mean + (tauv*positions[mi,i] + thetas[mi,i]*tauv^2/2)
        marg_mean[mi] = beta_mean/total_time
      }

    }
  }
  return(marg_mean)
}
c_mm2 <- function(pdmp, burnin = 1){
  ## Calc Marginal 2 mom
  d <- nrow(pdmp$positions)-1
  in_model <- (pdmp$positions[d+1, ] == 1 & pdmp$thetas[d+1, ] == 0)
  times <- pdmp$times
  positions <- pdmp$positions
  thetas <- pdmp$thetas

  maxIter <- length(times)
  marg_mean <- rep(0, d)

  for( mi in 1:d){
    beta_mean <- 0
    total_time <- 0
    for(i in (burnin):(maxIter-1)){
      if(in_model[i]){
        tv <- (times[i+1] - times[i])
        total_time <- total_time + tv
        x_t <- positions[mi,i]; v_t <- thetas[mi,i]
        beta_mean = beta_mean +
          (tv*x_t^2 + tv^2*v_t*x_t + tv^3*v_t^2/3)
        marg_mean[mi] = beta_mean/total_time
      }

    }
  }
  return(marg_mean)
}

load("GMM/gmm_data.Rdata")
exact_mom_1 <- rowMeans(mu_mat)
exact_mom_2 <- rowMeans(mu_mat^2 + dat$sigma[1])
num_reps = 20
alphas <- c(.1, .2, .3, .5, .7, .8)[6:1]
est_ZZ_all_2 <-est_ZZ_all_1 <- matrix(0, nrow = 8, ncol = num_reps)

for(j in 1:6){
  alpha <- alphas[j]
  est1_ZZ <- est2_ZZ <- est1_ZZ1 <-
    est2_ZZ1 <- matrix(0, nrow = 2, ncol = num_reps)

  n_iter <- rep(0, num_reps)

  stick_to_one <- rep(0,num_reps)
  comp_eff <- matrix(0, num_reps, 2)
i=1
for( i in 1:num_reps ){
  fe <- file.exists(paste0("~/Gits/ctzigzag_old/GMM/fix_alpha/gmm_alpha_",alpha,"_iter_",i,".Rdata"))

  if(fe){
    load(paste0("~/Gits/ctzigzag_old/GMM/fix_alpha/gmm_alpha_",alpha,"_iter_",i,".Rdata"))
    nr <- 2

    stick_to_one[i] = 1- rjpdmp::model_probabilities(zigzag_fit_ada$times, thetas = zigzag_fit_ada$thetas, marginals = 3)$marginal_prob
    est1_ZZ[,i] <- c_mm(zigzag_fit_ada); est2_ZZ[,i] <- c_mm2(zigzag_fit_ada)
    est_ZZ_all_1[j+1,i] <- est1_ZZ[1,i]
    est_ZZ_all_2[j+1,i] <- est2_ZZ[1,i]

    comp_eff[i,1] <- length(zigzag_fit_ada$times)/zigzag_fit_ada$nits

    load(paste0("~/Gits/ctzigzag_old/GMM/fix_alpha/gmm_alpha_",1,"_iter_",i,".Rdata"))
    est1_ZZ1[,i] <- c_mm(zigzag_fit_ada); est2_ZZ1[,i] <- c_mm2(zigzag_fit_ada)
    est_ZZ_all_1[1,i] <- est1_ZZ1[1,i]
    est_ZZ_all_2[1,i] <- est2_ZZ1[1,i]
    comp_eff[i,2] <- length(zigzag_fit_ada$times)/zigzag_fit_ada$nits

  } else{
    ind_rem <- c(ind_rem, i)
  }
}
prob_1 <- c(1,mean(stick_to_one))
comp <- colMeans(comp_eff)[2:1]
rmse_1 <- rbind(sqrt(rowMeans((est1_ZZ1 - exact_mom_1)^2)),
            sqrt(rowMeans((est1_ZZ - exact_mom_1)^2)))
rmse_2 <- rbind(sqrt(rowMeans((est2_ZZ1 - exact_mom_2)^2)),
                sqrt(rowMeans((est2_ZZ - exact_mom_2)^2)))
method <- c("Zig-Zag", "ctZig-Zag")
if(j ==1){
  row_add <- data.frame(method = method, alpha = c(1,alpha),
                        prob_1 = prob_1, Ex = rmse_1,
                        Ex2 = rmse_2, computation = comp)
} else {
  row_add <- rbind(row_add,
                   data.frame(method = method, alpha = c(1,alpha),
                        prob_1 = prob_1, Ex = rmse_1,
                        Ex2 = rmse_2, computation = comp)[-1,])
}
}
row_add
dim(zigzag_samples$samples)

## Add in the only tempering version (alpha = 0)....
alpha <- 0
est1_ZZ <- est2_ZZ <- est1_ZZ1 <-
  est2_ZZ1 <- matrix(0, nrow = 2, ncol = num_reps)

n_iter <- rep(0, num_reps)

stick_to_one <- rep(0,num_reps)
comp_eff <- matrix(0, num_reps, 2)
i=1
for( i in 1:num_reps ){
  fe <- file.exists(paste0("~/Gits/ctzigzag_old/GMM/fix_alpha/gmm_alpha_",alpha,"_iter_",i,".Rdata"))

  if(fe){
    load(paste0("~/Gits/ctzigzag_old/GMM/fix_alpha/gmm_alpha_",alpha,"_iter_",i,".Rdata"))
    nr <- 2

    stick_to_one[i] = 1- rjpdmp::model_probabilities(zigzag_fit_ada$times, thetas = zigzag_fit_ada$thetas, marginals = 3)$marginal_prob

    # est1_ZZ[,i] <- colSums(t(zigzag_samples$samples[1:nr,])*w_1_zz_norm)
    # est2_ZZ[,i] <- diag(crossprod(t(zigzag_samples$samples[1:nr,])*w_1_zz_norm, t(zigzag_samples$samples[1:nr,])))
    est1_ZZ[,i] <- colSums(t(zigzag_samplesl$samples[1:nr,])*w_1_zz_norml)
    est2_ZZ[,i] <- diag(crossprod(t(zigzag_samplesl$samples[1:nr,])*w_1_zz_norml, t(zigzag_samplesl$samples[1:nr,])))

    est_ZZ_all_1[8,i] <- est1_ZZ[1,i]
    est_ZZ_all_2[8,i] <- est2_ZZ[1,i]
    comp_eff[i,1] <- length(zigzag_fit_ada$times)/zigzag_fit_ada$nits

    load(paste0("~/Gits/ctzigzag_old/GMM/fix_alpha/gmm_alpha_",1,"_iter_",i,".Rdata"))
    est1_ZZ1[,i] <- c_mm(zigzag_fit_ada); est2_ZZ1[,i] <- c_mm2(zigzag_fit_ada)

    comp_eff[i,2] <- length(zigzag_fit_ada$times)/zigzag_fit_ada$nits

  } else{
    ind_rem <- c(ind_rem, i)
  }
}
par(mfrow = c(2,1), mar=c(2.5,2,2,2))
est <- t(est_ZZ_all_1)[,c(1,2,6,8)]
colnames(est) <- paste0("alpha = ", c(1,.8,.2,0))
boxplot(est,ylim =c(2.7,9.5))
title(TeX(r"(Recovery of $E\[X_1\]$)"))
abline(h = exact_mom_1[1])

est <- t(est_ZZ_all_2)[,c(1,2,6,8)]
colnames(est) <- paste0("alpha = ", c(1,.8,.2,0))
boxplot(est, ylim=c(5,55))
abline(h = exact_mom_2[1])
title(TeX(r"(Recovery of $E\[X_1^2\]$)"))

prob_1 <- c(1,mean(stick_to_one))
comp <- colMeans(comp_eff)[2:1]
rmse_1 <- rbind(sqrt(rowMeans((est1_ZZ1 - exact_mom_1)^2)),
                sqrt(rowMeans((est1_ZZ - exact_mom_1)^2)))
rmse_2 <- rbind(sqrt(rowMeans((est2_ZZ1 - exact_mom_2)^2)),
                sqrt(rowMeans((est2_ZZ - exact_mom_2)^2)))
method <- c("Zig-Zag", "ctZig-Zag")

row_add <- rbind(row_add,
                 data.frame(method = method, alpha = c(1,alpha),
                            prob_1 = prob_1, Ex = rmse_1,
                            Ex2 = rmse_2, computation = comp)[-1,])
stargazer::stargazer(row_add, summary = F)


boxplot(t(est_ZZ_all_1))



