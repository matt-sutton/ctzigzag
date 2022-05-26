## Base R files for processing Zig-Zag sampler

## Generate discrete samples from zigzag skeleton
gen_samples <- function(positions, times, theta,
                        nsample = 10^3, burn = 1){
  
  if(is.null(dim(positions))) positions <- matrix(positions, nrow = 1)
  
  positions <- positions[,burn:length(times), drop = F]
  theta <- theta[,burn:length(times), drop = F]
  times <- times[burn:length(times)] - times[burn]
  nsteps <- length(times)
  
  Tmax <- times[nsteps]
  dt <- Tmax/(nsample+2)
  t = dt
  t0 = times[1]
  x0 = positions[,1]
  samples <- matrix(0, nrow = length(x0), ncol = nsample)
  sample_times <- rep(0, nsample)
  n <- 0
  
  for(i in 2:nsteps){
    x1 = positions[,i]
    t1 = times[i]
    theta0 = theta[,i-1]
    while(t + dt < t1 && n < nsample){
      n <- n+1
      t <- t + dt
      x_s <- x0 + (t-t0)*theta0
      samples[,n] <- x_s
      sample_times[n] <- t
    }
    x0 = x1; t0 = t1
  }
  return(list(samples = samples, sample_times =sample_times))
}

## Calculate the marginal mean based on the pdmp trajectory 
path_marginal_mean <- function(pdmp, burnin = 1){
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
path_marginal_moment2 <- function(pdmp, burnin = 1){
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

## object for evaluating path quadrature estimate of logz(b)
phi <- function(xt){
  x <- xt[-length(xt)]; t <- xt[length(xt)]
  
  q <- target(x); q_0 <- temper(x)
  return(q$log_q - q_0$log_q)
}
## object for evaluating path quadrature
path_quad <- function(t, u, t_lower=0, t_upper=1){
  keep <- (t > t_lower) & (t < t_upper)
  if (sum(keep) > 0){
    t_uniq <- sort(unique(t[keep]))
    N_uniq <- length(t_uniq)
    u_bar <- rep(NA, N_uniq)
    for (i in 1:N_uniq){
      ok <- t == t_uniq[i]
      u_bar[i] <- mean(u[ok])
    }
    width <- c(t_uniq, t_upper) - c(t_lower, t_uniq)
    u_bar_avg <- (c(u_bar, u_bar[N_uniq]) + c(u_bar[1], u_bar))/ 2
    log_z <- c(0, cumsum(width*u_bar_avg))
    return(list(t=c(0, t_lower, t_uniq, t_upper, 1), log_z=c(log_z[1], log_z, log_z[length(log_z)])))
  }
  else
    return(list(t=c(0, 1), log_z=c(0, 0)))
}