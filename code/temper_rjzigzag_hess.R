library(ccpdmp)
source("code/utils.R")

## Code is set up for general tempering:
## log q(x,b) = bl(x) + sum(phi(x_i; bm, sigma^2) + d(x_i)) where l(x) is the likleihood

## Function returning the gradient for - log q(x, t)
d_nlogq <- function(x, t, theta, d_poly_coef, mus, sigmas){
  q <- target(x)
  which_nz_x <- which(abs(theta[1:length(x)])> 1e-16)

  q_0_log <- -sum(mus[which_nz_x]*(x[which_nz_x] - mus[which_nz_x]*t)/(sigmas[which_nz_x]^2))

  x_grads <- - t*q$d_log_q  + (x - mus*t)/sigmas^2

  t_grad <- pracma::polyval(d_poly_coef, t) - q$log_q + q_0_log

  grads <- c(x_grads, t_grad)
  return(grads)
}

## Return the rates for each dimension evaluated at times tau_grid ahead
return_rates_hess <- function(x, theta, tau_grid, d_poly_coef, mus, sigmas){
  l_tau <- length(tau_grid)
  nx <- length(x)
  rates_eval <- matrix(0, nx, l_tau)

  ## Find the a, b and c terms for the rates on x_j
  q <- target(x[-nx]);

  a_q <- -q$d_log_q*theta[-nx]; a_0 <- theta[-nx]*(x[-nx] - mus*x[nx])/sigmas^2
  b_q <- hess*abs(theta[-nx]); b_0 <- theta[-nx]*(theta[-nx] - mus*theta[nx])/sigmas^2

  t <- x[nx]
  a <- t*a_q + a_0;
  b <- theta[nx]*a_q + t*b_q + b_0
  c <- theta[nx]*b_q
  poly_coef_x <- rbind(c,b,a)

  ## Find the a, b and c terms for the rates on temp
  which_nz_x <- which(abs(theta[-nx])> 1e-16)
  n_nx <- length(which_nz_x)

  q_0_log <- -sum(mus[which_nz_x]*(x[which_nz_x] - mus[which_nz_x]*t)/(sigmas[which_nz_x]^2))

  a_q <- -theta[nx]*q$log_q; a_p <- theta[nx]*q_0_log
  b_q <- -theta[nx]*sum(theta[-nx]*q$d_log_q);  b_p <- -theta[nx]*sum((theta[which_nz_x]*mus[which_nz_x]- mus[which_nz_x]^2*theta[nx])/(sigmas[which_nz_x]^2))
  c_q <- 0.5*sum(hess*(theta[-nx])^2)*abs(theta[nx]); c_p <- 0
  a <- a_q + a_p; b <- b_q + b_p; c <- c_q + c_p
  poly_coef_t <- c(c,b,a)

  ## Evaluate the upper-bounding rate at the time points
  for( i in 1:(nx-1) ){
    rates_eval[i,] <- pracma::polyval(poly_coef_x[,i], tau_grid)
  }
  rates_eval[nx,] <- pracma::polyval(poly_coef_t, tau_grid) +
    theta[nx]*pracma::polyval(d_poly_coef, t + tau_grid*theta[nx])

  return(rates_eval)
}
## Return the rates to reintroduce
return_rates_reintro <- function(x, theta, tau_grid, prior_inc, mus, sigmas){
  
  l_tau <- length(tau_grid)
  nx <- length(x)
  which_sim <- which(abs(theta[-nx]) < 1e-15)
  nsim <- length(which_sim)
  rates_eval <- matrix(0, nsim, l_tau)

  w <- prior_inc[which_sim]

  ## Rate to reintro is convex so just return evaluation times
  for( i in 1:length(tau_grid) ){
    temp_s <- x[nx] + tau_grid[i]*theta[nx]
    ratios <- w/(1-w)
    rates_eval[,i] <- ratios*( 1/sqrt(2*pi*sigmas[which_sim]^2) )
  }
  return(rates_eval)
}


## calculate time until vector x hits bnd
time_to_bnd <- function(x, theta, bnd){
  times_bnd <- rep(Inf, length(x))

  towards_bnd <- which(sign(theta*(x-bnd)) == -1 )
  times_bnd[towards_bnd] <- (bnd-x[towards_bnd])/theta[towards_bnd]

  return(times_bnd)
}


zigzag_temp <- function(max_events, x0, theta0 = c(rep(1, length(x0)-1), 0.1),
                        alpha = 0.5, tau_max = 1, return_rates = return_rates_hess,
                        prior_prob_inc = rep(.5, length(x0)-1),
                        mus = rep(0, length(x0)-1),sigmas = rep(2, length(x0)-1),
                        poly_order = 2, echo = FALSE, poly_coef = c(0,0)){
  kappa_m <- function(t){
    return( exp( - pracma::polyval(poly_coef, t) ) )
  }

  poly_order_temp <- length(poly_coef)
  
  # derivative terms i.e. sum(i*a_it^{i-1})
  d_poly_coef <- c(poly_order_temp:1 - 1)*poly_coef
  d_poly_coef <- d_poly_coef[-poly_order_temp]

  ## Init
  t = 0; eps = 1e-10; nits <- 0
  x = x0; theta = theta0; nvel <- length(x)

  thetas <- positions <- matrix(0, nrow = nvel, ncol = max_events);
  times = rep(0,max_events); thetas[,1] <- theta; positions[,1] <- x;

  num_evts = 1
  event = FALSE

  ## Rates to jump from dirac on temperature
  if(abs(theta0[nvel]) < 1e-10){
    rate_jump_temperature <- (1-alpha)/(2*alpha)
  } else {
    rate_jump_temperature <- ((1-alpha))/(2*alpha)*abs(theta0[nvel])
  }

  ## Find time to hitting boundaries
  x_to_bnd <- time_to_bnd(x, theta, 0)
  t_to_bnd <- min(time_to_bnd(x[nvel], theta[nvel], 0), time_to_bnd(x[nvel], theta[nvel], 1))

  # Simulate times
  taus = rep(Inf, nvel);  u_s = rexp(nvel);  f_s = rep(Inf, nvel)
  tau_grid <- seq(0, to = min(c(tau_max,x_to_bnd, t_to_bnd)),  length.out = poly_order + 1)

  rates <- return_rates(x, theta, tau_grid, d_poly_coef, mus, sigmas)
  for( i in 1:nvel){
    tus = sim_rate_poly(eval_times = tau_grid, eval_rates = rates[i,], poly_order)
    taus[i] = tus$t
    u_s[i] = tus$u
    f_s[i] = tus$f_evall
  }

  # Simulate time for reintroducing variables
  sampling_pnt_mass <- abs(theta) < eps

  ## Reintroduce temperature variable
  if(sampling_pnt_mass[nvel]){
    taus[nvel] <- rexp(1)/rate_jump_temperature
  }

  ## Reintroduce the variables x
  sampling_pnt_mass_x <- sampling_pnt_mass[-nvel]
  if(any(sampling_pnt_mass_x)){
    reintro_x <- which(sampling_pnt_mass_x)
    rates <-  return_rates_reintro(x, theta, tau_grid, prior_prob_inc, mus, sigmas)
    for( i in 1:length(reintro_x)){
      tus = cc_sim(tau_grid, eval_rates = rbind(rates[i,],0,0))
      j <- reintro_x[i]
      taus[j] = tus$t
      u_s[j] = tus$u
      f_s[j] = tus$f_evall
    }
  }

  while(num_evts < max_events){

    mini_x <- which.min(taus)
    tau <- taus[mini_x]

    x = x + tau*theta
    t = t + tau

    ## If temp or x hits boundary it is an event
    bndry_x <- abs(x[-nvel] - rep(0,nvel-1)) < 1e-15
    bndry_temp <- abs(x[nvel] - c(0,1)) < 1e-15

    if(any(taus[which(bndry_x)] == tau & abs(theta[which(bndry_x)]) > 0)){
      mini_x <- which(abs(tau - x_to_bnd) < 1e-15)
    }

    update_temp <- (any(bndry_temp) & abs(tau - taus[nvel]) < eps) ||
      (abs(tau - taus[nvel]) < eps & u_s[nvel] < eps)

    if(update_temp){
      b_hit_1 <- abs(x[nvel]-1) < eps
      b_hit_0 <- abs(x[nvel]) < eps
      
      ## If hit 0 flip 
      if(b_hit_0){
        theta[nvel] = 1
        event = TRUE
      }
      ## If hit spike release or stick
      if( b_hit_1 ){
        if(sampling_pnt_mass[nvel]){
          theta[nvel] = -1
          sampling_pnt_mass[nvel] <- FALSE
        } else {
          theta[nvel] = 0
          sampling_pnt_mass[nvel] <- TRUE
        }
        event = TRUE
      }
      ## If regular event
      if(!b_hit_1 & !b_hit_0){
        ## If regular event
        if(u_s[nvel] < 1e-10){
          grad <- d_nlogq(x[-nvel], x[nvel], theta, d_poly_coef, mus, sigmas)
          rate <- grad[nvel]*theta[nvel]

          acc_prb <- rate/f_s[nvel]
          if(acc_prb > 1.0001){
            print(paste("Invalid thinning on inverse temp, thinning prob:",acc_prb))
          }
          if(runif(1) <= acc_prb){
            theta[nvel] = -theta[nvel]
            reintro_theta = theta[nvel]
            event = TRUE
          }
        }
      }

    } else {

      ## Check for variable selection
      if(abs(x[mini_x]) < 1e-15){
        if(sampling_pnt_mass[mini_x]){

          if(u_s[mini_x] < 1e-10){

            # Calculate rate to reintroduce
            or <- prior_prob_inc[mini_x]/(1-prior_prob_inc[mini_x])
            rate <- sqrt(1/(2*pi*sigmas[mini_x]^2))*or*
              exp( -mus[mini_x]^2*x[nvel]^2/(2*sigmas[mini_x]^2))

            acc_prb <- rate/f_s[mini_x]
            if(acc_prb > 1.001){
              print(paste("Invalid thinning on x, thinning prob:",acc_prb))
            }
            if(runif(1) <= acc_prb){
              ## Re-introduce velocity
              theta[mini_x] = sample(c(1,-1), 1)
              sampling_pnt_mass[mini_x] <- FALSE
              event = TRUE
            }
          }

        } else {
          if(prior_prob_inc[mini_x] < 1){
            theta[mini_x] = 0
            sampling_pnt_mass[mini_x] <- TRUE
            event = TRUE
          }
        }

      } else {
        ## Regular event
        if(u_s[mini_x] < 1e-10){
          grad <- d_nlogq(x[-nvel], x[nvel], theta, d_poly_coef, mus, sigmas)
          rate <- grad[mini_x]*theta[mini_x]

          acc_prb <- rate/f_s[mini_x]
          if(acc_prb > 1.001){
            print(paste("Invalid thinning on x, thinning prob:",acc_prb))
          }
          if(runif(1) <= acc_prb){
            theta[mini_x] = -theta[mini_x]
            event = TRUE
          }
        }
      }
    }

    if(event){
      # Store event info
      num_evts = num_evts + 1
      thetas[,num_evts] <- theta; positions[,num_evts] <- x;
      times[num_evts] = t

      event = FALSE
      if(echo & (num_evts %% 100 == 0)){
        print(num_evts)
      }

      # Simulate times
      x_to_bnd <- time_to_bnd(x, theta, 0)
      t_to_bnd <- min(time_to_bnd(x[nvel], theta[nvel], 0), time_to_bnd(x[nvel], theta[nvel], 1))
      tau_grid <- seq(0, to = min(c(tau_max,x_to_bnd, t_to_bnd)),  length.out = poly_order + 1)

      rates <- return_rates(x, theta, tau_grid, d_poly_coef, mus, sigmas)
      for( i in 1:nvel){
        tus = sim_rate_poly(eval_times = tau_grid, eval_rates = rates[i,], poly_order)
        taus[i] = tus$t
        u_s[i] = tus$u
        f_s[i] = tus$f_evall
      }

      ## Reintroduce temperature variable
      if(sampling_pnt_mass[nvel]){
        taus[nvel] <- rexp(1)/rate_jump_temperature
      }

      ## Reintroduce the variables x
      sampling_pnt_mass_x <- sampling_pnt_mass[-nvel]
      if(any(sampling_pnt_mass_x)){
        reintro_x <- which(sampling_pnt_mass_x)
        rates <-  return_rates_reintro(x, theta, tau_grid, prior_prob_inc, mus, sigmas)
        for( i in 1:length(reintro_x)){
          tus = cc_sim(tau_grid, eval_rates = rbind(rates[i,],0,0))
          j <- reintro_x[i]
          taus[j] = tus$t
          u_s[j] = tus$u
          f_s[j] = tus$f_evall
        }
      }

    } else {
      # If there was no event

      # Re-simulate times for all taus less than zero:
      x_to_bnd <- time_to_bnd(x, theta, 0)
      t_to_bnd <- min(time_to_bnd(x[nvel], theta[nvel], 0), time_to_bnd(x[nvel], theta[nvel], 1))
      tau_grid <- seq(0, to = min(c(tau_max,x_to_bnd, t_to_bnd)),  length.out = poly_order + 1)

      # Adjust simulated times
      taus <- taus-tau

      update_rates <- which(taus <= 0)

      rates <- return_rates(x, theta, tau_grid, d_poly_coef, mus, sigmas)
      for( j in update_rates){
        tus = sim_rate_poly(eval_times = tau_grid, eval_rates = rates[j,], poly_order)
        taus[j] = tus$t
        u_s[j] = tus$u
        f_s[j] = tus$f_evall
      }

      ## Reintroduce temperature variable
      if(sampling_pnt_mass[nvel]){
        taus[nvel] <- rexp(1)/rate_jump_temperature
      }

      ## Reintroduce the variables x
      sampling_pnt_mass_x <- sampling_pnt_mass[-nvel]
      if(any(sampling_pnt_mass_x)){
        reintro_x <- which(sampling_pnt_mass_x)
        rates <-  return_rates_reintro(x, theta, tau_grid, prior_prob_inc, mus, sigmas)
        for( i in 1:length(reintro_x)){
          tus = cc_sim(tau_grid, eval_rates = rbind(rates[i,],0,0))
          j <- reintro_x[i]
          taus[j] = tus$t
          u_s[j] = tus$u
          f_s[j] = tus$f_evall
        }
      }

    }
    nits = nits +1
  }
  return (list(positions=positions,thetas=thetas,times=times, nits = nits, 
               poly_coef = poly_coef, alpha = alpha))
}
