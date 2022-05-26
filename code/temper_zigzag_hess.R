library(ccpdmp)
source("code/utils.R")

## Function returning the gradient for - log q(x, t)
d_nlogq <- function(x, t, d_poly_coef){
  q <- target(x)
  q_0 <- temper(x)

  x_grads <- - t*q$d_log_q - (1-t)*q_0$d_log_q
  t_grad <- pracma::polyval(d_poly_coef, t) - q$log_q + q_0$log_q

  grads <- c(x_grads, t_grad)
  return(grads)
}

## Return the rates for each dimension evaluated at times tau_grid ahead
return_rates_hess <- function(x, theta, tau_grid, d_poly_coef){
  l_tau <- length(tau_grid)
  nx <- length(x)
  rates_eval <- matrix(0, nx, l_tau)

  ## Find the a, b and c terms for the rates on x_j
  q <- target(x[-nx]);  q_0 <- temper(x[-nx])

  a_q <- -q$d_log_q*theta[-nx]; a_q_0 <- -q_0$d_log_q*theta[-nx]
  b_q <- hess[2,]*abs(theta[-nx]); b_q_0 <- hess[1,]*abs(theta[-nx])

  t <- x[nx]
  a <- t*a_q + (1-t)*a_q_0;
  b <- theta[nx]*(a_q - a_q_0) + b_q*t + (1-t)*b_q_0
  c <- theta[nx]*(b_q - b_q_0)
  poly_coef_x <- rbind(c,b,a)

  ## Find the a, b and c terms for the rates on temp
  a_q <- -theta[nx]*q$log_q; a_q_0 <- theta[nx]*q_0$log_q
  b_q <- -theta[nx]*sum(theta[-nx]*q$d_log_q); b_q_0 <- theta[nx]*sum(theta[-nx]*q_0$d_log_q)
  c_q <- 0.5*sum(hess[2,]*(theta[-nx])^2)*abs(theta[nx]); c_q_0 <- 0.5*sum(hess[1,]*(theta[-nx])^2)*abs(theta[nx])
  a <- a_q + a_q_0; b <- b_q + b_q_0; c <- c_q + c_q_0
  poly_coef_t <- c(c,b,a)

  ## Evaluate the upper-bounding rate at the time points
  for( i in 1:(nx-1) ){
    rates_eval[i,] <- pracma::polyval(poly_coef_x[,i], tau_grid)
  }
  rates_eval[nx,] <- pracma::polyval(poly_coef_t, tau_grid) +
    theta[nx]*pracma::polyval(d_poly_coef, t + tau_grid*theta[nx])
  return(rates_eval)
}

## calculate time until inverse temp = 0 or = 1
time_to_bndry <- function(x, theta){
  if(theta == 0) {
    return(Inf)
  } else if( theta > 0){
    # if traveling to 1
    return( (1-x)/theta )
  } else {
    # if traveling to 0
    return(  -x/theta )
  }
}

## Tempered Zig-Zag algorithm
zigzag_temp <- function(max_events, x0, theta0 = c(rep(1, length(x0)-1), 0.1),
                        alpha = 0.2, tau_max = 1, 
                        return_rates = return_rates_hess,
                        poly_order = 2, echo = FALSE, 
                        poly_coef = c(0,0), nits_max = Inf){
  
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

  ## Rates to jump from dirac
  if(abs(theta0[nvel]) < 1e-10){
    rate_jump_from_spike <- (1-alpha)/(2*alpha)
    temp_scale <- 1
  } else {
    rate_jump_from_spike <- (1-alpha)/(2*alpha)*abs(theta0[nvel])
    temp_scale <- abs(theta0[nvel])
  }

  # Simulate times
  taus = rep(Inf, nvel);  u_s = rexp(nvel);  f_s = rep(Inf, nvel)
  temp_to_bndry <- time_to_bndry(x[nvel], theta[nvel])
  tau_grid <- seq(0, to = min(tau_max,temp_to_bndry),  length.out = poly_order + 1)

  ## Evaluates the (upper-bounded) rate 
  rates <- return_rates(x, theta, tau_grid, d_poly_coef)
  
  ## Simulate from the upper-bounding rate function using ccpdmp package 
  for( i in 1:nvel){
    tus = sim_rate_poly(eval_times = tau_grid, eval_rates = rates[i,], poly_order)
    taus[i] = tus$t
    u_s[i] = tus$u
    f_s[i] = tus$f_evall
  }

  # Simulate time to reintroduce inverse temp
  sampling_spike <- if(abs(theta[nvel]) < eps) TRUE else FALSE ## If b=1
  
  if(sampling_spike){
    taus[nvel] <- rexp(1)/rate_jump_from_spike
  }

  while(num_evts < max_events & nits < nits_max){

    mini_x <- which.min(taus)
    tau <- taus[mini_x]

    ## If temperature hits boundary it is an event
    update_temp <- if(abs(tau - taus[nvel]) < eps) TRUE else FALSE

    x = x + tau*theta
    t = t + tau

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
        if(sampling_spike){
          theta[nvel] = -1
          sampling_spike <- FALSE
        } else {
          theta[nvel] = 0
          sampling_spike <- TRUE
        }
        event = TRUE
      }
      ## If regular event
      if(!b_hit_1 & !b_hit_0){
        if(u_s[nvel] < 1e-10){
          
          ## Evaluate the rate
          grad <- d_nlogq(x[-nvel], x[nvel], d_poly_coef)
          rate <- grad[nvel]*theta[nvel]

          ## Thinning
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
      ## Proceed with regular event
      
      if(u_s[mini_x] < 1e-10){
        
        ## Calculate rate
        grad <- d_nlogq(x[-nvel], x[nvel], d_poly_coef)
        rate <- grad[mini_x]*theta[mini_x]

        acc_prb <- rate/f_s[mini_x]
        if(acc_prb > 1.0001){
          print(paste("Invalid thinning on x, thinning prob:",acc_prb))
        }
        if(runif(1) <= acc_prb){
          theta[mini_x] = -theta[mini_x]
          event = TRUE
        }
      } else{
        nits = nits - 1
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
      temp_to_bndry <- time_to_bndry(x[nvel], theta[nvel])
      tau_grid <- seq(0, to = min(tau_max,temp_to_bndry),  length.out = poly_order + 1)

      rates <- return_rates(x, theta, tau_grid, d_poly_coef)
      for( i in 1:nvel){
        tus = sim_rate_poly(eval_times = tau_grid, eval_rates = rates[i,], poly_order)
        taus[i] = tus$t
        u_s[i] = tus$u
        f_s[i] = tus$f_evall
      }

      # Simulate time for temp
      if(sampling_spike){
        taus[nvel] <- rexp(1)/rate_jump_from_spike
      }

    } else {
      # If there was no event

      # Re-simulate times for all taus less than zero:
      temp_to_bndry <- time_to_bndry(x[nvel], theta[nvel])
      tau_grid <- seq(0, to = min(tau_max,temp_to_bndry),  length.out = poly_order + 1)

      # Adjust simulated times
      taus <- taus-tau
      update_rates <- which(taus <= 0)

      rates <- return_rates(x, theta, tau_grid, d_poly_coef)
      for( j in update_rates){
        tus = sim_rate_poly(eval_times = tau_grid, eval_rates = rates[j,], poly_order)
        taus[j] = tus$t
        u_s[j] = tus$u
        f_s[j] = tus$f_evall
      }

      # Simulate time for temp
      if(sampling_spike){
        taus[nvel] <- rexp(1)/rate_jump_from_spike
      }
    }
    nits = nits +1
  }
  if(num_evts < max_events){
    return (list(positions=positions[,1:num_evts],thetas=thetas[,1:num_evts],times=times[1:num_evts],
                 nits = nits, poly_coef = poly_coef, alpha = alpha))
  } else{
    return (list(positions=positions,thetas=thetas,times=times, nits = nits,
                 poly_coef = poly_coef, alpha = alpha))
  }
}
