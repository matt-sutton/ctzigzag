//This file contains the code required to return gradients and log denisty evaluations 
// ctZigZag uses this information in the thinning

data{
  int<lower=1> K;         // number of mixture components
  int<lower=1> d;         // dimension
  row_vector[d] mu[K];     // Target mus
  vector[K] sigma;         // Target scales
  vector[d] mu0;           // Base mu
  vector[d] sigma0;        // Base precision
}

parameters{
  // Configuration state.
  vector[d] theta;
  real beta;
}

model {
  real ps[K];

 for(k in 1:K){
    ps[k] = log(1.0/K) + multi_normal_lpdf(theta|mu[k],diag_matrix(rep_vector(sigma[k],d)));
  }

  target += beta*log_sum_exp(ps);
  target += (1.0-beta)*multi_normal_lpdf(theta|mu0,diag_matrix(sigma0));
}
