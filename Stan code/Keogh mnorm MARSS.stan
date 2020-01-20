data { 
  int<lower=0> K; // number of species
  int<lower=0> N; // length of time-series
  vector[K] x[N];
  vector[K] y[N]; 
  vector[K] init_priors;
}

parameters { 
  vector<upper=0>[K] beta;
  vector[K] x0[N];
  cholesky_factor_corr[K] L_Omega_obs;
  cholesky_factor_corr[K] L_Omega_proc;  
  vector<lower=0>[K] L_sigma_obs;
  vector<lower=0>[K] L_sigma_proc;
}

transformed parameters{
  vector[K] mu[N];
  for (n in 1:N)
  {
    for(k in 1:K)
    {
      mu[n,k] = x[n,k]*beta[k] + x0[n,k];
    }
  }
}

model {
  beta ~ normal(-2e-3, 1);
  L_Omega_obs ~ lkj_corr_cholesky(1); 
  L_sigma_obs ~ student_t(3, 0, 2);
  L_Omega_proc ~ lkj_corr_cholesky(1); 
  L_sigma_proc ~ student_t(3, 0, 2);
  x0[1] ~ multi_normal_cholesky(init_priors,diag_pre_multiply(L_sigma_proc, L_Omega_proc));
  for(i in 2:N)
  {
    x0[i] ~ multi_normal_cholesky(x0[i-1],diag_pre_multiply(L_sigma_proc, L_Omega_proc));
  }
  y ~ multi_normal_cholesky(mu, diag_pre_multiply(L_sigma_obs, L_Omega_obs));
}
generated quantities {
  matrix[K, K] Omega;
  matrix[K, K] Sigma;
  matrix[K, K] Omega_proc;
  matrix[K, K] Sigma_proc;
  vector[N] log_lik;
  vector[K] y_ppd[N]; // recruitment productivity
  vector[K] R[N];
  Omega = multiply_lower_tri_self_transpose(L_Omega_obs);
  Sigma = quad_form_diag(Omega, L_sigma_obs);
  Omega_proc = multiply_lower_tri_self_transpose(L_Omega_proc);
  Sigma_proc = quad_form_diag(Omega_proc, L_sigma_proc);
  for(i in 1:N)
  {
    log_lik[i] = multi_normal_cholesky_lpdf(y[i] | mu[i],diag_pre_multiply(L_sigma_obs, L_Omega_obs));
    y_ppd[i] = multi_normal_cholesky_rng(mu[i],diag_pre_multiply(L_sigma_obs, L_Omega_obs));
    for(k in 1:K)
    {
      R[i,k] = exp(x0[i,k]) * x[i,k] * exp(beta[k]*x[i,k]);
    }
  }
}
