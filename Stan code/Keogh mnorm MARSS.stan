data { 
  int<lower=0> K; // number of species
  int<lower=0> N; // length of time-series
  vector[K] x[N];
  vector[K] y[N]; 
  vector[K] init_priors;
  int<lower=0> J1;
  matrix[J1,N] xx1;
  int<lower=0> J2;
  matrix[J2,N] xx2;
  int<lower=0> J3;
  matrix[J3,N] xx3;
  int<lower=0> J4;
  matrix[J4,N] xx4;
  int<lower=0> J5;
  matrix[J5,N] xx5;
}

parameters {
  vector[J1] beta_steel;
  vector[J2] beta_dolly;
  vector[J3] beta_cutt;
  vector[J4] beta_pink;
  vector[J5] beta_coho;
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
      if(k==1)
      {
        mu[n,k] = x[n,k]*beta[k] + x0[n,k] + xx1[n] * beta_steel;
      }
      if(k==2)
      {
        mu[n,k] = x[n,k]*beta[k] + x0[n,k] + xx2[n] * beta_dolly;
      }
      if(k==3)
      {
        mu[n,k] = x[n,k]*beta[k] + x0[n,k] + xx3[n] * beta_cutt;
      }
      if(k==4)
      {
        mu[n,k] = x[n,k]*beta[k] + x0[n,k] + xx4[n] * beta_pink;
      }
      if(k==5)
      {
        mu[n,k] = x[n,k]*beta[k] + x0[n,k] + xx5[n] * beta_coho;
      }
    }
  }
}

model {
  beta_steel ~ normal(0,1);
  beta_dolly ~ normal(0,1);
  beta_cutt ~ normal(0,1);
  beta_pink ~ normal(0,1);
  beta_coho ~ normal(0,1);
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
