data { 
  int<lower=0> K; // number of species
  int<lower=0> N; // length of time-series
  vector[K] x[N];
  vector[K] y[N]; 
  vector[K] init_priors;
  int<lower=0> J1;
  matrix[N,J1] xx1; // steelhead covariates
  int<lower=0> J2;
  matrix[N,J2] xx2; // dolly varden covariates
  int<lower=0> J3;
  matrix[N,J3] xx3; // cutthroat covariates
  int<lower=0> J4;
  matrix[N,J4] xx4; // pink salmon covariates
  int<lower=0> J5;
  matrix[N,J5] xx5; // coho covariates
  // steelhead life cycle below
  int<lower=0> N_obs; // missing steelhead survival year
  int M; // number of survival covariates
  matrix[N,M] x1; // what impacts steelhead survival
  int P; // number of steelhead adult covariates
  int Q; // number of steelhead run timing covariates
  matrix[N,Q] x2; // what impacts run times
  vector[N_obs] y1_obs; // the missing survival year
  vector[N] y2; // adult run time vector
  //vector[N] x3; // vector of time-varying covariates
  vector[N] juvCoh; // outmigrating smolt cohort size
  real init_s0; // prior for missing survival
}

parameters {
  vector[J1+1] beta_steel;
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
  // steelhead life cycle 
  vector[M] beta_surv;
  vector[Q+1] beta_run;
  vector[P] beta_adults;
  real<lower=0> sigma_surv_obs;
  real<lower=0> sigma_surv_pro;
  real<lower=0> sigma_run_obs;
  real<lower=0> sigma_run_pro;
  real<lower=0> sigma_adult_obs;
  real<lower=0> sigma_adult_pro;
  vector[N] s0;
  vector<lower=0>[N] r0;
  vector<lower=0>[N] a0;
  real y1_miss;
}

transformed parameters{
  vector[K] mu[N];
  vector[N] pred_surv;
  vector<lower=0>[N] pred_adults;
  vector[N] pred_run;
  vector[N] y1;
  
  y1[1] = y1_miss;
  y1[2:N] = y1_obs;
  
  for(i in 1:N)
  {
    pred_surv[i] = s0[i] + x1[i] * beta_surv;
  }
  
  for(i in 1:N)
  {
    pred_adults[i] = a0[i] + (pred_surv[i] - mean(pred_surv))/sd(pred_surv) * beta_adults[1] + juvCoh[i] * beta_adults[2];
  }
  
  for(i in 1:N)
  {
    pred_run[i] = r0[i] + (log(pred_adults[i])-mean(log(pred_adults)))/sd(log(pred_adults)) * beta_run[1] + x2[i] * beta_run[2:(Q+1)];
  }
  
  for (n in 1:N)
  {
    for(k in 1:K)
    {
      if(k==1)
      {
        mu[n,k] = x[n,k]*beta[k] + x0[n,k] + xx1[n] * beta_steel[1:J1] + (pred_run[n] - mean(pred_run))/sd(pred_run) * beta_steel[J1+1];
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
  // trend in marine survival
  y1_miss ~ normal(0,1);
  beta_surv ~ normal(0,5);
  s0[1] ~ normal(init_s0,sigma_surv_pro);
  for(i in 2:N)
  {
    s0[i] ~ normal(s0[i-1],sigma_surv_pro);
  }
  sigma_surv_pro ~ cauchy(0,5);
  sigma_surv_obs ~ cauchy(0,5);

  // trend in steelhead adults
  beta_adults ~ normal(0,20);
  a0[1] ~ normal(mean(x[1:N,1]),sigma_adult_pro);
  for(i in 2:N)
  {
    a0[i] ~ normal(a0[i-1],sigma_adult_pro);
  }
  sigma_adult_obs ~ cauchy(0,0.5);
  sigma_adult_pro ~ cauchy(0,15);
  
  // trend in adult run time
  beta_run ~ normal(0,15);
  r0[1] ~ normal(mean(y2),sigma_run_pro);
  for(i in 2:N)
  {
    r0[i] ~ normal(r0[i-1],sigma_run_pro);
  }
  sigma_run_pro ~ cauchy(0,50);
  sigma_run_obs ~ cauchy(0,50);
  
  // multivariate time-series for salmon community
  beta_steel ~ normal(0,5);
  beta_dolly ~ normal(0,5);
  beta_cutt ~ normal(0,5);
  beta_pink ~ normal(0,5);
  beta_coho ~ normal(0,5);
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
  // likelihood below
  y1 ~ normal(pred_surv,sigma_surv_obs);
  x[1:N,1] ~ lognormal(log(pred_adults),sigma_adult_obs);
  y2 ~ normal(pred_run,sigma_run_obs);
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
  vector[N] log_lik1;
  vector[N] log_lik2;
  vector[N] log_lik3;
  vector[N] y1_ppd;
  vector[N] x3_ppd;
  vector[N] y2_ppd;
  
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
      R[i,k] = exp(y_ppd[i,k])*x[i,k];
    }
    log_lik1[i] = normal_lpdf(y1[i] | pred_surv[i],sigma_surv_obs);
    log_lik2[i] = lognormal_lpdf(x[i,1] | log(pred_adults[i]),sigma_adult_obs);
    log_lik3[i] = normal_lpdf(y2[i] | pred_run[i],sigma_run_obs);
    y1_ppd[i] = inv_logit(normal_rng(pred_surv[i],sigma_surv_obs));
    x3_ppd[i] = lognormal_rng(log(pred_adults[i]),sigma_adult_obs);
    y2_ppd[i] = normal_rng(pred_run[i],sigma_run_obs);
  }
}
