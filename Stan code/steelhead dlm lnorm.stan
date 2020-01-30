functions {
  real normal_lb_ub_rng(real mu, real sigma, real lb, real ub) 
  {
    real p1 = normal_cdf(lb, mu, sigma);  // cdf with lower bound
    real p2 = normal_cdf(ub, mu, sigma);  // cdf with upper bound
    real u = uniform_rng(p1, p2);
    return (sigma * inv_Phi(u)) + mu;  // inverse cdf 
  }
}

data {
  int<lower=0> N;
  int<lower=0> N_obs;
  int K;
  matrix[N,K] x1;
  int J;
  matrix[N,J] x2;
  vector[N_obs] y1_obs;
  vector[N] y2;
  vector[N] y3;
  vector[N] x3; // vector of time-varying covariates
  int M;
  matrix[N,M] xx3;
  real init_s0;
}

parameters {
  
  vector[K] beta_surv;
  vector[J+1] beta_run;
  real beta_adults;
  
  real<lower=0> sigma_surv_obs;
  real<lower=0> sigma_surv_pro;
  real<lower=0> sigma_run_obs;
  real<lower=0> sigma_run_pro;
  real<lower=0> sigma_adult_obs;
  real<lower=0> sigma_adult_pro;
  vector[N] s0;
  vector<lower=0>[N] r0;
  vector<lower=0>[N] a0;
  vector<lower=0>[N] x0;
  real<upper=0> beta_rec;
  real<lower=0> sigma_rec_process;
  real<lower=0> sigma_rec_obs;
  vector[M+1] beta_rec_cov;
  real y1_miss;
}

transformed parameters {
  vector[N] pred_surv;
  vector[N] pred_adults;
  vector[N] pred_run;
  vector[N] pred_rec;
  vector[N] y1;
  
  y1[1] = y1_miss;
  y1[2:N] = y1_obs;
  
  for(i in 1:N)
  {
    pred_surv[i] = s0[i] + x1[i] * beta_surv;
  }
  
  for(i in 1:N)
  {
    pred_adults[i] = a0[i] + (pred_surv[i] - mean(pred_surv))/sd(pred_surv) * beta_adults;
  }
  
  for(i in 1:N)
  {
    pred_run[i] = r0[i] + log(x3[i]) * beta_run[1] + x2[i] * beta_run[2:(J+1)];
  }
  for(i in 1:N)
  {
    pred_rec[i] = x3[i] * beta_rec + x0[i] + (pred_run[i] - mean(pred_run))/sd(pred_run) * beta_rec_cov[1] + xx3[i] * beta_rec_cov[2:(M+1)];
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
  beta_adults ~ normal(0,15);
  a0[1] ~ normal(mean(x3),sigma_adult_pro);
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

  // trend in productivity
  beta_rec ~ normal(-2e-3,1);
  x0[1] ~ normal(mean(y3),sigma_rec_process);
  for(i in 2:N)
  {
    x0[i] ~ normal(x0[i-1],sigma_rec_process);
  }
  sigma_rec_process ~ cauchy(0,5);
  sigma_rec_obs ~ cauchy(0,5);
  beta_rec_cov ~ normal(0,5);
  
  // likelihood below
  y1 ~ normal(pred_surv,sigma_surv_obs);
  x3 ~ lognormal(log(pred_adults),sigma_adult_obs);
  y2 ~ normal(pred_run,sigma_run_obs);
  y3 ~ normal(pred_rec, sigma_rec_obs);
}

generated quantities {
  vector[N] log_lik1;
  vector[N] log_lik2;
  vector[N] log_lik3;
  vector[N] log_lik4;
  vector[N] y1_ppd;
  vector[N] x3_ppd;
  vector[N] y2_ppd;
  vector[N] y3_ppd;
  vector[N] R;
  for(i in 1:N) {
    log_lik1[i] = normal_lpdf(y1[i] | pred_surv[i],sigma_surv_obs);
    log_lik2[i] = lognormal_lpdf(x3[i] | log(pred_adults[i]),sigma_adult_obs);
    log_lik3[i] = normal_lpdf(y2[i] | pred_run[i],sigma_run_obs);
    log_lik4[i] = normal_lpdf(y3[i] | pred_rec[i], sigma_rec_obs);
    y1_ppd[i] = inv_logit(normal_rng(pred_surv[i],sigma_surv_obs));
    x3_ppd[i] = lognormal_rng(log(pred_adults[i]),sigma_adult_obs);
    y2_ppd[i] = normal_rng(pred_run[i],sigma_run_obs);
    y3_ppd[i] = normal_rng(pred_rec[i], sigma_rec_obs);
    R[i] = exp(pred_rec[i])*x3_ppd[i];
  }
}
