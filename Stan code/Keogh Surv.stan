data {
  int N;
  int K;
  matrix[N,K] X;
  vector[N] lSurv;
  int J;
  matrix[N,J] XX;
  vector[N] run_time;
}

parameters {
  vector[K] beta_surv;
  real pS0;
  vector[J] beta_run;
  real run0;
  real<lower=0> sigma_surv;
  real<lower=0> sigma_run;
}

transformed parameters{
  vector[N] mnSurv;
  vector[N] mnRun;
  mnSurv = pS0 + X * beta_surv;
  mnRun = run0 + XX * beta_run;
}

model {
  // declare priors on initial states
  run0 ~ normal(mean(run_time),50);
  pS0 ~ normal(0,5);
  // variances for process and observation error
  sigma_surv ~ cauchy(0,5);
  sigma_run ~ cauchy(0,5);
  // regression coefficients
  beta_surv ~ normal(0,5);
  beta_run ~ normal(0,5);
  // likelihood below
  lSurv ~ normal(mnSurv,sigma_surv);
  run_time ~ normal(mnRun,sigma_run);
  
}

generated quantities {
  vector[N] surv_new;
  vector[N] run_new;
  for(i in 1:N) {
    surv_new[i] = inv_logit(normal_rng(mnSurv[i],sigma_surv));
    run_new[i] = normal_rng(mnRun[i],sigma_run);
  }
}
