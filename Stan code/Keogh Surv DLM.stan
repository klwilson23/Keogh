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
  vector[J] beta_run;
  vector[N-1] pro_devS;
  vector[N-1] pro_devR;
  real pS0;
  real<lower=0> obs_sigma_surv;
  real<lower=0> pro_sigma_surv;
  real<lower=0> obs_sigma_run;
  real<lower=0> pro_sigma_run;
  real bSurv;
  real run0;
}

transformed parameters{
  vector[N] mnSurv;
  vector[N] obSurv;
  vector[N] mnRun;
  vector[N] obRun;

  mnSurv[1] = pS0;
  mnRun[1] = run0;
  for(t in 2:N) {
    mnSurv[t] = mnSurv[t-1] + pro_devS[t-1];
    mnRun[t] = mnRun[t-1] + pro_devR[t-1];
  }
  obSurv = mnSurv + X * beta_surv;
  obRun = mnRun + mnSurv * bSurv + XX * beta_run;
}

model {
  // declare priors on initial states
  run0 ~ normal(mean(run_time),25);
  pS0 ~ normal(0,5);
  // variances for process and observation error
  pro_sigma_surv ~ cauchy(0,10);
  obs_sigma_surv ~ cauchy(0,10);
  pro_sigma_run ~ cauchy(0,50);
  obs_sigma_run ~ cauchy(0,50);
  // declare time-varying process errors
  pro_devS ~ normal(0, pro_sigma_surv);
  pro_devR ~ normal(0, pro_sigma_run);
  // regression coefficients
  beta_surv ~ normal(0,5);
  beta_run ~ normal(0,5);
  bSurv ~ normal(0,5);
  // likelihood below
  lSurv ~ normal(obSurv,obs_sigma_surv);
  run_time ~ normal(obRun,obs_sigma_run);
  
}

generated quantities {
  vector[N] lsurv_new;
  vector[N] surv_new;
  vector[N] run_new;
  for(i in 1:N) {
    lsurv_new[i] = normal_rng(obSurv[i],obs_sigma_surv);
    surv_new[i] = inv_logit(lsurv_new[i]);
    run_new[i] = normal_rng(obRun[i],obs_sigma_run);
  }
}
