data {
  int N;
  int K;
  matrix[N,K] X;
  vector[N] lSurv;
  vector[N] run_time;
}

parameters {
  vector[K] beta;
  vector[N-1] pro_devS;
  real pS0;
  real obs_sigma_surv;
  real pro_sigma_surv;
  real obs_sigma_run;
  real bSurv;
  real run0;
}

transformed parameters{
  vector[N] mnSurv;
  vector[N] obSurv;
  vector[N] mnRun;

  mnSurv[1] = pS0;
  for(t in 2:N) {
    mnSurv[t] = mnSurv[t-1] + pro_devS[t-1];
  }
  obSurv = X * beta + mnSurv;
  mnRun = bSurv * obSurv + run0;
}

model {
  // declare priors on initial states
  run0 ~ normal(mean(run_time),10);
  pS0 ~ normal(0,5);
  // variances for process and observation error
  pro_sigma_surv ~ cauchy(0,5);
  obs_sigma_surv ~ cauchy(0,5);
  obs_sigma_run ~ cauchy(0,5);
  // declare time-varying process errors
  pro_devS ~ normal(0, pro_sigma_surv);
  // regression coefficients
  beta ~ normal(0,10);
  bSurv ~ normal(0,10);
  // likelihood below
  lSurv ~ normal(obSurv,obs_sigma_surv);
  run_time ~ normal(mnRun,obs_sigma_run);
}

generated quantities {
  vector[N] surv_new;
  vector[N] run_new;
  for(i in 1:N) {
    surv_new[i] = inv_logit(normal_rng(obSurv[i],obs_sigma_surv));
    run_new[i] = normal_rng(mnRun[i],obs_sigma_run);
  }
}
