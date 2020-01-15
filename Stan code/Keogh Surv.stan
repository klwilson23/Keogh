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
  real<lower=0> sigma_surv;
  real<lower=0> sigma_run;
  real<lower=-1,upper=1> phi_surv;
  real<lower=-1,upper=1> phi_run;
}

transformed parameters{
  vector[N] mnSurv;
  vector[N] mnRun;
  vector[N] epsilonSurv;
  vector[N] epsilonRun;
  real sigma_cor_surv;
  real sigma_cor_run;
  mnSurv[1] = X[1] * beta_surv;
  epsilonSurv[1] = lSurv[1] - mnSurv[1];
  mnRun[1] = XX[1] * beta_run;
  epsilonRun[1] = run_time[1] - mnRun[1];
  for(i in 2:N) {
    mnSurv[i] = X[i] * beta_surv + phi_surv*epsilonSurv[i-1];
    epsilonSurv[i] = (lSurv[i] - mnSurv[i]) - phi_surv*epsilonSurv[i-1];
    mnRun[i] = XX[i] * beta_run + phi_run*epsilonRun[i-1];
    epsilonRun[i] = (run_time[i] - mnRun[i]) - phi_run*epsilonRun[i-1];
  }
  sigma_cor_surv = sqrt( sigma_surv*sigma_surv * (1-phi_surv*phi_surv));
  sigma_cor_run = sqrt( sigma_run*sigma_run * (1-phi_run*phi_run));
}

model {
  // declare priors on autocorrelation
  phi_surv ~ normal(0,1);
  phi_run ~ normal(0,1);
  // regression coefficients
  beta_surv[1] ~ normal(mean(lSurv),sd(lSurv));
  beta_run[1] ~ normal(mean(run_time),sd(run_time));
  if(K>1) beta_surv[2:K] ~ normal(0,5);
  if(J>1) beta_run[2:J] ~ normal(0,5);
  // variances for process error
  sigma_surv ~ cauchy(0,sd(lSurv));
  sigma_run ~ cauchy(0,sd(run_time));
  // likelihood below
  lSurv ~ normal(mnSurv,sigma_cor_surv);
  run_time ~ normal(mnRun,sigma_cor_run);
}

generated quantities {
  vector[N] log_lik1;
  vector[N] log_lik2;
  vector[N] surv_new;
  vector[N] run_new;
  for(i in 1:N) {
    surv_new[i] = inv_logit(normal_rng(mnSurv[i],sigma_cor_surv));
    run_new[i] = normal_rng(mnRun[i],sigma_cor_run);
    log_lik1[i] = normal_lpdf(lSurv[i] | mnSurv[i],sigma_cor_surv);
    log_lik2[i] = normal_lpdf(run_time[i] | mnRun[i],sigma_cor_run);
  }
}
