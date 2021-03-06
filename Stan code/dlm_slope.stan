data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x; // vector of time-varying covariates
  int family; // 1 = normal, 2 = lognormal
}

parameters {
  real<lower=0> x0;
  vector<upper=0>[N] beta;
  real<lower=0> sigma_process;
  real<lower=0> sigma_obs;
}

transformed parameters {
  vector[N] pred;
  for(i in 1:N)
  {
    pred[i] = x[i]*beta[i] + x0;
  }
}

model {
  x0 ~ normal(mean(y),1);
  beta[1] ~ normal(-2e-3,sigma_process);
  for(i in 2:N)
  {
    beta[i] ~ normal(beta[i-1],sigma_process);
  }
  sigma_process ~ cauchy(0,5);
  sigma_obs ~ cauchy(0,10);
  if(family==1)
  {
    y ~ normal(pred, sigma_obs);
  }
  if(family==2)
  {
    y ~ lognormal(pred, sigma_obs);
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_ppd;
  vector[N] R;
  if(family==1)
  {
    for(i in 1:N) {
      log_lik[i] = normal_lpdf(y[i] | pred[i], sigma_obs);
      y_ppd[i] = normal_rng(pred[i], sigma_obs);
      R[i] = exp(x0)*x[i]*exp(beta[i]*x[i]);
    } 
  }
  if(family==2)
  {
    for(i in 1:N) {
      log_lik[i] = lognormal_lpdf(y[i] | pred[i], sigma_obs);
      y_ppd[i] = lognormal_rng(pred[i], sigma_obs);
      R[i] = exp(x0)*x[i]*exp(beta[i]*x[i]);
    } 
  }
}
