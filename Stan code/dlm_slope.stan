data {
  int<lower=0> N;
  vector[N] y;
  int<lower=0> K; // number of time-varying covariates
  matrix[N, K] x; // matrix of time-varying covariates
  int<lower=0> J; // number of time-constant covariates
  matrix[N,J] c;
  int family; // 1 = normal, 2 = gamma, 3 = lognormal
  int type; // 0 = time-constant intercept, 1 = time-varying intercept
}
parameters {
  real x0;
  vector[K] beta0;
  vector[K] pro_dev[N-1];
  real<lower=0> sigma_process[K];
  real<lower=0> sigma_obs;
  vector[J] Cvars;
}
transformed parameters {
  vector[N] pred;
  vector[K] beta[N]; // elements accessed [N,K]

  for(k in 1:K) {
   beta[1,k] = beta0[k];
   for(i in 2:N) {
    beta[i,k] = beta[i-1,k] + pro_dev[i-1,k];
   }
  }
  if(type==0)
  {
    for(i in 1:N) {
      pred[i] = x[i] * beta[i] + c[i] * Cvars + x0;
    }
  }
  if(type==1)
  {
    for(i in 1:N) {
      pred[i] = x[i] * beta[i] + c[i] * Cvars;
    }
  }
}
model {
  if(type==1) 
  {
    beta0[1] ~ normal(mean(y),sd(y));
   if(K>1)
   {
     beta0[2:K] ~ normal(0,1);
   }
  }
  if(type==0)
  {
    beta0 ~ normal(0,1);
    x0 ~ normal(mean(y),sd(y));
  }
  Cvars ~ normal(0,1);
  sigma_obs ~ cauchy(0,sd(y));
  sigma_process ~ cauchy(0,5);
  for(k in 1:K) {
    pro_dev[k] ~ normal(0, sigma_process[k]);
  }
  if(family==1) y ~ normal(pred, sigma_obs);
  if(family==2) y ~ gamma(sigma_obs, sigma_obs ./ exp(pred));
  if(family==3) y ~ lognormal(pred, sigma_obs);
}
generated quantities {
vector[N] log_lik;
  vector[N] y_ppd;
  if(family == 1) {
    for(i in 1:N) {
      log_lik[i] = normal_lpdf(y[i] | pred[i], sigma_obs);
      y_ppd[i] = normal_rng(pred[i], sigma_obs);
    }
  }
    if(family == 2) {
    for(i in 1:N) {
      log_lik[i] = gamma_lpdf(y[i] | sigma_obs , sigma_obs ./ exp(pred[i]));
      y_ppd[i] = gamma_rng(sigma_obs , sigma_obs ./ exp(pred[i]));
    }
  }
    if(family == 3) {
    for(i in 1:N) {
      log_lik[i] = lognormal_lpdf(y[i] | pred[i], sigma_obs);
      y_ppd[i] = lognormal_rng(pred[i], sigma_obs);
    }
  }
}
