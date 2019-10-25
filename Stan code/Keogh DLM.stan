data {
  int<lower=0> N; // number of years
  int<lower=0> n_obs; // number of datapoints
  int<lower=0> M; // number of species
  int<lower=0> K; // number of process covariates
  int<lower=0> J; // number of observation covariates
  int<lower=0> yearID[n_obs];
  int<lower=0> speciesID[n_obs];
  vector[n_obs] y; // data
  matrix[K,n_obs] procCovar; // design matrix for process covariates
  matrix[J,n_obs] obsCovar; // design matrix for observation covariates
  int slopeRand; // 1 = constant, 2 = time-varying
  int intRand; // 1 = constant, 2 = time-varying
  int multiVarObs; // 1 = independent deviates, 2 = multivariate deviates
  int multiVarProc; // 1 = independent deviates, 2 = multivariate deviates
}
parameters {
  vector[M] x0; // initial states
  vector[S] pro_dev[N-1];
  vector[n_trends] U;
  real<lower=0> sigma_process[S];
  real<lower=0> sigma_obs[n_obsvar];
}
transformed parameters {
  vector[M] pred[N];
  vector[S] x[N]; // elements accessed [N,K]
  matrix[M,N] intercept;
  matrix[M,N] slope;
  // random walk in intercept
  for(m in 1:M){
    if(intRand==2){
      intercept[m,1] = x0[m];
      for(i in 2:N) {
        intercept[m,i] = intercept[m,i-1] + int_pro_dev[m,i-1];
        }
    }
    if(intRand==1){
      for(i in 1:N){
        intercept[m,i] = x0[m];
      }
    }
  }
  
  // random walk in slope
  for(m in 1:M){
    if(slopeRand==2){
      slope[m,1] = slope0[m];
      for(i in 2:N) {
        slope[m,i] = slope[m,i-1] + slope_pro_dev[m,i-1];
        }
    }
    if(slopeRand==1){
      for(i in 1:N){
        slope[m,i] = slope0[m];
      }
    }
  }
  
  // map predicted states to time series
  for(m in 1:M){
    for(i in 2:N) {
      pred[m,i] = x[i] * slope[m,i] + intercept[m,i];
    }
  }
}
model {
  x0 ~ normal(0,10);
  for(i in 1:n_obsvar) {
    sigma_obs[i] ~ cauchy(0,5);
  }
  
  if(multiVarProc==2) {
    Sigma_proc <- diag_matrix(tauProc) * Omega_proc * diag_matrix(tauProc);
    pro_dev[s,i] ~ multi_normal(0, Sigma_proc);
  }
  
  if(multiVarProc==1) {
    for(s in 1:n_provar) {
      sigma_process[s] ~ cauchy(0,5); // process var
      }
    for(s in 1:S) {
      pro_dev[s] ~ normal(0, sigma_process[proVariances[s]]); // process deviations
      }
  }
  
  // likelihood
  for(i in 1:n_obs) {
    y[i] ~ normal(pred[speciesID[i], yearID[i]], sigma_obs[speciesID[i]]);
  }
}
generated quantities {
  vector[n_pos] log_lik;
  // regresssion example in loo() package
  for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | pred[col_indx_pos[n], row_indx_pos[n]], sigma_obs[obsVariances[row_indx_pos[n]]]);
}


