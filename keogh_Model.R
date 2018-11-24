model {

  for(j in 1:Nspecies) # the likelihood for recruitment dynamics
  {
    for(i in 1:Nyears)
    {
      # recruitment is a lag-1 process
      recruits[i,j] ~ dnorm(muRec[i,j],1/(muRec[i,j]*rec.cv[j])^2)T(0,)
      muRec[i,j] <- a_t[i,j]*obsStock[i,j]*exp(-b_t[i,j]*obsStock[i,j])
    }
  } # end the data portion


  for(j in 1:Nspecies) # the likelihood for stock dynamics
  {
    stock[1,j] ~ dnorm(muStock[1,j],1/(muStock[1,j]*stock.cv[j]))T(0,)
    muStock[1,j] <- obsRec[1,j]*marSurv_t[1,j]
    for(i in 2:Nyears)
    {
      stock[i,j] ~ dnorm(muStock[i,j],1/(muStock[i,j]*stock.cv[j]))T(0,)
      muStock[i,j] <- obsRec[i-1,j]*marSurv_t[i-1,j]
    }
  } # end the data portion

  for(j in 1:Nspecies) # the regressions
  {
    la_t[1,j] ~ dnorm(log(a[j]),tau.la[j])
    a_t[1,j] <- exp(la_t[1,j])
    b_t[1,j] ~ dnorm(b[j],tau.b[j])
    
    # marine survival from smolt to adult is a logit regression
    logit(mu_Mar[1,j]) <- -log(marSurv[j]/(1-marSurv[j]))+marTrend[j]*(0)
    mar.a[1,j] <- mu_Mar[1,j] / (sd.surv[j] * sd.surv[j]) # convert to JAGS dbeta shape parameters
    mar.b[1,j] <- (1-mu_Mar[1,j]) / (sd.surv[j] * sd.surv[j])
    
    # marine survival follows a beta distribution
    marSurv_t[1,j] ~ dbeta(mar.a[1,j], mar.b[1,j])
    
    for(i in 2:Nyears)
    {
      la_t[i,j] ~ dnorm(log(a[j])+alpha[j]*(i-1),tau.la[j])
      a_t[i,j] <- exp(la_t[i,j])
      b_t[i,j] ~ dnorm(b[j]+beta[j,1]*obsRec[i-1,1]+beta[j,2]*obsRec[i-1,2]+beta[j,3]*obsRec[i-1,3]+beta[j,4]*obsRec[i-1,4]+beta[j,5]*obsRec[i-1,5],tau.b[j])
      
      # marine survival from smolt to adult is a logit regression
      logit(mu_Mar[i,j]) <- -log(marSurv[j]/(1-marSurv[j]))+marTrend[j]*(i-1)
      mar.a[i,j] <- mu_Mar[i,j] / (sd.surv[j] * sd.surv[j]) # convert to JAGS dbeta shape parameters
      mar.b[i,j] <- (1-mu_Mar[i,j]) / (sd.surv[j] * sd.surv[j])

      # marine survival follows a beta distribution
      marSurv_t[i,j] ~ dbeta(mar.a[i,j], mar.b[i,j])
    }
  } # end the data portion

  for(j in 1:Nspecies)
  {
    sd.surv[j] ~ dgamma(1e-3,1e-3) # variance in marine survival
    marSurv[j] ~ dbeta(1e-3,1e-3) # baseline marine survival
    marTrend[j] ~ dnorm(0,1e-3) # slope in marine marine survival
    alpha[j] ~ dnorm(0,1e-3) # slope in average productivity
    tau.la[j] ~ dgamma(1e-3,1e-3) # variance in productivity
    tau.b[j] ~ dgamma(1e-3,1e-3) # variance in DD-survival slope
    stock.cv[j] ~ dgamma(1e-3,1e-3) # cv in recruitment
    rec.cv[j] ~ dgamma(1e-3,1e-3) # cv in recruitment
    b[j] ~ dnorm(0,1e-3) # DD-survival slope
    a[j] ~ dnorm(0,1e-3) # intercept
  }

  for(i in 1:Nspecies){
    for(j in 1:Nspecies){
      beta[i,j] ~ dnorm(0, 1e-3)
    }
  }
  

}

