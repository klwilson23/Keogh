keogh_recruitment <- "model {

  for(j in 1:Nspecies) # the likelihood for recruitment dynamics
  {
    # set the 1st year
    recruits[1,j] ~ dnorm(muRec[1,j],1/(muRec[1,j]*rec.cv[j])^2)T(0,)
    muRec[1,j] <- a_t[1,j]*stock[1,j]*exp(-b_t[1,j]*stock[1,j])

    for(i in 2:Nyears)
    {
      # recruitment is a lag-1 process
      recruits[i,j] ~ dnorm(muRec[i,j],1/(muRec[i,j]*rec.cv[j])^2)T(0,)
      muRec[i,j] <- a_t[i,j]*stock[i-1,j]*exp(-b_t[i,j]*stock[i-1,j])
    }
  } # end the data portion


  for(j in 1:Nspecies) # the likelihood for stock dynamics
  {
    stock[1,j] ~ dnorm(muStock[1,j],1/(muStock[1,j]*stock.cv[j]))T(0,)
    muStock[1,j] <- recruits[i,j]*marSurv_t[i,j]
    for(i in 2:Nyears)
    {
      stock[i,j] ~ dnorm(muStock[i,j],1/(muStock[i,j]*stock.cv[j]))T(0,)
      muStock[i,j] <- recruits[i-1,j]*marSurv_t[i-1,j]
    }
  } # end the data portion

  for(j in 1:Nspecies) # the regressions
  {
    for(i in 2:Nyears)
    {
      la_t[i,j]) ~ dnorm(log(a[j])+alpha[j]*i,tau.la)
      b_t[i,j] ~ dnorm(b[j]+beta[j,1]*recruits[i-1,1]+beta[j,2]*recruits[i-1,2]+beta[j,3]*recruits[i-1,3]+beta[j,4]*recruits[i-1,4]+beta[j,5]*recruits[i-1,5],tau.b)
      mar.a <- mu_Mar[j,i] / (sd.surv * sd.surv)
      mar.b <- (1-mu_Mar[j,i]) / (sd.surv * sd.surv)
      marSurv_t[i,j] ~ dbeta(marSurv[j] + marTrend[j]*i, tau.Surv)
    }
  } # end the data portion

}"

write(keogh_recruitment,"keoghRec.jags")
