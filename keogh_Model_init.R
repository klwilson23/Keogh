keogh_JAGS <- "model {
  
  for(j in 1:Nspecies) # the likelihood for recruitment dynamics
  { 
    for(i in 1:Nyears)
    {
      # recruitment is a lag-1 process
      recruits[i,j] ~ dnorm(recHat[i,j],1/pow(recHat[i,j]*rec.cv[j],2))T(0,)
    }
  } # end the observation portion

  for(j in 1:Nspecies) # the regressions
  {
    recHat[1,j] ~ dnorm(10000,1/pow(sys.sd[j],2))T(0,)
    for(i in 2:(Nyears))
    {
      a_t[i,j] ~ dlnorm(log(a_t[i-1,j]),tau.a[j])
      b_t[i,j] ~ dnorm(rho_b[j]*b_t[i-1,j],tau.b[j])T(0,)
      muRec[i,j] <-  a_t[i,j]*stock[i,j]*exp(-b_t[i,j]*stock[i,j])
      #muRec[i,j] <-  rho_Rec[j]*recHat[i-1,j] + a_t[i,j]*stock[i,j]*exp(-b_t[i,j]*stock[i,j])
      recHat[i,j] ~ dnorm(muRec[i,j],1/pow(sys.sd[j],2))T(0,)
    }
  } # end the process portion

  for(j in 1:Nspecies) {
    rho_a[j] ~ dunif(0,1) # autocorrelation
    rho_b[j] ~ dunif(0,1) # autocorrelation
    rho_Rec[j] ~ dunif(0,1) # autocorrelation
    sys.sd[j] ~ dgamma(1e-3,1e-3) # sd in system process
    rec.cv[j] ~ dgamma(1e-3,1e-3) # cv in recruitment
    a_t[1,j] ~ dlnorm(log(10),tau.a[j])
    b_t[1,j] ~ dnorm(0.00023,tau.b[j])T(0,)
    tau.b[j] ~ dgamma(1e-3,1e-3)
    tau.a[j] ~ dgamma(1e-3,1e-3)
  }
  
}
"
write(keogh_JAGS,"keogh_d1.jags")