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
    rho_a[j] ~ dunif(-1,1) # autocorrelation
    rho_b[j] ~ dunif(-1,1) # autocorrelation
    rho_Rec[j] ~ dnorm(0,1e-3)
    rec.cv[j] ~ dunif(0,0.6) # cv in recruitment
    a_t[1,j] ~ dnorm(10,1e-3)T(0,)
    b_t[1,j] ~ dnorm(0.00023,1e-1)T(0,)
    recHat[1,j] ~ dnorm(recruits[1,j],1/pow(recruits[1,j]*rec.cv[j],2))T(0,)
    for(i in 2:(Nyears))
    {
      log(a_t[i,j]) <- rho_a[j]*log(a_t[i-1,j])
      b_t[i,j] <- rho_b[j]*b_t[i-1,j]
      muRec[i,j] <-  rho_Rec[j]*recHat[i-1,j] + a_t[i,j]*stock[i,j]*exp(-b_t[i,j]*stock[i,j])
      recHat[i,j] ~ dnorm(muRec[i,j],1/pow(muRec[i,j]*rec.cv[j],2))T(0,)
    }
  } # end the process portion
  
}
"
write(keogh_JAGS,"keogh_d1.jags")