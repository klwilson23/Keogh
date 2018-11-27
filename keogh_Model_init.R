keogh_JAGS <- "model {
  
  for(i in 1:Nyears) # the likelihood for recruitment dynamics
  {
    for(j in 1:Nspecies)
    {
      # recruitment is a lag-1 process
      recruits[i,j] ~ dnorm(muRec[i,j],1/pow(muRec[i,j]*rec.cv[j],2))T(0,)
      muRec[i,j] <- a_t[i+1,j]*stock[i,j]*exp(-b_t[i+1,j]*stock[i,j])
    }
  } # end the observation portion
  sd.la ~ dhalfcauchy(25) # variance in productivity
  sd.b ~ dhalfcauchy(25) # variance in DD-survival slope
  
  for(j in 1:Nspecies) # the regressions
  {
    rho_a[j] ~ dunif(-1,1) # autocorrelation
    rho_b[j] ~ dunif(-1,1) # autocorrelation
    rec.cv[j] ~ dunif(0,0.6) # cv in recruitment
    la_t[1,j] <- log(10)
    a_t[1,j] <- exp(la_t[1,j])
    b_t[1,j] <- 0.0002302585
    for(i in 2:(Nyears+1))
    {
      la_t[i,j] ~ dnorm(rho_a[j]*la_t[i-1,j],1/(pow(sd.la,2)))
      a_t[i,j] <- exp(la_t[i,j])
      b_t[i,j] ~ dnorm(rho_b[j]*b_t[i-1,j],1/(pow(sd.b,2)))
    }
  } # end the process portion
  
}
"
write(keogh_JAGS,"keogh_d1.jags")