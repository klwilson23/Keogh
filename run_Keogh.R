#library(jagsUI)
#library(beepr)
library(coda)
library(rjags)
library(runjags)

dat <- readRDS("stockRec.rds")

# Set up MCMC -------------------------------------------------------------

nSamps <- 100
nBurnin <- 800
nThin <- 2
nChains <- 2
nAdapt <- 1000
parcomp <- FALSE

# Send to JAGS ------------------------------------------------------------

gen_ini <- function()
{
  sd.surv <- rgamma(Nspecies,1e-3,1e-3) # variance in marine survival
  marSurv <- rbeta(Nspecies,1,5) # baseline marine survival
  marTrend <- rnorm(Nspecies,0,1e-3) # slope in marine marine survival
  alpha <- rnorm(Nspecies,0,1e-3) # average productivity
  tau.la <- rgamma(Nspecies,1e-3,1e-3) # variance in productivity
  tau.b <- rgamma(Nspecies,1e-3,1e-3) # variance in DD-survival slope
  rec.cv <- rgamma(Nspecies,1e-3,1e-3) # cv in recruitment
  stock.cv <- rgamma(Nspecies,1e-3,1e-3) # cv in recruitment
  beta <- matrix(rnorm(Nspecies*Nspecies,0,1e-3),ncol=Nspecies,nrow=Nspecies)
  diag(beta) <- 0 # set the diagonals equal to 0
  b <- -1*abs(rnorm(Nspecies,0,1e-3)) # DD-survival slope
  a <- abs(rnorm(Nspecies,0,1e-3)) # intercept
  b_t <- matrix(-1*abs(rnorm(Nspecies,0,1e-3)),Nyears+1,Nspecies)
  la_t <- matrix(rnorm(Nspecies,0,1e-3),Nyears+1,Nspecies)
  marSurv_t <- matrix(rbeta(Nspecies,1,5),Nyears+1,Nspecies)
  return(list("sd.surv"=sd.surv,"marSurv"=marSurv,"marTrend"=marTrend,"alpha"=alpha,"tau.la"=tau.la,"tau.b"=tau.b,"rec.cv"=rec.cv,"stock.cv"=stock.cv,"beta"=beta,"b"=b,"a"=a,"b_t"=b_t,"la_t"=la_t,"marSurv_t"=marSurv_t))
}

ini <- lapply(1:nChains,FUN=function(x){gen_ini()})

mon_names <- names(ini[[1]])

results <- run.jags(model="keogh_Model.R",monitor=mon_names, data=dat, n.chains=nChains, method="rjags", inits=ini, plots=F,silent.jag=F, modules=c("glm","bugs"),sample=nSamps,adapt=nAdapt,burnin=nBurnin,thin=nThin,summarise=FALSE)


Nspecies <- 5
beta <- matrix(NA,Nspecies,Nspecies)

for(i in 1:Nspecies){
  beta[i,i] <- 0
}

for(i in 1:(Nspecies-1)){
  beta[i+1,i] <- 1
  beta[i-1,i] <- 1
  beta[i,Nspecies] <- 1
}

for(i in 2:(Nspecies-2)){
  beta[i+2,i] <- 1
  beta[i-2,i] <- 1
}
beta[3,1] <- 1
beta[4,1] <- 1
beta[5,1] <- 1
beta[5,2] <- 1
beta[1,4] <- 1
beta[2,4] <- 1
