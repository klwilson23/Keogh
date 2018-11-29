#library(jagsUI)
#library(beepr)
library(coda)
library(rjags)
library(runjags)

dat <- readRDS("stockRec.rds")

# Set up MCMC -------------------------------------------------------------

nSamps <- 1200
nThin <- 10
nBurnin <- nSamps*nThin
nChains <- 2
nAdapt <- round(0.1*nBurnin)

parcomp <- FALSE

Nspecies <- dat$Nspecies
Nyears <- dat$Nyears
# Send to JAGS ------------------------------------------------------------

gen_ini <- function()
{
  sd.surv <- rgamma(Nspecies,1,1e-1) # variance in marine survival
  marSurv <- rbeta(Nspecies,1,5) # baseline marine survival
  marTrend <- rnorm(Nspecies,0,1e-3) # slope in marine marine survival
  alpha <- rnorm(Nspecies,0,1e-3) # slope in productivity
  tau.la <- rgamma(Nspecies,100,100) # variance in productivity
  tau.b <- rgamma(Nspecies,100,100) # variance in DD-survival slope
  rec.cv <- runif(Nspecies,0,0.5) # cv in recruitment
  stock.cv <- runif(Nspecies,0,0.5) # cv in recruitment
  beta <- 1e-9*matrix(rnorm(Nspecies*Nspecies,0,1),ncol=Nspecies,nrow=Nspecies)
  diag(beta) <- 0 # set the diagonals equal to 0
  b <- abs(rnorm(Nspecies,0,1e-3)) # DD-survival slope
  a <- abs(rnorm(Nspecies,10,1)) # intercept
  b_t <- matrix(1e-3*abs(rnorm(Nspecies,0,1)),Nyears+1,Nspecies)
  la_t <- matrix(log(abs(rnorm(Nspecies,10,1))),Nyears+1,Nspecies)
  marSurv_t <- matrix(rbeta(Nspecies,1,5),Nyears+1,Nspecies)
  rho <- runif(Nspecies,-1,1)
  return(list("sd.surv"=sd.surv,"marSurv"=marSurv,"marTrend"=marTrend,"alpha"=alpha,"tau.la"=tau.la,"tau.b"=tau.b,"rec.cv"=rec.cv,"stock.cv"=stock.cv,"beta"=beta,"b"=b,"a"=a,"b_t"=b_t,"la_t"=la_t,"marSurv_t"=marSurv_t,"rho"=rho))
}

gen_ini <- function()
{
  rec.cv <- runif(Nspecies,0,0.5) # cv in recruitment
  sys.sd <- abs(rnorm(Nspecies,100,10)) # cv in recruitment
  
  rho_Rec <- runif(Nspecies,0,1)
  #b_t <- matrix(1e-3*abs(rnorm(Nspecies,0,1)),Nyears+1,Nspecies)
  #la_t <- matrix(log(abs(rnorm(Nspecies,10,1))),Nyears+1,Nspecies)
  rho_a <- runif(Nspecies,0,1)
  rho_b <- runif(Nspecies,0,1)
  tau.b <- 1e-3*rgamma(Nspecies,100,100) # variance in DD-survival slope
  tau.a <- 1e-3*rgamma(Nspecies,100,100) # variance in alpha
  
  sd.surv <- abs(runif(Nspecies,0,1)) # variance in marine survival
  marSurv <- rbeta(Nspecies,1,5) # baseline marine survival
  marTrend <- rnorm(Nspecies,0,1e-3) # slope in marine marine survival
  
  return(list("rho_Rec"=rho_Rec,"sys.sd"=sys.sd,"rec.cv"=rec.cv,"rho_a"=rho_a,"rho_b"=rho_b,"tau.b"=tau.b,"tau.a"=tau.a,"sd.surv"=sd.surv,"marSurv"=marSurv,"marTrend"=marTrend))
}


ini <- lapply(1:nChains,FUN=function(x){gen_ini()})

init <- ini[[1]]

sapply(1:Nyears,FUN=function(y){sapply(1:Nspecies,FUN=function(x){ricker_lin(dat$stock[1,x],init$a[x],init$b_t[x],init$alpha[x],init$beta[x],0,y)})})


mon_names <- names(ini[[1]])

results <- run.jags(model="keogh_d1.jags",monitor=mon_names, data=dat, n.chains=nChains, method="rjags", inits=ini, plots=F,silent.jag=F, modules=c("glm","bugs"),sample=nSamps,adapt=nAdapt,burnin=nBurnin,thin=nThin,summarise=FALSE)


summary(results)
