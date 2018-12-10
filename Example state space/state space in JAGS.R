# Author: Marie Auger-Methe
# Associated with manuscript: Auger-Methe et al. 2016 State-space models' dirty little secrets:
# even simple linear Gaussian models can have estimation problems. Scientific Reports


################################
# Simulation parameters
rho <- 0.7
sdObs <- 0.1
sdProSeq <- c(0.01,0.02,0.05,0.1,0.2,0.5,1)
sd0 <- 0.01 
n <- 100

# Fit model using JAGS code
library(rjags)

parameters <- c("sdObs", "sdPro", "rho", "x")  # The parameter(s) to be monitored
adaptSteps <- 50000 # Number of steps to "tune" the samplers
burnInSteps <- 50000 # Number of steps to "burn-in" the samplers
nChains <- 2 # Number of chains to run
nPerChain <- 50000 # Steps per chain
thinSteps <- 500 # Number of steps to "thin" (1=keep every step)


#############################
# Make a loop to look at the distribution of the parameter estimates
nSim <- 200

# Matrices that will keep track of information
parEst <- matrix(nrow=nSim*length(sdProSeq), ncol=4)
colnames(parEst) <- c("sdObsEst", "sdProEst", "rhoEst", "sdProSim")

rmseEst <- matrix(nrow=nSim*length(sdProSeq), ncol=2)
colnames(rmseEst) <- c("estPar", "simPar")

grEst <- matrix(nrow=nSim*length(sdProSeq), ncol=3)
colnames(grEst) <- c("sdObsGR", "sdProGR", "rhoGR")

set.seed(567)

for(k in seq_along(sdProSeq)){
  for(j in 1:nSim){
    indexSim <- (k-1)*nSim+j
    parEst[indexSim,4] <- sdProSeq[k]
    # Simulate the data
    x <- vector(length=n+1, mode="numeric") 
    # Assumes that x_1=0
    for(i in 2:(n+1)){
      x[i] <- rho*x[i-1] + rnorm(1,0,sdProSeq[k])
    }
    y <- x[-1] + rnorm(n, 0, sdObs)
    
    initFx <- function(){
      # Sampling from distribution close to the prior distribution
      # But a bit more centered to help convergence
      sdObs <- abs(rnorm(1, 0, 1)) # dnorm(0,0.0001)T(0,)
      sdPro <- abs(rnorm(1, 0, 1)) # dnorm(0,0.0001)T(0,)
      # for rho exactly the prior
      rho <-  runif(1, 0, 1)
      list(sdObs = sdObs, sdPro = sdPro, rho=rho) 
    }
    initVal <- list(initFx(),initFx())
    
    dataList <- list(y = y, N=length(y))
    
    # Create, initialize, and adapt the model:
    jagsModel = jags.model( "1DSSM.jags" , data=dataList , inits=initVal,
                            n.chains=nChains , n.adapt=adaptSteps )
    # Burn-in:
    update(jagsModel, n.iter=burnInSteps)
    # The saved MCMC chain:
    codaSamples <- coda.samples(jagsModel, variable.names=parameters, 
                                 n.iter=nPerChain, thin=thinSteps)
    
    # Parameter estimates
    parEst[indexSim,1:3] <- tryCatch(summary(codaSamples[,c(2,3,1)])[[1]][,1],
                                     error=function(e) rep(NA,3)) 

    # True state estimates
    xEst <- summary(codaSamples[,1:(n+1)+3])[[1]][,1]    
    rmseEst[indexSim,1] <- tryCatch(sqrt(sum((xEst - x)[-1]^2)/(n)),
                                    error=function(e) NA)
    
    # Convergence
    grEst[indexSim,] <- tryCatch(gelman.diag(codaSamples[,c(2,3,1)])[[1]][,1],
                                 error=function(e) NA)
    
    ##
    # Getting RMSE Doing the same with the simulated values
    dataList <- list(y = y, N=length(y), sdObs=sdObs, sdPro=sdProSeq[k], rho=rho)
    
    # Create, initialize, and adapt the model:
    jagsModelSimVal <- jags.model( "1DSSM_SimVal.jags", data=dataList,
                            n.chains=nChains, n.adapt=adaptSteps)
    # Burn-in:
    update(jagsModelSimVal, n.iter=burnInSteps)
    # The saved MCMC chain:
    codaSamplesSimVal <- coda.samples(jagsModelSimVal, variable.names=parameters, 
                                n.iter=nPerChain, thin=thinSteps)
    # True state estimates
    xEstSimVal <- summary(codaSamplesSimVal[,1:(n+1)+3])[[1]][,1]    
    rmseEst[indexSim,2] <- tryCatch(sqrt(sum((xEstSimVal - x)[-1]^2)/(n)),
                                    error=function(e) NA)
  }
}

# Figure to look at parameter estimates
quartz(width=6, height=6)
layout(matrix(1:(4*length(sdProSeq)),nrow=length(sdProSeq),byrow=TRUE))
par(mar=c(1.2,1.4,0.3,0.5), mgp=c(0.8,0.3,0), tck=-0.03, las=1, oma=c(1.5,2,0,0.3))
yleg <- 0.1
xleg <- -0.7
placl <- "topleft"
lt <- c(LETTERS, "AA", "BB")
for(i in 1:length(sdProSeq)){
  sdProIndex <- parEst[,4] == sdProSeq[i] 
  
  hh <- hist(parEst[sdProIndex,1],
             breaks=seq(0,max(parEst[,1],na.rm=TRUE), length.out=50), plot=FALSE)
  plot(hh, ylim =c(0,round(max(hh$counts)*1.24)),
       xlab="", main="", ylab="",
       border=FALSE, col="darkgrey")
  abline(v=sdObs)
  ys <- seq(0.982,0.04,-0.157)
  title(ylab="Frequency", line=0.4, cex.lab=1.2, outer=TRUE, 
        adj=ys[i])
  rat <- sdObs/sdProSeq[i]
  legend("top",
         legend = substitute(sigma[epsilon] == rat*~sigma[eta],list(rat=rat)), bty="n")
  box()
  legend(placl, lt[i*4-3], bty="n", x.intersp=xleg, y.intersp=yleg)
  
  hist(parEst[sdProIndex,3],
       breaks=seq(-1,1,length.out=50), xlab="", main="", 
       ylab="",
       border=FALSE, col="darkgrey")
  abline(v=rho)
  box()
  legend(placl, lt[i*4-2], bty="n", x.intersp=xleg, y.intersp=yleg)
  
  hist(parEst[sdProIndex,2],
       breaks=seq(0,max(parEst[sdProIndex,2],na.rm=TRUE),length.out=50), xlab="", main="",
       ylab="",
       border=FALSE, col="darkgrey")
  abline(v=unique(parEst[sdProIndex,4]), col="red")
  box()
  legend(placl, lt[i*4-1], bty="n", x.intersp=xleg, y.intersp=yleg)
  
  # RMSE
  hEst <- hist(rmseEst[sdProIndex,1],
               breaks=seq(min(rmseEst[sdProIndex,1:2], na.rm=TRUE),
                          max(rmseEst[sdProIndex,1:2], na.rm=TRUE), length.out=50), 
               plot=FALSE)
  hTrue <- hist(rmseEst[sdProIndex,2],
                breaks=seq(min(rmseEst[sdProIndex,1:2], na.rm=TRUE),
                           max(rmseEst[sdProIndex,1:2], na.rm=TRUE), length.out=50),
                plot=FALSE)
  plot(hEst, 
       ylim=c(0,round(max(hEst$counts,hTrue$counts)*1.21)),
       xlab="", main="", ylab="", border=FALSE, col=grey(0.5))
  plot(hTrue, border=FALSE, col=rgb(0,0,1,0.5), add=TRUE)
  box()
  legend(placl, lt[i*4], bty="n", x.intersp=xleg, y.intersp=yleg)
  
}
title(xlab=substitute(widehat(sigma)[epsilon]), line=0.3, cex.lab=1.2,outer=TRUE, 
      adj=0.13)
title(xlab=substitute(widehat(rho)), line=0.3, cex.lab=1.2,outer=TRUE, 
      adj=0.38)
title(xlab=substitute(widehat(sigma)[eta]), line=0.3, cex.lab=1.2,outer=TRUE, 
      adj=0.63)
title(xlab="RMSE", line=0.1, cex.lab=1.2,outer=TRUE, 
      adj=0.92)

warnings()

# Check whether the chains converged
# First check whether they are NA
any(is.na(grEst))
notConvIndex <- apply(grEst, 1, function(x) any(x>1.1))
grEst[notConvIndex,]
sum(notConvIndex)
sum(notConvIndex)/nrow(grEst)*100
cbind(sdProSeq, hist(parEst[notConvIndex,4], breaks=c(0,sdProSeq), plot=FALSE)$counts)
