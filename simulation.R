rm(list=ls(all=T))

get_beta <- function(mean,cv) #function that returns the alpha and beta shape parameters of a beta distribution, based on the mean and variation of a given beta distribution
{
  sd <- mean*cv
  alpha <- -((mean*(mean^2+sd^2-mean))/sd^2)
  beta <- alpha/mean-alpha
  return(list(alpha=alpha,beta=beta))
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

ricker <- function(stock,phiE0,R0,CR,beta,comp)
{
  a <- CR/phiE0
  b <- log(CR)/((R0+sum(beta*comp))*phiE0)
  recruits <- a*stock*exp(-b*stock)
  return(recruits)
}

Nspecies <- 5
Nyears <- 50

recCV <- rep(0.4,Nspecies)
phiE0 <- rep(1,Nspecies) # equilibrium eggs per recruit
R0 <- rep(1e4,Nspecies) # equilibrium recruitment
CR <- c(10,10,10,10,10) # recruitment compensation ratio for Ricker (Forrest et al. 2010 CJFAS)
mSurv <- rep(0.5,Nspecies)
mar.CV <- rep(0.01,Nspecies)

curve(ricker(x,phiE0[1],R0[1],CR[1],0,0),from=0,to=R0[1])


competition <- matrix(0,nrow=Nspecies,ncol=Nspecies,byrow=T,dimnames=list("victims"=1:Nspecies,"competitors"=1:Nspecies))
competition[1,2] <- 1
competition[2,1] <- -1
competition[4,3] <- 0
diag(competition) <- 0

a <- CR/phiE0
b <- log(CR)/(R0*phiE0)

competition <- competition*2e-1

stock <- matrix(NA,ncol=Nspecies,nrow=Nyears+1)
recruits <- matrix(NA,ncol=Nspecies,nrow=Nyears+1)

stock[1:2,] <- R0*mSurv
recruits[1:2,] <- sapply(1:Nspecies,function(x){ricker("stock"=stock[1,x],"phiE0"=phiE0[x],"R0"=R0[x],"CR"=CR[x],"beta"=rep(0,Nspecies),"comp"=rep(0,Nspecies))})

for(i in 3:(Nyears+1))
{
  shape <- get_beta(mSurv,mar.CV)
  mSurv.noise <- rbeta(Nspecies,shape$alpha,shape$beta)
  stock[i,] <- recruits[i-2,]*mSurv.noise
  mnRec <- sapply(1:Nspecies,function(x){ricker("stock"=stock[i,x],"phiE0"=phiE0[x],"R0"=R0[x],"CR"=CR[x],"beta"=competition[x,],"comp"=recruits[i-1,])})
  recruits[i,] <- pmax(0.1,rnorm(Nspecies,mean=mnRec,sd=mnRec*recCV))
}

species <- 2
c(R0[species],R0[species]+sum(competition[species,]*stock[i,]))

matplot(1:(Nyears+1),stock,type="l")
plot(stock[,1],recruits[,1])

pairs(recruits,upper.panel=panel.cor,lower.panel=panel.smooth)
pairs(stock,upper.panel=panel.cor,lower.panel=panel.smooth)

plot(recruits[1:Nyears,1],recruits[2:(Nyears+1),2])
