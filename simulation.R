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

ricker_lin <- function(stock,a,b,alpha,beta,comp,t)
{
  la_t <- log(a)+alpha*t
  b_t <- b+sum(beta*comp)
  recruits <- exp(la_t)*stock*exp(-b_t*stock)
  return(recruits)
}

Nspecies <- 5
Nyears <- 40

recCV <- rep(0.4,Nspecies)
phiE0 <- rep(1,Nspecies) # equilibrium eggs per recruit
R0 <- rep(1e4,Nspecies) # equilibrium recruitment
CR <- c(10,10,10,10,10) # recruitment compensation ratio for Ricker (Forrest et al. 2010 CJFAS)
mSurv <- rep(0.3,Nspecies)
mInits <- 1/(1+exp(-log(mSurv/(1-mSurv))))
mar.CV <- rep(0.1,Nspecies)

a <- CR/phiE0
b <- log(CR)/(R0*phiE0)

alpha <- c(0,0,0,-5e-2,-1e-1) # slope in alpha through time
marTrend <- c(0,0,-5e-3,0,0)

curve(ricker(x,phiE0[1],R0[1],CR[1],0,0),from=0,to=R0[1])


competition <- matrix(0,nrow=Nspecies,ncol=Nspecies,byrow=T,dimnames=list("victims"=1:Nspecies,"competitors"=1:Nspecies))
competition[1,2] <- 1
competition[2,1] <- -1
competition[4,3] <- 0
diag(competition) <- 0

competition <- competition
beta <- competition*1e-8 # scale competition coefficients for b

stock <- matrix(NA,ncol=Nspecies,nrow=Nyears+1)
recruits <- matrix(NA,ncol=Nspecies,nrow=Nyears+1)

stock[1,] <- R0*mSurv
recruits[1,] <- sapply(1:Nspecies,function(x){ricker_lin("stock"=stock[1,x],"a"=a[x],"b"=b[x],"alpha"=0,"beta"=rep(0,Nspecies),"comp"=rep(1,Nspecies),"t"=1)})

for(i in 2:(Nyears+1))
{
  # do logistic regression on marine survival
  shape <- get_beta(pmax(0,pmin(mSurv+marTrend*i,1)),mar.CV)
  mSurv.noise <- rbeta(Nspecies,shape$alpha,shape$beta)
  stock[i,] <- recruits[i-1,]*mSurv.noise
  #mnRec <- sapply(1:Nspecies,function(x){ricker("stock"=stock[i,x],"phiE0"=phiE0[x],"R0"=R0[x],"CR"=CR[x],"beta"=competition[x,],"comp"=recruits[i-1,])})
  
  mnRec <- sapply(1:Nspecies,function(x){ricker_lin("stock"=stock[1,x],"a"=a[x],"b"=b[x],"alpha"=alpha[x],"beta"=beta[x,],"comp"=recruits[i-1,],"t"=i)})
  
  recruits[i,] <- pmax(0.1,rnorm(Nspecies,mean=mnRec,sd=mnRec*recCV))
}

species <- 2
c(R0[species],R0[species]+sum(competition[species,]*stock[i,]))

matplot(1:(Nyears+1),stock,type="l",xlab="time")
plot(stock[1:Nyears,1],recruits[2:(Nyears+1),1],xlab="Stock",ylab="Recruits")

pairs(recruits,upper.panel=panel.cor,lower.panel=panel.smooth)
pairs(stock,upper.panel=panel.cor,lower.panel=panel.smooth)

layout(matrix(1:2,ncol=2))
par(mar=c(5,4,1,1))
plot(recruits[1:Nyears,1],recruits[2:(Nyears+1),2],xlab="Victims(t)",ylab="Competitors(t+1)",pch=21,bg="dodgerblue")
plot(recruits[1:Nyears,2],recruits[2:(Nyears+1),1],ylab="Competitors(t)",xlab="Victims(t+1)",pch=21,bg="dodgerblue")

layout(matrix(1))
plot(stock[2:(Nyears+1),3]/recruits[1:Nyears,3],ylab="Marine survival",pch=21,bg="grey50")

# model degrees of freedom
df <- sum(length(a),length(b),length(alpha),length(beta)-length(diag(beta)),length(),Nspecies) # number of parameters including variance terms
length(recruits)/df

stockRec <- list("stock"=stock,"recruits"=recruits,"Nyears"=nrow(recruits),"Nspecies"=ncol(recruits))

saveRDS(stockRec,"stockRec.rds")
