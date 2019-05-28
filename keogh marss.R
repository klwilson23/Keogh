# RShowDoc("Chapter_DFA.R",package="MARSS")
# C:\Users\kylel\R\win-library\3.5\MARSS\doc
getDFAfits <- function(MLEobj, alpha=0.05, covariates=NULL) {
  fits <- list()
  Ey <- MARSShatyt(MLEobj) # for var() calcs
  ZZ <- coef(MLEobj, type="matrix")$Z # estimated Z
  nn <- nrow(ZZ) # number of obs ts
  mm <- ncol(ZZ) # number of factors/states
  TT <- ncol(Ey$ytT)  # number of time steps
  ## check for covars
  if(!is.null(covariates)) {
    DD <- coef(MLEobj, type="matrix")$D
    cov_eff <- DD %*% covariates
  } else {
    cov_eff <- matrix(0, nn, TT)
  }
  ## model expectation
  fits$ex <- ZZ %*% MLEobj$states + cov_eff
  ## Var in model fits
  VtT <- MARSSkfss(MLEobj)$VtT
  VV <- NULL
  for(tt in 1:TT) {
    RZVZ <- coef(MLEobj, type="matrix")$R - ZZ%*%VtT[,,tt]%*%t(ZZ)
    SS <- Ey$yxtT[,,tt] - Ey$ytT[,tt,drop=FALSE] %*% t(MLEobj$states[,tt,drop=FALSE])
    VV <- cbind(VV,diag(RZVZ + SS%*%t(ZZ) + ZZ%*%t(SS)))
  }
  SE <- sqrt(VV)
  ## upper & lower (1-alpha)% CI
  fits$up <- qnorm(1-alpha/2)*SE + fits$ex
  fits$lo <- qnorm(alpha/2)*SE + fits$ex
  return(fits)
}


library(MARSS)
library(broom)
library(ggplot2)
SRdata <- read.csv("Keogh_StockRecruitment.csv",header=TRUE)
SRdata <- subset(SRdata,Year>=1976 & Year<=2013)
#SRdata <- SRdata[,-grep("ch_",colnames(SRdata))]
SRdata[,-1] <- SRdata[,-1]
sdScale <- attr(scale(SRdata[,-1],center=FALSE,scale=TRUE),"scaled:scale")
SRdata[,-1] <- scale(SRdata[,-1],center=FALSE,scale=TRUE)

environment <- read.csv("Data/keogh environmental covariates.csv",header=TRUE)
environment <- environment[environment$year>=min(SRdata$Year) & environment$year<=max(SRdata$Year),]
sdCovars <- attr(scale(environment,center=TRUE,scale=TRUE),"scaled:scale")
covarScale <- scale(environment,center=TRUE,scale=TRUE)
covarScale[is.na(covarScale)] <- 0

newDat <- SRdata[,-c(5,13)]

# run a DFA on the adult data
adultDat <- newDat[,c(2,4,5,7,9,11)]
adultDat <- t(as.matrix(adultDat))
colnames(adultDat) <- newDat$Year
ns <- nrow(adultDat)
B <- "diagonal and equal"
Q <- "unconstrained"
R <- diag(0.01,ns)
U <- "zero"
A <- "unequal"
x0 <- "zero"
mod.list = list(B=B, Q=Q, R=R, U=U, x0=x0, A=A,tinitx = 1)

# DFA model
nTrends <- 2
B <- matrix(list(0), nTrends, nTrends)
diag(B) <- paste("b",1:nTrends,sep="")

Q <- diag(1, nTrends)
R <- "diagonal and unequal"
#R <- diag(0.01,ns)
U <- "zero"
x0 <- "zero"
Z <- matrix(list(0), ns, nTrends)
Z[1:(ns * nTrends)] <- sapply(1:nTrends,function(x){paste0("z",x,1:ns)})
Z[upper.tri(Z)] <- 0
A <- "unequal"
#A <- "zero"
mod.list.dfa = list(B = B, Z = Z, Q = Q, R = R, U = U, A = A, 
                    x0 = x0)

m <- apply(adultDat, 1, mean, na.rm=TRUE)
fit <- MARSS(adultDat, model=mod.list, control=list(minit=200,maxit=50000+200), inits=list(A=matrix(m,ns,1)))
fit$AICc
Z.est = coef(fit, type="matrix")$Z
H.inv = 1
if(ncol(Z.est)>1) H.inv = varimax(coef(fit, type="matrix")$Z)$rotmat
# rotate factor loadings
Z.rot = Z.est %*% H.inv
# rotate trends
trends.rot = solve(H.inv) %*% fit$states

layout(1)
matplot(t(trends.rot),type="l")
matplot(t(fit$states),type="l")

fit.b = getDFAfits(fit)

d <- augment(fit, interval = "confidence")
d$Year <- d$t + 1975
d$Species <- d$.rownames

fit$states

adultFits <- matrix(sapply(1:nrow(d),function(x){d$.fitted[x]*sdScale[names(sdScale)%in%d$Species[x]]}),nrow=ncol(adultDat),ncol=nrow(adultDat),byrow=FALSE)

SRdata <- read.csv("Keogh_StockRecruitment.csv",header=TRUE)
SRdata <- subset(SRdata,Year>=1976 & Year<=2013)

ct_Adults.hat <- pmax(1e-4,ifelse(is.na(SRdata$ct_Adults),adultFits[,3],SRdata$ct_Adults))
dv_Adults.hat <- pmax(1e-4,ifelse(is.na(SRdata$dv_Adults),adultFits[,5],SRdata$dv_Adults))

cbind(SRdata$ct_Adults,adultFits[,3])

yy <- data.frame("Year"=newDat$Year,newDat[,c(2,4,5,7,9,11)])
yy <- reshape2::melt(yy, id.vars=c("Year"))
colnames(yy) <- c("Year","Species","Abundance")
p <- ggplot(data = d) + 
  geom_line(aes(Year, .fitted)) +
  geom_ribbon(aes(x=Year, ymin=.conf.low, ymax=.conf.up), linetype=2, alpha=0.5)
p <- p + geom_point(data=yy, mapping = aes(x=Year, y=Abundance))
p + facet_wrap(~Species) + xlab("") + ylab("Abundance")

layout(matrix(1:6,nrow=3,ncol=2))
apply(residuals(fit)$state.residuals[, 1:30], 1, acf)
mtext("State Residuals ACF", outer = TRUE, side = 3)

# run a DLM on stock-recruitment for steelhead only
SRdata <- read.csv("Keogh_StockRecruitment.csv",header=TRUE)
SRdata <- subset(SRdata,Year>=1976 & Year<=2013)
steelhead <- SRdata[,c("Year","sh_Adults","sh_Smolts")]
dolly <- SRdata[,c("Year","dv_Adults","dv_Smolts")]
cutty <- SRdata[,c("Year","ct_Adults","ct_Smolts")]


# including process and observation error
# adult and preciptation covariates only affect observation model
# NPGO affects process model
A <- U <- "zero"; Z <- "identity"
B <- "diagonal and unequal"
Q <- "equalvarcov"
D <- "unconstrained"
d <- t(matrix(c(steelhead$sh_Adults,covarScale[,"precip_1"]),ncol=2))
row.names(d) <- c("adults","precip")
C <- "unconstrained"
c <- t(matrix(covarScale[,"npgo"]))
row.names(c) <- c("npgo")
R <- "diagonal and unequal"
x0 <- "unequal"
tinitx <- 0
ln_RS <- steelhead$ln_RS <- log(steelhead$sh_Smolts/steelhead$sh_Adults)
model.list <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,D=D,d=d,C=C,c=c,x0=x0,tinitx=tinitx)
keoghSR <- MARSS(ln_RS, model=model.list)

keoghSR$model
keoghSR$states
mean(colSums(keoghSR$states))
MARSSparamCIs(keoghSR)
keoghSR$coef

library(broom)
library(ggplot2)
d <- augment(keoghSR, interval="confidence")
d$Year <- d$t + 1975
p <- ggplot(data = d) + 
  geom_line(aes(Year, .fitted)) +
  geom_ribbon(aes(x=Year, ymin=.conf.low, ymax=.conf.up), linetype=2, alpha=0.5)
p <- p + geom_point(data=steelhead, mapping = aes(x=Year, y=ln_RS))
p + xlab("") + ylab("ln_RS")


# including process and observation error
# adult and precipitation covariates only affect observation model
# NPGO affects process model
# time-varying affect on adults
C <- c <- A <- U <- "zero"
Nyears <- nrow(steelhead)
Nspecies <- 1
m <- 1
Z <- array(NA,c(1,m,Nyears),dimnames=list("Species"="steelhead","beta"=1:m,"year"=steelhead$Year))
Z[1,1,] <- steelhead$sh_Adults # time-varying beta
#Z[1,1,] <- rep(1,Nyears) # time-varying alpha
#Z[1,1,] <- covarScale[,"precip_1"] # time-varying precipitation

D <- "unconstrained"
d <- rbind(rep(1,Nyears),covarScale[,"precip_1"])
row.names(d) <- c("alpha","precip")

#C <- "unconstrained" # NPGO affects observation (time-varying component)
#c <- t(matrix(covarScale[,"npgo"]))
#row.names(c) <- c("npgo")

Q <- matrix(list(0),m,m)         ## 2x2; all 0 for now
diag(Q) <- c("q.adults") ## 2x2; diag = (q1,q2)
B <- diag(m)

R <- "diagonal and unequal"
x0 <- "unequal"
tinitx <- 0
ln_RS <- log(steelhead$sh_Smolts/steelhead$sh_Adults)
model.list <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,D=D,d=d,C=C,c=c,x0=x0,tinitx=tinitx)
model.listDLM <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,d=d,C=C,c=c,x0=x0,tinitx=tinitx)
keoghDLM <- MARSS(ln_RS, model=model.listDLM,control=list(maxit=2000,conv.test.slope.tol=0.1))

#From MARSSinfo()
#"This is not an error but rather an fyi.  Q or R is getting very small.  Because control$allow.degen=TRUE, the code is trying to set Q or R to 0, but in fact, the MLE Q or R is not 0 so setting to 0 is being blocked (correctly).  The code is warning you about this because when Q or R gets really small, you can have numerical problems in the algorithm.  Have you standardized the variances in your data and covariates (explanatory variables)?  In some types of models, that kind of mismatch can drive Q or R towards 0.  This is correct behavior, but you may want to standardize your data so that the variability is on similar scales.\n"

keoghDLM$model
keoghDLM$states
mean(colSums(keoghDLM$states))
MARSSparamCIs(keoghDLM)
keoghDLM$coef

fitDLM <- augment(keoghDLM, interval="confidence")
fitDLM$Year <- fitDLM$t + 1975
p <- ggplot(data = fitDLM) + 
  geom_line(aes(Year, .fitted)) +
  geom_ribbon(aes(x=Year, ymin=.conf.low, ymax=.conf.up), linetype=2, alpha=0.5)
p <- p + geom_point(data=steelhead, mapping = aes(x=Year, y=ln_RS))
p + xlab("") + ylab("ln_RS")

layout(1)
plot(sh_Smolts~sh_Adults,data=steelhead,pch=21,bg=ifelse(steelhead$Year>=1991,"orange","dodgerblue"))
curve(exp(keoghDLM$coef["D.(Y1,alpha)"]+keoghDLM$coef["D.(Y1,precip)"]*0)*x*exp(colSums(keoghDLM$states)[1]*x),add=TRUE,lwd=2,col="dodgerblue")
curve(exp(keoghDLM$coef["D.(Y1,alpha)"]+keoghDLM$coef["D.(Y1,precip)"]*0)*x*exp(colSums(keoghDLM$states)[Nyears]*x),add=TRUE,lwd=2,col="orange",lty=2)


ylabs <- c(expression("adults"~beta[t]))#, expression("pdo 2"~beta[t]),expression("pdo 3"~beta[t]))
colr <- c("darkgreen","dodgerblue","orange")
par(mfrow=c(m,1), mar=c(4,4,0.1,0), oma=c(0,0,2,0.5))
for(i in 1:m) {
  mn <- keoghDLM$states[i,]
  se <- keoghDLM$states.se[i,]
  plot(steelhead$Year,mn,xlab="",ylab=ylabs[i],bty="n",xaxt="n",type="n",
       ylim=c(min(mn-2*se),max(mn+2*se)))
  lines(steelhead$Year, rep(0,Nyears), lty="dashed")
  lines(steelhead$Year, mn, col=colr[i], lwd=3)
  lines(steelhead$Year, mn+2*se, col=colr[i])
  lines(steelhead$Year, mn-2*se, col=colr[i])
}
axis(1,at=seq(min(steelhead$Year),max(steelhead$Year),5))
mtext("Brood year", 1, line=3)
abline(v=1991)

# including process and observation error
# precipitation covariates only affect observation model
# time-varying beta & alpha

C <- c <- A <- U <- "zero"
Nyears <- nrow(steelhead)
Nspecies <- 1
m <- 2
Z <- array(NA,c(1,m,Nyears))
Z[1,1,] <- rep(1,Nyears) # time-varying alpha
Z[1,2,] <- steelhead$sh_Adults # time-varying beta
#Z[1,3,] <- covarScale[,"precip_1"] # time-varying precipitation

D <- "unconstrained"
d <- rbind(covarScale[,"precip_1"])
row.names(d) <- c("precip")

#C <- "unconstrained" # NPGO affects observation (time-varying component)
#c <- t(matrix(covarScale[,"npgo"]))
#row.names(c) <- c("npgo")

Q <- matrix(list(0),m,m)         ## 2x2; all 0 for now
diag(Q) <- c("q.alpha","q.adults") ## 2x2; diag = (q1,q2)
B <- diag(m)

R <- "diagonal and unequal"
x0 <- "unequal"
tinitx <- 0
inits.list <- list(x0=matrix(rep(0,m), nrow=m))
model.listDLMall <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,d=d,C=C,c=c)
keoghDLMall <- MARSS(ln_RS, model=model.listDLMall,control=list(maxit=2000,conv.test.slope.tol=0.1),inits=inits.list)


keoghDLMall$model
keoghDLMall$states
mean(colSums(keoghDLMall$states))
MARSSparamCIs(keoghDLMall)
keoghDLMall$coef

fitDLMall <- augment(keoghDLMall, interval="confidence")
fitDLMall$Year <- fitDLMall$t + 1975
p <- ggplot(data = fitDLMall) + 
  geom_line(aes(Year, .fitted)) +
  geom_ribbon(aes(x=Year, ymin=.conf.low, ymax=.conf.up), linetype=2, alpha=0.5)
p <- p + geom_point(data=steelhead, mapping = aes(x=Year, y=ln_RS))
p + xlab("") + ylab("ln_RS")

ylabs <- c(expression(alpha[t]),expression("adults"~beta[t]),expression(K[t]))#, expression("pdo 2"~beta[t]),expression("pdo 3"~beta[t]))
colr <- c("darkgreen","dodgerblue","orange")

layout(matrix(1:3,nrow=3,ncol=1))
par(mar=c(5,5,1,1))
for(i in 1:m) {
  mn <- keoghDLMall$states[i,]
  se <- keoghDLMall$states.se[i,]
  plot(steelhead$Year,mn,xlab="",ylab=ylabs[i],bty="n",xaxt="n",type="n",
       ylim=c(min(mn-2*se),max(mn+2*se)),cex.lab=2)
  lines(steelhead$Year, rep(0,Nyears), lty="dashed")
  lines(steelhead$Year, mn, col=colr[i], lwd=3)
  lines(steelhead$Year, mn+2*se, col=colr[i])
  lines(steelhead$Year, mn-2*se, col=colr[i])
  abline(v=1991)
}

mn <- -log(keoghDLMall$states[1,])/keoghDLMall$states[2,]
plot(steelhead$Year,mn,xlab="",ylab=ylabs[3],bty="n",xaxt="n",type="n",cex.lab=2)
lines(steelhead$Year, rep(0,Nyears), lty="dashed")
lines(steelhead$Year, mn, col=colr[3], lwd=3)
axis(1,at=seq(min(steelhead$Year),max(steelhead$Year),5),cex=2)
mtext("Brood year", 1, line=3,cex=2)
abline(v=1991)

AIC(keoghDLM,keoghDLMall)
keoghDLM$AICc
keoghDLMall$AICc

# dolly varden:
ln_RS_dv <- log(dolly$dv_Smolts/dolly$dv_Adults)

C <- c <- A <- U <- "zero"
Nyears <- nrow(steelhead)
Nspecies <- 1
m <- 2
Z <- array(NA,c(1,m,Nyears))
Z[1,1,] <- rep(1,Nyears) # time-varying alpha
Z[1,2,] <- dv_Adults.hat # time-varying beta
#Z[1,3,] <- covarScale[,"precip_1"] # time-varying precipitation

D <- "unconstrained"
d <- rbind(covarScale[,"precip_1"])
row.names(d) <- c("precip")

#C <- "unconstrained" # NPGO affects observation (time-varying component)
#c <- t(matrix(covarScale[,"npgo"]))
#row.names(c) <- c("npgo")

Q <- matrix(list(0),m,m)         ## 2x2; all 0 for now
diag(Q) <- c("q.alpha","q.adults") ## 2x2; diag = (q1,q2)
B <- diag(m)

R <- "diagonal and unequal"
x0 <- "unequal"
tinitx <- 0
inits.list <- list(x0=matrix(rep(0,m), nrow=m))
model.listDLMdv <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,d=d,C=C,c=c)
keoghDLMdv <- MARSS(ln_RS_dv, model=model.listDLMdv,control=list(maxit=2000,conv.test.slope.tol=0.1),inits=inits.list)


keoghDLMdv$model
keoghDLMdv$states
MARSSparamCIs(keoghDLMdv)
keoghDLMdv$coef

fitDLMdv <- augment(keoghDLMdv, interval="confidence")
fitDLMdv$Year <- fitDLMdv$t + 1975
p <- ggplot(data = fitDLMdv) + 
  geom_line(aes(Year, .fitted)) +
  geom_ribbon(aes(x=Year, ymin=.conf.low, ymax=.conf.up), linetype=2, alpha=0.5)
p <- p + geom_point(data=steelhead, mapping = aes(x=Year, y=ln_RS_dv))
p + xlab("") + ylab("ln_RS_dv")

ylabs <- c(expression(alpha[t]),expression("adults"~beta[t]),expression(K[t]))#, expression("pdo 2"~beta[t]),expression("pdo 3"~beta[t]))
colr <- c("darkgreen","dodgerblue","orange")

layout(matrix(1:3,nrow=3,ncol=1))
par(mar=c(5,5,1,1))
for(i in 1:m) {
  mn <- keoghDLMdv$states[i,]
  se <- keoghDLMdv$states.se[i,]
  plot(steelhead$Year,mn,xlab="",ylab=ylabs[i],bty="n",xaxt="n",type="n",
       ylim=c(min(mn-2*se),max(mn+2*se)),cex.lab=2)
  lines(steelhead$Year, rep(0,Nyears), lty="dashed")
  lines(steelhead$Year, mn, col=colr[i], lwd=3)
  lines(steelhead$Year, mn+2*se, col=colr[i])
  lines(steelhead$Year, mn-2*se, col=colr[i])
  abline(v=1991)
}

mn <- -log(keoghDLMdv$states[1,])/keoghDLMdv$states[2,]
plot(steelhead$Year,mn,xlab="",ylab=ylabs[3],bty="n",xaxt="n",type="n",cex.lab=2)
lines(steelhead$Year, rep(0,Nyears), lty="dashed")
lines(steelhead$Year, mn, col=colr[3], lwd=3)
axis(1,at=seq(min(steelhead$Year),max(steelhead$Year),5),cex=2)
mtext("Brood year", 1, line=3,cex=2)
abline(v=1991)

# cutthroat:
ln_RS_ct <- log(cutty$ct_Smolts/cutty$ct_Adults)

C <- c <- A <- U <- "zero"
Nyears <- nrow(steelhead)
Nspecies <- 1
m <- 2
Z <- array(NA,c(1,m,Nyears))
Z[1,1,] <- rep(1,Nyears) # time-varying alpha
Z[1,2,] <- ct_Adults.hat # time-varying beta
#Z[1,3,] <- covarScale[,"precip_1"] # time-varying precipitation

D <- "unconstrained"
d <- rbind(covarScale[,"precip_1"])
row.names(d) <- c("precip")

#C <- "unconstrained" # NPGO affects observation (time-varying component)
#c <- t(matrix(covarScale[,"npgo"]))
#row.names(c) <- c("npgo")

Q <- matrix(list(0),m,m)         ## 2x2; all 0 for now
diag(Q) <- c("q.alpha","q.adults") ## 2x2; diag = (q1,q2)
B <- diag(m)

R <- "diagonal and unequal"
x0 <- "unequal"
tinitx <- 0
inits.list <- list(x0=matrix(rep(0,m), nrow=m))
model.listDLMct <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,d=d,C=C,c=c)
keoghDLMct <- MARSS(ln_RS_ct, model=model.listDLMct,control=list(maxit=2000,conv.test.slope.tol=0.1),inits=inits.list)


keoghDLMct$model
keoghDLMct$states
MARSSparamCIs(keoghDLMct)
keoghDLMct$coef

fitDLMct <- augment(keoghDLMct, interval="confidence")
fitDLMct$Year <- fitDLMct$t + 1975
p <- ggplot(data = fitDLMct) + 
  geom_line(aes(Year, .fitted)) +
  geom_ribbon(aes(x=Year, ymin=.conf.low, ymax=.conf.up), linetype=2, alpha=0.5)
p <- p + geom_point(data=steelhead, mapping = aes(x=Year, y=ln_RS_ct))
p + xlab("") + ylab("ln_RS_ct")

ylabs <- c(expression(alpha[t]),expression("adults"~beta[t]),expression(K[t]))#, expression("pdo 2"~beta[t]),expression("pdo 3"~beta[t]))
colr <- c("darkgreen","dodgerblue","orange")

layout(matrix(1:3,nrow=3,ncol=1))
par(mar=c(5,5,1,1))
for(i in 1:m) {
  mn <- keoghDLMct$states[i,]
  se <- keoghDLMct$states.se[i,]
  plot(steelhead$Year,mn,xlab="",ylab=ylabs[i],bty="n",xaxt="n",type="n",
       ylim=c(min(mn-2*se),max(mn+2*se)),cex.lab=2)
  lines(steelhead$Year, rep(0,Nyears), lty="dashed")
  lines(steelhead$Year, mn, col=colr[i], lwd=3)
  lines(steelhead$Year, mn+2*se, col=colr[i])
  lines(steelhead$Year, mn-2*se, col=colr[i])
  abline(v=1991)
}

mn <- -log(keoghDLMct$states[1,])/keoghDLMct$states[2,]
plot(steelhead$Year,mn,xlab="",ylab=ylabs[3],bty="n",xaxt="n",type="n",cex.lab=2)
lines(steelhead$Year, rep(0,Nyears), lty="dashed")
lines(steelhead$Year, mn, col=colr[3], lwd=3)
axis(1,at=seq(min(steelhead$Year),max(steelhead$Year),5),cex=2)
mtext("Brood year", 1, line=3,cex=2)
abline(v=1991)

# multiple species: dolly varden & cutthroat trout
# including process and observation error
# precipitation covariates only affect observation model
# time-varying beta & alpha
# run a DLM on stock-recruitment for steelhead only

C <- c <- A <- U <- "zero"
Nyears <- nrow(steelhead)
Nspecies <- 3
m <- 2
Z <- array(0,c(Nspecies,m*Nspecies,Nyears))
Z[1,1,] <- rep(1,Nyears) # time-varying alpha
Z[1,2,] <- steelhead$sh_Adults # time-varying beta
Z[2,3,] <- rep(1,Nyears) # time-varying alpha
Z[2,4,] <- dv_Adults.hat # time-varying beta
Z[3,5,] <- rep(1,Nyears) # time-varying alpha
Z[3,6,] <- ct_Adults.hat # time-varying beta
#Z[1,3,] <- covarScale[,"precip_1"] # time-varying precipitation

D <- "unconstrained"
d <- rbind(covarScale[,"precip_1"])
row.names(d) <- c("precip")

#C <- "unconstrained" # NPGO affects observation (time-varying component)
#c <- t(matrix(covarScale[,"npgo"]))
#row.names(c) <- c("npgo")

Q <- matrix(list(0),m*Nspecies,m*Nspecies)        ## 2x2; all 0 for now
#diag(Q) <- rep(c("q.alpha","q.beta"),Nspecies)
diag(Q) <- c("q.sh_alpha","q.sh_adults","q.dv_alpha","q.dv_adults","q.ct_alpha","q.ct_adults") ## 2x2; diag = (q1,q2)
Q[1,3] <- "q.corr_alpha1"
Q[3,1] <- "q.corr_alpha1"
Q[3,5] <- "q.corr_alpha2"
Q[5,3] <- "q.corr_alpha2"
Q[1,5] <- "q.corr_alpha3"
Q[5,1] <- "q.corr_alpha3"

Q[2,4] <- "q.corr_beta1"
Q[4,2] <- "q.corr_beta1"
Q[4,6] <- "q.corr_beta2"
Q[6,4] <- "q.corr_beta2"
Q[2,6] <- "q.corr_beta3"
Q[6,2] <- "q.corr_beta3"

B <- diag(m*Nspecies)

R <- "unconstrained"
#R <- diag(0.01,Nspecies)
x0 <- "zero"
tinitx <- 0
model.listDLMspecies <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,d=d,C=C,c=c,x0=x0,tinitx=0)
ln_RS_sh <- log(steelhead$sh_Smolts/steelhead$sh_Adults)
ln_RS_dv <- log(dolly$dv_Smolts/dolly$dv_Adults)
ln_RS_ct <- log(cutty$ct_Smolts/cutty$ct_Adults)
dat <- rbind(ln_RS_sh,ln_RS_dv,ln_RS_ct)
a <- apply(dat,1,mean,na.rm=TRUE)
inits.list <- list(x0=array(matrix(rep(0,m), nrow=m),c(m,m,Nspecies)),
                   A=matrix(a,Nspecies,1))
keoghDLMspecies <- MARSS(dat, model=model.listDLMspecies,control=list(maxit=2000,conv.test.slope.tol=0.1),inits=inits.list)


keoghDLMspecies$model
keoghDLMspecies$states
mean(colSums(keoghDLMall$states))
MARSSparamCIs(keoghDLMall)
keoghDLMall$coef

coef(keoghDLMspecies)

keoghAllfit <- augment(keoghDLMspecies, interval="confidence")
keoghAllfit$Year <- keoghAllfit$t + 1975
keoghAllfit$Species <- keoghAllfit$.rownames
p <- ggplot(data = keoghAllfit) + 
  geom_line(aes(Year, .fitted)) +
  geom_ribbon(aes(x=Year, ymin=.conf.low, ymax=.conf.up), linetype=2, alpha=0.5)
p <- p + geom_point(data=keoghAllfit, mapping = aes(x=Year, y=y))
p + xlab("") + ylab("ln_RS") + facet_wrap(~Species)

layout(1)
plot(sh_Smolts~sh_Adults,data=steelhead,pch=21,bg=ifelse(steelhead$Year>=1991,"orange","dodgerblue"))
curve(exp(keoghDLM$coef["D.(Y1,alpha)"]+keoghDLM$coef["D.(Y1,precip)"]*0)*x*exp(colSums(keoghDLM$states)[1]*x),add=TRUE,lwd=2,col="dodgerblue")
curve(exp(keoghDLM$coef["D.(Y1,alpha)"]+keoghDLM$coef["D.(Y1,precip)"]*0)*x*exp(colSums(keoghDLM$states)[Nyears]*x),add=TRUE,lwd=2,col="orange",lty=2)


ylabs <- rep(c(expression(alpha[t]),expression(beta[t]~"adults")),Nspecies)#, expression("pdo 2"~beta[t]),expression("pdo 3"~beta[t]))

titles <- c("Steelhead","Dolly Varden","Cutthroat")

colr <- rep(c("darkgreen","dodgerblue","orange"),Nspecies)

layout(matrix(1:((m+1)*Nspecies),nrow=(m+1),ncol=Nspecies))
par(mar=c(5,5,2,1))
state <- 1
for(j in 1:Nspecies)
  {
  for(i in 1:m) 
  {
    mn <- keoghDLMspecies$states[state,]
    se <- keoghDLMspecies$states.se[state,]
    plot(steelhead$Year,mn,xlab="",ylab=ylabs[state],bty="n",xaxt="n",type="n",ylim=c(min(mn-2*se),max(mn+2*se)),cex.lab=1.5)
    if(i==1){
      title(titles[j])
    }
    
    lines(steelhead$Year, rep(0,Nyears), lty="dashed")
    lines(steelhead$Year, mn, col=colr[i], lwd=3)
    lines(steelhead$Year, mn+2*se, col=colr[i])
    lines(steelhead$Year, mn-2*se, col=colr[i])
    state <- state+1
  }
  mn <- pmax(0,-log(keoghDLMspecies$states[j*m-1,])/keoghDLMspecies$states[j*m,])
  plot(steelhead$Year,mn,xlab="",ylab=expression(~K[t]),bty="n",xaxt="n",type="n",cex.lab=2)
  lines(steelhead$Year, rep(0,Nyears), lty="dashed")
  lines(steelhead$Year, mn, col=colr[3], lwd=3)
  axis(1,at=seq(min(steelhead$Year),max(steelhead$Year),5),cex=2)
  mtext("Brood year", 1, line=3,cex=1)
  abline(v=1991)
}
