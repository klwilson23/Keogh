source("some functions.R")
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(MARSS)
library(broom)
keogh <- readRDS("Keogh_SR_enviro_long.rds")
keogh_long <- subset(keogh,Year<=2015 & Year>=1976)
keogh_long <- subset(keogh_long,Species!="Chum")

keogh <- subset(keogh_long,select = c(Year,Species,Stock,Recruits,juvCohort))
keogh <- reshape(keogh,direction = "wide",idvar="Year",timevar="Species")

environment <- subset(keogh_long,select = c(Year,Species,sumTemp,sumRain,winTemp,winRain,freshCoho,freshSteel,freshCutt,freshDolly,freshPink,seals,npgo,mei,oceanSalmon))
enviro <- reshape(environment,direction = "wide",idvar="Year",timevar="Species")

sdCovars <- attr(scale(enviro[,-1],center=TRUE,scale=TRUE),"scaled:scale")
covarScale <- scale(enviro[,-1],center=TRUE,scale=TRUE)
#covarScale[is.na(covarScale)] <- 0

# get estimates of missing data from DLM analysis
Nyears <- length(enviro$Year)
years <- enviro$Year
covars <- t(as.matrix(covarScale))
colnames(covars) <- enviro$Year
ns <- nrow(covars)
B <- "zero"
Q <- "unconstrained"
R <- diag(0.01, ns)
U <- "zero"
A <- "unequal"
x0 <- "zero"
m <- apply(covars, 1, mean, na.rm = TRUE)
mod.list.corr = list(B = B, Q = Q, R = R, U = U, x0 = x0, A = A, tinitx = 0)
fit.corr <- MARSS(covars, model = mod.list.corr, control = list(maxit = 5000), inits = list(A = matrix(m, ns, 1)))

d <- augment(fit.corr, interval = "confidence")
d$Year <- d$t + (min(years)-1)
d$covars <- d$.rownames

covarFits <- matrix(sapply(1:nrow(d),function(x){d$.fitted[x]*sdCovars[names(sdCovars)%in%d$covars[x]]}),nrow=ncol(covarScale),ncol=nrow(covarScale),byrow=FALSE)

sh_Adults.hat <- pmax(1e-4,ifelse(is.na(keogh$sh_Adults),adultFits[,1],keogh$sh_Adults))
ct_Adults.hat <- pmax(1e-4,ifelse(is.na(keogh$ct_Adults),adultFits[,3],keogh$ct_Adults))
dv_Adults.hat <- pmax(1e-4,ifelse(is.na(keogh$dv_Adults),adultFits[,2],keogh$dv_Adults))
co_Adults.hat <- pmax(1e-4,ifelse(is.na(keogh$co_Adults),adultFits[,5],keogh$co_Adults))
pk_Adults.hat <- pmax(1e-4,ifelse(is.na(keogh$pk_Adults),adultFits[,4],keogh$pk_Adults))


# get estimates of adult densities from DFA analysis:
# run a DFA on the standardized adult data
adultDat <- keogh[,c("Stock.Steelhead","Stock.Dolly Varden","Stock.Cutthroat","Stock.Pink","Stock.Coho")]
sdScale <- attr(scale(adultDat,center=FALSE,scale=TRUE),"scaled:scale")
adultDat <- scale(adultDat,center=FALSE,scale=TRUE)

adultDat <- t(as.matrix(adultDat))
colnames(adultDat) <- keogh$Year
ns <- nrow(adultDat)
B <- "diagonal and equal"
Q <- "unconstrained"
R <- diag(0.01,ns)
U <- "zero"
A <- "unequal"
x0 <- "zero"
mod.list = list(B=B, Q=Q, R=R, U=U, x0=x0, A=A,tinitx = 1)

mod.list.noTime = list(B = matrix(1), U = "zero", Q = "zero", 
                       Z = matrix(1,ns,1), A = "zero", R = "diagonal and equal", x0 = matrix("mu"))

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
mod.list.dfa = list(B = B, Z = Z, Q = Q, R = R, U = U, A = A, x0 = x0)

m <- apply(adultDat, 1, mean, na.rm=TRUE)
fit <- MARSS(adultDat, model=mod.list, control=list(minit=200,maxit=50000+200), inits=list(A=matrix(m,ns,1)))

fit <- MARSS(adultDat, model=mod.list.noTime, control=list(minit=200,maxit=50000+200))
summary(fit)
print(fit)
plot(fit)
d <- augment(fit, interval = "confidence")
d$Year <- d$t + 1975
d$Species <- d$.rownames

adultFits <- matrix(sapply(1:nrow(d),function(x){d$.fitted[x]*sdScale[names(sdScale)%in%d$Species[x]]}),nrow=ncol(adultDat),ncol=nrow(adultDat),byrow=FALSE)

sh_Adults.hat <- pmax(1e-4,ifelse(is.na(keogh$sh_Adults),adultFits[,1],keogh$sh_Adults))
ct_Adults.hat <- pmax(1e-4,ifelse(is.na(keogh$ct_Adults),adultFits[,3],keogh$ct_Adults))
dv_Adults.hat <- pmax(1e-4,ifelse(is.na(keogh$dv_Adults),adultFits[,2],keogh$dv_Adults))
co_Adults.hat <- pmax(1e-4,ifelse(is.na(keogh$co_Adults),adultFits[,5],keogh$co_Adults))
pk_Adults.hat <- pmax(1e-4,ifelse(is.na(keogh$pk_Adults),adultFits[,4],keogh$pk_Adults))

# multiple species: dolly varden, cutthroat trout, pink salmon, coho salmon
# including process and observation error
# precipitation covariates only affect observation model
# time-varying beta & alpha
# run a DLM on stock-recruitment for steelhead only

C <- c <- A <- U <- "zero"
Nyears <- nrow(steelhead)
Nspecies <- 5
m <- 2 # number of time-varying parameters
Z <- array(0,c(Nspecies,m*Nspecies,Nyears))
Z[1,1,] <- rep(1,Nyears) # time-varying alpha
Z[1,2,] <- steelhead$sh_Adults # time-varying beta
Z[2,3,] <- rep(1,Nyears) # time-varying alpha
Z[2,4,] <- dv_Adults.hat # time-varying beta
Z[3,5,] <- rep(1,Nyears) # time-varying alpha
Z[3,6,] <- ct_Adults.hat # time-varying beta
Z[4,7,] <- rep(1,Nyears) # time-varying alpha
Z[4,8,] <- ct_Adults.hat # time-varying beta
Z[5,9,] <- rep(1,Nyears) # time-varying alpha
Z[5,10,] <- ct_Adults.hat # time-varying beta
#Z[1,3,] <- covarScale[,"precip_1"] # time-varying precipitation

D <- "unconstrained"
d <- rbind(covarScale[,"precip_1"])
row.names(d) <- c("precip")

C <- "unconstrained" # NPGO, seals, pacific salmon, and mei affects observation (time-varying component)
c <- rbind(covarScale[,"npgo"],covarScale[,"seal_density"],covarScale[,"total"],covarScale[,"mei"])
row.names(c) <- c("npgo","seal_density","total","mei")

Q <- matrix(list(0),m*Nspecies,m*Nspecies)        ## 2x2; all 0 for now
#diag(Q) <- rep(c("q.alpha","q.beta"),Nspecies)
diag(Q) <- c("q.sh_alpha","q.sh_adults","q.dv_alpha","q.dv_adults","q.ct_alpha","q.ct_adults") ## 2x2; diag = (q1,q2)
Q[1,3] <- "q.corr_alpha1"
Q[3,1] <- "q.corr_alpha1"
Q[3,5] <- "q.corr_alpha2"
Q[5,3] <- "q.corr_alpha2"
Q[1,5] <- "q.corr_alpha3"
Q[5,1] <- "q.corr_alpha3"
Q[1,7] <- "q.corr_alpha4"
Q[7,1] <- "q.corr_alpha4"
Q[1,9] <- "q.corr_alpha5"
Q[9,1] <- "q.corr_alpha5"

Q[2,4] <- "q.corr_beta1"
Q[4,2] <- "q.corr_beta1"
Q[4,6] <- "q.corr_beta2"
Q[6,4] <- "q.corr_beta2"
Q[2,6] <- "q.corr_beta3"
Q[6,2] <- "q.corr_beta3"
Q[2,8] <- "q.corr_beta4"
Q[8,2] <- "q.corr_beta4"
Q[2,10] <- "q.corr_beta5"
Q[10,2] <- "q.corr_beta5"

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