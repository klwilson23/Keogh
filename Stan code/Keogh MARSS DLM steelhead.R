source("some functions.R")
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(MARSS)
library(broom)
library(wesanderson)
#library(devtools)
#devtools::install_github("nwfsc-timeseries/atsar")
library(atsar)

keogh <- readRDS("Keogh_newJuv_enviro.rds")
keogh_long <- subset(keogh,Year<=2015 & Year>=1976)
keogh_long <- subset(keogh_long,Species!="Chum")

adults <- subset(keogh_long,select = c(Year,Species,Stock))
adults <- reshape(adults,direction = "wide",idvar="Year",timevar="Species")
recruits <- subset(keogh_long,select = c(Year,Species,Recruits))
recruits <- reshape(recruits,direction = "wide",idvar="Year",timevar="Species")

marSurv <- subset(keogh_long,select = c(Year,Species,Stock,juvCohort))
marSurv$Surv <- pmax(1e-3,pmin(1-1e-3,marSurv$Stock/marSurv$juvCohort),na.rm=TRUE)
marSurv$logitSurv <- log(marSurv$Surv/(1-marSurv$Surv))
marSurv <- subset(marSurv,select = c(Year,Species,logitSurv))
ocean_survival <- reshape(marSurv,direction = "wide",idvar="Year",timevar="Species")

summary(lm(ocean_survival$logitSurv.Steelhead~keogh_long$seals[keogh_long$Species=="Steelhead"]+keogh_long$oceanSalmon[keogh_long$Species=="Steelhead"]))

juv_enviro <- subset(keogh_long,select = c(Year,Species,sumTemp,sumRain,winTemp,winRain,freshCoho))
fresh_enviro <- reshape(juv_enviro,direction = "wide",idvar="Year",timevar="Species")
#fresh_enviro <- fresh_enviro[,-match(c("freshSteel.Pink","freshDolly.Pink","freshCutt.Pink","freshPink.Pink","freshCoho.Pink"),colnames(fresh_enviro))]
adult_enviro <- subset(keogh_long,select = c(Year,Species,seals,npgo,oceanSalmon))
adult_enviro$oceanSalmon <- residuals(lm(oceanSalmon~seals:Species,data=adult_enviro))

ocean_enviro <- reshape(adult_enviro,direction = "wide",idvar="Year",timevar="Species")

freshEnviroNew <- fresh_enviro
sdCovarsFresh <- attr(scale(fresh_enviro[,-1],center=TRUE,scale=TRUE),"scaled:scale")
mnCovarsFresh <- attr(scale(fresh_enviro[,-1],center=TRUE,scale=TRUE),"scaled:center")
freshCovarScale <- scale(fresh_enviro[,-1],center=TRUE,scale=TRUE)
#covarScale[is.na(covarScale)] <- 0
oceanEnviroNew <- ocean_enviro
sdCovarsOcean <- attr(scale(ocean_enviro[,-1],center=TRUE,scale=TRUE),"scaled:scale")
mnCovarsOcean <- attr(scale(ocean_enviro[,-1],center=TRUE,scale=TRUE),"scaled:center")
oceanCovarScale <- scale(ocean_enviro[,-1],center=TRUE,scale=TRUE)

# all environmental predictors
all_enviro <- subset(keogh_long,select = c(Year,Species,sumTemp,winRain,freshCoho,seals,npgo,oceanSalmon))
all_enviro$seals <- residuals(lm(seals~oceanSalmon:Species,data=all_enviro))
#all_enviro$seals <- log(all_enviro$seals)
all_covars <- reshape(all_enviro,direction = "wide",idvar="Year",timevar="Species")
allEnviroNew <- all_covars
sdCovarsAll <- attr(scale(all_covars[,-1],center=TRUE,scale=TRUE),"scaled:scale")
mnCovarsAll <- attr(scale(all_covars[,-1],center=TRUE,scale=TRUE),"scaled:center")
allCovarScale <- scale(all_covars[,-1],center=TRUE,scale=TRUE)

# multiple species: dolly varden, cutthroat trout, pink salmon, coho salmon
# including process and observation error
# precipitation covariates only affect observation model
# time-varying beta & alpha
# run a DLM on stock-recruitment for steelhead only

# including process and observation error
# adult covariates only affect observation model
# time-varying affect on adults
C <- c <- A <- U <- "zero"
Nyears <- nrow(adults)
Nspecies <- 1
m <- 1
Z <- array(NA,c(1,m,Nyears),dimnames=list("Species"="steelhead","beta"=1:m,"year"=adults$Year))
Z[1,1,] <- adults$Stock.Steelhead # time-varying beta
D <- "unconstrained"
d <- rbind(rep(1,Nyears))
row.names(d) <- c("alpha")
Q <- matrix(list(0),m,m)         ## 2x2; all 0 for now
diag(Q) <- c("q.adults") ## 2x2; diag = (q1,q2)
B <- diag(m)
R <- "unconstrained"
x0 <- "zero"
initx <- 0
ln_RS <- log(recruits$Recruits.Steelhead/adults$Stock.Steelhead)
model.listDLM <- list(B=diag(Nspecies),U=U,Q="unconstrained",Z=Z,A=A,R=R,d=d,D=D,x0=x0,tinitx=initx)
keoghDLM <- MARSS(ln_RS, model=model.listDLM,control=list(maxit=2000,conv.test.slope.tol=0.1))

adults$ln_RS <- log(recruits$Recruits.Steelhead/adults$Stock.Steelhead)
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
p <- p + geom_point(data=adults, mapping = aes(x=Year, y=ln_RS))
p + xlab("") + ylab("ln_RS")

layout(1)
plot(adults$Stock.Steelhead,recruits$Recruits.Steelhead,pch=21,bg=ifelse(adults$Year>=1991,"orange","dodgerblue"))
curve(exp(keoghDLM$coef["D.D"])*x*exp(colSums(keoghDLM$states)[1]*x),add=TRUE,lwd=2,col="dodgerblue")
curve(exp(keoghDLM$coef["D.D"])*x*exp(colSums(keoghDLM$states)[Nyears]*x),add=TRUE,lwd=2,col="orange",lty=2)


ylabs <- c(expression("adults"~beta[t]))#, expression("pdo 2"~beta[t]),expression("pdo 3"~beta[t]))
colr <- c("darkgreen","dodgerblue","orange")
par(mfrow=c(m,1), mar=c(4,4,0.1,0), oma=c(0,0,2,0.5))
for(i in 1:m) {
  mn <- keoghDLM$states[i,]
  se <- keoghDLM$states.se[i,]
  plot(adults$Year,mn,xlab="",ylab=ylabs[i],bty="n",xaxt="n",type="n",
       ylim=c(min(mn-2*se),max(mn+2*se)))
  lines(adults$Year, rep(0,Nyears), lty="dashed")
  lines(adults$Year, mn, col=colr[i], lwd=3)
  lines(adults$Year, mn+2*se, col=colr[i])
  lines(adults$Year, mn-2*se, col=colr[i])
}
axis(1,at=seq(min(adults$Year),max(adults$Year),5))
mtext("Brood year", 1, line=3)
abline(v=1991)
