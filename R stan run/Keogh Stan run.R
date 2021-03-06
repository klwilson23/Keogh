source("some functions.R")
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(MARSS)
library(broom)
library(wesanderson)
library(rstan)
library(loo)
library(rethinking)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

keogh <- readRDS("Keogh_newJuv_enviro.rds")
run_time <- readRDS("Data/steelhead_run.rds")
sh_annual <- readRDS("Data/steelhead_run_annual.rds")
keogh_long <- subset(keogh,Year<=2015 & Year>=1976)
keogh_long <- subset(keogh_long,Species!="Chum")
keogh_long$Species <- factor(keogh_long$Species,levels=unique(keogh_long$Species))
keogh_long$marSurv <- keogh_long$Stock/keogh_long$juvCohort
keogh_long$logitSurv <- log(keogh_long$marSurv/(1-keogh_long$marSurv))
keogh_long$prod <- log(keogh_long$Recruits/keogh_long$Stock)

plot(sh_annual$year,sh_annual$sex,type="l")

boxplot(run~ifelse(sex>=0.6,"F",ifelse(sex<=0.4,"M","E")),data=sh_annual)

X <- model.matrix(~-1+Species+Stock:Species+sumRain:Species+sumTemp:Species+winRain:Species+winTemp:Species,data=keogh_long)

trends <- model.matrix(~-1+logitSurv:Species+seals:Species,data=keogh_long)

Xvars <- c("seals","npgo")
sdSurv_sh <- attr(scale(sh_annual[,Xvars],center=TRUE,scale=TRUE),"scaled:center")
mnSurv_sh <- attr(scale(sh_annual[,Xvars],center=TRUE,scale=TRUE),"scaled:center")
enviro <- scale(sh_annual[,Xvars],center=TRUE,scale=TRUE)
enviro <- data.frame(enviro)
sh_trends <- model.matrix(~-1+seals+npgo,data=enviro)

XXvars <- c("total_rain_run","mean_temp_run")
sdSurv_run <- attr(scale(sh_annual[,XXvars],center=TRUE,scale=TRUE),"scaled:center")
mnSurv_run <- attr(scale(sh_annual[,XXvars],center=TRUE,scale=TRUE),"scaled:center")
enviro_run <- scale(sh_annual[,XXvars],center=TRUE,scale=TRUE)
enviro_run <- data.frame(enviro_run)
run_trends <- model.matrix(~-1+total_rain_run+mean_temp_run,data=enviro_run)

dat <- list("N"=nrow(sh_trends),
            "K"=ncol(sh_trends),
            "X"=sh_trends,
            "lSurv"=sh_annual$logit_surv,
            "J"=ncol(run_trends),
            "XX"=run_trends,
            "run_time"=sh_annual$run)
trackPars <- c("beta_surv","beta_run","bSurv","pS0","run0","obs_sigma_surv","obs_sigma_run","pro_sigma_surv","pro_sigma_run","pro_devS","pro_devR","surv_new","run_new","mnSurv","mnRun")
fit <- stan(file = "Stan code/Keogh Surv DLM.stan", pars=trackPars,data=dat, iter=5000,chains=4,cores=4,control=list("adapt_delta"=0.9))

summary(fit, pars=trackPars[!grepl("new",trackPars)],probs=c(0.025,0.975))$summary
mypost <- as.data.frame(fit)
surv_ppd <- mypost[,grep("surv_new",colnames(mypost))]
mn_ppd <- apply(surv_ppd,2,mean)
ci_ppd <- apply(surv_ppd,2,HPDI,prob=0.89)
plot(1/(1+exp(-dat$lSurv)),mn_ppd,pch=21,bg="grey50",ylim=range(ci_ppd),xlim=range(ci_ppd), main = "Survival",xlab="Observed survival",ylab="Posterior predictive")
segments(x0=1/(1+exp(-dat$lSurv)),y0=ci_ppd[1,],y1=ci_ppd[2,],lwd=1)
abline(b=1,a=0,lwd=2,lty=1,col="red")


run_ppd <- mypost[,grep("run_new",colnames(mypost))]
mn_ppd <- apply(run_ppd,2,mean)
ci_ppd <- apply(run_ppd,2,HPDI,prob=0.89)
plot(dat$run_time,mn_ppd,pch=21,bg="grey50",ylim=range(ci_ppd),xlim=range(ci_ppd), main = "Survival",xlab="Observed survival",ylab="Posterior predictive")
segments(x0=dat$run_time,y0=ci_ppd[1,],y1=ci_ppd[2,],lwd=1)
abline(b=1,a=0,lwd=2,lty=1,col="red")

logit_surv <- colMeans(mypost[,grep("mnSurv",colnames(mypost))])
mn_surv <- exp(logit_surv)/(1+exp(logit_surv))

plot(colMeans(run_ppd),type="l")
plot(colMeans(surv_ppd),type="l")
#lines(mn_surv,col="red")

plot(dat$X[,1],mn_surv,xlab="seal densities",ylab="marine survival (%)")

