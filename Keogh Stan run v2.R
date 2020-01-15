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
#options(mc.cores = parallel::detectCores(logical=FALSE))
#Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

keogh <- readRDS("Keogh_newJuv_enviro.rds")
run_time <- readRDS("Data/steelhead_run.rds")
sh_annual <- readRDS("Data/steelhead_run_annual.rds")
sh_annual$time <- 1:length(sh_annual$year)
keogh_long <- subset(keogh,Year<=2015 & Year>=1976)
keogh_long <- subset(keogh_long,Species!="Chum")
keogh_long$Species <- factor(keogh_long$Species,levels=unique(keogh_long$Species))
keogh_long$marSurv <- keogh_long$Stock/keogh_long$juvCohort
keogh_long$logitSurv <- log(keogh_long$marSurv/(1-keogh_long$marSurv))
keogh_long$prod <- log(keogh_long$Recruits/keogh_long$Stock)

Xvars <- c("seals","npgo")
sdSurv_sh <- attr(scale(sh_annual[,Xvars],center=TRUE,scale=TRUE),"scaled:scale")
mnSurv_sh <- attr(scale(sh_annual[,Xvars],center=TRUE,scale=TRUE),"scaled:center")
enviro <- scale(sh_annual[,Xvars],center=TRUE,scale=TRUE)
enviro <- data.frame(Xvars=enviro)
colnames(enviro) <- Xvars
sh_trends <- model.matrix(~seals+npgo,data=enviro)

XXvars <- c("total_rain_run","mean_temp_run")
sdSurv_run <- attr(scale(sh_annual[,XXvars],center=TRUE,scale=TRUE),"scaled:scale")
mnSurv_run <- attr(scale(sh_annual[,XXvars],center=TRUE,scale=TRUE),"scaled:center")
enviro_run <- scale(sh_annual[,XXvars],center=TRUE,scale=TRUE)
enviro_run <- data.frame(enviro_run)
colnames(enviro_run) <- XXvars
run_trends <- model.matrix(~total_rain_run+mean_temp_run,data=enviro_run)

dat <- list("N"=nrow(sh_trends),
            "K"=ncol(sh_trends),
            "X"=sh_trends,
            "lSurv"=sh_annual$logit_surv,
            "J"=ncol(run_trends),
            "XX"=run_trends,
            "run_time"=sh_annual$run)

fit <- stan(file = "Stan code/Keogh Surv.stan", data=dat, iter=5000,chains=4,cores=4,control=list("adapt_delta"=0.8))
summary(fit, pars=c("beta_surv","beta_run","sigma_surv","sigma_run","phi_surv","phi_run"),probs=c(0.025,0.975))$summary

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

datNew <- list("N"=nrow(sh_trends),
               "K"=ncol(sh_trends),
               "X"=sh_trends,
               "lSurv"=sh_annual$logit_surv,
               "J"=ncol(run_trends),
               "XX"=run_trends,
               "run_time"=sh_annual$run)
trackPars <- c("beta_surv","beta_run","bSurv","pS0","run0","obs_sigma_surv","obs_sigma_run","pro_sigma_surv","pro_sigma_run","pro_devS","pro_devR","surv_new","run_new","mnSurv","mnRun","log_lik1","log_lik2")
fit2 <- stan(file = "Stan code/Keogh Surv DLM.stan", data=datNew,pars=trackPars, iter=5000,chains=4,cores=1,control=list("adapt_delta"=0.9))

summary(fit2, pars=c("beta_surv","beta_run","bSurv","pS0","run0","obs_sigma_surv","obs_sigma_run","pro_sigma_surv","pro_sigma_run","pro_devS","pro_devR"),probs=c(0.1,0.9))$summary

mypost <- as.data.frame(fit)
surv_ppd <- extract(fit)$surv_new
mn_ppd <- colMeans(surv_ppd)
ci_ppd <- apply(surv_ppd,2,HPDI,prob=0.89)
plot(1/(1+exp(-dat$lSurv)),mn_ppd,pch=21,bg="grey50",ylim=range(ci_ppd),xlim=range(1/(1+exp(-dat$lSurv))), main = "Survival",xlab="Observed survival",ylab="Posterior predictive")
segments(x0=1/(1+exp(-dat$lSurv)),y0=ci_ppd[1,],y1=ci_ppd[2,],lwd=1)
abline(b=1,a=0,lwd=2,lty=1,col="red")

run_ppd <- extract(fit)$run_new
mn_ppd <- colMeans(run_ppd)
ci_ppd <- apply(run_ppd,2,HPDI,prob=0.89)
plot(dat$run_time,mn_ppd,pch=21,bg="grey50",ylim=range(ci_ppd),xlim=range(dat$run_time), main = "Run time",xlab="Observed run time",ylab="Posterior predictive")
segments(x0=dat$run_time,y0=ci_ppd[1,],y1=ci_ppd[2,],lwd=1)
abline(b=1,a=0,lwd=2,lty=1,col="red")

run_ci <- apply(extract(fit)$run_new,2,HPDI,prob=0.95)
plot(colMeans(run_ppd),ylim=range(run_ci),lwd=2,type="l")
polygon(x=c(1:nrow(dat$X),rev(1:nrow(dat$X))),y=c(run_ci[1,],rev(run_ci[2,])),col=adjustcolor("grey",0.5))
lines(colMeans(run_ppd),lwd=2)
points(dat$run_time,pch=21,bg="red")

surv_ci <- apply(extract(fit)$surv_new,2,HPDI,prob=0.95)
plot(colMeans(extract(fit)$surv_new),ylim=range(surv_ci),lwd=2,type="l")
polygon(x=c(1:nrow(dat$X),rev(1:nrow(dat$X))),y=c(surv_ci[1,],rev(surv_ci[2,])),col=adjustcolor("grey",0.5))
lines(colMeans(extract(fit)$surv_new),lwd=2)
points(sh_annual$Stock/sh_annual$juvCohort,pch=21,bg="red")

# seals:
surv_ppd <- extract(fit)$surv_new[,order(sh_annual$seals)]
surv_ci <- apply(surv_ppd,2,HPDI,prob=0.95)
plot(sort(sh_annual$seals),colMeans(surv_ppd),ylim=range(surv_ci),lwd=2,type="l")
polygon(x=c(sort(sh_annual$seals),rev(sort(sh_annual$seals))),y=c(surv_ci[1,],rev(surv_ci[2,])),col=adjustcolor("grey",0.5))
lines(sort(sh_annual$seals),colMeans(surv_ppd),lwd=2)

plot(sh_annual$seals,colMeans(extract(fit)$surv_new))
plot(sh_annual$total_rain_run,colMeans(extract(fit)$run_new))
loo_1 <- loo(extract_log_lik(fit,parameter_name = c("log_lik1","log_lik2")),cores=1)
loo_2 <- loo(extract_log_lik(fit2,parameter_name = c("log_lik1","log_lik2")),cores=1)
loo_compare(loo_1,loo_2)
