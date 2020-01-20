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
options(mc.cores = parallel::detectCores(logical=FALSE))
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7')

keogh <- readRDS("Keogh_newJuv_enviro.rds")
run_time <- readRDS("Data/steelhead_run.rds")
sh_annual <- readRDS("Data/steelhead_run_annual.rds")
sh_annual$time <- 1:length(sh_annual$year)
sh_annual <- subset(sh_annual,year<=2015 & year>=1976)
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
c1 <- model.matrix(~seals+npgo,data=enviro)

XXvars <- c("total_rain_run","mean_temp_run")
sdSurv_run <- attr(scale(sh_annual[,XXvars],center=TRUE,scale=TRUE),"scaled:scale")
mnSurv_run <- attr(scale(sh_annual[,XXvars],center=TRUE,scale=TRUE),"scaled:center")
enviro_run <- scale(sh_annual[,XXvars],center=TRUE,scale=TRUE)
enviro_run <- data.frame(enviro_run)
colnames(enviro_run) <- XXvars
c2 <- model.matrix(~total_rain_run+mean_temp_run,data=enviro_run)

XXXvars <- c("mean_temp_egg", "total_rain_egg")
sdSurv_prod <- attr(scale(sh_annual[,XXXvars],center=TRUE,scale=TRUE),"scaled:scale")
mnSurv_prod <- attr(scale(sh_annual[,XXXvars],center=TRUE,scale=TRUE),"scaled:center")
enviro_prod <- scale(sh_annual[,XXXvars],center=TRUE,scale=TRUE)
enviro_prod <- data.frame(XXXvars=enviro_prod)
colnames(enviro_prod) <- XXXvars
c3 <- model.matrix(~1,data=enviro_prod)
x3 <- model.matrix(~-1+Stock,data=sh_annual)

dat <- list("N"=nrow(x3),
            "x"=as.numeric(x3),
            "y"=log(sh_annual$Recruits/sh_annual$Stock),
            "family"=1)
plot(y~x,data=dat)
fit <- stan(file="Stan code/dlm_slope.stan",data=dat, iter=10000,chains=4,cores=4,control=list("adapt_delta"=0.95))

summary(fit,pars=c("x0","beta","sigma_obs","sigma_process"),probs=c(0.1,0.9))$summary
predictions <- extract(fit)$pred
plot(colMeans(predictions),ylim=range(predictions),type="l",lwd=2)
points(dat$y,pch=21,bg="orange")

bs <- extract(fit)$beta
ci <- apply(bs,2,quantile,probs=c(0.05,0.95))
plot(colMeans(bs),ylim=range(c(0,ci)),type="l",lwd=2)
polygon(c(1:dat$N,rev(1:dat$N)),c(ci[1,],rev(ci[2,])),col=adjustcolor("grey50",0.5))

ppd <- extract(fit)$y_ppd
ci <- apply(ppd,2,quantile,probs=c(0.05,0.95))
plot(colMeans(ppd),ylim=range(ci),type="l",lwd=2)
polygon(c(1:dat$N,rev(1:dat$N)),c(ci[1,],rev(ci[2,])),col=adjustcolor("grey50",0.5))
points(dat$y,pch=21,bg="orange")

ppd <- extract(fit)$R
ci <- apply(ppd,2,quantile,probs=c(0.05,0.95))
plot(colMeans(ppd),ylim=range(ci),type="l",lwd=2)
polygon(c(1:dat$N,rev(1:dat$N)),c(ci[1,],rev(ci[2,])),col=adjustcolor("grey50",0.5))
points(sh_annual$Recruits,pch=21,bg="orange")

fit2 <- stan(file="Stan code/dlm_int.stan",data=dat, iter=10000,chains=4,cores=4,control=list("adapt_delta"=0.99))

summary(fit2,pars=c("x0","beta","sigma_process"),probs=c(0.1,0.9))$summary
predictions <- extract(fit2)$pred
ci <- apply(predictions,2,quantile,probs=c(0.05,0.95))
plot(colMeans(predictions),ylim=range(predictions),type="l",lwd=2)
polygon(c(1:dat$N,rev(1:dat$N)),c(ci[1,],rev(ci[2,])),col=adjustcolor("grey50",0.5))
points(dat$y,pch=21,bg="orange")

xs <- extract(fit2)$x0
ci <- apply(xs,2,quantile,probs=c(0.05,0.95))
plot(colMeans(xs),ylim=range(c(0,ci)),type="l",lwd=2)
polygon(c(1:dat$N,rev(1:dat$N)),c(ci[1,],rev(ci[2,])),col=adjustcolor("grey50",0.5))

ppd <- extract(fit2)$y_ppd
ci <- apply(ppd,2,quantile,probs=c(0.05,0.95))
plot(colMeans(ppd),ylim=range(ci),type="l",lwd=2)
polygon(c(1:dat$N,rev(1:dat$N)),c(ci[1,],rev(ci[2,])),col=adjustcolor("grey50",0.5))
points(dat$y,pch=21,bg="orange")

ppd <- extract(fit2)$R
ci <- apply(ppd,2,quantile,probs=c(0.05,0.95))
plot(colMeans(ppd),ylim=range(ci),type="l",lwd=2)
polygon(c(1:dat$N,rev(1:dat$N)),c(ci[1,],rev(ci[2,])),col=adjustcolor("grey50",0.5))
points(sh_annual$Recruits,pch=21,bg="orange")

loo_compare(loo(fit),loo(fit2))

post <- extract(fit2)
K <- sapply(1:nrow(post$x0),function(x){-post$x0[x,]/post$beta[x]})
ci <- apply(K,1,quantile,probs=c(0.05,0.95))
plot(rowMeans(K),ylim=range(ci,sh_annual$Recruits),type="l",lwd=2)
polygon(c(1:dat$N,rev(1:dat$N)),c(ci[1,],rev(ci[2,])),col=adjustcolor("grey50",0.5))
points(sh_annual$Recruits,pch=21,bg="red")

plot(sh_annual$logit_surv,exp(colMeans(extract(fit2)$x0)))

alphas <- exp(colMeans(extract(fit2)$x0))
beta <- mean(extract(fit2)$beta)

plot(Recruits~Stock,data=sh_annual,pch=21,bg=ifelse(sh_annual$year>=1990,"orange","dodgerblue"),xlim=range(c(0,rowMeans(K))),ylim=range(c(0,extract(fit2)$R)))
curve(mean(alphas[1:10])*x*exp(beta*x),add=TRUE,lwd=2,col="dodgerblue")
curve(mean(alphas[11:20])*x*exp(beta*x),add=TRUE,lwd=2,col="grey50")
curve(mean(alphas[21:30])*x*exp(beta*x),add=TRUE,lwd=2,col="orange")
curve(mean(alphas[31:40])*x*exp(beta*x),add=TRUE,lwd=2,col="orange")
curve(mean(alphas[regime])*x*exp(beta*x),add=TRUE,lwd=2,col="orange")

regime <- which.max(pareto_k_values(loo(fit2)))
plot(log(Recruits/Stock)~Stock,data=sh_annual,pch=21,bg=0,col=0,xlim=range(c(0,rowMeans(K))))
for(i in 1:dat$N)
{
  points(dat$x[i],dat$y[i],pch=21,bg="orange")
  curve(log(mean(alphas[i]))+beta*x,add=TRUE)
  Sys.sleep(1)
}
# all models
dat <- list("N"=nrow(x1),
            "K"=ncol(x1),
            "J"=ncol(x2),
            "M"=ncol(x3),
            "c1"=c1,
            "c2"=c2,
            "c3"=c3,
            "x1"=x1,
            "x2"=x2,
            "x3"=x3,
            "y1"=sh_annual$logit_surv,
            "y2"=sh_annual$run,
            "y3"=log(sh_annual$Recruits/sh_annual$Stock))
