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
c3 <- model.matrix(~-1+mean_temp_egg+total_rain_egg,data=enviro_prod)

x3 <- model.matrix(~-1+Stock,data=sh_annual)

dat <- list("N"=nrow(x3),
            "K"=ncol(x3),
            "J"=ncol(c3),
            "c"=c3,
            "x"=x3,
            "y"=log(sh_annual$Recruits/sh_annual$Stock),
            "family"=1,
            "type"=0)

fit <- stan(file="Stan code/dlm_slope.stan",data=dat, iter=5000,chains=4,cores=4,control=list("adapt_delta"=0.8))

predictions <- extract(fit)$pred
plot(colMeans(predictions),ylim=range(predictions))
points(dat$y,pch=21,bg="orange")

bs <- extract(fit)$beta
plot(colMeans(bs[,,1]))
plot(colMeans(bs[,,2]))

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
