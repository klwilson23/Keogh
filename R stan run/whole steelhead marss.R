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
sh_annual$logit_surv[1] <- NA

sh_annual$log_adults <- log(sh_annual$Stock)
Xvars <- c("seals","npgo","oceanSalmon")
sdSurv_sh <- attr(scale(sh_annual[,Xvars],center=TRUE,scale=TRUE),"scaled:scale")
mnSurv_sh <- attr(scale(sh_annual[,Xvars],center=TRUE,scale=TRUE),"scaled:center")
enviro <- scale(sh_annual[,Xvars],center=TRUE,scale=TRUE)
enviro <- data.frame(Xvars=enviro)
colnames(enviro) <- Xvars
x1 <- model.matrix(~-1+seals+npgo+oceanSalmon,data=enviro)

XXvars <- c("total_rain_run","mean_temp_run")
sdSurv_run <- attr(scale(sh_annual[,XXvars],center=TRUE,scale=TRUE),"scaled:scale")
mnSurv_run <- attr(scale(sh_annual[,XXvars],center=TRUE,scale=TRUE),"scaled:center")
enviro_run <- scale(sh_annual[,XXvars],center=TRUE,scale=TRUE)
enviro_run <- data.frame(enviro_run)
colnames(enviro_run) <- XXvars
x2 <- model.matrix(~-1+total_rain_run+mean_temp_run,data=enviro_run)

XXXvars <- c("mean_temp_egg", "total_rain_egg")
sdSurv_prod <- attr(scale(sh_annual[,XXXvars],center=TRUE,scale=TRUE),"scaled:scale")
mnSurv_prod <- attr(scale(sh_annual[,XXXvars],center=TRUE,scale=TRUE),"scaled:center")
enviro_prod <- scale(sh_annual[,XXXvars],center=TRUE,scale=TRUE)
enviro_prod <- data.frame(XXXvars=enviro_prod)
colnames(enviro_prod) <- XXXvars
xx3 <- model.matrix(~-1+mean_temp_egg+total_rain_egg,data=enviro_prod)
x3 <- model.matrix(~-1+Stock,data=sh_annual)

# all models
dat <- list("N"=nrow(x1),
            "N_obs"=sum(!is.na(sh_annual$logit_surv)),
            "K"=ncol(x1),
            "J"=ncol(x2),
            "x1"=x1,
            "x2"=x2,
            "x3"=as.numeric(x3),
            "y1_obs"=sh_annual$logit_surv[!is.na(sh_annual$logit_surv)],
            "y2"=sh_annual$run,
            "y3"=log(sh_annual$Recruits/sh_annual$Stock),
            "M"=ncol(xx3),
            "xx3"=xx3,
            "init_s0"=mean(sh_annual$logit_surv[1:10],na.rm=TRUE))

fit <- stan(file="Stan code/steelhead dlm lnorm.stan",data=dat, iter=2000,chains=4,cores=4,control=list("adapt_delta"=0.99,"max_treedepth"=15))

saveRDS(fit,file="~/Google Drive/SFU postdoc/Keogh river/Stan fits/keogh steelhead.rds")

summary(fit,pars=c("y1_miss","beta_surv","beta_adults","beta_run","beta_rec","beta_rec_cov"),probs=c(0.025,0.975))$summary

summary(fit,pars=c("sigma_surv_pro","sigma_surv_obs","sigma_adult_obs","sigma_adult_pro","sigma_run_obs","sigma_run_pro","sigma_rec_process","sigma_rec_obs"),probs=c(0.025,0.975))$summary

summary(fit,pars=c("s0","a0","r0","x0"),probs=c(0.025,0.975))$summary

#pairs(fit,pars=c("sigma_surv_pro","sigma_surv_obs","sigma_run_pro","sigma_run_obs","sigma_rec_process","sigma_rec_obs"))
# survival
jpeg(paste("Figures/Steelhead cycle marss ",Sys.Date(),".jpeg",sep=""),res=800,height=7,width=8,units="in")
regime <- which(sh_annual==1991)
matLayout <- matrix(0,nrow=15,ncol=15)
matLayout[1:3,1:11] <- 1
matLayout[1:3,12:15] <- 2

matLayout[4:6,1:11] <- 3
matLayout[4:6,12:15] <- 4

matLayout[7:9,1:11] <- 5
matLayout[7:9,12:15] <- 6

matLayout[10:12,1:11] <- 7
matLayout[10:12,12:15] <- 8

matLayout[13:15,1:11] <- 9
matLayout[13:15,12:15] <- 10

ptSize <- 0.7

layout(matLayout)
par(mar=c(3.5,4,0.5,0.5))
ppd <- extract(fit)$y1_ppd
ci <- apply(ppd,2,quantile,probs=c(0.05,0.95))
plot(colMeans(ppd),ylim=range(ci,quantile(1/(1+exp(-extract(fit)$y1_miss)),probs=c(0.05,0.95))),type="l",lwd=2,xlab="",ylab="Marine survival (yr-1)",xaxt="n")
axis(1,tick=TRUE,labels=FALSE)
polygon(c(1:dat$N,rev(1:dat$N)),c(ci[1,],rev(ci[2,])),col=adjustcolor("dodgerblue",0.3),border=NA)
polygon(c(regime:40,rev(regime:40)),c(ci[1,regime:40],rev(ci[2,regime:40])),col=adjustcolor("blue",0.3),border=NA)
ppx <- extract(fit)$y1
ci <- apply((1/(1+exp(-ppx))),2,quantile,probs=c(0.05,0.95))
segments(x0=1:ncol(ppx),y0=ci[1,],y1=ci[2,],lwd=2,col="darkorange")
points(colMeans(1/(1+exp(-ppx))),pch=21,bg="darkorange")

# some effects on survival
s0 <- extract(fit)$s0
betas <- extract(fit)$beta_surv
seal_seq <- seq(from=-2,to=2,by=0.01)
ppd <- sapply(1:nrow(betas),function(x){betas[x,1] * seal_seq})
ci <- apply(ppd,1,quantile,probs=c(0.05,0.95))
plot(seal_seq,rowMeans(ppd),xlab="",yaxt="n",type="l",ylab="",ylim=range(ci),lwd=2,xaxt="n")
axis(1,line=0)
axis(2,line=0)
mtext("\u0394 survival (logit)",side=2,line=2.5,cex=ptSize)
mtext("Seal densities",side=1,line=2.5,cex=ptSize)
polygon(c(seal_seq,rev(seal_seq)),c(ci[1,],rev(ci[2,])),col=adjustcolor("dodgerblue",0.3),border=NA)

# adult abundance
ppd <- extract(fit)$x3_ppd
ci <- apply(ppd,2,quantile,probs=c(0.05,0.95))
colMed <- apply(ppd,2,quantile,probs=0.5)
plot(colMed,ylim=range(ci),type="l",lwd=2,xlab="",ylab="Adult steelhead",xaxt="n")
axis(1,tick=TRUE,labels=FALSE)
polygon(c(1:dat$N,rev(1:dat$N)),c(ci[1,],rev(ci[2,])),col=adjustcolor("dodgerblue",0.3),border=NA)
polygon(c(regime:40,rev(regime:40)),c(ci[1,regime:40],rev(ci[2,regime:40])),col=adjustcolor("blue",0.3),border=NA)
points(sh_annual$Stock,pch=21,bg="orange")

# some effects on adult abundance
a0 <- extract(fit)$a0
betas <- extract(fit)$beta_adults
surv_seq <- seq(from=-2,to=2,by=0.01)
ppd <- sapply(1:length(betas),function(x){betas[x] * surv_seq})
ci <- apply(ppd,1,quantile,probs=c(0.05,0.95))
plot(surv_seq,rowMeans(ppd),xlab="",yaxt="n",type="l",ylab="",ylim=range(ci),lwd=2,xaxt="n")
axis(1,line=0)
axis(2,line=0)
mtext("\u0394 adults",side=2,line=2.5,cex=ptSize)
mtext("\u0394 survival (logit)",side=1,line=2.5,cex=ptSize)
polygon(c(surv_seq,rev(surv_seq)),c(ci[1,],rev(ci[2,])),col=adjustcolor("dodgerblue",0.3),border=NA)

# run time
ppd <- extract(fit)$y2_ppd
ci <- apply(ppd,2,quantile,probs=c(0.05,0.95))
plot(colMeans(ppd),ylim=range(ci),type="l",lwd=2,xlab="",ylab="Adult run (day of season)",xaxt="n")
axis(1,tick=TRUE,labels=FALSE)
polygon(c(1:dat$N,rev(1:dat$N)),c(ci[1,],rev(ci[2,])),col=adjustcolor("dodgerblue",0.3),border=NA)
polygon(c(regime:40,rev(regime:40)),c(ci[1,regime:40],rev(ci[2,regime:40])),col=adjustcolor("blue",0.3),border=NA)
points(dat$y2,pch=21,bg="orange")

# some effects on run time
r0 <- extract(fit)$r0
betas <- extract(fit)$beta_run
surv_seq <- seq(from=-2,to=2,length.out=25)
ppd <- sapply(1:nrow(betas),function(x){betas[x,2] * surv_seq})
ci <- apply(ppd,1,quantile,probs=c(0.05,0.95))
plot(surv_seq,rowMeans(ppd),xlab="",type="l",ylab="",ylim=range(ci),lwd=2,xaxt="n",yaxt="n")
axis(1,line=0)
axis(2,line=0)
mtext("\u0394 run time",side=2,line=2.5,cex=ptSize)
mtext("ln adult abundance",side=1,line=2.5,cex=ptSize)
polygon(c(surv_seq,rev(surv_seq)),c(ci[1,],rev(ci[2,])),col=adjustcolor("dodgerblue",0.3),border=NA)

# productivity
ppd <- extract(fit)$y3_ppd
ci <- apply(ppd,2,quantile,probs=c(0.05,0.95))
plot(colMeans(ppd),ylim=range(ci),type="l",lwd=2,xlab="",ylab="Smolt productivity (ln R/S)",xaxt="n")
axis(1,labels=FALSE,tick=TRUE)
polygon(c(1:dat$N,rev(1:dat$N)),c(ci[1,],rev(ci[2,])),col=adjustcolor("dodgerblue",0.3),border=NA)
polygon(c(regime:40,rev(regime:40)),c(ci[1,regime:40],rev(ci[2,regime:40])),col=adjustcolor("blue",0.3),border=NA)
points(dat$y3,pch=21,bg="orange")

# some effects on productivity
x0 <- extract(fit)$x0
betas <- extract(fit)$beta_rec_cov
rain_seq <- seq(from=-2,to=2,length.out=25)
ppd <- sapply(1:nrow(betas),function(x){betas[x,3] * rain_seq})
ci <- apply(ppd,1,quantile,probs=c(0.05,0.95))
plot(rain_seq,rowMeans(ppd),xlab="",type="l",ylab="",ylim=range(ci),lwd=2,xaxt="n",yaxt="n")
axis(1,line=0)
axis(2,line=0)
mtext("\u0394 productivity",side=2,line=2.5,cex=ptSize)
mtext("Rainfall during egg incubation",side=1,line=2,cex=ptSize)
polygon(c(rain_seq,rev(rain_seq)),c(ci[1,],rev(ci[2,])),col=adjustcolor("dodgerblue",0.3),border=NA)

# recruitment
ppd <- extract(fit)$R
ci <- apply(ppd,2,quantile,probs=c(0.05,0.95))
plot(colMeans(ppd),ylim=range(ci),type="l",lwd=2,xlab="",ylab="Smolt recruitment",xaxt="n")
axis(1,at=seq(from=0,to=dat$N,length=5),labels=seq(from=min(sh_annual$year)-1,to=max(sh_annual$year),length=5),tick=TRUE)
mtext("Year",side=1,line=2.5,cex=ptSize)
polygon(c(1:dat$N,rev(1:dat$N)),c(ci[1,],rev(ci[2,])),col=adjustcolor("dodgerblue",0.3),border=NA)
polygon(c(regime:40,rev(regime:40)),c(ci[1,regime:40],rev(ci[2,regime:40])),col=adjustcolor("blue",0.3),border=NA)
points(sh_annual$Recruits,pch=21,bg="orange")

# some effects on recruitment
betas <- extract(fit)$beta_rec
alphas <- extract(fit)$x0
adult_seq <- seq(from=0,to=max(dat$x3[1:(regime-1)]),length.out=25)
adult_seq2 <- seq(from=0,to=max(dat$x3[regime:length(dat$x3)]),length.out=25)
ppd <- sapply(1:nrow(betas),function(x){exp(mean(alphas[x,1:(regime-1)]))*adult_seq*exp(betas[x] * adult_seq)})
ci <- apply(ppd,1,quantile,probs=c(0.05,0.95))
ppd2 <- sapply(1:nrow(betas),function(x){exp(mean(alphas[x,regime:length(dat$x3)]))*adult_seq2*exp(betas[x] * adult_seq2)})
ci2 <- apply(ppd2,1,quantile,probs=c(0.05,0.95))
plot(adult_seq,rowMeans(ppd),xlab="",type="l",ylab="",ylim=range(ci,ci2),lwd=2,xaxt="n",yaxt="n")
axis(1,line=0)
axis(2,line=0)
mtext("Smolt recruitment",side=2,line=2.5,cex=ptSize)
mtext("Steelhead adults",side=1,line=2.5,cex=ptSize)
polygon(c(adult_seq,rev(adult_seq)),c(ci[1,],rev(ci[2,])),col=adjustcolor("dodgerblue",0.3),border=NA)
lines(adult_seq2,rowMeans(ppd2),lwd=2)
polygon(c(adult_seq2,rev(adult_seq2)),c(ci2[1,],rev(ci2[2,])),col=adjustcolor("blue",0.3),border=NA)
dev.off()

# panel plots for covariates
layout(matrix(1:8,nrow=4,ncol=2,byrow=TRUE))

# survival and seals
ppy <- extract(fit)$pred_surv
ppx <- dat$x1[,1]
ciy <- apply(ppy,2,quantile,probs=c(0.05,0.95))
plot(ppx,colMeans(ppy),ylim=range(ciy),pch=21,bg="dodgerblue")
segments(x0=ppx,y0=ciy[1,],y1=ciy[2,],lwd=0.75)
points(ppx,colMeans(ppy),pch=21,bg="dodgerblue")

# adult N and survival
ppy <- extract(fit)$pred_adults
ppx <- extract(fit)$pred_surv
cix <- apply(ppx,2,quantile,probs=c(0.05,0.95))
ciy <- apply(ppy,2,quantile,probs=c(0.05,0.95))
plot(colMeans(ppx),colMeans(ppy),ylim=range(ciy),xlim=range(cix),pch=21,bg="dodgerblue")
segments(x0=colMeans(ppx),y0=ciy[1,],y1=ciy[2,],lwd=0.75)
segments(x0=cix[1,],x1=cix[2,],y0=colMeans(ppy),lwd=0.75)
points(colMeans(ppx),colMeans(ppy),pch=21,bg="dodgerblue")

# run timing and adult N
ppy <- extract(fit)$pred_run
ppx <- dat$x2[,2]
ciy <- apply(ppy,2,quantile,probs=c(0.05,0.95))
plot(ppx,colMeans(ppy),ylim=range(ciy),pch=21,bg="dodgerblue")
segments(x0=ppx,y0=ciy[1,],y1=ciy[2,],lwd=0.75)
points(ppx,colMeans(ppy),pch=21,bg="dodgerblue")

# run timing and rainfall
ppy <- extract(fit)$pred_run
ppx <- dat$x2[,1]
ciy <- apply(ppy,2,quantile,probs=c(0.05,0.95))
plot(ppx,colMeans(ppy),ylim=range(ciy),pch=21,bg="dodgerblue")
segments(x0=ppx,y0=ciy[1,],y1=ciy[2,],lwd=0.75)
points(ppx,colMeans(ppy),pch=21,bg="dodgerblue")

# run timing and smolt production
ppy <- extract(fit)$pred_rec
ppx <- extract(fit)$pred_run
cix <- apply(ppx,2,quantile,probs=c(0.05,0.95))
ciy <- apply(ppy,2,quantile,probs=c(0.05,0.95))
plot(colMeans(ppx),colMeans(ppy),ylim=range(ciy),xlim=range(cix),pch=21,bg="dodgerblue")
segments(x0=colMeans(ppx),y0=ciy[1,],y1=ciy[2,],lwd=0.75)
segments(x0=cix[1,],x1=cix[2,],y0=colMeans(ppy),lwd=0.75)
points(colMeans(ppx),colMeans(ppy),pch=21,bg="dodgerblue")

# mean temperature and smolt production
ppy <- extract(fit)$pred_rec
ppx <- dat$xx3[,1]
ciy <- apply(ppy,2,quantile,probs=c(0.05,0.95))
plot(ppx,colMeans(ppy),ylim=range(ciy),pch=21,bg="dodgerblue")
segments(x0=ppx,y0=ciy[1,],y1=ciy[2,],lwd=0.75)
points(ppx,colMeans(ppy),pch=21,bg="dodgerblue")

# total rainfall and smolt production
ppy <- extract(fit)$pred_rec
ppx <- dat$xx3[,2]
ciy <- apply(ppy,2,quantile,probs=c(0.05,0.95))
plot(ppx,colMeans(ppy),ylim=range(ciy),pch=21,bg="dodgerblue")
segments(x0=ppx,y0=ciy[1,],y1=ciy[2,],lwd=0.75)
points(ppx,colMeans(ppy),pch=21,bg="dodgerblue")

# posterior frequencies for coefficients
betas <- extract(fit)$beta_surv
colSums(betas>0)/nrow(betas)
betas <- extract(fit)$beta_run
colSums(betas>0)/nrow(betas)
betas <- extract(fit)$beta_run_cov
sum(betas>0)/length(betas)
betas <- extract(fit)$beta_rec_cov
colSums(betas>0)/nrow(betas)
