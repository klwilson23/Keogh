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

keogh <- readRDS("Keogh_newRec_enviro.rds")
keogh_long <- subset(keogh,Year<=2015 & Year>=1976)
keogh_long <- subset(keogh_long,Species!="Chum")
keogh_long$Species <- factor(keogh_long$Species,levels=unique(keogh_long$Species))
keogh_long$marSurv <- keogh_long$Stock/keogh_long$juvCohort
keogh_long$logitSurv <- log(keogh_long$marSurv/(1-keogh_long$marSurv))
keogh_long$prod <- log(keogh_long$Recruits/keogh_long$Stock)

keogh_sel <- subset(keogh_long,select = c(Year,Species,Stock,prod))
keogh_wide <- reshape(keogh_sel,direction = "wide",idvar="Year",timevar="Species")

x <- as.matrix(keogh_wide[,grep("Stock",colnames(keogh_wide))])
y <- as.matrix(keogh_wide[,grep("prod",colnames(keogh_wide))])
y[which(is.na(y),arr.ind=TRUE)] <- colMeans(y[1:10,],na.rm=TRUE)[which(is.na(y),arr.ind=TRUE)[,2]]

xx1 <- scale(model.matrix(~-1+meanLogging+sumTemp+sumRain+winTemp+winRain+freshPink,data=keogh_long[keogh_long$Species=="Steelhead",]),center=TRUE,scale=TRUE)
xx2 <- scale(model.matrix(~-1+meanLogging+sumTemp+sumRain+winTemp+winRain+freshPink,data=keogh_long[keogh_long$Species=="Dolly Varden",]),center=TRUE,scale=TRUE)
xx3 <- scale(model.matrix(~-1+meanLogging+sumTemp+sumRain+winTemp+winRain+freshPink,data=keogh_long[keogh_long$Species=="Cutthroat",]),center=TRUE,scale=TRUE)
xx4 <- scale(model.matrix(~-1+meanLogging+sumTemp+sumRain+seals+oceanSalmon+npgo+winTemp+winRain,data=keogh_long[keogh_long$Species=="Pink",]),center=TRUE,scale=TRUE)
xx5 <- scale(model.matrix(~-1+meanLogging+sumTemp+sumRain+winTemp+winRain+freshPink,data=keogh_long[keogh_long$Species=="Coho",]),center=TRUE,scale=TRUE)


# all models
dat <- list("N"=nrow(x),
            "K"=ncol(x),
            "x"=x,
            "y"=y,
            "init_priors"=rep(-2e-3,ncol(x)),
            "J1"=ncol(xx1),
            "J2"=ncol(xx2),
            "J3"=ncol(xx3),
            "J4"=ncol(xx4),
            "J5"=ncol(xx5),
            "xx1"=xx1,
            "xx2"=xx2,
            "xx3"=xx3,
            "xx4"=xx4,
            "xx5"=xx5)
init_fx <- function(chain_id)
{
  list("beta_steel"=rep(0,dat$J1),
       "beta_dolly"=rep(0,dat$J2),
       "beta_cutt"=rep(0,dat$J3),
       "beta_pink"=rep(0,dat$J4),
       "beta_coho"=rep(0,dat$J5))
}

fit <- stan(file="Stan code/Keogh mnorm MARSS.stan",data=dat, iter=10000,chains=6,cores=6,control=list("adapt_delta"=0.95,"max_treedepth"=15),init=init_fx)

saveRDS(fit,file="~/Google Drive/SFU postdoc/Keogh river/Stan fits/keogh mvnorm covars.rds")

summary(fit,pars=c("beta","beta_steel","beta_dolly","beta_cutt","beta_pink","beta_coho","x0","L_Omega_obs","L_sigma_obs","L_Omega_proc","L_sigma_proc","Sigma_proc"),probs=c(0.025,0.975))$summary

pro_corr <- extract(fit)$Omega_proc
pro_corr_mn <- apply(pro_corr,c(3,2),mean)
colnames(pro_corr_mn) <- row.names(pro_corr_mn) <- unique(keogh_long$Species)

pro_cov <- extract(fit)$Sigma_proc
pro_cov_mn <- apply(pro_cov,c(3,2),mean)
colnames(pro_cov_mn) <- row.names(pro_cov_mn) <- unique(keogh_long$Species)

obs_corr <- extract(fit)$Omega
obs_corr_mn <- apply(obs_corr,c(3,2),mean)
colnames(obs_corr_mn) <- row.names(obs_corr_mn) <- unique(keogh_long$Species)

pdf("Figures/Keogh species interactions.pdf",width=6,height=6)
matLayout <- matrix(1:4,nrow=2,ncol=2,byrow=TRUE)
layout(matLayout)
par(mar=c(5,4,1,1))
ppd <- extract(fit)$x0
alpha_mns <- apply(ppd,c(3,2),mean)
ci <- apply(ppd,c(3,2),quantile,probs=c(0.025,0.975))
for(k in 1:dat$K)
{
  for(j in (1:dat$K)[(1:dat$K)!=k])
  {
    plot(alpha_mns[k,],alpha_mns[j,],xlab=paste("alpha",unique(keogh_long$Species)[k]),ylab=paste("alpha",unique(keogh_long$Species)[j]),pch=21,bg="orange",ylim=range(ci[,j,]),xlim=range(ci[,k,]))
    segments(x0=ci[1,k,],x1=ci[2,k,],y0=alpha_mns[j,],lwd=0.5)
    segments(x0=alpha_mns[k,],y0=ci[1,j,],y1=ci[2,j,],lwd=0.5)
    points(alpha_mns[k,],alpha_mns[j,],pch=21,bg="orange")
    Corner_text(paste(unique(keogh_long$Species)[k],"v.",unique(keogh_long$Species)[j]),"topleft")
  }
}
dev.off()

# productivity
jpeg("Figures/Keogh productivity stan marss.jpeg",width=8,height=8,units="in",res=800)
matLayout <- matrix(c(1:5,0),nrow=3,ncol=2,byrow=TRUE)
layout(matLayout)
par(mar=c(5,4,1,1))
ppd <- extract(fit)$y_ppd
mns <- apply(ppd,c(3,2),mean)
ci <- apply(ppd,c(3,2),quantile,probs=c(0.025,0.975))
for(k in 1:dat$K)
{
  plot(mns[k,],ylim=range(ci[,k,]),type="l",lwd=2,xlab="Year",ylab="Smolt productivity (ln R/S)",xaxt="n")
  axis(1,labels=NULL,tick=TRUE)
  polygon(c(1:dat$N,rev(1:dat$N)),c(ci[1,k,],rev(ci[2,k,])),col=adjustcolor("dodgerblue",0.3),border=NA)
  points(dat$y[,k],pch=21,bg="orange")
  Corner_text(unique(keogh_long$Species)[k],"topleft")
}
dev.off()

# productivity
jpeg("Figures/Keogh alpha stan marss.jpeg",width=8,height=8,units="in",res=800)
matLayout <- matrix(c(1:5,0),nrow=3,ncol=2,byrow=TRUE)
layout(matLayout)
par(mar=c(5,4,1,1))
ppd <- extract(fit)$x0
mns <- apply(ppd,c(3,2),mean)
ci <- apply(ppd,c(3,2),quantile,probs=c(0.025,0.975))
for(k in 1:dat$K)
{
  plot(mns[k,],ylim=range(ci[,k,]),type="l",lwd=2,xlab="Year",ylab="Productivity (ln \u03B1)",xaxt="n")
  axis(1,labels=NULL,tick=TRUE)
  polygon(c(1:dat$N,rev(1:dat$N)),c(ci[1,k,],rev(ci[2,k,])),col=adjustcolor("dodgerblue",0.3),border=NA)
  Corner_text(unique(keogh_long$Species)[k],"topleft")
}
dev.off()


# recruitment
jpeg("Figures/Keogh recruitment stan marss.jpeg",width=8,height=8,units="in",res=800)
matLayout <- matrix(c(1:5,0),nrow=3,ncol=2,byrow=TRUE)
layout(matLayout)
par(mar=c(5,4,1,1))
ppd <- extract(fit)$R
mns <- apply(ppd,c(3,2),mean)
ci <- apply(ppd,c(3,2),quantile,probs=c(0.025,0.975))
for(k in 1:dat$K)
{
  rec <- keogh_long[keogh_long$Species==unique(keogh_long$Species)[k],"Recruits"]
  plot(mns[k,],ylim=range(ci[,k,],rec,na.rm=TRUE),type="l",lwd=2,xlab="Year",ylab="Recruitment",xaxt="n")
  axis(1,labels=NULL,tick=TRUE)
  polygon(c(1:dat$N,rev(1:dat$N)),c(ci[1,k,],rev(ci[2,k,])),col=adjustcolor("dodgerblue",0.3),border=NA)
  points(rec,pch=21,bg="orange")
  Corner_text(unique(keogh_long$Species)[k],"topleft")
}
dev.off()

pdf("Figures/Keogh regression stan marss.pdf",width=5.5,height=8)
ppd <- extract(fit)$x0
mns <- apply(ppd,c(3,2),mean)
ci <- apply(ppd,c(3,2),quantile,probs=c(0.025,0.975))
mus <- extract(fit)$mu
for(k in 1:dat$K)
{
  xnames <- c("Cumulative logging (15-years lagged)","Summer air temperature","Summer rainfall","Winter air temperature","Winter rainfall","Pink salmon abundance in early life")
  if(k==1) { betas <- extract(fit)$beta_steel; xx <- xx1}
  if(k==2) { betas <- extract(fit)$beta_dolly; xx <- xx2}
  if(k==3) { betas <- extract(fit)$beta_cutt; xx <- xx3}
  if(k==4) { betas <- extract(fit)$beta_pink; xx <- xx4;
  xnames <- xnames <- c("Cumulative logging (15-years lagged)","Summer air temperature","Summer rainfall","Seal densities","Pacific salmon biomass","NPGO","Winter air temperature","Winter rainfall","Pink salmon abundance in early life")}
  if(k==5) { betas <- extract(fit)$beta_coho; xx <- xx5}
  matLayout <- matrix(c(1:ncol(xx),rep(0,8-ncol(xx))),nrow=4,ncol=2,byrow=TRUE)
  layout(matLayout)
  par(mar=c(5,4,1,1))
  for(j in 1:ncol(betas))
  {
    covar_seq <- seq(from=min(xx[,j]),to=max(xx[,j]),length=25)
    resid_alphas <- mus[,,k] - dat$x[,k]*extract(fit)$beta[,k] - ppd[,,k] - apply(xx[,-j],1,function(x){rowSums(x*betas[,-j])})
    pred <- apply(resid_alphas,2,quantile,probs=c(0.025,0.975))
    plot(xx[,j],colMeans(resid_alphas),ylim=range(pred,na.rm=TRUE),type="p",pch=21,xlab=xnames[j],ylab="Residual productivity (ln R/S)",xaxt="n",bg="grey50")
    segments(x0=xx[,j],y0=pred[1,],y1=pred[2,],lwd=0.5)
    regress <- sapply(covar_seq,function(x){x*betas[,j]})
    ci <- apply(regress,2,quantile,probs=c(0.025,0.975))
    axis(1,labels=NULL,tick=TRUE)
    polygon(c(covar_seq,rev(covar_seq)),c(ci[1,],rev(ci[2,])),col=adjustcolor("dodgerblue",0.3),border=NA)
    points(xx[,j],colMeans(resid_alphas),pch=21,bg="grey50")
    Corner_text(unique(keogh_long$Species)[k],"topleft")
  }
}
dev.off()