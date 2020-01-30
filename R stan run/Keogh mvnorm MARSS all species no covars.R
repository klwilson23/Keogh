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

xx1 <- scale(model.matrix(~-1+cumul_footprint+seals+oceanSalmon+freshPink+npgo+winTemp+winRain,data=keogh_long[keogh_long$Species=="Steelhead",]),center=TRUE,scale=TRUE)
xx2 <- scale(model.matrix(~-1+cumul_footprint+seals+oceanSalmon+freshPink+npgo+winTemp+winRain,data=keogh_long[keogh_long$Species=="Dolly Varden",]),center=TRUE,scale=TRUE)
xx3 <- scale(model.matrix(~-1+cumul_footprint+seals+oceanSalmon+freshPink+npgo+winTemp+winRain,data=keogh_long[keogh_long$Species=="Cutthroat",]),center=TRUE,scale=TRUE)
xx4 <- scale(model.matrix(~-1+cumul_footprint+seals+oceanSalmon+freshPink+npgo+winTemp+winRain,data=keogh_long[keogh_long$Species=="Pink",]),center=TRUE,scale=TRUE)
xx5 <- scale(model.matrix(~-1+cumul_footprint+seals+oceanSalmon+freshPink+npgo+winTemp+winRain,data=keogh_long[keogh_long$Species=="Coho",]),center=TRUE,scale=TRUE)


# all models
dat <- list("N"=nrow(x),
            "K"=ncol(x),
            "x"=x,
            "y"=y,
            "init_priors"=rep(-2e-3,ncol(x)))

fit <- stan(file="Stan code/Keogh mnorm MARSS no covars.stan",data=dat, iter=10000,chains=4,cores=4,control=list("adapt_delta"=0.9))

saveRDS(fit,file="~/Google Drive/SFU postdoc/Keogh river/Stan fits/keogh mvnorm.rds")

summary(fit,pars=c("beta","x0","L_Omega_obs","L_sigma_obs","L_Omega_proc","L_sigma_proc","Sigma_proc"),probs=c(0.025,0.975))$summary

pro_corr <- extract(fit)$Omega_proc
pro_corr_mn <- apply(pro_corr,c(3,2),mean)
colnames(pro_corr_mn) <- row.names(pro_corr_mn) <- unique(keogh_long$Species)

pro_cov <- extract(fit)$Sigma_proc
pro_cov_mn <- apply(pro_cov,c(3,2),mean)
colnames(pro_cov_mn) <- row.names(pro_cov_mn) <- unique(keogh_long$Species)

obs_corr <- extract(fit)$Omega
obs_corr_mn <- apply(obs_corr,c(3,2),mean)
colnames(obs_corr_mn) <- row.names(obs_corr_mn) <- unique(keogh_long$Species)

pdf("Figures/Keogh alpha productivity stan marss.pdf",width=6,height=6)
matLayout <- matrix(1:4,nrow=2,ncol=2,byrow=TRUE)
layout(matLayout)
par(mar=c(5,4,1,1))
ppd <- extract(fit)$x0
alpha_mns <- apply(ppd,c(3,2),mean)
for(k in 1:dat$K)
{
  for(j in (1:dat$K)[(1:dat$K)!=k])
  {
    plot(alpha_mns[k,],alpha_mns[j,],xlab=paste("alpha",unique(keogh_long$Species)[k]),ylab=paste("alpha",unique(keogh_long$Species)[j]),pch=21,bg="orange")
    Corner_text(paste(unique(keogh_long$Species)[k],"v.",unique(keogh_long$Species)[j]),"topleft")
  }
}
dev.off()

# productivity
jpeg("Figures/Keogh smolt productivity stan marss.jpeg",width=8,height=8,units="in",res=800)
matLayout <- matrix(c(1:5,0),nrow=3,ncol=2,byrow=TRUE)
layout(matLayout)
par(mar=c(5,4,1,1))
ppd <- extract(fit)$y_ppd
mns <- apply(ppd,c(3,2),mean)
ci <- apply(ppd,c(3,2),quantile,probs=c(0.05,0.95))
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
ci <- apply(ppd,c(3,2),quantile,probs=c(0.05,0.95))
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
ci <- apply(ppd,c(3,2),quantile,probs=c(0.05,0.95))
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