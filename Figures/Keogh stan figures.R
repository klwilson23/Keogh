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
library(corrplot)

fit <- readRDS(file="~/Google Drive/SFU postdoc/Keogh river/Stan fits/keogh mvnorm covars.rds")

keogh <- readRDS("Keogh_newRec_enviro.rds")
keogh_long <- subset(keogh,Year<=2015 & Year>=1976)
keogh_long <- subset(keogh_long,Species!="Chum")
keogh_long$Species <- factor(keogh_long$Species,levels=unique(keogh_long$Species))
keogh_long$marSurv <- keogh_long$Stock/keogh_long$juvCohort
keogh_long$logitSurv <- log(keogh_long$marSurv/(1-keogh_long$marSurv))
keogh_long$prod <- log(keogh_long$Recruits/keogh_long$Stock)

keogh_sel <- subset(keogh_long,select = c(Year,Species,Stock,prod))
keogh_wide <- reshape(keogh_sel,direction = "wide",idvar="Year",timevar="Species")

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

XXXvars <- c("meanLogging","total_rain_egg","mean_temp_egg","freshPink")
sdSurv_prod <- attr(scale(sh_annual[,XXXvars],center=TRUE,scale=TRUE),"scaled:scale")
mnSurv_prod <- attr(scale(sh_annual[,XXXvars],center=TRUE,scale=TRUE),"scaled:center")
enviro_prod <- scale(sh_annual[,XXXvars],center=TRUE,scale=TRUE)
enviro_prod <- data.frame(XXXvars=enviro_prod)
colnames(enviro_prod) <- XXXvars

x <- as.matrix(keogh_wide[,grep("Stock",colnames(keogh_wide))])
y <- as.matrix(keogh_wide[,grep("prod",colnames(keogh_wide))])
y[which(is.na(y),arr.ind=TRUE)] <- colMeans(y[1:10,],na.rm=TRUE)[which(is.na(y),arr.ind=TRUE)[,2]]

xx1 <- scale(model.matrix(~-1+meanLogging+total_rain_egg+mean_temp_egg+freshPink,data=enviro_prod),center=TRUE,scale=TRUE)
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
            "xx5"=xx5,
            "N_obs"=sum(!is.na(sh_annual$logit_surv)),
            "M"=ncol(x1),
            "Q"=ncol(x2),
            "P"=2,
            "x1"=x1,
            "x2"=x2,
            "juvCoh"=as.vector(scale(sh_annual$juvCohort)),
            #"x3"=as.numeric(x3),
            "y1_obs"=sh_annual$logit_surv[!is.na(sh_annual$logit_surv)],
            "y2"=sh_annual$run,
            "init_s0"=mean(sh_annual$logit_surv[1:10],na.rm=TRUE))

init_fx <- function(chain_id)
{
  list("beta_steel"=rep(0,dat$J1+1),
       "beta_dolly"=rep(0,dat$J2),
       "beta_cutt"=rep(0,dat$J3),
       "beta_pink"=rep(0,dat$J4),
       "beta_coho"=rep(0,dat$J5),
       "beta_adults"=rep(0,dat$P))
}

intervals <- c(0.1,0.9)

# posterior frequencies for coefficients
betas <- extract(fit)$beta_surv
colSums(betas>0)/nrow(betas)
betas <- extract(fit)$beta_adult
colSums(betas>0)/nrow(betas)
betas <- extract(fit)$beta_run
colSums(betas>0)/nrow(betas)
betas <- extract(fit)$beta_steel
colSums(betas>0)/nrow(betas)

pro_corr <- extract(fit)$Omega_proc
pro_corr_mn <- apply(pro_corr,c(3,2),mean)
colnames(pro_corr_mn) <- row.names(pro_corr_mn) <- unique(keogh_long$Species)

pro_cov <- extract(fit)$Sigma_proc
pro_cov_mn <- apply(pro_cov,c(3,2),mean)
colnames(pro_cov_mn) <- row.names(pro_cov_mn) <- unique(keogh_long$Species)

obs_cov <- extract(fit)$Sigma
obs_cov_mn <- apply(obs_cov,c(3,2),mean)
colnames(obs_cov_mn) <- row.names(obs_cov_mn) <- unique(keogh_long$Species)

obs_corr <- extract(fit)$Omega
obs_corr_mn <- apply(obs_corr,c(3,2),mean)
colnames(obs_corr_mn) <- row.names(obs_corr_mn) <- unique(keogh_long$Species)

pdf("Figures/Keogh species interactions.pdf",width=6,height=6)
matLayout <- matrix(1:4,nrow=2,ncol=2,byrow=TRUE)
layout(matLayout)
par(mar=c(5,4,1,1))
ppd <- extract(fit)$x0
alpha_mns <- apply(ppd,c(3,2),mean)
r <- matrix(NA,nrow=dat$K,ncol=dat$K,dimnames=list(unique(keogh_long$Species),unique(keogh_long$Species)))
diag(r) <- 1
ci <- apply(ppd,c(3,2),quantile,probs=intervals)
for(k in 1:dat$K)
{
  for(j in (1:dat$K)[(1:dat$K)!=k])
  {
    r[k,j] <- mean(apply(ppd,1,function(x){cor(x[,k],x[,j])}))
    plot(alpha_mns[k,],alpha_mns[j,],xlab=paste("alpha",unique(keogh_long$Species)[k]),ylab=paste("alpha",unique(keogh_long$Species)[j]),pch=21,bg="orange",ylim=range(ci[,j,]),xlim=range(ci[,k,]))
    segments(x0=ci[1,k,],x1=ci[2,k,],y0=alpha_mns[j,],lwd=0.5)
    segments(x0=alpha_mns[k,],y0=ci[1,j,],y1=ci[2,j,],lwd=0.5)
    points(alpha_mns[k,],alpha_mns[j,],pch=21,bg="orange")
    Corner_text(paste(unique(keogh_long$Species)[k],"v.",unique(keogh_long$Species)[j]),"topleft")
  }
}
dev.off()

jpeg("Figures/Keogh species correlations.jpeg",width=8,height=8,units="in",res=800)
layout(1)
par(mar=c(1,1,1,1),mai=c(0.2,0.2,0.2,0.2))
colFun <- colorRampPalette(c("tomato","orange","dodgerblue","darkblue"))
corrplot.mixed(r,upper="ellipse",lower.col="black",tl.col="black",upper.col = colFun(100))
dev.off()

# productivity
jpeg("Figures/Keogh productivity stan marss.jpeg",width=8,height=8,units="in",res=800)
matLayout <- matrix(c(1:5,0),nrow=3,ncol=2,byrow=TRUE)
layout(matLayout)
par(mar=c(5,4,1,1))
ppd <- extract(fit)$y_ppd
mns <- apply(ppd,c(3,2),median)
ci <- apply(ppd,c(3,2),quantile,probs=intervals)
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
regime <- which(sh_annual==1991)
regime2 <- which(sh_annual==2008)
matLayout <- matrix(0,nrow=10,ncol=10,byrow=TRUE)
matLayout[1:2,1:5] <- 1
matLayout[3:4,1:5] <- 2
matLayout[5:6,1:5] <- 3
matLayout[7:8,1:5] <- 4
matLayout[9:10,1:5] <- 5
matLayout[1:6,6:10] <- 6
matLayout[7:10,6:10] <- 7
layout(matLayout)
par(mar=c(5,4,1,1),mai=c(0.5,0.5,0.05,0.05))
ppd <- extract(fit)$x0
mns <- apply(ppd,c(3,2),median)
ci <- apply(ppd,c(3,2),quantile,probs=intervals)
for(k in 1:dat$K)
{
  plot(mns[k,],ylim=range(ci[,k,]),type="l",lwd=2,xlab="",ylab="",xaxt="n")
  axis(1,labels=FALSE,tick=TRUE)
  mtext("Productivity (ln \u03B1)",side=2,line=2.25,cex=ptSize)
  polygon(c(1:regime,rev(1:regime)),c(ci[1,k,1:regime],rev(ci[2,k,1:regime])),col=adjustcolor("orange",0.5),border=NA)
  polygon(c(regime:regime2,rev(regime:regime2)),c(ci[1,k,regime:regime2],rev(ci[2,k,regime:regime2])),col=adjustcolor("dodgerblue",0.5),border=NA)
  polygon(c(regime2:40,rev(regime2:40)),c(ci[1,k,regime2:40],rev(ci[2,k,regime2:40])),col=adjustcolor("darkblue",0.5),border=NA)
  lines(mns[k,],lwd=2,col="black")
  Corner_text(paste("(a) - ",unique(keogh_long$Species)[k],sep=""),"topleft")
}
mtext("Year",side=1,line=2,cex=ptSize)
axis(1,at=seq(from=0,to=dat$N,length=5),labels=seq(from=min(sh_annual$year)-1,to=max(sh_annual$year),length=5),tick=TRUE)
legend("bottomright",c("Early","Compensation","Decline"),title="Regime",pch=22,pt.bg=c(adjustcolor("orange",0.5),adjustcolor("dodgerblue",0.5),adjustcolor("darkblue",0.5)),bty="n")

# covariates
par(mar=c(5,2,1,1),mai=c(0.5,2,0.05,0.05))
plot(rep(0,sum(dat$J1+1+dat$J2+dat$J3+dat$J4+dat$J5)),1:sum(dat$J1+1+dat$J2+dat$J3+dat$J4+dat$J5),xlim=c(-4,4),type="p",pch=0,ylab="",xlab="",bg=0,col=0,yaxt="n")
axis(2,at=1:sum(dat$J1+1+dat$J2+dat$J3+dat$J4+dat$J5),labels=FALSE,tick=TRUE)
abline(h=1:sum(dat$J1+1+dat$J2+dat$J3+dat$J4+dat$J5),lwd=0.5,lty=3,col="grey85")
abline(v=0,lty=1,lwd=1,col="red")
mtext("Effect size",side=1,line=2,cex=ptSize)

counter <- sum(dat$J1+1+dat$J2+dat$J3+dat$J4+dat$J5)
for(k in 1:dat$K)
{
  xnames <- c("Cumulative logging","Summer air temperature","Summer rainfall","Winter air temperature","Winter rainfall","Pink salmon abundance in FW","Adult run time")
  if(k==1) { betas <- extract(fit)$beta_steel; xx <- data.frame(xx1,runtime);
  xnames <- c("Cumulative logging","Rainfall after run","Temp. after run","Pink salmon in FW","Spawn date")}
  if(k==2) { betas <- extract(fit)$beta_dolly; xx <- xx2}
  if(k==3) { betas <- extract(fit)$beta_cutt; xx <- xx3}
  if(k==4) { betas <- extract(fit)$beta_pink; xx <- xx4;
  xnames <- c("Cumulative logging","Summer air temperature","Summer rainfall","Seal densities","Pacific salmon biomass","NPGO","Winter air temperature","Winter rainfall")}
  if(k==5) { betas <- extract(fit)$beta_coho; xx <- xx5}
  
  par(mar=c(5,4,1,1))
  for(j in 1:ncol(betas))
  {
    if(j==1){text(-4,y=counter+0.5,paste("(b) - ", unique(keogh_long$Species)[k],sep=""),cex=ptSize,adj=0,xpd=NA)}
    points(mean(betas[,j])/sd(dat$y[,k]),counter,type="p",pch=21,bg="dodgerblue",xpd=NA)
    segments(x0=quantile(betas[,j]/sd(dat$y[,k]),probs=intervals[1]),x1=quantile(betas[,j]/sd(dat$y[,k]),probs=intervals[2]),y0=counter,lwd=2,xpd=NA)
    axis(2,at=counter,labels=xnames[j],tick=FALSE,las=1,line=-10.75,cex=ptSize)
    ifelse((quantile(betas[,j]/sd(dat$y[,k]),probs=intervals[1])*quantile(betas[,j]/sd(dat$y[,k]),probs=intervals[2]))>0,pch_color <- "red",pch_color<-"grey50")
    points(mean(betas[,j])/sd(dat$y[,k]),counter,pch=21,bg=pch_color,xpd=NA)
    counter <- counter - 1
  } 
}

par(mar=c(1,1,1,1),mai=c(0.3,0,0.2,0.05))
colFun <- colorRampPalette(c("tomato","orange","dodgerblue","darkblue"))
corrplot.mixed(r,upper="ellipse",lower.col="black",tl.col="black",upper.col = colFun(100),tl.cex = ptSize)
text("(c) - Species correlations",x=1,y=5.75,xpd=NA)
dev.off()


# recruitment
jpeg("Figures/Keogh recruitment stan marss.jpeg",width=8,height=8,units="in",res=800)
matLayout <- matrix(c(1:5,0),nrow=3,ncol=2,byrow=TRUE)
layout(matLayout)
par(mar=c(5,4,1,1))
ppd <- extract(fit)$R
mns <- apply(ppd,c(3,2),median)
ci <- apply(ppd,c(3,2),quantile,probs=intervals)
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
ci <- apply(ppd,c(3,2),quantile,probs=intervals)
mus <- extract(fit)$mu
runtime <- extract(fit)$pred_run
runtime <- apply(runtime,1,FUN=function(x){(x-mean(x))/sd(x)})
runtime <- rowMeans(runtime)
for(k in 1:dat$K)
{
  xnames <- c("Cumulative logging (15-years lagged)","Summer air temperature","Summer rainfall","Winter air temperature","Winter rainfall","Pink salmon abundance in early life","Adult run time")
  if(k==1) { betas <- extract(fit)$beta_steel; xx <- data.frame(xx1,runtime);
  xnames <- c("Cumulative logging (15-years lagged)","Rainfall after run","Temp. after run","Pink salmon in FW","Spawn date")}
  if(k==2) { betas <- extract(fit)$beta_dolly; xx <- xx2}
  if(k==3) { betas <- extract(fit)$beta_cutt; xx <- xx3}
  if(k==4) { betas <- extract(fit)$beta_pink; xx <- xx4;
  xnames <- c("Cumulative logging (15-years lagged)","Summer air temperature","Summer rainfall","Seal densities","Pacific salmon biomass","NPGO","Winter air temperature","Winter rainfall")}
  if(k==5) { betas <- extract(fit)$beta_coho; xx <- xx5}
  matLayout <- matrix(c(1:ncol(xx),rep(0,8-ncol(xx))),nrow=4,ncol=2,byrow=TRUE)
  layout(matLayout)
  par(mar=c(5,4,1,1))
  for(j in 1:ncol(betas))
  {
    covar_seq <- seq(from=min(xx[,j]),to=max(xx[,j]),length=25)
    resid_alphas <- mus[,,k] - dat$x[,k]*extract(fit)$beta[,k] - ppd[,,k] - apply(xx[,-j],1,function(x){rowSums(x*betas[,-j])})
    pred <- apply(resid_alphas,2,quantile,probs=intervals)
    plot(xx[,j],colMeans(resid_alphas),ylim=range(pred,na.rm=TRUE),type="p",pch=21,xlab=xnames[j],ylab="Residual productivity (ln R/S)",xaxt="n",bg="grey50")
    segments(x0=xx[,j],y0=pred[1,],y1=pred[2,],lwd=0.5)
    regress <- sapply(covar_seq,function(x){x*betas[,j]})
    ci <- apply(regress,2,quantile,probs=intervals)
    axis(1,labels=NULL,tick=TRUE)
    polygon(c(covar_seq,rev(covar_seq)),c(ci[1,],rev(ci[2,])),col=adjustcolor("dodgerblue",0.3),border=NA)
    points(xx[,j],colMeans(resid_alphas),pch=21,bg="grey50")
    Corner_text(unique(keogh_long$Species)[k],"topleft")
  }
}
dev.off()

# steelhead lifecycle
# survival
jpeg("Figures/Steelhead cycle marss.jpeg",res=800,height=8,width=8,units="in")
regime <- which(sh_annual==1991)
regime2 <- which(sh_annual==2008)
matLayout <- matrix(0,nrow=15,ncol=15)
matLayout[1:3,1:9] <- 1
matLayout[1:3,10:15] <- 2

matLayout[4:6,1:9] <- 3
matLayout[4:6,10:15] <- 4

matLayout[7:9,1:9] <- 5
matLayout[7:9,10:15] <- 6

matLayout[10:12,1:9] <- 7
matLayout[10:12,10:15] <- 8

matLayout[13:15,1:9] <- 9
matLayout[13:15,10:15] <- 10

ptSize <- 0.7

layout(matLayout)
par(mar=c(3.5,4,0,0.5),mai=c(0.45,0.55,0.05,0.05))
ppd <- extract(fit)$y1_ppd
ci <- apply(ppd,2,quantile,probs=intervals)
plot(colMeans(ppd),ylim=range(ci,quantile(1/(1+exp(-extract(fit)$y1_miss)),probs=intervals)),type="l",lwd=2,xlab="",ylab="",xaxt="n")
mtext("Marine survival (yr-1)",2,line=2.5,xpd=NA,cex=ptSize)
axis(1,tick=TRUE,labels=FALSE)
polygon(c(1:regime,rev(1:regime)),c(ci[1,1:regime],rev(ci[2,1:regime])),col=adjustcolor("orange",0.5),border=NA)
polygon(c(regime:regime2,rev(regime:regime2)),c(ci[1,regime:regime2],rev(ci[2,regime:regime2])),col=adjustcolor("dodgerblue",0.5),border=NA)
polygon(c(regime2:40,rev(regime2:40)),c(ci[1,regime2:40],rev(ci[2,regime2:40])),col=adjustcolor("darkblue",0.5),border=NA)
lines(colMeans(ppd),lwd=2)
ppx <- extract(fit)$y1
ci <- apply((1/(1+exp(-ppx))),2,quantile,probs=intervals)
segments(x0=1:ncol(ppx),y0=ci[1,],y1=ci[2,],lwd=2,col="black")
points(colMeans(1/(1+exp(-ppx))),pch=21,bg="grey50")
legend("topright",c("Early","Compensation","Decline"),title="Regime",pch=22,pt.bg=c(adjustcolor("orange",0.5),adjustcolor("dodgerblue",0.5),adjustcolor("darkblue",0.5)),bty="n",pt.cex = 1.2)
Corner_text("(a)","topleft")
# some effects on survival

betas <- extract(fit)$beta_surv
std_eff <- apply(betas,2,function(x){c(quantile(x,probs=intervals[1]),mean(x),quantile(x,probs=intervals[2]))})/sd(dat$y1_obs)
colnames(std_eff) <- c("Seal densities","NPGO","NP Salmon biomass")
plot(std_eff[2,],ylim=1.1*range(c(0,std_eff)),xaxt="n",xlim=range(0.5:(ncol(std_eff)+0.5)),yaxt="n",ylab="",col=0,xlab="")
axis(1,at=1:ncol(std_eff),labels=colnames(std_eff),cex.axis=1)
axis(2,line=0)
mtext("Effect size",side=2,cex=ptSize,line=2.5)
segments(x0=1:3,y0=std_eff[1,],y1=std_eff[3,],lwd=2)
points(std_eff[2,],pch=21,bg="dodgerblue",cex=1.5)
abline(h=0,lwd=1,lty=5)
Corner_text("(b)","topleft")

# adult abundance
ppd <- extract(fit)$x3_ppd
ci <- apply(ppd,2,quantile,probs=intervals)
colMed <- apply(ppd,2,quantile,probs=0.5)
plot(colMed,ylim=range(ci),type="l",lwd=2,xlab="",ylab="",xaxt="n")
mtext("Adult female returns",2,line=2.5,xpd=NA,cex=ptSize)
axis(1,tick=TRUE,labels=FALSE)
polygon(c(1:regime,rev(1:regime)),c(ci[1,1:regime],rev(ci[2,1:regime])),col=adjustcolor("orange",0.5),border=NA)
polygon(c(regime:regime2,rev(regime:regime2)),c(ci[1,regime:regime2],rev(ci[2,regime:regime2])),col=adjustcolor("dodgerblue",0.5),border=NA)
polygon(c(regime2:40,rev(regime2:40)),c(ci[1,regime2:40],rev(ci[2,regime2:40])),col=adjustcolor("darkblue",0.5),border=NA)
lines(colMed,lwd=2)
points(sh_annual$Stock,pch=21,bg="grey50")
Corner_text("(c)","topleft")

# some effects on adult abundance
betas <- extract(fit)$beta_adults
std_eff <- apply(betas,2,function(x){c(quantile(x,probs=intervals[1]),mean(x),quantile(x,probs=intervals[2]))})/log(sd(dat$x[,1]))
colnames(std_eff) <- c("Marine survival","Smolt cohort size")
plot(std_eff[2,],ylim=1.1*range(c(0,std_eff)),xaxt="n",xlim=range(0.5:(ncol(std_eff)+0.5)),yaxt="n",ylab="",col=0,xlab="")
axis(1,at=1:ncol(std_eff),labels=colnames(std_eff),cex.axis=1)
axis(2,line=0)
mtext("Effect size",side=2,cex=ptSize,line=2.5)
segments(x0=1:3,y0=std_eff[1,],y1=std_eff[3,],lwd=2)
points(std_eff[2,],pch=21,bg="dodgerblue",cex=1.5)
abline(h=0,lwd=1,lty=5)
Corner_text("(d)","topleft")

# run time
ppd <- extract(fit)$y2_ppd
ci <- apply(ppd,2,quantile,probs=intervals)
plot(colMeans(ppd),ylim=range(ci),type="l",lwd=2,xlab="",ylab="",xaxt="n")
mtext("Spawning run (day of season)",2,line=2.5,xpd=NA,cex=ptSize)
axis(1,tick=TRUE,labels=FALSE)
polygon(c(1:regime,rev(1:regime)),c(ci[1,1:regime],rev(ci[2,1:regime])),col=adjustcolor("orange",0.5),border=NA)
polygon(c(regime:regime2,rev(regime:regime2)),c(ci[1,regime:regime2],rev(ci[2,regime:regime2])),col=adjustcolor("dodgerblue",0.5),border=NA)
polygon(c(regime2:40,rev(regime2:40)),c(ci[1,regime2:40],rev(ci[2,regime2:40])),col=adjustcolor("darkblue",0.5),border=NA)
lines(colMeans(ppd),lwd=2)
points(dat$y2,pch=21,bg="grey50")
Corner_text("(e)","topleft")

# some effects on run time
betas <- extract(fit)$beta_run
std_eff <- apply(betas,2,function(x){c(quantile(x,probs=intervals[1]),mean(x),quantile(x,probs=intervals[2]))})/sd(dat$y2)
colnames(std_eff) <- c("ln(Adults)","Rainfall before run","Temp. before run")
plot(std_eff[2,],ylim=1.1*range(c(0,std_eff)),xaxt="n",xlim=range(0.5:(ncol(std_eff)+0.5)),yaxt="n",ylab="",col=0,xlab="")
axis(1,at=1:ncol(std_eff),labels=colnames(std_eff),cex.axis=0.8)
axis(2,line=0)
mtext("Effect size",side=2,cex=ptSize,line=2.5)
segments(x0=1:3,y0=std_eff[1,],y1=std_eff[3,],lwd=2)
points(std_eff[2,],pch=21,bg="dodgerblue",cex=1.5)
abline(h=0,lwd=1,lty=5)
Corner_text("(f)","topleft")

# productivity
ppd <- extract(fit)$y_ppd[,,1]
ci <- apply(ppd,2,quantile,probs=intervals)
plot(colMeans(ppd),ylim=range(ci),type="l",lwd=2,xlab="",ylab="",xaxt="n")
mtext("Smolt productivity (ln R/S)",2,line=2.5,xpd=NA,cex=ptSize)
axis(1,labels=FALSE,tick=TRUE)
polygon(c(1:regime,rev(1:regime)),c(ci[1,1:regime],rev(ci[2,1:regime])),col=adjustcolor("orange",0.5),border=NA)
polygon(c(regime:regime2,rev(regime:regime2)),c(ci[1,regime:regime2],rev(ci[2,regime:regime2])),col=adjustcolor("dodgerblue",0.5),border=NA)
polygon(c(regime2:40,rev(regime2:40)),c(ci[1,regime2:40],rev(ci[2,regime2:40])),col=adjustcolor("darkblue",0.5),border=NA)
lines(colMeans(ppd),lwd=2)
points(dat$y[,1],pch=21,bg="grey50")
Corner_text("(g)","topleft")

# some effects on productivity
betas <- extract(fit)$beta_steel
std_eff <- apply(betas,2,function(x){c(quantile(x,probs=intervals[1]),mean(x),quantile(x,probs=intervals[2]))})/sd(dat$y[,1])
#colnames(std_eff) <- c("Cumulative logging","Summer temp.","Summer rain","Winter temp.","Winter rain","Pink salmon in early life","Spawn date")
colnames(std_eff) <- c("Logging","Rainfall","Temper.","Pink salmon","Spawn date")
plot(std_eff[2,],ylim=1.1*range(c(0,std_eff)),xaxt="n",xlim=range(0.5:(ncol(std_eff)+0.5)),yaxt="n",ylab="",col=0,xlab="")
axis(1,at=1:ncol(std_eff),labels=colnames(std_eff),cex.axis=ptSize)
axis(2,line=0)
mtext("Effect size",side=2,cex=ptSize,line=2.5)
segments(x0=1:length(colnames(std_eff)),y0=std_eff[1,],y1=std_eff[3,],lwd=2)
points(std_eff[2,],pch=21,bg="dodgerblue",cex=1.5)
abline(h=0,lwd=1,lty=5)
Corner_text("(h)","topleft")

# recruitment
ppd <- extract(fit)$R[,,1]
ci <- apply(ppd,2,quantile,probs=intervals)
plot(colMeans(ppd),ylim=range(ci),type="l",lwd=2,xlab="",ylab="",xaxt="n")
mtext("Smolt recruitment",2,line=2.5,xpd=NA,cex=ptSize)
axis(1,at=seq(from=0,to=dat$N,length=5),labels=seq(from=min(sh_annual$year)-1,to=max(sh_annual$year),length=5),tick=TRUE)
mtext("Year",side=1,line=2.5,cex=ptSize)
polygon(c(1:regime,rev(1:regime)),c(ci[1,1:regime],rev(ci[2,1:regime])),col=adjustcolor("orange",0.5),border=NA)
polygon(c(regime:regime2,rev(regime:regime2)),c(ci[1,regime:regime2],rev(ci[2,regime:regime2])),col=adjustcolor("dodgerblue",0.5),border=NA)
polygon(c(regime2:40,rev(regime2:40)),c(ci[1,regime2:40],rev(ci[2,regime2:40])),col=adjustcolor("darkblue",0.5),border=NA)
lines(colMeans(ppd),lwd=2)
points(sh_annual$Recruits,pch=21,bg="grey50")
Corner_text("(i)","topleft")

# some effects on recruitment
betas <- extract(fit)$beta[,1]
alphas <- extract(fit)$x0[,,1]
runtime <- extract(fit)$pred_run
runtime <- t(apply(runtime,1,FUN=function(x){(x-mean(x))/sd(x)}))
runtime <- colMeans(runtime)
beta_cov <- colMeans(extract(fit)$beta_steel)
xx <- data.frame(xx1,runtime)

adult_seq <- seq(from=0,to=max(dat$x[1:(regime-1),1]),length.out=25)
adult_seq2 <- seq(from=0,to=max(dat$x[regime:nrow(dat$x),1]),length.out=25)
ppd <- sapply(1:length(betas),function(x){exp(mean(alphas[x,1:(regime-1)]))*adult_seq*exp(betas[x] * adult_seq + sum(beta_cov * colMeans(xx[1:(regime-1),])))})
ci <- apply(ppd,1,quantile,probs=intervals)
ppd2 <- sapply(1:length(betas),function(x){exp(mean(alphas[x,regime:nrow(dat$x)]))*adult_seq2*exp(betas[x] * adult_seq2 + sum(beta_cov * colMeans(xx[regime:nrow(dat$x),])))})
ci2 <- apply(ppd2,1,quantile,probs=intervals)
plot(adult_seq,rowMeans(ppd),xlab="",type="l",ylab="",ylim=range(ci,ci2),lwd=2,xaxt="n",yaxt="n")
axis(1,line=0)
axis(2,line=0)
mtext("Smolt recruitment",side=2,line=2.5,cex=ptSize)
mtext("Adult steelhead",side=1,line=2.5,cex=ptSize)
polygon(c(adult_seq,rev(adult_seq)),c(ci[1,],rev(ci[2,])),col=adjustcolor("orange",0.4),border=NA)
lines(adult_seq2,rowMeans(ppd2),lwd=2)
polygon(c(adult_seq2,rev(adult_seq2)),c(ci2[1,],rev(ci2[2,])),col=adjustcolor("dodgerblue",0.4),border=NA)
Corner_text("(j)","topleft")

dev.off()

sh_adults <- extract(fit)$pred_adults
ad <- colMeans(sh_adults)
plot(log(ad),dat$y2)
plot(log(x[,1]),dat$y2)
