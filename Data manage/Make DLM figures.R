source("some functions.R")
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(MARSS)
library(broom)
library(wesanderson)
keoghDLM <- readRDS("Results/keoghDLM.rds")
wesAnderson <- "Darjeeling1"

keogh <- readRDS("Keogh_newJuv_enviro.rds")
keogh_long <- subset(keogh,Year<=2015 & Year>=1976)
keogh_long <- subset(keogh_long,Species!="Chum")

adults <- subset(keogh_long,select = c(Year,Species,Stock))
adults <- reshape(adults,direction = "wide",idvar="Year",timevar="Species")
recruits <- subset(keogh_long,select = c(Year,Species,Recruits))
recruits <- reshape(recruits,direction = "wide",idvar="Year",timevar="Species")

juv_enviro <- subset(keogh_long,select = c(Year,Species,sumTemp,sumRain,winTemp,winRain,freshCoho,freshSteel,freshCutt,freshDolly,freshPink))
fresh_enviro <- reshape(juv_enviro,direction = "wide",idvar="Year",timevar="Species")
adult_enviro <- subset(keogh_long,select = c(Year,Species,seals,npgo,mei,oceanSalmon,juvCohort))
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

# multiple species: dolly varden, cutthroat trout, pink salmon, coho salmon
# including process and observation error
# precipitation covariates only affect observation model
# time-varying beta & alpha
# run a DLM on stock-recruitment for steelhead only

years <- recruits$Year
Nyears <- nrow(recruits)
Nspecies <- 5
colr <- wes_palette(wesAnderson,Nspecies,type=c("discrete"))

ln_RS_sh <- log(recruits$Recruits.Steelhead/adults$Stock.Steelhead)
ln_RS_dv <- log(recruits$`Recruits.Dolly Varden`/adults$`Stock.Dolly Varden`)
ln_RS_ct <- log(recruits$Recruits.Cutthroat/adults$Stock.Cutthroat)
ln_RS_pk <- log(recruits$Recruits.Pink/adults$Stock.Pink)
ln_RS_co <- log(recruits$Recruits.Coho/adults$Stock.Coho)
dat <- rbind(ln_RS_sh,ln_RS_dv,ln_RS_ct,ln_RS_pk,ln_RS_co)

keoghDLM$model
keoghDLM$states
mean(colSums(keoghDLM$states))
rowMeans(keoghDLM$states)
species <- unique(keogh_long$Species)
x <- species[4]
refYear <- 
alphas <- exp(keoghDLM$states[(1:10)%%2==1,15]) + sapply(species,function(x){sum(coef(keoghDLM)$D[grep(paste(x,"_",sep=""),row.names(coef(keoghDLM)$D))])})
betas <- keoghDLM$states[(1:10)%%2==0,15] + sapply(species,function(x){sum(coef(keoghDLM)$C[grep(x,row.names(coef(keoghDLM)$C))])})
names(alphas) <- names(betas) <- species

jpeg("Figures/recruitment.jpeg",width=7.5,height=5,units="in",res=800)
layout(matrix(1:6,nrow=2,ncol=3,byrow=TRUE))
par(mar=c(5,4,1,1))
for(i in 1:length(species))
{
  plot(keogh_long$Stock[keogh_long$Species==species[i]],keogh_long$Recruits[keogh_long$Species==species[i]],pch=21,bg=colr[i],xlab=paste(species[i],"adults",sep=" "),ylab=paste(species[i],"recruits",sep=" "),ylim=range(c(0,keogh_long$Recruits[keogh_long$Species==species[i]]),na.rm=TRUE),xlim=c(0,max(keogh_long$Stock[keogh_long$Species==species[i]],na.rm=TRUE)))
  for(t in 1:Nyears)
  {
    alphas <- exp(keoghDLM$states[(1:10)%%2==1,t][i]) + sum(coef(keoghDLM)$D[grep(paste(species[i],"_",sep=""),row.names(coef(keoghDLM)$D))])
    betas <- keoghDLM$states[(1:10)%%2==0,t][i] + sum(coef(keoghDLM)$C[grep(paste(species[i],"_",sep=""),row.names(coef(keoghDLM)$C))])
    curve(alphas*x*exp(betas*x),add=TRUE,xlab=paste(species[i],"adults",sep=" "),ylab=paste(species[i],"recruits",sep=" "),ylim=range(c(0,keogh_long$Recruits[keogh_long$Species==species[i]]),na.rm=TRUE),col=adjustcolor(colr[i],alpha=0.5))
  }
}
dev.off()

MARSSparamCIs(keoghDLM)
keoghDLM$coef

keoghAllfit <- augment(keoghDLM, interval="confidence")
keoghAllfit$Year <- keoghAllfit$t + 1975
keoghAllfit$Species <- keoghAllfit$.rownames
keoghAllfit$Species <- keogh_long$Species
margins <- c(0.5,0.5,0.5,1.1)
p <- ggplot(data = keoghAllfit) + 
  geom_line(aes(Year, .fitted)) +
  geom_ribbon(aes(x=Year, ymin=.conf.low, ymax=.conf.up), linetype=2, alpha=0.5) +
  geom_point(data=keoghAllfit, mapping = aes(x=Year, y=y,colour=Species)) +
  xlab("Year") + ylab("ln (recruits per spawner)") + facet_wrap(~Species,scales="free") +
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))
megaP <- ggarrange(p,ncol=1,nrow=1,legend="top",common.legend=TRUE)
pAnnotated <- annotate_figure(megaP,bottom=text_grob(wrapper("MARSS model fits to Keogh River stock-recruitment data",width=125),color="black",hjust=0,x=0.01,face="italic",size=10))

ggsave("Model fits.jpeg",plot=pAnnotated,units="in",height=5,width=7,dpi=800)

# temporal trends:
jpeg("Figures/temporal trends.jpeg",width=10,height=6,units="in",res=800)
m <- 2
ylabs <- rep(c(expression(log(alpha[t])),expression(beta[t]~"adults")),Nspecies)
titles <- c("Steelhead","Dolly Varden","Cutthroat","Pink","Coho")
layout(matrix(1:((m+1)*Nspecies),nrow=(m+1),ncol=Nspecies))
states <- (1:(Nspecies*m))[1:(Nspecies*m) %%2==1]
for(j in 1:Nspecies)
{
  state <- states[j]
  #temporal trends in alpha
  mn <- keoghDLM$states[state,]
  se <- keoghDLM$states.se[state,]
  par(mar=c(0,5,7,1))
  plot(years,mn,xlab="",ylab=ylabs[state],bty="n",xaxt="n",type="n",ylim=c(min(mn-2*se),max(mn+2*se)),cex.lab=1.5)
  title(titles[j])
  abline(v=1991,lwd=3)
  
  lines(years, rep(0,Nyears), lty="dashed")
  lines(years, mn, col=colr[j], lwd=3)
  lines(years, mn+2*se, col=colr[j])
  lines(years, mn-2*se, col=colr[j])
  # temporal trends in beta
  par(mar=c(3.5,5,3.5,1))
  mn <- keoghDLM$states[state+1,]
  se <- keoghDLM$states.se[state+1,]
  plot(years,mn,xlab="",ylab=ylabs[state+1],bty="n",xaxt="n",type="n",ylim=c(min(mn-2*se),max(mn+2*se)),cex.lab=1.5)

  lines(years, rep(0,Nyears), lty="dashed")
  lines(years, mn, col=colr[j], lwd=3)
  lines(years, mn+2*se, col=colr[j])
  lines(years, mn-2*se, col=colr[j])
  abline(v=1991,lwd=3)
  
  # temporal trends in carrying capacity
  par(mar=c(7,5,0,1))
  mn <- ifelse(-log(exp(keoghDLM$states[state,]))/keoghDLM$states[state+1,]<0,NA,-log(exp(keoghDLM$states[state,]))/keoghDLM$states[state+1,])
  plot(years,mn,xlab="",ylab=expression(~K[t]),bty="n",xaxt="n",type="n",cex.lab=2)
  lines(years, rep(0,Nyears), lty="dashed")
  lines(years, mn, col=colr[j], lwd=3)
  axis(1,at=seq(min(years),max(years),5),cex=2)
  mtext("Brood year", 1, line=3,cex=1)
  abline(v=1991,lwd=3)
}
dev.off()
matrix(coef(keoghDLM)$Q,nrow=10,byrow=FALSE)
coef(keoghDLM)$D[order(coef(keoghDLM)$D),1]
coef(keoghDLM)$C[order(coef(keoghDLM)$C),1]
keoghDLM$par$A
coef(keoghDLM)
coef(keoghDLM)$C
coef(keoghDLM,what="par.se")

CIs <- MARSSboot(keoghDLM,nboot = 100,param.gen="MLE")

d <- tidy(keoghDLM,type="states")
bootyMcBoot <- MARSSparamCIs(keoghDLM)

kf.out <- MARSSkfss(keoghDLM)

Phi <- kf.out$Vtt1
## obs variance; 1x1 matrix
R.est <- coef(keoghDLM, type="matrix")$R