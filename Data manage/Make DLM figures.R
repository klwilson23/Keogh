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
keogh_long$oceanSalmon <- residuals(lm(oceanSalmon~seals:Species,data=keogh_long))
keogh_long$seals <- log(keogh_long$seals)
Nspecies <- length(unique(keogh_long$Species))
adults <- subset(keogh_long,select = c(Year,Species,Stock))
adults <- reshape(adults,direction = "wide",idvar="Year",timevar="Species")
recruits <- subset(keogh_long,select = c(Year,Species,Recruits))
recruits <- reshape(recruits,direction = "wide",idvar="Year",timevar="Species")

years <- recruits$Year
Nyears <- nrow(recruits)
colr <- wes_palette(wesAnderson,Nspecies,type=c("discrete"))

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

# make plot of environment
margins <- c(0.5,0.5,0.5,1.1)
p1 <- ggplot(data = keogh_long) + 
  geom_line(aes(x=Year, y=seals,colour=Species)) +
  geom_smooth(data=keogh_long,aes(x=Year, y=seals))+
  xlab("Year") + ylab("Seal densities") +
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))
p2 <- ggplot(data = keogh_long) + 
  geom_line(aes(x=Year, y=oceanSalmon,colour=Species)) +
  geom_smooth(data=keogh_long,aes(x=Year, y=oceanSalmon))+
  xlab("Year") + ylab("North Pacific salmon (mt)") +
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))
p3 <- ggplot(data = keogh_long) + 
  geom_line(aes(x=Year, y=npgo,colour=Species)) +
  geom_smooth(data=keogh_long,aes(x=Year, y=npgo))+
  xlab("Year") + ylab("North Pacific gyre oscillation") +
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))
megaP <- ggarrange(p1,p2,p3,ncol=1,nrow=3,legend="top",common.legend=TRUE)
pAnnotated <- annotate_figure(megaP,bottom=text_grob(wrapper("Trends in marine and coastal conditions for Pacific salmonids in the Keogh River",width=125),color="black",hjust=0,x=0.01,face="italic",size=10))

ggsave("Figures/Seal and salmon rebuild.jpeg",plot=pAnnotated,units="in",height=7,width=6,dpi=800)

# make plot of environment
margins <- c(0.5,0.5,0.5,1.1)
p1 <- ggplot(data = keogh_long) + 
  geom_line(aes(x=Year, y=sumTemp,colour=Species)) +
  geom_smooth(data=keogh_long,aes(x=Year, y=sumTemp))+
  xlab("Year") + ylab("Summer air temperatures") +
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))
p2 <- ggplot(data = keogh_long) + 
  geom_line(aes(x=Year, y=winRain,colour=Species)) +
  geom_smooth(data=keogh_long,aes(x=Year, y=winRain))+
  xlab("Year") + ylab("Winter rainfall (mm)") +
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))
p3 <- ggplot(data = keogh_long) + 
  geom_line(aes(x=Year, y=freshCoho,colour=Species)) +
  geom_smooth(data=keogh_long,aes(x=Year, y=freshCoho))+
  xlab("Year") + ylab("Competing Coho salmon") +
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))
megaP <- ggarrange(p1,p2,p3,ncol=1,nrow=3,legend="top",common.legend=TRUE)
pAnnotated <- annotate_figure(megaP,bottom=text_grob(wrapper("Trends in freshwater conditions for Pacific salmonids in the Keogh River",width=125),color="black",hjust=0,x=0.01,face="italic",size=10))

ggsave("Figures/freshwater trends.jpeg",plot=pAnnotated,units="in",height=7,width=6,dpi=800)


# make Ricker plot
jpeg("Figures/Ricker model 1.jpeg",width=5,height=6,units="in",res=800)
layout(matrix(1:2,nrow=2))
par(mar=c(5,4,1,1))
alpha <- 4
beta <- -0.002
Smax <- 1/-beta
Rmax <- alpha*Smax*exp(beta*Smax)
K <- -log(alpha)/beta
curve(alpha*x*exp(beta*x),from=0,to=2*Smax,xlab="Spawner abundance",ylab="Recruits",lwd=2,col=colr[1],ylim=c(0,Rmax))
points(Smax,Rmax,pch=21,bg=colr[1])
curve(log(alpha)+beta*x,from=0,to=2*Smax,xlab="Spawner abundance",ylab="log(Recruits/Spawner)",lwd=2,col=colr[1],ylim=range(log(alpha)+beta*2*Smax,log(alpha)+beta*1e-3))
dev.off()

jpeg("Figures/Ricker model 2.jpeg",width=5,height=6,units="in",res=800)
layout(matrix(1:2,nrow=2))
par(mar=c(5,4,1,1))
alpha <- 4
beta <- -0.002
Smax <- 1/-beta
Rmax <- alpha*Smax*exp(beta*Smax)
K <- -log(alpha)/beta

alpha2 <- 4*0.5
beta2 <- -0.002
Smax2 <- 1/-beta2
Rmax2 <- alpha2*Smax2*exp(beta2*Smax2)
K2 <- -log(alpha2)/beta2

curve(alpha*x*exp(beta*x),from=0,to=2*Smax,xlab="Spawner abundance",ylab="Recruits",lwd=2,col=colr[1],ylim=c(0,Rmax))
curve(alpha2*x*exp(beta2*x),from=0,to=2*Smax,lwd=2,col=colr[2],ylim=c(0,Rmax2),add=TRUE)
points(Smax,Rmax,pch=21,bg=colr[1])
points(Smax2,Rmax2,pch=21,bg=colr[2])
curve(log(alpha)+beta*x,from=0,to=2*Smax,xlab="Spawner abundance",ylab="log(Recruits/Spawner)",lwd=2,col=colr[1],ylim=range(log(alpha)+beta*2*Smax,log(alpha)+beta*1e-3))
curve(log(alpha2)+beta2*x,from=0,to=2*Smax,lwd=2,col=colr[2],ylim=range(log(alpha2)+beta2*2*Smax2,log(alpha2)+beta2*1e-3),add=TRUE)

dev.off()

jpeg("Figures/Ricker model 3.jpeg",width=5,height=6,units="in",res=800)
layout(matrix(1:2,nrow=2))
par(mar=c(5,4,1,1))
alpha <- 4
beta <- -0.002
Smax <- 1/-beta
Rmax <- alpha*Smax*exp(beta*Smax)
K <- -log(alpha)/beta

alpha2 <- 4*0.5
beta2 <- -0.002
Smax2 <- 1/-beta2
Rmax2 <- alpha2*Smax2*exp(beta2*Smax2)
K2 <- -log(alpha2)/beta2

alpha3 <- 4
beta3 <- -0.002*2
Smax3 <- 1/-beta3
Rmax3 <- alpha3*Smax3*exp(beta3*Smax3)
K3 <- -log(alpha3)/beta3

curve(alpha*x*exp(beta*x),from=0,to=2*Smax,xlab="Spawner abundance",ylab="Recruits",lwd=2,col=colr[1],ylim=c(0,Rmax))
curve(alpha2*x*exp(beta2*x),from=0,to=2*Smax,lwd=2,col=colr[2],ylim=c(0,Rmax2),add=TRUE)
curve(alpha3*x*exp(beta3*x),from=0,to=2*Smax,lwd=2,col=colr[3],ylim=c(0,Rmax2),add=TRUE)
points(Smax,Rmax,pch=21,bg=colr[1])
points(Smax2,Rmax2,pch=21,bg=colr[2])
points(Smax3,Rmax3,pch=21,bg=colr[3])

curve(log(alpha)+beta*x,from=0,to=2*Smax,xlab="Spawner abundance",ylab="log(Recruits/Spawner)",lwd=2,col=colr[1],ylim=range(log(alpha)+beta*2*Smax,log(alpha)+beta*1e-3))
curve(log(alpha2)+beta2*x,from=0,to=2*Smax,lwd=2,col=colr[2],add=TRUE)
curve(log(alpha3)+beta3*x,from=0,to=2*Smax,lwd=2,col=colr[3],add=TRUE)
dev.off()

# multiple species: dolly varden, cutthroat trout, pink salmon, coho salmon
# including process and observation error
# precipitation covariates only affect observation model
# time-varying beta & alpha
# run a DLM on stock-recruitment for steelhead only

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
alphas <- exp(coef(keoghDLM)$D)
betas <- keoghDLM$states + sapply(species,function(x){sum(coef(keoghDLM)$C[grep(paste(x,"_",sep=""),row.names(coef(keoghDLM)$C))])})
names(alphas) <- names(betas) <- species

jpeg("Figures/recruitment.jpeg",width=7.5,height=5,units="in",res=800)
layout(matrix(1:6,nrow=2,ncol=3,byrow=TRUE))
par(mar=c(5,4,1,1))
for(i in 1:length(species))
{
  plot(keogh_long$Stock[keogh_long$Species==species[i]],keogh_long$Recruits[keogh_long$Species==species[i]],pch=21,bg=colr[i],xlab=paste(species[i],"adults",sep=" "),ylab=paste(species[i],"recruits",sep=" "),ylim=range(c(0,keogh_long$Recruits[keogh_long$Species==species[i]]),na.rm=TRUE),xlim=c(0,max(keogh_long$Stock[keogh_long$Species==species[i]],na.rm=TRUE)))
  for(t in 1:Nyears)
  {
    alpha <- alphas[i]
    betas <- keoghDLM$states[i,t]
    curve(alpha*x*exp(betas*x),add=TRUE,xlab=paste(species[i],"adults",sep=" "),ylab=paste(species[i],"recruits",sep=" "),ylim=range(c(0,keogh_long$Recruits[keogh_long$Species==species[i]]),na.rm=TRUE),col=adjustcolor(colr[i],alpha=0.5))
  }
}
dev.off()

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
jpeg("Figures/temporal trends.jpeg",width=10,height=4.5,units="in",res=800)
m <- 1
titles <- c("Steelhead","Dolly Varden","Cutthroat","Pink","Coho")
layout(matrix(1:((m+1)*Nspecies),nrow=(m+1),ncol=Nspecies))
states <- 1:5
for(j in 1:Nspecies)
{
  state <- states[j]
  # temporal trends in beta
  par(mar=c(2,5,5,1))
  mn <- keoghDLM$states[state,]
  se <- keoghDLM$states.se[state,]
  plot(years,mn,xlab="",ylab="Density-dependence",bty="n",xaxt="n",type="n",ylim=c(min(mn-2*se),max(mn+2*se)),cex.lab=1.1)

  lines(years, rep(0,Nyears), lty="dashed")
  lines(years, mn, col=colr[j], lwd=3)
  lines(years, mn+2*se, col=colr[j])
  lines(years, mn-2*se, col=colr[j])
  abline(v=1991,lwd=3)
  title(titles[j],font=2,cex=0.9,line=1)
  # temporal trends in carrying capacity
  Smax <- pmax(0,1/-(keoghDLM$states[state,]),na.rm=TRUE)
  Rmax <-  alphas[j,1]*Smax*exp(keoghDLM$states[state,]*Smax)
  Kt <- -log(alphas[j,1])/keoghDLM$states[state,]
  par(mar=c(5,5,2,1))
  mn <- ifelse(Rmax<0,NA,Rmax)
  plot(years,mn,xlab="",ylab="Freshwater capacity",bty="n",xaxt="n",type="n",cex.lab=1.1)
  lines(years, rep(0,Nyears), lty="dashed")
  lines(years, mn, col=colr[j], lwd=3)
  axis(1,at=seq(min(years),max(years),5),cex=2)
  mtext("Brood year", 1, line=3,cex=0.9)
  abline(v=1991,lwd=3)
}
dev.off()

# plot seals:
jpeg("Figures/effect of seals.jpeg",width=7,height=5.5,units="in",res=800)
species <- c("Steelhead","Dolly Varden","Cutthroat","Pink","Coho")
layout(matrix(1:6,nrow=3,ncol=2,byrow=TRUE))
par(mar=c(5,4,1,1))
for(i in 1:length(species)){
  plot(keogh_long$seals[keogh_long$Species==species[1]],keoghDLMspecies$states[1,],pch=21,bg=ifelse(years>1990,adjustcolor(colr[i],1),adjustcolor(colr[i],0.5)),xlab="Seal densities",ylab="Strength of density-dependence")
}
dev.off()

corr_mat <- coef(keoghDLMspecies,type="matrix")$Q

for(i in 1:nrow(coef(keoghDLMspecies,type="matrix")$Q)){
  for(j in 1:nrow(coef(keoghDLMspecies,type="matrix")$Q))
  {
    corr_mat[i,j] <- coef(keoghDLMspecies,type="matrix")$Q[i,j]/(sqrt(diag(coef(keoghDLMspecies,type="matrix")$Q)[i])*sqrt(diag(coef(keoghDLMspecies,type="matrix")$Q)[j]))
  }
}
