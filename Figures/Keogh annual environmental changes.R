source("some functions.R")
library(tidyverse)
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(MARSS)
library(broom)
library(wesanderson)
library(grid)
environ <- readRDS("environ_covars.rds")
environment <- read.csv("Data/keogh environmental covariates.csv",header=TRUE)
environ$seal_abundance[environ$year>=2014 & is.na(environ$seal_abundance)] <- environ$seal_abundance[environ$year==2014]
environ$seal_density[environ$year>=2014 & is.na(environ$seal_density)] <- environ$seal_density[environ$year==2014]
environ$total[environ$year>=2015 & is.na(environ$total)] <- environ$total[environ$year==2015]
forestry <- readRDS("Data/keogh_logging.rds")
forestry <- subset(forestry,forestry$Year>=min(environ$year))

environ$Year <- environ$year
environ <- environ[,-match("year",colnames(environ))]
environmental <- merge(environ,forestry,by="Year")
environmental <- environmental[environmental$Year>=1970,]
sc_enviro <- data.frame("Year"=environmental$Year,scale(environmental[,-match("Year",colnames(environmental))]))
layout(matrix(1:2,nrow=2,ncol=1))
par(mar=c(5,4,0.5,0.5),mai=c(0,1,1,0.1))
plot(npgo~Year,data=sc_enviro,col="blue",lty=1,type="l",ylim=c(-3,3),xaxt="n",ylab="Marine environment")
axis(1,labels=FALSE,tick=TRUE)
lines(sc_enviro$Year,sc_enviro$seal_density,lwd=1,col="blue",lty=2)
lines(sc_enviro$Year,sc_enviro$total,lwd=1,col="blue",lty=3)
par(mar=c(5,4,0.5,0.5),mai=c(1,1,0,0.1))
plot(win_mean_temp~Year,data=sc_enviro,col="forestgreen",lty=1,type="l",ylim=c(-3,3),ylab="Freshwater environment")
lines(sc_enviro$Year,sc_enviro$cumul_log,lwd=1,col="forestgreen",lty=2)
lines(sc_enviro$Year,sc_enviro$win_rain,lwd=1,col="forestgreen",lty=3)

df_mar <- sc_enviro %>%
  select(Year, seal_density, npgo, total,win_mean_temp, win_rain, cumul_log) %>%
  gather(key = "variable", value = "value", -Year)
df_mar$variable <- factor(df_mar$variable,levels=c("seal_density","npgo","total","win_mean_temp", "win_rain", "cumul_log"),labels=c("Seal densities","North Pacific Gyre Oscillations","North Pacific salmon biomass","Winter air temperatures","Winter rainfall","Cumulative logging (15 year window)"))
mar_labs <- c("Seal densities","North Pacific Gyre Oscillations","North Pacific salmon biomass","Winter air temperatures","Winter rainfall","Cumulative logging (15 year window)")
df_mar$ecosystem <- factor(ifelse(df_mar$variable %in% c("Seal densities","North Pacific Gyre Oscillations","North Pacific salmon biomass"),"Marine","Freshwater"),levels=c("Marine","Freshwater"),labels=c("(a) - Marine","(b) - Freshwater"))


margins <- c(0.5,0.5,0.5,1.1)
p1 <- ggplot(df_mar, aes(x = Year, y = value)) + 
  geom_line(aes(color=variable),size=1.2) +
  labs(x="Year",y="Covariate value",color="Covariate") +
  #scale_color_manual(values = c("dodgerblue", "steelblue","darkblue","darkgreen", "forestgreen","olivedrab4")) +
  scale_color_brewer(type="div",direction = -1) +
  facet_wrap(~ ecosystem,nrow=2,scales="free") +
  theme_minimal() +
  theme(legend.position="none",plot.margin=unit(margins,"line"),strip.text=element_text(hjust=0))
p1
megaPlot <- ggarrange(p1,ncol=1,nrow=1,legend="top",common.legend=TRUE,labels="")
megaPlot_Annotated <- annotate_figure(megaPlot,bottom=text_grob(wrapper("Changing marine and freshwater environments in the Keogh River",width=125),color="black",hjust=0,x=0.01,face="italic",size=10))
megaPlot_Annotated
ggsave("Changing marine and freshwater environments.jpeg",plot=megaPlot_Annotated,units="in",height=6,width=8,dpi=800)


# make plot of environment
wesAnderson <- "Darjeeling1"
keogh <- readRDS("Keogh_newRec_enviro.rds")
keogh_long <- subset(keogh,Year<=2015 & Year>=1976)
keogh_long <- subset(keogh_long,Species!="Chum")
keogh_long$Species <- factor(keogh_long$Species,levels=unique(keogh_long$Species))
keogh_long$marSurv <- keogh_long$Stock/keogh_long$juvCohort
keogh_long$logitSurv <- log(keogh_long$marSurv/(1-keogh_long$marSurv))
keogh_long$prod <- log(keogh_long$Recruits/keogh_long$Stock)
Nspecies <- length(unique(keogh_long$Species))

margins <- c(0.5,0.5,0.5,1.1)
p1 <- ggplot(data = keogh_long) + 
  geom_line(aes(x=Year, y=Recruits,colour=Species)) +
  geom_point(aes(x=Year, y=Recruits,colour=Species)) +
  #geom_smooth(data=keogh_long,aes(x=Year, y=Recruits))+
  xlab("Year") + ylab("Recruits") +
  facet_wrap(~Species,scales="free") +
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))
p2 <- ggplot(data = keogh_long) + 
  geom_line(aes(x=Year, y=Stock,colour=Species)) +
  geom_point(aes(x=Year, y=Stock,colour=Species)) +
  #geom_smooth(data=keogh_long,aes(x=Year, y=Stock))+
  xlab("Year") + ylab("Adults") +
  facet_wrap(~Species,scales="free") +
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))
megaP <- ggarrange(p1,p2,ncol=1,nrow=2,legend="top",common.legend=TRUE)
pAnnotated <- annotate_figure(megaP,bottom=text_grob(wrapper("Recruits and spawners since 1976 on the Keogh River",width=125),color="black",hjust=0,x=0.01,face="italic",size=10))

ggsave("Figures/stock and recruitment trends.jpeg",plot=pAnnotated,units="in",height=6,width=9,dpi=800)