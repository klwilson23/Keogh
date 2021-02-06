source("some functions.R")
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(MARSS)
library(broom)
library(wesanderson)
keogh <- readRDS("Data/Keogh_newJuv_enviro.rds")

keogh_long <- acast(sampSize,year~life_stage,value.var="number")
wesAnderson <- "Darjeeling1"

example <- reshape(keogh,direction = "long",varying = list(c("sh_Adults","dv_Adults","ct_Adults","pk_Adults","ch_Adults","co_Adults"),c("sh_Smolts","dv_Smolts","ct_Smolts","pk_Recruits","ch_Recruits","co_Smolts")),v.names=c("Stock","Recruits"),idvar="Species")

example$Species <- rep(c("Steelhead","Dolly Varden","Cutthroat","Pink","Chum","Coho"),each=nrow(keogh))
example$Species <- factor(example$Species,levels=c("Steelhead","Dolly Varden","Cutthroat","Pink","Chum","Coho"))

keogh_long <- example
keogh_long <- subset(keogh_long,keogh_long$Year<=2015)
head(keogh_long,1)
colfunc <- colorRampPalette(c("royalblue4","dodgerblue","lightblue","darkorange1","firebrick"))

keogh_long <- keogh
Nspecies <- length(unique(keogh_long$Species))
# http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
margins <- c(0.5,0.5,0.5,1.1)
p1 <- ggplot(data=keogh_long, aes(x=seals, y=Stock, color=Species))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  labs(x="Straight of Georgia seal densities (per km)",y="Adults",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

p2 <- ggplot(data=keogh_long, aes(x=npgo, y=Stock, color=Species,palette="Dark2"))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  labs(x="North Pacific Gyre Oscillation index",y="Adults",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

p3 <- ggplot(data=keogh_long, aes(x=oceanSalmon, y=Stock, color=Species))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  labs(x="North Pacific salmon biomass (thousand metric tonnes)",y="Adults",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

p4 <- ggplot(data=keogh_long, aes(x=juvCohort, y=Stock, color=Species))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  labs(x="Outmigrating cohort",y="Adults",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

#megaPlot <- grid.arrange(p1,p2,p3,p4,nrow=4)
megaPlot <- ggarrange(p1,p2,ncol=1,nrow=2,legend="top",common.legend=TRUE)
megaPlot_Annotated <- annotate_figure(megaPlot,bottom=text_grob(wrapper("Relationships between adult salmon and marine environment in the Keogh River",width=125),color="black",hjust=0,x=0.01,face="italic",size=10))

ggsave("adults and environment in Keogh.jpeg",plot=megaPlot_Annotated,units="in",height=6,width=8,dpi=800)

p1 <- ggplot(data=keogh_long, aes(x=sumTemp, y=Recruits, color=Species))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  labs(x="Summer air temperatures",y="Recruits",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

p2 <- ggplot(data=keogh_long, aes(x=sumRain, y=Recruits, color=Species,palette="Dark2"))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  labs(x="Summer rainfall (mm)",y="Recruits",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

p3 <- ggplot(data=keogh_long, aes(x=winTemp, y=Recruits, color=Species))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  labs(x="Winter air temperature",y="Recruits",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

p4 <- ggplot(data=keogh_long, aes(x=winRain, y=Recruits, color=Species))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  labs(x="Winter rainfall (mm)",y="Recruits",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

#megaPlot <- grid.arrange(p1,p2,p3,p4,nrow=4)
megaPlot <- ggarrange(p1,p2,p3,p4,ncol=1,nrow=4,legend="top",common.legend=TRUE)
megaPlot_Annotated <- annotate_figure(megaPlot,bottom=text_grob(wrapper("Relationships between salmon recruits and freshwater environment in the Keogh River.Summer season defined as April through October.Winter season defined as November through February of adjacent years.",width=125),color="black",hjust=0,x=0.01,face="italic",size=10))

ggsave("recruits and environment in Keogh.jpeg",plot=megaPlot_Annotated,units="in",height=12,width=10,dpi=800)


p1 <- ggplot(data=keogh_long, aes(x=sumTemp, y=Recruits, color=Species))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  labs(x="Summer air temperatures",y="Recruits",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))
pTemp <- ggarrange(p1,ncol=1,nrow=1,legend="top",common.legend=TRUE)

megaPlot_Annotated <- annotate_figure(pTemp,bottom=text_grob(wrapper("Relationships between salmon recruits and freshwater environment in the Keogh River.Summer season defined as April through October.",width=125),color="black",hjust=0,x=0.01,face="italic",size=10))
ggsave("recruits and temperature in Keogh.jpeg",plot=megaPlot_Annotated,units="in",height=5,width=7.5,dpi=800)



p1 <- ggplot(data=keogh_long, aes(x=freshPink, y=Recruits, color=Species))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  labs(x="Summer air temperatures",y="Recruits",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

p2 <- ggplot(data=keogh_long, aes(x=sumRain, y=Recruits, color=Species,palette="Dark2"))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  labs(x="Summer rainfall (mm)",y="Recruits",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

p3 <- ggplot(data=keogh_long, aes(x=winTemp, y=Recruits, color=Species))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  labs(x="Winter air temperature",y="Recruits",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

p4 <- ggplot(data=keogh_long, aes(x=winRain, y=Recruits, color=Species))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_colour_manual(values=wes_palette(n=Nspecies, name=wesAnderson)) +
  labs(x="Winter rainfall (mm)",y="Recruits",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

megaPlot <- ggarrange(p1,p2,p3,p4,ncol=1,nrow=4,legend="top",common.legend=TRUE)
megaPlot_Annotated <- annotate_figure(megaPlot,bottom=text_grob(wrapper("Relationships between salmon recruits and freshwater environment in the Keogh River.Summer season defined as April through October.Winter season defined as November through February of adjacent years.",width=125),color="black",hjust=0,x=0.01,face="italic",size=10))

ggsave("Figures/recruits and environment in Keogh.jpeg",plot=megaPlot_Annotated,units="in",height=12,width=10,dpi=800)