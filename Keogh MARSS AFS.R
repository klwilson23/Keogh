source("some functions.R")
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggplot2)
keogh <- readRDS("Keogh_stockRec_enviro.rds")
head(keogh,2)

keogh_long <- acast(sampSize,year~life_stage,value.var="number")

example <- reshape(keogh,direction = "long",varying = list(c("sh_Adults","dv_Adults","ct_Adults","pk_Adults","ch_Adults","co_Adults"),c("sh_Smolts","dv_Smolts","ct_Smolts","pk_Recruits","ch_Recruits","co_Smolts")),v.names=c("Stock","Recruits"),idvar="Species")

example$Species <- rep(c("Steelhead","Dolly Varden","Cutthroat","Pink","Chum","Coho"),each=nrow(keogh))
example$Species <- factor(example$Species,levels=c("Steelhead","Dolly Varden","Cutthroat","Pink","Chum","Coho"))

keogh_long <- example
keogh_long <- subset(keogh_long,keogh_long$Year<=2015)
head(keogh_long,1)
colfunc <- colorRampPalette(c("royalblue4","dodgerblue","lightblue","darkorange1","firebrick"))

# http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
margins <- c(0.5,0.5,0.5,1.1)
p1 <- ggplot(data=keogh_long, aes(x=seal_density, y=Stock, color=Species))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_color_brewer(palette="Dark2") +
  labs(x="Straight of Georgia seal densities (per km)",y="Adults",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

p2 <- ggplot(data=keogh_long, aes(x=npgo, y=Stock, color=Species,palette="Dark2"))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_color_brewer(palette="Dark2") +
  labs(x="North Pacific Gyre Oscillation index",y="Adults",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

p3 <- ggplot(data=keogh_long, aes(x=total, y=Stock, color=Species))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_color_brewer(palette="Dark2") +
  labs(x="North Pacific salmon biomass (thousand metric tonnes)",y="Adults",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

p4 <- ggplot(data=keogh_long, aes(x=mean_temp, y=Stock, color=Species))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_color_brewer(palette="Dark2") +
  labs(x="Mean summer air temperatures",y="Adults",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

#megaPlot <- grid.arrange(p1,p2,p3,p4,nrow=4)
megaPlot <- ggarrange(p1,p2,p3,p4,ncol=1,nrow=4,legend="top",common.legend=TRUE)
megaPlot_Annotated <- annotate_figure(megaPlot,bottom=text_grob(wrapper("Relationships between adult salmon and marine and freshwater environment in the Keogh River",width=125),color="black",hjust=0,x=0.01,face="italic",size=10))

ggsave("adults and environment in Keogh.jpeg",plot=megaPlot_Annotated,units="in",height=12,width=10,dpi=800)

p1 <- ggplot(data=keogh_long, aes(x=temp_3, y=Recruits, color=Species))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_color_brewer(palette="Dark2") +
  labs(x="Mean summer air temperatures (3-year running average)",y="Recruits",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

p2 <- ggplot(data=keogh_long, aes(x=precip_3, y=Recruits, color=Species,palette="Dark2"))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_color_brewer(palette="Dark2") +
  labs(x="Summer rainfall (mm, 3-year running average)",y="Recruits",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

p3 <- ggplot(data=keogh_long, aes(x=win_temp_3, y=Recruits, color=Species))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_color_brewer(palette="Dark2") +
  labs(x="Mean winter air temperature (3-year running average)",y="Recruits",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

p4 <- ggplot(data=keogh_long, aes(x=win_precip_3, y=Recruits, color=Species))+ 
  geom_point()+
  geom_smooth()+
  facet_wrap(~Species,scales="free")+
  scale_color_brewer(palette="Dark2") +
  labs(x="Winter rainfall (3-year running average)",y="Recruits",color="Species") +
  theme_minimal() +
  theme(legend.position="none",strip.text.x = element_blank(),plot.margin=unit(margins,"line"))

#megaPlot <- grid.arrange(p1,p2,p3,p4,nrow=4)
megaPlot <- ggarrange(p1,p2,p3,p4,ncol=1,nrow=4,legend="top",common.legend=TRUE)
megaPlot_Annotated <- annotate_figure(megaPlot,bottom=text_grob(wrapper("Relationships between salmon recruits and freshwater environment in the Keogh River.Summer season defined as April through October.Winter season defined as November through February of adjacent years.",width=125),color="black",hjust=0,x=0.01,face="italic",size=10))

ggsave("recruits and environment in Keogh.jpeg",plot=megaPlot_Annotated,units="in",height=12,width=10,dpi=800)


layout(1)

columns <- c("sh_Smolts","sh_Adults","dv_Smolts","dv_Adults","ct_Smolts","ct_Adults","pk_Recruits","pk_Adults","co_Smolts","co_Adults","precip_3","win_precip_3","temp_3","npgo","mei","total","seal_density","bakun")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue",log="xy")

# steelhead
columns <- c("sh_Smolts","sh_Adults","total_rain","precip_1","precip_2","precip_3","win_rain","win_precip_1","win_precip_2","win_precip_3","mean_temp","temp_1","temp_2","temp_3","npgo","mei","pink","seal_density","bakun")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue",log="xy")

# dolly varden
columns <- c("dv_Smolts","dv_Adults","total_rain","precip_1","precip_2","precip_3","win_rain","win_precip_1","win_precip_2","win_precip_3","mean_temp","temp_1","temp_2","temp_3","npgo","mei","total","seal_density","bakun")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue",log="xy")

# coastal cutthroat
columns <- c("ct_Smolts","ct_Adults","total_rain","precip_1","precip_2","precip_3","win_rain","win_precip_1","win_precip_2","win_precip_3","mean_temp","temp_1","temp_2","temp_3","npgo","mei","pink","seal_density","bakun")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue",log="xy")

#coho
columns <- c("co_Smolts","co_Adults","total_rain","precip_1","precip_2","precip_3","win_rain","win_precip_1","win_precip_2","win_precip_3","mean_temp","temp_1","temp_2","temp_3","npgo","mei","pink","seal_density","bakun")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue",log="xy")

#pink
columns <- c("pk_Recruits","pk_Adults","total_rain","precip_1","precip_2","precip_3","win_rain","win_precip_1","win_precip_2","win_precip_3","mean_temp","temp_1","temp_2","temp_3","npgo","mei","pink","seal_density","bakun")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue")
pairs(keogh[,match(columns,colnames(keogh))],upper.panel=panel.smooth,lower.panel=panel.cor,pch=".",bg="dodgerblue",log="xy")
