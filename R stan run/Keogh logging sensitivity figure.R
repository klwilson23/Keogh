library(ggplot2)
lag_seq <- c(5,10,15,30)
beta_logging <- readRDS(file="Results/forestry lag sensitivity.rds")
new_logging <- lapply(beta_logging,as.data.frame)
new_df <- Reduce(function(...) merge(... , all=T,sort=F),new_logging)
new_df$species <- factor(rep(c("Dolly Varden","Steelhead","Coastal Cutthroat","Coho","Pink"),length(lag_seq)),levels=c(c("Dolly Varden","Steelhead","Coastal Cutthroat","Coho","Pink")))
new_df$lag_years <- rep(lag_seq,each=length(unique(new_df$species)))
new_df$`2.5%`
Nspecies <- length(unique(new_df$species))
ggplot(new_df,aes(x=lag_years,y=mean,colour=species)) +
  geom_line() +
  geom_ribbon(aes(x=lag_years,ymin=`2.5%`,ymax=`97.5%`,fill=species),alpha=0.1) +
  geom_point(aes(fill=species),pch=21,colour='black') +
  facet_wrap(~species,scales='fixed') +
  geom_hline(yintercept=0)+
  xlab("Lagged years since logging") + ylab("Effect sizes on recruitment productivity") +
  scale_colour_brewer("Species",type="qual",palette = 2) +
  scale_fill_brewer("Species",type="qual",palette = 2) +
  theme_minimal() +
  theme(legend.position="top",strip.text.x=element_text(hjust=0),axis.text.x=element_text(size=7,angle=45,hjust=1))

ggsave(filename="Figures/Lagged forestry sensitivity.jpeg",width=7,height=6,dpi=800)
        