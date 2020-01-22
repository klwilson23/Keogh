library(ggplot2)
library(ggridges)
library(viridis)

# where adult run is a dataframe of individual fish where the columns are: (1) Year of sampling and (2) the day of the run (in this case day of spawning year, beginning November 15th)

adult_run$year <- factor(adult_run$year,levels=rev(unique(sort(adult_run$year))))
ggplot(adult_run, aes(x = run, y=as.factor(year),fill=..x..)) + 
  geom_density_ridges_gradient(scale=2.5,rel_min_height=1e-3) +
  scale_x_continuous(name="Run date (starting Nov. 15th)",expand=c(0.01,0)) +
  scale_y_discrete(expand=c(0.01,0)) +
  #scale_fill_viridis(name="Run time",option="C") +
  theme_ridges(font_size=13,grid=TRUE,center=TRUE) + 
  theme(axis.title.y=element_blank()) +
  scale_fill_gradient2(name="Run date",low="blue",high="orange",mid="dodgerblue",midpoint=75)

ggsave(filename="Figures/runtime v2.jpeg",units="in",height=7,width=6,dpi=800)

