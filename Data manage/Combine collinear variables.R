source("some functions.R")
keogh <- readRDS("Data/Keogh_newRec_enviro.rds")
ocean_interact <- ocean_interact_2 <- rep(NA,nrow(keogh))
counter <- 1
for(i in unique(keogh$Species))
{
  sub_keogh <- keogh[keogh$Species==i,]
  nsub <- counter+nrow(sub_keogh)-1
  oceanCovar <- cbind(scale(sub_keogh[,c("seals","oceanSalmon")],center=F))
  pca <- princomp(oceanCovar)
  summary(pca)
  pca$loadings
  ocean_interact[counter:nsub] <- pca$scores[,1]
  ocean_interact_2[counter:nsub] <- pca$scores[,2]
  counter <- nsub+1
}
keogh$ocean_interact <- ocean_interact
keogh$ocean_covar_2 <- ocean_interact_2

saveRDS(keogh,"Data/Keogh_collinear_enviro.rds")

# plot(ocean_interact~oceanSalmon,data=keogh)
# plot(ocean_interact~seals,data=keogh)
# plot(ocean_interact_2~oceanSalmon,data=keogh)
# plot(ocean_interact_2~seals,data=keogh)
# plot(ocean_interact)
