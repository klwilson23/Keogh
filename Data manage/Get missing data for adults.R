source("some functions.R")
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(MARSS)
library(broom)
keogh <- readRDS("Keogh_SR_enviro_new.rds")
keogh_long <- subset(keogh,Year<=2015 & Year>=1976)
keogh_long <- subset(keogh_long,Species!="Chum")

keogh_adults <- subset(keogh_long,select = c(Year,Species,Stock))

adults <- reshape(keogh_adults,direction = "wide",idvar="Year",timevar="Species")

adultNew <- adults
adultDat <- adults[,-1]
sdScale <- attr(scale(adultDat,center=FALSE,scale=TRUE),"scaled:scale")
adultDat <- scale(adultDat,center=FALSE,scale=TRUE)
adultDat <- t(as.matrix(adultDat))
colnames(adultDat) <- adults$Year

ns <- nrow(adultDat)
# AR1 model
B <- "diagonal and equal"
Q <- "unconstrained"
R <- diag(0.01, ns)
U <- "zero"
A <- "unequal"
x0 <- "unequal"
mod.list.ar1 = list(B = B, Q = Q, R = R, U = U, x0 = x0, A = A, tinitx = 1)

# DFA model
nTrends <- 1
B <- matrix(list(0), nTrends, nTrends)
diag(B) <- paste("b",1:nTrends,sep="")
Q <- diag(1, nTrends)
R <- "diagonal and unequal"
#R <- diag(0.01,ns)
U <- "zero"
x0 <- "zero"
Z <- matrix(list(0), ns, nTrends)
Z[1:(ns * nTrends)] <- sapply(1:nTrends,function(x){paste0("z",x,1:ns)})
Z[upper.tri(Z)] <- 0
A <- "unequal"
#A <- "zero"
mod.list.dfa = list(B = B, Z = Z, Q = Q, R = R, U = U, A = A, x0 = x0)

m <- apply(adultDat, 1, mean, na.rm=TRUE)
fit <- MARSS(adultDat, model=mod.list.dfa, control=list(minit=200,maxit=50000+200), inits=list(A=matrix(m,ns,1)))

d <- augment(fit, interval = "confidence")
d$Year <- d$t + 1975
spp <- matrix(unlist(strsplit(as.character(d$.rownames),".",fixed=TRUE)),ncol=2,byrow=TRUE)[,2]
d$Species <- factor(spp,levels=levels(keogh_adults$Species))

adultFits <- matrix(sapply(1:nrow(d),function(x){matches <- names(sdScale)%in%d$.rownames[x];
d$.fitted[x]*sdScale[matches]}),nrow=ncol(adultDat),ncol=nrow(adultDat),byrow=FALSE)

yy <- data.frame("Year"=adults$Year,t(adultDat))
yy <- reshape2::melt(yy, id.vars=c("Year"))
yy$Species <- d$Species
colnames(yy) <- c("Year","variable","Abundance","Species")
p <- ggplot(data = d) + 
  geom_line(aes(Year, .fitted)) + 
  geom_ribbon(aes(x = Year, ymin = .conf.low, ymax = .conf.up), linetype = 2, alpha = 0.5)
p <- p + geom_point(data = yy, mapping = aes(x = Year, y = Abundance))
p <- p + facet_wrap(~Species) + xlab("") + ylab("Abundance")
print(p)

adultEst <- data.frame("Year"=adultNew$Year,adultFits)
adultNew[which(is.na(adultNew),arr.ind=TRUE)] <- adultEst[which(is.na(adultNew),arr.ind=TRUE)]
adultNew[which(adultNew<=0,arr.ind=TRUE)] <- 1e-3
reAdult <- reshape(adultNew,direction="long",idvar="Year",timevar="Species")
colnames(reAdult) <- colnames(keogh_adults)

keogh_SR <- keogh_long[,!colnames(keogh_long)%in%"Stock"]
keogh_new <- data.frame(keogh_SR,"Stock"=reAdult$Stock)
plot(log(Recruits/Stock)~Stock,data=keogh_new[keogh_new$Species=="Coho",])
saveRDS(keogh_new,file="Keogh_newStock_enviro.rds")

