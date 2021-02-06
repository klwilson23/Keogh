source("some functions.R")
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(MARSS)
library(broom)
keogh <- readRDS("Data/Keogh_newJuv_enviro.rds")
keogh_long <- subset(keogh,Year<=2015 & Year>=1976)
keogh_long <- subset(keogh_long,Species!="Chum")

keogh_rec <- subset(keogh_long,select = c(Year,Species,Recruits))
recruits <- reshape(keogh_rec,direction = "wide",idvar="Year",timevar="Species")

recNew <- recruits
recDat <- log(recruits[,-1])
sdScale <- attr(scale(recDat,center=FALSE,scale=TRUE),"scaled:scale")
recDat <- scale(recDat,center=FALSE,scale=TRUE)
recDat <- t(as.matrix(recDat))
colnames(recDat) <- recruits$Year

ns <- nrow(recDat)
# AR1 model
B <- "diagonal and equal"
Q <- "unconstrained"
R <- diag(0.01, ns)
U <- "zero"
A <- "unequal"
x0 <- "unequal"
mod.list.ar1 = list(B = B, Q = Q, R = R, U = U, x0 = x0, A = A, tinitx = 1)

# DFA model
nTrends <- 2
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

m <- apply(recDat, 1, mean, na.rm=TRUE)
fit <- MARSS(recDat, model=mod.list.ar1, control=list(minit=200,maxit=50000+200), inits=list(A=matrix(m,ns,1)))

d <- augment(fit, interval = "confidence")
d$Year <- d$t + 1975
spp <- matrix(unlist(strsplit(as.character(d$.rownames),".",fixed=TRUE)),ncol=2,byrow=TRUE)[,2]
d$Species <- factor(spp,levels=levels(keogh_rec$Species))

recFits <- matrix(sapply(1:nrow(d),function(x){matches <- names(sdScale)%in%d$.rownames[x];
d$.fitted[x]*sdScale[matches]}),nrow=ncol(recDat),ncol=nrow(recDat),byrow=FALSE)

yy <- data.frame("Year"=recruits$Year,t(recDat))
yy <- reshape2::melt(yy, id.vars=c("Year"))
yy$Species <- d$Species
colnames(yy) <- c("Year","variable","Abundance","Species")
p <- ggplot(data = d) + 
  geom_line(aes(Year, .fitted)) + 
  geom_ribbon(aes(x = Year, ymin = .conf.low, ymax = .conf.up), linetype = 2, alpha = 0.5)
p <- p + geom_point(data = yy, mapping = aes(x = Year, y = Abundance))
p <- p + facet_wrap(~Species) + xlab("") + ylab("Abundance")
print(p)

recEst <- data.frame("Year"=recNew$Year,recFits)
recNew[which(is.na(recNew),arr.ind=TRUE)] <- exp(recEst[which(is.na(recNew),arr.ind=TRUE)])
recNew[which(recNew<=0,arr.ind=TRUE)] <- log(1e-3)
reRec <- reshape(recNew,direction="long",idvar="Year",timevar="Species")
colnames(reRec) <- colnames(keogh_rec)

keogh_SR <- keogh_long[,!colnames(keogh_long)%in%"Recruits"]
keogh_new <- data.frame(keogh_SR,"Recruits"=reRec$Recruits)
plot(Recruits~Year,data=keogh_new[keogh_new$Species=="Dolly Varden",])
plot(log(Recruits/Stock)~Year,data=keogh_new[keogh_new$Species=="Dolly Varden",])

saveRDS(keogh_new,file="Data/Keogh_newRec_enviro.rds")

