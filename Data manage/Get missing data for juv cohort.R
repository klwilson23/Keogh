source("some functions.R")
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(MARSS)
library(broom)
keogh <- readRDS("Keogh_newStock_enviro.rds")
keogh_long <- subset(keogh,Year<=2015 & Year>=1976)
keogh_long <- subset(keogh_long,Species!="Chum")

keogh_juv <- subset(keogh_long,select = c(Year,Species,juvCohort))
juvCohort <- reshape(keogh_juv,direction = "wide",idvar="Year",timevar="Species")

juvNew <- juvCohort
juvDat <- juvCohort[,-1]
sdScale <- attr(scale(juvDat,center=FALSE,scale=TRUE),"scaled:scale")
juvDat <- scale(juvDat,center=FALSE,scale=TRUE)
juvDat <- t(as.matrix(juvDat))
colnames(juvDat) <- juvCohort$Year

ns <- nrow(juvDat)
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

m <- apply(juvDat, 1, mean, na.rm=TRUE)
fit <- MARSS(juvDat, model=mod.list.ar1, control=list(minit=200,maxit=50000+200), inits=list(A=matrix(m,ns,1)))

d <- augment(fit, interval = "confidence")
d$Year <- d$t + 1975
spp <- matrix(unlist(strsplit(as.character(d$.rownames),".",fixed=TRUE)),ncol=2,byrow=TRUE)[,2]
d$Species <- factor(spp,levels=levels(keogh_juv$Species))

juvFits <- matrix(sapply(1:nrow(d),function(x){matches <- names(sdScale)%in%d$.rownames[x];
d$.fitted[x]*sdScale[matches]}),nrow=ncol(juvDat),ncol=nrow(juvDat),byrow=FALSE)

yy <- data.frame("Year"=juvCohort$Year,t(juvDat))
yy <- reshape2::melt(yy, id.vars=c("Year"))
yy$Species <- d$Species
colnames(yy) <- c("Year","variable","Abundance","Species")
p <- ggplot(data = d) + 
  geom_line(aes(Year, .fitted)) + 
  geom_ribbon(aes(x = Year, ymin = .conf.low, ymax = .conf.up), linetype = 2, alpha = 0.5)
p <- p + geom_point(data = yy, mapping = aes(x = Year, y = Abundance))
p <- p + facet_wrap(~Species) + xlab("") + ylab("Abundance")
print(p)

juvEst <- data.frame("Year"=juvNew$Year,juvFits)
juvNew[which(is.na(juvNew),arr.ind=TRUE)] <- juvEst[which(is.na(juvNew),arr.ind=TRUE)]
juvNew[which(juvNew<=0,arr.ind=TRUE)] <- 1e-3
reJuv <- reshape(juvNew,direction="long",idvar="Year",timevar="Species")
colnames(reJuv) <- colnames(keogh_juv)

keogh_SR <- keogh_long[,!colnames(keogh_long)%in%"juvCohort"]
keogh_new <- data.frame(keogh_SR,"juvCohort"=reJuv$juvCohort)
plot(Stock/juvCohort~Year,data=keogh_new[keogh_new$Species=="Coho",])
saveRDS(keogh_new,file="Keogh_newJuv_enviro.rds")

