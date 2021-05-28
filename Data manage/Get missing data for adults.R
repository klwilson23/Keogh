source("some functions.R")
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(MARSS)
library(broom)
keogh <- readRDS("Data/Keogh_SR_enviro_new.rds")
keogh_long <- subset(keogh,Year<=2015 & Year>=1976)
keogh_long <- subset(keogh_long,Species!="Chum")

keogh_adults <- subset(keogh_long,select = c(Year,Species,Stock))

coho <- read.csv("Data/Keogh coho adults.csv",stringsAsFactors = F,header=T)
coho <- coho[order(coho$Year),]
coho <- coho[coho$Year>=1976,]

adults <- reshape(keogh_adults,direction = "wide",idvar="Year",timevar="Species")
adults$Stock.Coho2 <- rep(NA,nrow(adults))
adults$Stock.Coho2[adults$Year < 1998] <- coho$Adults[coho$Year < 1998]
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
Z[nrow(Z),] <- Z[nrow(Z)-1,]
A <- "unequal"
#A <- "zero"
mod.list.dfa = list(B = B, Z = Z, Q = Q, R = R, U = U, A = A, x0 = x0)

m <- apply(adultDat, 1, mean, na.rm=TRUE)
fit <- MARSS(adultDat, model=mod.list.dfa, control=list(minit=200,maxit=5000+200), inits=list(A=matrix(m,ns,1)))

d <- fitted(fit,interval="confidence")
d$Year <- d$t + 1975
spp <- matrix(unlist(strsplit(as.character(d$.rownames),".",fixed=TRUE)),ncol=2,byrow=TRUE)[,2]
d$Species <- factor(spp,levels=c(levels(keogh_adults$Species),"Coho2"))

adultFits <- matrix(sapply(1:nrow(d),function(x){matches <- names(sdScale)%in%d$.rownames[x];
d$.fitted[x]*sdScale[matches]}),nrow=ncol(adultDat),ncol=nrow(adultDat),byrow=FALSE)

yy <- data.frame("Year"=adults$Year,t(adultDat))
yy <- reshape2::melt(yy, id.vars=c("Year"))
d$Species <- factor(d$Species,levels=levels(d$Species),labels=c("Steelhead","Dolly Varden","Cutthroat","Pink","Chum","Coho - Resistivity","Coho - NuSEDS"))
yy$Species <- factor(d$Species,levels=levels(d$Species),labels=c("Steelhead","Dolly Varden","Cutthroat","Pink","Chum","Coho - Resistivity","Coho - NuSEDS"))
colnames(yy) <- c("Year","variable","Abundance","Species")
p <- ggplot(data = d) + 
  geom_line(aes(Year, .fitted)) + 
  geom_ribbon(aes(x = Year, ymin = .conf.low, ymax = .conf.up), linetype = 2, alpha = 0.5)
p <- p + geom_point(data = yy, mapping = aes(x = Year, y = Abundance))
p <- p + facet_wrap(~Species) + xlab("") + ylab("Standardized abundance")
print(p)

ggsave("Figures/Adult missing data.jpeg",plot=p,dpi=800,units="in",height=5,width=7)

adultEst <- data.frame("Year"=adultNew$Year,adultFits)
adultNew[which(is.na(adultNew),arr.ind=TRUE)] <- adultEst[which(is.na(adultNew),arr.ind=TRUE)]
adultNew[which(adultNew<=0,arr.ind=TRUE)] <- 1e-3
#adultNew <- adultNew[,-(match("Stock.Coho2",colnames(adultNew)))]
reAdult <- reshape(adultNew,direction="long",idvar="Year",timevar="Species")
reAdult <- reAdult[,-match("Stock.Coho2",colnames(reAdult))]
colnames(reAdult) <- colnames(keogh_adults)
#reAdult <- reAdult[reAdult$Species,]
keogh_SR <- keogh_long[,!colnames(keogh_long)%in%"Stock"]
keogh_new <- data.frame(keogh_SR,"Stock"=reAdult$Stock)
plot(log(Recruits/Stock)~Stock,data=keogh_new[keogh_new$Species=="Coho",])

coho <- read.csv("Data/Keogh coho adults.csv",stringsAsFactors = F,header=T)
coho <- coho[order(coho$Year),]

plot(coho$Year,coho$Adults,pch=21,bg="orange",ylab="Coho adults",xlab="Year")
points(adults$Year,adults$Stock.Coho,pch=21,bg="dodgerblue")
points(adultEst$Year[is.na(adults$Stock.Coho)],adultEst$X5[is.na(adults$Stock.Coho)],type="l",lwd=2,col="black")
legend("topright",c("NuSeds Escapement","Keogh resistivity counter","Imputed from DFA"),pch=c(21,21,NA),lty=c(NA,NA,1),lwd=c(NA,NA,2),pt.bg=c("orange","dodgerblue",NA),col=c(1,1,1),bty="n")
saveRDS(keogh_new,file="Data/Keogh_newStock_enviro.rds")

