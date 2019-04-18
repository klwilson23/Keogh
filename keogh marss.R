# RShowDoc("Chapter_DFA.R",package="MARSS")
# C:\Users\kylel\R\win-library\3.5\MARSS\doc

library(MARSS)
library(broom)
library(ggplot2)
SRdata <- read.csv("Keogh_StockRecruitment.csv",header=TRUE)
SRdata <- subset(SRdata,Year>=1976 & Year<=2013)
#SRdata <- SRdata[,-grep("ch_",colnames(SRdata))]
sdScale <- attr(scale(SRdata[,-1],center=FALSE,scale=TRUE),"scaled:scale")
SRdata[,-1] <- scale(SRdata[,-1],center=FALSE,scale=TRUE)

newDat <- SRdata[,-c(5,13)]
adultDat <- newDat[,c(2,4,5,7,9,11)]
adultDat <- t(as.matrix(adultDat))
colnames(adultDat) <- newDat$Year
ns <- nrow(adultDat)
B <- "diagonal and equal"
Q <- "unconstrained"
R <- diag(0.01,ns)
U <- "zero"
A <- "unequal"
x0 <- "zero"
mod.list = list(B=B, Q=Q, R=R, U=U, x0=x0, A=A,tinitx = 1)

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
mod.list.dfa = list(B = B, Z = Z, Q = Q, R = R, U = U, A = A, 
                    x0 = x0)

m <- apply(adultDat, 1, mean, na.rm=TRUE)
fit <- MARSS(adultDat, model=mod.list, control=list(minit=200,maxit=50000+200), inits=list(A=matrix(m,ns,1)))
fit$AICc
Z.est = coef(fit, type="matrix")$Z
H.inv = 1
if(ncol(Z.est)>1) H.inv = varimax(coef(fit, type="matrix")$Z)$rotmat
# rotate factor loadings
Z.rot = Z.est %*% H.inv
# rotate trends
trends.rot = solve(H.inv) %*% fit$states

layout(1)
matplot(t(trends.rot),type="l")
matplot(t(fit$states),type="l")

fit.b = getDFAfits(fit)

d <- augment(fit, interval = "confidence")
d$Year <- d$t + 1975
d$Species <- d$.rownames

fit$states

yy <- data.frame("Year"=newDat$Year,newDat[,c(2,4,5,7,9,11)])
yy <- reshape2::melt(yy, id.vars=c("Year"))
colnames(yy) <- c("Year","Species","Abundance")
p <- ggplot(data = d) + 
  geom_line(aes(Year, .fitted)) +
  geom_ribbon(aes(x=Year, ymin=.conf.low, ymax=.conf.up), linetype=2, alpha=0.5)
p <- p + geom_point(data=yy, mapping = aes(x=Year, y=Abundance))
p + facet_wrap(~Species) + xlab("") + ylab("Abundance")

layout(matrix(1:6,nrow=3,ncol=2))
apply(residuals(fit)$state.residuals[, 1:30], 1, acf)
mtext("State Residuals ACF", outer = TRUE, side = 3)
