library(MARSS)
library(broom)
SRdata <- read.csv("Keogh_StockRecruitment.csv",header=TRUE)
SRdata <- subset(SRdata,Year>=1976 & Year<=2013)
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

B <- matrix(list(0), 2, 2)
B[1, 1] <- "b1"
B[2, 2] <- "b2"
Q <- diag(1, 2)
R <- "diagonal and unequal"
U <- "zero"
x0 <- "zero"
Z <- matrix(list(0), ns, 2)
Z[1:(ns * 2)] <- c(paste0("z1", 1:ns), paste0("z2", 1:ns))
Z[1, 2] <- 0
A <- "unequal"
mod.list.dfa = list(B = B, Z = Z, Q = Q, R = R, U = U, A = A, 
                    x0 = x0)

m <- apply(adultDat, 1, mean, na.rm=TRUE)
fit <- MARSS(adultDat, model=mod.list, control=list(maxit=5000), inits=list(A=matrix(m,ns,1)))

d <- augment(fit, interval = "confidence")
d$Year <- d$t + 1975
d$Species <- d$.rownames

fit$states[1,]

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
