source("some functions.R")
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(MARSS)
library(broom)
library(wesanderson)
keogh <- readRDS("Keogh_newJuv_enviro.rds")
keogh_long <- subset(keogh,Year<=2015 & Year>=1976)
keogh_long <- subset(keogh_long,Species!="Chum")

adults <- subset(keogh_long,select = c(Year,Species,Stock))
adults <- reshape(adults,direction = "wide",idvar="Year",timevar="Species")
recruits <- subset(keogh_long,select = c(Year,Species,Recruits))
recruits <- reshape(recruits,direction = "wide",idvar="Year",timevar="Species")

juv_enviro <- subset(keogh_long,select = c(Year,Species,sumTemp,sumRain,winTemp,winRain,freshCoho,freshSteel,freshCutt,freshDolly,freshPink))
fresh_enviro <- reshape(juv_enviro,direction = "wide",idvar="Year",timevar="Species")
adult_enviro <- subset(keogh_long,select = c(Year,Species,seals,npgo,mei,oceanSalmon,juvCohort))
ocean_enviro <- reshape(adult_enviro,direction = "wide",idvar="Year",timevar="Species")

freshEnviroNew <- fresh_enviro
sdCovarsFresh <- attr(scale(fresh_enviro[,-1],center=TRUE,scale=TRUE),"scaled:scale")
mnCovarsFresh <- attr(scale(fresh_enviro[,-1],center=TRUE,scale=TRUE),"scaled:center")
freshCovarScale <- scale(fresh_enviro[,-1],center=TRUE,scale=TRUE)
#covarScale[is.na(covarScale)] <- 0
oceanEnviroNew <- ocean_enviro
sdCovarsOcean <- attr(scale(ocean_enviro[,-1],center=TRUE,scale=TRUE),"scaled:scale")
mnCovarsOcean <- attr(scale(ocean_enviro[,-1],center=TRUE,scale=TRUE),"scaled:center")
oceanCovarScale <- scale(ocean_enviro[,-1],center=TRUE,scale=TRUE)

# multiple species: dolly varden, cutthroat trout, pink salmon, coho salmon
# including process and observation error
# precipitation covariates only affect observation model
# time-varying beta & alpha
# run a DLM on stock-recruitment for steelhead only

C <- c <- A <- U <- "zero"
years <- recruits$Year
Nyears <- nrow(recruits)
Nspecies <- 5
m <- 2 # number of time-varying parameters
Z <- array(0,c(Nspecies,m*Nspecies,Nyears))
Z[1,1,] <- rep(1,Nyears) # time-varying alpha
Z[1,2,] <- adults$Stock.Steelhead # time-varying beta
Z[2,3,] <- rep(1,Nyears) # time-varying alpha
Z[2,4,] <- adults$`Stock.Dolly Varden` # time-varying beta
Z[3,5,] <- rep(1,Nyears) # time-varying alpha
Z[3,6,] <- adults$Stock.Cutthroat # time-varying beta
Z[4,7,] <- rep(1,Nyears) # time-varying alpha
Z[4,8,] <- adults$Stock.Pink # time-varying beta
Z[5,9,] <- rep(1,Nyears) # time-varying alpha
Z[5,10,] <- adults$Stock.Coho # time-varying beta

# species competition, summer climate, and winter climate affects process model
d <- rbind(t(freshCovarScale))
D <- matrix(list(0),nrow=nrow(d),ncol=Nspecies)
coefNames <- paste("b",matrix(unlist(strsplit(row.names(d),"\\.")),ncol=2,byrow=TRUE)[,2],matrix(unlist(strsplit(row.names(d),"\\.")),ncol=2,byrow=TRUE)[,1],sep="_")
D[grep("Steelhead",coefNames),1] <- coefNames[grep("Steelhead",coefNames)]
D[grep("Dolly Varden",coefNames),2] <- coefNames[grep("Dolly Varden",coefNames)]
D[grep("Cutthroat",coefNames),3] <- coefNames[grep("Cutthroat",coefNames)]
D[grep("Pink_",coefNames),4] <- coefNames[grep("Pink_",coefNames)]
D[grep("Coho_",coefNames),5] <- coefNames[grep("Coho_",coefNames)]

# NPGO, seals, pacific salmon, and mei affects observation model
c <- rbind(t(oceanCovarScale))
C <- matrix(list(0),nrow=m*Nspecies,ncol=nrow(c))
alphaCoefNames <- paste("a",matrix(unlist(strsplit(row.names(c),"\\.")),ncol=2,byrow=TRUE)[,2],matrix(unlist(strsplit(row.names(c),"\\.")),ncol=2,byrow=TRUE)[,1],sep="_")
betaCoefNames <- paste("b",matrix(unlist(strsplit(row.names(c),"\\.")),ncol=2,byrow=TRUE)[,2],matrix(unlist(strsplit(row.names(c),"\\.")),ncol=2,byrow=TRUE)[,1],sep="_")
C[1:2,grep("Steelhead",adultCoefNames)] <- rbind(alphaCoefNames[grep("Steelhead",adultCoefNames)],
                                                 betaCoefNames[grep("Steelhead",adultCoefNames)])
C[3:4,grep("Dolly Varden",adultCoefNames)] <- rbind(alphaCoefNames[grep("Dolly Varden",adultCoefNames)],
                                                  betaCoefNames[grep("Dolly Varden",adultCoefNames)])
C[5:6,grep("Cutthroat",adultCoefNames)] <- rbind(alphaCoefNames[grep("Cutthroat",adultCoefNames)],
                                               betaCoefNames[grep("Cutthroat",adultCoefNames)])
C[7:8,grep("Pink_",adultCoefNames)] <- rbind(alphaCoefNames[grep("Pink_",adultCoefNames)],
                                           betaCoefNames[grep("Pink_",adultCoefNames)])
C[9:10,grep("Coho_",adultCoefNames)] <- rbind(alphaCoefNames[grep("Coho_",adultCoefNames)],
                                           betaCoefNames[grep("Coho_",adultCoefNames)])

c <- rbind(t(oceanCovarScale))
C2 <- matrix(list(0),nrow=m*Nspecies,ncol=nrow(c))
C2[2,grep("Steelhead",adultCoefNames)] <- betaCoefNames[grep("Steelhead",adultCoefNames)]
C2[4,grep("Dolly Varden",adultCoefNames)] <- betaCoefNames[grep("Dolly Varden",adultCoefNames)]
C2[6,grep("Cutthroat",adultCoefNames)] <- betaCoefNames[grep("Cutthroat",adultCoefNames)]
C2[8,grep("Pink_",adultCoefNames)] <- betaCoefNames[grep("Pink_",adultCoefNames)]
C2[10,grep("Coho_",adultCoefNames)] <- betaCoefNames[grep("Coho_",adultCoefNames)]

Q <- matrix(list(0),m*Nspecies,m*Nspecies)        ## 2x2; all 0 for now
#diag(Q) <- rep(c("q.alpha","q.beta"),Nspecies)
diag_Q <- c("sh_alpha","sh_adults","dv_alpha","dv_adults","ct_alpha","ct_adults","pk_alpha","pk_adults","co_alpha","co_adults") ## 2x2; diag = (q1,q2)
alphaCovariances <- paste("corr_alpha",0:((Nspecies-1)*(Nspecies-1)-1),sep="")
betaCovariances <- paste("corr_beta",0:((Nspecies-1)*(Nspecies-1)-1),sep="")

index = 1
Q[1,1:ncol(Q)%%2==1] <- alphaCovariances[1:sum(1:ncol(Q)%%2==1)]
Q[1:nrow(Q)%%2==1,1] <- alphaCovariances[1:sum(1:nrow(Q)%%2==1)]
index = index+sum(1:ncol(Q)%%2==1)
Nrow <- 3
coord <- (Nrow:ncol(Q))[(Nrow:ncol(Q))%%2==1]
Q[Nrow,coord] <- alphaCovariances[index:((index-1)+length(coord))]
Q[coord,Nrow] <- alphaCovariances[index:((index-1)+length(coord))]

index = index+sum(Nrow:ncol(Q)%%2==1)
Nrow <- 5
coord <- (Nrow:ncol(Q))[(Nrow:ncol(Q))%%2==1]
Q[Nrow,coord] <- alphaCovariances[index:((index-1)+length(coord))]
Q[coord,Nrow] <- alphaCovariances[index:((index-1)+length(coord))]

index = index+sum(Nrow:ncol(Q)%%2==1)
Nrow <- 7
coord <- (Nrow:ncol(Q))[(Nrow:ncol(Q))%%2==1]
Q[Nrow,coord] <- alphaCovariances[index:((index-1)+length(coord))]
Q[coord,Nrow] <- alphaCovariances[index:((index-1)+length(coord))]

index = index+sum(Nrow:ncol(Q)%%2==1)
Nrow <- 9
coord <- (Nrow:ncol(Q))[(Nrow:ncol(Q))%%2==1]
Q[Nrow,coord] <- alphaCovariances[index:((index-1)+length(coord))]
Q[coord,Nrow] <- alphaCovariances[index:((index-1)+length(coord))]

# covariances between beta
index = 1
Q[2,1:ncol(Q)%%2==0] <- betaCovariances[1:sum(1:ncol(Q)%%2==0)]
Q[1:nrow(Q)%%2==0,2] <- betaCovariances[1:sum(1:nrow(Q)%%2==0)]
index = index+sum(1:ncol(Q)%%2==1)
Nrow <- 4
coord <- (Nrow:ncol(Q))[(Nrow:ncol(Q))%%2==0]
Q[Nrow,coord] <- betaCovariances[index:((index-1)+length(coord))]
Q[coord,Nrow] <- betaCovariances[index:((index-1)+length(coord))]

index = index+sum(Nrow:ncol(Q)%%2==0)
Nrow <- 6
coord <- (Nrow:ncol(Q))[(Nrow:ncol(Q))%%2==0]
Q[Nrow,coord] <- betaCovariances[index:((index-1)+length(coord))]
Q[coord,Nrow] <- betaCovariances[index:((index-1)+length(coord))]

index = index+sum(Nrow:ncol(Q)%%2==0)
Nrow <- 8
coord <- (Nrow:ncol(Q))[(Nrow:ncol(Q))%%2==0]
Q[Nrow,coord] <- betaCovariances[index:((index-1)+length(coord))]
Q[coord,Nrow] <- betaCovariances[index:((index-1)+length(coord))]

index = index+sum(Nrow:ncol(Q)%%2==0)
Nrow <- 10
coord <- (Nrow:ncol(Q))[(Nrow:ncol(Q))%%2==0]
Q[Nrow,coord] <- betaCovariances[index:((index-1)+length(coord))]
Q[coord,Nrow] <- betaCovariances[index:((index-1)+length(coord))]
diag(Q) <- diag_Q
#Q <- "unconstrained"
#Q <- "diagonal and unequal"
B <- diag(m*Nspecies)
#B <- "unconstrained"
R <- "unconstrained"
#R <- diag(0.01,Nspecies)
x0 <- "zero"
initx <- 0
model.listDLMspecies <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,d=d,D=D,C=C2,c=c,x0=x0,tinitx=initx)
ln_RS_sh <- log(recruits$Recruits.Steelhead/adults$Stock.Steelhead)
ln_RS_dv <- log(recruits$`Recruits.Dolly Varden`/adults$`Stock.Dolly Varden`)
ln_RS_ct <- log(recruits$Recruits.Cutthroat/adults$Stock.Cutthroat)
ln_RS_pk <- log(recruits$Recruits.Pink/adults$Stock.Pink)
ln_RS_co <- log(recruits$Recruits.Coho/adults$Stock.Coho)
dat <- rbind(ln_RS_sh,ln_RS_dv,ln_RS_ct,ln_RS_pk,ln_RS_co)
a <- apply(dat,1,mean,na.rm=TRUE)
inits.list <- list(x0=array(matrix(rep(0,m), nrow=m),c(m,m,Nspecies)),A=matrix(a,Nspecies,1))
keoghDLMspecies <- MARSS(dat, model=model.listDLMspecies,control=list(maxit=250000,conv.test.slope.tol=0.5),inits=inits.list)
saveRDS(keoghDLMspecies,"Results/keoghDLM.rds")

keoghAllfit <- augment(keoghDLMspecies, interval="confidence")
keoghAllfit$Year <- keoghAllfit$t + 1975
keoghAllfit$Species <- keoghAllfit$.rownames
p <- ggplot(data = keoghAllfit) + 
  geom_line(aes(Year, .fitted)) +
  geom_ribbon(aes(x=Year, ymin=.conf.low, ymax=.conf.up), linetype=2, alpha=0.5)
p <- p + geom_point(data=keoghAllfit, mapping = aes(x=Year, y=y))
p + xlab("") + ylab("ln_RS") + facet_wrap(~Species)

