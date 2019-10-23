source("some functions.R")
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(MARSS)
library(broom)
library(wesanderson)
library(devtools)
install_github(repo="nwfsc-timeseries/atsar")
install.packages("atsar")
library(atsar)

keogh <- readRDS("Keogh_newJuv_enviro.rds")
keogh_long <- subset(keogh,Year<=2015 & Year>=1976)
keogh_long <- subset(keogh_long,Species!="Chum")

adults <- subset(keogh_long,select = c(Year,Species,Stock))
adults <- reshape(adults,direction = "wide",idvar="Year",timevar="Species")
recruits <- subset(keogh_long,select = c(Year,Species,Recruits))
recruits <- reshape(recruits,direction = "wide",idvar="Year",timevar="Species")

marSurv <- subset(keogh_long,select = c(Year,Species,Stock,juvCohort))
marSurv$Surv <- pmax(1e-3,pmin(1-1e-3,marSurv$Stock/marSurv$juvCohort),na.rm=TRUE)
marSurv$logitSurv <- log(marSurv$Surv/(1-marSurv$Surv))
marSurv <- subset(marSurv,select = c(Year,Species,logitSurv))
ocean_survival <- reshape(marSurv,direction = "wide",idvar="Year",timevar="Species")

summary(lm(ocean_survival$logitSurv.Steelhead~keogh_long$seals[keogh_long$Species=="Steelhead"]+keogh_long$oceanSalmon[keogh_long$Species=="Steelhead"]))

juv_enviro <- subset(keogh_long,select = c(Year,Species,sumTemp,sumRain,winTemp,winRain,freshCoho))
fresh_enviro <- reshape(juv_enviro,direction = "wide",idvar="Year",timevar="Species")
#fresh_enviro <- fresh_enviro[,-match(c("freshSteel.Pink","freshDolly.Pink","freshCutt.Pink","freshPink.Pink","freshCoho.Pink"),colnames(fresh_enviro))]
adult_enviro <- subset(keogh_long,select = c(Year,Species,seals,npgo,oceanSalmon))
adult_enviro$oceanSalmon <- residuals(lm(oceanSalmon~seals:Species,data=adult_enviro))

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

# all environmental predictors
all_enviro <- subset(keogh_long,select = c(Year,Species,sumTemp,winRain,freshCoho,seals,npgo,oceanSalmon))
all_enviro$seals <- residuals(lm(seals~oceanSalmon:Species,data=all_enviro))
#all_enviro$seals <- log(all_enviro$seals)
all_covars <- reshape(all_enviro,direction = "wide",idvar="Year",timevar="Species")
allEnviroNew <- all_covars
sdCovarsAll <- attr(scale(all_covars[,-1],center=TRUE,scale=TRUE),"scaled:scale")
mnCovarsAll <- attr(scale(all_covars[,-1],center=TRUE,scale=TRUE),"scaled:center")
allCovarScale <- scale(all_covars[,-1],center=TRUE,scale=TRUE)

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

Z1 <- array(0,c(Nspecies,Nspecies,Nyears))
Z1[1,1,] <- adults$Stock.Steelhead # time-varying beta
Z1[2,2,] <- adults$`Stock.Dolly Varden` # time-varying beta
Z1[3,3,] <- adults$Stock.Cutthroat # time-varying beta
Z1[4,4,] <- adults$Stock.Pink # time-varying beta
Z1[5,5,] <- adults$Stock.Coho # time-varying beta

Z2 <- array(0,c(Nspecies,Nspecies,Nyears))
Z2[1,1,] <- rep(1,Nyears) # time-varying alpha
Z2[2,2,] <- rep(1,Nyears) # time-varying alpha
Z2[3,3,] <- rep(1,Nyears) # time-varying alpha
Z2[4,4,] <- rep(1,Nyears) # time-varying alpha
Z2[5,5,] <- rep(1,Nyears) # time-varying alpha

# species competition, summer climate, and winter climate affects process model
d <- rbind(t(freshCovarScale))
D <- matrix(list(0),nrow=nrow(d),ncol=Nspecies)
coefNames <- paste("b",matrix(unlist(strsplit(row.names(d),"\\.")),ncol=2,byrow=TRUE)[,2],matrix(unlist(strsplit(row.names(d),"\\.")),ncol=2,byrow=TRUE)[,1],sep="_")
D[grep("Steelhead",coefNames),1] <- coefNames[grep("Steelhead",coefNames)]
D[grep("Dolly Varden",coefNames),2] <- coefNames[grep("Dolly Varden",coefNames)]
D[grep("Cutthroat",coefNames),3] <- coefNames[grep("Cutthroat",coefNames)]
D[grep("Pink_",coefNames),4] <- coefNames[grep("Pink_",coefNames)]
D[grep("Coho_",coefNames),5] <- coefNames[grep("Coho_",coefNames)]


d2 <- rbind(t(freshCovarScale))
D2 <- matrix(list(0),nrow=nrow(d2),ncol=Nspecies)
coefNames <- paste("b",matrix(unlist(strsplit(row.names(d2),"\\.")),ncol=2,byrow=TRUE)[,2],matrix(unlist(strsplit(row.names(d2),"\\.")),ncol=2,byrow=TRUE)[,1],sep="_")
D2[grep("Steelhead",coefNames),1] <- coefNames[grep("Steelhead",coefNames)]
D2[grep("Dolly Varden",coefNames),2] <- coefNames[grep("Dolly Varden",coefNames)]
D2[grep("Cutthroat",coefNames),3] <- coefNames[grep("Cutthroat",coefNames)]
D2[grep("Pink_",coefNames)[!grepl("fresh",coefNames[grep("Pink_",coefNames)])],4] <- coefNames[grep("Pink_",coefNames)][!grepl("fresh",coefNames[grep("Pink_",coefNames)])]
D2[grep("Coho_",coefNames),5] <- coefNames[grep("Coho_",coefNames)]

d3 <- rbind(adults$Stock.Steelhead,
            adults$`Stock.Dolly Varden`,
            adults$Stock.Cutthroat,
            adults$Stock.Pink,
            adults$Stock.Coho)
row.names(d3) <- c("stock.Steelhead","stock.Dolly Varden","stock.Cutthroat","stock.Pink","stock.Coho")
#d3 <- rbind(d3,d2)
D3 <- matrix(list(0),nrow=nrow(d3),ncol=Nspecies)
coefNames <- paste("b",matrix(unlist(strsplit(row.names(d3),"\\.")),ncol=2,byrow=TRUE)[,2],matrix(unlist(strsplit(row.names(d3),"\\.")),ncol=2,byrow=TRUE)[,1],sep="_")
D3[grep("Steelhead",coefNames),1] <- coefNames[grep("Steelhead",coefNames)]
D3[grep("Dolly Varden",coefNames),2] <- coefNames[grep("Dolly Varden",coefNames)]
D3[grep("Cutthroat",coefNames),3] <- coefNames[grep("Cutthroat",coefNames)]
D3[grep("Pink_",coefNames)[!grepl("fresh",coefNames[grep("Pink_",coefNames)])],4] <- coefNames[grep("Pink_",coefNames)][!grepl("fresh",coefNames[grep("Pink_",coefNames)])]
D3[grep("Coho_",coefNames),5] <- coefNames[grep("Coho_",coefNames)]

d4 <- rbind(rep(1,Nyears),
            rep(1,Nyears),
            rep(1,Nyears),
            rep(1,Nyears),
            rep(1,Nyears))
row.names(d4) <- c("alpha.Steelhead","alpha.Dolly Varden","alpha.Cutthroat","alpha.Pink","alpha.Coho")
#d4 <- rbind(d4,d2)
D4 <- matrix(list(0),nrow=nrow(d4),ncol=Nspecies)
coefNames <- paste("b",matrix(unlist(strsplit(row.names(d4),"\\.")),ncol=2,byrow=TRUE)[,2],matrix(unlist(strsplit(row.names(d4),"\\.")),ncol=2,byrow=TRUE)[,1],sep="_")
D4[grep("Steelhead",coefNames),1] <- coefNames[grep("Steelhead",coefNames)]
D4[grep("Dolly Varden",coefNames),2] <- coefNames[grep("Dolly Varden",coefNames)]
D4[grep("Cutthroat",coefNames),3] <- coefNames[grep("Cutthroat",coefNames)]
D4[grep("Pink_",coefNames)[!grepl("fresh",coefNames[grep("Pink_",coefNames)])],4] <- coefNames[grep("Pink_",coefNames)][!grepl("fresh",coefNames[grep("Pink_",coefNames)])]
D4[grep("Coho_",coefNames),5] <- coefNames[grep("Coho_",coefNames)]

# time-constant model
d5 <- rbind(rep(1,Nyears),
            adults$Stock.Steelhead,
            rep(1,Nyears),
            adults$`Stock.Dolly Varden`,
            rep(1,Nyears),
            adults$Stock.Cutthroat,
            rep(1,Nyears),
            adults$Stock.Pink,
            rep(1,Nyears),
            adults$Stock.Coho)
row.names(d5) <- c("alpha.Steelhead","stock.Steelhead",
                   "alpha.Dolly Varden","stock.Dolly Varden",
                   "alpha.Cutthroat","stock.Cutthroat",
                   "alpha.Pink","stock.Pink",
                   "alpha.Coho","stock.Coho")

#d3 <- rbind(d3,d2)
D5 <- matrix(list(0),nrow=nrow(d5),ncol=Nspecies)
coefNames <- paste("b",matrix(unlist(strsplit(row.names(d5),"\\.")),ncol=2,byrow=TRUE)[,2],matrix(unlist(strsplit(row.names(d5),"\\.")),ncol=2,byrow=TRUE)[,1],sep="_")
D5[grep("Steelhead",coefNames),1] <- coefNames[grep("Steelhead",coefNames)]
D5[grep("Dolly Varden",coefNames),2] <- coefNames[grep("Dolly Varden",coefNames)]
D5[grep("Cutthroat",coefNames),3] <- coefNames[grep("Cutthroat",coefNames)]
D5[grep("Pink_",coefNames)[!grepl("fresh",coefNames[grep("Pink_",coefNames)])],4] <- coefNames[grep("Pink_",coefNames)][!grepl("fresh",coefNames[grep("Pink_",coefNames)])]
D5[grep("Coho_",coefNames),5] <- coefNames[grep("Coho_",coefNames)]

# time-constant beta
# NPGO, seals, pacific salmon, and mei affects observation model
c <- rbind(t(oceanCovarScale))
C <- matrix(list(0),nrow=m*Nspecies,ncol=nrow(c))
alphaCoefNames <- paste("a",matrix(unlist(strsplit(row.names(c),"\\.")),ncol=2,byrow=TRUE)[,2],matrix(unlist(strsplit(row.names(c),"\\.")),ncol=2,byrow=TRUE)[,1],sep="_")
betaCoefNames <- paste("b",matrix(unlist(strsplit(row.names(c),"\\.")),ncol=2,byrow=TRUE)[,2],matrix(unlist(strsplit(row.names(c),"\\.")),ncol=2,byrow=TRUE)[,1],sep="_")
C[1:2,grep("Steelhead",alphaCoefNames)] <- rbind(alphaCoefNames[grep("Steelhead",alphaCoefNames)],
                                                 betaCoefNames[grep("Steelhead",betaCoefNames)])
C[3:4,grep("Dolly Varden",alphaCoefNames)] <- rbind(alphaCoefNames[grep("Dolly Varden",alphaCoefNames)],
                                                  betaCoefNames[grep("Dolly Varden",betaCoefNames)])
C[5:6,grep("Cutthroat",alphaCoefNames)] <- rbind(alphaCoefNames[grep("Cutthroat",alphaCoefNames)],
                                               betaCoefNames[grep("Cutthroat",betaCoefNames)])
C[7:8,grep("Pink_",alphaCoefNames)] <- rbind(alphaCoefNames[grep("Pink_",alphaCoefNames)],
                                           betaCoefNames[grep("Pink_",betaCoefNames)])
C[9:10,grep("Coho_",alphaCoefNames)] <- rbind(alphaCoefNames[grep("Coho_",alphaCoefNames)],
                                           betaCoefNames[grep("Coho_",betaCoefNames)])

c <- rbind(t(oceanCovarScale))
C2 <- matrix(list(0),nrow=m*Nspecies,ncol=nrow(c))
C2[2,grep("Steelhead",betaCoefNames)] <- betaCoefNames[grep("Steelhead",betaCoefNames)]
C2[4,grep("Dolly Varden",betaCoefNames)] <- betaCoefNames[grep("Dolly Varden",betaCoefNames)]
C2[6,grep("Cutthroat",betaCoefNames)] <- betaCoefNames[grep("Cutthroat",betaCoefNames)]
C2[8,grep("Pink_",betaCoefNames)] <- betaCoefNames[grep("Pink_",betaCoefNames)]
C2[10,grep("Coho_",betaCoefNames)] <- betaCoefNames[grep("Coho_",betaCoefNames)]

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

Q2 <- matrix(list(0),m*Nspecies,m*Nspecies)        ## 2x2; all 0 for now
a_bCovariances <- paste("corr_ab",unique(keogh_long$Species),sep="_")
Q2[1,2] <- a_bCovariances[1]
Q2[2,1] <- a_bCovariances[1]

Q2[3,4] <- a_bCovariances[2]
Q2[4,3] <- a_bCovariances[2]

Q2[5,6] <- a_bCovariances[3]
Q2[6,5] <- a_bCovariances[3]

Q2[7,8] <- a_bCovariances[4]
Q2[8,7] <- a_bCovariances[4]

Q2[9,10] <- a_bCovariances[5]
Q2[10,9] <- a_bCovariances[5]
diag(Q2) <- diag_Q
#Q <- "unconstrained"
#Q <- "diagonal and unequal"
B <- diag(m*Nspecies)
#B <- "unconstrained"
R <- "unconstrained"
#R <- diag(0.01,Nspecies)
x0 <- "zero"
initx <- 0
model.listDLMspecies <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,d=d,D=D,C=C2,c=c,x0=x0,tinitx=initx)
model.listDLMv2 <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,d=d2,D=D2,x0=x0,tinitx=initx)

model.listTVab <- list(B=diag(m*Nspecies),U=U,Q=Q,Z=Z,A=A,R=R,x0=x0,tinitx=initx)
model.listTVa <- list(B=diag(Nspecies),U=U,Q="unconstrained",Z=Z2,A=A,R=R,d=d3,D=D3,x0=x0,tinitx=initx)
model.listTVb <- list(B=diag(Nspecies),U=U,Q="unconstrained",Z=Z1,A=A,R=R,d=d4,D=D4,x0=x0,tinitx=initx)
model.listConab <- list(B="identity",Z="identity",U="zero",Q="zero",A=A,R=R,d=d5,D=D5,x0="zero",tinitx=initx)

ln_RS_sh <- log(recruits$Recruits.Steelhead/adults$Stock.Steelhead)
ln_RS_dv <- log(recruits$`Recruits.Dolly Varden`/adults$`Stock.Dolly Varden`)
ln_RS_ct <- log(recruits$Recruits.Cutthroat/adults$Stock.Cutthroat)
ln_RS_pk <- log(recruits$Recruits.Pink/adults$Stock.Pink)
ln_RS_co <- log(recruits$Recruits.Coho/adults$Stock.Coho)
dat <- rbind(ln_RS_sh,ln_RS_dv,ln_RS_ct,ln_RS_pk,ln_RS_co)

a <- apply(dat,1,mean,na.rm=TRUE)
inits.listTVab <- list(x0=array(matrix(rep(0,m), nrow=m),c(m,m,Nspecies)),A=matrix(a,Nspecies,1))
inits.listTVa <- list(x0=array(matrix(rep(0,1), nrow=1),c(1,1,Nspecies)),A=matrix(a,Nspecies,1))
inits.listTVb <- list(x0=array(matrix(rep(0,1), nrow=1),c(1,1,Nspecies)),A=matrix(a,Nspecies,1))
inits.listConab <- list(A=matrix(a,Nspecies,1))

begin <- Sys.time()
begin
keoghDLMTVab <- MARSS(dat, model=model.listTVab,control=list(maxit=10000,conv.test.slope.tol=0.3),inits=inits.listTVab)
end <- Sys.time()
time_elapse <- end-begin
time_elapse

begin <- Sys.time()
begin
keoghDLMTVa <- MARSS(dat, model=model.listTVa,control=list(maxit=10000,conv.test.slope.tol=0.3),inits=inits.listTVa)
end <- Sys.time()
time_elapse <- end-begin
time_elapse

begin <- Sys.time()
begin
keoghDLMTVb <- MARSS(dat, model=model.listTVb,control=list(maxit=10000,conv.test.slope.tol=0.3),inits=inits.listTVb)
end <- Sys.time()
time_elapse <- end-begin
time_elapse

begin <- Sys.time()
begin
keoghConab <- MARSS(dat, model=model.listConab,control=list(maxit=10000,conv.test.slope.tol=0.3),inits=inits.listConab)
end <- Sys.time()
time_elapse <- end-begin
time_elapse

AIC(keoghDLMTVab,keoghDLMTVa,keoghDLMTVb,keoghConab)
AICres <- data.frame("Model"=c("time-varying a,b","time-varying a","time-varying b","time-constant a,b"),"df"=c(keoghDLMTVab$num.params,keoghDLMTVa$num.params,keoghDLMTVb$num.params,keoghConab$num.params),"AICc"=c(keoghDLMTVab$AICc,keoghDLMTVa$AICc,keoghDLMTVb$AICc,keoghConab$AICc))
AICres$dAIC <- AICres$AICc-min(AICres$AICc)
AICres[order(AICres$dAIC),]
saveRDS(keoghDLMTVb,"Results/keoghDLM_tvar.rds")

# run time-varying beta model with covariates
dCovar <- rbind(rep(1,Nyears),
            rep(1,Nyears),
            rep(1,Nyears),
            rep(1,Nyears),
            rep(1,Nyears))
row.names(dCovar) <- c("alpha.Steelhead","alpha.Dolly Varden","alpha.Cutthroat","alpha.Pink","alpha.Coho")
dCovar <- rbind(dCovar[1,],d2[grep("\\.Steelhead",row.names(d2)),],
                dCovar[2,],d2[grep("\\.Dolly Varden",row.names(d2)),],
                dCovar[3,],d2[grep("\\.Cutthroat",row.names(d2)),],
                dCovar[4,],d2[grep("\\.Pink",row.names(d2)),],
                dCovar[5,],d2[grep("\\.Coho",row.names(d2)),])
row.names(dCovar)[row.names(dCovar)==""] <- c("alpha.Steelhead","alpha.Dolly Varden","alpha.Cutthroat","alpha.Pink","alpha.Coho")
Dcovar <- matrix(list(0),nrow=nrow(dCovar),ncol=Nspecies)
coefNames <- paste("b",matrix(unlist(strsplit(row.names(dCovar),"\\.")),ncol=2,byrow=TRUE)[,2],matrix(unlist(strsplit(row.names(dCovar),"\\.")),ncol=2,byrow=TRUE)[,1],sep="_")
Dcovar[grep("Steelhead",coefNames),1] <- coefNames[grep("Steelhead",coefNames)]
Dcovar[grep("Dolly Varden",coefNames),2] <- coefNames[grep("Dolly Varden",coefNames)]
Dcovar[grep("Cutthroat",coefNames),3] <- coefNames[grep("Cutthroat",coefNames)]
Dcovar[grep("Pink_",coefNames)[!grepl("fresh",coefNames[grep("Pink_",coefNames)])],4] <- coefNames[grep("Pink_",coefNames)][!grepl("fresh",coefNames[grep("Pink_",coefNames)])]
Dcovar[grep("Coho_",coefNames),5] <- coefNames[grep("Coho_",coefNames)]

# create a C matrix just for the time-varying beta model
c3 <- rbind(t(oceanCovarScale))
C3 <- matrix(list(0),nrow=Nspecies,ncol=nrow(c3))
coefNames <- paste("b",matrix(unlist(strsplit(row.names(c3),"\\.")),ncol=2,byrow=TRUE)[,2],matrix(unlist(strsplit(row.names(c3),"\\.")),ncol=2,byrow=TRUE)[,1],sep="_")
C3[1,grep("Steelhead",coefNames)] <- coefNames[grep("Steelhead",coefNames)]
C3[2,grep("Dolly Varden",coefNames)] <- coefNames[grep("Dolly Varden",coefNames)]
C3[3,grep("Cutthroat",coefNames)] <- coefNames[grep("Cutthroat",coefNames)]
C3[4,grep("Pink_",coefNames)] <- coefNames[grep("Pink_",coefNames)]
C3[5,grep("Coho_",coefNames)] <- coefNames[grep("Coho_",coefNames)]

# create a model matrix for the time-varying beta model with covariates
model.listTVbCovar <- list(B=diag(Nspecies),U=U,Q="unconstrained",Z=Z1,A=A,R=R,c=c3,C=C3,d=dCovar,D=Dcovar,x0=x0,tinitx=initx)
begin <- Sys.time()
begin
keoghDLMspecies <- MARSS(dat, model=model.listTVbCovar,control=list(maxit=10000,conv.test.slope.tol=0.3),inits=inits.listTVb)
end <- Sys.time()
time_elapse <- end-begin
time_elapse

saveRDS(keoghDLMspecies,"Results/keoghDLM.rds")

# get 95% CIS for model coefficients for bootstrap
begin <- Sys.time()
begin
CIs <- MARSSboot(keoghDLMspecies,nboot = 250,param.gen="MLE")
end <- Sys.time()
time_elapse <- end-begin
time_elapse
CI_95s <- apply(CIs$boot.params,1,quantile,probs=c(0.1,0.9))
CI_med <- apply(CIs$boot.params,1,median)
covarDD <- rbind(CI_med[grep("U.b_",colnames(CI_95s))],CI_95s[,grep("U.b_",colnames(CI_95s))])
sd_b <- sqrt(diag(coef(keoghDLMspecies,type="matrix")$Q))
sd_b <- sapply(1:Nspecies,function(x){sd(keoghDLMspecies$states[x,])})
names(sd_b) <- unique(keogh_long$Species)
std_beta <- sapply(1:ncol(covarDD),function(x){covarDD[,x]/rep(sd_b,times=c(3,3,3,3,3))[x]})
colnames(std_beta) <- colnames(covarDD)

varNames <- apply(matrix(unlist(strsplit(colnames(std_beta),"_")),ncol=3,byrow=TRUE)[,2:3],1,paste,collapse=" | ",sep="")
jpeg("effect sizes.jpeg",width=8,height=6,units="in",res=800)
layout(1)
par(mar=c(5,12,3,3))
plot(std_beta[1,],1:length(colnames(std_beta)),xlim=c(-1,1),pch=21,bg="dodgerblue",yaxt="n",xlab="Effect size",ylab="",cex=1.2)
#grid(NULL,NULL,lty=6,col="grey70")
abline(v=c(-4,-2,0,2,4),lty=6,col="grey70")
abline(h=c(3.5,6.5,9.5,12.5),lty=6,col="grey70")
segments(x0=std_beta[2,],y0=1:length(colnames(std_beta)),x1=std_beta[3,],col="black",lwd=2)
abline(v=0,lwd=2,lty=2,col="black")
sig <- sapply(1:ncol(std_beta),function(x){any(all(std_beta[,x] > 0) | all(std_beta[,x] < 0)) })
points(std_beta[1,!sig],(1:ncol(std_beta))[!sig],pch=21,bg="dodgerblue",cex=1.2)
points(std_beta[1,sig],(1:ncol(std_beta))[sig],pch=21,bg="tomato",cex=1.5)
axis(2,at=1:ncol(std_beta),labels=varNames,las=2,cex.axis=0.8)
text(0.5,16,"Reduced density-dependence",xpd=NA,cex=0.7,font=2)
text(-0.5,16,"Increased density-dependence",xpd=NA,cex=0.7,font=2)
dev.off()

covarDD <- rbind(CI_med[grep("A.b_",colnames(CI_95s))],CI_95s[,grep("A.b_",colnames(CI_95s))])
covarDD <- covarDD[,-grep("alpha",colnames(covarDD))]
sd_b <- sapply(1:Nspecies,function(x){sd(dat[x,],na.rm=TRUE)})
names(sd_b) <- unique(keogh_long$Species)
std_beta <- sapply(1:ncol(covarDD),function(x){covarDD[,x]/rep(sd_b,times=c(5,5,5,5,5))[x]})
colnames(std_beta) <- colnames(covarDD)

varNames <- apply(matrix(unlist(strsplit(colnames(std_beta),"_")),ncol=3,byrow=TRUE)[,2:3],1,paste,collapse=" | ",sep="")
jpeg("effect sizes on productivity.jpeg",width=8,height=6,units="in",res=800)
layout(1)
par(mar=c(5,12,3,3))
plot(std_beta[1,],1:length(colnames(std_beta)),xlim=c(-4,4),pch=21,bg="dodgerblue",yaxt="n",xlab="Effect size",ylab="",cex=1.2)
#grid(NULL,NULL,lty=6,col="grey70")
abline(v=c(-4,-2,0,2,4),lty=6,col="grey70")
abline(h=c(5.5,10.5,15.5,20.5),lty=6,col="grey70")
segments(x0=std_beta[2,],y0=1:length(colnames(std_beta)),x1=std_beta[3,],col="black",lwd=2)
abline(v=0,lwd=2,lty=2,col="black")
sig <- sapply(1:ncol(std_beta),function(x){any(all(std_beta[,x] > 0) | all(std_beta[,x] < 0)) })
points(std_beta[1,!sig],(1:ncol(std_beta))[!sig],pch=21,bg="dodgerblue",cex=1.2)
points(std_beta[1,sig],(1:ncol(std_beta))[sig],pch=21,bg="tomato",cex=1.5)
axis(2,at=1:ncol(std_beta),labels=varNames,las=2,cex.axis=0.8)
text(2,25.5,"Increased productivity",xpd=NA,cex=0.7,font=2)
text(-2,25.5,"Reduced productivity",xpd=NA,cex=0.7,font=2)
dev.off()

keoghDLMspecies$AICc
keoghDLMspecies$num.params

keoghAllfit <- augment(keoghDLMspecies, interval="confidence")
keoghAllfit$Year <- keoghAllfit$t + 1975
keoghAllfit$Species <- keoghAllfit$.rownames
p <- ggplot(data = keoghAllfit) + 
  geom_line(aes(Year, .fitted)) +
  geom_ribbon(aes(x=Year, ymin=.conf.low, ymax=.conf.up), linetype=2, alpha=0.5)
p <- p + geom_point(data=keoghAllfit, mapping = aes(x=Year, y=y))
p + xlab("") + ylab("ln_RS") + facet_wrap(~Species)

matplot(t(keoghDLMspecies$states))

plot(c3[4,],keoghDLMspecies$states[1,],pch=21,bg=ifelse(years>1990,"orange","dodgerblue"))
lines(c3[4,],c(0,keoghDLMspecies$states[1,1:(Nyears-1)]+coef(keoghDLMspecies)$C[3,1]*c3[4,2:Nyears]))
coef(keoghDLMspecies,type="matrix")$C

coef(keoghDLMspecies)$C
#
layout(matrix(1:6,nrow=3,ncol=2,byrow=TRUE))
par(mar=c(5,4,1,1))
plot(c3[grep("seals",row.names(c3))[1],],keoghDLMspecies$states[1,],pch=21,bg=ifelse(years>1990,"orange","dodgerblue"),xlab="Seals",ylab="Density-dependence")
plot(c3[grep("seals",row.names(c3))[2],],keoghDLMspecies$states[2,],pch=21,bg=ifelse(years>1990,"orange","dodgerblue"),xlab="Seals",ylab="Density-dependence")
plot(c3[grep("seals",row.names(c3))[3],],keoghDLMspecies$states[3,],pch=21,bg=ifelse(years>1990,"orange","dodgerblue"),xlab="Seals",ylab="Density-dependence")
plot(c3[grep("seals",row.names(c3))[4],],keoghDLMspecies$states[4,],pch=21,bg=ifelse(years>1990,"orange","dodgerblue"),xlab="Seals",ylab="Density-dependence")
plot(c3[grep("seals",row.names(c3))[5],],keoghDLMspecies$states[5,],pch=21,bg=ifelse(years>1990,"orange","dodgerblue"),xlab="Seals",ylab="Density-dependence")

plot(c3[1,],keoghDLMspecies$states[1,],pch=21,bg=ifelse(years>1990,"orange","dodgerblue"))
plot(c3[4,],keoghDLMspecies$states[1,],pch=21,bg=ifelse(years>1990,"orange","dodgerblue"))
plot(c3[6,],keoghDLMspecies$states[1,],pch=21,bg=ifelse(years>1990,"orange","dodgerblue"))

plot(c3[27,],keoghDLMspecies$states[5,],pch=21,bg=ifelse(years>1990,"orange","dodgerblue"))
plot(c3[29,],keoghDLMspecies$states[5,],pch=21,bg=ifelse(years>1990,"orange","dodgerblue"))

coef(keoghDLMspecies)$C
plot(c3[3,],keoghDLMspecies$states[2,],pch=21,bg=ifelse(years>1990,"orange","dodgerblue"))
plot(c3[10,],keoghDLMspecies$states[2,],pch=21,bg=ifelse(years>1990,"orange","dodgerblue"))
plot(c3[12,],keoghDLMspecies$states[2,],pch=21,bg=ifelse(years>1990,"orange","dodgerblue"))


plot(keoghDLMspecies$states[1,1]+coef(keoghDLMspecies)$C[4,1]*c3[4,2:Nyears],keoghDLMspecies$states[1,2:(Nyears)])
abline(b=1,a=0)

plot(coef(keoghDLMspecies)$C[4,1]*c3[4,2:Nyears],keoghDLMspecies$states[1,2:(Nyears)])

plot(coef(keoghDLMspecies)$D[1]+keoghDLMspecies$states[1,]*adults$Stock.Steelhead,dat[1,])
abline(b=1,a=0)
plot(seals~Year,data=keogh_long)
