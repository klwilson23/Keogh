source("some functions.R")
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggplot2)
library(MARSS)
library(broom)
keogh <- readRDS("Data/Keogh_SR_enviro_long.rds")
keogh_long <- subset(keogh,Year<=2015 & Year>=1976)
keogh_long <- subset(keogh_long,Species!="Chum")

keogh <- subset(keogh_long,select = c(Year,Species,Stock,Recruits,juvCohort))
keogh <- reshape(keogh,direction = "wide",idvar="Year",timevar="Species")

environment <- subset(keogh_long,select = c(Year,Species,sumTemp,sumRain,winTemp,winRain,freshCoho,freshSteel,freshCutt,freshDolly,freshPink,seals,npgo,mei,oceanSalmon,Logging,cumul_log,cumul_footprint,meanLogging,cumLogging,fertil))
enviro <- reshape(environment,direction = "wide",idvar="Year",timevar="Species")

enviroNew <- enviro
sdCovars <- attr(scale(enviro[,-1],center=TRUE,scale=TRUE),"scaled:scale")
mnCovars <- attr(scale(enviro[,-1],center=TRUE,scale=TRUE),"scaled:center")
covarScale <- scale(enviro[,-1],center=TRUE,scale=TRUE)
#covarScale[is.na(covarScale)] <- 0

# get estimates of missing data from DLM analysis
Nyears <- length(enviro$Year)
years <- enviro$Year
covarNames <- c("sumTemp","sumRain","winTemp","winRain","freshCoho","freshSteel","freshCutt","freshDolly","freshPink","seals","npgo","mei","oceanSalmon","Logging","cumul_log","cumul_footprint","meanLogging","cumLogging","fertil")

for(i in 1:length(covarNames))
{
  covarSub <- covarScale[,grep(covarNames[i],colnames(covarScale))]
  covars <- t(as.matrix(covarSub))
  colnames(covars) <- enviro$Year
  if(any(is.na(enviro[,grep(covarNames[i],colnames(enviro))])))
  {
    ns <- nrow(covars)
    nTrends <- 2 # one for trout and one for salmon
    B <- matrix(list(0), nTrends, nTrends)
    diag(B) <- paste("b",1:nTrends,sep="")
    Q <- diag(1, nTrends)
    #R <- "diagonal and unequal"
    R <- diag(0.01,ns)
    U <- "zero"
    x0 <- "zero"
    Z <- matrix(list(0), ns, nTrends)
    Z[1:(ns * nTrends)] <- sapply(1:nTrends,function(x){paste0("z",x,1:ns)})
    Z[upper.tri(Z)] <- 0
    A <- "unequal"
    mod.list.dfa = list(B = B, Z = Z, Q = Q, R = R, U = U, A = A, x0 = x0)
    m <- apply(covars, 1, mean, na.rm=TRUE)
    fit.dfa <- MARSS(covars, model = mod.list.dfa, control = list(maxit = 50000), inits = list(A = matrix(m, ns, 1)))
    
    d <- fitted(fit.dfa,interval="confidence")
    d$Year <- d$t + (min(years)-1)
    d$covars <- d$.rownames
    spp <- matrix(unlist(strsplit(as.character(d$.rownames),".",fixed=TRUE)),ncol=2,byrow=TRUE)[,2]
    d$Species <- factor(spp,levels=levels(keogh_long$Species))
  
    covarFits <- matrix(sapply(1:nrow(d),function(x){matches <- names(sdCovars)%in%d$covars[x];
      d$.fitted[x]*sdCovars[matches] + mnCovars[matches]
    }),nrow=ncol(covarSub),ncol=nrow(covarSub),byrow=FALSE)
    
    covarEst <- as.data.frame(t(covarFits))
    enviroNew[,grep(covarNames[i],colnames(enviroNew))][which(is.na(enviroNew[,grep(covarNames[i],colnames(enviroNew))]),arr.ind=TRUE)] <- covarEst[which(is.na(enviroNew[,grep(covarNames[i],colnames(enviroNew))]),arr.ind=TRUE)] 
  
    yy <- data.frame("Year"=enviro$Year,covarSub)
    yy <- reshape2::melt(yy, id.vars=c("Year"))
    yy$Species <- d$Species
    colnames(yy) <- c("Year","Covar","Covariate","Species")
    
    p <- ggplot(data = d) + 
      geom_line(aes(Year, .fitted)) + 
      geom_ribbon(aes(x = Year, ymin = .conf.low, ymax = .conf.up), linetype = 2, alpha = 0.5)
    p <- p + geom_point(data = yy, mapping = aes(x = Year, y = Covariate))
    p <- p + facet_wrap(~Species) + xlab("") + ylab(covarNames[i])
    print(p)
  }else{
    next
  }
}

enviroNew[,grep("fresh",colnames(enviroNew))] <- apply(enviroNew[,grep("fresh",colnames(enviroNew))],2,function(x){pmax(0,x,na.rm=TRUE)})
range(enviroNew[,grep("fresh",colnames(enviroNew))])

reEnviro <- reshape(enviroNew,direction="long",idvar="Year",timevar="Species")
colnames(reEnviro) <- colnames(environment)

keogh_SR <- subset(keogh_long,select = c(Year,Species,Stock,Recruits,juvCohort))
keogh_new <- data.frame(keogh_SR,reEnviro[,!colnames(reEnviro) %in% c("Year","Species")])

saveRDS(keogh_new,file="Data/Keogh_SR_enviro_new.rds")
