rm(list=ls(all=TRUE))
# data management sources
# part 1 - get forestry data into annual time-series
lag_seq <- c(5,10,15,30)
beta_logging <- list()
for(i in 1:length(lag_seq))
{
  lag_forest <- lag_seq[i]
  source("some functions.R")
  library(reshape2)
  library(gridExtra)
  library(ggpubr)
  library(ggplot2)
  library(MARSS)
  library(broom)
  library(wesanderson)
  library(rstan)
  library(loo)
  library(corrplot)
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores(logical=FALSE))
  #Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7')
  keogh <- readRDS(paste("Data/Keogh_collinear_enviro_forestlag_",lag_forest,".rds",sep=""))
  keogh_long <- subset(keogh,Year<=2015 & Year>=1976)
  keogh_long <- subset(keogh_long,Species!="Chum")
  keogh_long$Species <- factor(keogh_long$Species,levels=unique(keogh_long$Species))
  keogh_long$marSurv <- keogh_long$Stock/keogh_long$juvCohort
  keogh_long$logitSurv <- log(keogh_long$marSurv/(1-keogh_long$marSurv))
  keogh_long$prod <- log(keogh_long$Recruits/keogh_long$Stock)
  
  keogh_sel <- subset(keogh_long,select = c(Year,Species,Stock,prod))
  keogh_wide <- reshape(keogh_sel,direction = "wide",idvar="Year",timevar="Species")
  
  run_time <- readRDS("Data/steelhead_run.rds")
  sh_annual <- readRDS("Data/steelhead_run_annual.rds")
  sh_annual$time <- 1:length(sh_annual$year)
  sh_annual <- subset(sh_annual,year<=2015 & year>=1976)
  sh_annual$logit_surv[1] <- NA
  sh_annual$log_adults <- log(sh_annual$Stock)
  
  Xvars <- c("ocean_interact","npgo","ocean_covar_2")
  sdSurv_sh <- attr(scale(sh_annual[,Xvars],center=TRUE,scale=TRUE),"scaled:scale")
  
  mnSurv_sh <- attr(scale(sh_annual[,Xvars],center=TRUE,scale=TRUE),"scaled:center")
  enviro <- scale(sh_annual[,Xvars],center=TRUE,scale=TRUE)
  enviro <- data.frame(Xvars=enviro)
  colnames(enviro) <- Xvars
  x1 <- model.matrix(~-1+ocean_interact+npgo+ocean_covar_2,data=enviro)
  
  XXvars <- c("total_rain_run","mean_temp_run")
  sdSurv_run <- attr(scale(sh_annual[,XXvars],center=TRUE,scale=TRUE),"scaled:scale")
  mnSurv_run <- attr(scale(sh_annual[,XXvars],center=TRUE,scale=TRUE),"scaled:center")
  enviro_run <- scale(sh_annual[,XXvars],center=TRUE,scale=TRUE)
  enviro_run <- data.frame(enviro_run)
  colnames(enviro_run) <- XXvars
  x2 <- model.matrix(~-1+total_rain_run+mean_temp_run,data=enviro_run)
  
  XXXvars <- c("meanLogging","total_rain_egg","mean_temp_egg","freshPink","fertil")
  sdSurv_prod <- attr(scale(sh_annual[,XXXvars],center=TRUE,scale=TRUE),"scaled:scale")
  mnSurv_prod <- attr(scale(sh_annual[,XXXvars],center=TRUE,scale=TRUE),"scaled:center")
  enviro_prod <- scale(sh_annual[,XXXvars],center=TRUE,scale=TRUE)
  enviro_prod <- data.frame(XXXvars=enviro_prod)
  colnames(enviro_prod) <- XXXvars
  
  x <- as.matrix(keogh_wide[,grep("Stock",colnames(keogh_wide))])
  y <- as.matrix(keogh_wide[,grep("prod",colnames(keogh_wide))])
  y[which(is.na(y),arr.ind=TRUE)] <- colMeans(y[1:10,],na.rm=TRUE)[which(is.na(y),arr.ind=TRUE)[,2]]
  
  xx1 <- scale(model.matrix(~-1+meanLogging+total_rain_egg+mean_temp_egg+freshPink+fertil,data=enviro_prod),center=TRUE,scale=TRUE)
  xx2 <- scale(model.matrix(~-1+meanLogging+sumTemp+sumRain+winTemp+winRain+freshPink+fertil,data=keogh_long[keogh_long$Species=="Dolly Varden",]),center=TRUE,scale=TRUE)
  xx3 <- scale(model.matrix(~-1+meanLogging+sumTemp+sumRain+winTemp+winRain+freshPink+fertil,data=keogh_long[keogh_long$Species=="Cutthroat",]),center=TRUE,scale=TRUE)
  xx4 <- scale(model.matrix(~-1+meanLogging+sumTemp+sumRain+ocean_interact+ocean_covar_2+npgo+winTemp+winRain,data=keogh_long[keogh_long$Species=="Pink",]),center=TRUE,scale=TRUE)
  xx5 <- scale(model.matrix(~-1+meanLogging+sumTemp+sumRain+winTemp+winRain+freshPink+fertil,data=keogh_long[keogh_long$Species=="Coho",]),center=TRUE,scale=TRUE)
  
  
  # all models
  dat <- list("N"=nrow(x),
              "K"=ncol(x),
              "x"=x,
              "y"=y,
              "init_priors"=rep(-2e-3,ncol(x)),
              "J1"=ncol(xx1),
              "J2"=ncol(xx2),
              "J3"=ncol(xx3),
              "J4"=ncol(xx4),
              "J5"=ncol(xx5),
              "xx1"=xx1,
              "xx2"=xx2,
              "xx3"=xx3,
              "xx4"=xx4,
              "xx5"=xx5,
              "N_obs"=sum(!is.na(sh_annual$logit_surv)),
              "M"=ncol(x1),
              "Q"=ncol(x2),
              "P"=2,
              "x1"=x1,
              "x2"=x2,
              "juvCoh"=as.vector(scale(sh_annual$juvCohort)),
              #"x3"=as.numeric(x3),
              "y1_obs"=sh_annual$logit_surv[!is.na(sh_annual$logit_surv)],
              "y2"=sh_annual$run_date_corrected,
              "init_s0"=mean(sh_annual$logit_surv[1:10],na.rm=TRUE))
  
  init_fx <- function(chain_id)
  {
    list("beta_steel"=rep(0,dat$J1+1),
         "beta_dolly"=rep(0,dat$J2),
         "beta_cutt"=rep(0,dat$J3),
         "beta_pink"=rep(0,dat$J4),
         "beta_coho"=rep(0,dat$J5),
         "beta_adults"=rep(0,dat$P))
  }
  
  fit <- stan(file="Stan code/Keogh mnorm MARSS and steelhead.stan",data=dat, iter=2000,chains=6,cores=6,control=list("adapt_delta"=0.85,"max_treedepth"=15),init=init_fx)
  intervals <- c(0.025,0.975)
  summ <- rstan::summary(fit,pars=c("beta_dolly","beta_steel","beta_cutt","beta_coho","beta_pink"),probs=intervals)$summary
  beta_logging[[i]] <- summ[grep("\\b[1]",row.names(summ)),]
  rm(fit);rm(summ)
}
names(beta_logging) <- paste("lag time",lag_seq)
saveRDS(beta_logging,file="Results/forestry lag sensitivity.rds")


usdm::vifstep(as.data.frame(xx1))
usdm::vifstep(as.data.frame(x1))
usdm::vifstep(as.data.frame(x2))
usdm::vifstep(as.data.frame(xx2))
usdm::vifstep(as.data.frame(xx3))
usdm::vifstep(as.data.frame(xx4))
usdm::vifstep(as.data.frame(xx5))
