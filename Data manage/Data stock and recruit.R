library(reshape2)
library(igraph)
library(ggplot2)
library(ggridges)
library(viridis)
library(weathercan)
library(zoo)
source("some functions.R")

# download weather canada data
# Download them separately for the whole time range
# (NAs on the ends they will be trimmed, as you saw)
s1 <- weather_dl(station_ids = 202, start = "1975-01-01", end = "2018-12-31",
                 interval = "day")
s2 <- weather_dl(station_ids = 51319, start = "1975-01-01", end = "2018-12-31",
                 interval = "day")
portHardy_pg <- rbind(s1, s2)
duration <- 14 # how many days of rain/temperature before steelhead make run
portHardy_pg$mean_temp_run <- roll_mean(portHardy_pg$mean_temp,duration)
portHardy_pg$max_temp_run <- roll_mean(portHardy_pg$max_temp,duration)
portHardy_pg$min_temp_run <- roll_mean(portHardy_pg$min_temp,duration)
portHardy_pg$total_rain_run <- roll_mean(portHardy_pg$total_rain,duration)

incubation <- 30 # how many days after the adult run are eggs laid
portHardy_pg$total_rain_egg <- roll_mean_forward(portHardy_pg$total_rain,incubation)
portHardy_pg$mean_temp_egg <- roll_mean_forward(portHardy_pg$mean_temp,incubation)

ggplot(data = portHardy_pg, aes(x = date, y = mean_temp_run, colour = factor(station_id))) +
  geom_point()

data_check <- "new"

co_lag <- 1 # fixed freshwater residency for coho: 1 year on average from Wade & Irvine 2018 report from Keogh and Holtby et al. 1990 CJFAS paper from Carnation Creek
co_smolt_lag <- 2
ct_lag <- 2 # fixed freshwater residency for cutthroat: 2-4 Armstrong 1971 TAFS, Trotter 1989 suggests age 2 for estuarine/coastal:  Losee et al. (2018) suggests age 2 for coastal cutties
# Smith 1980 BC FLNRO suggests age 3 for Keogh River
ct_lags <- 2:4
ct_lags_prop <- c(0.32,0.45,0.20)
ct_lags_prop <- ct_lags_prop/sum(ct_lags_prop)
ct_smolt_lags <- 2:4
ct_smolt_lags_prop <- c(1,2,1)
ct_smolt_lags_prop <- ct_smolt_lags_prop/sum(ct_smolt_lags_prop)

dv_lag <- 4 # freshwater residence: Armstrong 1970 FRBC and Dolloff & Reeves 1990 CJFAS suggest ~ age 1-4 for dolly varden smoltification
dv_lags <- 2:4
dv_lags_prop <- c(1, 8, 1) # Smith and Slaney 1980
dv_lags_prop <- dv_lags_prop/sum(dv_lags_prop)
dv_smolt_lags <- 0:4
dv_smolt_lags_prop <- c(8,37,38,11,7)
dv_smolt_lags_prop <- dv_smolt_lags_prop/sum(dv_smolt_lags_prop)
ch_lag <- 4 # fixed 4 year life cycle for chum: Neave et al. 1952
pink_lag <- 2 # fixed 2 year life cycle

environ <- readRDS("environ_covars.rds")
environment <- read.csv("Data/keogh environmental covariates.csv",header=TRUE)
environ$seal_abundance[environ$year>=2014 & is.na(environ$seal_abundance)] <- environ$seal_abundance[environ$year==2014]
environ$seal_density[environ$year>=2014 & is.na(environ$seal_density)] <- environ$seal_density[environ$year==2014]
environ$total[environ$year>=2015 & is.na(environ$total)] <- environ$total[environ$year==2015]
forestry <- readRDS("Data/keogh_logging.rds")
forestry <- subset(forestry,forestry$Year>=min(environ$year))
#saveRDS(environ,"environ_covars.rds")
#environ <- readRDS("environ_covars.rds")

keogh <- read.csv("Data/Keogh_Database_Final_01oct18.csv",stringsAsFactors = F,header=T)
keogh[keogh$scale=="y " | keogh$scale=="Y","scale"] <- "y"
keogh_atlas <- read.csv("Data/Keogh sh smolts Atlas 2015.csv",stringsAsFactors=F,header=T)
keogh_instream <- read.csv("Data/Keogh smolts instream outmigration.csv",stringsAsFactors=F,header=T)
keogh_adults <- read.csv("Data/Keogh sh adults.csv",stringsAsFactors=F,header=T)

pinks <- read.csv("Data/Keogh_Pink_Salm_NuSeds.csv",stringsAsFactors = F,header=T)
pinks <- pinks[order(pinks$Year),]
pinks$Stock[pinks$Year<=1997] <- 2*pinks$Stock[pinks$Year<=1997] # from Bailey et al. 2018 - Pink salmon prior to 1997 sampled by stream walks and need to be doubled to correct for abundance patterns. After 1997 pink salmon counted by resistivitycounter.

clux_pinks <- read.csv("Data/Clux_Pink_Salm_NuSeds.csv",stringsAsFactors = F,header=T)
clux_pinks <- clux_pinks[order(clux_pinks$Year),]

pink_data <- merge(pinks,clux_pinks,by="Year",all=T)
pink_data <- pink_data[complete.cases(pink_data),]
missYears <- lm(log(pink_data$Stock)~log(pink_data$Clux_Stock))
summary(missYears)
pinks$Stock[is.na(pinks$Stock)] <- coef(missYears)[1]+coef(missYears)[2]*clux_pinks$Clux_Stock[is.na(pinks$Stock)]
pinks <- round(pinks,0)

coho <- read.csv("Data/Keogh coho adults.csv",stringsAsFactors = F,header=T)
coho <- coho[order(coho$Year),]
#coho$Adults[coho$Year<=1997] <- NA#2*coho$Adults[coho$Year<=1997] # from Bailey et al. 2018 - Pink salmon prior to 1997 sampled by stream walks and need to be doubled to correct for abundance patterns. After 1997 pink salmon counted by resistivitycounter.

chum <- read.csv("Data/Keogh chum adults.csv",stringsAsFactors = F,header=T)
chum <- chum[order(chum$Year),]
chum$Adults[chum$Adults==0] <- 0
#chum$Adults[chum$Year<=1997] <- 2*chum$Adults[chum$Year<=1997] # from Bailey et al. 2018 - Pink salmon prior to 1997 sampled by stream walks and need to be doubled to correct for abundance patterns. After 1997 pink salmon counted by resistivitycounter.

# read in coho smolt data from Tom Johnston dataset
coSm <- read.csv("Data/Keogh coho smolts.csv",stringsAsFactors = F,header=T)
coSm <- coSm[order(coSm$Year),]

keogh <- keogh[keogh$hatchery!=1|is.na(keogh$hatchery),]

keogh$life_stage[keogh$life_stage=="F"] <- "f"
keogh$life_stage[keogh$life_stage==" "] <- ""

keogh$total_age <- sapply(keogh$age_final,function(x){as.numeric(strsplit(x,split="\\.")[[1]][1])+sum(as.numeric(strsplit(gsub("[^[:digit:]]","",strsplit(x,split="\\.")[[1]][2]),split="")[[1]]))+nchar(gsub("[[:digit:]]","",strsplit(x,split="\\.")[[1]][2]))})

keogh$age2 <- as.numeric(keogh$fresh_age)
keogh$age_final <- ifelse(!is.na(keogh$total_age),keogh$age_final,keogh$X1_ager)
keogh$age <- sapply(keogh$age_final,function(x){as.numeric(strsplit(x,split="\\.")[[1]][1])})
keogh$OceanAge <- sapply(keogh$age_final,function(x){sum(as.numeric(strsplit(gsub("[^[:digit:]]","",strsplit(x,split="\\.")[[1]][2]),split="")[[1]]))+nchar(gsub("[[:digit:]]","",strsplit(x,split="\\.")[[1]][2]))})

keogh$AdultStrat <- sapply(keogh$age_final,function(x){strsplit(x,split="\\.")[[1]][2]})

# set ocean ages for all fish, smolts are given age 0 in the ocean

keogh$age_ocean <- as.numeric(keogh$OceanAge)
keogh$age_ocean[keogh$life_stage=="s"] <- 0

# assign brood year to fish as the total freshwater and ocean ages

keogh$hatch_year <- ifelse(keogh$life_stage=="s",keogh$year - keogh$age, keogh$year-(keogh$total_age)) # if adult, brood year is current year minus total age
keogh$smolt_year <- keogh$year-(keogh$age_ocean)

steelhead_data <- keogh[keogh$species=="sh",]
steelhead_data$sampAge <- steelhead_data$age+steelhead_data$age_ocean
individuals <- steelhead_data[complete.cases(steelhead_data$sampAge,steelhead_data$fork_length),]

adult_life <- individuals[individuals$life_stage=="a"|individuals$life_stage=="k",]
fresh_ocean_ages <- as.factor(paste(adult_life$age,adult_life$age_ocean,sep="|"))

annual_age <- aggregate(number~smolt_year+hatch_year+species+life_stage+age,data=keogh,FUN=sum,na.rm=T)
colnames(annual_age) <- c("Year","Hatch","Species","Stage","Age","Abundance")
annual_age <- subset(annual_age,Stage%in%c("s","a","k"))
annual_cohorts <- aggregate(Abundance~Hatch+Year+Species+Age,data=annual_age,FUN=sum,na.rm=T)

steel_age <- subset(annual_age,Species=="sh"&Stage=="s")

annual <- aggregate(number~year+species+life_stage,data=keogh,FUN=sum,na.rm=T)

colnames(annual) <- c("Year","Species","Stage","Abundance")
annual$Stage <- factor(annual$Stage,levels=c("f","p","s","j","a","k","","u"))
annual$Stage_Num <- as.numeric(annual$Stage)

size <- aggregate(fork_length~smolt_year+species+life_stage,data=keogh,FUN=mean,na.rm=T)
variance <- aggregate(fork_length~smolt_year+species+life_stage,data=keogh,FUN=sd,na.rm=T)
size$sigma <- variance$fork_length
colnames(size) <- c("Year","Species","Stage","Length","Sigma")

annual_cohort <- dcast(annual_cohorts,Year~Species+Hatch,value.var="Abundance",fun.aggregate=sum)
annual_cohort_sh <- annual_cohort[,grep("sh",colnames(annual_cohort))]
annual_cohort_prop <- annual_cohort_sh/rowSums(annual_cohort_sh,na.rm=T)
annual_prop_sh <- annual_cohort_prop[,grep("sh",colnames(annual_cohort_prop))]
row.names(annual_prop_sh) <- row.names(annual_cohort_prop) <- annual_cohort$Year
colnames(annual_prop_sh) <- min(annual_cohorts$Hatch):max(annual_cohorts$Hatch)
annual_abund <- dcast(annual,Year~Species+Stage,value.var="Abundance")
annual_size <- dcast(size,Year~Species+Stage,value.var="Length")
annual_sigma <- dcast(size,Year~Species+Stage,value.var="Sigma")

# adults with known ocean ages have a known smolt year: they should count as 'abundant' for their smolt year
annual_outmigAdults <- subset(annual_age,Stage%in%c("a","k") & Species%in%c("sh"))
annual_outmigAdults <- aggregate(Abundance~Hatch+Year+Species+Age,data=annual_outmigAdults,FUN=sum,na.rm=T)
annual_outmigAdults <- dcast(annual_outmigAdults,Year~Species+Hatch,value.var="Abundance",fun.aggregate=sum)
yrs <- annual_outmigAdults$Year
annual_outmigAdults <- rowSums(annual_outmigAdults[,grep("sh",colnames(annual_outmigAdults))])
names(annual_outmigAdults) <- yrs

# we want to find how many of the smolts sampled in Year X were from the spawners at Year Y
# to do that, we find the proportions aged and multiply the current year smolts by that proportion from hatch year Y
annual_sh <- subset(annual,Species=="sh")
annual_sh_abund <- dcast(annual_sh,Year~Species+Stage,value.var="Abundance")

annual_shs_abund <- merge(keogh_instream,annual_sh_abund[,c("Year","sh_s")],by="Year",all=TRUE)
annual_shs_abund$Smolts <- ifelse(!is.na(annual_shs_abund$Smolts),annual_shs_abund$Smolts,annual_shs_abund$sh_s)
keogh_smolt <- merge(data.frame("Year"=as.numeric(names(annual_outmigAdults)),"Smolts"=annual_outmigAdults),annual_shs_abund[,c("Year","Smolts")],by="Year",all=T)
keogh_smolt <- data.frame(keogh_smolt$Year,rowSums(keogh_smolt[,-1],na.rm=T))
colnames(keogh_smolt) <- c("Year","Smolts")

mn_smolt_age <- table(keogh$age[keogh$species=="sh" & (keogh$life_stage=="s" | keogh$life_stage=="k" | keogh$life_stage=="a")])
mn_smolt_age <- mn_smolt_age/sum(mn_smolt_age)

years <- min(keogh$hatch_year,na.rm=TRUE):max(keogh$year,na.rm=TRUE)
annual_smolts_sh <- rep(NA,length(years))
names(annual_smolts_sh) <- years
prop_smolts_year <- matrix(NA,nrow=length(years),ncol=length(mn_smolt_age),dimnames=list("Year"=years,"Ages"=1:length(mn_smolt_age)))
for(i in 1:length(annual_smolts_sh))
{
  IIyear <- as.numeric(names(annual_smolts_sh)[i])
  smolt_age <- table(keogh$age[keogh$species=="sh" & (keogh$life_stage=="s" | keogh$life_stage=="k" | keogh$life_stage=="a") & keogh$year==IIyear])
  if(sum(smolt_age,na.rm=TRUE)>0)
  {
    smolt_age <- smolt_age/sum(smolt_age)
  }else{
    smolt_age <- mn_smolt_age
  }
  prop_smolts_year[i,1:length(smolt_age)] <- smolt_age
  annual_smolts_sh[i] <- round(sum(unlist(sapply(1:length(smolt_age),function(x){keogh_smolt$Smolts[match(IIyear+as.numeric(names(smolt_age)[x]),keogh_smolt$Year,nomatch = 0)]*smolt_age[x]})),na.rm=TRUE),0)
}


sh_smolts <- annual_smolts_sh
stock_rec <- data.frame("Year"=as.numeric(names(sh_smolts)),"Smolts"=sh_smolts)
stock_rec$Smolts[stock_rec$Smolts==0] <- NA
stock_rec <- merge(keogh_adults,stock_rec,by="Year",all=T)

keogh_smolts <- merge(keogh_atlas[,c("Year","Total")],stock_rec,by="Year",all=T)
plot(keogh_smolts$Year,keogh_smolts$Total,type="l",xlab="Year",ylab="Smolt abundance",ylim=c(0,max(keogh_smolts,na.rm=T)))
points(keogh_smolts$Year,keogh_smolts$Smolts,pch=21,bg="dodgerblue",type="b",lty=2,lwd=1)
#points(keogh_sh_v2$Year,keogh_sh_v2$Smolts,pch=21,bg="orange",ylim=c(0,15000),type="b",lty=1,lwd=1,col="orange")
legend("topright",c("Tom Johnston","Current QA/QC"),col=c("black","dodgerblue"),pch=c(NA,21),pt.bg=c(NA,"dodgerblue"),lty=c(1,2),lwd=c(2,2),bty="n")

# find out the smolts that produced the steelhead adults
# we want to find how many of the adults sampled in Year X were from the smolts at Year Y
# to do that, we find the proportions aged and multiply the current year adults by that proportion from smolting year Y

mn_smolt_age <- table(keogh$age_ocean[keogh$species=="sh" & (keogh$life_stage=="a" | keogh$life_stage=="k")])
mn_smolt_age <- mn_smolt_age/sum(mn_smolt_age)

annual_smolts_sh <- rep(NA,nrow(stock_rec))
names(annual_smolts_sh) <- stock_rec$Year
prop_adults_year <- matrix(NA,nrow=nrow(stock_rec),ncol=length(mn_smolt_age),dimnames=list("Year"=stock_rec$Year,"Ages"=1:length(mn_smolt_age)))
for(i in 1:nrow(stock_rec))
{
  IIyear <- as.numeric(names(annual_smolts_sh)[i])
  smolt_age <- table(keogh$age_ocean[keogh$species=="sh" & (keogh$life_stage=="a" | keogh$life_stage=="k") & keogh$year==IIyear])
  if(sum(smolt_age,na.rm=TRUE)>0)
  {
    smolt_age <- smolt_age/sum(smolt_age)
  }else{
    smolt_age <- mn_smolt_age
  }
  prop_adults_year[i,1:length(smolt_age)] <- smolt_age
  annual_smolts_sh[i] <- round(sum(unlist(sapply(1:length(smolt_age),function(x){stock_rec$Smolts[match(IIyear-as.numeric(names(smolt_age)[x]),stock_rec$Year,nomatch = 0)]*smolt_age[x]})),na.rm=TRUE),0)
}

sh_smolting <- annual_smolts_sh

smolting_df <- data.frame("Year"=as.numeric(names(sh_smolting)),"juv_Cohort"=sh_smolting)
smolting_df$juv_Cohort[smolting_df$juv_Cohort==0] <- NA

stock_rec <- merge(stock_rec,smolting_df,by="Year",all=T)
plot(log(Smolts/Adults)~Adults,data=stock_rec,pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))

keogh_smolts <- merge(keogh_atlas[,c("Year","Total")],stock_rec,by="Year",all=T)

# compile cutthroat trout data

ct_R <- aggregate(number~year+life_stage,data=keogh[keogh$species=="ct",],FUN=sum,na.rm=T)
ct_age <- aggregate(number~year+fresh_age+life_stage,data=keogh[keogh$species=="ct",],FUN=sum,na.rm=T)
CT_annual <- dcast(ct_R,year~life_stage,value.var="number")
colnames(CT_annual)[colnames(CT_annual)=="year"] <- "Year"
CT_annual <- merge(stock_rec,CT_annual,by="Year",all=TRUE)
CT_annual <- CT_annual[,-match(c("Adults","Smolts"),colnames(CT_annual))]

#ct_SR <- data.frame("Year"=CT_annual$Year[1:(length(CT_annual$Year)-ct_lag)],"Stock"=CT_annual$a[1:(length(CT_annual$Year)-ct_lag)],"Recruits"=CT_annual$s[(ct_lag+1):length(CT_annual$Year)])

ct_SR <- data.frame("Year"=CT_annual$Year,"Stock"=CT_annual$a)
ct_SR$Recruits <- round(rowSums(sapply(1:length(ct_lags),function(x){c(CT_annual$s[(ct_lags[x]+1):length(CT_annual$Year)]*ct_lags_prop[x],rep(NA,ct_lags[x]))}),na.rm=TRUE),0)

ct_SR$Juv_cohort <- round(rowSums(sapply(1:length(ct_smolt_lags),function(x){c(rep(NA,ct_smolt_lags[x]),ct_SR$Recruits[1:(length(CT_annual$Year)-(ct_smolt_lags[x]))]*ct_smolt_lags_prop[x])}),na.rm=TRUE),0)
ct_SR[ct_SR==0] <- NA
plot(log(Recruits/Stock)~Stock,data=ct_SR,pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))

# dolly vardens
dv_R <- aggregate(number~year+life_stage,data=keogh[keogh$species=="dv",],FUN=sum,na.rm=T)
dv_age <- aggregate(number~year+fresh_age+life_stage,data=keogh[keogh$species=="dv",],FUN=sum,na.rm=T)
DV_annual <- dcast(dv_R,year~life_stage,value.var="number")
colnames(DV_annual)[colnames(DV_annual)=="year"] <- "Year"
DV_annual <- merge(stock_rec,DV_annual,by="Year",all=TRUE)
DV_annual <- DV_annual[,-match(c("Adults","Smolts"),colnames(DV_annual))]

dv_SR <- data.frame("Year"=DV_annual$Year,"Stock"=DV_annual$a)
dv_SR$Recruits <- round(rowSums(sapply(1:length(dv_lags),function(x){c(DV_annual$s[(dv_lags[x]+1):length(DV_annual$Year)]*dv_lags_prop[x],rep(NA,dv_lags[x]))}),na.rm=TRUE),0)
dv_SR$Juv_cohort <- round(rowSums(sapply(1:length(dv_smolt_lags),function(x){c(rep(NA,dv_smolt_lags[x]),dv_SR$Recruits[1:(length(DV_annual$Year)-(dv_smolt_lags[x]))]*dv_smolt_lags_prop[x])}),na.rm=TRUE),0)
dv_SR[dv_SR==0] <- NA
plot(log(Recruits/Stock)~Stock,data=dv_SR,pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))

# compile data on coho salmon
co_R <- aggregate(number~year+life_stage,data=keogh[keogh$species=="co",],FUN=sum,na.rm=T)
co_age <- aggregate(number~year+fresh_age+life_stage,data=keogh[keogh$species=="co",],FUN=sum,na.rm=T)
co_annual <- dcast(co_R,year~life_stage,value.var="number")

co_S <- rep(NA,length(coho$Year))
names(co_S) <- coho$Year
co_S <- coho$Adults

co_R <- rep(NA,length(coho$Year))
names(co_R) <- coho$Year
co_R[match(coSm$Year,coho$Year,nomatch=0)] <- coSm$Smolts[match(coho$Year,coSm$Year,nomatch=0)]

co_SR <- data.frame("Year"=coho$Year[1:(length(coho$Year)-co_lag)],"Stock"=co_S[1:(length(coho$Year)-co_lag)],"Recruits"=co_R[(co_lag+1):length(coho$Year)])
co_SR$Juv_cohort <- c(rep(NA,co_smolt_lag),co_SR$Recruits[1:(length(co_SR$Year)-(co_smolt_lag))])

# to include or not include pre-1997 coho
co_SR$Stock[co_SR$Year<=1997] <- NA
plot(log(Recruits/Stock)~Stock,data=co_SR,pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))

# compile data on chum salmon
ch_SR <- data.frame("Year"=chum$Year,"Stock"=chum$Adults)
ch_SR$Recruits <- c(chum$Adult[(ch_lag+1):length(chum$Year)],rep(NA,ch_lag))
ch_SR$Juv_cohort <- c(rep(NA,ch_lag),ch_SR$Stock[1:(length(ch_SR$Year)-(ch_lag))])
plot(log(Recruits/Stock)~Stock,data=ch_SR,pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))

# read in and calculate pink salmon stock-recruitment
pinks <- data.frame("Year"=pinks$Year,"Stock"=pinks$Stock)
pinks$Recruits <- c(pinks$Stock[(pink_lag+1):length(pinks$Year)],rep(NA,pink_lag))
pinks$Juv_cohort <- c(rep(NA,pink_lag),pinks$Stock[1:(length(pinks$Year)-(pink_lag))])
plot(log(Recruits/Stock)~Stock,data=pinks,pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))

# compile information for all salmon species
colnames(stock_rec) <- c("Year","sh_Adults","sh_Smolts","sh_juv_Cohort")
colnames(pinks) <- c("Year","pk_Adults","pk_Recruits","pk_juv_Cohort")
colnames(ct_SR) <- c("Year","ct_Adults","ct_Smolts","ct_juv_Cohort")
colnames(co_SR) <- c("Year","co_Adults","co_Smolts","co_juv_Cohort")
colnames(dv_SR) <- c("Year","dv_Adults","dv_Smolts","dv_juv_Cohort")
colnames(ch_SR) <- c("Year","ch_Adults","ch_Recruits","ch_juv_Cohort")

keogh_StockRec <- round(merge(stock_rec,merge(pinks,merge(ct_SR,merge(co_SR,merge(dv_SR,ch_SR,by="Year",all=T),by="Year",all=T),by="Year",all=T),by="Year",all=T),by="Year",all=T),0)

colnames(environ)[colnames(environ)=="year"] <- "Year"
keogh_SR <- merge(keogh_StockRec,environ,by="Year")
keogh_SR <- merge(keogh_SR,forestry,by="Year")

example <- reshape(keogh_SR,direction = "long",varying = list(c("sh_Adults","dv_Adults","ct_Adults","pk_Adults","ch_Adults","co_Adults"),c("sh_Smolts","dv_Smolts","ct_Smolts","pk_Recruits","ch_Recruits","co_Smolts"),c("sh_juv_Cohort","dv_juv_Cohort","ct_juv_Cohort","pk_juv_Cohort","ch_juv_Cohort","co_juv_Cohort")),v.names=c("Stock","Recruits","juvCohort"),idvar="Species")

example$Species <- rep(c("Steelhead","Dolly Varden","Cutthroat","Pink","Chum","Coho"),each=nrow(keogh_SR))
example$Species <- factor(example$Species,levels=c("Steelhead","Dolly Varden","Cutthroat","Pink","Chum","Coho"))

keogh_long <- example

seals <- rep(NA,nrow(keogh_long))
names(seals) <- keogh_long$Year
cumul_logging <- mean_logging <- sumTemp <- sumRain <- winTemp <- winRain <- freshCoho <- freshSteel <- freshCutt <- freshDolly <- freshPink <- npgo <- oceanSalmon <- mei <- seals

for(i in 1:nrow(keogh_long))
{
  year_matches <- NA
  age_matches <- NA
  pink_matches <- NA
  spp <- keogh_long$Species[i]
  if(spp=="Steelhead"){
    # lag freshwater conditions forwards from hatch year
    if(any(as.numeric(row.names(prop_smolts_year)) %in% keogh_long$Year[i])){
      year_matches <- match(keogh_long$Year[i] + as.integer(colnames(prop_smolts_year)),keogh_long$Year,nomatch=0)
      year_0 <- year_matches!=0
      year_matches <- sort(year_matches[year_matches!=0])
      age_matches <- match(as.numeric(row.names(prop_smolts_year)),keogh_long$Year[i],nomatch=0)
      props <- prop_smolts_year[age_matches,year_0]/sum(prop_smolts_year[age_matches,year_0])
      pink_matches <- match(keogh_long$Year[year_matches],keogh_long$Year[keogh_long$Species=="Pink"],nomatch=0)
    }else{
      year_matches <- match(keogh_long$Year[i] + as.integer(colnames(prop_smolts_year)),keogh_long$Year,nomatch=0)
      year_0 <- year_matches!=0
      year_matches <- sort(year_matches[year_matches!=0])
      age_matches <- match(as.numeric(row.names(prop_smolts_year)),keogh_long$Year[i],nomatch=0)
      props <- prop_smolts_year[1,year_0]/sum(prop_smolts_year[1,year_0])
      pink_matches <- match(keogh_long$Year[year_matches],keogh_long$Year[keogh_long$Species=="Pink"],nomatch=0)
    }
    lags <- 1:length(props)
    mean_temp <- sapply(lags,function(x){mean(keogh_long$mean_temp[year_matches[1:x]],na.rm=TRUE)})
    total_rain <- sapply(lags,function(x){mean(keogh_long$total_rain[year_matches[1:x]],na.rm=TRUE)})
    win_mean_temp <- sapply(lags,function(x){mean(keogh_long$win_mean_temp[year_matches[1:x]],na.rm=TRUE)})
    win_rain <- sapply(lags,function(x){mean(keogh_long$win_rain[year_matches[1:x]],na.rm=TRUE)})
    cum_log <- sapply(lags,function(x){mean(keogh_long$cumul_footprint[year_matches[1:x]],na.rm=TRUE)})
    mn_log <- sapply(lags,function(x){mean(keogh_long$cumul_log[year_matches[1:x]],na.rm=TRUE)})
    
    Coho <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Coho"][pink_matches[1:x]],na.rm=TRUE)})
    Steel <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Steelhead"][pink_matches[1:x]],na.rm=TRUE)})
    Cutty <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Cutthroat"][pink_matches[1:x]],na.rm=TRUE)})
    Dolly <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Dolly Varden"][pink_matches[1:x]],na.rm=TRUE)})
    Pink <- sapply(lags,function(x){mean(keogh_long$Stock[keogh_long$Species=="Pink"][pink_matches[1:x]],na.rm=TRUE)})
    
    sumTemp[i] <- sum(props * mean_temp,na.rm=TRUE)
    sumRain[i]<- sum(props * total_rain,na.rm=TRUE)
    winTemp[i]<- sum(props * win_mean_temp,na.rm=TRUE)
    winRain[i]<- sum(props * win_rain,na.rm=TRUE)
    freshCoho[i] <- sum(props * Coho,na.rm=TRUE)
    freshSteel[i] <- sum(props * Steel,na.rm=TRUE)
    freshCutt[i] <- sum(props * Cutty,na.rm=TRUE)
    freshDolly[i] <- sum(props * Dolly,na.rm=TRUE)
    freshPink[i]<- sum(props * Pink,na.rm=TRUE)
    cumul_logging[i] <- sum(props * cum_log,na.rm=TRUE)
    mean_logging[i] <- sum(props * mn_log,na.rm=TRUE)
    seals[i] <- sum(props * keogh_long$seal_density[i],na.rm=TRUE)
    
    # lag ocean condition backwards from spawning year
    if(any(as.numeric(row.names(prop_adults_year)) %in% keogh_long$Year[i])){
      year_matches <- match(keogh_long$Year[i] - as.integer(colnames(prop_adults_year)),keogh_long$Year,nomatch=0)
      year_0 <- year_matches!=0
      year_matches <- sort(year_matches[year_matches!=0])
      age_matches <- match(as.numeric(row.names(prop_smolts_year)),keogh_long$Year[i],nomatch=0)
      props <- prop_adults_year[age_matches,year_0]/sum(prop_adults_year[age_matches,year_0])
      pink_matches <- match(keogh_long$Year[year_matches],keogh_long$Year[keogh_long$Species=="Pink"],nomatch=0)
    }else{
      year_matches <- match(keogh_long$Year[i] - as.integer(colnames(prop_adults_year)),keogh_long$Year,nomatch=0)
      year_0 <- year_matches!=0
      year_matches <- sort(year_matches[year_matches!=0])
      age_matches <- match(as.numeric(row.names(prop_smolts_year)),keogh_long$Year[i],nomatch=0)
      props <- prop_adults_year[1,year_0]/sum(prop_adults_year[1,year_0])
      pink_matches <- match(keogh_long$Year[year_matches],keogh_long$Year[keogh_long$Species=="Pink"],nomatch=0)
    }
    lags <- 1:length(props)
    npgo_lag <- sapply(lags,function(x){mean(keogh_long$npgo[year_matches[1:x]],na.rm=TRUE)})
    oceanSalmon_lag <- sapply(lags,function(x){mean(keogh_long$total[year_matches[1:x]],na.rm=TRUE)})
    mei_lag <- sapply(lags,function(x){mean(keogh_long$mei[year_matches[1:x]],na.rm=TRUE)})
    
    npgo[i] <- sum(props * npgo_lag,na.rm=TRUE)
    oceanSalmon[i] <- sum(props * oceanSalmon_lag,na.rm=TRUE)
    mei[i] <- sum(props * mei_lag,na.rm=TRUE)
    seals[i] <- seals[i] + sum(props * keogh_long$seal_density[year_matches],na.rm=TRUE)
    
  }
  if(spp=="Pink"){
    juvLag <- 0
    adultLag <- pink_lag
    
    # lag freshwater conditions forwards from hatch year
    year_matches <- match(keogh_long$Year[i] + juvLag,keogh_long$Year,nomatch=0)
    year_0 <- year_matches!=0
    year_matches <- sort(year_matches[year_matches!=0])
    pink_matches <- match(keogh_long$Year[year_matches],keogh_long$Year[keogh_long$Species=="Pink"],nomatch=0)
    props <- 1
    lags <- length(props)
    mean_temp <- sapply(lags,function(x){mean(keogh_long$mean_temp[year_matches[1:x]],na.rm=TRUE)})
    total_rain <- sapply(lags,function(x){mean(keogh_long$total_rain[year_matches[1:x]],na.rm=TRUE)})
    win_mean_temp <- sapply(lags,function(x){mean(keogh_long$win_mean_temp[year_matches[1:x]],na.rm=TRUE)})
    win_rain <- sapply(lags,function(x){mean(keogh_long$win_rain[year_matches[1:x]],na.rm=TRUE)})
    cum_log <- sapply(lags,function(x){mean(keogh_long$cumul_footprint[year_matches[1:x]],na.rm=TRUE)})
    mn_log <- sapply(lags,function(x){mean(keogh_long$cumul_log[year_matches[1:x]],na.rm=TRUE)})
    Coho <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Coho"][pink_matches[1:x]],na.rm=TRUE)})
    Steel <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Steelhead"][pink_matches[1:x]],na.rm=TRUE)})
    Cutty <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Cutthroat"][pink_matches[1:x]],na.rm=TRUE)})
    Dolly <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Dolly Varden"][pink_matches[1:x]],na.rm=TRUE)})
    Pink <- sapply(lags,function(x){mean(keogh_long$Stock[keogh_long$Species=="Pink"][pink_matches[1:x]],na.rm=TRUE)})
    
    sumTemp[i] <- sum(props * mean_temp,na.rm=TRUE)
    sumRain[i]<- sum(props * total_rain,na.rm=TRUE)
    winTemp[i]<- sum(props * win_mean_temp,na.rm=TRUE)
    winRain[i]<- sum(props * win_rain,na.rm=TRUE)
    freshCoho[i] <- sum(props * Coho,na.rm=TRUE)
    freshSteel[i] <- sum(props * Steel,na.rm=TRUE)
    freshCutt[i] <- sum(props * Cutty,na.rm=TRUE)
    freshDolly[i] <- sum(props * Dolly,na.rm=TRUE)
    freshPink[i]<- sum(props * Pink,na.rm=TRUE)
    cumul_logging[i] <- sum(props * cum_log,na.rm=TRUE)
    mean_logging[i] <- sum(props * mn_log,na.rm=TRUE)

    # lag ocean condition backwards from spawning year
    year_matches <- match(keogh_long$Year[i] - adultLag,keogh_long$Year,nomatch=0)
    year_0 <- year_matches!=0
    year_matches <- sort(year_matches[year_matches!=0])
    props <- 1
    lags <- 1:length(props)
    npgo_lag <- sapply(lags,function(x){mean(keogh_long$npgo[year_matches[1:x]],na.rm=TRUE)})
    oceanSalmon_lag <- sapply(lags,function(x){mean(keogh_long$total[year_matches[1:x]],na.rm=TRUE)})
    mei_lag <- sapply(lags,function(x){mean(keogh_long$mei[year_matches[1:x]],na.rm=TRUE)})
    npgo[i] <- sum(props * npgo_lag,na.rm=TRUE)
    oceanSalmon[i] <- sum(props * oceanSalmon_lag,na.rm=TRUE)
    mei[i] <- sum(props * mei_lag,na.rm=TRUE)
    seals[i] <- sum(props * keogh_long$seal_density[i],na.rm=TRUE)
  }
  if(spp=="Dolly Varden"){
    
    # lag freshwater conditions forwards from hatch year
    year_matches <- match(keogh_long$Year[i] + dv_lags,keogh_long$Year,nomatch=0)
    year_0 <- year_matches!=0
    year_matches <- sort(year_matches[year_matches!=0])
    props <- dv_lags_prop[year_0]/sum(dv_lags_prop[year_0])
    pink_matches <- match(keogh_long$Year[year_matches],keogh_long$Year[keogh_long$Species=="Pink"],nomatch=0)
    
    lags <- 1:length(props)
    mean_temp <- sapply(lags,function(x){mean(keogh_long$mean_temp[year_matches[1:x]],na.rm=TRUE)})
    total_rain <- sapply(lags,function(x){mean(keogh_long$total_rain[year_matches[1:x]],na.rm=TRUE)})
    win_mean_temp <- sapply(lags,function(x){mean(keogh_long$win_mean_temp[year_matches[1:x]],na.rm=TRUE)})
    win_rain <- sapply(lags,function(x){mean(keogh_long$win_rain[year_matches[1:x]],na.rm=TRUE)})
    cum_log <- sapply(lags,function(x){mean(keogh_long$cumul_footprint[year_matches[1:x]],na.rm=TRUE)})
    mn_log <- sapply(lags,function(x){mean(keogh_long$cumul_log[year_matches[1:x]],na.rm=TRUE)})
    Coho <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Coho"][pink_matches[1:x]],na.rm=TRUE)})
    Steel <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Steelhead"][pink_matches[1:x]],na.rm=TRUE)})
    Cutty <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Cutthroat"][pink_matches[1:x]],na.rm=TRUE)})
    Dolly <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Dolly Varden"][pink_matches[1:x]],na.rm=TRUE)})
    Pink <- sapply(lags,function(x){mean(keogh_long$Stock[keogh_long$Species=="Pink"][pink_matches[1:x]],na.rm=TRUE)})
    sumTemp[i] <- sum(props * mean_temp,na.rm=TRUE)
    sumRain[i]<- sum(props * total_rain,na.rm=TRUE)
    winTemp[i]<- sum(props * win_mean_temp,na.rm=TRUE)
    winRain[i]<- sum(props * win_rain,na.rm=TRUE)
    freshCoho[i] <- sum(props * Coho,na.rm=TRUE)
    freshSteel[i] <- sum(props * Steel,na.rm=TRUE)
    freshCutt[i] <- sum(props * Cutty,na.rm=TRUE)
    freshDolly[i] <- sum(props * Dolly,na.rm=TRUE)
    freshPink[i]<- sum(props * Pink,na.rm=TRUE)
    cumul_logging[i] <- sum(props * cum_log,na.rm=TRUE)
    mean_logging[i] <- sum(props * mn_log,na.rm=TRUE)
    seals[i] <- sum(props * keogh_long$seal_density[i],na.rm=TRUE)
    
    # lag ocean condition backwards from spawning year
    year_matches <- match(keogh_long$Year[i] - dv_smolt_lags,keogh_long$Year,nomatch=0)
    year_0 <- year_matches!=0
    year_matches <- sort(year_matches[year_matches!=0])
    props <- dv_smolt_lags_prop[year_0]/sum(dv_smolt_lags_prop[year_0])
    lags <- 1:length(props)
    npgo_lag <- sapply(lags,function(x){mean(keogh_long$npgo[year_matches[1:x]],na.rm=TRUE)})
    oceanSalmon_lag <- sapply(lags,function(x){mean(keogh_long$total[year_matches[1:x]],na.rm=TRUE)})
    mei_lag <- sapply(lags,function(x){mean(keogh_long$mei[year_matches[1:x]],na.rm=TRUE)})
    npgo[i] <- sum(props * npgo_lag,na.rm=TRUE)
    oceanSalmon[i] <- sum(props * oceanSalmon_lag,na.rm=TRUE)
    mei[i] <- sum(props * mei_lag,na.rm=TRUE)
    seals[i] <- seals[i] + sum(props * keogh_long$seal_density[year_matches],na.rm=TRUE)
  }
  if(spp=="Cutthroat"){
    # lag freshwater conditions forwards from hatch year
    year_matches <- match(keogh_long$Year[i] + ct_lags,keogh_long$Year,nomatch=0)
    year_0 <- year_matches!=0
    year_matches <- sort(year_matches[year_matches!=0])
    props <- ct_lags_prop[year_0]/sum(ct_lags_prop[year_0])
    pink_matches <- match(keogh_long$Year[year_matches],keogh_long$Year[keogh_long$Species=="Pink"],nomatch=0)
    lags <- 1:length(props)
    mean_temp <- sapply(lags,function(x){mean(keogh_long$mean_temp[year_matches[1:x]],na.rm=TRUE)})
    total_rain <- sapply(lags,function(x){mean(keogh_long$total_rain[year_matches[1:x]],na.rm=TRUE)})
    win_mean_temp <- sapply(lags,function(x){mean(keogh_long$win_mean_temp[year_matches[1:x]],na.rm=TRUE)})
    win_rain <- sapply(lags,function(x){mean(keogh_long$win_rain[year_matches[1:x]],na.rm=TRUE)})
    cum_log <- sapply(lags,function(x){mean(keogh_long$cumul_footprint[year_matches[1:x]],na.rm=TRUE)})
    mn_log <- sapply(lags,function(x){mean(keogh_long$cumul_log[year_matches[1:x]],na.rm=TRUE)})
    Coho <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Coho"][pink_matches[1:x]],na.rm=TRUE)})
    Steel <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Steelhead"][pink_matches[1:x]],na.rm=TRUE)})
    Cutty <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Cutthroat"][pink_matches[1:x]],na.rm=TRUE)})
    Dolly <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Dolly Varden"][pink_matches[1:x]],na.rm=TRUE)})
    Pink <- sapply(lags,function(x){mean(keogh_long$Stock[keogh_long$Species=="Pink"][pink_matches[1:x]],na.rm=TRUE)})
    
    sumTemp[i] <- sum(props * mean_temp,na.rm=TRUE)
    sumRain[i]<- sum(props * total_rain,na.rm=TRUE)
    winTemp[i]<- sum(props * win_mean_temp,na.rm=TRUE)
    winRain[i]<- sum(props * win_rain,na.rm=TRUE)
    freshCoho[i] <- sum(props * Coho,na.rm=TRUE)
    freshSteel[i] <- sum(props * Steel,na.rm=TRUE)
    freshCutt[i] <- sum(props * Cutty,na.rm=TRUE)
    freshDolly[i] <- sum(props * Dolly,na.rm=TRUE)
    freshPink[i]<- sum(props * Pink,na.rm=TRUE)
    cumul_logging[i] <- sum(props * cum_log,na.rm=TRUE)
    mean_logging[i] <- sum(props * mn_log,na.rm=TRUE)
    seals[i] <- sum(props * keogh_long$seal_density[i],na.rm=TRUE)
    
    # lag ocean condition backwards from spawning year
    year_matches <- match(keogh_long$Year[i] - ct_smolt_lags,keogh_long$Year,nomatch=0)
    year_0 <- year_matches!=0
    year_matches <- sort(year_matches[year_matches!=0])
    props <- ct_smolt_lags_prop[year_0]/sum(ct_smolt_lags_prop[year_0])
    lags <- 1:length(props)
    npgo_lag <- sapply(lags,function(x){mean(keogh_long$npgo[year_matches[1:x]],na.rm=TRUE)})
    oceanSalmon_lag <- sapply(lags,function(x){mean(keogh_long$total[year_matches[1:x]],na.rm=TRUE)})
    mei_lag <- sapply(lags,function(x){mean(keogh_long$mei[year_matches[1:x]],na.rm=TRUE)})
    
    npgo[i] <- sum(props * npgo_lag,na.rm=TRUE)
    oceanSalmon[i] <- sum(props * oceanSalmon_lag,na.rm=TRUE)
    mei[i] <- sum(props * mei_lag,na.rm=TRUE)
    seals[i] <- seals[i] + sum(props * keogh_long$seal_density[year_matches],na.rm=TRUE)
  }
  if(spp=="Coho"){
    juvLag <- co_lag
    adultLag <- co_smolt_lag
    
    # lag freshwater conditions forwards from hatch year
    year_matches <- match(keogh_long$Year[i] + juvLag,keogh_long$Year,nomatch=0)
    year_0 <- year_matches!=0
    year_matches <- sort(year_matches[year_matches!=0])
    pink_matches <- match(keogh_long$Year[year_matches],keogh_long$Year[keogh_long$Species=="Pink"],nomatch=0)
    props <- 1
    lags <- length(props)
    mean_temp <- sapply(lags,function(x){mean(keogh_long$mean_temp[year_matches[1:x]],na.rm=TRUE)})
    total_rain <- sapply(lags,function(x){mean(keogh_long$total_rain[year_matches[1:x]],na.rm=TRUE)})
    win_mean_temp <- sapply(lags,function(x){mean(keogh_long$win_mean_temp[year_matches[1:x]],na.rm=TRUE)})
    win_rain <- sapply(lags,function(x){mean(keogh_long$win_rain[year_matches[1:x]],na.rm=TRUE)})
    cum_log <- sapply(lags,function(x){mean(keogh_long$cumul_footprint[year_matches[1:x]],na.rm=TRUE)})
    mn_log <- sapply(lags,function(x){mean(keogh_long$cumul_log[year_matches[1:x]],na.rm=TRUE)})
    Coho <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Coho"][pink_matches[1:x]],na.rm=TRUE)})
    Steel <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Steelhead"][pink_matches[1:x]],na.rm=TRUE)})
    Cutty <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Cutthroat"][pink_matches[1:x]],na.rm=TRUE)})
    Dolly <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Dolly Varden"][pink_matches[1:x]],na.rm=TRUE)})
    Pink <- sapply(lags,function(x){mean(keogh_long$Stock[keogh_long$Species=="Pink"][pink_matches[1:x]],na.rm=TRUE)})
    
    sumTemp[i] <- sum(props * mean_temp,na.rm=TRUE)
    sumRain[i]<- sum(props * total_rain,na.rm=TRUE)
    winTemp[i]<- sum(props * win_mean_temp,na.rm=TRUE)
    winRain[i]<- sum(props * win_rain,na.rm=TRUE)
    freshCoho[i] <- sum(props * Coho,na.rm=TRUE)
    freshSteel[i] <- sum(props * Steel,na.rm=TRUE)
    freshCutt[i] <- sum(props * Cutty,na.rm=TRUE)
    freshDolly[i] <- sum(props * Dolly,na.rm=TRUE)
    freshPink[i]<- sum(props * Pink,na.rm=TRUE)
    cumul_logging[i] <- sum(props * cum_log,na.rm=TRUE)
    mean_logging[i] <- sum(props * mn_log,na.rm=TRUE)
    
    # lag ocean condition backwards from spawning year
    year_matches <- match(keogh_long$Year[i] - adultLag,keogh_long$Year,nomatch=0)
    year_0 <- year_matches!=0
    year_matches <- sort(year_matches[year_matches!=0])
    props <- 1
    lags <- 1:length(props)
    npgo_lag <- sapply(lags,function(x){mean(keogh_long$npgo[year_matches[1:x]],na.rm=TRUE)})
    oceanSalmon_lag <- sapply(lags,function(x){mean(keogh_long$total[year_matches[1:x]],na.rm=TRUE)})
    mei_lag <- sapply(lags,function(x){mean(keogh_long$mei[year_matches[1:x]],na.rm=TRUE)})
    npgo[i] <- sum(props * npgo_lag,na.rm=TRUE)
    oceanSalmon[i] <- sum(props * oceanSalmon_lag,na.rm=TRUE)
    mei[i] <- sum(props * mei_lag,na.rm=TRUE)
    seals[i] <- sum(props * keogh_long$seal_density[i],na.rm=TRUE)
  }
  if(spp=="Chum"){
    juvLag <- 0
    adultLag <- ch_lag
    # lag freshwater conditions forwards from hatch year
    year_matches <- match(keogh_long$Year[i] + juvLag,keogh_long$Year,nomatch=0)
    year_0 <- year_matches!=0
    year_matches <- sort(year_matches[year_matches!=0])
    pink_matches <- match(keogh_long$Year[year_matches],keogh_long$Year[keogh_long$Species=="Pink"],nomatch=0)
    props <- 1
    lags <- length(props)
    mean_temp <- sapply(lags,function(x){mean(keogh_long$mean_temp[year_matches[1:x]],na.rm=TRUE)})
    total_rain <- sapply(lags,function(x){mean(keogh_long$total_rain[year_matches[1:x]],na.rm=TRUE)})
    win_mean_temp <- sapply(lags,function(x){mean(keogh_long$win_mean_temp[year_matches[1:x]],na.rm=TRUE)})
    win_rain <- sapply(lags,function(x){mean(keogh_long$win_rain[year_matches[1:x]],na.rm=TRUE)})
    cum_log <- sapply(lags,function(x){mean(keogh_long$cumul_footprint[year_matches[1:x]],na.rm=TRUE)})
    mn_log <- sapply(lags,function(x){mean(keogh_long$cumul_log[year_matches[1:x]],na.rm=TRUE)})
    Coho <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Coho"][pink_matches[1:x]],na.rm=TRUE)})
    Steel <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Steelhead"][pink_matches[1:x]],na.rm=TRUE)})
    Cutty <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Cutthroat"][pink_matches[1:x]],na.rm=TRUE)})
    Dolly <- sapply(lags,function(x){mean(keogh_long$Recruits[keogh_long$Species=="Dolly Varden"][pink_matches[1:x]],na.rm=TRUE)})
    Pink <- sapply(lags,function(x){mean(keogh_long$Stock[keogh_long$Species=="Pink"][pink_matches[1:x]],na.rm=TRUE)})
    
    sumTemp[i] <- sum(props * mean_temp,na.rm=TRUE)
    sumRain[i]<- sum(props * total_rain,na.rm=TRUE)
    winTemp[i]<- sum(props * win_mean_temp,na.rm=TRUE)
    winRain[i]<- sum(props * win_rain,na.rm=TRUE)
    freshCoho[i] <- sum(props * Coho,na.rm=TRUE)
    freshSteel[i] <- sum(props * Steel,na.rm=TRUE)
    freshCutt[i] <- sum(props * Cutty,na.rm=TRUE)
    freshDolly[i] <- sum(props * Dolly,na.rm=TRUE)
    freshPink[i]<- sum(props * Pink,na.rm=TRUE)
    cumul_logging[i] <- sum(props * cum_log,na.rm=TRUE)
    mean_logging[i] <- sum(props * mn_log,na.rm=TRUE)
    
    # lag ocean condition backwards from spawning year
    year_matches <- match(keogh_long$Year[i] - adultLag,keogh_long$Year,nomatch=0)
    year_0 <- year_matches!=0
    year_matches <- sort(year_matches[year_matches!=0])
    props <- 1
    lags <- 1:length(props)
    npgo_lag <- sapply(lags,function(x){mean(keogh_long$npgo[year_matches[1:x]],na.rm=TRUE)})
    oceanSalmon_lag <- sapply(lags,function(x){mean(keogh_long$total[year_matches[1:x]],na.rm=TRUE)})
    mei_lag <- sapply(lags,function(x){mean(keogh_long$mei[year_matches[1:x]],na.rm=TRUE)})
    npgo[i] <- sum(props * npgo_lag,na.rm=TRUE)
    oceanSalmon[i] <- sum(props * oceanSalmon_lag,na.rm=TRUE)
    mei[i] <- sum(props * mei_lag,na.rm=TRUE)
    seals[i] <- sum(props * keogh_long$seal_density[i],na.rm=TRUE)
  }
}
keogh_long$sumTemp <- ifelse(sumTemp==0,NA,sumTemp)
keogh_long$sumRain <- ifelse(sumRain==0,NA,sumRain)
keogh_long$winTemp <- ifelse(winTemp==0,NA,winTemp)
keogh_long$winRain <- ifelse(winRain==0,NA,winRain)
keogh_long$freshCoho <- ifelse(freshCoho==0,NA,freshCoho)
keogh_long$freshSteel <- ifelse(freshSteel==0,NA,freshSteel)
keogh_long$freshCutt <- ifelse(freshCutt==0,NA,freshCutt)
keogh_long$freshDolly <- ifelse(freshDolly==0,NA,freshDolly)
keogh_long$freshPink <- ifelse(freshPink==0,NA,freshPink)
keogh_long$seals <- ifelse(seals==0,NA,seals)
keogh_long$npgo <- ifelse(npgo==0,NA,npgo)
keogh_long$mei <- ifelse(mei==0,NA,mei)
keogh_long$oceanSalmon <- ifelse(oceanSalmon==0,NA,oceanSalmon)
keogh_long$meanLogging <- ifelse(mean_logging==0,NA,mean_logging)
keogh_long$cumLogging <- ifelse(cumul_logging==0,NA,cumul_logging)


plot(log(Recruits/Stock)~log(cumul_footprint),data=keogh_long[keogh_long$Species=="Dolly Varden",],pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))

plot(Recruits~sumTemp,data=keogh_long[keogh_long$Species=="Dolly Varden",],pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))
plot(Recruits~sumRain,data=keogh_long[keogh_long$Species=="Dolly Varden",],pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))
plot(Recruits~winTemp,data=keogh_long[keogh_long$Species=="Dolly Varden",],pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))
plot(Recruits~winRain,data=keogh_long[keogh_long$Species=="Dolly Varden",],pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))
plot(Recruits~meanLogging,data=keogh_long[keogh_long$Species=="Dolly Varden",],pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))
plot(Recruits~meanLogging,data=keogh_long[keogh_long$Species=="Cutthroat",],pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))


plot(Stock~seals,data=keogh_long[keogh_long$Species=="Dolly Varden",],pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))
plot(Stock~npgo,data=keogh_long[keogh_long$Species=="Dolly Varden",],pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))
plot(Stock~mei,data=keogh_long[keogh_long$Species=="Dolly Varden",],pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))
plot(Stock~oceanSalmon,data=keogh_long[keogh_long$Species=="Dolly Varden",],pch=21,bg=ifelse(Year>=1990,"orange","dodgerblue"))


refYear <- 1991
pre89 <- keogh_long[keogh_long$Year<refYear & keogh_long$Species=="Steelhead",]
post89 <- keogh_long[keogh_long$Year>=refYear & keogh_long$Species=="Steelhead",]
ricker_early <- lm(log(Recruits/Stock)~Stock,data=pre89)
ricker_late <- lm(log(Recruits/Stock)~Stock,data=post89)

jpeg("steelhead Ricker breakpoint.jpeg",width=7,height=7,units="in",res=800)
layout(1)
plot(log(Recruits/Stock)~Stock,data=keogh_long[keogh_long$Species=="Steelhead",],pch=21,bg=ifelse(keogh_long$Year[keogh_long$Species=="Steelhead"]<refYear,"dodgerblue","orange"))

with(keogh_long[keogh_long$Species=="Steelhead",], text(log(Recruits/Stock)~Stock, labels = Year,pos=4,font=5,cex=0.75))

abline(ricker_early,col="dodgerblue")
abline(ricker_late,col="orange")

dev.off()
saveRDS(keogh_SR,"Keogh_stockRec_enviro.rds")
saveRDS(keogh_long,"Keogh_SR_enviro_long.rds")
write.csv(keogh_long,"Keogh_StockRecruitment.csv",row.names=F)

# keogh river adult run timing:

keogh <- read.csv("Data/Keogh_Database_Final_01oct18.csv",stringsAsFactors = F,header=T)
keogh[keogh$scale=="y " | keogh$scale=="Y","scale"] <- "y"

keogh$life_stage[keogh$life_stage=="F"] <- "f"
keogh$life_stage[keogh$life_stage==" "] <- ""

keogh$total_age <- sapply(keogh$age_final,function(x){as.numeric(strsplit(x,split="\\.")[[1]][1])+sum(as.numeric(strsplit(gsub("[^[:digit:]]","",strsplit(x,split="\\.")[[1]][2]),split="")[[1]]))+nchar(gsub("[[:digit:]]","",strsplit(x,split="\\.")[[1]][2]))})

keogh$age2 <- as.numeric(keogh$fresh_age)
keogh$age_final <- ifelse(!is.na(keogh$total_age),keogh$age_final,keogh$X1_ager)
keogh$age <- sapply(keogh$age_final,function(x){as.numeric(strsplit(x,split="\\.")[[1]][1])})
keogh$OceanAge <- sapply(keogh$age_final,function(x){sum(as.numeric(strsplit(gsub("[^[:digit:]]","",strsplit(x,split="\\.")[[1]][2]),split="")[[1]]))+nchar(gsub("[[:digit:]]","",strsplit(x,split="\\.")[[1]][2]))})

keogh$AdultStrat <- sapply(keogh$age_final,function(x){strsplit(x,split="\\.")[[1]][2]})
# set ocean ages for all fish, smolts are given age 0 in the ocean
keogh$age_ocean <- as.numeric(keogh$OceanAge)
keogh$age_ocean[keogh$life_stage=="s"] <- 0
# assign brood year to fish as the total freshwater and ocean ages
keogh$hatch_year <- ifelse(keogh$life_stage=="s",keogh$year - keogh$age, keogh$year-(keogh$total_age)) # if adult, brood year is current year minus total age
keogh$smolt_year <- keogh$year-(keogh$age_ocean)

keogh$sampAge <- keogh$age+keogh$age_ocean
adult_life <- individuals[individuals$life_stage=="a"|individuals$life_stage=="k",]
fresh_ocean_ages <- as.factor(paste(adult_life$age,adult_life$age_ocean,sep="|"))

sh_adults <- keogh[keogh$species=="sh" & keogh$life_stage=="a" & keogh$direction!="d" & !is.na(keogh$julian) & (keogh$method=="trap"|keogh$method=="angled"|is.na(keogh$method)|keogh$method==""),]

sh_adults$date <- as.Date(format(strptime(sh_adults$date, format = "%d-%b-%y"), "%Y-%m-%d"))
sh_adults$number <- ifelse(is.na(sh_adults$number),1,sh_adults$number)
base_date <- 319 #julian date of year where run time could begin
sh_adults <- sh_adults[!(sh_adults$julian<=base_date & sh_adults$julian>=250),]
sh_adults$run_year <- ifelse(sh_adults$julian>=base_date,sh_adults$year+1,sh_adults$year)
sh_adults <- sh_adults[sh_adults$run_year>=1976,]
sh_adults$run_date <- ifelse((sh_adults$julian-base_date)<=0,sh_adults$julian-base_date+366,sh_adults$julian-base_date)

sub_keogh <- keogh[!is.na(keogh$julian),]
keogh_run_year <- ifelse(sub_keogh$julian>=base_date,sub_keogh$year+1,sub_keogh$year)
keogh_base_dates <- ifelse((sub_keogh$julian-base_date)<=0,sub_keogh$julian-base_date+366,sub_keogh$julian-base_date)
years <- sort(unique(keogh$year))
keogh_samp_start <- sapply(years,function(x){min(keogh_base_dates[which(keogh_run_year==x & sub_keogh$trap_fishing=="TRUE")],na.rm=TRUE)})
names(keogh_samp_start) <- years
keogh_samp_start["1985"] <- min(keogh_base_dates[which(keogh_run_year==1985 & sub_keogh$species_code=="sha")],na.rm=TRUE)
prop_run <- sapply(years,function(x){sum(sh_adults$run_date[sh_adults$year==x]>=keogh_samp_start[names(keogh_samp_start)==x])/length((sh_adults$run_date[sh_adults$year==x]))})
sh_adults$sex_est <- ifelse(grepl("f",sh_adults$sex),1,ifelse(grepl("m",sh_adults$sex),0,NA))
# create new dataframe repeating number of adults per run date:
adult_run <- data.frame("year"=rep(sh_adults$run_year,sh_adults$number),"date"=rep(sh_adults$date,sh_adults$number),"run"=rep(sh_adults$run_date,sh_adults$number),"size"=rep(sh_adults$fork_length,sh_adults$number),"age"=rep(sh_adults$age_ocean,sh_adults$number),"sex"=rep(sh_adults$sex_est,sh_adults$number))
adult_run <- adult_run[order(adult_run$date),]

# Bind the rows together
adult_run <- merge(adult_run,portHardy_pg[,c("date","mean_temp_run","max_temp_run","min_temp_run","total_rain_run","mean_temp_egg","total_rain_egg")],by=c("date"))
adult_run$cum_sum <- unlist(sapply(unique(adult_run$year),function(x){(1:sum(adult_run$year==x))/sum(adult_run$year==x)}))
adult_run$peak <- unlist(sapply(unique(adult_run$year),function(x){rep(table(adult_run$run[adult_run$year==x])/max(table(adult_run$run[adult_run$year==x])),table(adult_run$run[adult_run$year==x]))}))
#sh_adults$run_date <- sh_adults$run_date-min(sh_adults$run_date,na.rm=TRUE)

sapply(c(1977,1987,1997,2007),function(x){sum(adult_run$year==x & adult_run$run<=100)/sum(adult_run$year==x)})
sapply(c(1977,1987,1997,2007),function(x){sum(adult_run$year==x & adult_run$run<=100)})
sapply(c(1977,1987,1997,2007),function(x){median(adult_run$run[adult_run$year==x],na.rm=TRUE)})


run_time <- aggregate(cbind(run,size,age,peak,sex,mean_temp_run,max_temp_run,min_temp_run,total_rain_run,mean_temp_egg,total_rain_egg)~year,data=adult_run,FUN=mean,na.rm=TRUE)
run_time$run <- aggregate(run~year,data=adult_run,FUN=median,na.rm=TRUE)$run
run_time$N <- aggregate(run~year,data=adult_run,FUN=length)$run
run_time$run_f <- aggregate(run~year,data=adult_run[adult_run$sex==1,],FUN=median,na.rm=TRUE)$run
run_time$run_m <- aggregate(run~year,data=adult_run[adult_run$sex==0,],FUN=median,na.rm=TRUE)$run
run_time$samp_start <- keogh_samp_start[match(run_time$year,names(keogh_samp_start))]
run_time$since_start <- run_time$run-run_time$samp_start
run_time$prop_after_start <- prop_run[match(run_time$year,names(keogh_samp_start))]

run_time$fence_samps <- run_time$prop_after_start*run_time$N
plot(run_time$fence_samps)
layout(matrix(1))
jpeg("Figures/correcting run date through time.jpeg",width=5.5,height=7.5,units="in",res=800)
layout(matrix(1:3,nrow=3))
par(mar=c(5,5,0.1,0.1))
plot(prop_after_start~samp_start,data=run_time[run_time$year>=1997 | run_time$year<1978,],ylim=c(0,1),xlab="Fence install date (starting Nov. 15th)",ylab="Proportion of run\n(after fence install)")
m1 <- glm(prop_after_start~samp_start,data=run_time[run_time$year>=1997 | run_time$year<1978,],family=binomial(link="logit"),weights=run_time$N[run_time$year>=1997 | run_time$year<1978])
curve(exp(coef(m1)[1]+coef(m1)[2]*x)/(1+exp(coef(m1)[1]+coef(m1)[2]*x)),add=TRUE)
summary(m1)
x1=-(coef(m1)[1]/coef(m1)[2])
odds <- 0.75
ybar <- (1/odds-1)
x80=(-log(ybar)-(coef(m1)[1]))/coef(m1)[2]
abline(v=x80)
abline(h=odds)

correction <- predict(m1,newdata = run_time,type="response")

plot(since_start~samp_start,data=run_time[run_time$year>=1997 | run_time$year<1978,],xlab="Fence install date (starting Nov. 15th)",ylab="\u0394 median run and install date")
m0 <- lm(since_start~samp_start,data=run_time[run_time$year>=1997 | run_time$year<1978,],weights=run_time$N[run_time$year>=1997 | run_time$year<1978])
summary(m0)
abline(m0)
xcorrect <- predict(m0,newdata=run_time)
run_time$run_date_corrected <- (correction*run_time$run)+(1-correction)*(xcorrect+run_time$samp_start)

plot(run_time$year,run_time$run_date_corrected,ylim=c(0,190),type="b",pch=21,bg="grey50",xlab="Year",ylab="Median run date (since Nov. 15th)")
lines(run_time$year,run_time$run,type="b",pch=21,bg="dodgerblue")
legend("bottomright",c("Corrected","Observed"),pch=21,lty=1,col=c("black","black"),pt.bg=c("grey50","dodgerblue"),bty="n")
dev.off()

plot(samp_start~year,data=run_time)
plot(run_time$since_start,run_time$run)

jpeg("Figures/run date through time.jpeg",width=6,height=5,units="in",res=800)
layout(1)
par(mar=c(5,4,1,1))
plot(run~year,data=adult_run,pch=21,bg="grey80",xlab="Brood year",ylab="Run date (starting Nov. 15th)")
lines(run_time$year,run_time$run,lwd=3,col="black")
lines(run_time$year,run_time$run_f,lwd=3,col="orange")
lines(run_time$year,run_time$run_m,lwd=3,col="dodgerblue")

dev.off()

sh_SR <- keogh_long[keogh_long$Species=="Steelhead" & keogh_long$Year>=1976 & keogh_long$Year<=2017,]
colnames(sh_SR)[colnames(sh_SR)=="Year"] <- "year"
run_time <- merge(run_time,sh_SR,by="year")

run_time$logit_surv <- log((run_time$Stock/run_time$juvCohort)/(1-(run_time$Stock/run_time$juvCohort)))

jpeg("Figures/drivers of steelhead run date and smolt production.jpeg",width=6,height=5,units="in",res=800)
layout(matrix(1:6,ncol=2,nrow=3,byrow=TRUE))
par(mar=c(5,4,1,1))
plot(run_time$seals,run_time$logit_surv,ylab="Marine survival (logit)",xlab="Seal density")
abline(lm(logit_surv~seals,data=run_time),lwd=2,col="red")

plot(run_time$oceanSalmon,run_time$logit_surv,ylab="Marine survival (logit)",xlab="Pacific salmon biomass")
abline(lm(logit_surv~oceanSalmon,data=run_time),lwd=2,col="red")

plot(run_time$logit_surv,run_time$run_date_corrected,xlab="Marine survival (logit)",ylab="Average run date")
m1 <- lm(run~logit_surv,data=run_time)
m2 <- nls(run~Asym/(1 + exp((xmid - logit_surv))),data=run_time,start=list(Asym=140,xmid=-4))
summary(m2)
AIC(m1,m2)
curve(coef(m2)["Asym"]/(1 + exp((coef(m2)["xmid"] - x))),add=TRUE,lwd=2,col="red")

plot(run_time$total_rain_run,run_time$run_date_corrected,xlab="Mean rainfall 14 days before run (mm)",ylab="Average run date")
lines(smooth.spline(run_time$total_rain_run,run_time$run),lwd=2,col="red")

plot(run_time$run_date_corrected,run_time$Recruits,xlab="Average run date",ylab="Smolts")
lines(smooth.spline(run_time$run_date_corrected[!is.na(run_time$Recruits)],run_time$Recruits[!is.na(run_time$Recruits)],cv=TRUE),lwd=2,col="red")

plot(run_time$total_rain_egg,run_time$Recruits,xlab="Mean rainfall 30 days after run (mm)",ylab="Smolts")
lines(smooth.spline(run_time$total_rain_egg[!is.na(run_time$Recruits)],run_time$Recruits[!is.na(run_time$Recruits)],cv=TRUE),lwd=2,col="red")

dev.off()

yr_run <- aggregate(run~year,data=adult_run,FUN=mean)
yr_run$sd <- aggregate(run~year,data=adult_run,FUN=function(x){sd(x)/mean(x)})$run
plot(sd~year,data=yr_run)
plot(sd~run,data=yr_run)
adult_run$year <- factor(adult_run$year,levels=rev(unique(sort(adult_run$year))))
run_time_year <- factor(run_time$year,levels=rev(unique(sort(run_time$year))))
ggplot(adult_run, aes(x = run, y=as.factor(year),fill=..x..)) + 
  geom_density_ridges_gradient(scale=2.5,rel_min_height=1e-3) +
  scale_x_continuous(name="Run date (starting Nov. 15th)",expand=c(0.01,0)) +
  scale_y_discrete(expand=c(0.01,0)) +
  geom_segment(data = run_time, aes(x = samp_start, xend = samp_start, y = as.numeric(run_time_year),yend = as.numeric(run_time_year) + 1.5),color = "black",lwd=0.5) +
  #scale_fill_viridis(name="Run time",option="C") +
  theme_ridges(font_size=13,grid=TRUE,center=TRUE) + 
  theme(axis.title.y=element_blank()) +
  scale_fill_gradient2(name="Run date",low="blue",high="orange",mid="dodgerblue",midpoint=75)

ggsave(filename="Figures/runtime v2.jpeg",units="in",height=7,width=6,dpi=800)

oceanCovar <- cbind(scale(run_time[,c("seals","oceanSalmon")],center=T))
pca <- princomp(oceanCovar)
summary(pca)
pca$loadings
run_time$ocean_interact <- pca$scores[,1]
run_time$ocean_covar_2 <- pca$scores[,2]

saveRDS(adult_run,"Data/steelhead_run.rds")
saveRDS(run_time,"Data/steelhead_run_annual.rds")
