library(reshape2)
source("some functions.R")
data_check <- "new"

co_lag <- 1 # fixed freshwater residency for coho: 1 year on average from Wade & Irvine 2018 report from Keogh and Holtby et al. 1990 CJFAS paper from Carnation Creek
co_smolt_lag <- 2
ct_lag <- 2 # fixed freshwater residency for cutthroat: 2-4 Armstrong 1971 TAFS, Trotter 1989 suggests age 2 for estuarine/coastal:  Losee et al. (2018) suggests age 2 for coastal cutties
# Smith 1980 BC FLNRO suggests age 3 for Keogh River
ct_lags <- 1:4
ct_lags_prop <- c(0.32,0.45,0.20,0.03)
ct_smolt_lags <- 2:4
ct_smolt_lags_prop <- c(1/3/2,1/3/2,2/3)
dv_lag <- 4 # freshwater residence: Armstrong 1970 FRBC and Dolloff & Reeves 1990 CJFAS suggest ~ age 1-4 for dolly varden smoltification
dv_lags <- 2:4
dv_lags_prop <- c(0.11,0.80,0.09) # Smith and Slaney 1980
dv_smolt_lags <- 0:4
dv_smolt_lags_prop <- c(8,37,38,11,7)/sum(c(8,37,38,11,7))
ch_lag <- 4 # fixed 4 year life cycle for chum: Neave et al. 1952
pink_lag <- 2 # fixed 2 year life cycle

environ <- readRDS("environ_covars.rds")

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

# adults with known ocean ages have a known smolt year: they should count as 'abundant' for their smolt year
row.names(annual_prop_sh) <- row.names(annual_cohort_prop) <- annual_cohort$Year
colnames(annual_prop_sh) <- min(annual_cohorts$Hatch):max(annual_cohorts$Hatch)
annual_abund <- dcast(annual,Year~Species+Stage,value.var="Abundance")
annual_size <- dcast(size,Year~Species+Stage,value.var="Length")
annual_sigma <- dcast(size,Year~Species+Stage,value.var="Sigma")

# we want to find how many of the smolts sampled in Year X were from the spawners at Year Y
# to do that, we find the proportions aged and multiply the current year smolts by that proportion from hatch year Y

mn_smolt_age <- table(keogh$age[keogh$species=="sh" & (keogh$life_stage=="s" | keogh$life_stage=="k" | keogh$life_stage=="a")])
mn_smolt_age <- mn_smolt_age/sum(mn_smolt_age)

years <- min(keogh$hatch_year,na.rm=TRUE):max(keogh$year,na.rm=TRUE)
annual_smolts_sh <- rep(NA,length(years))
names(annual_smolts_sh) <- years
prop_smolts_year <- matrix(NA,nrow=nrow(stock_rec),ncol=length(mn_smolt_age),dimnames=list("Year"=stock_rec$Year,"Ages"=1:length(mn_smolt_age)))
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

keogh_smolts <- merge(keogh_atlas[,c("Year","Total")],stock_rec,by="Year",all=T)

# compile cutthroat trout data

ct_R <- aggregate(number~year+life_stage,data=keogh[keogh$species=="ct",],FUN=sum,na.rm=T)
ct_age <- aggregate(number~year+fresh_age+life_stage,data=keogh[keogh$species=="ct",],FUN=sum,na.rm=T)
CT_annual <- dcast(ct_R,year~life_stage,value.var="number")
colnames(CT_annual)[colnames(CT_annual)=="year"] <- "Year"
CT_annual <- merge(stock_rec,CT_annual,by="Year",all=TRUE)
CT_annual <- CT_annual[,-match(c("sh_Adults","sh_Smolts"),colnames(CT_annual))]

#ct_SR <- data.frame("Year"=CT_annual$Year[1:(length(CT_annual$Year)-ct_lag)],"Stock"=CT_annual$a[1:(length(CT_annual$Year)-ct_lag)],"Recruits"=CT_annual$s[(ct_lag+1):length(CT_annual$Year)])

ct_SR <- data.frame("Year"=CT_annual$Year,"Stock"=CT_annual$a,"Recruits"=c(CT_annual$s[(ct_lag+1):length(CT_annual$Year)],rep(NA,ct_lag)))
ct_SR$Juv_cohort <- c(rep(NA,5),ct_SR$Recruits[1:(length(CT_annual$Year)-(5))])

ct_SR$Recruits <- round(rowSums(sapply(1:length(ct_lags),function(x){c(CT_annual$s[(ct_lags[x]+1):length(CT_annual$Year)]*ct_lags_prop[x],rep(NA,ct_lags[x]))}),na.rm=TRUE),0)

ct_SR$Juv_cohort <- round(rowSums(sapply(1:length(ct_smolt_lags),function(x){c(rep(NA,ct_smolt_lags[x]),ct_SR$Recruits[1:(length(CT_annual$Year)-(ct_smolt_lags[x]))]*ct_smolt_lags_prop[x])}),na.rm=TRUE),0)
ct_SR[ct_SR==0] <- NA
plot(ct_SR$Stock,ct_SR$Recruits)

# dolly vardens
dv_R <- aggregate(number~year+life_stage,data=keogh[keogh$species=="dv",],FUN=sum,na.rm=T)
dv_age <- aggregate(number~year+fresh_age+life_stage,data=keogh[keogh$species=="dv",],FUN=sum,na.rm=T)
DV_annual <- dcast(dv_R,year~life_stage,value.var="number")
colnames(DV_annual)[colnames(DV_annual)=="year"] <- "Year"
DV_annual <- merge(stock_rec,DV_annual,by="Year",all=TRUE)
DV_annual <- DV_annual[,-match(c("sh_Adults","sh_Smolts"),colnames(DV_annual))]

dv_SR <- data.frame("Year"=DV_annual$Year,"Stock"=DV_annual$a,"Recruits"=c(DV_annual$s[(dv_lag+1):length(DV_annual$Year)],rep(NA,dv_lag)))

dv_SR$Recruits <- round(rowSums(sapply(1:length(dv_lags),function(x){c(DV_annual$s[(dv_lags[x]+1):length(DV_annual$Year)]*dv_lags_prop[x],rep(NA,dv_lags[x]))}),na.rm=TRUE),0)

dv_SR$Juv_cohort <- round(rowSums(sapply(1:length(dv_smolt_lags),function(x){c(rep(NA,dv_smolt_lags[x]),dv_SR$Recruits[1:(length(DV_annual$Year)-(dv_smolt_lags[x]))]*dv_smolt_lags_prop[x])}),na.rm=TRUE),0)
dv_SR[dv_SR==0] <- NA

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


# compile data on chum salmon
ch_SR <- data.frame("Year"=chum$Year,"Stock"=chum$Adults)
ch_SR$Recruits <- c(chum$Adult[(ch_lag+1):length(chum$Year)],rep(NA,ch_lag))
ch_SR$Juv_cohort <- c(rep(NA,ch_lag),ch_SR$Stock[1:(length(ch_SR$Year)-(ch_lag))])

# read in and calculate pink salmon stock-recruitment
pinks <- data.frame("Year"=pinks$Year[-((length(pinks$Stock)-pink_lag+1):length(pinks$Stock))],"Stock"=pinks$Stock[-((length(pinks$Stock)-pink_lag+1):length(pinks$Stock))],"Recruits"=pinks$Stock[-(1:pink_lag)])

pinks <- data.frame("Year"=pinks$Year,"Stock"=pinks$Stock)
pinks$Recruits <- c(pinks$Stock[(pink_lag+1):length(pinks$Year)],rep(NA,pink_lag))
pinks$Juv_cohort <- c(rep(NA,pink_lag),pinks$Stock[1:(length(pinks$Year)-(pink_lag))])

# to include or not include pre-1997 coho
co_SR$Stock[co_SR$Year<=1997] <- NA

# compile data for Jordan: two data frames
# First: Annual time-series for brood year stock-recruit
# Second: Annual time-series for outmigrating year, smolt abundance, smolt body size, smolt age-structure, spawners
smolt_size <- aggregate(fork_length~smolt_year,data=keogh[keogh$species=="sh" & keogh$life_stage=="s",],FUN=mean,na.rm=T)
smolt_sd <- aggregate(fork_length~smolt_year,data=keogh[keogh$species=="sh" & keogh$life_stage=="s",],FUN=sd,na.rm=T)
smolt_size$sigma <- smolt_sd$fork_length
colnames(smolt_size) <- c("Year","Fork length (mm)","Std. dev.")

smolt_outmigrate <- merge(steelDailyMax,smolt_size,by="Year",all=T)
steel_outmigrate <- merge(keogh_adults,smolt_outmigrate,by="Year",all=T)
steel_age <- dcast(steel_age[steel_age$Species=="sh" & steel_age$Stage=="s",],Year~Age,value.var="Abundance")

steel_outmigrate <- merge(steel_outmigrate,steel_age,by="Year",all=T)
colnames(steel_outmigrate) <- c("Year","Spawner abundance (mark-recap)","Smolt abundance (fence count)","Smolt fork length (mm)","Smolt FL sd","Age 1","Age 2","Age 3","Age 4","Age 5")

steel_outmigrate <- steel_outmigrate[steel_outmigrate$Year>=1975,]

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

example <- reshape(keogh_SR,direction = "long",varying = list(c("sh_Adults","dv_Adults","ct_Adults","pk_Adults","ch_Adults","co_Adults"),c("sh_Smolts","dv_Smolts","ct_Smolts","pk_Recruits","ch_Recruits","co_Smolts"),c("sh_juv_Cohort","dv_juv_Cohort","ct_juv_Cohort","pk_juv_Cohort","ch_juv_Cohort","co_juv_Cohort")),v.names=c("Stock","Recruits","juvCohort"),idvar="Species")

example$Species <- rep(c("Steelhead","Dolly Varden","Cutthroat","Pink","Chum","Coho"),each=nrow(keogh_SR))
example$Species <- factor(example$Species,levels=c("Steelhead","Dolly Varden","Cutthroat","Pink","Chum","Coho"))

keogh_long <- example

seals <- rep(NA,nrow(keogh_long))
names(seals) <- keogh_long$Year
sumTemp <- sumRain <- winTemp <- winRain <- freshCoho <- freshSteel <- freshCutt <- freshDolly <- freshPink <- npgo <- oceanSalmon <- mei <- seals

for(i in 1:nrow(keogh_long))
{
  spp <- keogh_long$Species[i]
  if(spp=="Steelhead"){
    # lag freshwater conditions forwards from hatch year
    if(any(as.numeric(row.names(prop_smolts_year)) %in% keogh_long$Year[i])){
      year_matches <- sort(match(keogh_long$Year[i] + as.integer(colnames(prop_smolts_year)),keogh_long$Year,nomatch=0))
      age_matches <- match(as.numeric(row.names(prop_smolts_year)),keogh_long$Year[i],nomatch=0)
      pink_matches <- match(keogh_long$Year[year_matches],keogh_long$Year[keogh_long$Species=="Pink"],nomatch=0)
      sumTemp[i] <- sum(prop_smolts_year[age_matches,] * keogh_long$max_temp[year_matches],na.rm=TRUE)
      sumRain[i]<- sum(prop_smolts_year[age_matches,] * keogh_long$total_rain[year_matches],na.rm=TRUE)
      winTemp[i]<- sum(prop_smolts_year[age_matches,] * keogh_long$win_mean_temp[year_matches],na.rm=TRUE)
      winRain[i]<- sum(prop_smolts_year[age_matches,] * keogh_long$win_rain[year_matches],na.rm=TRUE)
      freshCoho[i] <- sum(prop_smolts_year[age_matches,] * keogh_long$Recruits[keogh_long$Species=="Coho"][pink_matches],na.rm=TRUE)
      freshSteel[i] <- sum(prop_smolts_year[age_matches,] * keogh_long$Recruits[keogh_long$Species=="Steelhead"][pink_matches],na.rm=TRUE)
      freshCutt[i] <- sum(prop_smolts_year[age_matches,] * keogh_long$Recruits[keogh_long$Species=="Cutthroat"][pink_matches],na.rm=TRUE)
      freshDolly[i] <- sum(prop_smolts_year[age_matches,] * keogh_long$Recruits[keogh_long$Species=="Dolly Varden"][pink_matches],na.rm=TRUE)
      freshPink[i]<- sum(prop_smolts_year[age_matches,] * keogh_long$Recruits[keogh_long$Species=="Pink"][pink_matches],na.rm=TRUE)
    }
    # lag ocean condition backwards from spawning year
    if(any(as.numeric(row.names(prop_adults_year)) %in% keogh_long$Year[i])){
      year_matches <- sort(match(keogh_long$Year[i] - as.integer(colnames(prop_adults_year)),keogh_long$Year,nomatch=0))
      age_matches <- match(as.numeric(row.names(prop_adults_year)),keogh_long$Year[i],nomatch=0)
      pink_matches <- match(keogh_long$Year[year_matches],keogh_long$Year[keogh_long$Species=="Pink"],nomatch=0)
      
      npgo[i] <- sum(prop_adults_year[age_matches,] * keogh_long$npgo[year_matches],na.rm=TRUE)
      oceanSalmon[i] <- sum(prop_adults_year[age_matches,] * keogh_long$total[year_matches],na.rm=TRUE)
      mei[i] <- sum(prop_adults_year[age_matches,] * keogh_long$mei[year_matches],na.rm=TRUE)
      seals[i] <- sum(prop_adults_year[age_matches,] * keogh_long$seal_density[year_matches],na.rm=TRUE)
    }
    
  }
  if(spp=="Pink"){
    juvLag <- 0
    adultLag <- pink_lag
    sumTemp[i]
    sumRain[i]
    winTemp[i]
    winRain[i]
    freshComp[i]
    freshPink[i]
    npgo[i]
    oceanSalmon[i]
    mei[i]
    seals[i]
  }
  if(spp=="Dolly Varden"){
    sumTemp[i]
    sumRain[i]
    winTemp[i]
    winRain[i]
    freshComp[i]
    freshPink[i]
    npgo[i]
    oceanSalmon[i]
    mei[i]
    seals[i]
  }
  if(spp=="Cutthroat"){
    sumTemp[i]
    sumRain[i]
    winTemp[i]
    winRain[i]
    freshComp[i]
    freshPink[i]
    npgo[i]
    oceanSalmon[i]
    mei[i]
    seals[i]
  }
  if(spp=="Coho"){
    sumTemp[i]
    sumRain[i]
    winTemp[i]
    winRain[i]
    freshComp[i]
    freshPink[i]
    npgo[i]
    oceanSalmon[i]
    mei[i]
    seals[i]
  }
  if(spp=="Chum"){
    sumTemp[i]
    sumRain[i]
    winTemp[i]
    winRain[i]
    freshComp[i]
    freshPink[i]
    npgo[i]
    oceanSalmon[i]
    mei[i]
    seals[i]
  }
}

saveRDS(keogh_SR,"Keogh_stockRec_enviro.rds")
write.csv(keogh_StockRec,"Keogh_StockRecruitment.csv",row.names=F)
