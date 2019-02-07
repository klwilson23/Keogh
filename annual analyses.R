library(reshape2)

Corner_text <- function(text, location="topright",...)
{
  legend(location,legend=text, bty ="n", pch=NA,...)
}

keogh <- read.csv("Data/Keogh_Database_Final_01oct18.csv",stringsAsFactors = F,header=T)
#keogh <- read.csv("Data/Keogh_Database_Final_19nov18.csv",stringsAsFactors = F,header=T)
keogh[keogh$scale=="y " | keogh$scale=="Y","scale"] <- "y"
keogh_atlas <- read.csv("Data/Keogh sh smolts Atlas 2015.csv",stringsAsFactors=F,header=T)
keogh_instream <- read.csv("Data/Keogh smolts instream outmigration.csv",stringsAsFactors=F,header=T)

keogh_adults <- read.csv("Data/Keogh sh adults.csv",stringsAsFactors=F,header=T)

pinks <- read.csv("Data/Keogh_Pink_Salm.csv",stringsAsFactors = F,header=T)
pinks <- pinks[order(pinks$Year),]
pinks$Stock[pinks$Year<=1997] <- 2*pinks$Stock[pinks$Year<=1997] # from Bailey et al. 2018 - Pink salmon prior to 1997 sampled by stream walks and need to be doubled to correct for abundance patterns. After 1997 pink salmon counted by resistivitycounter.
coho <- read.csv("Data/Keogh coho adults.csv",stringsAsFactors = F,header=T)
coho <- coho[order(coho$Year),]
coho$Adults[coho$Year<=1997] <- NA#2*coho$Adults[coho$Year<=1997] # from Bailey et al. 2018 - Pink salmon prior to 1997 sampled by stream walks and need to be doubled to correct for abundance patterns. After 1997 pink salmon counted by resistivitycounter.

chum <- read.csv("Data/Keogh chum adults.csv",stringsAsFactors = F,header=T)
chum <- chum[order(chum$Year),]
chum$Adults[chum$Adults==0] <- NA
chum$Adults[chum$Year<=1997] <- 2*chum$Adults[chum$Year<=1997] # from Bailey et al. 2018 - Pink salmon prior to 1997 sampled by stream walks and need to be doubled to correct for abundance patterns. After 1997 pink salmon counted by resistivitycounter.

# read in coho smolt data from Tom Johnston dataset
coSm <- read.csv("Data/Keogh coho smolts.csv",stringsAsFactors = F,header=T)
coSm <- coSm[order(coSm$Year),]

# read in and calculate pink salmon stock-recruitment
pinks <- data.frame("Year"=pinks$Year,"Stock"=c(pinks$Stock[-1],NA),"Recruits"=pinks$Stock)
plot(pinks$Stock,pinks$Recruits)
keogh <- keogh[keogh$hatchery!=1|is.na(keogh$hatchery),]

keogh$life_stage[keogh$life_stage=="F"] <- "f"
keogh$life_stage[keogh$life_stage==" "] <- ""

keogh$total_age <- sapply(keogh$age_final,function(x){as.numeric(strsplit(x,split="\\.")[[1]][1])+sum(as.numeric(strsplit(gsub("[^[:digit:]]","",strsplit(x,split="\\.")[[1]][2]),split="")[[1]]))+nchar(gsub("[[:digit:]]","",strsplit(x,split="\\.")[[1]][2]))})

keogh$age2 <- as.numeric(keogh$fresh_age)
keogh$age_final <- ifelse(!is.na(keogh$total_age),keogh$age_final,keogh$X1_ager)
keogh$age <- sapply(keogh$age_final,function(x){as.numeric(strsplit(x,split="\\.")[[1]][1])})
keogh$OceanAge <- sapply(keogh$age_final,function(x){sum(as.numeric(strsplit(gsub("[^[:digit:]]","",strsplit(x,split="\\.")[[1]][2]),split="")[[1]]))+nchar(gsub("[[:digit:]]","",strsplit(x,split="\\.")[[1]][2]))})

# example code for validating the string splits for age data when age reads a character sequence like 3.1s2s

#keogh[keogh$year==1985 & keogh$species=="sh" & keogh$life_stage=="a",c("age_final","age","ocean_age","OceanAge")]

#as.numeric(strsplit(keogh$age_final[98301],split="\\.")[[1]][2])

#sum(as.numeric(strsplit(gsub("[^[:digit:]]","",strsplit(keogh$X1_ager[185116],split="\\.")[[1]][2]),split="")[[1]]))+nchar(gsub("[[:digit:]]","",strsplit(keogh$X1_ager[185116],split="\\.")[[1]][2]))

#keogh <- keogh[!keogh$species_code%in%c("cf","ctt","ctr","ks","shweres","shlgbres","cot"),]

# find out sample size for age data by year & life stage

sampSize <- aggregate(number~year+life_stage,data=keogh[keogh$species=="sh"&!is.na(keogh$age),],FUN=function(x){sum(x)})

keogh$scaleLog <- ifelse(keogh$scale=="y",1,0)
keogh$scaleN <- keogh$scaleLog*keogh$number
sampSizev2 <- aggregate(scaleLog~year+life_stage,data=keogh[keogh$species=="sh"&!is.na(keogh$age),],FUN=function(x){sum(x)})

ageEstimates <- dcast(sampSize,year~life_stage,value.var="number")
ageEstimatesv2 <- dcast(sampSizev2,year~life_stage,value.var="scaleLog")

cbind(ageEstimates$year,ageEstimates$s,ageEstimatesv2$s)

# set ocean ages for all fish, smolts are given age 0 in the ocean

keogh$age_ocean <- as.numeric(keogh$OceanAge)
keogh$age_ocean[keogh$life_stage=="s"] <- 0

# assign brood year to fish as the total freshwater and ocean ages

keogh$hatch_year <- ifelse(keogh$life_stage=="s",keogh$year - keogh$age, keogh$year-(keogh$total_age)) # if adult, brood year is current year minus total age
keogh$smolt_year <- keogh$year-(keogh$age_ocean)

annual_age <- aggregate(number~smolt_year+hatch_year+species+life_stage+age,data=keogh,FUN=sum,na.rm=T)
colnames(annual_age) <- c("Year","Hatch","Species","Stage","Age","Abundance")
annual_age <- subset(annual_age,Stage%in%c("s","a","k"))
annual_cohorts <- aggregate(Abundance~Hatch+Year+Species+Age,data=annual_age,FUN=sum,na.rm=T)

steel_age <- subset(annual_age,Species=="sh"&Stage=="s")

dailyMax <- NULL

for(i in unique(keogh$year))
{
  subKeogh <- keogh[keogh$year==i & keogh$species=="sh" & keogh$life_stage=="s" & keogh$alive!="false",]
  for(j in unique(subKeogh$julian))
  {
    subDaily <- subKeogh[subKeogh$julian==j,]
    dailyMax <- rbind(dailyMax,subDaily[which(subDaily$number==max(subDaily$number,na.rm=T)),])
  }
}

steelDailyMax <- aggregate(number~year,dailyMax,FUN=sum,na.rm=T)
colnames(steelDailyMax) <- c("Year","DailyMax")

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

# adults with known ocean ages have a known smolt year: they should count as 'abundant' for their smolt year
annual_outmigAdults <- subset(annual_age,Stage%in%c("a") & Species%in%c("sh"))
annual_outmigAdults <- aggregate(Abundance~Hatch+Year+Species+Age,data=annual_outmigAdults,FUN=sum,na.rm=T)
annual_outmigAdults <- dcast(annual_outmigAdults,Year~Species+Hatch,value.var="Abundance",fun.aggregate=sum)
yrs <- annual_outmigAdults$Year
annual_outmigAdults <- rowSums(annual_outmigAdults[,grep("sh",colnames(annual_outmigAdults))])

names(annual_outmigAdults) <- yrs

row.names(annual_prop_sh) <- row.names(annual_cohort_prop) <- annual_cohort$Year
colnames(annual_prop_sh) <- min(annual_cohorts$Hatch):max(annual_cohorts$Hatch)
annual_abund <- dcast(annual,Year~Species+Stage,value.var="Abundance")
annual_size <- dcast(size,Year~Species+Stage,value.var="Length")
annual_sigma <- dcast(size,Year~Species+Stage,value.var="Sigma")

# we want to find how many of the smolts sampled in Year X were from the spawners at Year Y
# to do that, we find the proportions aged and multiply the current year smolts by that proportion from hatch year Y
annual_sh <- subset(annual,Species=="sh")
annual_sh_abund <- dcast(annual_sh,Year~Species+Stage,value.var="Abundance")

keogh_smolt <- merge(data.frame("Year"=as.numeric(names(annual_outmigAdults)),"Smolts"=annual_outmigAdults),annual_sh_abund[,c("Year","sh_s")],by="Year",all=T)
keogh_smolt <- data.frame(keogh_smolt$Year,rowSums(keogh_smolt[,-1],na.rm=T))
colnames(keogh_smolt) <- c("Year","Smolts")

keogh_smoltv2 <- merge(data.frame("Year"=as.numeric(names(annual_outmigAdults)),"Smolts"=annual_outmigAdults),steelDailyMax,by="Year",all=T)
keogh_smoltv2 <- data.frame(keogh_smoltv2$Year,rowSums(keogh_smoltv2[,-1],na.rm=T))
colnames(keogh_smoltv2) <- c("Year","Smolts")

annual_smolts_sh <- annual_smolts_shv2 <- rep(NA,length(colnames(annual_prop_sh)))
names(annual_smolts_sh) <- names(annual_smolts_shv2) <- colnames(annual_prop_sh)
for(i in as.numeric(colnames(annual_prop_sh)))
{
  IIbrood <- which(as.numeric(colnames(annual_prop_sh))==i)
  IIsample <- grep(paste("^",row.names(annual_prop_sh)[annual_prop_sh[,IIbrood]>0],"$",collapse="|",sep=""),keogh_smolt$Year)
  IIage <- which(annual_prop_sh[,IIbrood]>0)

  annual_smolts_sh[IIbrood] <- sum(keogh_smolt$Smolts[IIsample]*annual_prop_sh[IIage,IIbrood],na.rm=T)
  annual_smolts_shv2[IIbrood] <- sum(keogh_smoltv2$Smolts[IIsample]*annual_prop_sh[IIage,IIbrood],na.rm=T)
}

stock_rec <- data.frame("Year"=as.numeric(names(annual_smolts_shv2)),"Smolts"=annual_smolts_shv2)
stock_rec$Smolts[stock_rec$Smolts==0] <- NA

stock_rec <- merge(keogh_adults,stock_rec,by="Year",all=T)

saveRDS(stock_rec,"steelhead_stockRec.rds")

keogh_sh <- readRDS("steelhead_stockRec.rds")
keogh_sh_v2 <- readRDS("steelhead_stockRec_v2.rds")

keogh_smolts <- merge(keogh_atlas[,c("Year","Total")],stock_rec,by="Year",all=T)

plot(keogh_smolts$Year,keogh_smolts$Total,type="l",xlab="Year",ylab="Smolt abundance",ylim=c(0,max(keogh_smolts,na.rm=T)))
points(keogh_smolts$Year,keogh_smolts$Smolts,pch=21,bg="dodgerblue",type="b",lty=2,lwd=1)
#points(keogh_sh_v2$Year,keogh_sh_v2$Smolts,pch=21,bg="orange",ylim=c(0,15000),type="b",lty=1,lwd=1,col="orange")
legend("topright",c("Tom Johnston","Current QA/QC"),col=c("black","dodgerblue"),pch=c(NA,21),pt.bg=c(NA,"dodgerblue"),lty=c(1,2),lwd=c(2,2),bty="n")

plot(stock_rec$Adults,log(stock_rec$Smolts/stock_rec$Adults))
sh_SR <- lm(log(stock_rec$Smolts/stock_rec$Adults)~stock_rec$Adults)
abline(sh_SR)

plot(stock_rec$Adults,stock_rec$Smolts,xlab="Adults",ylab="Smolts",ylim=c(0,max(stock_rec$Smolts,na.rm=T)))
lines(smooth.spline(stock_rec$Adults[complete.cases(stock_rec)],stock_rec$Smolts[complete.cases(stock_rec)],cv=T),lwd=2,lty=2,col="dodgerblue")
curve(exp(coef(sh_SR)[1])*x*exp(coef(sh_SR)[2]*x),add=T,from=0,to=max(stock_rec$Adults,na.rm=T))

# analysis on pink salmon
plot(pinks$Stock,pinks$Recruits)
pk_SR <- lm(log(pinks$Recruits/pinks$Stock)~pinks$Stock)
curve(exp(coef(pk_SR)[1])*x*exp(coef(pk_SR)[2]*x),add=T,from=0,to=max(pinks$Stock,na.rm=T))

plot(stock_rec$Year,stock_rec$Adults/stock_rec$Smolts,xlab="Year",ylab="Marin survival (%)",type="b")


matplot(stock_rec[,2:3],type="l")

# compile cutthroat trout data

ct_R <- aggregate(number~year+life_stage,data=keogh[keogh$species=="ct",],FUN=sum,na.rm=T)
ct_age <- aggregate(number~year+fresh_age+life_stage,data=keogh[keogh$species=="ct",],FUN=sum,na.rm=T)
CT_annual <- dcast(ct_R,year~life_stage,value.var="number")

ct_lag <- 2 # range is 2-4 Armstrong 1971 TAFS, Trotter 1989 suggests age 2 for estuarine/coastal:  Losee et al. (2018) suggests age 2 for coastal cutties
# Smith 1980 BC FLNRO suggests age 3 for Keogh River
ct_SR <- data.frame("Year"=CT_annual$year[1:(length(CT_annual$year)-ct_lag)],"Stock"=CT_annual$a[1:(length(CT_annual$year)-ct_lag)],"Recruits"=CT_annual$s[(ct_lag+1):length(CT_annual$year)])

plot(ct_SR$Stock,ct_SR$Recruits)
ct_Ricker <- lm(log(ct_SR$Recruits/ct_SR$Stock)~ct_SR$Stock)
curve(exp(coef(ct_Ricker)[1])*x*exp(coef(ct_Ricker)[2]*x),add=T,from=0,to=max(ct_SR$Stock,na.rm=T))

# dolly vardens
dv_R <- aggregate(number~year+life_stage,data=keogh[keogh$species=="dv",],FUN=sum,na.rm=T)
dv_age <- aggregate(number~year+fresh_age+life_stage,data=keogh[keogh$species=="dv",],FUN=sum,na.rm=T)
DV_annual <- dcast(dv_R,year~life_stage,value.var="number")

dv_lag <- 2 # Armstrong 1970 FRBC and Dolloff & Reeves 1990 CJFAS suggest ~ age 1-4 for dolly varden smoltification
dv_SR <- data.frame("Year"=DV_annual$year[1:(length(DV_annual$year)-dv_lag)],"Stock"=DV_annual$a[1:(length(DV_annual$year)-dv_lag)],"Recruits"=DV_annual$s[(dv_lag+1):length(DV_annual$year)])

plot(dv_SR$Stock,dv_SR$Recruits)
dv_Ricker <- lm(log(dv_SR$Recruits/dv_SR$Stock)~dv_SR$Stock)
curve(exp(coef(dv_Ricker)[1])*x*exp(coef(dv_Ricker)[2]*x),add=T,from=0,to=max(dv_SR$Stock,na.rm=T))


# compile data on coho salmon
co_R <- aggregate(number~year+life_stage,data=keogh[keogh$species=="co",],FUN=sum,na.rm=T)
co_age <- aggregate(number~year+fresh_age+life_stage,data=keogh[keogh$species=="co",],FUN=sum,na.rm=T)
co_annual <- dcast(co_R,year~life_stage,value.var="number")

co_S <- rep(NA,length(coho$Year))
names(co_S) <- coho$Year
co_S <- coho$Adults

co_R <- rep(NA,length(coho$Year))
names(co_R) <- coho$Year
co_R[coho$Year%in%coSm$Year] <- coSm$Smolts

# comparison between Coho smolt data from database v. Wade & Irvine report
#plot(cbind(co_R,coSm$Smolts))
#abline(b=1,a=0)

co_lag <- 1 # fixed lag time for coho: 1 year on average from Wade & Irvine 2018 report from Keogh and Holtby et al. 1990 CJFAS paper from Carnation Creek
co_SR <- data.frame("Year"=coho$Year[1:(length(coho$Year)-co_lag)],"Stock"=co_S[1:(length(coho$Year)-co_lag)],"Recruits"=co_R[(co_lag+1):length(coho$Year)])

plot(co_SR$Stock,co_SR$Recruits)
co_Ricker <- lm(log(co_SR$Recruits/co_SR$Stock)~co_SR$Stock)
curve(exp(coef(co_Ricker)[1])*x*exp(coef(co_Ricker)[2]*x),add=T,from=0,to=max(co_SR$Stock,na.rm=T))

# compile data on chum salmon
ch_lag <- 4 # fixed lag time for chum: Neave et al. 1952
ch_SR <- data.frame("Year"=chum$Year[1:(length(chum$Year)-ch_lag)],"Stock"=chum$Adults[1:(length(chum$Year)-ch_lag)],"Recruits"=chum$Adults[(ch_lag+1):length(chum$Year)])

plot(ch_SR$Stock,ch_SR$Recruits)
ch_Ricker <- lm(log(ch_SR$Recruits/ch_SR$Stock)~ch_SR$Stock)
curve(exp(coef(ch_Ricker)[1])*x*exp(coef(ch_Ricker)[2]*x),add=T,from=0,to=max(ch_SR$Stock,na.rm=T))


#tiff("keogh ricker.tiff",compression="lzw",units="in",height=7,width=8,res=800)
layout(matrix(1:6,nrow=3,ncol=2,byrow=T))
par(mar=c(4,4,1,1))
# panel a
plot(stock_rec$Adults,stock_rec$Smolts,xlab="steelhead adults",ylab="steelhead smolts",ylim=c(0,max(stock_rec$Smolts,na.rm=T)))
curve(exp(coef(sh_SR)[1])*x*exp(coef(sh_SR)[2]*x),add=T,from=0,to=max(stock_rec$Adults,na.rm=T))
Corner_text("a)","topleft")

# panel b
plot(pinks$Stock,pinks$Recruits,xlab="pink adults (t)",ylab="pink adults (t+2)")
pk_SR <- lm(log(pinks$Recruits/pinks$Stock)~pinks$Stock)
curve(exp(coef(pk_SR)[1])*x*exp(coef(pk_SR)[2]*x),add=T,from=0,to=max(pinks$Stock,na.rm=T))
Corner_text("b)","topleft")

# panel c
plot(ct_SR$Stock,ct_SR$Recruits,xlab="cutthroat adults (t)",ylab="cutthroat smolts (t+2)")
ct_Ricker <- lm(log(ct_SR$Recruits/ct_SR$Stock)~ct_SR$Stock)
curve(exp(coef(ct_Ricker)[1])*x*exp(coef(ct_Ricker)[2]*x),add=T,from=0,to=max(ct_SR$Stock,na.rm=T))
Corner_text("c)","topleft")

#panel d
plot(dv_SR$Stock,dv_SR$Recruits,xlab="dolly varden adults (t)",ylab="dolly varden smolts (t+4)")
dv_Ricker <- lm(log(dv_SR$Recruits/dv_SR$Stock)~dv_SR$Stock)
curve(exp(coef(dv_Ricker)[1])*x*exp(coef(dv_Ricker)[2]*x),add=T,from=0,to=max(dv_SR$Stock,na.rm=T))
Corner_text("d)","topleft")

#panel e
plot(co_SR$Stock,co_SR$Recruits,xlab="coho adults (t)",ylab="coho smolts (t+2)",xlim=c(0,max(co_SR$Stock[!is.na(co_SR$Recruits)],na.rm=T)))
co_Ricker <- lm(log(co_SR$Recruits/co_SR$Stock)~co_SR$Stock)
curve(exp(coef(co_Ricker)[1])*x*exp(coef(co_Ricker)[2]*x),add=T,from=0,to=max(co_SR$Stock[!is.na(co_SR$Recruits)],na.rm=T))
Corner_text("e)","topleft")

#panel f
plot(ch_SR$Stock,ch_SR$Recruits,xlab="chum adults (t)",ylab="chum adults (t+4)")
ch_Ricker <- lm(log(Recruits/Stock)~Stock,ch_SR[complete.cases(ch_SR),])
curve(exp(coef(ch_Ricker)[1])*x*exp(coef(ch_Ricker)[2]*x),add=T,from=0,to=max(ch_SR$Stock,na.rm=T))
Corner_text("f)","topleft")

#dev.off()

mat1 <- matrix(sapply(1:18,FUN=function(x){rep(x,3)}),nrow=18,ncol=3,byrow=F)
mat1 <- matrix(sapply(1:ncol(mat1),FUN=function(x){rbind(rep(mat1[,x],3))}),nrow=18,ncol=9,byrow=F)
mat1 <- rbind(mat1,rep(0,ncol(mat1)))
mat1 <- rbind(rep(0,ncol(mat1)),mat1)

regimePlot <- function(){c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey90",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA))}

#tiff("keogh time series since 1953.tiff",compression="lzw",units="in",height=9,width=9,res=800)

layout(mat1)
par(mar=c(1,4,1,1))
YrRange <- c(1953,2018)

salmoYrsReg1 <- c(1953,1975)
salmoYrsReg2 <- c(1991,2021)

plot(stock_rec$Year,stock_rec$Adults,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Steelhead adults",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
text("Early",x=1960,y=580)
text("Peak",x=1985,y=580)
text("Collapsed",x=2010,y=580)

title("Spawners",line=1,xpd=NA)

plot(pinks$Year,pinks$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Pink adults (t)",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
plot(ct_SR$Year,ct_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Cutthroat adults (t)",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
plot(dv_SR$Year,dv_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Dolly varden adults (t)",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
plot(co_SR$Year,co_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Coho adults (t)",panel.first=regimePlot())
axis(1,tick=T,labels=TRUE)

plot(ch_SR$Year,ch_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Chum adults (t)",panel.first=regimePlot())
axis(1,tick=T,labels=TRUE)

mtext("Year",side=1,line=2,cex=0.8)

# plot smolts
plot(stock_rec$Year,stock_rec$Smolts,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Steelhead smolts",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
title("Recruits",line=1,xpd=NA)

plot(pinks$Year,pinks$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Pink adults (t+2)",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)

plot(ct_SR$Year,ct_SR$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Cutthroat smolts (t+2)",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)

plot(dv_SR$Year,dv_SR$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Dolly varden smolts (t+4)",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)

plot(co_SR$Year,co_SR$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Coho smolts (t+1)",panel.first=regimePlot())
axis(1,tick=T,labels=TRUE)

plot(ch_SR$Year,ch_SR$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Chum adults (t+4)",panel.first=regimePlot())
axis(1,tick=T,labels=TRUE)

mtext("Year",side=1,line=2,cex=0.8)

# plot smolt production
plot(stock_rec$Year,stock_rec$Smolts/stock_rec$Adults,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Steelhead",panel.first=regimePlot())

axis(1,tick=T,labels=FALSE)
title("Recruits per spawner",line=1,xpd=NA)
plot(pinks$Year,pinks$Recruits/pinks$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Pink",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
plot(ct_SR$Year,ct_SR$Recruits/ct_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Cutthroat",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
plot(dv_SR$Year,dv_SR$Recruits/dv_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Dolly varden",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
plot(co_SR$Year,co_SR$Recruits/co_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Coho",panel.first=regimePlot())
axis(1,tick=T,labels=TRUE)

plot(ch_SR$Year,ch_SR$Recruits/ch_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Chum",panel.first=regimePlot())
axis(1,tick=T,labels=TRUE)

mtext("Year",side=1,line=2,cex=0.8)
#dev.off()

#tiff("keogh time series since 1975.tiff",compression="lzw",units="in",height=9,width=9,res=800)

layout(mat1)
par(mar=c(1,4,1,1))
YrRange <- c(1975,2018)

salmoYrsReg1 <- c(1953,1975)
salmoYrsReg2 <- c(1990,2021)

plot(stock_rec$Year,stock_rec$Adults,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Steelhead adults",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
text("Early",x=1960,y=580)
text("Peak",x=1985,y=580)
text("Collapsed",x=2010,y=580)

title("Spawners",line=1,xpd=NA)

plot(pinks$Year,pinks$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Pink adults (t)",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
plot(ct_SR$Year,ct_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Cutthroat adults (t)",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
plot(dv_SR$Year,dv_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Dolly varden adults (t)",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
plot(co_SR$Year,co_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Coho adults (t)",panel.first=regimePlot())
axis(1,tick=T,labels=TRUE)

plot(ch_SR$Year,ch_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Chum adults (t)",panel.first=regimePlot())
axis(1,tick=T,labels=TRUE)

mtext("Year",side=1,line=2,cex=0.8)

# plot smolts
plot(stock_rec$Year,stock_rec$Smolts,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Steelhead smolts",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
title("Recruits",line=1,xpd=NA)

plot(pinks$Year,pinks$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Pink adults (t+2)",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)

plot(ct_SR$Year,ct_SR$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Cutthroat smolts (t+2)",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)

plot(dv_SR$Year,dv_SR$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Dolly varden smolts (t+4)",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)

plot(co_SR$Year,co_SR$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Coho smolts (t+1)",panel.first=regimePlot())
axis(1,tick=T,labels=TRUE)

plot(ch_SR$Year,ch_SR$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Chum adults (t+4)",panel.first=regimePlot())
axis(1,tick=T,labels=TRUE)

mtext("Year",side=1,line=2,cex=0.8)

# plot smolt production
plot(stock_rec$Year,stock_rec$Smolts/stock_rec$Adults,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Steelhead",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
title("Recruits per spawner",line=1,xpd=NA)
plot(pinks$Year,pinks$Recruits/pinks$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Pink",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
plot(ct_SR$Year,ct_SR$Recruits/ct_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Cutthroat",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
plot(dv_SR$Year,dv_SR$Recruits/dv_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Dolly varden",panel.first=regimePlot())
axis(1,tick=T,labels=FALSE)
plot(co_SR$Year,co_SR$Recruits/co_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Coho",panel.first=regimePlot())
axis(1,tick=T,labels=TRUE)

plot(ch_SR$Year,ch_SR$Recruits/ch_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Chum",panel.first=regimePlot())
axis(1,tick=T,labels=TRUE)

mtext("Year",side=1,line=2,cex=0.8)
#dev.off()


smolts_compare <- merge(annual_sh_abund[,c("Year","sh_s")],keogh_instream[,c("Year","Smolts")],by="Year",all=T)
smolts_compare <- merge(smolts_compare,steelDailyMax,by="Year",all=T)

jpeg("smolt abundance by sample year.jpeg",units="in",res=800,width=7,height=4)
layout(matrix(1:2,ncol=2,byrow=F))
par(mar=c(5,4,1,1))
matplot(smolts_compare$Year,smolts_compare[,-1],type="l",lty=c(1,1,1),lwd=2,col=c("black","dodgerblue","orange"),xlab="Outmigrating year",ylab="Smolt abundance")
legend("topright",c("Current database","Instream","Daily maximum"),bty="n",lwd=2,lty=1,col=c("black","dodgerblue","orange"))


matplot(smolts_compare$Year,sapply(2:4,function(x){(smolts_compare[,x]-smolts_compare[,3])/sd(smolts_compare[,3],na.rm=T)}),type="l",lty=1,lwd=2,col=c("black","dodgerblue","orange"),xlab="Outmigrating year",ylab="Bias (standardized residuals)")
abline(h=0,lty=3,lwd=1,col="red")

dev.off()


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

write.csv(steel_outmigrate,"steelhead_outmigrate.csv",row.names=F)
write.csv(stock_rec,"steelhead_brood_SR.csv",row.names=F)
