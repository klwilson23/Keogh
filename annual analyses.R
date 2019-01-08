library(reshape2)

Corner_text <- function(text, location="topright",...)
{
  legend(location,legend=text, bty ="n", pch=NA,...)
}

keogh <- read.csv("Data/Keogh_Database_Final_01oct18.csv",stringsAsFactors = F,header=T)
#keogh <- read.csv("Data/Keogh_Database_Final_19nov18.csv",stringsAsFactors = F,header=T)
keogh_atlas <- read.csv("Data/Keogh sh smolts Atlas 2015.csv",stringsAsFactors=F,header=T)
pinks <- read.csv("Data/Keogh_Pink_Salm.csv",stringsAsFactors = F,header=T)
pinks <- pinks[order(pinks$Year),]

coho <- read.csv("Data/Keogh coho adults.csv",stringsAsFactors = F,header=T)
coho <- coho[order(coho$Year),]

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


#keogh[keogh$year==1985 & keogh$species=="sh" & keogh$life_stage=="a",c("age_final","age","ocean_age","OceanAge")]

#as.numeric(strsplit(keogh$age_final[98301],split="\\.")[[1]][2])

#sum(as.numeric(strsplit(gsub("[^[:digit:]]","",strsplit(keogh$X1_ager[185116],split="\\.")[[1]][2]),split="")[[1]]))+nchar(gsub("[[:digit:]]","",strsplit(keogh$X1_ager[185116],split="\\.")[[1]][2]))


#keogh <- keogh[!keogh$species_code%in%c("cf","ctt","ctr","ks","shweres","shlgbres","cot"),]
sampSize <- aggregate(number~year+life_stage,data=keogh[keogh$species=="sh"&!is.na(keogh$age),],FUN=function(x){sum(x)})
ageEstimates <- dcast(sampSize,year~life_stage,value.var="number")

keogh$age_ocean <- as.numeric(keogh$OceanAge)
keogh$age_ocean[keogh$life_stage=="s"] <- 0
keogh$hatch_year <- ifelse(keogh$life_stage=="s",keogh$year - keogh$age, keogh$year-(keogh$total_age))
keogh$smolt_year <- keogh$year-(keogh$age_ocean)

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
annual_cohort_prop <- annual_cohort_sh[,-1]/rowSums(annual_cohort_sh[,-1],na.rm=T)

annual_prop_sh <- annual_cohort_prop[,grep("sh",colnames(annual_cohort_prop))]
row.names(annual_prop_sh) <- row.names(annual_cohort_prop) <- annual_cohort$Year

annual_abund <- dcast(annual,Year~Species+Stage,value.var="Abundance")
annual_size <- dcast(size,Year~Species+Stage,value.var="Length")
annual_sigma <- dcast(size,Year~Species+Stage,value.var="Sigma")

# we want to find how many of the smolts sampled in Year X were from the spawners at Year Y
# to do that, we find the proportions aged and multiply the current year smolts by that proportion from hatch year Y
annual_sh <- subset(annual,Species=="sh")
annual_sh_abund <- dcast(annual_sh,Year~Species+Stage,value.var="Abundance")
annual_smolts_sh <- rep(NA,length(annual_sh_abund$Year))
names(annual_smolts_sh) <- annual_sh_abund$Year
for(i in annual_sh_abund$Year)
{
  annual_smolts_sh[which(annual_sh_abund$Year==i)] <- sum(annual_sh_abund$sh_s[grep(paste("^",row.names(annual_prop_sh)[annual_prop_sh[,grep(i,colnames(annual_prop_sh))]>0],"$",collapse="|",sep=""),annual_sh_abund$Year)]*annual_prop_sh[,grep(i,colnames(annual_prop_sh))][annual_prop_sh[,grep(i,colnames(annual_prop_sh))]>0],na.rm=T)
}

c(0.037,0.585,0.719,0.022,0.0237)*annual_sh_abund$sh_s[3:7]

smolts <- colSums(apply(annual_prop_sh[row.names(annual_prop_sh)%in%annual_sh_abund$Year,],2,FUN=function(x){x*annual_sh_abund$sh_s}))

stock_rec <- data.frame("Year"=annual_sh_abund$Year,"Adults"=annual_sh_abund$sh_a,"Smolts"=annual_smolts_sh)
stock_rec$Smolts[stock_rec$Smolts==0] <- NA

saveRDS(stock_rec,"steelhead_stockRec.rds")

keogh_sh <- readRDS("steelhead_stockRec.rds")
keogh_sh_v2 <- readRDS("steelhead_stockRec_v2.rds")

plot(keogh_atlas$Year,keogh_atlas$Total,ylim=c(0,max(keogh_sh$Smolts,keogh_atlas$Total,na.rm=T)),type="l",xlab="Year",ylab="Smolts (assigned to brood year)")
points(keogh_sh$Year,keogh_sh$Smolts,pch=21,bg="dodgerblue",type="b",lty=2,lwd=1)
#points(keogh_sh_v2$Year,keogh_sh_v2$Smolts,pch=21,bg="orange",ylim=c(0,15000),type="b",lty=1,lwd=1,col="orange")
legend("topright",c("Tom Johnston","Current QA/QC"),col=c("black","dodgerblue"),pch=c(NA,21),pt.bg=c(NA,"dodgerblue"),lty=c(1,2),lwd=c(2,2),bty="n")


plot(stock_rec$Adults,log(stock_rec$Smolts/stock_rec$Adults))
sh_SR <- lm(log(stock_rec$Smolts/stock_rec$Adults)~stock_rec$Adults)
abline(sh_SR)

plot(stock_rec$Adults,stock_rec$Smolts,xlab="Adults",ylab="Smolts",ylim=c(0,max(stock_rec$Smolts,na.rm=T)))
lines(smooth.spline(stock_rec$Adults[!is.na(stock_rec$Smolts)],stock_rec$Smolts[!is.na(stock_rec$Smolts)],cv=T),lwd=2,lty=2,col="dodgerblue")
curve(exp(coef(sh_SR)[1])*x*exp(coef(sh_SR)[2]*x),add=T,from=0,to=max(stock_rec$Adults))

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

co_S <- rep(NA,length(co_annual$year))
names(co_S) <- co_annual$year
co_S[co_annual$year%in%coho$Year] <- coho$Adults

co_lag <- 1 # fixed lag time for coho: 1 year on average from Wade & Irvine 2018 report from Keogh and Holtby et al. 1990 CJFAS paper from Carnation Creek
co_SR <- data.frame("Year"=co_annual$year[1:(length(co_annual$year)-co_lag)],"Stock"=co_S[1:(length(co_annual$year)-co_lag)],"Recruits"=co_annual$s[(co_lag+1):length(co_annual$year)])

plot(co_SR$Stock,co_SR$Recruits)
co_Ricker <- lm(log(co_SR$Recruits/co_SR$Stock)~co_SR$Stock)
curve(exp(coef(co_Ricker)[1])*x*exp(coef(co_Ricker)[2]*x),add=T,from=0,to=max(co_SR$Stock,na.rm=T))


tiff("keogh ricker.tiff",compression="lzw",units="in",height=7,width=8,res=800)
layout(matrix(1:6,nrow=3,ncol=2,byrow=T))
par(mar=c(4,4,1,1))
# panel a
plot(stock_rec$Adults,stock_rec$Smolts,xlab="steelhead adults",ylab="steelhead smolts",ylim=c(0,max(stock_rec$Smolts,na.rm=T)))
curve(exp(coef(sh_SR)[1])*x*exp(coef(sh_SR)[2]*x),add=T,from=0,to=max(stock_rec$Adults))
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
plot(co_SR$Stock,co_SR$Recruits,xlab="coho adults (t)",ylab="coho smolts (t+2)")
co_Ricker <- lm(log(co_SR$Recruits/co_SR$Stock)~co_SR$Stock)
curve(exp(coef(co_Ricker)[1])*x*exp(coef(co_Ricker)[2]*x),add=T,from=0,to=max(co_SR$Stock,na.rm=T))
Corner_text("d)","topleft")
dev.off()

mat1 <- matrix(sapply(1:15,FUN=function(x){rep(x,3)}),nrow=15,ncol=3,byrow=F)
mat1 <- matrix(sapply(1:ncol(mat1),FUN=function(x){rbind(rep(mat1[,x],3))}),nrow=15,ncol=9,byrow=F)
mat1 <- rbind(mat1,rep(0,ncol(mat1)))
mat1 <- rbind(rep(0,ncol(mat1)),mat1)


layout(mat1)
par(mar=c(1,4,1,1))
YrRange <- c(1953,2018)

salmoYrsReg1 <- c(1953,1975)
salmoYrsReg2 <- c(1990,2021)

tiff("keogh time series since 1953.tiff",compression="lzw",units="in",height=7,width=8,res=800)

plot(stock_rec$Year,stock_rec$Adults,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Steelhead adults",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
title("Spawners",line=1,xpd=NA)

plot(pinks$Year,pinks$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Pink adults (t)",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(ct_SR$Year,ct_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Cutthroat adults (t)",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(dv_SR$Year,dv_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Dolly varden adults (t)",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(co_SR$Year,co_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Coho adults (t)",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=TRUE)
mtext("Year",side=1,line=2,cex=0.8)

# plot smolts
plot(stock_rec$Year,stock_rec$Smolts,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Steelhead smolts",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
title("Recruits",line=1,xpd=NA)
plot(pinks$Year,pinks$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Pink adults (t+2)",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(ct_SR$Year,ct_SR$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Cutthroat smolts (t+2)",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(dv_SR$Year,dv_SR$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Dolly varden smolts (t+4)",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(co_SR$Year,co_SR$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Coho smolts (t+1)",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=TRUE)
mtext("Year",side=1,line=2,cex=0.8)

# plot smolt production
plot(stock_rec$Year,stock_rec$Smolts/stock_rec$Adults,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Steelhead",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
title("Recruits per spawner",line=1,xpd=NA)
plot(pinks$Year,pinks$Recruits/pinks$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Pink",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(ct_SR$Year,ct_SR$Recruits/ct_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Cutthroat",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(dv_SR$Year,dv_SR$Recruits/dv_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Dolly varden",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(co_SR$Year,co_SR$Recruits/co_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Coho",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=TRUE)
mtext("Year",side=1,line=2,cex=0.8)
dev.off()

layout(mat1)
par(mar=c(1,4,1,1))
YrRange <- c(1975,2018)

salmoYrsReg1 <- c(1953,1975)
salmoYrsReg2 <- c(1990,2021)

tiff("keogh time series since 1975.tiff",compression="lzw",units="in",height=7,width=8,res=800)

plot(stock_rec$Year,stock_rec$Adults,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Steelhead adults",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
title("Spawners",line=1,xpd=NA)

plot(pinks$Year,pinks$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Pink adults (t)",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(ct_SR$Year,ct_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Cutthroat adults (t)",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(dv_SR$Year,dv_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Dolly varden adults (t)",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(co_SR$Year,co_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Coho adults (t)",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=TRUE)
mtext("Year",side=1,line=2,cex=0.8)

# plot smolts
plot(stock_rec$Year,stock_rec$Smolts,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Steelhead smolts",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
title("Recruits",line=1,xpd=NA)
plot(pinks$Year,pinks$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Pink adults (t+2)",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(ct_SR$Year,ct_SR$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Cutthroat smolts (t+2)",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(dv_SR$Year,dv_SR$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Dolly varden smolts (t+4)",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(co_SR$Year,co_SR$Recruits,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Coho smolts (t+1)",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=TRUE)
mtext("Year",side=1,line=2,cex=0.8)

# plot smolt production
plot(stock_rec$Year,stock_rec$Smolts/stock_rec$Adults,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Steelhead",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
title("Recruits per spawner",line=1,xpd=NA)
plot(pinks$Year,pinks$Recruits/pinks$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Pink",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(ct_SR$Year,ct_SR$Recruits/ct_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Cutthroat",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(dv_SR$Year,dv_SR$Recruits/dv_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Dolly varden",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=FALSE)
plot(co_SR$Year,co_SR$Recruits/co_SR$Stock,xlim=YrRange,type="l",lwd=2,col="dodgerblue",xaxt="n",xlab="",ylab="Coho",panel.first=c(rect(xleft=salmoYrsReg1[1],xright=salmoYrsReg1[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA),rect(xleft=salmoYrsReg2[1],xright=salmoYrsReg2[2],ybottom=-1e6,ytop=1e6,col="grey70",border=NA)))
axis(1,tick=T,labels=TRUE)
mtext("Year",side=1,line=2,cex=0.8)
dev.off()
