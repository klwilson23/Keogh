
Corner_text <- function(text, location="topright",...)
{
  legend(location,legend=text, bty ="n", pch=NA,...)
}

keogh <- read.csv("Data/Keogh_Database_Final_01oct18.csv",stringsAsFactors = F,header=T)
pinks <- read.csv("Data/Keogh_Pink_Salm.csv",stringsAsFactors = F,header=T)
pinks <- pinks[order(pinks$Year),]

# read in and calculate pink salmon stock-recruitment
pinks <- data.frame("Year"=pinks$Year,"Stock"=c(pinks$Stock[-1],NA),"Recruits"=pinks$Stock)
plot(pinks$Stock,pinks$Recruits)
#keogh <- keogh[keogh$hatchery!=1|is.na(keogh$hatchery),]

keogh <- keogh[!keogh$species_code%in%c("cf","ctt","ctr","ks","shweres","shlgbres","cot"),]
keogh$life_stage[keogh$life_stage=="F"] <- "f"
keogh$life_stage[keogh$life_stage==" "] <- ""
keogh$age <- as.numeric(keogh$fresh_age)
keogh$age_ocean <- as.numeric(keogh$ocean_age)
keogh$age_ocean[keogh$life_stage=="s"] <- 0
keogh$hatch_year <- keogh$year-(keogh$age-keogh$age_ocean)

annual_age <- aggregate(number~year+hatch_year+species+life_stage+age+age_ocean,data=keogh,FUN=sum,na.rm=T)
colnames(annual_age) <- c("Year","Hatch","Species","Stage","Age","OceanAge","Abundance")
annual_age <- subset(annual_age,Stage%in%c("s","a","k"))
annual_age$Hatch <- annual_age$Year-(annual_age$Age+annual_age$OceanAge)
annual_cohorts <- aggregate(Abundance~Hatch+Year+Species+Age,data=annual_age,FUN=sum,na.rm=T)

steel_age <- subset(annual_age,Species=="sh"&Stage=="s")

annual <- aggregate(number~year+species+life_stage,data=keogh,FUN=sum,na.rm=T)

colnames(annual) <- c("Year","Species","Stage","Abundance")
annual$Stage <- factor(annual$Stage,levels=c("f","p","s","j","a","k","","u"))
annual$Stage_Num <- as.numeric(annual$Stage)

size <- aggregate(fork_length~year+species+life_stage,data=keogh,FUN=mean,na.rm=T)
variance <- aggregate(fork_length~year+species+life_stage,data=keogh,FUN=sd,na.rm=T)
size$sigma <- variance$fork_length
colnames(size) <- c("Year","Species","Stage","Length","Sigma")


library(reshape2)
annual_cohort <- dcast(annual_cohorts,Year~Species+Hatch,value.var="Abundance")

annual_cohort_prop <- annual_cohort[,-1]/rowSums(annual_cohort[,-1],na.rm=T)

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

stock_rec <- data.frame("Year"=annual_sh_abund$Year,"Adults"=annual_sh_abund$sh_a,"Smolts"=annual_smolts_sh)
stock_rec$Smolts[stock_rec$Smolts==0] <- NA

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

ct_lag <- 3 # range is 2-4 Armstrong 1971 TAFS, Trotter 1989 suggests age 2 for estuarine/coastal
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


# coho salmon
co_R <- aggregate(number~year+life_stage,data=keogh[keogh$species=="co",],FUN=sum,na.rm=T)
co_age <- aggregate(number~year+fresh_age+life_stage,data=keogh[keogh$species=="co",],FUN=sum,na.rm=T)
co_annual <- dcast(co_R,year~life_stage,value.var="number")

co_lag <- 1 # Armstrong 1970 FRBC and Dolloff & Reeves 1990 CJFAS suggest ~ age 1-4 for dolly varden smoltification
co_SR <- data.frame("Year"=co_annual$year[1:(length(co_annual$year)-co_lag)],"Stock"=co_annual$a[1:(length(co_annual$year)-co_lag)],"Recruits"=co_annual$s[(co_lag+1):length(co_annual$year)])

plot(co_SR$Stock,co_SR$Recruits)
co_Ricker <- lm(log(co_SR$Recruits/co_SR$Stock)~co_SR$Stock)
curve(exp(coef(co_Ricker)[1])*x*exp(coef(co_Ricker)[2]*x),add=T,from=0,to=max(co_SR$Stock,na.rm=T))


#tiff("keogh ricker.tiff",compression="lzw",units="in",height=7,width=8,res=800)
layout(matrix(1:4,nrow=2,ncol=2,byrow=T))
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

#dev.off()

plot(annual_smolts_sh)
