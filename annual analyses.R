
keogh <- read.csv("Data/Keogh_Database_Final_01oct18.csv",stringsAsFactors = F,header=T)
pinks <- read.csv("Data/Keogh_Pink_Salm.csv",stringsAsFactors = F,header=T)


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

plot(stock_rec$Adults,stock_rec$Smolts,xlab="Adults",ylab="Smolts")
lines(smooth.spline(stock_rec$Adults[!is.na(stock_rec$Smolts)],stock_rec$Smolts[!is.na(stock_rec$Smolts)],cv=T),lwd=2,lty=2,col="dodgerblue")


matplot(stock_rec[,2:3],type="l")

steel <- subset(annual,annual$Species=="sh")
matplot(steel$Year,steel$Abundance,type="l",lwd=steel$Stage,col=0,xlab="Year",ylab="Abundance")
for(i in unique(steel$Stage_Num)){
  lines(steel$Year[steel$Stage_Num==i],steel$Abundance[steel$Stage_Num==i],lwd=2,col=i)
}
lines(annual_cohort$Hatch,annual_cohort$sh,lwd=2,lty=1,col="dodgerblue")

plot(annual_cohort$Hatch,annual_cohort$sh/max(annual_cohort$sh,na.rm=T),xlim=c(1970,2018),ylim=c(0,1),pch=21,bg="grey50",xlab="Year",ylab="Relative abundance")
lines(steel$Year[steel$Stage_Num==3],steel$Abundance[steel$Stage_Num==3]/max(steel$Abundance[steel$Stage_Num==3],na.rm=T))
legend("topright",c("Hatch","Smolts"),lty=c(NA,2),pch=c(21,NA),pt.bg=c("grey50",NA),col="black",bty="n")


layout(matrix(1:4,nrow=2,byrow=T))
par(mar=c(5,4,1,1))
sh_rec <- rep(NA, length(annual_abund$sh_s))

sh_rec[grep(paste(annual_cohort$Hatch,collapse="|",sep=""),annual_abund$Year)] <- annual_cohort$sh[grep(paste(annual_abund$Year,collapse="|",sep=""),annual_cohort$Hatch)]

plot(sh_s~sh_a,annual_abund,pch=21,bg="grey50")
lines(smooth.spline(annual_abund$sh_a[!is.na(annual_abund$sh_s) & !is.na(annual_abund$sh_a)],annual_abund$sh_s[!is.na(annual_abund$sh_s) & !is.na(annual_abund$sh_a)],cv=T),lwd=2,lty=2,col="dodgerblue")

plot(dv_s~dv_a,data=annual_abund,pch=21,bg="grey50")
lines(smooth.spline(annual_abund$dv_a[!is.na(annual_abund$dv_s) & !is.na(annual_abund$dv_a)],annual_abund$dv_s[!is.na(annual_abund$dv_s) & !is.na(annual_abund$dv_a)],cv=F),lwd=2,lty=2,col="dodgerblue")

plot(ct_s~ct_a,data=annual_abund,pch=21,bg="grey50")
lines(smooth.spline(annual_abund$ct_a[!is.na(annual_abund$ct_s) & !is.na(annual_abund$ct_a)],annual_abund$ct_s[!is.na(annual_abund$ct_s) & !is.na(annual_abund$ct_a)],cv=F),lwd=2,lty=2,col="dodgerblue")

plot(co_s~co_a,data=annual_abund,pch=21,bg="grey50")

layout(matrix(1,nrow=1,ncol=1))

cbind(annual_abund$Year,annual_abund$sh_a,sh_rec)
plot(annual_abund$sh_a,sh_rec)
