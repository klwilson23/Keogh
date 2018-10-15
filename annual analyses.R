
keogh <- read.csv("Data/Keogh_Database_Final_01oct18.csv",stringsAsFactors = F,header=T)

keogh <- keogh[!keogh$species_code%in%c("cf","ctt","ctr","ks","shweres","shlgbres","cot"),]
keogh$life_stage[keogh$life_stage=="F"] <- "f"
keogh$life_stage[keogh$life_stage==" "] <- ""
unique(keogh$fresh_age)
keogh$age <- as.numeric(keogh$fresh_age)

annual_age <- aggregate(number~year+species+life_stage+age,data=keogh,FUN=sum,na.rm=T)
colnames(annual_age) <- c("Year","Species","Stage","Age","Abundance")
annual_age$Hatch <- annual_age$Year-annual_age$Age
annual_cohorts <- aggregate(Abundance~Hatch+Species+Stage,data=annual_cohorts,FUN=sum,na.rm=T)


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
annual_cohort <- dcast(annual_cohorts[annual_cohorts$Stage=="s",],Hatch~Species,value.var="Abundance")

annual_abund <- dcast(annual,Year~Species+Stage,value.var="Abundance")
annual_size <- dcast(size,Year~Species+Stage,value.var="Length")
annual_sigma <- dcast(size,Year~Species+Stage,value.var="Sigma")

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

