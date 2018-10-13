
keogh <- read.csv("Data/Keogh_Database_Final_01oct18.csv",stringsAsFactors = F,header=T)

keogh <- keogh[!keogh$species_code%in%c("cf","ctt","ctr","ks","shweres","shlgbres","cot"),]
keogh$life_stage[keogh$life_stage=="F"] <- "f"
keogh$life_stage[keogh$life_stage==" "] <- ""

annual <- aggregate(number~year+species+life_stage,data=keogh,FUN=sum,na.rm=T)
colnames(annual) <- c("Year","Species","Stage","Abundance")
annual$Stage <- factor(annual$Stage,levels=c("f","p","s","j","a","k","","u"))
annual$Stage_Num <- as.numeric(annual$Stage)

size <- aggregate(fork_length~year+species+life_stage,data=keogh,FUN=mean,na.rm=T)
variance <- aggregate(fork_length~year+species+life_stage,data=keogh,FUN=sd,na.rm=T)
size$sigma <- variance$fork_length
colnames(size) <- c("Year","Species","Stage","Length","Sigma")


library(reshape2)

annual_abund <- dcast(annual,Year~Species+Stage,value.var="Abundance")
annual_size <- dcast(size,Year~Species+Stage,value.var="Length")
annual_sigma <- dcast(size,Year~Species+Stage,value.var="Sigma")

steel <- subset(annual,annual$Species=="sh")
matplot(steel$Year,steel$Abundance,type="l",lwd=steel$Stage,col=0,xlab="Year",ylab="Abundance")
for(i in unique(steel$Stage_Num)){
  lines(steel$Year[steel$Stage_Num==i],steel$Abundance[steel$Stage_Num==i],lwd=2,col=i)
}

layout(matrix(1:4,nrow=2,byrow=T))
par(mar=c(5,4,1,1))
plot(sh_s~sh_a,data=annual_abund,pch=21,bg="grey50")
lines(smooth.spline(annual_abund$sh_a[!is.na(annual_abund$sh_s) & !is.na(annual_abund$sh_a)],annual_abund$sh_s[!is.na(annual_abund$sh_s) & !is.na(annual_abund$sh_a)],cv=T),lwd=2,lty=2,col="dodgerblue")

plot(dv_s~dv_a,data=annual_abund,pch=21,bg="grey50")
lines(smooth.spline(annual_abund$dv_a[!is.na(annual_abund$dv_s) & !is.na(annual_abund$dv_a)],annual_abund$dv_s[!is.na(annual_abund$dv_s) & !is.na(annual_abund$dv_a)],cv=F),lwd=2,lty=2,col="dodgerblue")

plot(ct_s~ct_a,data=annual_abund,pch=21,bg="grey50")
lines(smooth.spline(annual_abund$ct_a[!is.na(annual_abund$ct_s) & !is.na(annual_abund$ct_a)],annual_abund$ct_s[!is.na(annual_abund$ct_s) & !is.na(annual_abund$ct_a)],cv=F),lwd=2,lty=2,col="dodgerblue")

plot(co_s~co_a,data=annual_abund,pch=21,bg="grey50")
