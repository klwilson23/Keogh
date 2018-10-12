
keogh <- read.csv("Data/Keogh_Database_Final_01oct18.csv",stringsAsFactors = F,header=T)

keogh <- keogh[!keogh$species_code%in%c("cf","ctt","ctr","ks","shweres","shlgbres","cot"),]
keogh$life_stage[keogh$life_stage=="F"] <- "f"
keogh$life_stage[keogh$life_stage==" "] <- ""
unique(keogh$life_stage)
annual <- aggregate(fork_length~year+species+life_stage,data=keogh,FUN=length)
colnames(annual) <- c("Year","Species","Stage","Abundance")
annual$Stage <- as.numeric(factor(annual$Stage,levels=c("f","p","s","j","a","k","","u")))
steel <- subset(annual,annual$Species=="sh")
matplot(steel$Year,steel$Abundance,type="l",lwd=steel$Stage,col=steel$Stage)

plot()