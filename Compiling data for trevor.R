library(reshape2)

keogh <- read.csv("Data/Keogh_Database_Final_01oct18.csv",stringsAsFactors = F,header=T)
#keogh <- read.csv("Data/Keogh_Database_Final_19nov18.csv",stringsAsFactors = F,header=T)
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

# read in and calculate pink salmon stock-recruitment
pinks <- data.frame("Year"=pinks$Year[-((length(pinks$Stock)-pink_lag+1):length(pinks$Stock))],"Stock"=pinks$Stock[-((length(pinks$Stock)-pink_lag+1):length(pinks$Stock))],"Recruits"=pinks$Stock[-(1:pink_lag)])
keogh <- keogh[keogh$hatchery!=1|is.na(keogh$hatchery),]

keogh$life_stage[keogh$life_stage=="F"] <- "f"
keogh$life_stage[keogh$life_stage==" "] <- ""

keogh$total_age <- sapply(keogh$age_final,function(x){as.numeric(strsplit(x,split="\\.")[[1]][1])+sum(as.numeric(strsplit(gsub("[^[:digit:]]","",strsplit(x,split="\\.")[[1]][2]),split="")[[1]]))+nchar(gsub("[[:digit:]]","",strsplit(x,split="\\.")[[1]][2]))})

keogh$age2 <- as.numeric(keogh$fresh_age)
keogh$age_final <- ifelse(!is.na(keogh$total_age),keogh$age_final,keogh$X1_ager)
keogh$age <- sapply(keogh$age_final,function(x){as.numeric(strsplit(x,split="\\.")[[1]][1])})
keogh$OceanAge <- sapply(keogh$age_final,function(x){sum(as.numeric(strsplit(gsub("[^[:digit:]]","",strsplit(x,split="\\.")[[1]][2]),split="")[[1]]))+nchar(gsub("[[:digit:]]","",strsplit(x,split="\\.")[[1]][2]))})

keogh$age_ocean <- as.numeric(keogh$OceanAge)
keogh$age_ocean[keogh$life_stage=="s"] <- 0


keogh$sampAge <- keogh$age+keogh$age_ocean

data <- keogh[keogh$year>=2011 & keogh$species%in%c("sh","co"),]
annual_age <- aggregate(number~year+species+life_stage+sampAge,data=data,FUN=sum,na.rm=T)
colnames(annual_age) <- c("Year","Species","Stage","Age","Abundance")
annual_age <- subset(annual_age,Stage%in%c("s","a","k"))

annual_cohort <- dcast(annual_age,Year~Species+Stage+Age,value.var="Abundance",fun.aggregate=sum)
annual_cohort_sh <- annual_cohort[,grep("sh",colnames(annual_cohort))]
annual_cohort_prop <- annual_cohort_sh/rowSums(annual_cohort_sh,na.rm=T)

annual_cohort_prop <- sapply(c("co_s","sh_s","sh_a","sh_k"),
                             function(x){
                               annual_cohort[,grep(x,colnames(annual_cohort))]/rowSums(annual_cohort[,grep(x,colnames(annual_cohort))],na.rm=T)
                             })
age_prop <- do.call(cbind, annual_cohort_prop)
colnames(age_prop) <- sapply(colnames(age_prop),function(x){paste("Prop.(",strsplit(x,split="\\.")[[1]][2],")",sep="")})

age_prop <- cbind(annual_cohort$Year,age_prop)
colnames(age_prop)[1] <- "Year"

annual_N <- aggregate(number~year+species+life_stage,data=data[data$life_stage%in%c("s","a","k"),],FUN=sum,na.rm=T)
colnames(annual_N) <- c("Year","Species","Stage","Abundance")
annual_N <- dcast(annual_N,Year~Species+Stage,value.var="Abundance",fun.aggregate=sum)

annual_N$co_a[annual_N$Year%in%coho$Year] <- coho$Adults[coho$Year %in% annual_N$Year]

trevor_data <- merge(annual_N,round(age_prop,2),by="Year")

write.csv(trevor_data,"Data/data for trevor.csv",row.names=F)
