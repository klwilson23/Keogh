#install.packages("weathercan")
library(weathercan)
library(igraph)
library(reshape2)
#vignette(package = "weathercan")
#vignette("weathercan", package = "weathercan")
#vignette("glossary", package = "weathercan")

stations_search("Port Hardy",interval="day")
stations_search(coords=c(50.67,-127.3),dist=50,interval="day")

small_data <- weather_dl(station_ids = c(201,202, 51319), start = "1975-01-01", end = "1976-01-01",interval = "day")

months <- c(paste(0,1:9,sep=""),10:12)
summer_months <- months[5:8]
winter_months <- months[-(3:10)]
portHardy_pg <- weather_dl(station_ids = c(202, 51319), start = "1950-01-01", end = "2018-12-31",interval = "day")

portHardy_pg$brood_year <- ifelse(portHardy_pg$month%in%c("01","02"),as.numeric(portHardy_pg$year)-1,as.numeric(portHardy_pg$year))

mean_min_max <- function(x){c("mean"=mean(x),"min"=min(x),"max"=max(x))}
temp <- do.call(data.frame,aggregate(mean_temp~brood_year,data=portHardy_pg[portHardy_pg$month%in%summer_months,],FUN=function(y){c("mean"=mean(y),"min"=min(y),"max"=max(y))}))
colnames(temp) <- c("year","mean_temp","min_mean_temp","max_mean_temp")
rain <- do.call(data.frame,aggregate(total_rain~brood_year,data=portHardy_pg[portHardy_pg$month%in%summer_months,],FUN=sum))
colnames(rain) <- c("year","total_rain")
heating <- do.call(data.frame,aggregate(max_temp~brood_year,data=portHardy_pg[portHardy_pg$month%in%summer_months,],FUN=mean))
colnames(heating) <- c("year","max_temp")

win_temp <- do.call(data.frame,aggregate(mean_temp~brood_year,data=portHardy_pg[portHardy_pg$month%in%winter_months,],FUN=function(y){c("mean"=mean(y),"min"=min(y),"max"=max(y))}))
colnames(win_temp) <- c("year","win_mean_temp","win_min_mean_temp","win_max_mean_temp")
win_rain <- do.call(data.frame,aggregate(total_rain~brood_year,data=portHardy_pg[portHardy_pg$month%in%winter_months,],FUN=sum))
colnames(win_rain) <- c("year","win_rain")
win_heating <- do.call(data.frame,aggregate(max_temp~brood_year,data=portHardy_pg[portHardy_pg$month%in%winter_months,],FUN=mean))
colnames(win_heating) <- c("year","win_heat")

climate <- merge(temp,rain,by="year")
climate <- merge(climate,heating,by="year")
climate <- merge(climate,win_temp,by="year")
climate <- merge(climate,win_rain,by="year")
climate <- merge(climate,win_heating,by="year")
climate$year <- as.numeric(as.character(climate$year))
climate <- subset(climate,climate$year!=2019)

environment <- read.csv("Data/keogh environmental covariates.csv",header=TRUE)

environ_covars <- merge(environment,climate,by="year")
environ_covars <- merge(environment[,!colnames(environment)%in%c("precip_1","precip_2","precip_3","precip_4","temp_2","temp_3","temp_4")],climate,by="year")

bakun <- read.csv("Data/Bakun Index.csv",header=TRUE)
bakun_Spr <- do.call(data.frame,aggregate(Bakun_Index_.51N_131W~Year,data=bakun,FUN=mean))
colnames(bakun_Spr) <- c("year","bakun")
environ_covars <- merge(environ_covars,bakun_Spr,by="year",all=TRUE)
saveRDS(environ_covars,"Data/environ_covars.rds")