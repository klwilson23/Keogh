#install.packages("weathercan")
library(weathercan)
library(igraph)
vignette(package = "weathercan")
vignette("weathercan", package = "weathercan")
vignette("glossary", package = "weathercan")

stations[stations$station_name=="PortHardyA",]
stations_search("Port Hardy",interval="day")
stations_search(coords=c(50.67,-127.3),dist=50,interval="day")

small_data <- weather_dl(station_ids = c(201,202, 51319), start = "1975-01-01", end = "1976-01-01",interval = "day")

months <- c(paste(0,1:9,sep=""),10:12)
summer_months <- months[4:10]
winter_months <- months[-(3:10)]
portHardy_pg <- weather_dl(station_ids = c(202, 51319), start = "1950-01-01", end = "2019-01-01",interval = "day")

head(portHardy_pg,10)
colnames(portHardy_pg)
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

climate$rain_range <- climate$total_rain/climate$win_rain

environment <- read.csv("Data/keogh environmental covariates.csv",header=TRUE)

environ_covars <- merge(environment,climate,by="year")
environ_covars <- merge(environment[,!colnames(environment)%in%c("precip_1","precip_2","precip_3","precip_4","temp_2","temp_3","temp_4")],climate,by="year")

lag <- 1
environ_covars$temp_1 <- c(rep(NA,lag),environ_covars$mean_temp[-nrow(environ_covars)])
environ_covars$precip_1 <- c(rep(NA,lag),environ_covars$total_rain[-nrow(environ_covars)])
environ_covars$win_precip_1 <- c(rep(NA,lag),environ_covars$win_rain[-nrow(environ_covars)])
environ_covars$rain_range_1 <- environ_covars$precip_1/environ_covars$win_precip_1


lag <- 2
environ_covars$temp_2 <- c(rep(NA,lag),running_mean(environ_covars$mean_temp,lag)[-((nrow(environ_covars)-(lag-1)):nrow(environ_covars))])
environ_covars$precip_2 <- c(rep(NA,lag),running_mean(environ_covars$total_rain,lag)[-((nrow(environ_covars)-(lag-1)):nrow(environ_covars))])
environ_covars$win_precip_2 <- c(rep(NA,lag),running_mean(environ_covars$win_rain,lag)[-((nrow(environ_covars)-(lag-1)):nrow(environ_covars))])
environ_covars$rain_range_2 <- environ_covars$precip_2/environ_covars$win_precip_2

lag <- 3
environ_covars$temp_3 <- c(rep(NA,lag),running_mean(environ_covars$mean_temp,lag)[-((nrow(environ_covars)-(lag-1)):nrow(environ_covars))])
environ_covars$precip_3 <- c(rep(NA,lag),running_mean(environ_covars$total_rain,lag)[-((nrow(environ_covars)-(lag-1)):nrow(environ_covars))])
environ_covars$win_precip_3 <- c(rep(NA,lag),running_mean(environ_covars$win_rain,lag)[-((nrow(environ_covars)-(lag-1)):nrow(environ_covars))])
environ_covars$rain_range_3 <- environ_covars$precip_3/environ_covars$win_precip_3

plot(rain_range_2~year,data=environ_covars,type="l")
plot(rain_range_3~year,data=environ_covars,type="l")
plot(win_precip_3~year,data=environ_covars,type="l")
plot(rain_range~year,data=environ_covars,type="l")

bakun <- read.csv("Data/Bakun Index.csv",header=TRUE)
bakun_Spr <- do.call(data.frame,aggregate(Bakun_Index_.51N_131W~Year,data=bakun[bakun$Month%in%c(1,2,3,4),],FUN=mean))
colnames(bakun_Spr) <- c("year","bakun")
environ_covars <- merge(environ_covars,bakun_Spr,by="year",all=TRUE)
head(environ_covars)
saveRDS(environ_covars,"environ_covars.rds")

plot(mean_temp~year,data=climate,type="l")
plot(min_mean_temp~year,data=climate,type="l")
plot(max_mean_temp~year,data=climate,type="l")
plot(total_rain~year,data=climate,type="l")
plot(win_mean_temp~year,data=climate,type="l")
plot(win_min_mean_temp~year,data=climate,type="l")
plot(win_max_mean_temp~year,data=climate,type="l")
plot(max_mean_temp~year,data=climate,type="l")
plot(win_rain~year,data=climate,type="l")
plot(win_heat~year,data=climate,type="l")
plot(rain_range~year,data=climate,type="l")
