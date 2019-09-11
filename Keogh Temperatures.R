
library(weathercan)
vignette(package = "weathercan")

vignette("weathercan", package = "weathercan")
vignette("glossary", package = "weathercan")

head(stations)
stations[stations$station_name=="PortHardyA",]
stations_search("Port Hardy",interval="day")

small_data <- weather_dl(station_ids = c(201,202, 51319), start = "1975-01-01", end = "1976-01-01",interval = "day")

months <- c(paste(0,1:9,sep=""),10:12)
summer_months <- months[4:10]
winter_months <- months[-(3:10)]
portHardy_pg <- weather_dl(station_ids = c(202, 51319), start = "1950-01-01", end = "2019-01-01",interval = "day")
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

environ_covars <- merge(environment[,!colnames(environment)%in%c("precip_1","precip_2","precip_3","precip_4","temp_2","temp_3","temp_4")],climate,by="year")

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
