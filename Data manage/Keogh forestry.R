library(sp)
library(sf)

library(raster)
str_name<-'~/Google Drive/SFU postdoc/Keogh river/BC_disturbance/logging_ageclass2012/age_class.tif' 

imported_raster <- raster(str_name)
str(imported_raster)
imported_raster@data
imported_raster@crs
proj4string(imported_raster) <- CRS("+init=EPSG:3005")
layout(1)
plot(imported_raster)

keogh <- data.frame("river"="Keogh","lat"=50.674333,"long"=-127.351641)
coordinates(keogh) <- c("long","lat")
proj4string(keogh) <- CRS("+proj=longlat +datum=WGS84")
spKeogh <- spTransform(keogh, imported_raster@crs)
points(spKeogh,pch=21,cex=2)
log_hist <- extract(imported_raster,             # raster layer
                    spKeogh,   # SPDF with centroids for buffer
                    buffer = 25000,     # buffer size, units depend on CRS
                    fun=hist,         # what to value to extract
                    df=TRUE)         # return a dataframe? 

log_max <- extract(imported_raster,             # raster layer
                   spKeogh,   # SPDF with centroids for buffer
                   buffer = 25000,     # buffer size, units depend on CRS
                   fun=median,         # what to value to extract
                   df=TRUE)         # return a dataframe? 
log_max

str_name<-'~/Google Drive/SFU postdoc/Keogh river/BC_disturbance/logging_ageclass2012/logging_year.tif' 

log_year <- raster(str_name)
str(log_year)
log_year@data
log_year@crs
proj4string(log_year) <- CRS("+init=EPSG:3005")
layout(1)
plot(log_year)
points(spKeogh,pch=21,cex=2,bg="tomato")

log_max <- extract(log_year,             # raster layer
                   spKeogh,   # SPDF with centroids for buffer
                   buffer = 40000,     # buffer size, units depend on CRS
                   fun=unique,         # what to value to extract
                   df=TRUE)         # return a dataframe? 

lu <- spPolygons( rbind(c(-180,-20), c(-160,5), c(-60, 0), 
                        c(-160,-60), c(-180,-20)),rbind(c(80,0), 
                                                        c(100,60), c(120,0), c(120,-55), c(80,0)))
lu <- SpatialPolygonsDataFrame(lu, data.frame(class=c("urban","ag")))
?reclassify
raster::subset(log_year,log_year@data@attributes[[1]]$ID==log_max[1])
ext <- extent(spKeogh)
ext@xmin <- ext@xmin-1.5e4
ext@xmax <- ext@xmax+1.5e4
ext@ymin <- ext@ymin-4e4
ext@ymax <- ext@ymax+1000
new_rast <- crop(log_year,ext)
plot(new_rast)
points(spKeogh,pch=21,cex=2,bg="tomato")


xscale <- (extent(log_year)@xmax-extent(log_year)@xmin)/log_year@ncols # meters per pixel
yscale <- (extent(log_year)@ymax-extent(log_year)@ymin)/log_year@nrows # meters per pixel
years <- 1752:2017

logging_hist <- table(new_rast@data@values)*xscale
sum_tot <- sum(logging_hist)
forestry <- data.frame("Year"=years,"Logging"=as.numeric(logging_hist[match(years,names(logging_hist))]))
forestry$Logging[is.na(forestry$Logging)] <- 0
lag <- 15
forestry$cumul_log <- c(rep(0,lag),sapply((lag+1):nrow(forestry),function(x){sum(forestry$Logging[(x-lag):x])}))

forestry$cumul_footprint <- cumsum(forestry$Logging)

layout(matrix(1:2,ncol=2))
par(mar=c(5,4,1,1))
plot(forestry$Year,forestry$cumul_log,type="l",lwd=2,xlab="Year",ylab="10-year moving window logged area (m2)")
abline(v=1991,lwd=2,lty=2,col="red")
plot(forestry$Year,forestry$cumul_footprint,type="l",lwd=2,xlab="Year",ylab="Cumulative logged area (m2)")
abline(v=1991,lwd=2,lty=2,col="red")

saveRDS(forestry,file="Data/keogh_logging.rds")
