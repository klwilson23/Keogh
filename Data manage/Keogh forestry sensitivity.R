library(bcdata)
library(dplyr)
library(bcmaps)
library(ggplot2)
library(ggspatial) # this is for adding the arrows
library(rgdal) #use this for data conversion
library(ggrepel) # to offset labels using geom_sf_label_repel  --> Not done here
library(riverdist) # to snap points to River --> Not done here
library(bcmapsdata)
library(sf)
library(sp)

# Set plot box.  Set Port McNeil as centre and box is 25km to each side
plot_area1 <- bc_cities() %>%   #Extract datafram of all bc cities location
  filter(NAME == "Port McNeill") %>%
  #filter(NAME == "Port Hardy") %>%
  st_buffer(dist = 25000)%>%    # 25000 is a decent zoom in one that covers the Keogh
  st_bbox() %>%                 # Turn into a square
  st_as_sfc()                   # converts to sfc object

###--------------    Read in my data ###########
# data1 <- read.csv('sample_location_2018.csv', stringsAsFactors = FALSE) %>% filter(Habitat !="Flow") %>% mutate(UTM_W=as.numeric(UTM_W),UTM_N=as.numeric(UTM_N))

### -- These are the data from the data1 files one using 'dput()'
data1<- structure(list(SiteID = c("1", "2", "3", "4", "5", "6", "7",  "7", "8", "9", "10", "10", "11", "11", "12", "13", "14", "15",  "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26",  "27"), UTM_W = c(616176, 616194, 616231, 616253, 616268, 621447, 621441, 621441, 621464, 621522, 620025, 620025, 620040, 620105, 620146, 620190, 620273, 620277, 625585, 625593, 625585, 625599, 625695, 625712, 629764, 629780, 629787, 629795, 629813, 629835), UTM_N = c(5614158, 5614125, 5614036, 5614014, 5613975, 5608302,    5608290, 5608290, 5608257, 5608138, 5607035, 5607020, 5606991, 5606927, 5606921, 5606916, 5606857, 5606846, 5600222, 5600166, 5600086, 5600079, 5600082, 5600058, 5596856, 5596846, 5596849, 5596853, 5596850, 5596835), Reach = c("W", "W", "W", "W", "W", "W", "W", "W", "W", "W", "X", "X", "X", "X", "X", "X", "X", "X", "Y", "Y", "Y", "Y", "Y", "Y", "Z", "Z", "Z", "Z", "Z", "Z"), Location = c("Pumphouse", "Pumphouse", "Pumphouse", "Pumphouse", "Pumphouse", "Rupert main", "Rupert main", "Rupert main", "Rupert main", "Rupert main", "Highway", "Highway", "Highway",  "Highway", "Highway", "Highway", "Highway", "Highway", "West main",  "West main", "West main", "West main", "West main", "West main",  "Wolfe creek", "Wolfe creek", "Wolfe creek", "Wolfe creek", "Wolfe creek", "Wolfe creek"), Habitat = c("Riffle", "Glide", "Riffle", "Pool", "Riffle", "Run", "Riffle", "Riffle", "Riffle",  "Pool", "Riffle", "Riffle", "Glide", "Glide", "Pool", "Riffle", "Run", "Riffle", "Riffle", "Glide", "Pool", "Riffle", "Run", "Riffle", "Run", "Pool", "Riffle", "Pool", "Glide", "Riffle"), Within_reach_ID = c("3", "1", "2", "1", "1", "1", "2", "2", "1", "1", "3", "3", "1", "1", "1", "2", "1", "1", "3", "1", "1", "2", "1", "1", "1", "2", "2", "1", "1", "1"), Activity = c("active", "active", "active", "active", "active", "active", "active","old", "active", "active", "old", "active", "active", "active","active", "active", "active", "active", "active", "active", "active", "active", "active", "active", "active", "active", "active", "active", "active", "active"), comments = c("", "", "", "", "", "", "2b", "", "", "", "", "3b", "Q-site", "EF", "", "", "", "", "3b", "", "", "", "", "", "", "", "", "", "", "")), class = "data.frame", row.names = c(NA, -30L))

## ---- Convert data: Change coordinates to UTM site dataset to the correct projection.
## ---- My data are UTM 9U and all the mapping layers 'm'.
data2 <- data1  %>% dplyr::select(UTM_W, UTM_N) #%>% mutate(UTM_W=as.numeric(UTM_W),UTM_N=as.numeric(UTM_N)) #
sputm <- SpatialPoints(data2, proj4string=CRS("+proj=utm +zone=9U =+datum=WGS84"))
spgeo <- spTransform(sputm, CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
data_points <- data1 %>% mutate(lat=spgeo$UTM_N,lon=spgeo$UTM_W) %>%
  st_as_sf(coords=c("lon","lat"),crs=3005) %>%  # This chantes the projection
  st_intersection(plot_area1) # I needeed to add this intersection point in the event that points fall outside the plotting region as ggplot will still put them outside the map

## Subset the points data to be a single point for each location (for labelling Locations)
data_point_labels <- data_points %>% group_by(Location) %>% slice(1)

## --- Collect Mapping layers
# below grabs the coastline data within plot box
coast_line <- bc_bound_hres() %>% st_intersection(plot_area1)

## coastlines taken from freshwater atlas
## This looks the same as above but takes a lot longer to run
#coastline_k <- bcdc_query_geodata("freshwater-atlas-coastlines") %>% collect() %>%  st_intersection(square_round_ph)

# Locations of all cities within plot box
city_df<- bc_cities() %>% st_intersection(plot_area1)

#Get a dataset of streams of order 3,4,5 within plot box
#
rivers_in_plot_area <- bcdc_query_geodata("92344413-8035-4c08-b996-65a9b3f62fca") %>%
  filter(STREAM_ORDER %in% c(3,4,5)) %>%  #Defines as only streams order 3,4,5 (too many including 1 &2)
  filter(INTERSECTS(plot_area1)) %>%      # not sure about this line
  collect() %>%                           #Extracts the data
  st_intersection(plot_area1)             #Where it intersects with plot line

# Extract only the keogh river data from the plot box
keogh_r <- rivers_in_plot_area %>%  filter(GNIS_NAME =="Keogh River")
keogh_points <- data.frame(sf::st_coordinates(keogh_r))
#sp::spTransform(keogh_points,CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
# load forestry data
library(raster)
str_name<-'~/Google Drive/SFU postdoc/Keogh river/BC_disturbance/logging_ageclass2012/logging_year.tif' 
log_year <- raster(str_name)
proj4string(log_year) <- CRS("+init=EPSG:3005")
logging_year <- data.frame(log_year@data@attributes[[1]])
log_max <- extract(log_year,             # raster layer
                   keogh_points[1,1:2],   # SPDF with centroids for buffer
                   buffer = 10000,     # buffer size, units depend on CRS
                   fun=table,         # what to value to extract
                   df=TRUE)         # return a dataframe?

lu <- spPolygons( rbind(c(-180,-20), c(-160,5), c(-60, 0), 
                        c(-160,-60), c(-180,-20)),rbind(c(80,0), 
                                                        c(100,60), c(120,0), c(120,-55), c(80,0)))
lu <- SpatialPolygonsDataFrame(lu, data.frame(class=c("urban","ag")))
#raster::subset(log_year,log_year@data@attributes[[1]]$ID==log_max[1])
ext <- extent(keogh_r)
ext@xmin <- ext@xmin-25000
ext@xmax <- ext@xmax+25000
ext@ymin <- ext@ymin-25000
ext@ymax <- ext@ymax+1000
new_rast <- crop(log_year,ext)

xscale <- (extent(log_year)@xmax-extent(log_year)@xmin)/log_year@ncols # meters per pixel
yscale <- (extent(log_year)@ymax-extent(log_year)@ymin)/log_year@nrows # meters per pixel
years <- 1752:2017

logging_hist <- table(new_rast@data@values)*xscale
sum_tot <- sum(logging_hist)
forestry <- data.frame("Year"=years,"Logging"=as.numeric(logging_hist[match(years,names(logging_hist))]))
forestry$Logging[forestry$Year>2012] <- forestry$Logging[forestry$Year==2012]
forestry$Logging[is.na(forestry$Logging)] <- 0
lag <- lag_forest
forestry$cumul_log <- c(rep(0,lag),sapply((lag+1):nrow(forestry),function(x){sum(forestry$Logging[(x-lag):x])}))

forestry$cumul_footprint <- cumsum(forestry$Logging)

layout(matrix(1:2,ncol=2))
par(mar=c(5,4,1,1))
plot(forestry$Year,forestry$cumul_log,type="l",lwd=2,xlab="Year",ylab="10-year moving window logged area (m2)")
abline(v=1991,lwd=2,lty=2,col="red")
plot(forestry$Year,forestry$cumul_footprint,type="l",lwd=2,xlab="Year",ylab="Cumulative logged area (m2)")
abline(v=1991,lwd=2,lty=2,col="red")

saveRDS(forestry,file="Data/keogh_logging.rds")

