# See: https://github.com/bcgov/bcmaps
# see: https://cran.r-project.org/web/packages/bcmaps/bcmaps.pdf
# availble layers in bcmaps:
# https://gist.github.com/ateucher/86674e12ffba66ecce87cceb7bffbf41
# https://github.com/poissonconsulting/fwabc#<https://github.com/poissonconsulting/fwabc>
# https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html

library(bcdata)
library(dplyr)
library(bcmaps)
library(sf)
library(sp)
library(ggplot2)
library(ggspatial) # this is for adding the arrows
library(rgdal) #use this for data conversion
library(ggrepel) # to offset labels using geom_sf_label_repel  --> Not done here
library(riverdist) # to snap points to River --> Not done here
library(bcmapsdata)
library(viridis)
# Set plot box.  Set Port Hardy as centre and box is 20km to each side
plot_area1 <- bc_cities() %>%   #Extract datafram of all bc cities location
  filter(NAME == "Port Hardy") %>%
  #filter(NAME == "Port Hardy") %>%
  st_buffer(dist = 20000)%>%    # 20000 is a decent zoom in one that covers the Keogh
  st_bbox() %>%                 # Turn into a square
  st_as_sfc()                   # converts to sfc object

rivers_in_plot_area <- bcdc_query_geodata("92344413-8035-4c08-b996-65a9b3f62fca") %>%
  filter(STREAM_ORDER %in% c(3,4,5)) %>%  #Defines as only streams order 3,4,5 (too many including 1 &2)
  filter(INTERSECTS(plot_area1)) %>%      # not sure about this line
  collect() %>%                           #Extracts the data
  st_intersection(plot_area1)             #Where it intersects with plot line

# Extract only the keogh river data from the plot box
keogh_r <- rivers_in_plot_area %>%  filter(GNIS_NAME =="Keogh River")
# filter the watershed down to just that of Keogh and tribs: 9208669
keogh_watersh <- rivers_in_plot_area %>%  filter(grepl("9208669",WATERSHED_CODE_50K))
keogh_watersh <- keogh_watersh %>%  filter(grepl(unique(keogh_r$WATERSHED_GROUP_ID),WATERSHED_GROUP_ID))

# filter out the keogh river itself, so we have the tribs as one set to be plotted and the mainstem as another

keogh_watersh <- keogh_watersh %>% filter(LEFT_RIGHT_TRIBUTARY!="NONE")

# Reset plot box.  Set the Keogh Watershed as centre and box is 8.5km to each side
plot_area1 <- keogh_watersh %>%   #Extract datafram of all bc cities location
  #filter(NAME == "Port Hardy") %>%
  #filter(NAME == "Port Hardy") %>%
  st_buffer(dist = 8500)%>%    # 8500 is a decent zoom in one that covers the Keogh
  st_bbox() %>%                 # Turn into a square
  st_as_sfc()                   # converts to sfc object

###--------------    Read in my data ###########
# data1 <- read.csv('sample_location_2018.csv', stringsAsFactors = FALSE) %>% filter(Habitat !="Flow") %>% mutate(UTM_W=as.numeric(UTM_W),UTM_N=as.numeric(UTM_N))

### -- These are the data from the data1 files one using 'dput()'
data1<- structure(list(SiteID = c("1", "2", "3", "4", "5", "6", "7",  "7", "8", "9", "10", "10", "11", "11", "12", "13", "14", "15",  "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26",  "27"), UTM_W = c(616176, 616194, 616231, 616253, 616268, 621447, 621441, 621441, 621464, 621522, 620025, 620025, 620040, 620105, 620146, 620190, 620273, 620277, 625585, 625593, 625585, 625599, 625695, 625712, 629764, 629780, 629787, 629795, 629813, 629835), UTM_N = c(5614158, 5614125, 5614036, 5614014, 5613975, 5608302,    5608290, 5608290, 5608257, 5608138, 5607035, 5607020, 5606991, 5606927, 5606921, 5606916, 5606857, 5606846, 5600222, 5600166, 5600086, 5600079, 5600082, 5600058, 5596856, 5596846, 5596849, 5596853, 5596850, 5596835), Reach = c("W", "W", "W", "W", "W", "W", "W", "W", "W", "W", "X", "X", "X", "X", "X", "X", "X", "X", "Y", "Y", "Y", "Y", "Y", "Y", "Z", "Z", "Z", "Z", "Z", "Z"), Location = c("Pumphouse", "Pumphouse", "Pumphouse", "Pumphouse", "Pumphouse", "Rupert main", "Rupert main", "Rupert main", "Rupert main", "Rupert main", "Highway", "Highway", "Highway",  "Highway", "Highway", "Highway", "Highway", "Highway", "West main",  "West main", "West main", "West main", "West main", "West main",  "Wolfe creek", "Wolfe creek", "Wolfe creek", "Wolfe creek", "Wolfe creek", "Wolfe creek"), Habitat = c("Riffle", "Glide", "Riffle", "Pool", "Riffle", "Run", "Riffle", "Riffle", "Riffle",  "Pool", "Riffle", "Riffle", "Glide", "Glide", "Pool", "Riffle", "Run", "Riffle", "Riffle", "Glide", "Pool", "Riffle", "Run", "Riffle", "Run", "Pool", "Riffle", "Pool", "Glide", "Riffle"), Within_reach_ID = c("3", "1", "2", "1", "1", "1", "2", "2", "1", "1", "3", "3", "1", "1", "1", "2", "1", "1", "3", "1", "1", "2", "1", "1", "1", "2", "2", "1", "1", "1"), Activity = c("active", "active", "active", "active", "active", "active", "active","old", "active", "active", "old", "active", "active", "active","active", "active", "active", "active", "active", "active", "active", "active", "active", "active", "active", "active", "active", "active", "active", "active"), comments = c("", "", "", "", "", "", "2b", "", "", "", "", "3b", "Q-site", "EF", "", "", "", "", "3b", "", "", "", "", "", "", "", "", "", "", "")), class = "data.frame", row.names = c(NA, -30L))

## ---- Convert data: Change coordinates to UTM site dataset to the correct projection.
## ---- My data are UTM 9U and all the mapping layers 'm'.
data2 <- data1  %>% select(UTM_W, UTM_N) #%>% mutate(UTM_W=as.numeric(UTM_W),UTM_N=as.numeric(UTM_N)) #
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

#Below is to get a dataset of lakes in plot box
lakes_in_plot_area <- bcdc_query_geodata("freshwater-atlas-lakes") %>%
  filter(INTERSECTS(plot_area1)) %>%
  collect() %>%
  st_intersection(plot_area1)

# Ocean colouring - this doesn't work well  because the resolution isn't the same
# ocean_colour <- bc_neighbours() %>%
# filter(type=="Ocean") %>%
# st_intersection(plot_area1)
###########

# load forestry data
library(raster)
str_name<-'~/Google Drive/SFU postdoc/Keogh river/BC_disturbance/logging_ageclass2012/logging_year.tif' 
log_year <- raster(str_name)
# project forestry data into map projections
proj4string(log_year) <- CRS("+init=EPSG:3005")
# find the map extents from the raster file, and "crop" or subset the raster, otherwise its a huge file
ext <- extent(plot_area1[[1]][1][[1]][1:4,])
log_year_crop <- crop(log_year,ext)
log_year_crop@data@attributes[[1]]$Rowid <- log_year_crop@data@attributes[[1]]$ID
logging_df <- as.data.frame(log_year_crop,xy=TRUE)
colnames(logging_df) <- c("x", "y","Year","count")
logging_df$Year <- as.numeric(logging_df$Year)
logging_df <- logging_df[complete.cases(logging_df),]
# maniuplate logging data to bin last logged year by 10 year bins
logging_df$Year <- round(logging_df$Year/10)*10
# anything earlier than 1950 was last logged "Pre-1950"
logging_df$Year <- ifelse(logging_df$Year<1950,"Pre-1950",logging_df$Year)
logging_df$Year <- factor(logging_df$Year,levels=c("Pre-1950",sort(unique(logging_df$Year)[!grepl("Pre",unique(logging_df$Year))])))

## The fun part! we get to make a map
## Plot map -- Order of plotting MATTERS
keogh_map <- ggplot() +
              geom_sf(data = coast_line, fill = "grey55") +                 #Plot coastline
              geom_sf(data = plot_area1, alpha = 0,colour='black') +        #Plot area box
              geom_tile(data=logging_df, aes(x=x, y=y, fill=Year), alpha=0.8) +
              #scale_fill_viridis(name = "Last logging year",option="viridis",direction=-1) +
  #scale_fill_gradient2(name="Last logged year",midpoint = 1920, low="darkgreen",mid="olivedrab3",high="goldenrod1",space="Lab") +
              scale_fill_brewer(name="Last logged year",palette = "BrBG",direction=-1) +
              geom_sf(data = rivers_in_plot_area, colour = "lightblue3") +  #Plot Rivers
              geom_sf(data = keogh_watersh, colour = "steelblue2", lwd=1) +
              geom_sf(data = lakes_in_plot_area, fill = "lightblue1") +     #Plot Lakes
              geom_sf(data = keogh_r, colour = "steelblue2",lwd=1.5)+              #Plot Keogh R in RED
              #geom_sf_label(data=data_point_labels[data_point_labels$Location=="Highway",],label=data_point_labels$Location[data_point_labels$Location=="Highway"],size=2)+ #Add labels
              coord_sf(expand = FALSE) +                                    #Expands box to axes
              #geom_sf(data=ocean_colour,fill='red') +                      #Plot Ocean (not working)
              geom_sf_label(data=city_df[city_df$NAME=="Port Hardy",],label=city_df$NAME[city_df$NAME=="Port Hardy"],label.size = 0.15) +            #Labels
              geom_sf_label(data=keogh_r[1,],label=keogh_r$GNIS_NAME[1],nudge_y=-3000,nudge_x=4500,label.size = 0.15) +            #Labels
              geom_sf_label(data=lakes_in_plot_area[grepl("Keogh",lakes_in_plot_area$GNIS_NAME_1),],label=lakes_in_plot_area$GNIS_NAME_1[grepl("Keogh",lakes_in_plot_area$GNIS_NAME_1)],nudge_y=-500,nudge_x=4500,label.size=0.15) + 
              xlab('Longitude') + ylab('Latitude') +                        #Axis titles
              annotation_scale(location = "tr", width_hint = 0.5) +         #Rose Compass
              annotation_north_arrow(location = "tr", which_north = "true",
                                     pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in"),
                                     style = north_arrow_fancy_orienteering,
                                     height = unit(1,"cm"), width = unit(1, "cm"))+
              theme(panel.background = element_rect('lightblue1')           #Make empty space blue to colour ocean
                    , panel.grid.major = element_line('lightblue1'),
                    legend.position="top",legend.text = element_text(size=6),legend.key.width=unit(0.35, "in"))

keogh_map
#Save the plot
ggsave('Figures/keogh_map.pdf',plot=keogh_map,width = 6, height = 6,units='in',dpi=800)
ggsave('Figures/keogh_map.jpeg',plot=keogh_map,width = 6, height = 6,units='in',dpi=800)
