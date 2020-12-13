# LIBRARIES
library(rjson)
library(tidycensus)
library(tidyverse)
library(sf)
library(spdep)
library(caret)
library(ckanr)
library(FNN)
library(grid)
library(gridExtra)
library(ggcorrplot)
library(jtools)  
library(viridis)

library(kableExtra)
library(rlist)
library(dplyr)
library(osmdata)
library(geosphere)
library(fastDummies)
library(FNN)
library(viridis)
library(stargazer)
library(pscl)
library(pROC)
library(plotROC)
library(RANN)
library(riem)
options(scipen=999)
options(tigris_class = "sf")

# THEMES AND FUNCTIONS
mapTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(size = 14,colour = "black"),
    plot.subtitle=element_text(face="italic"),
    plot.caption=element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),axis.title = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2)
  )
}

plotTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(size = 14,colour = "black"),
    plot.subtitle = element_text(face="italic"),
    plot.caption = element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_line("grey80", size = 0.1),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    strip.background = element_rect(fill = "grey80", color = "white"),
    strip.text = element_text(size=12),
    axis.title = element_text(size=12),
    axis.text = element_text(size=10),
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(colour = "black", face = "italic"),
    legend.text = element_text(colour = "black", face = "italic"),
    strip.text.x = element_text(size = 14)
  )
}

# CHANGE THESE
palette5 <- c("#25CB10", "#5AB60C", "#8FA108",   "#C48C04", "#FA7800")
palette2 <- c("#981FAC","#FF006A")

qBr <- function(df, variable, rnd) {
  if (missing(rnd)) {
    as.character(quantile(round(df[[variable]],0),
                          c(.01,.2,.4,.6,.8), na.rm=T))
  } else if (rnd == FALSE | rnd == F) {
    as.character(formatC(quantile(df[[variable]]), digits = 3),
                 c(.01,.2,.4,.6,.8), na.rm=T)
  }
}

q5 <- function(variable) {as.factor(ntile(variable, 5))}

# FUNCTIONS
nn_function <- function(measureFrom,measureTo,k) {
  measureFrom_Matrix <- as.matrix(measureFrom)
  measureTo_Matrix <- as.matrix(measureTo)
  nn <-   
    get.knnx(measureTo, measureFrom, k)$nn.dist
  output <-
    as.data.frame(nn) %>%
    rownames_to_column(var = "thisPoint") %>%
    gather(points, point_distance, V1:ncol(.)) %>%
    arrange(as.numeric(thisPoint)) %>%
    group_by(thisPoint) %>%
    summarize(pointDistance = mean(point_distance)) %>%
    arrange(as.numeric(thisPoint)) %>% 
    dplyr::select(-thisPoint) %>%
    pull()
  
  return(output)  
}

#Function MultipleRingBuffer
multipleRingBuffer <- function(inputPolygon, maxDistance, interval) 
{
  #create a list of distances that we'll iterate through to create each ring
  distances <- seq(0, maxDistance, interval)
  #we'll start with the second value in that list - the first is '0'
  distancesCounter <- 2
  #total number of rings we're going to create
  numberOfRings <- floor(maxDistance / interval)
  #a counter of number of rings
  numberOfRingsCounter <- 1
  #initialize an otuput data frame (that is not an sf)
  allRings <- data.frame()
  
  #while number of rings  counteris less than the specified nubmer of rings
  while (numberOfRingsCounter <= numberOfRings) 
  {
    #if we're interested in a negative buffer and this is the first buffer
    #(ie. not distance = '0' in the distances list)
    if(distances[distancesCounter] < 0 & distancesCounter == 2)
    {
      #buffer the input by the first distance
      buffer1 <- st_buffer(inputPolygon, distances[distancesCounter])
      #different that buffer from the input polygon to get the first ring
      buffer1_ <- st_difference(inputPolygon, buffer1)
      #cast this sf as a polygon geometry type
      thisRing <- st_cast(buffer1_, "POLYGON")
      #take the last column which is 'geometry'
      thisRing <- as.data.frame(thisRing[,ncol(thisRing)])
      #add a new field, 'distance' so we know how far the distance is for a give ring
      thisRing$distance <- distances[distancesCounter]
    }
    
    
    #otherwise, if this is the second or more ring (and a negative buffer)
    else if(distances[distancesCounter] < 0 & distancesCounter > 2) 
    {
      #buffer by a specific distance
      buffer1 <- st_buffer(inputPolygon, distances[distancesCounter])
      #create the next smallest buffer
      buffer2 <- st_buffer(inputPolygon, distances[distancesCounter-1])
      #This can then be used to difference out a buffer running from 660 to 1320
      #This works because differencing 1320ft by 660ft = a buffer between 660 & 1320.
      #bc the area after 660ft in buffer2 = NA.
      thisRing <- st_difference(buffer2,buffer1)
      #cast as apolygon
      thisRing <- st_cast(thisRing, "POLYGON")
      #get the last field
      thisRing <- as.data.frame(thisRing$geometry)
      #create the distance field
      thisRing$distance <- distances[distancesCounter]
    }
    
    #Otherwise, if its a positive buffer
    else 
    {
      #Create a positive buffer
      buffer1 <- st_buffer(inputPolygon, distances[distancesCounter])
      #create a positive buffer that is one distance smaller. So if its the first buffer
      #distance, buffer1_ will = 0. 
      buffer1_ <- st_buffer(inputPolygon, distances[distancesCounter-1])
      #difference the two buffers
      thisRing <- st_difference(buffer1,buffer1_)
      #cast as a polygon
      thisRing <- st_cast(thisRing, "POLYGON")
      #geometry column as a data frame
      thisRing <- as.data.frame(thisRing[,ncol(thisRing)])
      #add the distance
      thisRing$distance <- distances[distancesCounter]
    }  
    
    #rbind this ring to the rest of the rings
    allRings <- rbind(allRings, thisRing)
    #iterate the distance counter
    distancesCounter <- distancesCounter + 1
    #iterate the number of rings counter
    numberOfRingsCounter <- numberOfRingsCounter + 1
  }
  
  #convert the allRings data frame to an sf data frame
  allRings <- st_as_sf(allRings)
}

## READ IN DATA

all_counties <- st_read("C:/Users/owner160829a/Desktop/Graduate School/Penn/Courses/Fall 20/MUSA 508/Final Project/Geoprocessing/counties.shp") %>%
  st_transform('EPSG:2225')

fire_pt <- st_read("https://services1.arcgis.com/jUJYIo9tSA7EHvfZ/arcgis/rest/services/California_Fire_Perimeters/FeatureServer/0/query?outFields=*&where=1%3D1&f=geojson")%>%
  st_transform('EPSG:2225')

fire_suppression_facilities <- st_read("C:/Users/owner160829a/Desktop/Graduate School/Penn/Courses/Fall 20/MUSA 508/Final Project/Geoprocessing/fire_suppression_facilities.shp")

fishnet_unclipped <- st_read("C:/Users/owner160829a/Desktop/Graduate School/Penn/Courses/Fall 20/MUSA 508/Final Project/Geoprocessing/Fishnet with Joined Data/fishnet_halfmile_joins.shp") %>%
  st_transform('EPSG:2225')

selected_counties <- st_read("https://opendata.arcgis.com/datasets/a61c138d0a6946da8d1ebb8d1c9db13a_0.geojson") %>%
  filter(COUNTY_NAME == 'Del Norte' | COUNTY_NAME == 'Siskiyou' | COUNTY_NAME == 'Humboldt' | COUNTY_NAME == 'Trinity' |
           COUNTY_NAME == 'Shasta' | COUNTY_NAME == 'Tehama' | COUNTY_NAME == 'Mendocino' | COUNTY_NAME == 'Glenn' |
           COUNTY_NAME == 'Lake' | COUNTY_NAME == 'Colusa' | COUNTY_NAME == 'Sonoma' |COUNTY_NAME == 'Napa' | COUNTY_NAME == 'Yolo') %>%
  st_transform('EPSG:2225')

fishnet_clipped <- st_intersection(fishnet_unclipped,selected_counties)

fishnet_clipped <- fishnet_clipped %>% dplyr::select(WUI_MAJORI,FVEG_MAJOR,ELEVATION_,
                                                     SLOPE_MEAN,COVER_MAJ,JUL1819_ME,
                                                     AUG1819_ME, SEP1819_ME, OCT1819_ME,
                                                     COUNTY_NAME, geometry) %>%
  rename (WUI_MAJ=WUI_MAJORI, 
          FVEG_MAJ=FVEG_MAJOR,
          ELEVATION_AV=ELEVATION_,
          JULY_AVTEMP=JUL1819_ME,
          AUG_AVTEMP=AUG1819_ME,
          SEP_AVTEMP=SEP1819_ME,
          OCT_AVTEMP=OCT1819_ME)

# Replacing NAs with the mean, remove?
fishnet_clipped$JULY_AVTEMP <- ifelse(is.na(fishnet_clipped$JULY_AVTEMP), 15188.16, fishnet_clipped$JULY_AVTEMP) 

fishnet_clipped$AUG_AVTEMP <- ifelse(is.na(fishnet_clipped$AUG_AVTEMP), 15188.16, fishnet_clipped$AUG_AVTEMP) 

fishnet_clipped$SEP_AVTEMP <- ifelse(is.na(fishnet_clipped$SEP_AVTEMP), 15071.67, fishnet_clipped$SEP_AVTEMP) 

fishnet_clipped$OCT_AVTEMP <- ifelse(is.na(fishnet_clipped$OCT_AVTEMP), 14748.27, fishnet_clipped$OCT_AVTEMP) 

# Adding Unique IDs for each cell
fishnet_clipped$ID <-  seq.int(nrow(fishnet_clipped))

## Joining fire data to fishnets 
###2014-18
fire_perimeter1418 <-
  fire_pt %>%
  filter(YEAR_  == '2014' | YEAR_  == '2015' | YEAR_  == '2016' | YEAR_ =='2017' | YEAR_  == '2018') %>%
  st_transform('EPSG:2225')

clip1418 <- 
  st_intersection(st_make_valid(fire_perimeter1418),st_make_valid(fishnet_clipped)) %>%
  select(ID) %>%
  st_drop_geometry() %>%
  mutate(Fire1418 = 1) %>%
  distinct()

fishnet_clipped <-
  fishnet_clipped %>%
  left_join(., clip1418, on= 'ID') 

fishnet_clipped$Fire1418 <- ifelse(is.na(fishnet_clipped$Fire1418),0, fishnet_clipped$Fire1418)

###2019
fire_perimeter19 <-
  fire_pt %>%
  filter(YEAR_ =='2019') %>%
  st_transform('EPSG:2225')

clip19 <- 
  st_intersection(st_make_valid(fire_perimeter19),st_make_valid(fishnet_clipped)) %>%
  select(ID) %>%
  st_drop_geometry() %>%
  mutate(Fire19 = 1) %>%
  distinct()

fishnet_clipped <-
  fishnet_clipped %>%
  left_join(., clip19, on= 'ID') 

fishnet_clipped$Fire19 <- ifelse(is.na(fishnet_clipped$Fire19),0, fishnet_clipped$Fire19)

## WEATHER DATA

# vector 1 - of southern california station ids
weather_station_ids <- c("SIY", "CEC", "MHS", "O86", "ACV", "FOT", "RDD", "RBL", "CIC", "OVE",
                         "UKI", "MYV", "STS", "O69", "DVO", "APC", "SUU", "VCB", "EDU", "SMF", "LHM", "MYV")

# df - stations with lat/lon and name info (in addition to ids)
asos_socal_stations <- riem_stations("CA_ASOS") %>% filter(str_detect(id, paste(weather_station_ids, collapse="|")))
asos_socal_stations$weather_station_id <- asos_socal_stations$id
asos_socal_stations <-  st_as_sf(asos_socal_stations, coords = c("lon","lat"), crs = 4326, agr = "constant") %>% st_transform('EPSG:2225')
asos_socal_stations$weather_ID <-  seq.int(nrow(asos_socal_stations))

## Finding closest station
weather_coords <- 
  asos_socal_stations %>%
  select(geometry)

fishnet_coords <- 
  fishnet_clipped %>%
  select(geometry)

closest_weather_station_to_fishnet <- nn2(weather_coords, fishnet_coords, k = 1)$nn.idx

fishnet_clipped$weather_ID <- closest_weather_station_to_fishnet

ggplot()+
  geom_sf(data = selected_counties)+
  geom_sf(data = weather_coords)

# function and loop -- tried it with the vector of ids and the df of all info, but neither worked. This version below uses just the vector of station ids
get_weather_features_by_station <- function(weather_station_ids, start_year, end_year){
  
  year_vec <- seq(start_year, end_year)
  i <- 1
  weather_data_list <- list()
  for(station_id in weather_station_ids){
    print(paste("Processing station", station_id))
    for(year in year_vec){
      start_date = paste0(year, "07-01")
      end_date = paste0(year, "10-31")
      weather_data <- riem_measures(station = station_id, date_start = start_date, date_end = end_date) %>% 
        dplyr::summarise(weather_station_id = station_id,
                         year = year,
                         Max_Temp = max(tmpf, na.rm = TRUE),
                         Mean_Temp = mean(tmpf, na.rm = TRUE),
                         Mean_Precipitation = mean(p01i, na.rm = TRUE),
                         Mean_Humidity = mean(relh, na.rm = TRUE),
                         Mean_Wind_Speed = mean(sknt, na.rm = TRUE),
        ) 
      weather_data_list[[i]] <- weather_data
      i <- i + 1
    }
  }
  
  do.call("rbind", weather_data_list) 
}

weather_data2014 <- get_weather_features_by_station(weather_station_ids, 2014, 2014) %>%
  rename(Max_Temp14 = Max_Temp,
         Mean_Temp14 = Mean_Temp,
         Mean_Precipitation14 = Mean_Precipitation,
         Mean_Humidity14 = Mean_Humidity,
         Mean_Wind_Speed14 = Mean_Wind_Speed)

weather_2014 <- left_join(weather_data2014, asos_socal_stations, on = 'weather_station_id') %>%
  select (-weather_station_id, -id, -name, -year, -geometry) %>%
  distinct() 

fishnet_clipped <- left_join(fishnet_clipped, weather_2014, on = "weather_ID")

weather_data2015 <- get_weather_features_by_station(weather_station_ids, 2015, 2015) %>%
  rename(Max_Temp15 = Max_Temp,
         Mean_Temp15 = Mean_Temp,
         Mean_Precipitation15 = Mean_Precipitation,
         Mean_Humidity15 = Mean_Humidity,
         Mean_Wind_Speed15 = Mean_Wind_Speed)

weather_2015 <- left_join(weather_data2015, asos_socal_stations, on = 'weather_station_id') %>%
  select (-weather_station_id, -id, -name, -year, -geometry) %>%
  distinct()

fishnet_clipped <- left_join(fishnet_clipped, weather_2015, on = "weather_ID")

weather_data2016 <- get_weather_features_by_station(weather_station_ids, 2016, 2016) %>%
  rename(Max_Temp16 = Max_Temp,
         Mean_Temp16 = Mean_Temp,
         Mean_Precipitation16 = Mean_Precipitation,
         Mean_Humidity16 = Mean_Humidity,
         Mean_Wind_Speed16 = Mean_Wind_Speed)

weather_2016 <- left_join(weather_data2016, asos_socal_stations, on = 'weather_station_id') %>%
  select (-weather_station_id, -id, -name, -year, -geometry)%>%
  distinct() 

fishnet_clipped <- left_join(fishnet_clipped, weather_2016, on = "weather_ID")

weather_data2017 <- get_weather_features_by_station(weather_station_ids, 2017, 2017) %>%
  rename(Max_Temp17 = Max_Temp,
         Mean_Temp17 = Mean_Temp,
         Mean_Precipitation17 = Mean_Precipitation,
         Mean_Humidity17 = Mean_Humidity,
         Mean_Wind_Speed17 = Mean_Wind_Speed)

weather_2017 <- left_join(weather_data2017, asos_socal_stations, on = 'weather_station_id') %>%
  select (-weather_station_id, -id, -name, -year, -geometry)%>%
  distinct() 

fishnet_clipped <- left_join(fishnet_clipped, weather_2017, on = "weather_ID")

weather_data2018 <- get_weather_features_by_station(weather_station_ids, 2018, 2018) %>%
  rename(Max_Temp18 = Max_Temp,
         Mean_Temp18 = Mean_Temp,
         Mean_Precipitation18 = Mean_Precipitation,
         Mean_Humidity18 = Mean_Humidity,
         Mean_Wind_Speed18 = Mean_Wind_Speed)

weather_2018 <- left_join(weather_data2018, asos_socal_stations, on = 'weather_station_id') %>%
  select (-weather_station_id, -id, -name, -year, -geometry)%>%
  distinct() 

fishnet_clipped <- left_join(fishnet_clipped, weather_2018, on = "weather_ID")

weather_data2019 <- get_weather_features_by_station(weather_station_ids, 2019, 2019) %>%
  rename(Max_Temp19 = Max_Temp,
         Mean_Temp19 = Mean_Temp,
         Mean_Precipitation19 = Mean_Precipitation,
         Mean_Humidity19 = Mean_Humidity,
         Mean_Wind_Speed19 = Mean_Wind_Speed)

weather_2019 <- left_join(weather_data2019, asos_socal_stations, on = 'weather_station_id') %>%
  select (-weather_station_id, -id, -name, -year, -geometry)%>%
  distinct() 

fishnet_clipped <- left_join(fishnet_clipped, weather_2019, on = "weather_ID")

# EXPLORATORY ANALYSIS

fire_perimeter1019 <- fire_pt %>% filter(YEAR_=="2010"|YEAR_=="2011"|YEAR_=="2012"|YEAR_=="2013"|YEAR_=="2014"|YEAR_=="2015"|YEAR_=="2016"|YEAR_=="2017"|YEAR_=="2018"|YEAR_=="2019")

ggplot() +
  geom_sf(data = fire_perimeter1019, fill="orange")+
 geom_sf(data=all_counties, fill="transparent")+ 
  labs(title="California Wildfires",
       subtitle="Years 2010-2019")+
  mapTheme()

## Showing our selected counties
ggplot() +
  geom_sf(data = fire_perimeter1019, fill="orange", color="transparent")+
  geom_sf(data=all_counties, fill="transparent")+ 
  geom_sf(data=selected_counties, fill="transparent", color="blue", size=1)+
  labs(title="California Wildfires",
       subtitle="Years 2010-2019; Selected Counties in Blue")+
  mapTheme()

# FEATURE ENGINEERING
# Finding Means of Weather data
fishnet_clipped$MeanTemp1418 = (fishnet_clipped$Mean_Temp14+
                                  fishnet_clipped$Mean_Temp15+
                                  fishnet_clipped$Mean_Temp16+
                                  fishnet_clipped$Mean_Temp17+
                                  fishnet_clipped$Mean_Temp18)/5

fishnet_clipped$MeanHumidity1418 = (fishnet_clipped$Mean_Humidity14+
                                      fishnet_clipped$Mean_Humidity15+
                                      fishnet_clipped$Mean_Humidity16+
                                      fishnet_clipped$Mean_Humidity17+
                                      fishnet_clipped$Mean_Humidity18)/5

fishnet_clipped$MeanPrecip1418 = (fishnet_clipped$Mean_Precipitation14+
                                    fishnet_clipped$Mean_Precipitation15+
                                    fishnet_clipped$Mean_Precipitation16+
                                    fishnet_clipped$Mean_Precipitation17+
                                    fishnet_clipped$Mean_Precipitation18)/5

fishnet_clipped$MeanWindSpeed1418 = (fishnet_clipped$Mean_Wind_Speed14+
                                       fishnet_clipped$Mean_Wind_Speed15+
                                       fishnet_clipped$Mean_Wind_Speed16+
                                       fishnet_clipped$Mean_Wind_Speed17+
                                       fishnet_clipped$Mean_Wind_Speed18)/5

fishnet_clipped$MeanMaxTemp = (fishnet_clipped$Max_Temp14+
                                 fishnet_clipped$Max_Temp15+
                                 fishnet_clipped$Max_Temp16+
                                 fishnet_clipped$Max_Temp17+
                                 fishnet_clipped$Max_Temp18)/5

# Remove Unnecessary Variables
fishnet_clipped <- fishnet_clipped %>% dplyr::select (-JULY_AVTEMP,-AUG_AVTEMP,
                                                      -SEP_AVTEMP,-OCT_AVTEMP,
                                                      -Mean_Temp14,-Mean_Temp15,-Mean_Temp16,-Mean_Temp17,
                                                      -Mean_Temp18,-Mean_Humidity14,-Mean_Humidity15,
                                                      -Mean_Humidity16,-Mean_Humidity17,-Mean_Humidity18,
                                                      -Mean_Precipitation14,-Mean_Precipitation15,
                                                      -Mean_Precipitation16, -Mean_Precipitation17,
                                                      -Mean_Precipitation18,-Mean_Wind_Speed14,
                                                      -Mean_Wind_Speed15,-Mean_Wind_Speed16,
                                                      -Mean_Wind_Speed17,-Mean_Wind_Speed18,
                                                      -Max_Temp14,-Max_Temp15,-Max_Temp16,-Max_Temp17,-Max_Temp18)

## Changing Integers to Characters

fishnet_clipped$WUI_MAJ <- as.factor(fishnet_clipped$WUI_MAJ)

fishnet_clipped$FVEG_MAJ <- as.factor(fishnet_clipped$FVEG_MAJ)

fishnet_clipped$COVER_MAJ <- as.factor(fishnet_clipped$COVER_MAJ)

## Fire in last 3 years
fire_perimeter1013 <-
  fire_pt %>%
  filter(YEAR_  == '2010' | YEAR_ =='2011'| YEAR_ =='2012'| YEAR_ =='2013') %>%
  st_transform('EPSG:2225')

clip1013 <- 
  st_intersection(st_make_valid(fire_perimeter1013),st_make_valid(fishnet_clipped)) %>%
  select(ID) %>%
  st_drop_geometry() %>%
  mutate(Fire1013 = 1) %>%
  distinct()

fishnet_clipped <-
  fishnet_clipped %>%
  left_join(., clip1013, on= 'ID') 

fishnet_clipped$Fire1013 <- ifelse(is.na(fishnet_clipped$Fire1013),0, fishnet_clipped$Fire1013)

## Historical fire
##intersections of fire perimeters with each fishnet cell.
fishnet_clipped <- 
  fishnet_clipped %>% 
  mutate(n_fires_intersections = lengths(st_intersects(st_make_valid(fishnet_clipped), st_make_valid(fire_pt))))

## adding a column for y/n for historical fire presence
#fishnet_clipped <-
 # fishnet_clipped %>%
  #mutate(prev_fire = ifelse(fishnet_clipped$n_fires_intersections > 0, "1", "0"))

## Categorical Features
fishnet_clipped <- fishnet_clipped %>% mutate(CoverCat = case_when(fishnet_clipped$COVER_MAJ=="1"|fishnet_clipped$COVER_MAJ=="2"|fishnet_clipped$COVER_MAJ=="4" ~ "forest",
                                              fishnet_clipped$COVER_MAJ=="6"|fishnet_clipped$COVER_MAJ=="7" ~ "shrubland",
                                              fishnet_clipped$COVER_MAJ=="8"|fishnet_clipped$COVER_MAJ=="9"|fishnet_clipped$COVER_MAJ=="10"~ "savanna_grassland",
                                              fishnet_clipped$COVER_MAJ=="11"|fishnet_clipped$COVER_MAJ=="15"|fishnet_clipped$COVER_MAJ=="17"~ "wet",
                                              fishnet_clipped$COVER_MAJ=="14"|fishnet_clipped$COVER_MAJ=="12" ~ "cropland",
                                              fishnet_clipped$COVER_MAJ=="13" ~ "urban",
                                              fishnet_clipped$COVER_MAJ=="17" ~ "barren"))


fishnet_clipped <- fishnet_clipped %>% mutate(SlopeCat = case_when(fishnet_clipped$SLOPE_MEAN < 5 ~ "low",
                                                                   fishnet_clipped$SLOPE_MEAN >=5|fishnet_clipped$SLOPE_MEAN <15 ~ "medium",
                                                                   fishnet_clipped$SLOPE_MEAN >=15 ~ "high" ))

## Replacing Land Cover type 15 (perm snow and ice) with type 17 (waterbodies)
fishnet_clipped$COVER_MAJ <- ifelse(fishnet_clipped$COVER_MAJ=="15","17",fishnet_clipped$COVER_MAJ)

## Dummy Feature
fishnet_clipped <- fishnet_clipped %>% mutate (ElevationBi = if_else(fishnet_clipped$ELEVATION_AV>3000,1,0))

## Nearest Neighbor Features
conifer_points <- fishnet_clipped %>% filter(FVEG_MAJ=="1") %>% st_centroid()

shrub_points<- fishnet_clipped %>% filter(FVEG_MAJ=="2") %>% st_centroid()

hardwood_points <- fishnet_clipped %>% filter(FVEG_MAJ=="6") %>% st_centroid()

wui_points <- fishnet_clipped %>% filter(WUI_MAJ=="4") %>% st_centroid()

fire1013_points <- fishnet_clipped %>% filter(Fire1013=="1") %>% st_centroid()

fishnet_clipped <- fishnet_clipped %>%
  mutate(
    Conifer.nn =
      nn_function(st_coordinates(st_centroid(fishnet_clipped)), st_coordinates(conifer_points),1),
    Shrub.nn=
      nn_function(st_coordinates(st_centroid(fishnet_clipped)), st_coordinates(shrub_points),1),
    Hardwood.nn=
      nn_function(st_coordinates(st_centroid(fishnet_clipped)), st_coordinates(hardwood_points),1),
    Facilities.nn=
      nn_function(st_coordinates(st_centroid(fishnet_clipped)), st_coordinates(fire_suppression_facilities),3),
    WUI.nn=
      nn_function(st_coordinates(st_centroid(fishnet_clipped)), st_coordinates(wui_points),1),
    Fire.nn=
      nn_function(st_coordinates(st_centroid(fishnet_clipped)), st_coordinates(fire1013_points),5))

# DATA VISUALIZATIONS
##continuous variables
### Need to fix legend here
fishnet_clipped$FIRE <- ifelse(fishnet_clipped$Fire1418==1,"Fire","No_Fire")

fishnet_clipped %>% st_drop_geometry() %>%
  dplyr::select(FIRE, ELEVATION_AV, SLOPE_MEAN,MeanMaxTemp, MeanTemp1418,MeanHumidity1418,
                MeanPrecip1418,MeanWindSpeed1418, n_fires_intersections, Conifer.nn, Shrub.nn, Hardwood.nn, 
                Facilities.nn, WUI.nn) %>%
  rename("Elevation" = ELEVATION_AV, "Slope" = SLOPE_MEAN) %>%
  gather(Variable, value, -FIRE) %>%
  ggplot(aes(FIRE, value, fill=FIRE)) + 
  geom_bar(position = "dodge", stat = "summary", fun = "mean") + 
  facet_wrap(~Variable, scales = "free") +
  #scale_fill_manual(values = palette2) +
  labs(x="Fire", y="Mean", 
       title = "Feature associations with the likelihood of Wildfire",
       subtitle = "(continous outcomes)") +
  theme(legend.position = "none")


# Identifying Colinearity 
numericVars <- select_if(fishnet_clipped, is.numeric) %>% na.omit() %>% st_drop_geometry() %>%
  dplyr::select(ELEVATION_AV,SLOPE_MEAN,MeanTemp1418,MeanHumidity1418,
                MeanPrecip1418,MeanWindSpeed1418,MeanMaxTemp, n_fires_intersections,Conifer.nn, Shrub.nn, Hardwood.nn, 
                Facilities.nn, WUI.nn)

ggcorrplot(
  round(cor(numericVars), 1), 
  p.mat = cor_pmat(numericVars),
  colors = c("#25CB10", "white", "#FA7800"),
  type="lower",
  insig = "blank") +  
  labs(title = "Correlation across Characteristics") 

# LOGISTIC MODEL
set.seed(1214)
#fireTrain1 <- fishnet_clipped %>% filter(Fire1418==1) %>% sample_n()
#fireTrain0 <- fishnet_clipped %>% filter(Fire1418==0) %>% sample_n()
#fireTrain <- rbind(fireTrain1,fireTrain0)
#fireTest <- fishnet_clipped[!(fishnet_clipped$ID %in% fireTrain$ID),]

## Getting warning message that may be affecting levels
trainIndex <- createDataPartition(fishnet_clipped$Fire1418, p = .65, 
                                 y = paste(fishnet_clipped$COVER_MAJ),
                                  list = FALSE,
                                  times = 1)

fireTrain <- fishnet_clipped[ trainIndex,] %>% st_drop_geometry()
fireTest  <- fishnet_clipped[-trainIndex,] %>% st_drop_geometry()

# MODEL
fireModel <- glm(Fire1418 ~ .,
                    data=fireTrain %>% 
                      dplyr::select(-ID,-Fire19,-weather_ID, -WUI_MAJ, -WUI.nn, 
                                    -MeanMaxTemp,-FIRE, -CoverCat, -ElevationBi,
                                    -Fire.nn, -SlopeCat, -MeanPrecip1418, -Max_Temp19, 
                                    -Mean_Temp19, -Mean_Precipitation19, 
                                    -Mean_Humidity19, -Mean_Wind_Speed19, -Fire1013, -Shrub.nn),
                    family="binomial" (link="logit"))

summary(fireModel)

## Adding Coefficients
x <- fireModel$coefficients
exp(x)


## Fit metrics
pR2(fireModel)

## Prediction
testProbs <- data.frame(Outcome = as.factor(fireTest$Fire1418),
                        Probs = predict(fireModel, fireTest, type= "response"))

# Replace NAs with average prob
#testProbs$Probs <- ifelse(is.na(testProbs$Probs), 0.1043699, testProbs$Probs) 

#testProbskitchensink <- data.frame(Outcome = as.factor(housingTest$y_numeric),
                                   #Probs = predict(kitchensink, housingTest, type= "response"))

#Here we want more of a hump in the bottom plot around 1 to indicate that the reg is predictive
ggplot(testProbs, aes(x = Probs, fill = as.factor(Outcome))) + 
  geom_density() +
  facet_grid(Outcome ~ .) +
  scale_fill_manual(values = palette2) +
  labs(x = "Fire", y = "Density of probabilities",
       title = "Distribution of predicted probabilities by observed outcome",
       subtitle = "First Model") +
  theme(strip.text.x = element_text(size = 18),
        legend.position = "none")

## Confusion matrix
### Might want to change this threshold, here a probability >50% if being predicted as takes credit
testProbs <- 
  testProbs %>%
  mutate(predOutcome  = as.factor(ifelse(testProbs$Probs > 0.2 , 1, 0)))

caret::confusionMatrix(testProbs$predOutcome, testProbs$Outcome, 
                       positive = "1")

# ROC Curve
## This us a goodness of fit measure, 1 would be a perfect fit, .5 is a coin toss
auc(testProbs$Outcome, testProbs$Probs)

ggplot(testProbs, aes(d = as.numeric(testProbs$Outcome), m = Probs)) +
  geom_roc(n.cuts = 50, labels = FALSE, colour = "#FE9900") +
  style_roc(theme = theme_grey) +
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = 'grey') +
  labs(title = "ROC Curve - Model with Feature Engineering")

# Testing model on 2019
fireModel19 <- glm(Fire19 ~ .,
                 data=fishnet_clipped %>% st_drop_geometry()%>%
                   dplyr::select(-ID,-weather_ID, -Fire1418, -WUI_MAJ, -WUI.nn, 
                                 -MeanMaxTemp,-FIRE, -CoverCat, -ElevationBi,
                                 -Fire.nn, -SlopeCat, -MeanPrecip1418, -MeanMaxTemp, 
                                 -MeanTemp1418, -Mean_Precipitation19, 
                                 -MeanHumidity1418, -MeanWindSpeed1418, -Fire1013),
                 family="binomial" (link="logit"))

summary(fireModel19)

testProbs19 <- data.frame(Outcome = as.factor(fishnet_clipped$Fire19),
                        Probs = predict(fireModel19, fishnet_clipped, type= "response"))

testProbs19nofire <- testProbs19 %>% filter(Outcome==0)
testProbs19fire <- testProbs19 %>% filter (Outcome==1)

hist(testProbs19nofire$Probs)
hist(testProbs19fire$Probs)

##0.005986742
mean(testProbs19nofire$Probs)

##0.2573636
mean(testProbs19fire$Probs)

## Confusion Matrix for 2019 at .25 threshold
testProbs19 <- 
  testProbs19 %>%
  mutate(predOutcome  = as.factor(ifelse(testProbs19$Probs > 0.25 , 1, 0)))

caret::confusionMatrix(testProbs19$predOutcome, testProbs$Outcome, 
                       positive = "1")

# K Fold Model Validation
ctrl <- trainControl(method = "cv", number = 100, classProbs=TRUE, summaryFunction=twoClassSummary)

cvFit <- train(Fire1418 ~ ., data = fishnet_clipped %>% st_drop_geometry() %>%
                 dplyr::select(
                   -ID,-Fire19,-weather_ID, -WUI_MAJ,
                   -MeanMaxTemp,-FIRE, -CoverCat, -ElevationBi,
                   -Fire.nn, -SlopeCat, -MeanPrecip1418, -Max_Temp19, 
                   -Mean_Temp19, -Mean_Precipitation19, 
                   -Mean_Humidity19, -Mean_Wind_Speed19, -Fire1013)%>%
                 dplyr::mutate(Fire1418=ifelse(Fire1418==1,"Fire","No_Fire")),
               method="glm", family="binomial",
               metric="ROC", trControl = ctrl)

cvFit

dplyr::select(cvFit$resample, -Resample) %>%
  gather(metric, value) %>%
  left_join(gather(cvFit$results[2:4], metric, mean)) %>%
  ggplot(aes(value)) + 
  geom_histogram(bins=35, fill = "#FF006A") +
  facet_wrap(~metric) +
  geom_vline(aes(xintercept = mean), colour = "#981FAC", linetype = 3, size = 1.5) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x="Goodness of Fit", y="Count", title="CV Goodness of Fit Metrics",
       subtitle = "Across-fold mean reprented as dotted lines") +
  plotTheme()


# Model Validation
crossValidate <- function(dataset, id, dependentVariable, indVariables) {
  
  allPredictions <- data.frame()
  cvID_list <- unique(dataset[[id]])
  
  for (i in cvID_list) {
    
    thisFold <- i
    cat("This hold out fold is", thisFold, "\n")
    
    fold.train <- filter(dataset, dataset[[id]] != thisFold) %>% as.data.frame() %>% 
      dplyr::select(id, geometry, indVariables, dependentVariable)
    fold.test  <- filter(dataset, dataset[[id]] == thisFold) %>% as.data.frame() %>% 
      dplyr::select(id, geometry, indVariables, dependentVariable)
    
    regression <-
      glm(Fire1418 ~ ., family = "binomial", 
          data = fold.train %>% 
            dplyr::select(-geometry, -id))
    
    thisPrediction <- 
      mutate(fold.test, Prediction = predict(regression, fold.test, type = "response"))
    
    allPredictions <-
      rbind(allPredictions, thisPrediction)
    
  }
  return(st_sf(allPredictions))
}

reg.vars <- c("FVEG_MAJ","ELEVATION_AV","SLOPE_MEAN", "COVER_MAJ","COUNTY_NAME",
              "MeanTemp1418","MeanHumidity1418","MeanWindSpeed1418",
              "n_fires_intersections","Conifer.nn","Shrub.nn","Hardwood.nn",
              "Facilities.nn","WUI.nn")

reg.spatialCV <- crossValidate(
  dataset = fishnet_clipped,
  id = "COUNTY_NAME",
  dependentVariable = "Fire1418",
  indVariables = reg.vars) 

reg.spatialcv <-
  reg.spatialCV %>%
  dplyr::select(cvID = COUNTY_NAME, Fire1418, Prediction, geometry)

ggplot() +
  geom_sf(data = reg.spatialcv, aes(fill = Prediction), color = "transparent")+
  geom_sf(data = fire_perimeter1418, fill = "transparent", color = "red")

reg.spatialcv <-
  reg.spatialcv %>%
  mutate(Error = Prediction - Fire1418)

ggplot() +
  geom_sf(data = reg.spatialcv, aes(fill = Error), color = "transparent")

# Visualize goodness of fit metrics?

# Cost-Benefit Table (from book)
cost_benefit_table <-
  testProbs %>%
  count(predOutcome, Outcome) %>%
  summarize(True_Negative = sum(n[predOutcome==0 & Outcome==0]),
            True_Positive = sum(n[predOutcome==1 & Outcome==1]),
            False_Negative = sum(n[predOutcome==0 & Outcome==1]),
            False_Positive = sum(n[predOutcome==1 & Outcome==0])) %>%
  gather(Variable, Count) %>%
  mutate(Revenue =
           case_when(Variable == "True_Negative"  ~ Count * 30,
                     Variable == "True_Positive"  ~ ((30 - 8) * (Count * .50)) + 
                       (-32 * (Count * .50)),
                     Variable == "False_Negative" ~ (-30) * Count,
                     Variable == "False_Positive" ~ (30 - 8) * Count)) %>%
  bind_cols(data.frame(Description = c(
    "We predicted no churn and did not send a mailer",
    "We predicted churn and sent the mailer",
    "We predicted no churn and the customer churned",
    "We predicted churn and the customer did not churn")))

# Finding Optimal Threshold

iterateThresholds <- function(data, observedClass, predictedProbs, group) {
  #This function takes as its inputs, a data frame with an observed binomial class (1 or 0); a vector of predicted probabilities; and optionally a group indicator like race. It returns accuracy plus counts and rates of confusion matrix outcomes. It's a bit verbose because of the if (missing(group)). I don't know another way to make an optional parameter.
  observedClass <- enquo(observedClass)
  predictedProbs <- enquo(predictedProbs)
  group <- enquo(group)
  x = .01
  all_prediction <- data.frame()
  
  if (missing(group)) {
    
    while (x <= 1) {
      this_prediction <- data.frame()
      
      this_prediction <-
        data %>%
        mutate(predclass = ifelse(!!predictedProbs > x, 1,0)) %>%
        count(predclass, !!observedClass) %>%
        summarize(Count_TN = sum(n[predclass==0 & !!observedClass==0]),
                  Count_TP = sum(n[predclass==1 & !!observedClass==1]),
                  Count_FN = sum(n[predclass==0 & !!observedClass==1]),
                  Count_FP = sum(n[predclass==1 & !!observedClass==0]),
                  Rate_TP = Count_TP / (Count_TP + Count_FN),
                  Rate_FP = Count_FP / (Count_FP + Count_TN),
                  Rate_FN = Count_FN / (Count_FN + Count_TP),
                  Rate_TN = Count_TN / (Count_TN + Count_FP),
                  Accuracy = (Count_TP + Count_TN) / 
                    (Count_TP + Count_TN + Count_FN + Count_FP)) %>%
        mutate(Threshold = round(x,2))
      
      all_prediction <- rbind(all_prediction,this_prediction)
      x <- x + .01
    }
    return(all_prediction)
  }
  else if (!missing(group)) { 
    while (x <= 1) {
      this_prediction <- data.frame()
      
      this_prediction <-
        data %>%
        mutate(predclass = ifelse(!!predictedProbs > x, 1,0)) %>%
        group_by(!!group) %>%
        count(predclass, !!observedClass) %>%
        summarize(Count_TN = sum(n[predclass==0 & !!observedClass==0]),
                  Count_TP = sum(n[predclass==1 & !!observedClass==1]),
                  Count_FN = sum(n[predclass==0 & !!observedClass==1]),
                  Count_FP = sum(n[predclass==1 & !!observedClass==0]),
                  Rate_TP = Count_TP / (Count_TP + Count_FN),
                  Rate_FP = Count_FP / (Count_FP + Count_TN),
                  Rate_FN = Count_FN / (Count_FN + Count_TP),
                  Rate_TN = Count_TN / (Count_TN + Count_FP),
                  Accuracy = (Count_TP + Count_TN) / 
                    (Count_TP + Count_TN + Count_FN + Count_FP)) %>%
        mutate(Threshold = round(x,2))
      
      all_prediction <- rbind(all_prediction,this_prediction)
      x <- x + .01
    }
    return(all_prediction)
  }
}

whichThreshold <- 
  iterateThresholds(
    data=testProbs, observedClass = Outcome, predictedProbs = Probs)

whichThreshold <- 
  whichThreshold %>%
  dplyr::select(starts_with("Count"), Threshold) %>%
  gather(Variable, Count, -Threshold) %>%
  mutate(Revenue =
           case_when(Variable == "Count_TN"  ~ Count * 30,
                     Variable == "Count_TP"  ~ ((30 - 8) * (Count * .50)) + 
                       (-32 * (Count * .50)),
                     Variable == "Count_FN"  ~ (-30) * Count,
                     Variable == "Count_FP"  ~ (30 - 8) * Count))

whichThreshold %>%
  ggplot(.,aes(Threshold, Revenue, colour = Variable)) +
  geom_point() +
  scale_colour_manual(values = palette5[c(5, 1:3)]) +    
  labs(title = "Revenue by confusion matrix type and threshold",
       y = "Revenue") +
  plotTheme() +
  guides(colour=guide_legend(title = "Confusion Matrix"))

whichThreshold_revenue <- 
  whichThreshold %>% 
  mutate(actualChurn = ifelse(Variable == "Count_TP", (Count * .5),
                              ifelse(Variable == "Count_FN", Count, 0))) %>% 
  group_by(Threshold) %>% 
  summarize(Revenue = sum(Revenue),
            Actual_Churn_Rate = sum(actualChurn) / sum(Count),
            Actual_Churn_Revenue_Loss =  sum(actualChurn * 30),
            Revenue_Next_Period = Revenue - Actual_Churn_Revenue_Loss) 

# Need to provide additional maps and data visualizations to show that the model is useful
