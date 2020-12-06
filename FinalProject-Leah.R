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
palette5 <- c("#25CB10", "#5AB60C", "#8FA108",   "#C48C04", "#FA7800")
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

fire_perimeters <- st_read("C:/Users/owner160829a/Desktop/Graduate School/Penn/Courses/Fall 20/MUSA 508/Final Project/Geoprocessing/fire_perimeters.shp") %>%
  filter (YEAR_=="2018" | YEAR_=="2019") %>% st_transform('EPSG:2225')

fire_suppression_facilities <- st_read("C:/Users/owner160829a/Desktop/Graduate School/Penn/Courses/Fall 20/MUSA 508/Final Project/Geoprocessing/fire_suppression_facilities.shp")

fishnet_unclipped <- st_read("C:/Users/owner160829a/Desktop/Graduate School/Penn/Courses/Fall 20/MUSA 508/Final Project/Geoprocessing/Fishnet with Joined Data/fishnet_halfmile_joins.shp") %>%
  st_transform('EPSG:2225')

selected_counties <- st_read("C:/Users/owner160829a/Desktop/Graduate School/Penn/Courses/Fall 20/MUSA 508/Final Project/Geoprocessing/Selected Counties/selected_counties.shp") %>%
  st_transform('EPSG:2225')

fishnet_clipped <- st_intersection(fishnet_unclipped,selected_counties)

fishnet_clipped <- fishnet_clipped %>% dplyr::select(WUI_MAJORI,FVEG_MAJOR,ELEVATION_,
                                                     SLOPE_MEAN,COVER_MAJ,JUL1819_ME,
                                                     AUG1819_ME, SEP1819_ME, OCT1819_ME,
                                                     COUNTY_NAM,COUNTY_ABB,COUNTY_NUM,
                                                     COUNTY_COD, COUNTY_FIP,Shape_Leng,
                                                     Shape_Area, geometry)

# Adding Unique IDs for each cell
fishnet_clipped$ID <-  seq.int(nrow(fishnet_clipped))
                               
# Joining Fire Perimeters to Fishnet

fishnet_fires <- st_intersection(fire_perimeters,fishnet_clipped)
fishnet_fires2 <- fire_perimeters[fishnet_clipped,]
fishnet_fires3 <- st_intersection(fire_perimeters,fishnet_clipped)

## Creating a fishnet grid

#fishnet <- 
# st_make_grid(counties, cellsize = 5280) %>%
#st_sf() %>%
#mutate(uniqueID = rownames(.))

# EXPLORATORY ANALYSIS

# FEATURE ENGINEERING
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

fishnet_clipped <- fishnet_clipped %>% mutate (ElevationBi = if_else(fishnet_clipped$ELEVATION_>3000,"high","low"))

conifer_points <- fishnet_clipped %>% filter(FVEG_MAJOR=="1") %>% st_centroid()

shrub_points<- fishnet_clipped %>% filter(FVEG_MAJOR=="2") %>% st_centroid()

hardwood_points <- fishnet_clipped %>% filter(FVEG_MAJOR=="6") %>% st_centroid()

fishnet_clipped <- fishnet_clipped %>%
  mutate(
    Conifer.nn =
      nn_function(st_coordinates(st_centroid(fishnet_clipped)), st_coordinates(conifer_points),1),
    Shrub.nn=
      nn_function(st_coordinates(st_centroid(fishnet_clipped)), st_coordinates(shrub_points),1),
    Hardwood.nn=
      nn_function(st_coordinates(st_centroid(fishnet_clipped)), st_coordinates(hardwood_points),1),
    Facilities.nn=
      nn_function(st_coordinates(st_centroid(fishnet_clipped)), st_coordinates(fire_suppression_facilities),3))



rowMeans(data[ , c(1,2)], na.rm=TRUE)
# LOCAL MORAN's I
# Join nn features to our fishnet
## important to drop the geometry from joining features
final_net <-
  left_join(crime_net, st_drop_geometry(vars_net), by="uniqueID")

final_net <-
  st_centroid(final_net) %>%
  st_join(dplyr::select(Neighborhoods, NAME), by = "uniqueID") %>%
  st_join(dplyr::select(PoliceDistricts, ID), by = "uniqueID") %>%
  st_drop_geometry() %>%
  left_join(dplyr::select(final_net, geometry, uniqueID)) %>%
  st_sf() %>%
  na.omit()

## Local Moran's I for fishnet grid cells
## generates warnings from PROJ issues
## {spdep} to make polygon to neighborhoods... 
final_net.nb <- poly2nb(as_Spatial(final_net), queen=TRUE)
## ... and neighborhoods to list of weigths
final_net.weights <- nb2listw(final_net.nb, style="W", zero.policy=TRUE)

final_net.localMorans <- 
  cbind(
    as.data.frame(localmoran(final_net$countViolations, final_net.weights)),
    as.data.frame(final_net)) %>% 
  st_sf() %>%
  dplyr::select(Violations_Count = countViolations, 
                Local_Morans_I = Ii, 
                P_Value = `Pr(z > 0)`) %>%
  mutate(Significant_Hotspots = ifelse(P_Value <= 0.05, 1, 0)) %>%
  gather(Variable, Value, -geometry)

vars <- unique(final_net.localMorans$Variable)
varList <- list()

for(i in vars){
  varList[[i]] <- 
    ggplot() +
    geom_sf(data = filter(final_net.localMorans, Variable == i), 
            aes(fill = Value), colour=NA) +
    scale_fill_viridis(name="") +
    labs(title=i) +
    mapTheme() + theme(legend.position="bottom")}

do.call(grid.arrange,c(varList, ncol = 4, top = "Local Morans I statistics, Drug Violations"))

final_net <-
  final_net %>% 
  mutate(drugs.isSig = 
           ifelse(localmoran(final_net$countViolations, 
                             final_net.weights)[,5] <= 0.0000001, 1, 0)) %>%
  mutate(drugs.isSig.dist = 
           nn_function(st_coordinates(st_centroid(final_net)),
                       st_coordinates(st_centroid(
                         filter(final_net, drugs.isSig == 1))), 1))


# DATA VISUALIZATIONS
##continuous variables
house_subsidy %>%
  dplyr::select(y,unemploy_rate, spent_on_repairs, age, campaign, 
                previous,cons.price.idx,cons.conf.idx) %>%
  rename("Unemployment Rate" = unemploy_rate, "$ Spent on Repairs" = spent_on_repairs, "Age of Homeowner"=age, "# of contacts"=campaign, "# of previous contacts"=previous, "Cons. Price Index"=cons.price.idx, "Cons. Conf. Index"=cons.conf.idx) %>%
  gather(Variable, value, -y) %>%
  ggplot(aes(y, value, fill=y)) + 
  geom_bar(position = "dodge", stat = "summary", fun = "mean") + 
  facet_wrap(~Variable, scales = "free") +
  scale_fill_manual(values = palette2) +
  labs(x="y", y="Value", 
       title = "Feature associations with the likelihood of taking credit",
       subtitle = "(continous outcomes)") +
  theme(legend.position = "none")

# CORRELATIONS
numericVars1 <- 
  select_if(house_subsidy, is.numeric) %>% na.omit() %>%
  dplyr::select(age, unemploy_rate, cons.price.idx, cons.conf.idx, inflation_rate, spent_on_repairs,y_numeric)

ggcorrplot(
  round(cor(numericVars1), 1), 
  p.mat = cor_pmat(numericVars1),
  colors = c("#25CB10", "white", "#FA7800"),
  type="lower",
  insig = "blank") +  
  labs(title = "Correlation across Characteristics") 

correlation.long <-
  st_drop_geometry(final_net) %>%
  dplyr::select(-uniqueID, -cvID, -NAME, -ID) %>%
  gather(Variable, Value, -countViolations)

correlation.cor <-
  correlation.long %>%
  group_by(Variable) %>%
  summarize(correlation = cor(Value, countViolations, use = "complete.obs"))

ggplot(filter(correlation.long, Variable=="Abandoned Vehicles"), aes (x=Value, y=countViolations))+  
  geom_point(size = 0.1) +
  geom_text(data = filter(correlation.cor,Variable=="Abandoned Vehicles"), check_overlap=TRUE, aes(label = paste("r = ", round(correlation, 2))),
            x=-Inf, y=Inf, vjust = 1.5, hjust = -.1) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  xlab("Abandoned Vehicles") + ylab("Drug Violations") +
  plotTheme()


# LOGISTIC MODEL
set.seed(3456)
trainIndex <- createDataPartition(house_subsidy$y, p = .65, 
                                  y = paste(house_subsidy$Education_group,house_subsidy$Age_group,house_subsidy$Season,house_subsidy$Employment,house_subsidy$taxLien),
                                  list = FALSE,
                                  times = 1)
housingTrain <- house_subsidy[ trainIndex,]
housingTest  <- house_subsidy[-trainIndex,]

housingModel <- glm(y_numeric ~ .,
                    data=housingTrain %>% 
                      dplyr::select(-y,-X, -education, -age, -Season, -job, -inflation_rate,-Unemployed,-Pdays_group,-Day),
                    family="binomial" (link="logit"))

summary(housingModel)


## Adding Coefficients
x <- housingModel$coefficients
exp(x)


## Fit metrics
pR2(kitchensink)
pR2(housingModel)

## Prediction
testProbs <- data.frame(Outcome = as.factor(housingTest$y_numeric),
                        Probs = predict(housingModel, housingTest, type= "response"))

# Replace NAs with average prob
testProbs$Probs <- ifelse(is.na(testProbs$Probs), 0.1043699, testProbs$Probs) 

testProbskitchensink <- data.frame(Outcome = as.factor(housingTest$y_numeric),
                                   Probs = predict(kitchensink, housingTest, type= "response"))

#Here we want more of a hump in the bottom plot around 1 to indicate that the reg is predictive
ggplot(testProbskitchensink, aes(x = Probs, fill = as.factor(Outcome))) + 
  geom_density() +
  facet_grid(Outcome ~ .) +
  scale_fill_manual(values = palette2) +
  labs(x = "Click", y = "Density of probabilities",
       title = "Distribution of predicted probabilities by observed outcome",
       subtitle = "Kitchen Sink Model") +
  theme(strip.text.x = element_text(size = 18),
        legend.position = "none")

## Confusion matrix
### Might want to change this threshold, here a probability >50% if being predicted as takes credit
testProbskitchensink <- 
  testProbskitchensink %>%
  mutate(predOutcome  = as.factor(ifelse(testProbskitchensink$Probs > 0.5 , 1, 0)))

caret::confusionMatrix(testProbskitchensink$predOutcome, testProbskitchensink$Outcome, 
                       positive = "1")

# ROC Curve
## This us a goodness of fit measure, 1 would be a perfect fit, .5 is a coin toss
auc(testProbs$Outcome, testProbs$Probs)

ggplot(testProbs, aes(d = as.numeric(testProbs$Outcome), m = Probs)) +
  geom_roc(n.cuts = 50, labels = FALSE, colour = "#FE9900") +
  style_roc(theme = theme_grey) +
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = 'grey') +
  labs(title = "ROC Curve - Model with Feature Engineering")