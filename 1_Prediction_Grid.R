# File: 1_Prediction_Grid.R
# Author: Guillermo Martin 
# Template Created: Thu Mar 18 09:56:59 2021
# ---------------------------------------------------------------------------
# Description:
# Creating a grid to predict using the INLA models developed for the 
# standardization of Crab (Cancer pagurus) in the offshore NW fleet of Ireland
# ---------------------------------------------------------------------------
rm(list = ls())


DataDir<-file.path("C:/Users/ggonzales/Desktop/gmartin_work_folder/",
                   "Inshore_Data_compilations/Brown Crab/Data")


# Libraries
library(sf)
library(raster)

# Loading clean and filtered the data -------------------------------------
load(file.path(DataDir,"CRE_offshore_clean_9006_withMask.RData"))


# Log of effort as offset
CRE_data14$Log.E<-log(CRE_data14$Gear_Units)


# Standardize continuous covariates
MyStd <- function(x) { (x - mean(x)) / sd(x)}

CRE_data14[,c("Log.E.std",
              "Soak.std",
              "Gear_Saturation5K3D.std")]<-apply(CRE_data14[,c("Log.E",
                                                               "Soak_Time",
                                                               "Gear_Saturation5K3D")],
                                                 2,MyStd)

# Data preparation --------------------------------------------------------
CRE_data14$Soak_Time<-as.integer(CRE_data14$Soak_Time)
CRE_data14$Target_Units_Landed<-as.numeric(CRE_data14$Target_Units_Landed)

#Categorical variables as factors
CRE_data14$fYear<-factor(CRE_data14$fYear)
CRE_data14$fVESSElID<-as.factor(CRE_data14$fVESSElID)


# Stack for prediction in 1*1 degree --------------------------------------
coords<-data.frame(Lat=round(CRE_data14$Latitude_Start,1),# you may use 0.1*0.1 degree cells or 1*1 degree cells 
                   lon=round(CRE_data14$Longitude_Start,1))
coords$Latlon<-paste0(coords$Lat ,'/', coords$lon)

space = data.frame(Latlon=coords$Latlon[!duplicated(coords$Latlon)])
space = data.frame(Latlon=space[order(space[,1]),])


predDat2 = data.frame(expand.grid(fYear=sort(unique(CRE_data14$fYear)), 
                                  Latlon=as.character(space$Latlon)))


predDat2$Lat =  sapply(predDat2$Latlon, function(x) 
  as.numeric(unlist(strsplit(as.character(x), split='/')))[1])
predDat2$lon =  sapply(predDat2$Latlon, function(x) 
  as.numeric(unlist(strsplit(as.character(x), '/')))[2])

predDat2<-subset(predDat2,select =  -c(Latlon))

# First convert lat/long to UTM
cord.dec<-st_as_sf(predDat2,coords=c("lon","Lat"))
st_crs(cord.dec)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
par(mfrow=c(2,2))
plot(cord.dec[cord.dec$fYear=="2000",])

# Transforming coordinate to UTM 29 (Ireland)
cord.UTM <- st_transform(cord.dec, CRS("+proj=utm +zone=29N +units=km")) # Ireland UTM
cord.UTM
plot(cord.UTM,col="black")

predDat2$XKm<-st_coordinates(cord.UTM)[,1]
predDat2$YKm<-st_coordinates(cord.UTM)[,2]

predDat2<-subset(predDat2,select =  -c(Lat,lon))

#save(predDat2,file=file.path("C:/Users/ggonzales/Desktop/gmartin_work_folder/",
#                             "Inshore_Data_compilations/Brown Crab/Data/",
#                             "Prediction_Grid","Prediction_grid.RData"))
