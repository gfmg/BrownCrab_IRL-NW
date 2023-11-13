# File: CRE_NW_SPiCT_Scenarios.R
# Author: Guillermo Martin 
# Template Created: Thu Apr 29 14:40:14 2021
# ---------------------------------------------------------------------------
# Description:
# Running several different scenario for SPiCT models for the NW Crab data
# ---------------------------------------------------------------------------

rm(list = ls())

library(spict)
library(ggplot2);theme_set(theme_bw())
library(ggpubr)
library(corrplot)


dataDir<-file.path("")

# Load data
load(file.path(dataDir,"1_Landings_DataPrep/2022_Assessment_Landings",
               "CRE_NW_ObsC_new.RData")) #NW CRE landings

#Reconstructed NW CRE landings based on lm. Applied to NI from 1990-2001
load(file.path(dataDir,"1_Landings_DataPrep/2022_Assessment_Landings",
               "CRE_NW_ObsC_Reconstructed.RData")) 

#load(file.path(dataDir,"CRE_NW_ObsI_gam1_Allvessels.RData")) # NW CRE index from gam_1
load(file.path(dataDir,"3_GAM_Standarization/2022_Assessment_Index",
               "CRE_NW_ObsI_gam1_Above8m_21.RData")) # NW CRE index from gam_1
#load(file.path(dataDir,"CRE_NW_ObsI_gam1_quarter_Target.RData")) #NW_CRE_Quarterly based

#Load INLA seasonal index:
load(file.path(dataDir,"3_GAM_Standarization","CRE_NW_I4_Offshore_ar1_index.RData")) #NW_Crab_index from offshore vessels. Model I3

#Additional plots from Paul Boch to check uncertainty
source(file.path("",
                 "production_uncertainty_pb.R"))

#Additional function to compare production curves:
source(file.path("",
                 "extractspict.production.R"))

# Function to plot landings overlaid with the index observations
source(file.path("",
                 "plot_LandingsIndex.R"))
# Data prep ---------------------------------------------------------------

# Sum all landings by countries
#NW_Landings<-subset(NW_Landings,!Year %in% 2021) #Remove 2021 data
CRE_NW_Landings<-aggregate(Landings_ton_CRE~Year,data = NW_Landings,FUN = sum)

# Sum all landings by countries in reconstructed dataset
#NW_LandingsL<-subset(NW_LandingsL,!Year %in% 2020) #Remove 2020 data
CRE_NW_LandingsR<-aggregate(Landings_ton_CRE~Year,data = NW_LandingsL,FUN = sum)

#Add 2020 and 2021 to the reconstructed dataset
CRE_NW_LandingsR<-rbind(CRE_NW_LandingsR,
                        CRE_NW_Landings[CRE_NW_Landings$Year %in% c(2020,2021),])

# Catch series Year to numeric
CRE_NW_Landings$Year<-as.integer(as.character(CRE_NW_Landings$Year))
CRE_NW_LandingsR$Year<-as.integer(as.character(CRE_NW_LandingsR$Year))


# Scale indices of LPUE
Obs_I2$mu.std<-Obs_I2$mu.r/mean(Obs_I2$mu.r)
Obs_I2$ci.low.std<-Obs_I2$ci.low.r/mean(Obs_I2$mu.r)
Obs_I2$ci.up.std<-Obs_I2$ci.up.r/mean(Obs_I2$mu.r)

#Rename index
Obs_I_gam1_Above8m<-Obs_I2

#Year as numeric, reorder data.frame
Obs_I_gam1_Above8m$Year<-as.integer(as.character(Obs_I_gam1_Above8m$fYear))
Obs_I_gam1_Above8m<-Obs_I_gam1_Above8m[order(Obs_I_gam1_Above8m$Year),]

# Quarterly index. 
#Obs_I1_Q$YQ.dec<-ifelse(Obs_I1_Q$fQuarter==1,0,
#                    ifelse(Obs_I1_Q$fQuarter==2,0.25,
#                           ifelse(Obs_I1_Q$fQuarter==3,0.50,
#                                  ifelse(Obs_I1_Q$fQuarter==4,0.75,NA))))

#Obs_I1_Q$YQ<-as.integer(as.character(Obs_I1_Q$fYear))+Obs_I1_Q$YQ.dec

#Scale quarterly index:
#Obs_I1_Q$mu.std<-Obs_I1_Q$mu.r/mean(Obs_I1_Q$mu.r)
#Obs_I1_Q$ci.low.std<-Obs_I1_Q$ci.low.r/mean(Obs_I1_Q$mu.r)
#Obs_I1_Q$ci.up.std<-Obs_I1_Q$ci.up.r/mean(Obs_I1_Q$mu.r)

#Rename index
#Obs_I_gam1_Above8m_Q<-Obs_I1_Q

#Obs_I_gam1_Above8m_Q$Year<-as.integer(as.character(Obs_I_gam1_Above8m_Q$fYear))
#Obs_I_gam1_Above8m_Q<-Obs_I_gam1_Above8m_Q[order(Obs_I_gam1_Above8m_Q$Year),]

# Offshore index year as numeric
I4_index$Year<-as.integer(as.character(I4_index$fYear))

#Loop inputs
End_Year_Landings<-which(CRE_NW_Landings$Year==2021)#Last year to include landings and index data


# SPiCT function ----------------------------------------------------------
fit.SPICT.Scenarios<-function(Scenarios) {
  inp <- list()
  
  if(Scenarios=="F.1"){
    scenario<-"F.1"
    # Description: Landings from 2002 onward. Only index from the SVP fleet removing the 
    # outlier year 2008.
    
    #Landings: 
    inp$obsC <- CRE_NW_Landings$Landings_ton_CRE[23:End_Year_Landings] 
    inp$timeC <- CRE_NW_Landings$Year[23:End_Year_Landings] 
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    

  }
  if(Scenarios=="F.2"){
    scenario<-"F.2"
    
    # Description: Same as F.1 but fixing the surplus production curve to be symmetrical
    
    #Landings: 
    inp$obsC <- CRE_NW_Landings$Landings_ton_CRE[23:End_Year_Landings] 
    inp$timeC <- CRE_NW_Landings$Year[23:End_Year_Landings] 
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    
    #fix the surplus production curve to be symmetrical - optional
    inp$ini$logn <- log(2)
    inp$phases$logn <- -1
    
  }
  if(Scenarios=="F.3"){
    scenario<-"F.3"
    # Description: Same as F.2 but including a prior for the biomass depletion level
    
    #Landings: 
    inp$obsC <- CRE_NW_Landings$Landings_ton_CRE[23:End_Year_Landings] 
    inp$timeC <- CRE_NW_Landings$Year[23:End_Year_Landings] 
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    
    #fix the surplus production curve to be symmetrical - optional
    inp$ini$logn <- log(2)
    inp$phases$logn <- -1
    
    ### a prior for the initial biomass depletion - can help if guessable 
    inp$priors$logbkfrac <- c(log(0.7), 0.5, 1) #Low or no exploitation before beginning of available data
    
  }
  if(Scenarios=="F.4"){
    scenario<-"F.4"
    # Description: Same as F.3 but increasing uncertainty in the first years of landings data
    
    #Landings: 
    inp$obsC <- CRE_NW_Landings$Landings_ton_CRE[23:End_Year_Landings] 
    inp$timeC <- CRE_NW_Landings$Year[23:End_Year_Landings] 
    
    inp$stdevfacC <- rep(1, length(inp$obsC)) 
    inp$stdevfacC[1:5] <- 2 #Extra uncertainty in the first 5 years of landings
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    
    #fix the surplus production curve to be symmetrical - optional
    inp$ini$logn <- log(2)
    inp$phases$logn <- -1
    
    ### a prior for the initial biomass depletion - can help if guessable 
    inp$priors$logbkfrac <- c(log(0.7), 0.5, 1) #Low or no exploitation before beginning of available data
    
  }
  if(Scenarios=="F.5"){
    scenario<-"F.5"
    # Description: Same as F.4 but increasing uncertainty in the time series of landings
    
    #Landings: 
    inp$obsC <- CRE_NW_Landings$Landings_ton_CRE[23:End_Year_Landings] 
    inp$timeC <- CRE_NW_Landings$Year[23:End_Year_Landings] 
    
    inp$stdevfacC <- rep(2, length(inp$obsC)) #Extra uncertainty in the whole time series of landings
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    
    #fix the surplus production curve to be symmetrical - optional
    inp$ini$logn <- log(2)
    inp$phases$logn <- -1
    
    ### a prior for the initial biomass depletion - can help if guessable 
    inp$priors$logbkfrac <- c(log(0.7), 0.5, 1) #Low or no exploitation before beginning of available data
    
  }
  if(Scenarios=="F.6"){
    scenario<-"F.6"
    # Description: Same as F.5 but including the standarized index from 1990-2006 from the offshore fishery
    
    #Landings: 
    inp$obsC <- CRE_NW_Landings$Landings_ton_CRE[23:End_Year_Landings] 
    inp$timeC <- CRE_NW_Landings$Year[23:End_Year_Landings] 
    
    inp$stdevfacC <- rep(2, length(inp$obsC)) #Extra uncertainty in the whole time series of landings
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    inp$timeI[[2]] <- I4_index$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    inp$obsI[[2]] <- I4_index$mu.std
    
    #fix the surplus production curve to be symmetrical - optional
    inp$ini$logn <- log(2)
    inp$phases$logn <- -1
    
    ### a prior for the initial biomass depletion - can help if guessable 
    inp$priors$logbkfrac <- c(log(0.7), 0.5, 1) #Low or no exploitation before beginning of available data
    
  }
  if(Scenarios=="F.7"){
    scenario<-"F.7"
    # Description: Same as F.6 but including the standardized index from 1990-2006 from the offshore fishery
    # only from 2002 onwards
    
    #Landings: 
    inp$obsC <- CRE_NW_Landings$Landings_ton_CRE[23:End_Year_Landings] 
    inp$timeC <- CRE_NW_Landings$Year[23:End_Year_Landings] 
    
    inp$stdevfacC <- rep(2, length(inp$obsC)) #Extra uncertainty in the whole time series of landings
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    inp$timeI[[2]] <- I4_index$Year[12:16]+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    inp$obsI[[2]] <- I4_index$mu.std[12:16]
    
    #fix the surplus production curve to be symmetrical - optional
    inp$ini$logn <- log(2)
    inp$phases$logn <- -1
    
    ### a prior for the initial biomass depletion - can help if guessable 
    inp$priors$logbkfrac <- c(log(0.7), 0.5, 1) #Low or no exploitation before beginning of available data
    
  }
  if(Scenarios=="R.8"){
    scenario<-"R.8"
    # Description: using the reconstructed landings in NI from 1990 and only the SVP index. 
    # Remaining parameters as F.1
    
    #Landings: 
    inp$obsC <- CRE_NW_LandingsR$Landings_ton_CRE[11:End_Year_Landings] 
    inp$timeC <- CRE_NW_LandingsR$Year[11:End_Year_Landings] 
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    
    
  }
  if(Scenarios=="R.9"){
    scenario<-"R.9"
    # Description: using the reconstructed landings in NI from 1990 and both the SVP and
    # offshore vivier index. Remaining parameters as F.1
    
    #Landings: 
    inp$obsC <- CRE_NW_LandingsR$Landings_ton_CRE[11:End_Year_Landings] 
    inp$timeC <- CRE_NW_LandingsR$Year[11:End_Year_Landings] 
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    inp$timeI[[2]] <- I4_index$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    inp$obsI[[2]] <- I4_index$mu.std
    
  }
  if(Scenarios=="R.10"){
    scenario<-"R.10"
    # Description: Same as R.9 with extra uncertainty in the landings for the first years
    
    #Landings: 
    inp$obsC <- CRE_NW_LandingsR$Landings_ton_CRE[11:End_Year_Landings] 
    inp$timeC <- CRE_NW_LandingsR$Year[11:End_Year_Landings] 
    
    
    inp$stdevfacC <- rep(1, length(inp$obsC)) #Extra uncertainty in the whole time series of landings
    inp$stdevfacC[1:13] <- 2 #Extra uncertainty from 1990-2002
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    inp$timeI[[2]] <- I4_index$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    inp$obsI[[2]] <- I4_index$mu.std
    

    
  }
  if(Scenarios=="R.11"){
    scenario<-"R.11"
    # Description: Same as R.10 with extra uncertainty in the landings, but fixing production curve
    # and prior for depletion level. 
    
    #Landings: 
    inp$obsC <- CRE_NW_LandingsR$Landings_ton_CRE[11:End_Year_Landings] 
    inp$timeC <- CRE_NW_LandingsR$Year[11:End_Year_Landings] 
    
    
    inp$stdevfacC <- rep(1, length(inp$obsC)) #Extra uncertainty in the whole time series of landings
    inp$stdevfacC[1:13] <- 2 #Extra uncertainty from 1990-2002
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    inp$timeI[[2]] <- I4_index$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    inp$obsI[[2]] <- I4_index$mu.std
    
    #fix the surplus production curve to be symmetrical - optional
    inp$ini$logn <- log(2)
    inp$phases$logn <- -1
    
    ### a prior for the initial biomass depletion - can help if guessable 
    inp$priors$logbkfrac <- c(log(0.9), 0.5, 1) #Low or no exploitation before beginning of available data
    
    
    
  }
  if(Scenarios=="R.12"){
    scenario<-"R.12"
    # Description: Same as R.11 with extra uncertainty in the landings for all time series
    
    #Landings: 
    inp$obsC <- CRE_NW_LandingsR$Landings_ton_CRE[11:End_Year_Landings] 
    inp$timeC <- CRE_NW_LandingsR$Year[11:End_Year_Landings] 
    
    
    inp$stdevfacC <- rep(2, length(inp$obsC)) #Extra uncertainty in the whole time series of landings
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    inp$timeI[[2]] <- I4_index$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    inp$obsI[[2]] <- I4_index$mu.std
    
    #fix the surplus production curve to be symmetrical - optional
    inp$ini$logn <- log(2)
    inp$phases$logn <- -1
    
    ### a prior for the initial biomass depletion - can help if guessable 
    inp$priors$logbkfrac <- c(log(0.9), 0.5, 1) #Low or no exploitation before beginning of available data
    
    
  }
  if(Scenarios=="R.13"){
    scenario<-"R.13"
    # Description: Same as R.12 with extra uncertainty in the landings for all time series and
    #at the beggining of the time series
    
    #Landings: 
    inp$obsC <- CRE_NW_LandingsR$Landings_ton_CRE[11:End_Year_Landings] 
    inp$timeC <- CRE_NW_LandingsR$Year[11:End_Year_Landings] 
    
    
    inp$stdevfacC <- rep(2, length(inp$obsC)) #Extra uncertainty in the whole time series of landings
    inp$stdevfacC[1:13] <- 4
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    inp$timeI[[2]] <- I4_index$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    inp$obsI[[2]] <- I4_index$mu.std
    
    #fix the surplus production curve to be symmetrical - optional
    inp$ini$logn <- log(2)
    inp$phases$logn <- -1
    
    ### a prior for the initial biomass depletion - can help if guessable 
    inp$priors$logbkfrac <- c(log(0.9), 0.5, 1) #Low or no exploitation before beginning of available data
    
    #Management interval from 2021-2022
    inp$maninterval<-c(CRE_NW_Landings$Year[End_Year_Landings]+1,
                       CRE_NW_Landings$Year[End_Year_Landings]+2)
  }
  if(Scenarios=="G.1"){
    scenario<-"G.1"
    # Description: Same as R.13 but including a very informative prior for the intrinsic
    #growth rate
    
    #Landings: 
    inp$obsC <- CRE_NW_LandingsR$Landings_ton_CRE[11:End_Year_Landings] 
    inp$timeC <- CRE_NW_LandingsR$Year[11:End_Year_Landings] 
    
    
    inp$stdevfacC <- rep(2, length(inp$obsC)) #Extra uncertainty in the whole time series of landings
    inp$stdevfacC[1:13] <- 4
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    inp$timeI[[2]] <- I4_index$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    inp$obsI[[2]] <- I4_index$mu.std
    
    #fix the surplus production curve to be symmetrical - optional
    inp$ini$logn <- log(2)
    inp$phases$logn <- -1
    
    ### a prior for the initial biomass depletion - can help if guessable 
    inp$priors$logbkfrac <- c(log(0.9), 0.5, 1) #Low or no exploitation before beginning of available data
  
    # Prior for the intrinsic growth rate
    inp$priors$logr <- c(log(0.3), 0.2)
  }
  
  if(Scenarios=="G.2"){
    scenario<-"G.2"
    # Description: Same as G.1 but relaxing the prior for the growth rate.
    
    #Landings: 
    inp$obsC <- CRE_NW_LandingsR$Landings_ton_CRE[11:End_Year_Landings] 
    inp$timeC <- CRE_NW_LandingsR$Year[11:End_Year_Landings] 
    
    
    inp$stdevfacC <- rep(2, length(inp$obsC)) #Extra uncertainty in the whole time series of landings
    inp$stdevfacC[1:13] <- 4
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    inp$timeI[[2]] <- I4_index$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    inp$obsI[[2]] <- I4_index$mu.std
    
    #fix the surplus production curve to be symmetrical - optional
    inp$ini$logn <- log(2)
    inp$phases$logn <- -1
    
    ### a prior for the initial biomass depletion - can help if guessable 
    inp$priors$logbkfrac <- c(log(0.9), 0.5, 1) #Low or no exploitation before beginning of available data
    
    # Prior for the intrinsic growth rate
    inp$priors$logr <- c(log(0.3), 0.5)
    
    #Management interval from 2021-2022
    #inp$maninterval<-c(CRE_NW_Landings$Year[End_Year_Landings]+1,
    #                   CRE_NW_Landings$Year[End_Year_Landings]+2)
  }
  
  res <- fit.spict(inp)
  
  return(list(scenario, inp, res))
}


Scenarios.F<-c(paste0("F.",seq(1,7,by=1)))#paste0("F.",seq(1,9,by=1))
Scenario.R<-c(paste0("R.",seq(8,13,by=1)))
Scenario.G<-c(paste0("G.",seq(1,2,by=1)))
Scenarios<-c(Scenarios.F,Scenario.R,Scenario.G)
n.sce<-length(Scenarios)

res.pos <- vector("list", n.sce) 

# Fit several scenarios
for (s in 1:n.sce) {
  res.pos[[s]] <- fit.SPICT.Scenarios(Scenarios = Scenarios[s])
}
save(res.pos,file=file.path(".","2022_Assessment_SPiCT","/SPiCT_Scenarios.RData"))


# SPiCT diagnoses and Outputs ---------------------------------------------
#laod spict fits
#load("./SPiCT.fit.RData")
# Model diagnosis: 

# Convergence
res.pos[[15]][[3]]


#Normality/ Autocorrelation
plotspict.diagnostic(calc.osa.resid(res.pos[[15]][[3]]))

#Retrospective 
fit<-retro(res.pos[[15]][[3]])
plotspict.retro(fit)
plotspict.retro.fixed(fit)
mohns_rho(fit, what = c("FFmsy", "BBmsy"))

#Variance parameters finite
all(is.finite(res.pos[[15]]$sd))

#Realistic production curve
calc.bmsyk(res.pos[[15]][[3]])

# Magnitude difference
calc.om(res.pos[[15]][[3]])

#Sensitivity
set.seed(123)
fit <- check.ini(res.pos[[15]][[3]],ntrials = 30)
fit$check.ini$resmat


#AIC
get.AIC(res.pos[[3]][[3]])

# Results
plotspict.data(res.pos[[15]][[2]])
plot(res.pos[[15]][[3]])
res.pos[[2]][[3]]
plotspict.production(res.pos[[15]][[3]])
plot_production(res.pos[[15]][[3]], plot_it = T)  


# Paul Bouch plots --------------------------------------------------------
colp <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", 
                               "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                               "#4393C3", "#2166AC", "#053061"))) ## intiutively think cold is negative and blue

# Plot production curves from different scenarios one against each other: 
par(mfrow = c(round(n.sce/2), round(n.sce/2)))
for (i in 1:n.sce) {
  plot_production(res.pos[[i]][[3]], plot_it = T)
  
}
par(mfrow = c(1, 1))

corrplot(cov2cor(res.pos[[3]][[3]]$cov.fixed), method = "ellipse", 
         type = "upper", col = colp(200), addCoef.col = "black", diag = FALSE)

res.pos[[1]][[3]]$diag.cov.random

get.cov(res.pos[[1]][[3]], q, n)



# Extracting reference points ---------------------------------------------
# Extract MSY
MSYd<-array(NA,dim=c(n.sce,3))
MSYs<-array(NA,dim=c(n.sce,3))

for (s in 1:n.sce) {
  if (res.pos[[s]][[3]]$opt$convergence == 1) {
    MSYd[s,]<- NA
    MSYs[s,]<- NA
  }
  else {
    MSYd[s,]<- as.numeric(sumspict.drefpoints(res.pos[[s]][[3]])[3,1:3])
    MSYs[s,]<- as.numeric(sumspict.srefpoints(res.pos[[s]][[3]])[3,1:3])
  }
}
MSYd.df <- data.frame(MSYd)
colnames(MSYd.df) <- c("estimate","cilow","ciupp")
MSYd.df$Scenario <- Scenarios
MSYd.df$Ref.Point<-"MSY"
MSYd.df$Model<-"Deterministic"

MSYs.df <- data.frame(MSYs)
colnames(MSYs.df) <- c("estimate","cilow","ciupp")
MSYs.df$Scenario <- Scenarios
MSYs.df$Ref.Point<-"MSY"
MSYs.df$Model<-"Stochastic"

# Bind Stochastic and deterministic MSY
MSY<-rbind(MSYd.df,MSYs.df)


# Extract Bmsy
Bmsyd<-array(NA,dim=c(n.sce,3))
Bmsys<-array(NA,dim=c(n.sce,3))

for (s in 1:n.sce) {
  if (res.pos[[s]][[3]]$opt$convergence == 1) {
    Bmsyd[s,]<- NA
    Bmsys[s,]<- NA
  }
  else{
    Bmsyd[s,]<- as.numeric(sumspict.drefpoints(res.pos[[s]][[3]])[1,1:3])
    Bmsys[s,]<- as.numeric(sumspict.srefpoints(res.pos[[s]][[3]])[1,1:3])
  }
}
Bmsyd.df <- data.frame(Bmsyd)
colnames(Bmsyd.df) <- c("estimate","cilow","ciupp")
Bmsyd.df$Scenario <- Scenarios
Bmsyd.df$Ref.Point<-"Bmsy"
Bmsyd.df$Model<-"Deterministic"

Bmsys.df <- data.frame(Bmsys)
colnames(Bmsys.df) <- c("estimate","cilow","ciupp")
Bmsys.df$Scenario <- Scenarios
Bmsys.df$Ref.Point<-"Bmsy"
Bmsys.df$Model<-"Stochastic"

# Bind Stochastic and deterministic Bmsy
Bmsy<-rbind(Bmsyd.df,Bmsys.df)

# Extract Fmsy
Fmsyd<-array(NA,dim=c(n.sce,3))
Fmsys<-array(NA,dim=c(n.sce,3))

for (s in 1:n.sce) {
  if (res.pos[[s]][[3]]$opt$convergence == 1) {
    Fmsyd[s,]<- NA
    Fmsys[s,]<- NA
  }
  else{
    Fmsyd[s,]<- as.numeric(sumspict.drefpoints(res.pos[[s]][[3]])[2,1:3])
    Fmsys[s,]<- as.numeric(sumspict.srefpoints(res.pos[[s]][[3]])[2,1:3])
  }
}
Fmsyd.df <- data.frame(Fmsyd)
colnames(Fmsyd.df) <- c("estimate","cilow","ciupp")
Fmsyd.df$Scenario <- Scenarios
Fmsyd.df$Ref.Point<-"Fmsy"
Fmsyd.df$Model<-"Deterministic"

Fmsys.df <- data.frame(Fmsys)
colnames(Fmsys.df) <- c("estimate","cilow","ciupp")
Fmsys.df$Scenario <- Scenarios
Fmsys.df$Ref.Point<-"Fmsy"
Fmsys.df$Model<-"Stochastic"

# Bind Stochastic and deterministic Fmsy
Fmsy<-rbind(Fmsyd.df,Fmsys.df)

# Extract B/Bmsy and F/Fmsy
Bbmsy<-array(NA,dim=c(n.sce,3))
Ffmsy<-array(NA,dim=c(n.sce,3))

for (s in 1:n.sce) {
  if (res.pos[[s]][[3]]$opt$convergence == 1) {
    Bbmsy[s,]<- NA
    Ffmsy[s,]<- NA
  }
  else{
    Bbmsy[s,]<- as.numeric(calc.om(res.pos[[s]][[3]])[1,1:3])
    Ffmsy[s,]<- as.numeric(calc.om(res.pos[[s]][[3]])[2,1:3])
  }
}

Ffmsy.df <- data.frame(Ffmsy)
colnames(Ffmsy.df) <- c("cilow","estimate","ciupp")
Ffmsy.df$Scenario <- Scenarios
Ffmsy.df$Ref.Point<-"Ffmsy"
Ffmsy.df$Model<-NA

Bbmsy.df <- data.frame(Bbmsy)
colnames(Bbmsy.df) <- c("cilow","estimate","ciupp")
Bbmsy.df$Scenario <- Scenarios
Bbmsy.df$Ref.Point<-"Bbmsy"
Bbmsy.df$Model<-NA


#Intrinsic growth rate
r<-array(NA,dim=c(n.sce,3))

for (s in 1:n.sce) {
  if (res.pos[[s]][[3]]$opt$convergence == 1) {
    r[s,]<- NA
  } else {
    r[s,]<-get.par("r",res.pos[[s]][[3]],exp=FALSE,CI=.95)[,c("est","ll","ul")]
    #v<-sumspict.parest(res.pos[[s]][[3]])
    #v[rownames(v) == c("r  "),]                   
    #r[s,]<- as.numeric(v[rownames(v) == c("r  "),1:3])
  }
}

r.df <- data.frame(r)
colnames(r.df) <- c("estimate","cilow","ciupp")
r.df$Scenario <- Scenarios
r.df$Ref.Point<-"Intrinsic Growth Rate"
r.df$Model<-NA

# logsdb    logsdf
logsdb<-array(NA,dim=c(n.sce,3))
logsdf<-array(NA,dim=c(n.sce,3))

for (s in 1:n.sce) {
  if (res.pos[[s]][[3]]$opt$convergence == 1) {
    logsdb[s,]<- NA
    logsdf[s,]<- NA
  } else {
    logsdb[s,]<-get.par("logsdb",res.pos[[s]][[3]],exp=FALSE,CI=.95)[,c("est","ll","ul")]
    logsdf[s,]<-get.par("logsdf",res.pos[[s]][[3]],exp=FALSE,CI=.95)[,c("est","ll","ul")]
  }
}
logsdb.df <- data.frame(logsdb)
colnames(logsdb.df) <- c("estimate","cilow","ciupp")
logsdb.df$Scenario <- Scenarios
logsdb.df$Ref.Point<-"logsdb"
logsdb.df$Model<-NA

logsdf.df <- data.frame(logsdf)
colnames(logsdf.df) <- c("estimate","cilow","ciupp")
logsdf.df$Scenario <- Scenarios
logsdf.df$Ref.Point<-"logsdf"
logsdf.df$Model<-NA


# Bind outputs together
out<-rbind(Bmsy,Fmsy,MSY,
           Bbmsy.df, Ffmsy.df,r.df,logsdb.df,logsdf.df)

x11()
ggplot(subset(out,Scenario %in% c("F.4","F.6","R.12","R.13","G.1","G.2") & 
                Ref.Point %in% c("Bbmsy","Ffmsy","Intrinsic Growth Rate",
                                 "logsdb","logsdf")),
       aes(x =  as.factor(Scenario),
           y = estimate,
           group = Model,
           colour = Model,
           fill = Model)) +
  geom_point(position = position_dodge(w = 0.7), size = 2.5,shape=21) +
  geom_errorbar(position = position_dodge(w = 0.7),
                aes(ymin = cilow, ymax = ciupp), size = 0.3,
                width = 0.2) +
  geom_errorbar(position = position_dodge(w = 0.7),
                aes(ymin = cilow, ymax = ciupp), size = .5,
                width = 0.2) +
  labs(x = NULL) + 
  coord_flip()+
  facet_wrap(.~Ref.Point,scales="free_x",ncol=2) +
  scale_color_grey(start=0, end=0.7,guide = guide_legend(reverse=TRUE))+
  scale_fill_grey(start=0, end=0.7,guide = guide_legend(reverse=TRUE))+
  labs(y="")+
  theme_bw()+
  theme(panel.grid.major  = element_blank(),
        axis.title.x = element_text(size=rel(1.4),
                                    margin = margin(t = 25, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=rel(1.4),
                                    margin = margin(t = 0, r = 25, b = 0, l = 0)),
        axis.text.y = element_text(size=rel(1.4)),
        strip.text = element_text(size=14))


# Plot production curves on top of each other: 
for (s in 1:n.sce) {
  if (res.pos[[s]][[3]]$opt$convergence == 1) {
    test<- NA
  }else {
    test<-extractspict.production(res.pos[[s]][[3]],res.pos[[s]][[1]]) #Sourced function 
  }
  if(s==1) {test.f<-test}
  else {test.f<-rbind(test.f,test)}
}

X11()
ggplot(subset(test.f,Scenario %in% c("F.4","F.6","R.12","R.13","G.1","G.2")),
       aes(x=Bplot/Kest,y=Pst/Pstscal,colour=Scenario))+
  geom_line(size=1.2)+
  labs(y=paste0(unique(test.f$ylab)))+
  geom_hline(yintercept = 0,linetype=2)+
  geom_vline(xintercept = .5,linetype=2)+
  ggsci::scale_colour_lancet()


# Extrac Ffmsy and Bbmsy over the years
sce.c<-c(4,6,12,13,14,15)

for (s in 1:length(sce.c)) {
  rep<-res.pos[[sce.c[[s]]]][[3]]
  inp<-res.pos[[sce.c[[s]]]][[3]]$inp
  CI<-.95
  manflag <- any(names(rep) == "man")
  indxmax <- which(inp$time == ifelse(manflag, inp$timerange[2], max(inp$time)))
  BB<-get.par('logBBmsy', rep, exp=TRUE, CI = CI)[1:indxmax,]

  df<-data.frame(BB)
  df$time <- as.numeric(rownames(df))
  df$Scenario<-res.pos[[sce.c[[s]]]][[1]]
  
  if(s==1) {df.f<-df}
  else {df.f<-rbind(df.f,df)}
}

ggplot(subset(df.f,Scenario %in% c("F.4","F.6","R.12","R.13","G.1","G.2")),
              aes(x=time,y=est,colour=Scenario,grooup=Scenario))+
  geom_line(size=1)+
  #geom_ribbon(aes(ymin=ll,ymax=ul,fill=Scenario,colour=Scenario),alpha=.18)+
  scale_x_continuous(breaks = seq(min(df.f$time),max(df.f$time),by=1))+
  labs(y="Bbmsy")+
  ggsci::scale_color_lancet()+
  ggsci::scale_fill_lancet()


for (s in 1:length(sce.c)) {
  rep<-res.pos[[sce.c[[s]]]][[3]]
  inp<-res.pos[[sce.c[[s]]]][[3]]$inp
  CI<-.95
  manflag <- any(names(rep) == "man")
  indxmax <- which(inp$time == ifelse(manflag, inp$timerange[2], max(inp$time)))
  FF<- get.par('logFFmsynotS', rep, exp=TRUE, CI = CI)
  
  df<-data.frame(FF)
  df$time <- inp$time
  df$Scenario<-res.pos[[sce.c[[s]]]][[1]]
  
  if(s==1) {df.f<-df}
  else {df.f<-rbind(df.f,df)}
}

ggplot(df.f,aes(x=time,y=est,colour=Scenario,grooup=Scenario))+
  geom_line(size=1)+
  #geom_ribbon(aes(ymin=ll,ymax=ul,fill=Scenario,colour=Scenario),alpha=.18)+
  scale_x_continuous(breaks = seq(min(df.f$time),max(df.f$time),by=1))+
  labs(y="Ffmsy")+
  ggsci::scale_color_lancet()+
  ggsci::scale_fill_lancet()



# Below from Paul Bauch: 
#for (s in 1:n.sce) {
#  test<-plot_production(res.pos[[s]][[3]],plot_it = F)$pred
#  test$Scenario<-Scenarios[s]
#  
#  if(s==1) {test.f<-test}
#  else {test.f<-rbind(test.f,test)}
#}

#ggplot(test.f,aes(x=B,y=P,colour=Scenario))+
#  geom_line(size=1.2)+
#  scale_y_continuous(limits = c(0,max(test.f$P)))+
#  labs(y="Surplus Production",x="Biomass")+
#  ggsci::scale_colour_lancet()

#TAC
# Hockey-stick MSY rule with the 35th percentile
man.timeline(res.pos[[15]][[3]])

get.TAC(res.pos[[15]][[3]],
        fractiles = list(catch=0.35, bbmsy=0.35, ffmsy=0.35), breakpointB=0.5)
rep<-manage(res.pos[[15]][[3]],
            maninterval = c(2022,2023))
plot2(rep)
plotspict.bbmsy(rep)
sumspict.manage(rep)
plotspict.fb(rep,rel.axes = TRUE)

# What were landings in 2020? 
#IRL: 2001.5607
#SCO: 1564.028
#NI: 703.7378
#Tot:  4269.327

repIntPer <- manage(res.pos[[15]][[3]], 
                     maninterval = c(2021,2022), 
                     #maneval=2021.5,
                     intermediatePeriodCatch = 4270) 

repIntPerC <- add.man.scenario(repIntPer, scenarioTitle = "ices",
                              breakpointB = 0.5,
                              fractiles = list(catch=0.35, bbmsy=0.35, ffmsy=0.35),
                              intermediatePeriodCatch = 4270,
                              maninterval = c(2021,2022))


man.timeline(repIntPerC)
plot2(repIntPerC)
plotspict.fb(repIntPerC,rel.axes = TRUE,man.legend = FALSE)
plotspict.bbmsy(repIntPerC,ylim =c(0,3))
sumspict.manage(repIntPerC)

get.TAC(res.pos[[15]][[3]],breakpointB = 0.5,
        fractiles = list(catch=0.35), #bbmsy=0.35, ffmsy=0.35
        intermediatePeriodCatch = 4270,
        maninterval = c(2021,2022))


# Applying harvest control rules (2/3 rule)

r_NW<-mean(Obs_I_gam1_Above8m$mu.std[Obs_I_gam1_Above8m$Year %in% c(2020,2021)])/mean(Obs_I_gam1_Above8m$mu.std[Obs_I_gam1_Above8m$Year %in% c(2017,2018,2019)])
#0.5668599

#Cy=Cy-1*r
CRE_NW_Landings$Landings_ton_CRE[CRE_NW_Landings$Year %in% 2021]*r_NW
4521*r_NW #4270 are catches in 2020

b<-Obs_I_gam1_Above8m$mu.std[Obs_I_gam1_Above8m$Year %in% 2019]/1.4*(Obs_I_gam1_Above8m$mu.std[Obs_I_gam1_Above8m$Year %in% 2019])
#0.18

#f<- 
#m<-0.95

SVP_nom_index<-read.csv(file = file.path("",
                                         "SVP_nominal_index_AllStocks.csv"))
SVP_nom_index<-subset(SVP_nom_index,CatchType %in% "LPUE")

# For the SW Irish Stock (non-standardized index)
SW_index<-subset(SVP_nom_index,Stock %in% "SW Ireland" &
                   Target %in% "Target")

r_SW<-mean(SW_index$Value[SW_index$Year %in% c(2020,2021)])/mean(SW_index$Value[SW_index$Year %in% c(2017,2018,2019)])
#0.5593662

#For the Celtic Sea Stock (non-standardized index)
CS_index<-subset(SVP_nom_index,Stock %in% "Celtic Sea" &
                   Target %in% "Target")
CS_index$Value[CS_index$Year == 2020]<-mean(CS_index$Value[CS_index$Year == 2019],CS_index$Value[CS_index$Year == 2018])

r_CS<-mean(CS_index$Value[CS_index$Year %in% c(2020,2021)])/mean(CS_index$Value[CS_index$Year %in% c(2017,2018,20119)])
#0.7753303


