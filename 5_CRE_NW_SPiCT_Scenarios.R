# File: CRE_NW_SPiCT_Scenarios.R
# Author: Guillermo Martin 
# Template Created: Thu Apr 29 14:40:14 2021
# ---------------------------------------------------------------------------
# Description:
# Running several different scenario for SPiCT models for the NW Crab data
# ---------------------------------------------------------------------------

rm(list = ls())

library(spict)
library(ggplot2)



dataDir<-"C:/Users/ggonzales/Desktop/gmartin_work_folder/Stock_Assessment/CRE/NorthWest_SPiCT/Data/SPiCT"

# Load data
load(file.path(dataDir,"CRE_NW_ObsC.RData")) #NW CRE landings

#load(file.path(dataDir,"CRE_NW_ObsI_gam1_Allvessels.RData")) # NW CRE index from gam_1
load(file.path(dataDir,"CRE_NW_ObsI_gam1_Above8m.RData")) # NW CRE index from gam_1
load(file.path(dataDir,"CRE_NW_ObsI_gam1_quarter_Above8m.RData")) #NW_CRE_Quarterly based

#Load INLA seasonal index:
load(file.path(dataDir,"CRE_NW_I4_Offshore_ar1_index.RData")) #NW_Crab_index from offshore vessels. Model I3

#Additional plots from Paul Boch to check uncertainty
source(file.path("C:/Users/ggonzales/Desktop/gmartin_work_folder/",
                 "Stock_Assessment/CRE/NorthWest_SPiCT/SPiCT",
                 "cockle spict trial v2.R"))

# Data prep ---------------------------------------------------------------

# Sum all landings by countries
NW_Landings<-subset(NW_Landings,!Year %in% 2020) #Remove 2020 data
CRE_NW_Landings<-aggregate(Landings_ton_CRE~Year,data = NW_Landings,FUN = sum)

# Catch series Year to numeric
CRE_NW_Landings$Year<-as.integer(as.character(CRE_NW_Landings$Year))

# Scale indices of LPUE
Obs_I1$mu.std<-Obs_I1$mu.r/mean(Obs_I1$mu.r)
Obs_I1$ci.low.std<-Obs_I1$ci.low.r/mean(Obs_I1$mu.r)
Obs_I1$ci.up.std<-Obs_I1$ci.up.r/mean(Obs_I1$mu.r)

#Rename index
Obs_I_gam1_Above8m<-Obs_I1

#Year as numeric, reorder data.frame
Obs_I_gam1_Above8m$Year<-as.integer(as.character(Obs_I_gam1_Above8m$fYear))
Obs_I_gam1_Above8m<-Obs_I_gam1_Above8m[order(Obs_I_gam1_Above8m$Year),]

# Quarterly index. 
Obs_I1_Q$YQ.dec<-ifelse(Obs_I1_Q$fQuarter==1,0,
                    ifelse(Obs_I1_Q$fQuarter==2,0.25,
                           ifelse(Obs_I1_Q$fQuarter==3,0.50,
                                  ifelse(Obs_I1_Q$fQuarter==4,0.75,NA))))

Obs_I1_Q$YQ<-as.integer(as.character(Obs_I1_Q$fYear))+Obs_I1_Q$YQ.dec

#Scale quarterly index:
Obs_I1_Q$mu.std<-Obs_I1_Q$mu.r/mean(Obs_I1_Q$mu.r)
Obs_I1_Q$ci.low.std<-Obs_I1_Q$ci.low.r/mean(Obs_I1_Q$mu.r)
Obs_I1_Q$ci.up.std<-Obs_I1_Q$ci.up.r/mean(Obs_I1_Q$mu.r)

#Rename index
Obs_I_gam1_Above8m_Q<-Obs_I1_Q

Obs_I_gam1_Above8m_Q$Year<-as.integer(as.character(Obs_I_gam1_Above8m_Q$fYear))
Obs_I_gam1_Above8m_Q<-Obs_I_gam1_Above8m_Q[order(Obs_I_gam1_Above8m_Q$Year),]

# Ohsshore index year as numeric
I4_index$Year<-as.integer(as.character(I4_index$fYear))

# SPiCT function ----------------------------------------------------------
fit.SPICT.Scenarios<-function(Scenarios) {
  inp <- list()
  
  #Landings 
  inp$obsC <- CRE_NW_Landings$Landings_ton_CRE 
  inp$timeC <- CRE_NW_Landings$Year
  
  if(Scenarios=="F.1"){
    scenario<-"F.1"
    
    # Remove 2008 from time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
  }
  if(Scenarios=="F.2"){
    scenario<-"F.2"
    
    # Only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
  }
  if(Scenarios=="F.3"){
    scenario<-"F.3"
    #Adding offshore index at start of time series

    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    inp$timeI[[2]] <- I4_index$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    inp$obsI[[2]] <- I4_index$mu.std
  }
  if(Scenarios=="F.4"){
    scenario<-"F.4"
    #Using quarterly index
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m_Q$YQ
    inp$timeI[[2]] <- I4_index$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m_Q$mu.std
    inp$obsI[[2]] <- I4_index$mu.std
  }
  if(Scenarios=="F.5"){
    scenario<-"F.5"
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    inp$timeI[[2]] <- I4_index$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    inp$obsI[[2]] <- I4_index$mu.std
    
    ### a prior for the initial biomass depletion - can help if guessable 
    inp$priors$logbkfrac <- c(log(0.8), 0.5, 1) #Low or no exploitation before beginning of available data
    
  }
  if(Scenarios=="F.6"){
    scenario<-"F.6"
    #Adding offshore index at start of time series
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    inp$timeI[[2]] <- I4_index$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    inp$obsI[[2]] <- I4_index$mu.std
    
    ### a prior for the initial biomass depletion - can help if guessable 
    inp$priors$logbkfrac <- c(log(0.8), 0.5, 1) #Low or no exploitation before beginning of available data
    
    #fix the surplus production curve to be symmetrical - optional
    inp$ini$logn <- log(2)
    inp$phases$logn <- -1
  }
  if(Scenarios=="F.7"){
    scenario<-"F.7"
    #Adding offshore index at start of time series
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    inp$timeI[[2]] <- I4_index$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    inp$obsI[[2]] <- I4_index$mu.std
    
    ### a prior for the initial biomass depletion - can help if guessable 
    inp$priors$logbkfrac <- c(log(0.5), 0.5, 1) #Low or no exploitation before beginning of available data
    
    #fix the surplus production curve to be symmetrical - optional
    inp$ini$logn <- log(2)
    inp$phases$logn <- -1
  }
  if(Scenarios=="F.8"){
    scenario<-"F.8"
    #Adding offshore index at start of time series
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    inp$timeI[[2]] <- I4_index$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    inp$obsI[[2]] <- I4_index$mu.std
    
    ### a prior for the initial biomass depletion - can help if guessable 
    inp$priors$logbkfrac <- c(log(0.9), 0.25, 1) #Low or no exploitation before beginning of available data
    
    #fix the surplus production curve to be symmetrical - optional
    inp$ini$logn <- log(2)
    inp$phases$logn <- -1
  }
  if(Scenarios=="F.9"){
    scenario<-"F.9"
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)

    inp$stdevfacC <- rep(1, length(inp$obsC)) 
    inp$stdevfacC[1:15] <- 5 #Extra uncertainty in the first 15 years of landings
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    
    ### a prior for the initial biomass depletion - can help if guessable 
    inp$priors$logbkfrac <- c(log(0.8), 0.5, 1) #Low or no exploitation before beginning of available data
    
    #fix the surplus production curve to be symmetrical - optional
    inp$ini$logn <- log(2)
    inp$phases$logn <- -1
  }
  if(Scenarios=="F.10"){
    scenario<-"F.10"
    
    # only index from 2006 (more reliable)
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year >= 2005)
    
    # Remove 2008 from SVP time series as seems unreliable
    Obs_I_gam1_Above8m<-subset(Obs_I_gam1_Above8m,Year != 2008)
    
    inp$stdevfacC <- rep(1, length(inp$obsC)) 
    inp$stdevfacC[1:15] <- 5 #Extra uncertainty in the first 15 years of landings
    
    inp$timeI <- list()
    inp$timeI[[1]] <- Obs_I_gam1_Above8m$Year+0.5
    inp$timeI[[2]] <- I4_index$Year+0.5
    
    inp$obsI <- list()
    inp$obsI[[1]] <- Obs_I_gam1_Above8m$mu.std
    inp$obsI[[2]] <- I4_index$mu.std
    
    ### a prior for the initial biomass depletion - can help if guessable 
    inp$priors$logbkfrac <- c(log(0.8), 0.5, 1) #Low or no exploitation before beginning of available data
    
    #fix the surplus production curve to be symmetrical - optional
    inp$ini$logn <- log(2)
    inp$phases$logn <- -1
  }
  

  
  # Prior for the intrinsic growth rate?
  
  
  res <- fit.spict(inp)
  
  return(list(scenario, inp, res))
}


Scenarios<-c(paste0("F.",seq(1,10,by=1)))#paste0("F.",seq(1,9,by=1))
n.sce<-length(Scenarios)

res.pos <- vector("list", n.sce) 

# Fit several scenarios
for (s in 1:n.sce) {
  res.pos[[s]] <- fit.SPICT.Scenarios(Scenarios = Scenarios[s])
}



# SPiCT diagnoses and Outputs ---------------------------------------------
# Model diagnosis: 

# Convergence
res.pos[[10]][[3]]


#Normality/ Autocorrelation
plotspict.diagnostic(calc.osa.resid(res.pos[[10]][[3]]))

#Retrospective 
fit<-retro(res.pos[[10]][[3]])
plotspict.retro(fit)
mohns_rho(fit, what = c("FFmsy", "BBmsy"))

#Variance parameters finite
all(is.finite(res.pos[[10]]$sd))

#Realistic production curve
calc.bmsyk(res.pos[[9]][[3]])

# Magnitude difference
calc.om(res.pos[[10]][[3]])

#Sensitivity
set.seed(123)
fit <- check.ini(res.pos[[10]][[3]],ntrials = 30)
fit$check.ini$resmat


# Results
plotspict.data(res.pos[[8]][[2]])
plot(res.pos[[10]][[3]])
res.pos[[9]][[3]]
plot_production(res.pos[[1]][[3]], plot_it = T)  



corrplot(cov2cor(res.pos[[5]][[3]]$cov.fixed), 
         method = "ellipse", type = "upper", col = colp(200), 
         addCoef.col = "black", diag = FALSE)



# Extracting reference points ---------------------------------------------

# Extract MSY
MSYd<-array(NA,dim=c(n.sce,3))
MSYs<-array(NA,dim=c(n.sce,3))

for (s in 1:n.sce) {
  MSYd[s,]<- as.numeric(sumspict.drefpoints(res.pos[[s]][[3]])[3,1:3])
  MSYs[s,]<- as.numeric(sumspict.srefpoints(res.pos[[s]][[3]])[3,1:3])
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


# Extract Bbmsy
Bmsyd<-array(NA,dim=c(n.sce,3))
Bmsys<-array(NA,dim=c(n.sce,3))

for (s in 1:n.sce) {
  Bmsyd[s,]<- as.numeric(sumspict.drefpoints(res.pos[[s]][[3]])[1,1:3])
  Bmsys[s,]<- as.numeric(sumspict.srefpoints(res.pos[[s]][[3]])[1,1:3])
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
  Fmsyd[s,]<- as.numeric(sumspict.drefpoints(res.pos[[s]][[3]])[2,1:3])
  Fmsys[s,]<- as.numeric(sumspict.srefpoints(res.pos[[s]][[3]])[2,1:3])
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
  Bbmsy[s,]<- as.numeric(calc.om(res.pos[[s]][[3]])[1,1:3])
  Ffmsy[s,]<- as.numeric(calc.om(res.pos[[s]][[3]])[2,1:3])
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

# Bind reference points together
out<-rbind(Bmsy,Fmsy,MSY,
           Bbmsy.df, Ffmsy.df)
out$Scenario<-factor(out$Scenario,levels = paste0("F.",seq(1,10,by=1)))

ggplot(subset(out,Scenario %in% c("F.1","F.3","F.4","F.5","F.6","F.7","F.8",
                                  "F.9","F.10")),
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
        axis.text.y = element_text(size=rel(1.4)))


