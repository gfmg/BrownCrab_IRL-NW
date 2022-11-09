# File: 3_CRE_NW_modelling.R
# Author: Guillermo Martin 
# Template Created: Tue Apr 13 12:55:04 2021
# ---------------------------------------------------------------------------
# Description:
# Modeling the SVP CRE_NW data for index standarization
# ---------------------------------------------------------------------------
rm(list = ls())

library(ggplot2)
library(mgcv)
library(DHARMa)
library(glmmTMB)
library(ggsci)

# Load data
dataDir<-file.path("")


load(file.path(dataDir,
               "CRE_NW_clean_Above8m_21.Rdata"))

#Turn vessel Registration to ID
CRE_NW_11$fVesselID<-as.factor(as.numeric(CRE_NW_11$Vessel_Registration))
CRE_NW_11$Log.Effort<-log(CRE_NW_11$Effort)

#Removing NA in Retained_Kg
CRE_NW_12<-CRE_NW_11[which(!is.na(CRE_NW_11$Retained_Kg)),]

# Make a dataset wihiouth NA on SoakDays
CRE_NW_13<-CRE_NW_12[which(!is.na(CRE_NW_12$SoakDays)),]

#Remove 2021
# CRE_NW_14<-subset(CRE_NW_13,!Year %in% 2021)

#Selecting only the variable of interest
dat<-CRE_NW_13[, which(names(CRE_NW_13) %in% c("EventStartDate",
                                               "Retained_Kg",
                                               "Discarded_Kg",
                                               "CatchUnits_new",
                                               "Effort",
                                               "Log.Effort",
                                               "SoakDays",
                                               "Year",
                                               "fYear",
                                               "fQuarter",
                                               "fVesselID"
                                               ))]
dat$LPUE<-dat$Retained_Kg/dat$Effort


# GAM fixed year ----------------------------------------------------------

# LPUE as response variable
gam_1_LPUE<-gam(LPUE ~ 1 + fYear + 
                  s(SoakDays,bs="cs",k=8) + 
                  s(fVesselID,bs="re"),
                data = dat,
                family=tw(),#Gamma(link="log"),
                method="REML")
par(mfrow=c(2,2))
gam.check(gam_1_LPUE)
acf(residuals(gam_1_LPUE))
pacf(residuals(gam_1_LPUE))
simulationOutput <- simulateResiduals(fittedModel = gam_1_LPUE, plot = F)
plot(simulationOutput)
plot.gam(gam_1_LPUE)

#Remove vessel ID
gam_1.1<-gam(LPUE ~ 1 + fYear + #Retained_Kg
             #  offset(Log.Effort) + 
               s(SoakDays,bs="cs",k=8),
             data = dat,
             family=tw(),#Gamma(link="log"),
             method="REML")

# Remove soak time smooth
gam_1.2<-gam(LPUE ~ 1 + fYear + #Retained_Kg
            # offset(Log.Effort) + 
             s(fVesselID,bs="re"),
           data = dat,
           family=tw(),#Gamma(link="log"),
           method="REML")

AIC(gam_1_LPUE,gam_1.1,gam_1.2)
# Keep gam_1

# First model checks
plot.gam(gam_1_LPUE)
qq.gam(gam_1_LPUE)
par(mfrow=c(2,2))
gam.check(gam_1_LPUE)

# Residuals vs fitted values
dat$E_1<-resid(gam_1_LPUE)
dat$F_1<-fitted(gam_1_LPUE)

ggplot(dat,aes(x=F_1,y=E_1)) + 
  geom_point()+
  geom_hline(yintercept = 0)

#Residuals vs covariates
ggplot(dat,aes(y=E_1,x=fYear))+
  geom_boxplot() # Little heterogeneity on year

ggplot(dat,aes(y=E_1,x=fQuarter))+
  geom_boxplot() #Some heterogeneity in the quarter. 


# GAM fixed Quarter -------------------------------------------------------
# Try quarter model from 2009 onwards
dat_Q<-subset(dat,Year >= 2009)
gam_quarter<-gam(LPUE ~ 1 + fYear*fQuarter + #Retained_Kg
                   #offset(Log.Effort) + 
                   s(SoakDays,bs="cs",k=8) + 
                   s(fVesselID,bs="re"),
                 data = dat_Q,
                 family=tw(),#Gamma(link="log"),
                 method="REML")
 
# Drop interaction
gam_quarter.1<-gam(LPUE ~ 1 + fYear+fQuarter +
                   #offset(Log.Effort) + 
                   s(SoakDays,bs="cs",k=8) + 
                   s(fVesselID,bs="re"),
                 data = dat_Q,
                 family=tw(),#Gamma(link="log"),
                 method="REML")

AIC(gam_quarter,gam_quarter.1) #Keep the quarter interaction

# First model checks
plot.gam(gam_quarter)
qq.gam(gam_quarter)
par(mfrow=c(2,2))
gam.check(gam_quarter)


# Residuals vs fitted values
dat_Q$E_1<-resid(gam_quarter)
dat_Q$F_1<-fitted(gam_quarter)

ggplot(dat_Q,aes(x=F_1,y=E_1)) + 
  geom_point()+
  geom_hline(yintercept = 0)

#Residuals vs covariates
ggplot(dat_Q,aes(y=E_1,x=fYear))+
  geom_boxplot() # Heterogeneity reduced

ggplot(dat_Q,aes(y=E_1,x=fQuarter))+
  geom_boxplot() #Heterogeneity reduced



# Predictions and standarized index ---------------------------------------


# Dataset for predictions Year level --------------------------------------

predDat2<-data.frame(fYear=unique(dat$fYear),
                     SoakDays=2,
                     fVesselID="17")


preds2<-predict.gam(gam_1_LPUE,newdata = predDat2,
                    se.fit=TRUE,
                    exclude = "s(fVesselID)",
                    type = "link")

predDat2$mu<-preds2$fit
predDat2$se<-preds2$se.fit
predDat2$ci.low<-predDat2$mu - (2 * predDat2$se) #~95% pointwise CI
predDat2$ci.up<-predDat2$mu + (2 * predDat2$se) #~95% pointwise CI

#Response scale: 
predDat2$mu.r<-gam_1_LPUE$family$linkinv(predDat2$mu)
predDat2$ci.low.r<-gam_1_LPUE$family$linkinv(predDat2$ci.low)
predDat2$ci.up.r<-gam_1_LPUE$family$linkinv(predDat2$ci.up)

#Index: 
Obs_I2<-predDat2[,c("fYear","mu.r","ci.low.r","ci.up.r")]
Obs_I2$Response<-"LPUE"

ggplot(Obs_I2,aes(x=fYear,y=mu.r))+
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin=ci.low.r,ymax=ci.up.r),alpha=.5) + 
  theme_bw()+
  facet_wrap(.~Response,scales = "free")


# Dataset for predictions Quarter level -----------------------------------
predDatQ<-expand.grid(fYear=unique(dat_Q$fYear),
                      fQuarter=unique(dat_Q$fQuarter),
                      SoakDays=2,
                      #Log.Effort=mean(dat_Q$Log.Effort),
                      fVesselID="17")


predsQ<-predict.gam(gam_quarter,newdata = predDatQ,
                    se.fit=TRUE,
                    exclude = "s(fVesselID)",
                    type = "link")
                    
predDatQ$mu<-predsQ$fit
predDatQ$se<-predsQ$se.fit
predDatQ$ci.low<-predDatQ$mu - (2 * predDatQ$se) #~95% pointwise CI
predDatQ$ci.up<-predDatQ$mu + (2 * predDatQ$se) #~95% pointwise CI

#Response scale: 
predDatQ$mu.r<-gam_quarter$family$linkinv(predDatQ$mu)
predDatQ$ci.low.r<-gam_quarter$family$linkinv(predDatQ$ci.low)
predDatQ$ci.up.r<-gam_quarter$family$linkinv(predDatQ$ci.up)

predDatQ$YQ<-paste0(predDatQ$fYear,predDatQ$fQuarter)

# Index
Obs_I1_Q<-predDatQ[,c("fYear","fQuarter","YQ","mu.r","ci.low.r","ci.up.r")]
ggplot(Obs_I1_Q,aes(x=YQ,y=mu.r,group=1))+geom_point()+geom_line()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

