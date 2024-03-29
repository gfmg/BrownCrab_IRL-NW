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

# Implement the most complex GAM. Retained Kg as response variable
gam_1<-gam(Retained_Kg ~ 1 + fYear + 
             offset(Log.Effort) + 
             s(SoakDays,bs="cs",k=8) + 
             s(fVesselID,bs="re"),
           data = dat,
           family=tw(),#Gamma(link="log"), Using the Tweedie instead of the Gamma
           method="REML")
par(mfrow=c(2,2))
gam.check(gam_1)
simulationOutput <- simulateResiduals(fittedModel = gam_1, plot = F)
plot(simulationOutput)
plot.gam(gam_1)

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

#Residuals seem decreased in the model with LPUE and Tweedie distribution. 
# In terms of model predictions, we demonstrate later that are almost the same

#AR1 correlations
gamm_1_LPUE<-gamm(LPUE ~ 1 + Year + 
                  s(SoakDays,bs="cs",k=8), #+ s(fVesselID,bs="re"),
                  correlation = corAR1(form = ~ time | fVesselID),
                data = dat,
                family=Gamma(link="log"),
                method="REML")
plot(gamm_1_LPUE)

#glmmTMB test
glmm1<-glmmTMB(LPUE ~ 1 + fYear + 
          #offset(Log.Effort) + 
          #SoakDays + 
          ar1(fYear|fVesselID),
        data = dat,
        family=tweedie(link = "log"))
simulationOutput <- simulateResiduals(fittedModel = glmm1, plot = F)
plot(simulationOutput)

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
predDat1<-data.frame(fYear=unique(dat$fYear),
                    SoakDays=2,
                    Log.Effort=mean(dat$Log.Effort),
                    fVesselID="17")


preds1<-predict.gam(gam_1,newdata = predDat1,
                   se.fit=TRUE,
                   exclude = "s(fVesselID)",
                   type = "link")

predDat1$mu<-preds1$fit
predDat1$se<-preds1$se.fit
predDat1$ci.low<-predDat1$mu - (2 * predDat1$se) #~95% pointwise CI
predDat1$ci.up<-predDat1$mu + (2 * predDat1$se) #~95% pointwise CI

#Response scale: 
predDat1$mu.r<-gam_1$family$linkinv(predDat1$mu)
predDat1$ci.low.r<-gam_1$family$linkinv(predDat1$ci.low)
predDat1$ci.up.r<-gam_1$family$linkinv(predDat1$ci.up)

# And for the model with LPUE rather than Retained_Kg
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
Obs_I1<-predDat1[,c("fYear","mu.r","ci.low.r","ci.up.r")]
Obs_I1$Response<-"Retained"
Obs_I2<-predDat2[,c("fYear","mu.r","ci.low.r","ci.up.r")]
Obs_I2$Response<-"LPUE"

Obs_I<-rbind(Obs_I1,Obs_I2)

#save(Obs_I2,file =file.path("C:/Users/ggonzales/Desktop/gmartin_work_folder/",
#                            "Stock_Assessment/CRE/NorthWest_SPiCT/Code/",
#                            "3_GAM_Standarization/2022_Assessment_Index", 
#                            "CRE_NW_ObsI_gam1_Above8m_21.RData"))


# Dataset for predictions Quarter level -----------------------------------
predDatQ<-expand.grid(fYear=unique(dat_Q$fYear),
                      fQuarter=unique(dat_Q$fQuarter),
                      SoakDays=2,
                      Log.Effort=mean(dat_Q$Log.Effort),
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


#save(Obs_I1_Q,file =file.path(dataDir, 
#                          "CRE_NW_ObsI_gam1_quarter_Above8m.RData"))




# GAM Year smooth ---------------------------------------------------------
#Gam with year as a smoother
dat$ndate<-as.numeric(dat$EventStartDate)
gam_2<-gam(LPUE~ 1 + 
            # offset(Log.Effort) +
             s(ndate,bs="cs") +
             s(SoakDays,bs="cs",k=5) +
             s(fVesselID,bs="re"),
           data = dat,
           family=tw(),#Gamma(link="log"),
           method="REML")
plot.gam(gam_2)
par(mfrow=c(2,2))
gam.check(gam_2)

#Gam_2
predDat2<-data.frame(Year=unique(dat$Year),
                     SoakDays=2,
                     Log.Effort=mean(dat$Log.Effort),
                     fVesselID="17")


preds2<-predict.gam(gam_2,newdata = predDat2,
                    se.fit=TRUE,
                    exclude = "s(fVesselID)",
                    type = "link")
predDat2$mu<-preds2$fit
predDat2$se<-preds2$se.fit
predDat2$ci.low<-predDat2$mu - (2 * predDat2$se)
predDat2$ci.up<-predDat2$mu + (2 * predDat2$se)

#Response scale: 
predDat2$mu.r<-gam_1$family$linkinv(predDat2$mu)
predDat2$ci.low.r<-gam_1$family$linkinv(predDat2$ci.low)
predDat2$ci.up.r<-gam_1$family$linkinv(predDat2$ci.up)


Obs_I2<-predDat2[,c("Year","mu.r","ci.low.r","ci.up.r")]



# Plot standardized indexes ------------------------------------------------

# Scale index of LPUE fixed
Obs_I1$mu.std<-Obs_I1$mu.r/mean(Obs_I1$mu.r)
Obs_I1$ci.low.std<-Obs_I1$ci.low.r/mean(Obs_I1$mu.r)
Obs_I1$ci.up.std<-Obs_I1$ci.up.r/mean(Obs_I1$mu.r)
Obs_I1$model<-"gam_fixed"
names(Obs_I1)[names(Obs_I1) %in% "fYear"]<-"Year"

ggplot(Obs_I1,aes(x=Year,y=mu.std,group=Response,colour=Response))+
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin=ci.low.std,ymax=ci.up.std),alpha=.5) + 
  theme_bw()+
  theme(legend.position = "none")

# Scale index of LPUE smooth
Obs_I2$mu.std<-Obs_I2$mu.r/mean(Obs_I2$mu.r)
Obs_I2$ci.low.std<-Obs_I2$ci.low.r/mean(Obs_I2$mu.r)
Obs_I2$ci.up.std<-Obs_I2$ci.up.r/mean(Obs_I2$mu.r)
Obs_I2$model<-"gam_smooth"

Obs_I<-rbind(Obs_I1,Obs_I2)

# Raw LPUE
dat$LPUE<-dat$Retained_Kg/dat$Effort
raw_index<-aggregate(LPUE~fYear,dat,FUN=mean)
raw_index$mu.std<-raw_index$LPUE/mean(raw_index$LPUE)
raw_index$model<-"nominal"
names(raw_index)[names(raw_index) %in% "fYear"]<-"Year"


index<-ggplot(data = Obs_I,aes(x=Year,y=mu.std,group=model,
                            colour=model,
                            fill=model))+
  geom_ribbon(aes(ymin=ci.low.std,ymax=ci.up.std),alpha=.1)+
  geom_point(aes(shape=model),size=2)+
  geom_line(size=1)+
  scale_fill_lancet()+
  scale_colour_lancet()

index+
  geom_point(data = raw_index,aes(x=Year,y=mu.std),
             size=2,shape=7,colour="black")+
  labs(y="Standarized index",x= "Year", fill="Model", colour="Model",
       shape="Model")+
  theme_bw()+
  scale_y_continuous(limits = c(0,2))+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=11),
        axis.text.x = element_text(angle=45,hjust = 1),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        axis.title.x = element_blank())

aggregate(fVesselID~fYear,dat,FUN=function(x)length(unique(x)))
# Only 3 boats in 2008...
# Below 4 boats between 1996-2004


# Marginal of SoakTime ----------------------------------------------------
x11()
P.smooths<-plot(gam_1,scale=0,seWithMean = TRUE)

Soakupr<-P.smooths[[1]]$fit + (P.smooths[[1]]$se)
Soaklwr<-P.smooths[[1]]$fit - (P.smooths[[1]]$se)


SoakDays.plot<-ggplot()+
  geom_line(aes(x=P.smooths[[1]]$x, y=P.smooths[[1]]$fit))+
  geom_ribbon(aes(ymin=Soaklwr,ymax=Soakupr,x=P.smooths[[1]]$x),fill="grey",alpha=0.4)+
  xlab('Soak Days')+
  ylab('Marginal effect SoakDays')+
  theme_bw()+
  scale_y_continuous(limits = c(-0.35,0.2))+
  theme(axis.title.y=element_text(margin=margin(0,20,0,0)),
        axis.title.x=element_text(margin=margin(20,0,0,0)),
        axis.title=element_text(size=12,face="bold"),
        panel.grid=element_blank())

# Adding rug data
SoakDays.plot + 
  geom_rug(data=dat,aes(x=SoakDays,y=0),sides="b",
           inherit.aes = F,position = "jitter")
