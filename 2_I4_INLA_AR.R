# File: 2_I4_INLA_AR.R
# Author: Guillermo Martin 
# Template Created: Tue Jun 08 14:10:34 2021
# ---------------------------------------------------------------------------
# Description:
# Standardizing the time series (1990-2005) of vivier ofshore Brown Crab vessels
# LPUE's
#
# Using an Spatio-Temporal model within INLA. 
# 
# Descriptions of the standardization method based on the work by: 
# Zhou et al., 2019. https://doi.org/10.1093/icesjms/fsz034
# 
# Some of the code for the implementation of the Autoregressive INLA model
# from: 
# Zuur et al.,2017 
# https://www.highstat.com/index.php/beginner-s-guide-to-regression-models-with-spatial-and-temporal-correlation

# ---------------------------------------------------------------------------

rm(list=ls())

# Setting Directory and Loading data --------------------------------------

setwd(file.path("C:/Users/ggonzales/Desktop/gmartin_work_folder/Inshore_Data_compilations/Brown Crab/Modelling"))


scriptDir<-file.path("C:/Users/ggonzales/Desktop/gmartin_work_folder/Inshore_Data_compilations/Brown Crab/Modelling")

OutDir<-file.path("C:/Users/ggonzales/Desktop/gmartin_work_folder/Inshore_Data_compilations/Brown Crab/Outputs")

DataDir<-file.path("C:/Users/ggonzales/Desktop/gmartin_work_folder/Inshore_Data_compilations/Brown Crab/Data")


## Libraries
library(INLA)
library(sf)
library(colorRamps) 

## Data
load(file.path(paste(DataDir,"/","CRE_offshore_clean_9006_withMask.RData",sep="")))


#Log of effort
CRE_data14$Log.E<-log(CRE_data14$Gear_Units)

#Standarized continuous covariates
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

CRE_data14$fYear<-factor(CRE_data14$fYear)
CRE_data14$fVESSElID<-as.factor(CRE_data14$fVESSElID)

# Spatiotemporal model ----------------------------------------------------

# First we need to define the mesh for the spatial-temporl correlation
Loc<-cbind(CRE_data14$XKm,CRE_data14$YKm)
D<-dist(Loc)

par(mfrow = c(1, 1))
hist(D,freq=T,breaks = 50,main="",
     xlab="Distance between Hauls (Km)",
     ylab = "Frequency")
axis(side=1, at=seq(0,400, 10))

plot(x = sort(D), 
     y = (1:length(D))/length(D), 
     type = "l",
     xlab = "Distance between sites (km)",
     ylab = "Cumulative proportion")
text(10, 1, "B", cex = 1.5)

#Second figure indicates that 50% of the distances between sites are less than
# ~60Km apart. 

# Since this is our initial analysis, we will set the range to 80Km. 
#However, it seems is approx. 25 Km
RangeGuess<-80
MaxEdge<-RangeGuess/5
ConvHull<-inla.nonconvex.hull(Loc,convex = -0.05)
mesh1<-inla.mesh.2d(boundary = ConvHull,
                    max.edge = c(1,5)* MaxEdge,
                    cutoff = MaxEdge/5,
                    crs = CRS("+proj=utm +zone=29N +units=km"))

plot(mesh1,asp=1)
points(Loc,col=2,pch=16,cex=1)

# A mesh with a finer scale? 
RangeGuess2<-35
MaxEdge2<-RangeGuess2/5
ConvHull<-inla.nonconvex.hull(Loc,convex = -0.05)
mesh2<-inla.mesh.2d(boundary = ConvHull,
                    max.edge = c(1,5)* MaxEdge2,
                    cutoff = MaxEdge2/5)

plot(mesh2,asp=1)
points(Loc,col=2,pch=16,cex=1)

#Number of vertices  
inla.mesh.components(mesh1)
inla.mesh.components(mesh2)

# We need to group data by years

#Repl<-as.numeric(CRE_data14$fYear)
Group3<-as.integer(as.character(CRE_data14$fYear))-1990
NGroups<-length(unique(Group3))
#We have 14 years

#Defining weighting factors of the mesh
A1<-inla.spde.make.A(mesh1,
                     group = Group3,
                     loc = Loc)
dim(A1) # 17946 sampling locations and 602 vertices multiplied by 14 years!

#Defining the SPDE, default and with pcpriors. We said the range was approximately
#at 25Km
spdeAr<-inla.spde2.matern(mesh1,alpha = 2)

spde.pc.1<-inla.spde2.pcmatern(mesh1,
                               prior.range = c(20,0.05),
                               prior.sigma = c(3,0.05)
                               )



# Finer mesh
A2<-inla.spde.make.A(mesh2,
                     group = Group3,
                     loc = Loc)
dim(A2) # 19251 sampling locations and 882 vertices multiplied by 14 years!

spdeAr2<-inla.spde2.matern(mesh2,alpha = 2)

#Define the spatial random field
w1.index<-inla.spde.make.index(name="w",
                               n.spde=mesh1$n,
                               n.group = NGroups)

w2.index<-inla.spde.make.index(name="w",
                               n.spde=mesh2$n,
                               n.group = NGroups)


# Define rw2 smoother for Soak Time ---------------------------------------
## Hyperparameter for the rw2
U<-0.5
hyper.rw<-list(theta=list(prior="pc.prec",
                          param=c(U,0.05)))
# (0.5,0.05 for the smoother looks ok)



# Defining priors for hyperparameters -------------------------------------

# Penalised complexity (PC) prior for Vessel ID
M1 <- lm(Target_Units_Landed ~ 1, data = CRE_data14)

# sigma  of the gaussian distribution is estimated as 
summary(M1)$sigma
priorpc <- list(prec = list(prior = "pc.prec", param = c(700, 0.05)))


# Define the stack of the data --------------------------------------------
# Next, we define the Stack of the data. It links the spatial field and covariates to the
#observed data. We will still use the smoother for Soak time

#Categorical covarites need to be defined as dummy variables with zeros and ones
XfYear<- model.matrix(~fYear,data=CRE_data14)

N<-nrow(CRE_data14)
Covariates<-data.frame(Intercept=rep(1,N),
                       #fYear=CRE_data14$fYear,
                       fYear1992=XfYear[,"fYear1992"],
                       fYear1993=XfYear[,"fYear1993"],
                       fYear1994=XfYear[,"fYear1994"],
                       fYear1995=XfYear[,"fYear1995"],
                       fYear1996=XfYear[,"fYear1996"],
                       fYear1997=XfYear[,"fYear1997"],
                       fYear1998=XfYear[,"fYear1998"],
                       fYear1999=XfYear[,"fYear1999"],
                       fYear2000=XfYear[,"fYear2000"],
                       fYear2001=XfYear[,"fYear2001"],
                       fYear2002=XfYear[,"fYear2002"],
                       fYear2003=XfYear[,"fYear2003"],
                       fYear2004=XfYear[,"fYear2004"],
                       fYear2005=XfYear[,"fYear2005"],
                       fYear2006=XfYear[,"fYear2006"],
                       #fQuarter=CRE_data14$fQuarter, We will keep it simple Withouth
                       #Quarter
                       Soak_Time=CRE_data14$Soak_Time,
                       Gear_Saturation5K3D.std=CRE_data14$Gear_Saturation5K3D.std,
                       fVesselID=CRE_data14$fVESSElID,
                       Log.E=CRE_data14$Log.E)


#D=Combine the response variable, projection matrix, covariates 
#and spatial random field
Stk.est<- inla.stack(tag="Fit",
                     data=list(y=CRE_data14$Target_Units_Landed),
                     A=list(1,A1),
                     effects=list(Covariates=Covariates,
                                  w= w1.index))


# Stack for prediction in 1*1 degree --------------------------------------
load(file.path("C:/Users/ggonzales/Desktop/gmartin_work_folder/",
               "Stock_Assessment/CRE/NorthWest_SPiCT/GitHub_repo/BrownCrab_IRL-NW",
               "Prediction_grid.RData"))


Group4<-as.integer(as.character(predDat2$fYear))-1990
NGroups<-length(unique(Group4))

A.pred2 <- inla.spde.make.A(mesh=mesh1, 
                            loc = cbind(predDat2$XKm, predDat2$YKm), 
                            group=Group4, 
                            n.group=16) # 

# Covariates.pred
XfYear.Pred<- model.matrix(~fYear,data=predDat2)
N.pred<-nrow(predDat2)
Covariates.pred2<-data.frame(Intercept=rep(1,N.pred),
                             #fYear=predDat2$fYear,
                             fYear1992=XfYear.Pred[,"fYear1992"],
                             fYear1993=XfYear.Pred[,"fYear1993"],
                             fYear1994=XfYear.Pred[,"fYear1994"],
                             fYear1995=XfYear.Pred[,"fYear1995"],
                             fYear1996=XfYear.Pred[,"fYear1996"],
                             fYear1997=XfYear.Pred[,"fYear1997"],
                             fYear1998=XfYear.Pred[,"fYear1998"],
                             fYear1999=XfYear.Pred[,"fYear1999"],
                             fYear2000=XfYear.Pred[,"fYear2000"],
                             fYear2001=XfYear.Pred[,"fYear2001"],
                             fYear2002=XfYear.Pred[,"fYear2002"],
                             fYear2003=XfYear.Pred[,"fYear2003"],
                             fYear2004=XfYear.Pred[,"fYear2004"],
                             fYear2005=XfYear.Pred[,"fYear2005"],
                             fYear2006=XfYear.Pred[,"fYear2006"],
                             Soak_Time=rep(2,N.pred),
                             Gear_Saturation5K3D.std=rep(0,N.pred),
                             fVesselID=NA,
                             Log.E=rep(mean(CRE_data14$Log.E),N.pred)
)


Stk.pred2<- inla.stack(tag="Pred",
                       data=list(y=NA),
                       A=list(1,A.pred2),
                       effects=list(Covariates=Covariates.pred2,
                                    w= w1.index))



# Running the model -------------------------------------------------------

# Combine prediction and estimates stack
stck.all=inla.stack(Stk.est,Stk.pred2)

#Defininf prior for the overdispersion parameter
#load(file.path(paste(scriptDir,"/","Theta_spatiotemporal_mesh1.Rdata",sep="")))
#Hyper.NB <- list(size = list(initial = theta.pm.I2.nb.2.SPATIO, fixed = TRUE))

#Define formula
f_INLA_AR<- y ~ -1 + Intercept + fYear1992 + 
  fYear1993 + fYear1994 + fYear1995 + fYear1996 + fYear1997 + fYear1998 +
  fYear1999 + fYear2000 + fYear2001 + fYear2002 + fYear2003 +
  fYear2004 + fYear2005 + fYear2006 +
  Gear_Saturation5K3D.std +
  f(Soak_Time,model="rw2",hyper=hyper.rw) +  
  offset(Log.E) + 
  f(fVesselID,model = "iid",hyper = priorpc)+
  f(w, model = spde.pc.1, group= w.group, control.group=list(model="ar1"))

## Executing the model
I4<-inla(f_INLA_AR,
             family = "gamma",
             data = inla.stack.data(stck.all),
             control.compute = list(dic = TRUE, waic = TRUE , cpo=TRUE),
             control.predictor = list(A=inla.stack.A(stck.all),link=1),
             #control.family = list(hyper = Hyper.NB),
             #control.inla = list(strategy = "adaptive"),
             #num.threads = 2,
             quantiles = c(0.025, 0.975),
             verbose=T)



# Extract predicted values and fitted:
idx.fit<- inla.stack.index(stck.all, 'Fit')$data
idx.pred<- inla.stack.index(stck.all, 'Pred')$data # Indices of the predicted values
# in the stack


# Ploting the spatial random field ----------------------------------------
proj <- inla.mesh.projector(mesh1, 
                            xlim = range(mesh1$loc[,1]), 
                            ylim = range(mesh1$loc[,2]), 
                            dims = c(300, 300))


w     <- I4$summary.random$w$mean
Years <- c(rep(1991:2006,by=1))

df.f <- data.frame(x=numeric(),
                 y=numeric(), 
                 Year=factor(), 
                 mean_s=numeric()) 

Loc.Y.f <- data.frame(X1=numeric(),
                   X2=numeric(), 
                   fVessel=factor(), 
                   Year=factor()) 

for (i in 1:length(Years)){
  w.pm <- w[w1.index$w.group == i]
  
  Loc.Y<-data.frame(Loc[CRE_data14$fYear==Years[i],])
  Loc.Y$fVessel<-CRE_data14$fVESSElID[CRE_data14$fYear==Years[i]]
  Loc.Y$Year<-Years[i]
  
  field.proj <- inla.mesh.project(proj, w.pm)
  
  df<-expand.grid(x=proj$x,y=proj$y)
  df$Year<-Years[i]
  df$mean_s <- as.vector(field.proj)

    if(Years==1){
      df.f=df 
      Loc.Y.f<-Loc.Y
    }
  if(Years>1){
    df.f=rbind(df.f,df )
    Loc.Y.f =rbind(Loc.Y.f,Loc.Y)
  }
  
}
  

coastline_UTM<-spTransform(coastline, CRS("+proj=utm +zone=29N +units=km"))
UK_UTM<-spTransform(UK, crs(coastline_UTM))
Counties_UTM<-spTransform(Counties, crs(coastline_UTM))


p<-ggplot() +
  geom_raster(df.f, mapping= aes(x = x, y = y, fill = mean_s))+
  geom_polygon(data = coastline_UTM,aes(x=long,y=lat,group=group),
               fill='antiquewhite', colour="black",size=1)+
  geom_polygon(data = UK_UTM,aes(x=long,y=lat,group=group),
               fill='antiquewhite', colour="black",size=1)+
  geom_polygon(data = Counties_UTM,aes(x=long,y=lat,group=group),
               fill='antiquewhite', colour="black",size=1)+
  scale_fill_gradientn(colours = matlab.like(500),na.value = 'transparent',
                       limits=c(min(w),max(w)))+
  coord_cartesian(xlim = range(mesh1$loc[,1]),
                  ylim = range(mesh1$loc[,2]))+
  labs(x ="Easting (Km)", y = "Northing (Km)", fill="Mean\nSpatial\nrandom field")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "aliceblue"), axis.line = element_line(colour = "black"),
        plot.title = element_text(color="black", size=rel(1.5), face="bold",hjust = 0.5),
        axis.title.x = element_text(size=rel(1.4), face="bold",
                                    margin = margin(t = 25, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=rel(1.4), face="bold",
                                    margin = margin(t = 0, r = 25, b = 0, l = 0)), 
        axis.text = element_text(color="black", size=rel(1.3)),
        axis.text.y  = element_text(hjust=1),
        legend.title=element_text(size=rel(1.4),face="bold"),
        strip.text = element_text(size=rel(1.4)))+
  ggtitle("{frame_time}") +
  transition_time(Year) +
  ease_aes("linear") +
  enter_fade() +
  exit_fade()

animate(p, fps = 5, width = 750, height = 450)
anim_save("SpatialRandomfield4.gif")


GMRF_mean<-ggplot() +
  geom_raster(df.f, mapping= aes(x = x, y = y, fill = mean_s))+
  geom_polygon(data = coastline_UTM,aes(x=long,y=lat,group=group),
               fill='antiquewhite', colour="black",size=1)+
  geom_polygon(data = UK_UTM,aes(x=long,y=lat,group=group),
               fill='antiquewhite', colour="black",size=1)+
  geom_polygon(data = Counties_UTM,aes(x=long,y=lat,group=group),
               fill='antiquewhite', colour="black",size=1)+
  scale_fill_gradientn(colours = matlab.like(500),na.value = 'transparent',
                       limits=c(min(w),max(w)))+
  coord_cartesian(xlim = range(mesh1$loc[,1]),
                  ylim = range(mesh1$loc[,2]))+
  facet_wrap(.~Year)+
  labs(x ="Easting (Km)", y = "Northing (Km)", fill="Mean\nGMRF")+
  theme(panel.spacing.x = unit(.5, "lines"),
        panel.spacing.y = unit(.5, "lines"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = gray(.5),linetype = "dashed", size = 0.5),
        panel.margin = unit(0, "lines"),
        panel.background = element_rect(fill = "aliceblue"),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x = element_text(size=rel(1.4), face="bold",
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=rel(1.4), face="bold",
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)))

GMRF_mean

# Soak_Time smoother plot -------------------------------------------------
# Visualizing the smoother
##Plotting rw smooth
Smoother<-I4$summary.random$Soak_Time
Smoother<-Smoother %>%
  rename(SeLo="0.025quant",
         SeUp="0.975quant")

ggplot(data=Smoother,aes(x=ID,y=mean))+
  geom_line(size=.8)+
  geom_ribbon(data=Smoother,aes(x=ID,ymax=SeLo,ymin=SeUp),alpha=0.6)+
  geom_rug(data=CRE_data14,aes(x=Soak_Time,y=0),sides="b",
           inherit.aes = F,position = "jitter")+
  geom_hline(yintercept = 0,linetype="dashed")+
  scale_y_continuous(limits = c(-0.3,0.1))+
  scale_x_continuous(breaks = seq(1,11,by=1))+
  xlab("Days")+ylab("Smoother")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.title.x = element_text(size=rel(1.4),
                                    margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=rel(1.4),
                                    margin = margin(t = 5, r = 0, b = 0, l = 0)),                  
        axis.line = element_line(colour = "black"))
 


# Value of the autoregressive components ----------------------------------
#Value of phi
#Make sure the 4 is phi!!!!
phi.mean<-I4$summary.hyperpar[6,"mean"]
phi.marginal<-data.frame(I4$marginals.hyperpar[[6]])
phi.selo<-I4$summary.hyperpar[6,"0.025quant"]
phi.seup<-I4$summary.hyperpar[6,"0.975quant"]


ggplot(data = phi.marginal,
       aes(x=x,y=y))+
  geom_line(size=1.2,colour="skyblue")+
  geom_vline(xintercept = phi.mean,colour="green",size=1) +
  geom_vline(xintercept = phi.selo,colour="darkgreen",size=1,linetype="dashed") +
  geom_vline(xintercept = phi.seup,colour="darkgreen",size=1,linetype="dashed") +
  labs(x =expression(phi), y = "Density")+
  expand_limits(x = 0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(color="black", size=rel(1.5), face="bold",hjust = 0.5),
        axis.title.x = element_text(size=rel(1.4), face="bold",
                                    margin = margin(t = 25, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=rel(1.4), face="bold",
                                    margin = margin(t = 0, r = 25, b = 0, l = 0)), 
        axis.text = element_text(color="black", size=rel(1.3)),
        axis.text.y  = element_text(hjust=1),
        legend.title=element_text(size=rel(1.4),face="bold"),
        strip.text = element_text(size=rel(1.4)))



# Value of the spatial correlation ----------------------------------------
# Hyperparameters spatial random field
summary(I4)

MySpatialParams <- function(Model, ThisSpde) {
  SpFi <- inla.spde2.result(inla = Model, name = "w", spde = ThisSpde, do.transfer = TRUE) 
  Kappa  <- inla.emarginal(function(x) x, SpFi$marginals.kappa[[1]] )
  sigmau <- inla.emarginal(function(x) sqrt(x),SpFi$marginals.variance.nominal[[1]] )
  Range <- inla.emarginal(function(x) x, SpFi$marginals.range.nominal[[1]] )
  Out <- c(Kappa, sigmau, Range)
  names(Out) <- c("Kappa", "Sigma_u", "Range")
  Out
}

MySpatialParams(I4,spde.pc.1)

#IMPOSED correlation
Kappa<- 0.1057216
sigmau <- 0.2465466
r<-26.9570763   #Distance at whic dependecy < than 0.1


D<-as.matrix(dist(mesh1$loc[,1:2]))
d.vec<-seq(0,max(D),length=100)
Cor.M<-(Kappa*d.vec)*besselK(Kappa*d.vec,1)
Cor.M[1]<-1


range.correlation<-data.frame(d.vec=d.vec,
                              Cor.M=Cor.M)

ggplot()+
  geom_line(data=range.correlation,aes(x=d.vec,y=Cor.M),size=1)+
  geom_hline(yintercept = 0.1,colour="red",linetype = "dashed",size=.8)+
  labs(x ="Distance (Km)", y = "Correlation")+
  scale_x_continuous(limits = c(0,100))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(color="black", size=rel(1.5), face="bold",hjust = 0.5),
        axis.title.x = element_text(size=rel(1.4), face="bold",
                                    margin = margin(t = 25, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=rel(1.4), face="bold",
                                    margin = margin(t = 0, r = 25, b = 0, l = 0)), 
        axis.text = element_text(color="black", size=rel(1.3)),
        axis.text.y  = element_text(hjust=1),
        legend.title=element_text(size=rel(1.4),face="bold"),
        strip.text = element_text(size=rel(1.4)))


# Model Validation --------------------------------------------------------
N     <- nrow(CRE_data14)
mu.G.2   <- I4$summary.fitted.values[1:N,"mean"] 

#Parameter phi from the Gamma variance
phi.gamma    <- I4$summary.hyperpar[1, "mean"]
VarY1 <- mu.G.2^2 / phi.gamma
E.G.2    <- (CRE_data14$Target_Units_Landed - mu.G.2) / sqrt(VarY1)

CRE_data14$E2.2<-E.G.2
CRE_data14$mu.G.2<-mu.G.2

## Plot residuals vs fitted values
par(mfrow = c(1, 1))
plot(x = mu.G.2, 
     y = E.G.2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2, col = 1)


##
ggplot(data=CRE_data14,aes(x=mu.G.2,y=Target_Units_Landed))+
  geom_point()+
  geom_smooth(colour="red")+
  facet_wrap(.~fYear)+
  xlab("Fitted values")+ylab("Observed values")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.margin = unit(0, "lines"),
        axis.title.x = element_text(size=rel(1.4), face="bold",
                                    margin = margin(t = 25, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=rel(1.4), face="bold",
                                    margin = margin(t = 0, r = 25, b = 0, l = 0)))



## Plot PIT
hist(I4$cpo$pit,breaks=10)

##Plot posterio predictive p-values
predicted.p.value<-c()
n<-length(CRE_data14[,1])
for (i in 1:n) {
  predicted.p.value[i]<-inla.pmarginal(q=CRE_data14$Target_Units_Landed[i],
                                       marginal=I4$marginals.fitted.values[[i]])
}

hist(predicted.p.value,main="",xlab = "Posterior predictive p-value")


## Plot residuals vs each covariate
MyVar <- c("Soak_Time","Gear_Saturation5K3D.std","Log.E")
MyMultipanel.ggp2(Z = CRE_data14, 
                  varx = MyVar, 
                  vary = "E2.2", 
                  ylab = "Pearson residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)

Myxyplot(CRE_data14, MyVar, "E2.2", MyYlab = "Pearson residuals")

# Residuals vs categorical covariates
# Pearson residuals against categorical covariates
par(mfrow = c(1,2))
boxplot(E2.2 ~ fYear, data = CRE_data14) #Ok
boxplot(E2.2 ~ fVESSElID, data = CRE_data14) 
boxplot(E2.2 ~ fQuarter, data = CRE_data14) 


# Spatial dependency? Anysotropy? -----------------------------------------
years=unique(levels(factor(CRE_data14$fYear)))
Vario.all<-data.frame()

for (i in 1:length(years)) {
  ind_y=which(CRE_data14$fYear==years[i])
  CRE_data14_Year1=CRE_data14[ind_y,]
  
  mydata<-data.frame(CRE_data14_Year1$E2.2,
                     CRE_data14_Year1$XKm,
                     CRE_data14_Year1$YKm)
  coordinates(mydata)<-c("CRE_data14_Year1.XKm","CRE_data14_Year1.YKm")
  Vario<-variogram(object = CRE_data14_Year1$E2.2 ~ 1,
                   data = mydata,
                   cressie=TRUE,
                   cutoff=150, # Defining a a distance limit of 150Km for semivariance calculation
                   width=2)
  
  Vario$fYear<-years[i]
  
  if(i==1){
    Vario.all=Vario 
  }
  if(i>1){
    Vario.all=rbind(Vario.all,Vario) 
  }
}

p<-ggplot(data=Vario.all,aes(x=dist,y=gamma))
p<-p+geom_point()
p<-p+geom_smooth(method = "gam",formula = y ~ s(x,bs="cs"),colour="black")
#p<-p+ylim(0,1)
p<-p+theme(text=element_text(size=15))
p<-p+xlab("Distance (Km)")
p<-p+ylab("Sample variogram")
p<-p+facet_wrap(.~fYear)
p<-p+scale_y_continuous(limits = c(0,1.5))
p



# Constructing CPUE index -------------------------------------------------

I4_pred2<-I4$summary.fitted.values[idx.pred,]

I4_predI<-data.frame(fYear=predDat2$fYear,
                   mu=I4_pred2[,"mean"],
                   selow=I4_pred2[,"0.025quant"],
                   seup=I4_pred2[,"0.975quant"])
#Year index:
I4_index<-aggregate(I4_predI[c("mu", "selow","seup")],
                    by=list(I4_predI$fYear),FUN=mean)
colnames(I4_index)[1]<-"fYear"

# Scale indices of LPUE
I4_index$mu.std<-I4_index$mu/mean(I4_index$mu)
I4_index$selow.std<-I4_index$selow/mean(I4_index$mu)
I4_index$seup.std<-I4_index$seup/mean(I4_index$mu)

I4_index<-I4_index[,c("fYear","mu.std", "selow.std", "seup.std")]


ggplot(data=I4_index,aes(x=fYear,y=mu.std,group=1))+
  geom_point()+
  geom_line()+
  geom_ribbon(data=I4_index,aes(x=fYear,ymin=selow.std,ymax=seup.std),alpha=.3)

I4_index<-subset(I4_index,select =  -c(Est.Y,selow.Y,seup.Y,mu.Year.mean))


## Plotting the predictions
predictions<-cbind(predDat2,pred2)
predictions$mean_New<-predictions$mean

ggplot(predictions,aes(x=XKm, y= YKm,colour=mean_New))+
  geom_polygon(data = coastline_UTM,aes(x=long,y=lat,group=group),
               fill='grey', colour="black",size=1)+
  geom_point(shape=15,size=1)+
  scale_colour_gradientn(colours = matlab.like(500),na.value = 'transparent')+
  scale_y_continuous(limits = c(range(predictions$YKm)[1],range(predictions$YKm)[2]))+
  facet_wrap(.~fYear)+
  labs(x ="Easting (Km)", y = "Northing (Km)", fill="Predicted\nCPUE",colour="Predicted\nCPUE")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(color="black", size=rel(1.5), face="bold",hjust = 0.5),
        axis.title.x = element_text(size=rel(1.4), face="bold",
                                    margin = margin(t = 25, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=rel(1.4), face="bold",
                                    margin = margin(t = 0, r = 25, b = 0, l = 0)), 
        axis.text = element_text(color="black", size=rel(1.3)),
        axis.text.y  = element_text(hjust=1),
        legend.title=element_text(size=rel(1.4),face="bold"),
        strip.text = element_text(size=rel(1.4)))



# Plotting the grid of predictions with the mesh
coastline.sf<-st_as_sf(coastline)
UK.sf<-st_as_sf(UK)
Counties.sf<-st_as_sf(Counties)
cord.UTM.sf<-st_as_sf(cord.UTM)
points.Sf <- st_as_sf(CRE_data14, coords = c("Longitude_Start", "Latitude_Start"), 
                      crs = st_crs(coastline.sf))

coastline.sf<-st_transform(coastline.sf,
                           "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
cord.UTM.sf<-st_transform(cord.UTM.sf,crs = st_crs(coastline.sf))
Counties.sf<-st_transform(Counties.sf,crs = st_crs(coastline.sf))
points.Sf<-st_transform(points.Sf,crs = st_crs(coastline.sf))
UK.sf<-st_transform(UK.sf,crs = st_crs(coastline.sf))

ggplot()+
  gg(mesh1,size = 40,edge.color = "azure4",int.color = "blue",ext.color = "black",
     crs = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")) +
  geom_sf(data = coastline.sf,colour="black",size=1)+
  geom_sf(data = Counties.sf,colour="black",fill= "antiquewhite")+
  geom_sf(data = UK.sf,colour="black",fill= "antiquewhite")+
  geom_sf(data=cord.UTM.sf,fill="Red",colour="black",alpha=0.9,pch=22,size=2.5)+
  xlab("Longitude (°)")+ylab("Latitutde (°)")+
  annotation_scale(location = "tl", width_hint = 0.5) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(2950000 ,3400000),
           ylim = c(3650000,4000000))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = gray(.5),linetype = "dashed", size = 0.5),
        panel.margin = unit(0, "lines"),
        panel.background = element_rect(fill = "aliceblue"),
        axis.title.x = element_text(size=rel(1.4), face="bold",
                                    margin = margin(t = 25, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=rel(1.4), face="bold",
                                    margin = margin(t = 0, r = 25, b = 0, l = 0)))


#Plotting the INLA mesh with the points
ggplot()+
  gg(mesh1,size = 40,edge.color = "azure4",int.color = "red",ext.color = "black",
     crs = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")) +
  geom_sf(data = coastline.sf,colour="black",size=1)+
  geom_sf(data = Counties.sf,colour="black",fill= "antiquewhite")+
  geom_sf(data=points.Sf,fill="yellow",colour="black",alpha=0.9,pch=21)+
  xlab("Longitude (°)")+ylab("Latitude (°)")+
  annotation_scale(location = "tl", width_hint = 0.5) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(2900000 ,3400000),
           ylim = c(3550000,4100000))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = gray(.5),linetype = "dashed", size = 0.5),
        panel.margin = unit(0, "lines"),
        panel.background = element_rect(fill = "aliceblue"),
        axis.title.x = element_text(size=rel(1.4), face="bold",
                                    margin = margin(t = 25, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=rel(1.4), face="bold",
                                    margin = margin(t = 0, r = 25, b = 0, l = 0)))


