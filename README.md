# Brown Crab_IRL-NW :crab:

This repository contains the analysis carried for the stock assessment of Brown Crab (*Cancer pagurus*) in the North West of Ireland. 

The scripts include the process of: 
1. Data collation / cleaning 

2. Development of two independent standarized indices of abundance: 
     * <ins>1_Prediction_Grid.R</ins>: script to create a grid (1x1 degree resolution) in which to make the predictions of the INLA model to obtain the standarized index of abundance. 
    * <ins>2_I4_INLA_AR.R</ins>: shows the development of an standarized LPUE index of abundance using the Integrated Nested Laplace Approximation (INLA). Specifically, an autoregressive spatio-temporal model is implemented. The construction of the model and the development of the standarized index of abundance is based on the work by  https://doi.org/10.1093/icesjms/fsz034.
    * <ins>3_CRE_NW_GAM.R</ins>: shows the standarization of the LPUE index from the Sentinel Vessel Fleet programme (a self-sampling commercial scheme) using Generalized Additive Models (GAM's). 
 
3. <ins>5_CRE_NW_SPiCT_Scenarios.R</ins>: Implementation of different scenario of the Stochastic Surplus Production Model in Continuous Time (SPiCT) (https://github.com/DTUAqua/spict)

Several scripts still in progress
