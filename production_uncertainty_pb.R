##---------------------------
## Function to plot uncertainty in estimated
## production function and empirical raised
## surplus production values on top
## CM, PB: 18/02/2020
##---------------------------

library(mvtnorm)
library (spict)

plot_production <- function(fit, nsamp = 1e4, nB = 100, plot_it = TRUE){
  ## fit is a spict fit object
  ## nsamp is the number of replicate samples from m and K and possibly n
  ## nB is the number of biomass points to calculate production at
  if("logn" %in% names(fit$par.fixed)){
    pars <- c("logm", "logK", "logn")
  }else{
    pars <- c("logm", "logK")
  }
  mu <- fit$par.fixed[pars]
  S <- fit$cov.fixed[pars, pars]
  mhat <- exp(mu["logm"])
  Khat <- exp(mu["logK"])
  if("logn" %in% names(fit$par.fixed)){
    nhat <- exp(mu["logn"])
  }else{
    nhat <- 2
  }
  gammahat <- nhat^(nhat/(nhat-1))/(nhat-1)
  ## get the upper bounds on the parameters for use later
  summ <- sumspict.parest(fit, ndigits = 7)
  mupr <- summ["m", "ciupp"]
  kest <- summ["K", "estimate"]
  Kupr <- summ["K", "ciupp"]
  Klwr <- summ["K", "cilow"]
  ## generate samples
  samp <- as.data.frame(rmvnorm(nsamp, mu, S))
  samp$m <- exp(samp$logm)
  samp$K <- exp(samp$logK)
  if("logn" %in% names(fit$par.fixed)){
    samp$n <- exp(samp$logn)
  }else{
    samp$n <- 2
  }
  samp$gamma <- with(samp, n^(n/(n-1))/(n-1))
  ##with(samp, plot(K, m, col = "grey"))
  ##points(Khat, mhat, pch = 19, col = "red")
  ## calculate the empirical surplus production for each index
  qhat <- summ[grep("^(q[0-9]|q)", rownames(summ)), "estimate"]
  ni <- length(fit$inp$obsI)## number of indices
  C_df <- data.frame(year = fit$inp$timeC,
                     C = fit$inp$obsC)
  sp_df <- NULL
  for(j in 1:ni){
    I_df <- data.frame(year = floor(fit$inp$timeI[[j]]),
                       I = fit$inp$obsI[[j]],
                       index = j)
    I_df$Bt <- I_df$I / qhat[j]
    ## next year biomass
    I_df$Bt1 <- c(I_df$Bt[-1], NA)
    tmp <- merge(I_df, C_df)
    tmp$sp <- with(tmp, Bt1 - Bt + C)
    sp_df <- rbind(sp_df, tmp)
  }
  
  ylim <- range(sp_df$sp, na.rm = TRUE)
  # max sample is either 1.1 times max sp from index, or is k
  maxx <- ifelse(max(sp_df$Bt, na.rm = TRUE) < kest, kest, max(sp_df$Bt, na.rm = TRUE))
  xlim <- c(0, 1.3 * maxx)
  
  
  ##### pb edit
  ## was sampling far too high for some stocks (WhgNS, kupper 8.412696e+10)
  B <- seq(0, xlim[2], length = nB)
  
  sp_mat <- matrix(NA, nrow = nsamp, ncol = nB)
  for(i in 1:nsamp){   
    sp_mat[i, ] <- samp$gamma[i] * samp$m[i] * (B/samp$K[i]) - samp$gamma[i] * samp$m[i] * (B/samp$K[i])^samp$n[i]
  }
  ## 95% CI
  ci <- apply(sp_mat, 2, quantile, p = c(0.025, 0.975))
  ## median from estimates
  P <- gammahat * mhat * (B/Khat) - gammahat * mhat * (B/Khat)^nhat
  ##sp_mean <- apply(sp_mat, 2, mean)
  sp_median <- apply(sp_mat, 2, median)
  
  # plot
  if(plot_it){
    matplot(B, t(ci), type = "n", lty = 2, col = 1, ylim = ylim, xlim = xlim, 
            bty = "l", xlab = "Biomass", ylab = "Surplus production", xaxs = "i")
    polygon(c(B, rev(B)), c(ci[1,], rev(ci[2,])), border = NA, col = "grey")
    lines(B, P)
    ##points(B, sp_mean)
    ##points(B, sp_median, col = "red") ## slight bias here
    abline(v = c(Klwr, Kupr), lty = 2)
    abline(h = 0, lty = 2)
    ## add in empirical surplus production
    for(j in 1:ni){
      tmp <- subset(sp_df, index == j)
      ##with(tmp, lines(Bt, sp, type = "o", pch = j, lwd = 0.5, cex = 0.3))
      with(tmp, points(Bt, sp, pch = j, lwd = 0.5, cex = 1))
    }
  }    
  ## output
  res <- list()
  ## predicted production
  tmp <- as.data.frame(cbind(B, P, t(ci)))
  names(tmp)[3:4] <- c("Plwr", "Pupr")    
  res$pred <- tmp
  ## empirical SP
  res$emp <- sp_df
  ## K bounds
  res$Kbounds <- c(kest, Klwr, Kupr)
  return(res)
}