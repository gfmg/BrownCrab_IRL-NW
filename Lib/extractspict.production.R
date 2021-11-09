# File: extractspict.production.R
# Author: Guillermo Martin 
# Template Created: Tue Aug 17 09:20:49 2021
# ---------------------------------------------------------------------------
# Description:
# Function to extract the data used to plot the production curve in SPiCT 
# ---------------------------------------------------------------------------

extractspict.production <- function(rep,Scenario, main='Production curve',
                                     CI = 0.95){
  
  df<-NULL
  
 # namespace:spict::check.rep(rep)
  if (!'sderr' %in% names(rep)){
    inp <- rep$inp
    tvgflag <- rep$inp$timevaryinggrowth | rep$inp$logmcovflag
    Kest <- get.par('logK', rep, exp=TRUE, CI = CI)
    mest <- get.par('logm', rep, exp=TRUE, CI = CI)
    nr <- dim(mest)[1]
    gamma <- get.par('gamma', rep, CI = CI)
    n <- get.par('logn', rep, exp=TRUE, CI = CI)
    Pest <- get.par('P', rep, CI = CI)
    binds <- inp$ic[1:dim(Pest)[1]]
    Bmsy <- get.par('logBmsy', rep, exp=TRUE, CI = CI)
    Bmsy <- c(1,1)
    if (tvgflag){
      yscal <- get.par('logMSYvec', rep, exp=TRUE, CI = CI)[binds, 2]
    } else {
      yscal <- rep(1, length(binds))
    }
    nBplot <- 200
    Bplot <- seq(0.5*1e-8, Kest[2], length=nBplot)
    # Calculate production curve (Pst)
    pfun <- function(gamma, m, K, n, B) gamma*m/K*B*(1 - (B/K)^(n-1))
    Pst <- list()
    for (i in 1:nr){
      Pst[[i]] <- pfun(gamma[2], mest[i,2], Kest[2], n[2], Bplot)
    }
    Pstscal <- ifelse(tvgflag, max(unlist(Pst)), 1)
    ylim <- c(0, max(unlist(Pst)/Pstscal, na.rm=TRUE))
    if (inp$reportall){
      Best <- get.par('logB', rep, exp=TRUE, CI = CI)
      Bplot <- seq(0.5*min(c(1e-8, Best[, 2])), 1*max(c(Kest[2], Best[, 2])), length=nBplot)
      for (i in 1:nr){
        Pst[[i]] <- pfun(gamma[2], mest[i,2], Kest[2], n[2], Bplot)
      }
      
      Bvec <- Best[binds, 2]
      xlim <- range(Bvec/Kest[2], 0, 1)
      ylim <- c(min(0, Pest[,2]/yscal), max(Pest[,2]/yscal, unlist(Pst)/Pstscal, na.rm=TRUE))
    } else {
      xlim <- range(Bplot/Kest[2], na.rm=TRUE)
    }
    dt <- inp$dt[-1]
    inde <- inp$indest[-length(inp$indest)]
    indp <- inp$indpred[-1]-1
    ylab <- 'Production'
    if (tvgflag){
      ylab <- paste(ylab, '(normalised)')
    } else {
      ylab <- paste(ylab, inp$catchunit, sep=" ")
    }
  }
  df$Bplot<-Bplot
  df$Kest<-Kest[2]
  df$Pst<-Pst[[nr]]
  df$Pstscal<-Pstscal
  df$Scenario<-Scenario
  df$ylab<-ylab
  
  df <- data.frame(df)
  colnames(df) <- c("Bplot","Kest","Pst","Pstscal","Scenario","ylab")
  
  return(df)
}
