


################################################################################
Transf_Gauss_lognorm = function(E, stdev){
################################################################################
  #^* GOAL: provides mean and standard deviation of the Log-Normal distribution
  #^*       given mean and standard deviation of Normal distribution
  #^****************************************************************************
  #^* PROGRAMMER
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [real] E --> mean of the Gaussian distribution
  #^*    2. [real] stdev --> stand.dev. of the Gaussian distribution
  #^* OUT
  #^*    list()
  #^*    1. [real] mu --> mean of the Lognormal distribution
  #^*    2. [real] sd --> standard deviation of the lognormal distribution
  #^****************************************************************************
  #^* REF.: 
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  #^****************************************************************************
  mu = log((E^2)/(stdev^2+E^2)^0.5)
  sd = (log(1 +(stdev^2)/(E^2)))^0.5
  return(list(mu = mu,
              sd = sd))
}






################################################################################
find_peaks <- function (x, m){
################################################################################
  #^* GOAL: find all most important peaks (local maximum values) in a time
  #^*       series within a time interval m
  #^****************************************************************************
  #^* PROGRAMMER:
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [array of reals] x --> observations series
  #^*    2. [integer] m --> interval where to find maximum value
  #^* OUT
  #^*    1. [real] pks --> peaks values
  #^****************************************************************************
  #^* REF.: 
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  #^****************************************************************************  
  shape <- diff(sign(diff(x, na.pad=F)))
  pks <- sapply(which(shape<0), FUN=function(i){
    z <- i-m+1
    z <- ifelse(z>0,z,1)
    w <- i+m+1
    w <- ifelse(w<length(x),w,length(x))
    if(all(x[c(z:i, (i+2):w)] <= x[i+1])){
      return(i+1) 
    } else {
      return(numeric(0))
    }
  })
  pks<-unlist(pks)
  pks
}






