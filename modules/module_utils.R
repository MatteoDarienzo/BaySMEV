


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



