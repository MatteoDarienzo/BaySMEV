# needed packages
library(evd) # with functions for GEV distr


################################################################################
logLike=function(theta,y, t, thr){ 
################################################################################
  #^* GOAL: log-likelihood for full-nonstationary model, for mcmc algorithm
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [integer]
  #^*    2. [integer]
  #^*    3. [character]
  #^* OUT
  #^*    1. [real]
  #^****************************************************************************
  #^* REF.: takes hints from open source codes of Francesco Marra (University of
  #^*       Padua) for SMEV method and left-censoring 
  #^*       (https://zenodo.org/records/11934843) and non-stationary SMEV 
  #^*       (https://zenodo.org/records/15047817)
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************  
  # Non-stationary Weibull parameters (SMEV)
  a_w=theta[1]
  b_w=theta[2]
  a_C=theta[3]
  b_C=theta[4]
  t0=t[y<thr]
  shapes0=a_w+b_w*t0 # linear relation
  scales0=a_C+b_C*t0 # linear relation
  y1=y[y>=thr]
  t1=t[y>=thr]
  shapes1=a_w+b_w*t1 # linear relation
  scales1=a_C+b_C*t1 # linear relation
  
  # compute likelihood:
  #set.seed(1)
  if(any(scales1<=0)|any(scales0<=0)|any(shapes0<=0)|any(shapes1<=0)){
    # impossible values of the parameters
    out=-Inf
  } else { 
    #out=sum(log(wblcdf(thr,scales0,shapes0)))+sum(log(wblpdf(x1,scales1,shapes1)))
    out=sum(pweibull(q=thr,shape=shapes0,scale=scales0,log.p=T))+ 
        sum(dweibull(x=y1,shape=shapes1,scale=scales1,log=T))
  }
  return(out)
}






################################################################################
logLike_loo =function(thetat, y, t, thr){ 
################################################################################
  #^* GOAL: log-likelihood for full-nonstationary model, for LOO
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [array of reals] thetat
  #^*    2. [array of reals] y
  #^*    3. [array of reals]
  #^* OUT
  #^*    1. [real]
  #^****************************************************************************
  #^* REF.: takes hints from open source codes of Francesco Marra (University of
  #^*       Padua) for SMEV method and left-censoring 
  #^*       (https://zenodo.org/records/11934843) and non-stationary SMEV 
  #^*       (https://zenodo.org/records/15047817)
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************    
  # Non-stationary Weibull parameters (SMEV)
  a_w=thetat[1]
  b_w=thetat[2]
  a_C=thetat[3]
  b_C=thetat[4]
  out=c()
  for (d in 1:length(y)){
    shapes=a_w+b_w*t[d]
    scales=a_C+b_C*t[d]
    if (any(shapes<=0)|any(scales<=0)){
      out[d]=-Inf # impossible values of the parameters
    } else {
      if (y[d]>=thr){
        out[d]=dweibull(x=y[d],shape=shapes,scale=scales,log=T)
      } else {
        out[d]=pweibull(q=thr,shape=shapes,scale=scales,log.p=T)
      }
    }
  }
  return(out)
}







################################################################################
logLike1=function(x, y, t, thr){ 
################################################################################
  #^* GOAL: log-likelihood for full-nonstationary model, for fminsearch function 
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [integer]
  #^*    2. [integer]
  #^*    3. [character]
  #^* OUT
  #^*    1. [real]
  #^****************************************************************************
  #^* REF.: takes hints from open source codes of Francesco Marra (University of
  #^*       Padua) for SMEV method and left-censoring 
  #^*       (https://zenodo.org/records/11934843) and non-stationary SMEV 
  #^*       (https://zenodo.org/records/15047817)
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************  
  # Non-stationary Weibull parameters (SMEV)
  a_w=x[1]
  b_w=x[2]
  a_C=x[3]
  b_C=x[4]
  t0=t[y<thr]
  shapes0=a_w+b_w*t0
  scales0=a_C+b_C*t0
  y1=y[y>=thr]
  t1=t[y>=thr]
  shapes1=a_w+b_w*t1
  scales1=a_C+b_C*t1
  if(any(scales1<=0)|any(scales0<=0)|any(shapes0<=0)|any(shapes1<=0)){
    out=-Inf
  } else { 
    out=sum(pweibull(q=thr, shape=shapes0, scale=scales0, log.p=T)) + 
        sum(dweibull(x=y1,  shape=shapes1, scale=scales1, log=T))
  }
  return(-out)
}









################################################################################
logLike_partns=function(theta, y, t, thr){ 
################################################################################
  #^* GOAL: log-likelihood for Partial-nonstationary model, for mcmc algorithm
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [integer]
  #^*    2. [integer]
  #^*    3. [character]
  #^* OUT
  #^*    1. [real]
  #^****************************************************************************
  #^* REF.: takes hints from open source codes of Francesco Marra (University of
  #^*       Padua) for SMEV method and left-censoring 
  #^*       (https://zenodo.org/records/11934843) and non-stationary SMEV 
  #^*       (https://zenodo.org/records/15047817)
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************  
  # Partial non-stationary Weibull parameters (SMEV)  
  a_w=theta[1]
  a_C=theta[2]
  b_C=theta[3]
  t0=t[y<thr]
  shapes0=a_w*rep(1,length(t0))
  scales0=a_C+b_C*t0
  y1=y[y>=thr]
  t1=t[y>=thr]
  shapes1=a_w*rep(1,length(t1))
  scales1=a_C+b_C*t1
  #set.seed(1)
  if(any(scales1<=0)|any(scales0<=0)|any(shapes0<=0)|any(shapes1<=0)){
    # impossible values of the parameters
    out=-Inf
  } else { 
    # out=sum(log(wblcdf(thr,scales0,shapes0)))+sum(log(wblpdf(x1,scales1,shapes1)))
    out=sum(pweibull(q=thr, shape=shapes0, scale=scales0, log.p=T)) + 
        sum(dweibull(x=y1,  shape=shapes1, scale=scales1, log=T))
  }
  return(out)
}








################################################################################
logLike1_partns =function(x, y, t, thr){
################################################################################
  #^* GOAL: Log-Likelihood (Partial non-stationary) 
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [array of reals] x --> vector of 3 Weib param (shape int, scale int,
  #^*                              scale slope) 
  #^*    2. [array of reals] y --> array of observations 
  #^*    3. [array of reals] t --> array of observational times (years)
  #^*    4. [real] thr --> left-censoring threshold value (e.g., 90% quantile)
  #^* OUT
  #^*    1. [real] out --> Log-Likelihood (changed of sign) for fminsearch() 
  #^****************************************************************************
  #^* REF.: takes hints from open source codes of Francesco Marra (University of
  #^*       Padua) for SMEV method and left-censoring 
  #^*       (https://zenodo.org/records/11934843) and non-stationary SMEV 
  #^*       (https://zenodo.org/records/15047817)
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  #^****************************************************************************  
  # log-likelihood for partial-nonstationary model, for fminsearch function (MLE) 
  a_w=x[1]
  a_C=x[2]
  b_C=x[3]
  t0=t[y<thr]
  shapes0=a_w*rep(1,length(t0))
  scales0=a_C+b_C*t0
  y1=y[y>=thr]
  t1=t[y>=thr]
  shapes1=a_w*rep(1,length(t1))
  scales1=a_C+b_C*t1
  if(any(scales1<=0)|any(scales0<=0)|any(shapes0<=0)|any(shapes1<=0)){
    out=-Inf
  } else { 
    out=sum(pweibull(q=thr,shape=shapes0,scale=scales0,log.p=T))+ 
        sum(dweibull(x=y1,shape=shapes1,scale=scales1,log=T))
  }
  return(-out)
}






################################################################################
logLike_stat=function(theta, y, t, thr){
################################################################################
  #^* GOAL: Log-Likelihood (Stationary) 
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [array of reals] theta --> vector of 2 Weib parameters (shape int,
  #^*                                  scale int) 
  #^*    2. [array of reals] y --> array of observations 
  #^*    3. [array of reals] t --> array of observational times (years)
  #^*    4. [real] thr --> left-censoring threshold value (e.g., 90% quantile)
  #^* OUT
  #^*    1. [real] out --> Log-Likelihood
  #^****************************************************************************
  #^* REF.: takes hints from open source codes of Francesco Marra (University of
  #^*       Padua) for SMEV method and left-censoring 
  #^*       (https://zenodo.org/records/11934843) and non-stationary SMEV 
  #^*       (https://zenodo.org/records/15047817) 
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  #^****************************************************************************  
  shapes=theta[1]
  scales=theta[2]
  y1=y[y>=thr]
  if(any(scales<=0)|any(shapes<=0)){
    # impossible values of the parameters
    out=-Inf
  } else { 
    out=sum(y<thr)*pweibull(q=thr,shape=shapes,scale=scales,log.p=T)+ 
        sum(dweibull(x=y1,shape=shapes,scale=scales,log=T))
    # out=sum(pweibull(q=thr,shape=shapes,scale=scales,log.p=T))+ 
    #     sum(dweibull(x=y1,shape=shapes,scale=scales,log=T))
  }
  return(out)
}





################################################################################
logLike_stat_loo =function(thetat, y, t, thr){
################################################################################
  #^* GOAL: log-likelihood for stationary model, for loo package  
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [array of reals] thetat --> vector of 2 Weib parameters (shape int,
  #^*                                   scale int) 
  #^*    2. [array of reals] y --> array of observations 
  #^*    3. [array of reals] t --> array of observational times (years)
  #^*    4. [real] thr --> left-censoring threshold value (e.g., 90% quantile)
  #^* OUT
  #^*    1. [real] out --> Log-Likelihood for Loo
  #^****************************************************************************
  #^* REF.: takes hints from open source codes of Francesco Marra (University of
  #^*       Padua) for SMEV method and left-censoring 
  #^*       (https://zenodo.org/records/11934843) and non-stationary SMEV 
  #^*       (https://zenodo.org/records/15047817)
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  #^****************************************************************************  
  shapes=thetat[1]
  scales=thetat[2]
  out=c()
  if (any(shapes<=0) | any(scales<=0)){
    out[1:length(y)]=-Inf # impossible values of the parameters
  }
  for (d in 1:length(y)){
    if (y[d] >= thr){
      out[d]=dweibull(x=y[d],shape=shapes,scale=scales,log=T)
    } else {
      out[d]=pweibull(q=thr,shape=shapes,scale=scales,log.p=T)
    }
  }
  return(out)
}






################################################################################
logLike1_stat=function(x, y, t, thr){
################################################################################
  #^* GOAL: Log-Likelihood (Stationary SMEV)
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [array of reals] x --> vector of 2 Weib param (shape int, scale int) 
  #^*    2. [array of reals] y --> array of observations
  #^*    3. [array of reals] t --> array of observational times (years)
  #^*    4. [real] thr --> left-censoring threshold value (e.g., 90% quantile)
  #^* OUT
  #^*    1. [real] out --> Log-Likelihood (changed of sign) for fminsearch() 
  #^****************************************************************************
  #^* REF.: takes hints from open source codes of Francesco Marra (University of
  #^*       Padua) for SMEV method and left-censoring 
  #^*       (https://zenodo.org/records/11934843) and non-stationary SMEV 
  #^*       (https://zenodo.org/records/15047817)
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  #^****************************************************************************  
  # log-likelihood for stationary model, for fminsearch function  
  shapes=x[1]
  scales=x[2]
  y1=y[y>=thr]
  if(any(scales<=0)|any(shapes<=0)){
    out=-Inf # impossible values of the parameters
  } else {
    #out=sum(pweibull(q=thr,shape=shapes,scale=scales,log.p=T))+ 
    #    sum(dweibull(x=y1,shape=shapes,scale=scales,log=T))
    out=sum(y< thr)*pweibull(q=thr,shape=shapes,scale=scales,log.p=T)+
        sum(dweibull(x=y1,shape=shapes,scale=scales,log=T))
  }
  # return output with opposite sign
  return(-out)
}






################################################################################
logLike_stat_gev=function(theta, y){
################################################################################
  #^* GOAL: Log-Likelihood (Stationary) -- GEV
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [array of reals] theta --> vector of 3 GEV parameters (loc, scale,
  #^*                                  shape). 
  #^*    2. [array of reals] y --> array of observations (AMAX).
  #^* OUT
  #^*    1. [real] out --> Log-Likelihood
  #^****************************************************************************
  #^* REF.: use "evd" package to compute GEV
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  #^****************************************************************************  
  # For all i, Y_i ~ GEV(loc, scale, shape)
  loc=theta[1]
  scal=theta[2]
  shap=theta[3]
  if(scal<=0){
    out=-Inf
  } else { 
    out=sum(evd::dgev(y,loc=loc,scale=scal,shape=shap,log=T))
  }
  return(out)
}




################################################################################
logLike_stat_gev_loo =function(thetat, AMS){
################################################################################
  #^* GOAL: log-likelihood for stationary GEV model, for loo package
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [array of reals] thetat --> vector of 3 GEV parameters (loc, scale,
  #^*                                   shape). 
  #^*    2. [array of reals] AMS --> array of observations (AMAX).
  #^* OUT
  #^*    1. [real] out --> Log-Likelihood
  #^****************************************************************************
  #^* REF.: use "evd" package to compute GEV
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  #^****************************************************************************  
  # in this case "y" is not the array with ordinary events,
  # but it is the array with annual maxima AMS!!!
  # For all i,   AMS_i ~ GEV(loc, scale, shape)
  loc=thetat[1]
  scal=thetat[2]
  shap=thetat[3]
  out=c()
  if (scal<=0){
    out[1:length(AMS)]=-Inf  # Impossible values of the parameters!
  }
  for (d in 1:length(AMS)){
    out[d]=evd::dgev(AMS[d],loc=loc,scale=scal,shape=shap,log=T)
  }
  return(out)
}




