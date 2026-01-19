


################################################################################
logPrior = function(theta, priors){
################################################################################
  #^* GOAL: Log-Prior (full nonstationary weibull model SMEV)                  #
  #^***************************************************************************#
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy                  #
  #^***************************************************************************#
  #^* CREATED/MODIFIED: 2025/11/28                                             #
  #^***************************************************************************#
  #^* IN                                                                       #
  #^*   1. [array of reals] theta --> set of SMEV nonstat params for posterior # 
  #^*   2. [list of lists] priors --> list of priors for SMEV parameters,      #
  #^*                                 independently if stat or nonstat model   # 
  #^*                                 (shape intercept, shape slope,           #
  #^*                                 scale intercpet, scale slope).           #
  #^*                                 one list for each parameter, including   #
  #^*                                 distrib.name,mean,stdev,init value.      #
  #^* OUT                                                                      #
  #^*    1. [real] out --> value of log-prior                                  #
  #^****************************************************************************
  #^* REF.: Marra et al., 2020
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************
  # Parameters of non-stationary Weibull distribution:
  shape_intercept =theta[1]
  shape_slope     =theta[2]
  scale_intercept =theta[3]
  scale_slope     =theta[4]
  #set.seed(1)
  
  # shape intercept:
  if (priors$shape_intercept[[1]]=='uniform'){
    d1=dunif(x=shape_intercept, 
             min=priors$shape_intercept[[2]],
             max=priors$shape_intercept[[3]],log=T)
  } else if (priors$shape_intercept[[1]]=='normal'){
    d1=dnorm(x=shape_intercept, 
             mean=priors$shape_intercept[[2]], 
             sd=priors$shape_intercept[[3]],log=T)
  } else if (priors$shape_intercept[[1]]=='lognormal'){
    d1=dlnorm(x=shape_intercept, 
              meanlog=priors$shape_intercept[[2]],
              sdlog=priors$shape_intercept[[3]], 
              log=T)
  } else if (priors$shape_intercept[[1]]=='flatprior+'){
    d1=0
  } else {
    stop("prior distribution not available for shape intercept !!!")
  }
  
  # shape slope:
  if (priors$shape_slope[[1]]=='uniform'){
    d2=dunif(x=shape_slope,
             min=priors$shape_slope[[2]], 
             max=priors$shape_slope[[3]],log=T)
  } else if (priors$shape_slope[[1]]=='normal'){
    d2=dnorm(x=shape_slope, 
             mean=priors$shape_slope[[2]], 
             sd=priors$shape_slope[[3]],log=T)
  } else if (priors$shape_slope[[1]]=='lognormal'){
    d2=dlnorm(x=shape_slope, 
              meanlog=round(Transf_Gauss_lognorm(E=priors$shape_slope[[2]], 
                                                 stdev=priors$shape_slope[[3]])$mu,digits=4),
              sdlog=round(Transf_Gauss_lognorm(E=priors$shape_slope[[2]],
                                               stdev=priors$shape_slope[[3]])$sd,digits=4),
              log=T)
  } else if (priors$shape_slope[[1]]=='flatprior+'){
    d2=0
  } else {
    stop("prior distribution not available for shape intercept !!!")
  }
  
  # scale intercept:
  if (priors$scale_intercept[[1]]=='uniform'){
    d3=dunif(x=scale_intercept, 
             min=priors$scale_intercept[[2]], 
             max=priors$scale_intercept[[3]],log=T)
  } else if (priors$scale_intercept[[1]]=='normal'){
    d3=dnorm(x=scale_intercept, 
             mean=priors$scale_intercept[[2]], 
             sd=priors$scale_intercept[[3]], 
             log=T)
  } else if (priors$scale_intercept[[1]]=='lognormal'){
    d3=dlnorm(x=scale_intercept, 
              meanlog=priors$scale_intercept[[2]],
              sdlog=priors$scale_intercept[[3]],
              log=T)
  } else if (priors$scale_intercept[[1]]=='flatprior+'){
    d3=0
  } else {
    stop("prior distribution not available for shape intercept !!!")
  }
  
  # scale slope:
  if (priors$scale_slope[[1]]=='uniform'){
    d4=dunif(x=scale_slope, 
             min=priors$scale_slope[[2]],
             max=priors$scale_slope[[3]],log=T)
  } else if (priors$scale_slope[[1]]=='normal'){
    d4=dnorm(x=scale_slope, 
             mean=priors$scale_slope[[2]], 
             sd=priors$scale_slope[[3]],log=T)
  } else if (priors$scale_slope[[1]]=='lognormal'){
    d4=dlnorm(x=scale_slope, 
              meanlog=priors$scale_intercept[[2]],
              sdlog=priors$scale_intercept[[3]],
              log=T)
  } else if (priors$scale_slope[[1]]=='flatprior+'){
    d4=0
  } else {
    stop("prior distribution not available for shape intercept !!!")
  }
  
  # Compute Prior:
  out=d1+d2+d3+d4
  return(out)
}










################################################################################
logPrior_partns = function(theta, priors){
################################################################################
  #^* GOAL: Log-Prior (partial nonstationary weibull model SMEV)
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
  #^* REF.:  Marra et al., 2020
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************
  # Parameters of Partial non-stationary Weibull distribution:
  shape_intercept=theta[1]
  scale_intercept=theta[2]
  scale_slope=theta[3]
  #set.seed(1)
  
  # Shape intercept:
  if (priors$shape_intercept[[1]]=='Unif'){
    d1=dunif(x=shape_intercept, 
             min=priors$shape_intercept[[2]], 
             max=priors$shape_intercept[[3]],log=T)
  } else if (priors$shape_intercept[[1]]=='Gaussian'){
    d1=dnorm(x=shape_intercept, 
             mean=priors$shape_intercept[[2]], 
             sd=priors$shape_intercept[[3]],log=T)
  } else if (priors$shape_intercept[[1]]=='LogNormal'){
    d1=dlnorm(x=shape_intercept, 
              meanlog=round(Transf_Gauss_lognorm(E=priors$shape_intercept[[2]], 
                                                 stdev=priors$shape_intercept[[3]])$mu,digits=4),
              sdlog=round(Transf_Gauss_lognorm(E=priors$shape_intercept[[2]], 
                                               stdev=priors$shape_intercept[[3]])$sd,digits=4),
              log=T)
  } else if (priors$shape_intercept[[1]]=='FlatPrior+'){
    d1=0
  } else {
    stop("prior distribution not available for shape intercept !!!")
  }
  
  
  # scale intercept:
  if (priors$scale_intercept[[1]]=='Unif'){
    d3=dunif(x=scale_intercept,
             min=priors$scale_intercept[[2]],
             max=priors$scale_intercept[[3]],log=T)
  } else if (priors$scale_intercept[[1]]=='Gaussian'){
    d3=dnorm(x=scale_intercept, 
             mean=priors$scale_intercept[[2]],
             sd=priors$scale_intercept[[3]],log=T)
  } else if (priors$scale_intercept[[1]]=='LogNormal'){
    d3=dlnorm(x=scale_intercept, 
              meanlog=round(Transf_Gauss_lognorm(E=priors$scale_intercept[[2]],
                                                 stdev=priors$scale_intercept[[3]])$mu,digits=4),
              sdlog = round(Transf_Gauss_lognorm(E=priors$scale_intercept[[2]],
                                                 stdev= priors$scale_intercept[[3]])$sd,digits=4),
              log=T)
  } else if (priors$scale_intercept[[1]]=='FlatPrior+'){
    d3=0
  } else {
    stop("prior distribution not available for scale intercept !!!")
  }
  
  
  # Scale slope:
  if (priors$scale_slope[[1]]=='Unif'){
    d4=dunif(x=scale_slope, 
             min=priors$scale_slope[[2]], 
             max=priors$scale_slope[[3]],log=T)
  } else if (priors$scale_slope[[1]]=='Gaussian'){
    d4=dnorm(x=scale_slope, 
             mean=priors$scale_slope[[2]], 
             sd=priors$scale_slope[[3]],log=T)
  } else if (priors$scale_slope[[1]]=='LogNormal'){
    d4=dlnorm(x=scale_slope, 
              meanlog=round(Transf_Gauss_lognorm(E=priors$scale_slope[[2]], 
                                                 stdev=priors$scale_slope[[3]])$mu,digits=4),
              sdlog=round(Transf_Gauss_lognorm(E=priors$scale_slope[[2]], 
                                               stdev=priors$scale_slope[[3]])$sd,digits=4),
              log=T)
  } else if (priors$scale_slope[[1]]=='FlatPrior+'){
    d4=0
  } else {
    stop("prior distribution not available for scale slope intercept !!!")
  }
  
  # Compute Prior:
  out=d1+d3+d4
  return(out)
}











################################################################################
logPrior_stat = function(theta, priors){
################################################################################
  #^* GOAL: Log-Prior (stationary weibull model SMEV)
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
  #^* REF.:  Marra et al., 2020
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************
  # Parameters stationary Weibull distribution:
  shape_intercept = theta[1]
  scale_intercept = theta[2]
  
  # Shape intercept:
  if (priors$shape_intercept[[1]]=='Unif'){
    d1=dunif(x=shape_intercept, 
             min=priors$shape_intercept[[2]], 
             max=priors$shape_intercept[[3]],log=T)
  } else if (priors$shape_intercept[[1]]=='Gaussian'){
    d1=dnorm(x=shape_intercept, 
             mean=priors$shape_intercept[[2]], 
             sd=priors$shape_intercept[[3]],log=T)
  } else if (priors$shape_intercept[[1]]=='LogNormal'){
    d1=dlnorm(x=shape_intercept,
              meanlog=round(Transf_Gauss_lognorm(E=priors$shape_intercept[[2]], 
                                                 stdev=priors$shape_intercept[[3]])$mu,digits=4),
              sdlog=round(Transf_Gauss_lognorm(E=priors$shape_intercept[[2]], 
                                               stdev=priors$shape_intercept[[3]])$sd,digits=4),
              log=T)
  } else if (priors$shape_intercept[[1]]=='FlatPrior+'){
    d1=0
  } else {
    stop("prior distribution not available for shape intercept !!!")
  }
  
  
  # Scale intercept:
  if (priors$scale_intercept[[1]]=='Unif'){
    d3=dunif(x=scale_intercept, 
             min=priors$scale_intercept[[2]], 
             max=priors$scale_intercept[[3]],log=T)
  } else if (priors$scale_intercept[[1]]=='Gaussian'){
    d3=dnorm(x=scale_intercept,
             mean=priors$scale_intercept[[2]], 
             sd=priors$scale_intercept[[3]],log=T)
  } else if (priors$scale_intercept[[1]]=='LogNormal'){
    d3=dlnorm(x=scale_intercept, 
              meanlog=round(Transf_Gauss_lognorm(E=priors$scale_intercept[[2]],
                                                 stdev=priors$scale_intercept[[3]])$mu,digits=4),
              sdlog=round(Transf_Gauss_lognorm(E=priors$scale_intercept[[2]],
                                               stdev=priors$scale_intercept[[3]])$sd,digits=4),
              log=T)
  } else if (priors$scale_intercept[[1]]=='FlatPrior+'){
    d3=0
  } else {
    stop("prior distribution not available for scale intercept !!!")
  }
  
  # Compute Prior:
  out=d1+d3
  return(out)
}








################################################################################
logPrior_stat_gev = function(theta, priors_gev){
################################################################################
  #^* GOAL: Log-Prior (stationary GEV model)
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
  #^* REF.: 
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************
  # gev parameters (location, scale, shape):
  loc=theta[1]
  scal=theta[2]
  shap=theta[3]
    
  # location:
  if (priors_gev$loc[[1]]=='uniform'){
    d1=dunif(x=loc, 
             min=priors_gev$loc[[2]], 
             max=priors_gev$loc[[3]],log=T)
  } else if (priors_gev$loc[[1]]=='normal'){
    d1=dnorm(x=loc, 
             mean=priors_gev$loc[[2]], 
             sd=priors_gev$loc[[3]],log=T)
  } else if (priors_gev$loc[[1]]=='lognormal'){
    d1=dlnorm(x=loc,
              meanlog=round(Transf_Gauss_lognorm(E=priors_gev$loc[[2]],
                                                 stdev=priors_gev$loc[[3]])$mu,digits=4),
              sdlog=round(Transf_Gauss_lognorm(E=priors_gev$loc[[2]], 
                                               stdev=priors_gev$loc[[3]])$sd,digits=4),
              log=T)
  } else if (priors_gev$loc[[1]]=='flatprior+'){
    d1=0
  } else {
    stop("prior distribution not available for location !!!")
  }
  
  # Scale:
  if (priors_gev$scale[[1]]=='uniform'){
    d2=dunif(x=scal,
             min=priors_gev$scale[[2]],
             max=priors_gev$scale[[3]],log=T)
  } else if (priors_gev$scale[[1]]=='normal'){
    d2=dnorm(x=scal,
             mean=priors_gev$scale[[2]],
             sd=priors_gev$scale[[3]],log=T)
  } else if (priors_gev$scale[[1]]=='lognormal'){
    d2=dlnorm(x=scal, 
              meanlog=round(Transf_Gauss_lognorm(E=priors_gev$scale[[2]],
                                                 stdev=priors_gev$scale[[3]])$mu,digits=4),
              sdlog=round(Transf_Gauss_lognorm(E=priors_gev$scale[[2]], 
                                               stdev=priors_gev$scale[[3]])$sd,digits=4),
              log=T)
  } else if (priors_gev$scale[[1]]=='flatprior+'){
    d2=0
  } else {
    stop("prior distribution not available for location !!!")
  }
  
  # shape:
  if (priors_gev$shape[[1]]=='uniform'){
    d3=dunif(x=shap,
             min=priors_gev$shape[[2]], 
             max=priors_gev$shape[[3]],log=T)
  } else if (priors_gev$shape[[1]]=='normal'){
    d3=dnorm(x=shap, 
             mean=priors_gev$shape[[2]], 
             sd=priors_gev$shape[[3]],log=T)
  } else if (priors_gev$shape[[1]]=='lognormal'){
    d3=dlnorm(x=shap, 
                meanlog=round(Transf_Gauss_lognorm(E=priors_gev$shape[[2]],
                                                   stdev=priors_gev$shape[[3]])$mu,digits=4),
                sdlog=round(Transf_Gauss_lognorm(E=priors_gev$shape[[2]],
                                                 stdev=priors_gev$shape[[3]])$sd,digits=4),
                log=T)
  } else if (priors_gev$shape[[1]]=='flatprior+'){
    d3=0
  } else {
    stop("prior distribution not available for shape intercept !!!")
  }
  
  # Compute Prior:
  out=d1+d2+d3
  return(out)
}

