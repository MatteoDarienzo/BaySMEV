



################################################################################
logPost=function(theta, y, t, thr, priors){
################################################################################
  #^* GOAL: Log-Posterior of full nonstationary SMEV model (without constant)  #
  #^***************************************************************************#
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy                  #
  #^***************************************************************************#
  #^* CREATED/MODIFIED: 2025/11/28                                             #
  #^***************************************************************************#
  #^* IN                                                                       #
  #^*   1. [array of reals] theta --> set of SMEV nonstat params for posterior # 
  #^*   2. [array of reals] y --> array or ordinary events values              #
  #^*   3. [array of reals] t --> array of ordinary events years (from 0)      #
  #^*   4. [real] thr -->  threshold percentile for left censoring (0.9)       #
  #^*   5. [list of lists] priors --> list of priors for SMEV parameters,      #
  #^*                                 independently if stat or nonstat model   # 
  #^*                                 (shape intercept, shape slope,           #
  #^*                                 scale intercpet, scale slope).           #
  #^*                                 one list for each parameter, including   #
  #^*                                 distrib.name,mean,stdev,init value.      #
  #^* OUT                                                                      #
  #^*    1. [real] out --> value of log-posterior                              #
  #^***************************************************************************#
  #^* REF.:                                                                    #
  #^***************************************************************************#
  #^* to do:                                                                   #
  #^***************************************************************************#
  #^* COMMENTS:                                                                #
  #^***************************************************************************#
  out= logLike(theta, y, t, thr) + logPrior(theta, priors)
  return(out)
}








################################################################################
logPost_partns=function(theta, y, t, thr, priors){
################################################################################
  #^* GOAL: Log-Posterior of partial nonstationary SMEV model (without constant)
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN                                                                       #
  #^*   1. [array of reals] theta --> SMEV part.nonstat params for posterior   # 
  #^*   2. [array of reals] y --> array or ordinary events values [mm]         #
  #^*   3. [array of reals] t --> array of ordinary events years (from 0)      #
  #^*   4. [real] thr -->  threshold percentile for left censoring (0.9)       #
  #^*   5. [list of lists] priors --> list of priors for SMEV parameters,      #
  #^*                                 independently if stat or nonstat model   # 
  #^*                                 (shape intercept, shape slope,           #
  #^*                                 scale intercpet, scale slope).           #
  #^*                                 one list for each parameter, including   #
  #^*                                 distrib.name,mean,stdev,init value.      #
  #^* OUT                                                                      #
  #^*    1. [real] out --> value of log-posterior                              #
  #^***************************************************************************#
  #^* REF.:                                                                    #
  #^***************************************************************************#
  #^* to do:                                                                   #
  #^***************************************************************************#
  #^* COMMENTS:                                                                #
  #^***************************************************************************#
  out= logLike_partns(theta, y, t, thr) + logPrior_partns(theta, priors)
  return(out)
}







################################################################################
logPost_stat=function(theta, y, t, thr, priors){
################################################################################
  #^* GOAL: Log-Posterior of stationary SMEV model (without constant)          #
  #^***************************************************************************#
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy                  #
  #^***************************************************************************#
  #^* CREATED/MODIFIED: 2025/11/28                                             #
  #^***************************************************************************#
  #^* IN                                                                       #
  #^*   1. [array of reals] theta --> set of SMEV stat params for posterior    # 
  #^*   2. [array of reals] y --> array or ordinary events values [mm]         #
  #^*   3. [array of reals] t --> array of ordinary events years (from 0)      #
  #^*   4. [real] thr -->  threshold percentile for left censoring (0.9)       #
  #^*   5. [list of lists] priors --> list of priors for SMEV parameters,      #
  #^*                                 independently if stat or nonstat model   # 
  #^*                                 (shape intercept, shape slope,           #
  #^*                                 scale intercpet, scale slope).           #
  #^*                                 one list for each parameter, including   #
  #^*                                 distrib.name,mean,stdev,init value.      #
  #^* OUT                                                                      #
  #^*    1. [real] out --> value of log-posterior                              #
  #^***************************************************************************#
  #^* REF.:                                                                    #
  #^***************************************************************************#
  #^* to do:                                                                   #
  #^***************************************************************************#
  #^* COMMENTS:                                                                #
  #^***************************************************************************#
  out = logLike_stat(theta, y, t, thr) + logPrior_stat(theta, priors)
  return(out)
}








################################################################################
logPost_stat_gev=function(theta, y, priors_gev){
################################################################################
  #^* GOAL: Log-Posterior of stationary GEV model (without constant)           #
  #^***************************************************************************#
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy                  #
  #^***************************************************************************#
  #^* CREATED/MODIFIED: 2025/11/28                                             #
  #^***************************************************************************#
  #^* IN                                                                       #
  #^*   1. [array of reals] theta --> set of GEV stat params for posterior     # 
  #^*   2. [array of reals] y --> array or annual maximum values               #
  #^*   3. [list of lists] priors_gev --> list of priors for GEV parameters,   #
  #^*                                     (location, scale, shape):            #
  #^*                                   as distribname, mean, stdev, init value#
  #^* OUT                                                                      #
  #^*    1. [real] out --> value of log-posterior stat GEV                     #
  #^***************************************************************************#
  #^* REF.:                                                                    #
  #^***************************************************************************#
  #^* to do:                                                                   #
  #^***************************************************************************#
  #^* COMMENTS:                                                                #
  #^***************************************************************************#
  out = logLike_stat_gev(theta, y) + logPrior_stat_gev(theta, priors_gev)
  return(out)
}

