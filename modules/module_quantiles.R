library(evd) # with the GEV pdf, cdf, ...



################################################################################
weibull_inv <- function(X, scale, shape, n_ordev){
################################################################################
  #^* GOAL: inverts the 2 param Weibull to find quantiles at given prob 
  # target return periods 
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [real] X --> given probability
  #^*    2. [real] scale  --> Weibull scale parameter value
  #^*    3. [real] shape  --> Weibull shape parameter value
  #^*    4. [integer] n_ordev --> number of ordinary events per year
  #^* OUT
  #^*    1. [real] qnt --> quantile at the given probability X
  #^****************************************************************************
  #^* REF.:
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************
  # compute quantile of SMEV (Weibull):
  qnt=scale*(-log( 1-X^(1/n_ordev)))^(1/shape)
  return(qnt)
}






################################################################################
mevcdf <- function(x, theta, n_coeff, M, p){
################################################################################
  #^* GOAL: MEV cdf
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [real] x --> given probability
  #^*    2. [real] theta  --> Nonstationiry SMEV Weibull parameters values (4)
  #^*                         shape int, shape slope, scale int, scale slope
  #^*    3. [array of reals] n_coeff -->intercept and slope coeff of linear
  #^*                                   relation of "n" (number of events/year)
  #^*    4. [integer] M --> number of years in the studied period
  #^*    5. [real] p --> given probability
  #^* OUT
  #^*    1. [real] qnt --> quantile at the given probability X
  #^****************************************************************************
  #^* REF.: Marani & Ignaccolo 2015.
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************
  # ncount  years
  # ord_eve_13= y[which(t==13)]
  cdf=0; out=0
  # compute the mean over the period:
  for (year in 1:M){
    scale=theta[3]+theta[4]*year
    shape=theta[1]+theta[2]*year
    n=coeff[[1]]+coeff[[2]]*year
    cdf=cdf+(1-exp(-(x/scale)^shape))^n
  }
  cdf=cdf/M 
  out=cdf-p
  #cdf=mean((1-exp(-(y/phat[,1])^phat[,2]))^phat[,3])
  return(out)
}










################################################################################
compute_quantiles_SMEV<-function(s                = list(),
                                 estim_algorithm  = "cmdstan", 
                                 target_T         = c(2,5,10,50,100),
                                 flg_save         = T,
                                 id               = 1, 
                                 filter_spaghetti = 10,
                                 flag_smev_stat   = T,
                                 flag_smev_partns = F,
                                 flag_smev_ns     = F,
                                 flag_gev_stat    = F){
################################################################################
  #^* GOAL: Compute quantiles SMEV models for given set of return periods
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [list of lists] s --> object with all smev/GEV parameter estimates
  #^*    2. [character] estim_algorithm --> algorithm chosen for parameter estim:
  #^*                                 "cmdstan","stanr","adaptMCMC","fminsearch"
  #^*    3. [array of integers] target_T --> given return periods c(2,10,50,100)
  #^*    4. [logical] flg_save --> if T it saves all the input return periods, 
  #^*                              otherwise it computes quantiles only for a 
  #^*                              limited set (2,5,10,50,100 years).
  #^*    5. [integer] id --> duration index 
  #^*    6. [integer] filter_spaghetti --> filter spghetti to obtain a reduced
  #^*                                      ensemble for plotting reasons
  #^* OUT
  #^*    1. [list of lists] s_tmp --> list with all quantiles results
  #^****************************************************************************
  #^* REF.:
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  ##############################################################################
  # Compute quantiles curves with uncertainty:
  rm(s_tmp)
  s_tmp=c(); s_tmp$Tgrid=c()
  
  ##############################################################################
  if (flg_save==T){
  ##############################################################################
    # grid for prob:
    # pgrid=seq(0.01,0.99,0.1) 
    s_tmp$Tgrid=target_T
    # prob associated to the target_return_periods:
    s_tmp$pr=1-1/s_tmp$Tgrid  
    # grid of return periods
    # Tgrid=1/(1-pgrid)    
    
    #^**************************************************************************
    if (estim_algorithm=="cmdstan" | estim_algorithm=="stanr" | 
        estim_algorithm=="adaptMCMC"){
    #^**************************************************************************
      # initialization of arrays:
      if (flag_smev_ns){ 
        if (!is.null(s$s_ns_tmp[[id]]$zPost)){
          s_tmp$qCurve_Lik=matrix(NA,NROW(s$s_ns_tmp[[id]]$zLik),
                                length(s_tmp$Tgrid))
          s_tmp$qCurve_Post=matrix(NA,NROW(s$s_ns_tmp[[id]]$zPost),
                                 length(s_tmp$Tgrid))
          Nsim=NROW(s$s_ns_tmp[[id]]$zPost)
        } else {
          message("Posterior from Non-Stationary SMEV is missing!")
        }
      }
      if (flag_smev_partns){
        if (!is.null(s$s_partns_tmp[[id]]$zPost)){
          s_tmp$qCurve_Lik_partns=matrix(NA,NROW(s$s_partns_tmp[[id]]$zLik),
                                         length(s_tmp$Tgrid))
          s_tmp$qCurve_Post_partns=matrix(NA,NROW(s$s_partns_tmp[[id]]$zPost),
                                          length(s_tmp$Tgrid))
          Nsim=NROW(s$s_partns_tmp[[id]]$zPost)
        } else {
          message("Posterior from Partial-Non-Stationary SMEV is missing!")
        }
      }
      if (flag_smev_stat){
        if (!is.null(s$s_stat_tmp[[id]]$zPost)){
          s_tmp$qCurve_Lik_stat=matrix(NA,NROW(s$s_stat_tmp[[id]]$zLik),
                                       length(s_tmp$Tgrid))
          s_tmp$qCurve_Post_stat=matrix(NA,NROW(s$s_stat_tmp[[id]]$zPost),
                                        length(s_tmp$Tgrid))
          Nsim=NROW(s$s_stat_tmp[[id]]$zPost)
        } else {
          message("Posterior from Stationary SMEV is missing!")
        }
      }
      if (flag_gev_stat){
        if (!is.null(s$s_gev_tmp[[id]]$zPost)){
          s_tmp$qCurve_PostGEV=matrix(NA,NROW(s$s_gev_tmp[[id]]$zPost),
                                      length(s_tmp$Tgrid))
          Nsim=NROW(s$s_gev_tmp[[id]]$zPost)
        } else {
          message("Posterior sample from Stationary GEV is missing!")
        }
      }
      if (is.null(Nsim)){message("No models have been inferred. Please check!")}
      
      # All spaghetti:
      for(i in 1:Nsim){
        # NON-STATIONARY SMEV (Likelihood estimates):
        if (!is.null(s$s_ns_tmp[[id]]$zLik)){
          # compute shape at halfperiod year (with LIKE estimates):
          shape_t_Lik=s$s_ns_tmp[[id]]$zLik[i,1]+
            s$s_ns_tmp[[id]]$zLik[i,2]*s$halfperiod
          # compute scale at halfperiod year (with LIKE estimates):
          scale_t_Lik=s$s_ns_tmp[[id]]$zLik[i,3]+
            s$s_ns_tmp[[id]]$zLik[i,4]*s$halfperiod
          # quantiles at halfperiod year (with LIKE estimates):
          s_tmp$qCurve_Lik[i,]=sapply(scale=scale_t_Lik,
                                      shape=shape_t_Lik, 
                                      n_ordev=s$n_ordev[[id]], 
                                      FUN=weibull_inv,
                                      X=s_tmp$pr)
        }
        # NON-STATIONARY SMEV (Posterior estimates):
        if (!is.null(s$s_ns_tmp[[id]]$zPost)){
          # compute shape at halfperiod year (with POST estimates):
          shape_t_Post=s$s_ns_tmp[[id]]$zPost[i,1]+
            s$s_ns_tmp[[id]]$zPost[i,2]*s$halfperiod
          # compute scale at halfperiod year (with POST estimates):
          scale_t_Post=s$s_ns_tmp[[id]]$zPost[i,3]+
            s$s_ns_tmp[[id]]$zPost[i,4]*s$halfperiod
          # quantiles at halfperiod year (with POST estimates):
          s_tmp$qCurve_Post[i,]=sapply(scale=scale_t_Post, 
                                       shape=shape_t_Post,
                                       n_ordev=s$n_ordev[[id]],
                                       FUN=weibull_inv, 
                                       X=s_tmp$pr)
        }
        # PARTIAL-NON-STATIONARY SMEV (Likelihood estimates):
        if (!is.null(s$s_partns_tmp[[id]]$zLik)){
          # compute shape at halfperiod year (with LIKE estimates):
          shape_t_Lik_partns=s$s_partns_tmp[[id]]$zLik[i,1]
          # compute scale at halfperiod year (with LIKE estimates):
          scale_t_Lik_partns=s$s_partns_tmp[[id]]$zLik[i,2]+
            s$s_partns_tmp[[id]]$zLik[i,3]*s$halfperiod
          # quantiles at halfperiod year (with LIKE estimates):
          s_tmp$qCurve_Lik_partns[i,]=sapply(scale=scale_t_Lik_partns, 
                                             shape=shape_t_Lik_partns, 
                                             n_ordev=s$n_ordev[[id]], 
                                             FUN=weibull_inv, 
                                             X=s_tmp$pr)
        }
        # PARTIAL-NON-STATIONARY SMEV (Posterior estimates):
        if (!is.null(s$s_partns_tmp[[id]]$zPost)){
          # compute shape at halfperiod year (with POST estimates):
          shape_t_Post_partns=s$s_partns_tmp[[id]]$zPost[i,1]
          # compute shape at halfperiod year (with POST estimates):
          scale_t_Post_partns=s$s_partns_tmp[[id]]$zPost[i,2]+
            s$s_partns_tmp[[id]]$zPost[i,3]*s$halfperiod
          # quantiles at halfperiod year (with POST estimates):
          s_tmp$qCurve_Post_partns[i,]=sapply(scale=scale_t_Post_partns,
                                              shape=shape_t_Post_partns,
                                              n_ordev=s$n_ordev[[id]], 
                                              FUN=weibull_inv, 
                                              X=s_tmp$pr)
        }
        # STATIONARY SMEV (Posterior likelihood):
        if (!is.null(s$s_stat_tmp[[id]]$zLik)){
          shape_t_Lik_stat=s$s_stat_tmp[[id]]$zLik[i,1]
          scale_t_Lik_stat=s$s_stat_tmp[[id]]$zLik[i,2]
          s_tmp$qCurve_Lik_stat[i,]=sapply(scale=scale_t_Lik_stat, 
                                           shape=shape_t_Lik_stat, 
                                           n_ordev=s$n_ordev[[id]],
                                           FUN=weibull_inv, 
                                           X=s_tmp$pr)
        }
        # STATIONARY SMEV (Posterior estimates):
        if (!is.null(s$s_stat_tmp[[id]]$zPost)){
          shape_t_Post_stat=s$s_stat_tmp[[id]]$zPost[i,1]
          scale_t_Post_stat=s$s_stat_tmp[[id]]$zPost[i,2]
          s_tmp$qCurve_Post_stat[i,]=sapply(scale=scale_t_Post_stat,
                                            shape=shape_t_Post_stat,
                                            n_ordev=s$n_ordev[[id]], 
                                            FUN=weibull_inv, 
                                            X=s_tmp$pr)
        }
        # STATIONARY GEV (Posterior estimates):
        if (!is.null(s$s_gev_tmp[[id]]$zPost)){
          loc_t_PostGEV=s$s_gev_tmp[[id]]$zPost[i,1]
          shape_t_PostGEV=s$s_gev_tmp[[id]]$zPost[i,3]
          scale_t_PostGEV=s$s_gev_tmp[[id]]$zPost[i,2]
          s_tmp$qCurve_PostGEV[i,]=sapply(loc=loc_t_PostGEV, 
                                          scale=scale_t_PostGEV, 
                                          shape=shape_t_PostGEV, 
                                          FUN=qgev,
                                          X=s_tmp$pr)
        }
      }
      
      
      #############
      # MAP and CI:
      #############
      
      #^************************************************************************
      #                     NONSTATIONARY (MAXPOST):
      #^************************************************************************
      if (!is.null(s$s_ns_tmp[[id]]$theta_maxpost)){
        s_tmp$shape_t_maxpost_halfper=s$s_ns_tmp[[id]]$theta_maxpost[1]+
          s$s_ns_tmp[[id]]$theta_maxpost[2]*s$halfperiod
        s_tmp$scale_t_maxpost_halfper=s$s_ns_tmp[[id]]$theta_maxpost[3]+
          s$s_ns_tmp[[id]]$theta_maxpost[4]*s$halfperiod
        
        # s_tmp$shape_t_maxpost_0y=s$s_ns_tmp[[id]]$theta_maxpost[1]   # firs year
        # s_tmp$scale_t_maxpost_0y=s$s_ns_tmp[[id]]$theta_maxpost[3]   # firs year
        # s_tmp$shape_t_maxpost_10y=s$s_ns_tmp[[id]]$theta_maxpost[1]+ 
        #   s$s_ns_tmp[[id]]$theta_maxpost[2]*10                       # decade
        # s_tmp$scale_t_maxpost_10y=s$s_ns_tmp[[id]]$theta_maxpost[3]+
        #   s$s_ns_tmp[[id]]$theta_maxpost[4]*10                       # decade
        # s_tmp$shape_t_maxpost_lasty=s$s_ns_tmp[[id]]$theta_maxpost[1]+ 
        #   s$s_ns_tmp[[id]]$theta_maxpost[2]*tail(Snn[[1]]$years,1)   # last year
        # s_tmp$scale_t_maxpost_lasty=s$s_ns_tmp[[id]]$theta_maxpost[3]+ 
        #   s$s_ns_tmp[[id]]$theta_maxpost[4]*tail(Snn[[1]]$years,1)   # last year
        # s_tmp$n_decades=tail(S_tmp$years,1) / 10 # n. of decades in the studied period
        # s_tmp$n_0y[id]=Snn[[jj]]$lmfit[[dur]]$coefficients[[1]]      # firs year
        # s_tmp$n_10y[id]=Snn[[jj]]$lmfit[[dur]]$coefficients[[1]]+
        #   Snn[[jj]]$lmfit[[dur]]$coefficients[[2]]*10                # decade
        # s_tmp$n_lasty[id]=Snn[[jj]]$lmfit[[dur]]$coefficients[[1]]+ 
        #   Snn[[jj]]$lmfit[[dur]]$coefficients[[2]]* tail(Snn[[1]]$years,1) # last year
        
        s_tmp$qCurve_maxpost_halfper=sapply(scale=s_tmp$scale_t_maxpost, 
                                            shape=s_tmp$shape_t_maxpost, 
                                            n_ordev=s$n_ordev[[id]],  # "n" at halfperiod 
                                            FUN=weibull_inv, 
                                            X=s_tmp$pr)
        # s_tmp$qCurve_maxpost_0y = sapply(scale  = s_tmp$scale_t_maxpost, 
        #                                   shape  = s_tmp$shape_t_maxpost, 
        #                                   n_ordev= s$n_0y[[id]], 
        #                                   FUN    = weibull_inv, 
        #                                   X      = s_tmp$pr)
        # s_tmp$qCurve_maxpost_10y = sapply(scale  = s_tmp$scale_t_maxpost, 
        #                                   shape  = s_tmp$shape_t_maxpost, 
        #                                   n_ordev= s$n_10y[[id]], 
        #                                   FUN    = weibull_inv, 
        #                                   X      = s_tmp$pr)
        # s_tmp$qCurve_maxpost_lasty = sapply(scale  = s_tmp$scale_t_maxpost, 
        #                                     shape  = s_tmp$shape_t_maxpost, 
        #                                     n_ordev= s$n_lasty[[id]], 
        #                                     FUN    = weibull_inv, 
        #                                     X      = s_tmp$pr)
        # 
        
        # Median curve + credibility interval at 90% and at 95%:
        med_ns_post=apply(s_tmp$qCurve_Post,2,median)
        low_ns_post=apply(s_tmp$qCurve_Post,2,quantile,0.05,na.rm=T)
        high_ns_post=apply(s_tmp$qCurve_Post,2,quantile,0.95,na.rm=T)
        low2.5_ns_post=apply(s_tmp$qCurve_Post,2,quantile,0.025,na.rm=T)
        high97.5_ns_post=apply(s_tmp$qCurve_Post,2,quantile,0.975,na.rm=T)
        s_tmp$spaghetti_qnt_df_ns_post <- data.frame(
          x=s_tmp$Tgrid,
          low_ns_post=low_ns_post,
          high_ns_post=high_ns_post,
          low2.5_ns_post=low2.5_ns_post,
          high97.5_ns_post=high97.5_ns_post,
          med_ns_post=med_ns_post,
          MAP_ns_post=s_tmp$qCurve_maxpost)
        # spaghetti (filtered for plots, e.g., filter_spaghetti=10):
        s_tmp$spaghetti_df_ns_post=data.frame(
          grid=s_tmp$Tgrid,
          spaghetti=list(t(s_tmp$qCurve_Post[seq(1,Nsim,filter_spaghetti),])))
        s_tmp$spaghetti_df_melt_ns_post<-melt(s_tmp$spaghetti_df_ns_post, 
                                              id.vars='grid', 
                                              variable.name='series') 
      }
      #^************************************************************************
      #                    NONSTATIONARY (MAXLIKE):
      #^************************************************************************
      if (!is.null(s$s_ns_tmp[[id]]$theta_maxlik_mcmc)){
        s_tmp$shape_t_maxlik=s$s_ns_tmp[[id]]$theta_maxlik_mcmc[1]+
          s$s_ns_tmp[[id]]$theta_maxlik_mcmc[2]*s$halfperiod
        s_tmp$scale_t_maxlik=s$s_ns_tmp[[id]]$theta_maxlik_mcmc[3]+
          s$s_ns_tmp[[id]]$theta_maxlik_mcmc[4]*s$halfperiod
        s_tmp$qCurve_maxlik=sapply(scale=s_tmp$scale_t_maxlik,  
                                   shape=s_tmp$shape_t_maxlik,  
                                   n_ordev=s$n_ordev[[id]], 
                                   FUN=weibull_inv, 
                                   X=s_tmp$pr)
        
        # Median curve + credibility interval at 90% and at 95%:
        med_ns_lik=apply(s_tmp$qCurve_Lik,2,median)
        low_ns_lik=apply(s_tmp$qCurve_Lik,2,quantile,0.05,na.rm=T)
        high_ns_lik=apply(s_tmp$qCurve_Lik,2,quantile,0.95,na.rm=T)
        low2.5_ns_lik=apply(s_tmp$qCurve_Lik,2,quantile,0.025,na.rm=T)
        high97.5_ns_lik=apply(s_tmp$qCurve_Lik,2,quantile,0.975,na.rm=T)
        s_tmp$spaghetti_qnt_df_ns_lik<-data.frame(x=s_tmp$Tgrid,
                                                  low_ns_lik=low_ns_lik,
                                                  high_ns_lik=high_ns_lik,
                                                  low2.5_ns_lik=low2.5_ns_lik,
                                                  high97.5_ns_lik=high97.5_ns_lik,
                                                  med_ns_lik=med_ns_lik,
                                                  MAP_ns_lik=s_tmp$qCurve_maxlik)
        # spaghetti (filtered for plots, e.g., filter_spaghetti=10):
        s_tmp$spaghetti_df_ns_lik=data.frame(grid=s_tmp$Tgrid,
                                             spaghetti=list(t(s_tmp$qCurve_Lik[
                                               seq(1,Nsim,filter_spaghetti),])))
        s_tmp$spaghetti_df_melt_ns_lik<-melt(s_tmp$spaghetti_df_ns_lik, 
                                             id.vars='grid', 
                                             variable.name='series') 
      }
      #^************************************************************************  
      #             PARTIAL-NON-STATIONARY (MAXPOST):
      #^************************************************************************
      if (!is.null(s$s_partns_tmp[[id]]$theta_maxpost)){
        s_tmp$shape_t_maxpost_partns= s$s_partns_tmp[[id]]$theta_maxpost[1]
        s_tmp$scale_t_maxpost_partns= s$s_partns_tmp[[id]]$theta_maxpost[2] +
          s$s_partns_tmp[[id]]$theta_maxpost[3]*s$halfperiod
        s_tmp$qCurve_maxpost_partns=sapply(scale=s_tmp$scale_t_maxpost_partns, 
                                           shape=s_tmp$shape_t_maxpost_partns, 
                                           n_ordev=s$n_ordev[[id]], 
                                           FUN=weibull_inv, 
                                           X=s_tmp$pr)
        # Median curve + credibility interval at 90% and at 95%:
        med_partns_post=apply(s_tmp$qCurve_Post_partns,2,median)
        low_partns_post=apply(s_tmp$qCurve_Post_partns,2,quantile,0.05,na.rm=T)
        high_partns_post=apply(s_tmp$qCurve_Post_partns,2,quantile,0.95,na.rm=T)
        low2.5_partns_post=apply(s_tmp$qCurve_Post_partns,2,quantile,0.025,na.rm=T)
        high97.5_partns_post=apply(s_tmp$qCurve_Post_partns,2,quantile,0.975,na.rm=T)
        s_tmp$spaghetti_qnt_df_partns_post <- data.frame(
          x= s_tmp$Tgrid,
          low_partns_post=low_partns_post,
          high_partns_post=high_partns_post,
          low2.5_partns_post=low2.5_partns_post,
          high97.5_partns_post=high97.5_partns_post,
          med_partns_post=med_partns_post,
          MAP_partns_post=s_tmp$qCurve_maxpost_partns)
        # spaghetti (filtered for plots, e.g.,filter_spaghetti=10):
        s_tmp$spaghetti_df_partns_post=data.frame(
          grid=s_tmp$Tgrid,
          spaghetti=list(t(s_tmp$qCurve_Post_partns[seq(1,Nsim,filter_spaghetti),])))
        s_tmp$spaghetti_df_melt_partns_post <- melt(s_tmp$spaghetti_df_partns_post , 
                                               id.vars='grid',
                                               variable.name='series') 
      }
      #^************************************************************************    
      #             PARTIAL-NON-STATIONARY (MAXLIKE):
      #^************************************************************************
      if (!is.null(s$s_partns_tmp[[id]]$theta_maxlik_mcmc)){
        s_tmp$shape_t_maxlik_partns=s$s_partns_tmp[[id]]$theta_maxlik_mcmc[1]
        s_tmp$scale_t_maxlik_partns=s$s_partns_tmp[[id]]$theta_maxlik_mcmc[2]+
          s$s_partns_tmp[[id]]$theta_maxlik_mcmc[3]*s$halfperiod
        s_tmp$qCurve_maxlik_partns=sapply(scale=s_tmp$scale_t_maxlik_partns,  
                                            shape=s_tmp$shape_t_maxlik_partns,  
                                            n_ordev=s$n_ordev[[id]], 
                                            FUN=weibull_inv,
                                            X=s_tmp$pr)
        # Median curve + credibility interval at 90% and at 95%:
        med_partns_lik=apply(s_tmp$qCurve_Lik_partns,2,median)
        low_partns_lik=apply(s_tmp$qCurve_Lik_partns,2,quantile,0.05,na.rm=T)
        high_partns_lik=apply(s_tmp$qCurve_Lik_partns,2,quantile,0.95,na.rm=T)
        low2.5_partns_lik=apply(s_tmp$qCurve_Lik_partns,2,quantile,0.025,na.rm=T)
        high97.5_partns_lik=apply(s_tmp$qCurve_Lik_partns,2,quantile,0.975,na.rm=T)
        s_tmp$spaghetti_qnt_df_partns_lik<-data.frame(x=s_tmp$Tgrid,
                                                      low_partns_lik=low_partns_lik,
                                                      high_partns_lik=high_partns_lik,
                                                      low2.5_partns_lik=low2.5_partns_lik,
                                                      high97.5_partns_lik=high97.5_partns_lik,
                                                      med_partns_lik=med_partns_lik,
                                                      MAP_partns_lik=s_tmp$qCurve_maxlik_partns)
        # spaghetti (filtered for plots, e.g., filter_spaghetti=10):
        s_tmp$spaghetti_df_partns_lik=data.frame(grid=s_tmp$Tgrid,
                                                 spaghetti=list(t(s_tmp$qCurve_Lik_partns[
                                                   seq(1,Nsim,filter_spaghetti),])))
        s_tmp$spaghetti_df_melt_partns_lik<-melt(s_tmp$spaghetti_df_partns_lik, 
                                                 id.vars='grid', 
                                                 variable.name='series')
      }
      #^************************************************************************    
      #                       STATIONARY SMEV (MAXPOST):
      #^************************************************************************
      if (!is.null(s$s_stat_tmp[[id]]$theta_maxpost)){
        s_tmp$shape_t_maxpost_stat=s$s_stat_tmp[[id]]$theta_maxpost[1]
        s_tmp$scale_t_maxpost_stat=s$s_stat_tmp[[id]]$theta_maxpost[2]
        s_tmp$qCurve_maxpost_stat=sapply(scale=s_tmp$scale_t_maxpost_stat, 
                                         shape=s_tmp$shape_t_maxpost_stat, 
                                         n_ordev=s$n_ordev[[id]],
                                         FUN=weibull_inv, 
                                         X=s_tmp$pr)
        # Median curve + credibility interval at 90% and at 95%:
        med_stat_post=apply(s_tmp$qCurve_Post_stat,2,median,na.rm=T)
        low_stat_post=apply(s_tmp$qCurve_Post_stat,2,quantile,0.05,na.rm=T)
        high_stat_post=apply(s_tmp$qCurve_Post_stat,2,quantile,0.95,na.rm=T)
        low2.5_stat_post=apply(s_tmp$qCurve_Post_stat,2,quantile,0.025,na.rm=T)
        high97.5_stat_post=apply(s_tmp$qCurve_Post_stat,2,quantile,0.975,na.rm=T)
        s_tmp$spaghetti_qnt_df_stat_post<-data.frame(x=s_tmp$Tgrid,
                                                     low_stat_post=low_stat_post,
                                                     high_stat_post=high_stat_post,
                                                     low2.5_stat_post=low2.5_stat_post,
                                                     high97.5_stat_post=high97.5_stat_post,
                                                     med_stat_post=med_stat_post,
                                                     MAP_stat_post=s_tmp$qCurve_maxpost_stat)
        # spaghetti (filtered for plots, e.g., filter_spaghetti=10):
        s_tmp$spaghetti_df_stat_post = data.frame(
          grid=s_tmp$Tgrid,
          spaghetti=list(t(s_tmp$qCurve_Post_stat[seq(1,Nsim,filter_spaghetti),])))
        s_tmp$spaghetti_df_melt_stat_post<-melt(s_tmp$spaghetti_df_stat_post, 
                                                id.vars='grid', 
                                                variable.name='series')
      }
      #^************************************************************************   
      #                       STATIONARY SMEV (MAXLIKE):
      #^************************************************************************
      if (!is.null(s$s_stat_tmp[[id]]$theta_maxlik_mcmc)){
        s_tmp$shape_t_maxlik_stat        = s$s_stat_tmp[[id]]$theta_maxlik_mcmc[1]
        s_tmp$scale_t_maxlik_stat        = s$s_stat_tmp[[id]]$theta_maxlik_mcmc[2]
        s_tmp$qCurve_maxlik_stat = sapply(scale  = s_tmp$scale_t_maxlik_stat,  
                                          shape  = s_tmp$shape_t_maxlik_stat,  
                                          n_ordev= s$n_ordev[[id]], 
                                          FUN    = weibull_inv, 
                                          X      = s_tmp$pr)
        # Median curve + credibility interval at 90% and at 95%:
        med_stat_lik     = apply(s_tmp$qCurve_Lik_stat, 2, median)
        low_stat_lik     = apply(s_tmp$qCurve_Lik_stat, 2, quantile,0.05, na.rm=T)
        high_stat_lik    = apply(s_tmp$qCurve_Lik_stat, 2, quantile,0.95, na.rm=T)
        low2.5_stat_lik  = apply(s_tmp$qCurve_Lik_stat, 2, quantile,0.025, na.rm=T)
        high97.5_stat_lik= apply(s_tmp$qCurve_Lik_stat, 2, quantile,0.975, na.rm=T)
        s_tmp$spaghetti_qnt_df_stat_lik <- data.frame(x                = s_tmp$Tgrid,
                                                      low_stat_lik     = low_stat_lik,
                                                      high_stat_lik    = high_stat_lik,
                                                      low2.5_stat_lik  = low2.5_stat_lik,
                                                      high97.5_stat_lik= high97.5_stat_lik,
                                                      med_stat_lik     = med_stat_lik,
                                                      MAP_stat_lik     = s_tmp$qCurve_maxlik_stat)
        # spaghetti (filtered for plots, e.g., filter_spaghetti=10):
        s_tmp$spaghetti_df_stat_lik=data.frame(
          grid=s_tmp$Tgrid,
          spaghetti=list(t(s_tmp$qCurve_Lik_stat[seq(1,Nsim,filter_spaghetti),])))
        s_tmp$spaghetti_df_melt_stat_lik <- melt(s_tmp$spaghetti_df_stat_lik , 
                                                 id.vars='grid', 
                                                 variable.name='series')
      }
      #^************************************************************************   
      #                   NON STATIONARY SMEV (MAXLIKE fminsearch):
      #^************************************************************************
      if (!is.null(s$s_ns_tmp[[id]]$theta_maxlik_fmins)){
        # shape
        s_tmp$shape_t_maxlik_fmins_ns=s$s_ns_tmp[[id]]$theta_maxlik_fmins[1]+ 
          s$s_ns_tmp[[id]]$theta_maxlik_fmins[2]*s$halfperiod
        # scale
        s_tmp$scale_t_maxlik_fmins_ns=s$s_ns_tmp[[id]]$theta_maxlik_fmins[3]+ 
          s$s_ns_tmp[[id]]$theta_maxlik_fmins[4]*s$halfperiod
        # quantile
        s_tmp$qCurve_maxlik_fmins_ns=sapply(scale  = s_tmp$scale_t_maxlik_fmins_ns,
                                            shape  = s_tmp$shape_t_maxlik_fmins_ns, 
                                            n_ordev= s$n_ordev[[id]], 
                                            FUN    = weibull_inv, 
                                            X      = s_tmp$pr)
      }
      #^************************************************************************   
      #                PARTIAL NON STATIONARY SMEV (MAXLIKE fminsearch):
      #^************************************************************************
      if (!is.null(s$s_partns_tmp[[id]]$theta_maxlik_fmins)){
        s_tmp$shape_t_maxlik_fmins_partns= s$s_partns_tmp[[id]]$theta_maxlik_fmins[1]
        s_tmp$scale_t_maxlik_fmins_partns= s$s_partns_tmp[[id]]$theta_maxlik_fmins[2]+
          s$s_partns_tmp[[id]]$theta_maxlik_fmins[3]*s$halfperiod
        s_tmp$qCurve_maxlik_fmins_partns = sapply(scale  = s_tmp$scale_t_maxlik_fmins_partns, 
                                                  shape  = s_tmp$shape_t_maxlik_fmins_partns, 
                                                  n_ordev= s$n_ordev[[id]], 
                                                  FUN    = weibull_inv, 
                                                  X      = s_tmp$pr)
      }
      #^************************************************************************   
      #                   STATIONARY SMEV (MAXLIKE fminsearch):
      #^************************************************************************
      if (!is.null(s$s_stat_tmp[[id]]$theta_maxlik_fmins)){
        s_tmp$shape_t_maxlik_fmins_stat  = s$s_stat_tmp[[id]]$theta_maxlik_fmins[1]
        s_tmp$scale_t_maxlik_fmins_stat  = s$s_stat_tmp[[id]]$theta_maxlik_fmins[2]
        s_tmp$qCurve_maxlik_fmins_stat = sapply(scale  = s_tmp$scale_t_maxlik_fmins_stat, 
                                                shape  = s_tmp$shape_t_maxlik_fmins_stat, 
                                                n_ordev= s$n_ordev[[id]], 
                                                FUN    = weibull_inv, 
                                                X      = s_tmp$pr)
      }
      #^************************************************************************    
      #                      STATIONARY GEV (MAXPOST):
      #^************************************************************************
      if (!is.null(s$s_gev_tmp[[id]]$theta_maxpost)){
        s_tmp$loc_t_maxpostGEV           = s$s_gev_tmp[[id]]$theta_maxpost[1]
        s_tmp$shape_t_maxpostGEV         = s$s_gev_tmp[[id]]$theta_maxpost[3]
        s_tmp$scale_t_maxpostGEV         = s$s_gev_tmp[[id]]$theta_maxpost[2] 
        s_tmp$qCurve_maxpostGEV = sapply(loc   = s_tmp$loc_t_maxpostGEV, 
                                         scale = s_tmp$scale_t_maxpostGEV, 
                                         shape = s_tmp$shape_t_maxpostGEV, 
                                         FUN   = qgev, 
                                         X=s_tmp$pr)
        # Median curve + credibility interval at 90% and at 95%:
        med_gev_post=apply(s_tmp$qCurve_PostGEV,2,median)
        low_gev_post=apply(s_tmp$qCurve_PostGEV,2,quantile,0.05,na.rm=T)
        high_gev_post=apply(s_tmp$qCurve_PostGEV,2,quantile,0.95,na.rm=T)
        low2.5_gev_post=apply(s_tmp$qCurve_PostGEV,2,quantile,0.025,na.rm=T)
        high97.5_gev_post=apply(s_tmp$qCurve_PostGEV,2,quantile,0.975,na.rm=T)
        s_tmp$spaghetti_qnt_df_gev_post <- data.frame(
          x=s_tmp$Tgrid,
          low_gev_post=low_gev_post,
          high_gev_post=high_gev_post,
          low2.5_gev_post=low2.5_gev_post,
          high97.5_gev_post=high97.5_gev_post,
          med_gev_post=med_gev_post,
          MAP_gev_post=s_tmp$qCurve_maxpostGEV)
        # spaghetti (filtered for plots, e.g., filter_spaghetti=10):
        s_tmp$spaghetti_df_gev_post = data.frame(
          grid= s_tmp$Tgrid,
          spaghetti=list(t(s_tmp$qCurve_PostGEV[seq(1,Nsim,filter_spaghetti),])))
        s_tmp$spaghetti_df_melt_gev_post <- melt(s_tmp$spaghetti_df_gev_post, 
                                                 id.vars='grid', 
                                                 variable.name='series')
      }
      
    #^**************************************************************************
    } else if (estim_algorithm=="fminsearch"){
    #^**************************************************************************
      
      #^************************************************************************   
      #                   NON STATIONARY SMEV (MAXLIKE fminsearch):
      #^************************************************************************
      if (!is.null(s$s_ns_tmp[[id]]$theta_maxlik_fmins)){
        # shape
        s_tmp$shape_t_maxlik_fmins_ns=s$s_ns_tmp[[id]]$theta_maxlik_fmins[1]+ 
          s$s_ns_tmp[[id]]$theta_maxlik_fmins[2]*s$halfperiod
        # scale
        s_tmp$scale_t_maxlik_fmins_ns=s$s_ns_tmp[[id]]$theta_maxlik_fmins[3]+ 
          s$s_ns_tmp[[id]]$theta_maxlik_fmins[4]*s$halfperiod
        # quantile
        s_tmp$qCurve_maxlik_fmins_ns=sapply(scale  = s_tmp$scale_t_maxlik_fmins_ns,
                                            shape  = s_tmp$shape_t_maxlik_fmins_ns, 
                                            n_ordev= s$n_ordev[[id]], 
                                            FUN    = weibull_inv, 
                                            X      = s_tmp$pr)
      }
      #^************************************************************************   
      #                PARTIAL NON STATIONARY SMEV (MAXLIKE fminsearch):
      #^************************************************************************
      if (!is.null(s$s_partns_tmp[[id]]$theta_maxlik_fmins)){
        s_tmp$shape_t_maxlik_fmins_partns= s$s_partns_tmp[[id]]$theta_maxlik_fmins[1]
        s_tmp$scale_t_maxlik_fmins_partns= s$s_partns_tmp[[id]]$theta_maxlik_fmins[2]+
          s$s_partns_tmp[[id]]$theta_maxlik_fmins[3]*s$halfperiod
        s_tmp$qCurve_maxlik_fmins_partns = sapply(scale  = s_tmp$scale_t_maxlik_fmins_partns, 
                                                  shape  = s_tmp$shape_t_maxlik_fmins_partns, 
                                                  n_ordev= s$n_ordev[[id]], 
                                                  FUN    = weibull_inv, 
                                                  X      = s_tmp$pr)
      }
      #^************************************************************************   
      #                   STATIONARY SMEV (MAXLIKE fminsearch):
      #^************************************************************************
      if (!is.null(s$s_stat_tmp[[id]]$theta_maxlik_fmins)){
        s_tmp$shape_t_maxlik_fmins_stat  = s$s_stat_tmp[[id]]$theta_maxlik_fmins[1]
        s_tmp$scale_t_maxlik_fmins_stat  = s$s_stat_tmp[[id]]$theta_maxlik_fmins[2]
        s_tmp$qCurve_maxlik_fmins_stat = sapply(scale  = s_tmp$scale_t_maxlik_fmins_stat, 
                                                shape  = s_tmp$shape_t_maxlik_fmins_stat, 
                                                n_ordev= s$n_ordev[[id]], 
                                                FUN    = weibull_inv, 
                                                X      = s_tmp$pr)
      }
    }       
    
    # # save all statistics in one dataframe:
    # ######################################
    # s_tmp$spaghetti_qnt_df <- data.frame(
    #   x                 = s_tmp$Tgrid, 
    #   # ns:
    #   low               = low,
    #   high              = high,
    #   low2.5            = low2.5,
    #   high97.5          = high97.5,
    #   med               = med,
    #   MAP               = s_tmp$qCurve_maxpost,
    #   maxlike           = s_tmp$qCurve_maxlik,
    #   maxlike_fmins_ns  = s_tmp$qCurve_maxlik_fmins_ns,
    #   
    #   # part ns:
    #   low_partns        = low_partns,
    #   high_partns       = high_partns,
    #   low2.5_partns     = low2.5_partns,
    #   high97.5_partns   = high97.5_partns,
    #   med_partns        = med_partns,
    #   MAP_partns        = s_tmp$qCurve_maxpost_partns,
    #   maxlike_partns    = s_tmp$qCurve_maxlik_partns,
    #   maxlike_fmins_partns=s_tmp$qCurve_maxlik_fmins_partns,
    #   
    #   # stat:
    #   low_stat          = low_stat,
    #   high_stat         = high_stat,
    #   low2.5_stat       = low2.5_stat,
    #   high97.5_stat     = high97.5_stat,
    #   med_stat          = med_stat,
    #   MAP_stat          = s_tmp$qCurve_maxpost_stat,
    #   maxlike_stat      = s_tmp$qCurve_maxlik_stat,
    #   maxlike_fmins_stat= s_tmp$qCurve_maxlik_fmins_stat,
    #   
    #   # GEV stat:
    #   MAPgev            = s_tmp$qCurve_maxpostGEV)

    
  
    
  ##############################################################################
  } else {
  ##############################################################################
    
    # if flag_save is FALSE then compute quantiles only for a few return 
    # periods: 2, 5, 10, 50, 100 years.
    s_tmp$Tgrid=c(2,5,10,50,100) # just a few (HARDCODED)! and no MaxLikelihood
    s_tmp$pr=1-1/s_tmp$Tgrid # prob associated to the target_return_periods
    
    #^**************************************************************************
    if (estim_algorithm=="cmdstan" | estim_algorithm=="stanr" | 
        estim_algorithm=="adaptMCMC"){
    #^**************************************************************************
      # initialization of arrays:
      if (!is.null(s$s_ns_tmp[[id]]$zPost)){
        s_tmp$qCurve_Post=matrix(NA,NROW(s$s_ns_tmp[[id]]$zPost),
                                        length(s_tmp$Tgrid))
        s_tmp$qCurve_Post_partns=matrix(NA,NROW(s$s_partns_tmp[[id]]$zPost),
                                        length(s_tmp$Tgrid))
        s_tmp$qCurve_Post_stat=matrix(NA,NROW(s$s_stat_tmp[[id]]$zPost),
                                      length(s_tmp$Tgrid))
        s_tmp$qCurve_PostGEV=matrix(NA,NROW(s$s_gev_tmp[[id]]$zPost),
                                    length(s_tmp$Tgrid))
      } else {
        message("posterior sample from ns smev is missing! please, check!!! ")
      }
      
      # for all spaghetti:
      for(i in 1:NROW(s$s_ns_tmp[[id]]$zPost)){
        if (!is.null(s$s_ns_tmp[[id]]$zPost)){   # FULL NONSTAT SMEV
          shape_t_Post=s$s_ns_tmp[[id]]$zPost[i,1]+
            s$s_ns_tmp[[id]]$zPost[i,2]*s$halfperiod
          scale_t_Post=s$s_ns_tmp[[id]]$zPost[i,3]+
            s$s_ns_tmp[[id]]$zPost[i,4]*s$halfperiod
          s_tmp$qCurve_Post[i,]=sapply(scale=scale_t_Post, 
                                       shape=shape_t_Post,
                                       n_ordev=s$n_ordev[[id]],
                                       FUN=weibull_inv, 
                                       X=s_tmp$pr)
        }
        if (!is.null(s$s_partns_tmp[[id]]$zPost)){   # PART NONSTAT SMEV
          shape_t_Post_partns=s$s_partns_tmp[[id]]$zPost[i,1]
          scale_t_Post_partns=s$s_partns_tmp[[id]]$zPost[i,2]+
            s$s_partns_tmp[[id]]$zPost[i,3]*s$halfperiod
          s_tmp$qCurve_Post_partns[i,]=sapply(scale=scale_t_Post_partns,
                                              shape=shape_t_Post_partns,
                                              n_ordev=s$n_ordev[[id]], 
                                              FUN=weibull_inv, 
                                              X=s_tmp$pr)
        }
        if (!is.null(s$s_stat_tmp[[id]]$zPost)){   # STAT SMEV
          shape_t_Post_stat   = s$s_stat_tmp[[id]]$zPost[i,1]
          scale_t_Post_stat   = s$s_stat_tmp[[id]]$zPost[i,2]
          s_tmp$qCurve_Post_stat[i,]  = sapply(scale=scale_t_Post_stat,
                                               shape=shape_t_Post_stat,
                                               n_ordev=s$n_ordev[[id]], 
                                               FUN=weibull_inv, 
                                               X=s_tmp$pr)
        }
        if (!is.null(s$s_gev_tmp[[id]]$zPost)){    # GEV
          loc_t_PostGEV       = s$s_gev_tmp[[id]]$zPost[i,1]
          shape_t_PostGEV     = s$s_gev_tmp[[id]]$zPost[i,3]
          scale_t_PostGEV     = s$s_gev_tmp[[id]]$zPost[i,2]
          s_tmp$qCurve_PostGEV[i,]    = sapply(loc=loc_t_PostGEV, 
                                               scale=scale_t_PostGEV, 
                                               shape=shape_t_PostGEV, 
                                               FUN=qgev,
                                               X=s_tmp$pr)
        }
      }
      # MAP and CI:
      #############
      if (!is.null(s$s_ns_tmp[[id]]$theta_maxpost)){
        s_tmp$shape_t_maxpost=s$s_ns_tmp[[id]]$theta_maxpost[1]+
          s$s_ns_tmp[[id]]$theta_maxpost[2]*s$halfperiod
        s_tmp$scale_t_maxpost=s$s_ns_tmp[[id]]$theta_maxpost[3]+
          s$s_ns_tmp[[id]]$theta_maxpost[4]*s$halfperiod
        s_tmp$qCurve_maxpost=sapply(scale=s_tmp$scale_t_maxpost, 
                                    shape=s_tmp$shape_t_maxpost, 
                                    n_ordev=s$n_ordev[[id]], 
                                    FUN=weibull_inv, 
                                    X=s_tmp$pr)
        # Median curve + credibility interval at 90% and at 95%:
        med_ns_post     =apply(s_tmp$qCurve_Post,2,median)
        low_ns_post     =apply(s_tmp$qCurve_Post,2,quantile,0.05,na.rm=T)
        high_ns_post    =apply(s_tmp$qCurve_Post,2,quantile,0.95,na.rm=T)
        low2.5_ns_post  =apply(s_tmp$qCurve_Post,2,quantile,0.025,na.rm=T)
        high97.5_ns_post=apply(s_tmp$qCurve_Post,2,quantile,0.975,na.rm=T)
        s_tmp$spaghetti_qnt_df_ns_post<-data.frame(x=s_tmp$Tgrid,
                                                   low_ns_post=low_ns_post,
                                                   high_ns_post=high_ns_post,
                                                   low2.5_ns_post=low2.5_ns_post,
                                                   high97.5_ns_post=high97.5_ns_post,
                                                   med_ns_post=med_ns_post,
                                                   MAP_ns_post=s_tmp$qCurve_maxpost)
        # spaghetti (filtered for plots, e.g., filter_spaghetti=10):
        s_tmp$spaghetti_df_ns_post=data.frame(
          grid=s_tmp$Tgrid,
          spaghetti=list(t(s_tmp$qCurve_Post[seq(1,NROW(s_tmp$qCurve_Post), 
                                                 filter_spaghetti),])))
        s_tmp$spaghetti_df_melt_ns_post<-melt(s_tmp$spaghetti_df_ns_post, 
                                              id.vars='grid', 
                                              variable.name ='series') 
      }
      
      if (!is.null(s$s_partns_tmp[[id]]$theta_maxpost)){
        s_tmp$shape_t_maxpost_partns=s$s_partns_tmp[[id]]$theta_maxpost[1]
        s_tmp$scale_t_maxpost_partns=s$s_partns_tmp[[id]]$theta_maxpost[2] +
          s$s_partns_tmp[[id]]$theta_maxpost[3]*s$halfperiod
        s_tmp$qCurve_maxpost_partns = sapply(scale  = s_tmp$scale_t_maxpost_partns, 
                                             shape  = s_tmp$shape_t_maxpost_partns, 
                                             n_ordev= s$n_ordev[[id]], 
                                             FUN    = weibull_inv, 
                                             X      = s_tmp$pr)
        # Median curve + credibility interval at 90% and at 95%:
        med_partns_post=apply(s_tmp$qCurve_Post_partns,2,median)
        low_partns_post=apply(s_tmp$qCurve_Post_partns,2,quantile,0.05,na.rm=T)
        high_partns_post=apply(s_tmp$qCurve_Post_partns,2,quantile,0.95,na.rm=T)
        low2.5_partns_post=apply(s_tmp$qCurve_Post_partns,2,quantile,0.025,na.rm=T)
        high97.5_partns_post=apply(s_tmp$qCurve_Post_partns,2,quantile,0.975,na.rm=T)
        s_tmp$spaghetti_qnt_df_partns_post <- data.frame(
          x                   = s_tmp$Tgrid,
          low_partns_post     = low_partns_post,
          high_partns_post    = high_partns_post,
          low2.5_partns_post  = low2.5_partns_post,
          high97.5_partns_post= high97.5_partns_post,
          med_partns_post     = med_partns_post,
          MAP_partns_post     = s_tmp$qCurve_maxpost_partns)
        
        # spaghetti (filtered for plots, e.g., filter_spaghetti=10):
        s_tmp$spaghetti_df_partns_post = data.frame(
          grid=s_tmp$Tgrid,
          spaghetti=list(t(s_tmp$qCurve_Post_partns[seq(1,NROW(s_tmp$qCurve_Post_partns), 
                                                        filter_spaghetti),])))
        s_tmp$spaghetti_df_melt_partns_post<-melt(s_tmp$spaghetti_df_partns_post,
                                                  id.vars='grid', 
                                                  variable.name='series') 
      }
      
      if (!is.null(s$s_stat_tmp[[id]]$theta_maxpost)){
        s_tmp$shape_t_maxpost_stat=s$s_stat_tmp[[id]]$theta_maxpost[1]
        s_tmp$scale_t_maxpost_stat=s$s_stat_tmp[[id]]$theta_maxpost[2]
        s_tmp$qCurve_maxpost_stat=sapply(scale=s_tmp$scale_t_maxpost_stat, 
                                         shape=s_tmp$shape_t_maxpost_stat, 
                                         n_ordev=s$n_ordev[[id]],
                                         FUN=weibull_inv, 
                                         X=s_tmp$pr)
        # Median curve + credibility interval at 90% and at 95%:
        med_stat_post     =apply(s_tmp$qCurve_Post_stat,2,median)
        low_stat_post     =apply(s_tmp$qCurve_Post_stat,2,quantile,0.05,na.rm=T)
        high_stat_post    =apply(s_tmp$qCurve_Post_stat,2,quantile,0.95,na.rm=T)
        low2.5_stat_post  =apply(s_tmp$qCurve_Post_stat,2,quantile,0.025,na.rm=T)
        high97.5_stat_post=apply(s_tmp$qCurve_Post_stat,2,quantile,0.975,na.rm=T)
        s_tmp$spaghetti_qnt_df_stat_post <- data.frame(x                 =s_tmp$Tgrid,
                                                       low_stat_post     =low_stat_post,
                                                       high_stat_post    =high_stat_post,
                                                       low2.5_stat_post  =low2.5_stat_post,
                                                       high97.5_stat_post=high97.5_stat_post,
                                                       med_stat_post     =med_stat_post,
                                                       MAP_stat_post     =s_tmp$qCurve_maxpost_stat)
        # spaghetti (filtered for plots, e.g., filter_spaghetti=10):
        s_tmp$spaghetti_df_stat_post=data.frame(
          grid=s_tmp$Tgrid,
          spaghetti=list(t(s_tmp$qCurve_Post_stat[seq(1,NROW(s_tmp$qCurve_Post_stat), 
                                                      filter_spaghetti),])))
        s_tmp$spaghetti_df_melt_stat_post<-melt(s_tmp$spaghetti_df_stat_post,
                                                id.vars='grid', 
                                                variable.name='series')
      }
      
      
      if (!is.null(s$s_gev_tmp[[id]]$theta_maxpost)){
        s_tmp$loc_t_maxpostGEV=s$s_gev_tmp[[id]]$theta_maxpost[1]
        s_tmp$shape_t_maxpostGEV=s$s_gev_tmp[[id]]$theta_maxpost[3]
        s_tmp$scale_t_maxpostGEV= s$s_gev_tmp[[id]]$theta_maxpost[2] 
        s_tmp$qCurve_maxpostGEV=sapply(loc=s_tmp$loc_t_maxpostGEV, 
                                       scale=s_tmp$scale_t_maxpostGEV, 
                                       shape=s_tmp$shape_t_maxpostGEV, 
                                       FUN=qgev, 
                                       X=s_tmp$pr)
        # Median curve + credibility interval at 90% and at 95%:
        med_gev_post     =apply(s_tmp$qCurve_PostGEV,2,median)
        low_gev_post     =apply(s_tmp$qCurve_PostGEV,2,quantile,0.05,na.rm=T)
        high_gev_post    =apply(s_tmp$qCurve_PostGEV,2,quantile,0.95,na.rm=T)
        low2.5_gev_post  =apply(s_tmp$qCurve_PostGEV,2,quantile,0.025,na.rm=T)
        high97.5_gev_post=apply(s_tmp$qCurve_PostGEV,2,quantile,0.975,na.rm=T)
        s_tmp$spaghetti_qnt_df_gev_post <- data.frame(x                =s_tmp$Tgrid,
                                                      low_gev_post     =low_gev_post,
                                                      high_gev_post    =high_gev_post,
                                                      low2.5_gev_post  =low2.5_gev_post,
                                                      high97.5_gev_post=high97.5_gev_post,
                                                      med_gev_post     =med_gev_post,
                                                      MAP_gev_post     =s_tmp$qCurve_maxpostGEV)
        # spaghetti (filtered for plots, e.g., filter_spaghetti=10):
        s_tmp$spaghetti_df_gev_post = data.frame(
          grid= s_tmp$Tgrid,
          spaghetti=list(t(s_tmp$qCurve_PostGEV[seq(1,NROW(s_tmp$qCurve_PostGEV), 
                                                    filter_spaghetti),])))
        s_tmp$spaghetti_df_melt_gev_post <- melt(s_tmp$spaghetti_df_gev_post , 
                                                 id.vars='grid', 
                                                 variable.name='series')
      }
    
    #^**************************************************************************  
    } else if (estim_algorithm=="fminsearch"){
    #^**************************************************************************
    
      if (!is.null(s$s_ns_tmp[[id]]$theta_maxlik_fmins)){
        s_tmp$shape_t_maxlik_fmins_ns=s$s_ns_tmp[[id]]$theta_maxlik_fmins[1]+ 
          s$s_ns_tmp[[id]]$theta_maxlik_fmins[2]*s$halfperiod
        s_tmp$scale_t_maxlik_fmins_ns=s$s_ns_tmp[[id]]$theta_maxlik_fmins[3]+
          s$s_ns_tmp[[id]]$theta_maxlik_fmins[4]*s$halfperiod
        s_tmp$qCurve_maxlik_fmins_ns=sapply(scale  = s_tmp$scale_t_maxlik_fmins_ns,
                                            shape  = s_tmp$shape_t_maxlik_fmins_ns, 
                                            n_ordev= s$n_ordev[[id]], 
                                            FUN    = weibull_inv, 
                                            X      = s_tmp$pr)
        s_tmp$spaghetti_qnt_df_fmins_ns <- data.frame(
          x=s_tmp$Tgrid,
          qCurve_maxlik_fmins_ns=s_tmp$qCurve_maxlik_fmins_ns)
      }
      
      if (!is.null(s$s_partns_tmp[[id]]$theta_maxlik_fmins)){
        s_tmp$shape_t_maxlik_fmins_partns=s$s_partns_tmp[[id]]$theta_maxlik_fmins[1]
        s_tmp$scale_t_maxlik_fmins_partns=s$s_partns_tmp[[id]]$theta_maxlik_fmins[2]+
          s$s_partns_tmp[[id]]$theta_maxlik_fmins[3]*s$halfperiod
        s_tmp$qCurve_maxlik_fmins_partns=sapply(scale=s_tmp$scale_t_maxlik_fmins_partns, 
                                                shape=s_tmp$shape_t_maxlik_fmins_partns, 
                                                n_ordev=s$n_ordev[[id]], 
                                                FUN=weibull_inv, 
                                                X=s_tmp$pr)
        s_tmp$spaghetti_qnt_df_fmins_partns <- data.frame(
          x= s_tmp$Tgrid,
          qCurve_maxlik_fmins_partns=s_tmp$qCurve_maxlik_fmins_partns)
      }
      
      if (!is.null(s$s_stat_tmp[[id]]$theta_maxlik_fmins)){
        s_tmp$shape_t_maxlik_fmins_stat=s$s_stat_tmp[[id]]$theta_maxlik_fmins[1]
        s_tmp$scale_t_maxlik_fmins_stat=s$s_stat_tmp[[id]]$theta_maxlik_fmins[2]
        s_tmp$qCurve_maxlik_fmins_stat=sapply(scale=s_tmp$scale_t_maxlik_fmins_stat, 
                                              shape=s_tmp$shape_t_maxlik_fmins_stat, 
                                              n_ordev=s$n_ordev[[id]], 
                                              FUN=weibull_inv, 
                                              X=s_tmp$pr)
        s_tmp$spaghetti_qnt_df_fmins_stat <- data.frame(
          x=s_tmp$Tgrid,
          qCurve_maxlik_fmins_stat=s_tmp$qCurve_maxlik_fmins_stat)
      }
    }
    
  }
  # return main object with all quantiles results
  return(s_tmp)
}















################################################################################
compute_quantiles_MEV<-function(flg_save=T, 
                                s,
                                id=1,
                                start_val,
                                filter_spaghetti=10){
################################################################################
  #^* GOAL: Compute quantiles SMEV model with MEV formulation for given set of 
  #^*       return periods
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [list of lists] s --> 
  #^*    2. [array of integers] target_T --> given return periods c(2,10,50,100)
  #^*    3. [logical] flg_save --> if T it saves all the input return periods, 
  #^*                              otherwise it computes quantiles only for a 
  #^*                              limited set (2,5,10,50,100 years).
  #^*    4. [integer] id --> duration index 
  #^*    5. [integer] filter_spaghetti --> filter spghetti to obtain a reduced
  #^*                                      ensemble for plotting reasons
  #^* OUT
  #^*    1. [list of lists] s_tmp --> list with all quantiles results
  #^****************************************************************************
  #^* REF.:
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  ##############################################################################
  # Compute quantiles curves with uncertainty:
  s_tmp=c()
  if (flg_save==T){
    # define a limited set of return periods:
    Tgrid_onlyfew = sort(c(exp(seq(log(1.1), log(50),0.1)),2,5,10,20,50,100))
    qnt = matrix(NA, length(seq(1,NROW(s$s_ns_tmp[[id]]$zPost), 
                                filter_spaghetti)),length(Tgrid_onlyfew))
    k=0
    for(i in seq(1,NROW(s$s_ns_tmp[[id]]$zPost), filter_spaghetti)){
      #print(i)
      k=k+1
      for (Tr in 1:length(Tgrid_onlyfew)){
        p = 1-1/Tgrid_onlyfew[Tr]
        # perform inversion looking for the zeros with the fzero R function.
        # starting point (use smev result as starting value)
        x0 = sapply(scale   = start_val[1], 
                    shape   = start_val[2], 
                    n_ordev = s$n_ordev[[id]], 
                    FUN     = weibull_inv, 
                    X       = p)
        
        if ((s$s_ns_tmp[[id]]$zPost[i,1]+s$s_ns_tmp[[id]]$zPost[i,2]*s$M >0) & 
            (s$s_ns_tmp[[id]]$zPost[i,3]+s$s_ns_tmp[[id]]$zPost[i,4]*s$M >0)){
          # apply fzero function to solve the inversion:  
          tmp = pracma::fzero(f      =mevcdf,
                              x      =c(x0/10,x0*10),
                              theta  =s$s_ns_tmp[[id]]$zPost[i,],
                              coeff  =s$lmfit[[id]]$coefficients,
                              maxiter=500,
                              M      =s$M,
                              p      =p) # s$ord_events[[id]])
          # print(tail(tmp$x,1))
          qnt[k,Tr]=tmp$x
        } else {
          qnt[k,Tr]=NA
        }
      }
      # print(qnt[k,])
      # print(qnt[i,])
    }
    s_tmp$spaghetti_mev_df=data.frame(grid=Tgrid_onlyfew, 
                                      spaghetti=t(qnt))
    # t(qCurve_Post[seq(1,NROW(qCurve_Post),100),])))
    s_tmp$qnt=qnt
    s_tmp$spaghetti_mev_df_melt<-melt(s_tmp$spaghetti_mev_df,  
                                      id.vars='grid',
                                      variable.name='series')
    
    # median curve + credibility interval at 90%
    med_mev=apply(qnt,2,median,na.rm=T)
    low_mev=apply(qnt,2,quantile,0.05,na.rm=T)
    high_mev=apply(qnt,2,quantile,0.95,na.rm=T)
    low2.5_mev=apply( qnt,2,quantile,0.025,na.rm=T)
    high97.5_mev=apply( qnt,2,quantile,0.975,na.rm=T)
    
    # dataframe:
    s_tmp$spaghetti_mev_qnt_df<-data.frame(x=Tgrid_onlyfew,
                                           med=med_mev,
                                           low=low_mev,
                                           high=high_mev,
                                           low2.5=low2.5_mev,
                                           high97.5 =high97.5_mev)
  } else {
    s_tmp=NULL
  }
  return(s_tmp)
}











