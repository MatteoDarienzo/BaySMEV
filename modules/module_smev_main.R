################################################################################
smev_analysis_main = function(rain.df,       # rainfall timeseries (dataframe)
                              name_project      = "TEST",
                              dir_results       = "",
                              date_initial      = NULL,
                              date_final        = NULL,
                              min_rain          = 0.1,
                              time_resolution   = 60,
                              min_ev_duration   = 30,
                              separation_in_min = 1440,
                              duration          = c(1,3,6,24)*60, # minutes!
                              target_T          = c(2,5,10,20,25,50,100),
                              thr_leftcens      = 0.9,
                              estim_algorithm   = "cmdstan", # "fminsearch"
                              priors            = list(),                               
                              priors_gev        = list(loc=list("uniform",0,100,1),
                                                       scale=list("uniform",0,100,1),
                                                       shape=list("normal",0,0.3,0.1)),                           
                              Nmcmc             = 4000,                        
                              Nmcmc_max         = 40000,                     
                              acc.rate          = 0.234,                       
                              Nslim             = 1,                         
                              nburn             = 0.5,                     
                              scales            = c(0.02,0.002,0.1,0.005),  
                              scales_gev        = c(0.1,0.1,0.1),
                              pix_row           = NULL, 
                              pix_col           = NULL,
                              sat_prod          = NULL,
                              flag_like         = F,                         
                              flag_saveALL      = T,                        
                              flag_gev_stat     = F,
                              flag_smev_stat    = T,
                              flag_smev_ns      = T,
                              flag_smev_partns  = F,
                              flag_model_compar = T,
                              flag_save2rds     = T
                              ){
################################################################################
  #^* GOAL: Simplified Metastatistical Extreme Value (SMEV) method for rainfall#
  #^*       extreme analysis.                                                  #
  #^*       Both non-stationary and stationary models are available.           #
  #^*       Both MLE (fminsearch) and Bayesian (STAN) methods are available.   #
  #^***************************************************************************#
  #^* PROGRAMMERS: Matteo Darienzo, University of Padua, Italy                 #
  #^*              (matteo.darienzo@unipd.it)                                  #
  #^*              Francesco Marra, University of Padua, Italy                 #
  #^*              (francesco.marra@unipd.it)                                  #
  #^***************************************************************************#
  #^* CREATED/MODIFIED: 2025/11/28                                             #
  #^***************************************************************************#
  #^* VERSION: v. 1.0                                                          #
  #^***************************************************************************#
  #^* FUNDS:                                                                   #
  #^* supported by the INTENSE project (raINfall exTremEs and their impacts:   #
  #^* from the local to the National ScalE) funded by European Union-NextGenEU #
  #^* within PRIN (Progetti di ricerca di Rilevante Interesse Nazionale)       #
  #^* programme (grant 2022ZC2522).                                            #
  #^***************************************************************************#
  #^*                            DESCRIPTION:                                  #
  #^***************************************************************************#
  #^* This program performs the SMEV method (Marra et al., 2020) for rainfall  #
  #^* extreme analysis, which is a simplified version of the MEV approach      #
  #^* (Marani and Ignaccolo,2015). This method makes use of all ordinary events#
  #^* not only of a small sample of extremes (e.g., annual maxima).            #
  #^* More specifically, it uses a 2-parameter Weibull distrib. with the       #
  #^* left-Censoring approach to model the tail of ordinary events.            #
  #^* Different parameter estimation methods are available:                    #
  #^*    -  Maximum Likelihood estimation method with the fminsearch()         #
  #^*       minimisation function, as in the original matlab version of the    #
  #^*       method (https://zenodo.org/records/11934843).                      #
  #^*    -  Maximum likelihood or Posterior sampling within a Bayesian         #
  #^*       framework using an adaptive Metropolis algorithm ("adaptMCMC"      #
  #^*       R package of A. Scheidegger doi:10.1007/s11222-011-9269-5)         #
  #^*    -  Maximum likelihood or Posterior sampling within a Bayesian         #
  #^*       framework using an Hamiltonian MCMC algorithm with packages "Rstan"#
  #^*       "cmdstanr" R packages.                                             #
  #^* Moreover, both stationary and non-stationary SMEV frameworks available.  #
  #^* For comparison, also a stationary GEV model is available.                # 
  #^***************************************************************************#
  #^*                          ACKNOWLEDGEMENTS:                               #
  #^***************************************************************************#
  #^* For the extreme analysis it uses methods: MEV (Marani and Ignaccolo,2015)#
  #^* and its simplified version SMEV (Marra et al., 2019. 2020).              #
  #^* It uses open source codes of Francesco Marra (University of Padua) for   #
  #^* storm separation, identification of ordinary events, SMEV method and     #
  #^* left-censoring approach, and non-stationary model freely available at:   #
  #^* https://zenodo.org/records/11934843 (SMEV methodology)                   #
  #^* https://zenodo.org/records/15047817 (non-stationary SMEV)                #
  #^*                                                                          #
  #^* For the Bayesian and MCMC algorithms it uses some ideas from:            #
  #^* - RsTooDS of Benjamin Renard (INRAE): https://zenodo.org/records/5075760 #
  #^* - BaM of Benjamin Renard (INRAE): https://github.com/BaM-tools/RBaM      #
  #^* - BayDERS of Matteo Darienzo(INRAE):                                     #
  #^*   https://github.com/MatteoDarienzo/BayDERS                              #
  #^* - packages "rstan" and "cmdstanr" for stan Bayesian inference + MCMC alg.#
  #^* - package "adaptMCMC" for adaptive metropolis algorithm.                 #
  #^* For the trend analysis it uses package 'trend' of Thorsten Pohlert:      #
  #^* (https://cran.r-project.org/web/packages/trend/index.html)               #
  #^***************************************************************************#
  #^*                              DISCLAIMER:                                 #
  #^***************************************************************************#
  #^* Please, notice that is an experimental code. Further analysis is required#
  #^* for validating the proposed tools and scripts. We are not responsible for#
  #^* any loss (e.g.,of data, profits, business, customers, orders, other costs#
  #^* or disturbances) derived by their use in the operational practice.       #
  #^* The authorized user accepts the risks associated with using the software #
  #^* given its nature of free software. It is reserved to expert Users        #
  #^* (developers or professionals) having knowledge in computer science,      #
  #^* hydrology, extreme analysis and statistics.                              #
  #^***************************************************************************#
  #^*                               READ ME:                                   #
  #^***************************************************************************#
  #^* This is a main function, which makes use of several other functions /    #
  #^* modules which must be priorly loaded.                                    #
  #^* Moreover, please make sure you have installed the necessary packages:    #
  #^*          "methods","lattice","gridExtra","reshape","reshape2","ggplot2", #
  #^*          "extrafont","grid","gtable","chron","RColorBrewer","cowplot",   #
  #^*          "psych","mosaicData","tidyr","lubridate","parallel","latex2exp",#
  #^*          "scales","dplyr","png","pracma","ncdf4","CFtime","rnaturalearth"#
  #^*          "raster","maptools","Kendall","trend","tidyverse","ggpubr",     #
  #^*          "sf","coda","mcmc","evd","stats","R.matlab","extRemes","imager",#
  #^*          "mmand","adaptMCMC","PEIP","posterior","loo","bayesplot","rstan"#
  #^***************************************************************************#
  #^*                               INPUTS                                     #
  #^***************************************************************************#
  #^*   1. [dataframe] rain.df --> input rainfall timeseries.                  #
  #^*                              One column 'timestamps' with the dates in   #
  #^*                              timestamp, and a column named 'y' with      #
  #^*                              rainfall values [mm]. Time step is defined  #
  #^*                              in 'time_resolution' setting.               #
  #^*   2. [character] name_project --> name of the folder for the results     #
  #^*   3. [character] dir_results --> path of the folder with results         #
  #^*   4. [timestamp] date_initial --> [optional] starting date of the period #
  #^*   5. [timestamp] date_final --> [optional] ending date of the period     #
  #^*   6. [real] min_rain --> minimum precipitation value to define wet values#
  #^*   7. [integer] time_resolution --> time step of the input prec timeseries#
  #^*                                    in minutes (e.g., 60 for hourly data).#
  #^*   8. [integer] min_ev_duration --> min duration [min] of an event, ex:30 #
  #^*   9. [integer] separation_in_min --> min time [min] for event separation,#
  #^*                                      ex:1440 for 1 day of dry period     # 
  #^*   10. [array of integers] duration --> durations (in minutes) to analyse #
  #^*                                        =c(1,3,6,24)*60                   #
  #^*   11. [array of integers] target_T --> return times for the quantiles:   #
  #^*                                        =c(2,5,10,20,25,50, 100)          #
  #^*   12. [real] thr_leftcens --> threshold percentile for the left censoring#
  #^*                               (e.g., 0.9 to censor 90% of ordin.events). #
  #^*   13. [character] estim_algorithm --> choose algorithm for estimation of #
  #^*                                       SMEV parameters:                   #
  #^*                                       "cmdstan", "rstan", "adaptMCMC",   #
  #^*                                       "fminsearch".                      #
  #^*                                       "fminsearch" (MLE) is by default.  #
  #^*   14. [list of lists] priors --> list of priors for SMEV parameters,     #
  #^*                                  independently if stat or nonstat model  # 
  #^*                                  (shape intercept, shape slope,          #
  #^*                                   scale intercpet, scale slope).         #
  #^*                                  one list for each selected duration,and #
  #^*                                  one list for each parameter, including  #
  #^*                                  distrib.name,mean,stdev,startvalue, ex: #
  #^*                                  priors=list(); priors[[1]]=list(        #
  #^*                                  shape_intercept=list("lognormal",       #
  #^*                                                       log(0.7),0.4,0.7), #
  #^*                                  shape_slope=list("normal",0,0.05,0),    #
  #^*                                  scale_intercept=list(...),              #    
  #^*                                  scale_slope=list(...))                  #
  #^*   15. [list of array] priors_gev --> priors for Bayesian estimation of   #
  #^*                                      stationary GEV (if flag_gev_stat=T) # 
  #^*                                      one list for each parameter, as:    #
  #^*                                      distrib.name,mean,stdev,init value. #
  #^*   16. [integer] Nmcmc --> tot number of MCMC simulations, e.g., 4000     # 
  #^*   17. [integer] Nmcmc_max --> (only if "adaptMCMC" algorithm is selected)#
  #^*                               max total number of MCMC sim, e.g., 30000  # 
  #^*   18. [real] acc.rate --> (only if "adaptMCMC" algorithm is selected)    #
  #^*                           targeted acceptance rate for mcmc iter. (0.234)#
  #^*   19. [integer] Nslim --> (only if "adaptMCMC" algorithm is selected)    #
  #^*                           slim mcmc for posterior statistics (to reduce  #
  #^*                           lag-autocorrelation)                           # 
  #^*   20. [real] nburn --> burn first fraction of the MCMC iterations to     #
  #^*                        reduce impact of start value on mcmc convergence. #                                    
  #^*   21. [array of reals] start_values --> (only if "adaptMCMC" algorithm   #
  #^*                                         is selected).                    #
  #^*                                         Values for each model parameters,#
  #^*                                         (shape intercept, shape slope,   #
  #^*                                         scale intercept, scale slope),ex:#     
  #^*                                         =c(0.7, 0, 5, 0)                 #
  #^*   22. [array of reals] scales -->(only if "adaptMCMC" algorithm is chosen#
  #^*                                    scaling setting of mcmc jumps to      #
  #^*                                    explore SMEV parameters distribution: #
  #^*                                    =c(0.02, 0.002, 0.1, 0.005)           #
  #^*   23. [array of reals] scales_gev --> (only if "adaptMCMC" algorithm is  #
  #^*                                    selected).                            #
  #^*                                    scaling setting of mcmc jumps to      #
  #^*                                    explore GEV parameters distribution:  #
  #^*                                    =c(0.1,0.1,0.1), for loc,scale,shape  #
  #^*   24. [integer] pix_row --> if gridded dataset, row index                #
  #^*   25. [integer] pix_col --> if gridded dataset, column index             #
  #^*   26. [character] sat_prod --> if any, name of the satellite product     # 
  #^*   27. [logical] flag_like --> flag to activate maximum likelihood method #
  #^*   28. [logical] flag_saveALL --> flag to compute all quantiles and       #
  #^*                                  comparison between the chosen models or #
  #^*                                  just a few quantiles (2,5,10,20,50,100) #
  #^*   29. [logical] flag_gev_stat --> flag to infer stationary GEV model     #
  #^*   30. [logical] flag_smev_stat --> flag to infer stationary SMEV model   #
  #^*   31. [logical] flag_smev_ns --> flag to infer nonstationary SMEV model  #
  #^*   32. [logical] flag_smev_partns --> flag to infer partial nonstationary #
  #^*                                      SMEV model                          #
  #^*   33. [logical] flag_model_compar --> flag to compute comparison (plot)  #
  #^*                                       among all selected models.         #   
  #^*   34. [logical] flag_save2rds --> flag to save all results into .rds file#
  #^*                                                                          #
  #^***************************************************************************#
  #^*                               OUTPUTS                                    #
  #^***************************************************************************#
  #^*   1. [list] s --> list object containing all results of the SMEV analysis#
  #^*                                                                          #
  #^***************************************************************************#
  #^*                               to do:                                     #
  #^***************************************************************************#
  #^* Implementation of classical MEV method                                   #
  #^* Implementation of MEV and SMEV methods including different storm types.  #
  #^* Implementation of nonstationarity with other covariate (e.g. temperature)#
  #^* Implementation of different non-stationary relation (e.g., exponential)  #
  #^* Implementation of other estimation methods (e.g., L-moments) and other   #
  #^* Bayesian MCMC algorithms.                                                #
  #^* Computation of Marginal Likelihood for computing Bayes Factor and other  #
  #^* Bayesian criteria.                                                       #
  ##############################################################################
  message("- Ordinary events separation ...")
  start_t= Sys.time()
  # initialize object "s" that will contain all results of the analysis
  rm(s)
  s<-list()
  
  #^***************************************************************************#
  #  CHECK AND ARRANGE INPUT PRECIPITATION TIMESERIES
  #^***************************************************************************#
  # reinitialize s by adding input rainfall time series to object s, as "vals"
  s<-list(vals=rain.df$y)
  # Fix negative (or below min threshold, "min_rain") values
  s$vals[is.na(s$vals)]=0
  s$vals[s$vals<min_rain]=0
  # Get time variable (within selected time interval, if available):
  if ((!is.null(date_initial)) & (!is.null(date_final))){
    s$time=rain.df$timestamps[which(
                              rain.df$timestamps==date_initial):which(
                              rain.df$timestamps==date_final)]
  } else {
    s$time=rain.df$timestamps
  }
  
  # add to object "s" a few extra information:
  s$time_resolution=time_resolution
  s$durations=durations 
  s$name_project=name_project
  s$sat_prod=sat_prod
  s$row=pix_row
  s$col=pix_col
  
  # initialize title and filename of plots to be saved.
  if ((!is.null(s$row))|(!is.null(s$col))|(!is.null(s$sat_prod))){
    nfplot=paste0("_",s$sat_prod,'_',s$row,'_',s$col,".png")
    title_plot=paste0('(',s$sat_prod,', ',year(s$time[1]),'-',
                      year(tail(s$time,1)),', pix:',s$row,'-',s$col,')')
  } else {
    nfplot=paste0("_", s$name_project,".png")
    title_plot=paste0('(', s$name_project,', ', 
                      year(s$time[1]),' - ',year(tail(s$time,1)),')')
  }
  
  if (flag_saveALL){
    # plot precipitation timeseries:
    plot_rain_series(time=s$time, 
                     vals=s$vals, 
                     title_x="Year",
                     title_y="Rainfall (mm/h)",
                     main_title=paste0("Rainfall time series (removed values <",
                                        min_rain,"mm)\n",title_plot),
                     path2save=paste0(dir_results,"/rain_timeseries",nfplot))
  }  
  
  #^****************************************************************************
  #^* Just verification: aggregate to daily rain and plot time series and AMAX
  #^****************************************************************************
  flag_check_daily=F        # FOR DEBUGGING ONLY !!!!
  if (flag_check_daily){
    message('checking AMAX series on daily aggregates...')
    precip_boulder=data.frame(date=s$time,
                              values=s$vals)
    daily_sum_precip<-precip_boulder %>%
      mutate(day=as.POSIXct(date,format="%Y-%m-%d")) %>%
      group_by(day) %>% # group by the day column
      # compute SUM of precipitation for each day
      summarise(total_precip=sum(values)) %>%
      na.omit()
    df_daily_p=as.data.frame(daily_sum_precip)
    tmp_array_aggr_daily=df_daily_p$total_precip
    timestamps_ALL_daily=df_daily_p$day
    plot_rain_series(time=timestamps_ALL_daily,
                     vals=tmp_array_aggr_daily,
                     title_x="Year",
                     title_y="Rainfall (mm/day)",
                     main_title=paste0("Daily Rainfall time series mm/day \n",
                                       title_plot),
                     path2save=paste0(dir_results,"/rain_timeseries_daily_",
                                      nfplot))
    # Amax:
    Pmax=c(); years=unique(substring(timestamps_ALL_daily,1,4))
    for (i in 1:length(years)){
      days=which(substring(timestamps_ALL_daily,1,4)==years[i])
      Pdays=tmp_array_aggr_daily[days]
      Pmax[i]=max(Pdays,na.rm=T)
    }
    plot_rain_series(time=as.Date(years,format="%Y"),
                     vals=Pmax,
                     title_x="Year",
                     title_y="Rainfall (mm/day)",
                     main_title=paste0("Daily AMAX Rainfall mm/day \n",
                                       title_plot),
                     path2save=paste0(dir_results,"/rain_AMAX_daily_",nfplot))
  }
  
  
  
  
  
  ##############################################################################
  # SEPARATION defined by 'separation' DRY time intervals:
  ##############################################################################
  # Use function from F. Marra (Unipd, June 2024).
  # Defines start and end of all events separated by at least tot minutes 
  # ('separation'), 
  # with a filter to remove ordinary events across non-subsequent years
  # it provides also s$yr (the covariate in years for the ns smev analysis)  
  s=event_separation_dry_spell(s=s,
                               separation=separation_in_min/s$time_resolution, 
                               min_ev_duration=min_ev_duration)
  #^****************************************************************************
  # Defines ordinary events and AMS for each duration:
  #^****************************************************************************
  # define the blocks (years in this case)
  blocks_id=s$yr
  # analysed years of available data
  blocks=unique(blocks_id)
  # analysed years of all data (including nan)
  blocks_all=unique(year(s$time))
  # number of years
  M=length(blocks_all) 
  # Initialize Annual Maxima with zeros:
  AMS=matrix(0,M,length(durations)) # zeros(M,length(durations)) 
  # add M and blocks to object s
  s$M=M
  s$blocks=blocks
  
  # intialization of all results lists for current station/cell and all durations:
  s_tmp_SMEV=s_tmp_MEV=s$ord_events=list()
  ncount=lmfit=n_ordev=NULL
  quantileplot_smev_ALL=quantileplot_mev_ALL=spaghettiplot_smev_ALL=list()
  spaghettiplot_mev_ALL=bayesfactorplot_ALL=list()
  bayesfactorplot_ALL_bis=modelselectionplot_ALL=list()
  trendAMAXplot_ALL=trend_n_plot_ALL<-list()

  
  # Cycle on all selected durations (e.g.,  1h, 3h, 6h, 24h):
  ##############################################################################
  for (dur in 1:length(durations)){
  ##############################################################################
    # rainfall duration
    message('')
    print('##################')
    print(paste0('Duration ',durations[dur]/60,'h'))
    print('##################')
    # create directory for duration "dur"
    dir.create(paste0(dir_results,"/",durations[dur]/60,"h"),recursive=T)
    dir.res_pix=paste0(dir_results,"/",durations[dur]/60,"h")
    
    #^**************************************************************************
    # get all ordinary events for duration "dur":
    #^**************************************************************************
    s$ord_events[[dur]]=rep(NA,length(s$from))
    iid=durations[dur]/s$time_resolution
    for (ev in 1:length(s$from)){
      event=pracma::conv(s$vals[s$from[ev]:s$to[ev]],rep(1,iid))  #,'same')
      s$ord_events[[dur]][ev]=max(event)
      # ib=which.max(event)
      # s$ord_time[[dur]][ev]=s$time[s$from[ev]+ib-1]
    }
    
    # convert rainfall values to [mm/h] from [mm]
    s$ord_events[[dur]]=s$ord_events[[dur]]*60/durations[dur]
  
    # Compute the annual maxima for duration "dur":
    for (yy in 1:M){   
      if(any(blocks_id==blocks_all[yy])){
        AMS[yy,dur]=max(s$ord_events[[dur]][blocks_id==blocks_all[yy]]) # AMAX
      } else {
        AMS[yy,dur]=NA
        message(paste0('ATTENTION !!! year ',blocks_all[yy],' is missing !!!'))
      }
    }
    
    # define t,y variables for smev:
    t=0;y=0;
    t=s$yr-min(s$yr) # to start time from 0.
    y=s$ord_events[[dur]]
  
    # Apply "left-censoring" threshold for the Weibull distr. (e.g., 90%):
    # We are mostly interested in the right tail of the distribution:  
    # quantile fun (we use type 5 to obtain value consistent with matlab fun)
    thr=quantile(y,thr_leftcens,type=5) 
    thr=thr[[1]]
    # update object s:
    s$thr[[dur]]=thr
    s$AMS=AMS
    s$t=t

    
    #^**************************************************************************
    # Apply trend analysis on the AMAX series (MK and Sens's Slope tests):
    #^**************************************************************************
    # if stat. significant trend is estimated then update and switch on the 
    # flag for more in-depth analysis and saving all results:
    if (flag_saveALL){
      trendAMAX=trend_AMAX(s=s, 
                           alpha=0.05, # significance level 
                           dir.res_trend=dir.res_pix,
                           dur=dur)
      # save trend results to object "s":
      s$signif_trend[[dur]] = trendAMAX$signif_trend
      s$res_MK[[dur]]       = trendAMAX$res_MK
      s$senss[[dur]]        = trendAMAX$senss
      s$sensslope[[dur]]    = trendAMAX$senss$estimates[["Sen's slope"]]
      s$flag_saveALL[[dur]] = trendAMAX$flag_save
      s$flag_like[[dur]]    = trendAMAX$flag_like
      trendAMAXplot_ALL     = c(trendAMAXplot_ALL,list(trendAMAX$trendAMAXplot))
    }
    
    
    #^**************************************************************************
    # BAYESIAN INFERENCE (using adaptive metropolis):
    #^**************************************************************************
    message(paste0("Statistical modelling:"))
    # 4 models currently implemented:
    # - Full nonstationary SMEV model
    # - Partial nonstationary SMEV model
    # - Stationary SMEV model
    # - Stationary GEV model
    #
    # 4 differnt parameters estimation methods currently implemented:
    # - "cmdstanr" Bayesian estimation with Hamiltonian MCMC of STAN package
    # - "stanr" Bayesian estimation with Hamiltonian MCMC of STAN package (old)
    # - "adaptMCMC" Bayesian estimation with adaptive Metropolis MCMC
    # - "fminsearch" Maximum Likelihood Method (minimization function) 
    #^**************************************************************************
    # FULL NON-STATIONARY SMEV MODEL:
    #^**************************************************************************
    if (flag_smev_ns){
      message(paste0("\n--> NON-STATIONARY smev model:"))
      # available algorithms: "cmdstan", "rstan", "adaptMCMC"", "fminsearch" 
      # call function from "module_mod_smev_ns.R"
      s$s_ns_tmp[[dur]]= smev_ns(algorithm  = estim_algorithm,  
                                 Nmcmc      = Nmcmc,
                                 Nmcmc_max  = Nmcmc_max, # needed only for "adaptMCMC" method
                                 nburn      = nburn,
                                 Nslim      = Nslim, 
                                 scales     = scales, # needed only for "adaptMCMC" method
                                 priors_id  = priors[[dur]],
                                 acc.rate   = acc.rate, # needed only for "adaptMCMC" method
                                 y          = y, 
                                 t          = t, 
                                 thr        = thr,
                                 dir.res_id = paste0(dir.res_pix,'/ns/'),
                                 dur        = dur,
                                 flag_save  = flag_saveALL)
    }
  
    
    if (flag_smev_partns){ 
      #^************************************************************************
      # PARTIAL NON-STATIONARY SMEV MODEL:
      #^************************************************************************
      # SMEV with only scale dependent on time
      message(paste0("\n--> PART-NONSTATIONARY smev model:"))
      s$s_partns_tmp[[dur]]= smev_partns(algorithm  = estim_algorithm,
                                         Nmcmc      = Nmcmc,
                                         Nmcmc_max  = Nmcmc_max,
                                         nburn      = nburn,
                                         Nslim      = Nslim,
                                         scales     = scales,
                                         priors_id  = priors[[dur]],
                                         acc.rate   = acc.rate,
                                         y          = y,
                                         t          = t,
                                         thr        = thr,
                                         dir.res_id = paste0(dir.res_pix,'/partial_ns/'),
                                         dur        = dur,
                                         flag_save  = flag_saveALL)
    }
    
    
    #^**************************************************************************
    # STATIONARY SMEV MODEL:
    #^**************************************************************************
    if (flag_smev_stat){
      message(paste0("\n--> STATIONARY smev model:"))
      # call function from "module_mod_smev_stat.R"
      s$s_stat_tmp[[dur]] = smev_stat(algorithm = estim_algorithm,  
                                      Nmcmc     = Nmcmc,
                                      Nmcmc_max = Nmcmc_max,
                                      nburn     = nburn,
                                      Nslim     = Nslim, 
                                      scales    = scales,
                                      priors_id = priors[[dur]],
                                      acc.rate  = acc.rate,
                                      y         = y, 
                                      t         = t, 
                                      thr       = thr,
                                      dir.res_id= paste0(dir.res_pix,'/stat/'),
                                      dur       = dur,
                                      flag_save = flag_saveALL)
    }
        
    
    #^**************************************************************************
    # STATIONARY GEV MODEL:
    #^**************************************************************************
    if (flag_gev_stat){
      message(paste0("\n--> STATIONARY gev model:"))
      # call function from "module_mod_gev_stat.R"
      s$s_gev_tmp[[dur]]= gev_stat(algorithm    = estim_algorithm,  
                                   Nmcmc        = Nmcmc,
                                   Nmcmc_max    = Nmcmc_max,
                                   nburn        = nburn,
                                   Nslim        = Nslim,
                                   scales       = scales_gev,
                                   priors_gev   = priors_gev,
                                   acc.rate     = acc.rate,
                                   y            = na.omit(AMS[,dur]),
                                   dir.res_id   = paste0(dir.res_pix,'/gev/'),
                                   dur          = dur,
                                   flag_save    = flag_saveALL)
    }
    
    
    
    
    
    #^**************************************************************************
    # STATISTICAL SIGNIFICANCE / MODEL COMPARISON:
    #^**************************************************************************
    if (flag_model_compar){
      criteria = model_comparison(s               = s, 
                                  flag_smev_ns    = flag_smev_ns, 
                                  flag_smev_stat  = flag_smev_stat,
                                  flag_smev_partns= flag_smev_partns,
                                  flg_save        = T,
                                  dur             = dur, 
                                  sat_prod        = s$sat_prod,
                                  pix_row         = s$row,
                                  pix_col         = s$col,
                                  dir.res_pix     = dir.res_pix,
                                  nfplot          = nfplot,
                                  title_plot      = title_plot)
      modelselectionplot_ALL=c(modelselectionplot_ALL,
                               list(criteria$modelselectionplot))
      bayesfactorplot_ALL=c(bayesfactorplot_ALL, 
                            list(criteria$bayesfactorplot[[1]]))
      #bayesfactorplot_ALL_bis=c(bayesfactorplot_ALL_bis, 
      #                         list(criteria$bayesfactorplot[[2]]))
    }  
  
    
    
    #^**************************************************************************
    # Count number of ordinary events in each year:
    #^**************************************************************************
    t.df=data.frame(t=s$t)
    s$agg_tbl[[dur]]<-t.df %>% group_by(t) %>% 
                      summarise(total_count=n(),.groups='drop')
    s$ncount[[dur]]=s$agg_tbl[[dur]]$total_count
    s$years=unique(t)
    s$years_unique=blocks  # unique(year(s$time))        

    
    
    #^**************************************************************************
    # MK trend and Linear regression of n:
    #^**************************************************************************
    trendn=trend_n(s_obj         = s, 
                   alpha         = 0.05, # Assign significance level of the test
                   dir.res_trend = dir.res_pix,
                   dur           = dur,
                   nfplot        = nfplot,
                   title_plot    = title_plot)
    
    # save trend in n to object "s":
    s$signif_trend_n[[dur]] = trendn$signif_trend_n
    s$n_MK[[dur]]           = trendn$n_MK
    s$n_senss[[dur]]        = trendn$n_senss
    s$n_sensslope[[dur]]    = trendn$n_sensslope
    s$n_lmfit[[dur]]        = trendn$n_lmfit
    trend_n_plot_ALL        = c(trend_n_plot_ALL, list(trendn$trend_n_plot))
  
    
    
    #^**************************************************************************
    # GET HALF PERIOD ORDINARY EVENTS:
    #^**************************************************************************
    s$halfperiod=ceil(tail(s$years,1)/2)
    # half period year
    s$year_halfperiod=year(s$time[1])+s$halfperiod 
    # n at half period (e.g., 13 year for period 1998-2023):
    s$n_ordev[[dur]]=s$n_lmfit[[dur]]$coefficients[[1]]+
               s$n_lmfit[[dur]]$coefficients[[2]]*s$halfperiod
    # n at half period
    s$n_ordev_halfperiod[[dur]]=s$n_ordev[[dur]]
    # all ordinary events of the half period year (e.g., 13 year for 1998-2023):
    s$ord_eve_halfperiod[[dur]]=y[which(s$t==s$halfperiod)]
    # to do next (case if halfperiod year is nan ?!?????????)
  
    
    
    #^**************************************************************************
    # Compute Empirical Quantiles
    #^**************************************************************************
    # PPF= i/(N+1)   for all distributions.
    # more in general the fomrula is = (i-b)/(N+1-2b)
    
    # # find exceedance probability
    precip_exceed=(1:length(na.omit(s$AMS[,dur])))/(
      length(na.omit(s$AMS[,dur]))+1)
    precip_period=1/(1-precip_exceed)
    # dataframe with AMS and T:
    s$empir_qnt_df[[dur]]=data.frame(x=precip_period, 
                                     y=sort(na.omit(s$AMS[,dur])))
    # # points(precip_period, sort(ord_eve_13),pch=19)
    # # Hazen formula for the non-exceedance empirical frequency:
    # freq=((1:length(s$AMS[,dur])) - 0.5)/length(s$AMS[,dur]) 
    # # Empirical Return periods
    # Temp=1/(1-freq)   
    # points(Temp,sort(AMS),pch=19)
    
    
    
    #^**************************************************************************
    # compute SMEV quantiles for the target return periods (for a given year)
    #^**************************************************************************
    # NB:
    # - All inferred models must have same number of posterior draws:
    # - At least posterior draws from smev ns are needed!
    print("compute quantiles...")
    s_tmp_SMEV[[dur]] = compute_quantiles_SMEV(s                = s,
                                               estim_algorithm  = estim_algorithm,
                                               target_T         = target_T,
                                               flg_save         = T,
                                               id               = dur,
                                               filter_spaghetti = 10,
                                               flag_smev_stat   = flag_smev_stat,
                                               flag_smev_partns = flag_smev_partns,
                                               flag_smev_ns     = flag_smev_ns,
                                               flag_gev_stat    = flag_gev_stat)
    s_tmp_SMEV[[dur]]$duration=paste0(durations[dur]/60, 'h')
    
    # spaghetti plot (SMEV halfperiod):
    spaghettiplot_smev=plot_spaghetti_smev(flg_save    = T,
                                           s           = s,
                                           quant       = s_tmp_SMEV[[dur]],
                                           id          = dur,
                                           dir.res_pix = dir.res_pix)        
    spaghettiplot_smev_ALL=c(spaghettiplot_smev_ALL, list(spaghettiplot_smev))
    
    if (estim_algorithm != "fminsearch"){
      # 90% uncertainty plot (SMEV halfperiod):
      quantileplot_smev=plot_quantiles_smev(flg_save    = T,
                                            s           = s,
                                            quant       = s_tmp_SMEV[[dur]],
                                            id          = dur,
                                            dir.res_pix = dir.res_pix)
      quantileplot_smev_ALL= c(quantileplot_smev_ALL, list(quantileplot_smev))
    }
    # end loop on durations !!!
  }
  # save the AMS array
  s$AMS = AMS
  # update object s with quantiles for each duration:
  s=c(s, 
      quantiles_smev=list(s_tmp_SMEV), 
      quantiles_mev=list(s_tmp_MEV))
  
  
  

  
  ##############################################################################
  #                    Multi-plots with all durations:
  ##############################################################################
  
  #^****************************************************************************
  # BAYES FACTOR:
  #^****************************************************************************
  if (length(bayesfactorplot_ALL)>0){
    # plot BF values:
    plot_quantiles_gridplot(plotlist=bayesfactorplot_ALL,
                            title=paste0('Bayes Factor on SMEV models \n',
                                         title_plot),
                            pathout=paste0(dir_results,"/bayesfactor_smev",
                                           nfplot))
  }
  
  #^****************************************************************************
  # BIC/DIC/...:
  #^****************************************************************************
  if (length(modelselectionplot_ALL)>0){
    plot_quantiles_gridplot(plotlist=modelselectionplot_ALL,
                            title=paste0('Model selection criteria on SMEV \n',
                                         title_plot),
                            pathout=paste0(dir_results,"/criteria_smev",nfplot))
  }
  
  #^****************************************************************************
  # TREND ON AMAX:
  #^****************************************************************************
  if (flag_saveALL){
    plot_quantiles_gridplot(plotlist=trendAMAXplot_ALL,
                            title=paste0('Trends on AMAX rainfall \n',title_plot),
                            pathout=paste0(dir_results,"/trends_AMAX",nfplot))
  }
  #^****************************************************************************
  # TREND ON n:
  #^****************************************************************************
  plot_quantiles_gridplot(plotlist=trend_n_plot_ALL,
                          title=paste0('Trends on number of ordinary events, n\n',
                                       title_plot),
                          pathout=paste0(dir_results,"/trends_n",nfplot))
  
  #^****************************************************************************
  # QUANTILES:
  #^****************************************************************************
  if (flag_saveALL==T){
    if (length(spaghettiplot_smev_ALL)>=4){
    
      # SMEV
      #######
      # quantiles spaghetti:
      plot_quantiles_gridplot(plotlist=spaghettiplot_smev_ALL,
                              title=paste0('Quantiles non-stationary SMEV (',
                                           s$year_halfperiod, ') - spaghetti \n',
                                           title_plot),
                              pathout=paste0(dir_results,"/SMEV_spagh_quant_", 
                                             s$year_halfperiod, nfplot))
      # quantiles with uncertainty ribbon:
      plot_quantiles_gridplot(plotlist=quantileplot_smev_ALL,
                              title=paste0('Quantiles non-stationary SMEV (', 
                                           s$year_halfperiod, ') \n',title_plot),
                              pathout=paste0(dir_results,"/SMEV_uncert_quant_", 
                                             s$year_halfperiod, nfplot))
      
      # # MEV (NOT AVAILABLE YET):
      # ##########################
      # # spaghetti
      # title1=paste0('Quantiles non-stationary MEV-spaghetti\n',title_plot) 
      # pathout1=paste0(dir_results,"/MEV_spaghetti_post_quant_",nfplot)
      # plot_quantiles_gridplot(plotlist=spaghettiplot_mev_ALL,
      #                         title=title1,
      #                         pathout=pathout1)
      # # Quantiles with uncertainty ribbon:
      # plot_quantiles_gridplot(plotlist=quantileplot_mev_ALL,
      #                         title=paste0('Quantiles non-stationary MEV \n', 
      #                                      title_plot),
      #                         pathout=paste0(dir_results,"/MEV_uncert_quant_",
      #                                        nfplot))
    } else {
      print("not all durations have flag_save=T")
    }
  }      
    
  
  ##############################################################################
  #                     SAVING RESULTS INTO .RDS FILE:
  ##############################################################################
  # store object s in the list Sn:
  # Sn[[idx]] = s
  if (flag_save2rds){
    if ((!is.null(s$row))|(!is.null(s$col))|(!is.null(s$sat_prod))){
      # save object s to rds file for pixel idx:
      saveRDS(s,paste0(dir_results,"/s_",s$sat_prod,"_pix_",row,"_",col,".rds"))
    } else {
      saveRDS(s,paste0(dir_results,"/s_",s$name_project,".rds"))
    }
  }  
  # remove object s from memory 
  # rm(s)
  # compute computation time
  end_t=Sys.time()
  print(paste0("comp.time=",end_t-start_t))
  # return the main object "s" containing all results:
  return(s=s)
}





