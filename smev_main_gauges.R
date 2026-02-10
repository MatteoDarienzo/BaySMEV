################################################################################
#                               BaySMEV                                        # 
################################################################################
#                 R version of SMEV method for extreme analysis         	     #
################################################################################
# Main authors: Matteo Darienzo                                                #
# Contributors: Francesco Marra                                                #
# Contacts: matteo.darienzo@unipd.it, francesco.marra@unipd.it                 #
# Institute: University of Padova (Italy), Project INTENSE                     #
# Year:  2025                                                                  #
# version: v. 1.0                                                              #
# last update: 09/07/2025		                                                   #
#------------------------------------------------------------------------------#
#                        Main objective of the program                         #
#------------------------------------------------------------------------------#
# This program performs the SMEV method (Marra et al., 2020) for rainfall      #
# extreme analysis, which is a simplified version of the MEV approach (Marani  #
# and Ignaccolo,2015). This method makes use of all ordinary events, not only  #
# a small sample of extremes (e.g., annual maxima).                            #
# More specifically, it uses a 2-parameter Weibull distribution with the       #
# left-Censoring approach.                                                     #
# Different parameter estimation approaches are available here:                #
# - Maximum Likelihood estimation method with the fminsearch() minimisation    #
#   function as in the original matlab version of the method (please, see      #
#   https://zenodo.org/records/11934843).                                      #
# - Maximum likelihood or Posterior sampling within a Bayesian framework using #
#   an adaptive Metropolis algorithm.                                          #
# - Maximum likelihood or Posterior sampling within a Bayesian framework using #
#   an Hamiltonian MCMC algorithm with Rstan/cmdstanr R packages.              #
# Moreover, both stationary and non-stationary frameworks are available.       #
#                                                                              #
#------------------------------------------------------------------------------#
#                               Acknowledgements                               #
#------------------------------------------------------------------------------#
# For the extreme analysis Uses methods:                                       #
# - MEV (Marani and Ignaccolo, 2015)                                           #
# - SMEV (Marra et al., 2019. 2020)                                            # 
# - Uses open source codes of F. Marra (University of Padua) for SMEV method   #
#   and left-censoring, https://zenodo.org/records/11934843 (SMEV methodology) #
#   and https://zenodo.org/records/15047817 (non-stationary framework)         #
#                                                                              #
#                                                                              #  
# For the Bayesian and MCMC algoryhtm Uses ideas from methods:                 #
# - RsTooDS of Benjamin Renard (INRAE), https://zenodo.org/records/5075760     #
# - BaM of Benjamin Renard (INRAE), https://github.com/BaM-tools/RBaM          #
# - BayDERS of Matteo Darienzo(INRAE),https://github.com/MatteoDarienzo/BayDERS#
# - package "rstan" and "cmdstanr" for stan Bayesian inference and diagnostics.#
# - package "adaptMCMC" for adaptive metropolis algorithm.                     # 
#                                                                              #
#                                                                              #
# Funds:                                                                       #
# Supported by the INTENSE project (raINfall exTremEs and their impacts:       #
# from the local to the National ScalE) funded by the European Union -NextGenEU#
# within PRIN (Progetti di ricerca di Rilevante Interesse Nazionale) programme #
# (grant 2022ZC2522).                                                          #
#------------------------------------------------------------------------------#
#                              Disclaimer                                      #
#------------------------------------------------------------------------------#
# Please, notice that is an experimental code. Further analysis is required for#
# validating the proposed tools. We are not responsible for any loss (e.g., of #
# data, profits, business, customers, orders, other costs or disturbances)     #
# derived by their use in the operational practice. The authorized user accepts#
# the risks associated with using the software given its nature of free        #
# software. It is reserved to expert Users (developers or professionals) having#
# knowledge in computer science, hydrology, extreme analysis and statistics.   #
#------------------------------------------------------------------------------#
#                                READ ME!!!                                    #
#------------------------------------------------------------------------------#
# This is a main launcher, which makes use of several other functions/modules  #
# which must be priorly imported/compiled. The main function is "smev_main()". #
# Moreover, please make sure you installed the following required packages:    #
# "methods", "lattice", "gridExtra",                                           #
# "reshape","reshape2", "ggplot2", "extrafont", "grid", "gtable", "chron",     #
# "RColorBrewer", "cowplot", "psych",  "mosaicData", "tidyr", "lubridate",     #
# "parallel", "latex2exp", "scales", "dplyr", "png","pracma", "ncdf4","CFtime",#
# "rnaturalearth", "raster", "maptools", "Kendall", "trend", "tidyverse",      #
# "ggpubr", "sf","coda", "mcmc","evd","stats", "R.matlab", "extRemes","imager",#
# "mmand", "adaptMCMC", "PEIP", "posterior", "loo", "bayesplot", "rstan")      #
################################################################################





################################################################################
#             PREAMBLE (install packages and load modules)                     #
################################################################################
# Get the main directory folder automatically:         
if(!is.null(dev.list())) dev.off() # Clear plots
cat("\014") # Clear console
rm(list=ls())# Clean workspace
#Libraries used  (it automatically detects and installs missing packages):
pack = c("rstudioapi",
         "methods", "lattice", "gridExtra", "reshape","reshape2", "ggplot2",
         "extrafont", "grid", "gtable", "chron",  "RColorBrewer", "cowplot",
         "psych",  "mosaicData", "tidyr", "lubridate", "changepoint", "parallel",
         "latex2exp", "scales", "dplyr", "png","pracma", 
         "ncdf4", "CFtime", "rnaturalearth", "raster", "maptools", 
         "Kendall", "trend", "tidyverse", "ggpubr", "sf", "coda",
         "mcmc", "evd", "stats", "R.matlab", "extRemes", "imager",  "mmand",
         "adaptMCMC", "PEIP", "posterior", "loo", "bayesplot", "rstan")
install_pack <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE )
    }
  }
}

#  Then try/install packages:
install_pack( pack ) 
# set the directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dir_code <- getwd()
#dir_code <- '/home/idrologia/share/MATTEO/BayDERS/'
setwd(dir_code)
# define the directory with the modules files:
dir.modules = paste0(dir_code, "/modules")
# Include all modules for computation:
sapply(paste0(dir.modules,"/", 
              list.files(paste0(dir_code, "/modules"))), source, chdir=F)





################################################################################
#                       inputs and settings
################################################################################
case_study='Italy' # case study name (the same name of the case study folder !!!): 
name_project='BiasCorrection' # name of results subfolder
time_resolut_prod='Hourly' # time resolution of input rainfall series
min_rain=.1 # [mm] minimum value to define a time interval as 'wet'
separation_in_min=1440 # [min] length of the dry periods to separate storms (e.g. 360)
min_ev_duration=30 # [min] minimum duration of a storm
target_T=sort(c(exp(seq(log(1.1),log(50),.1)),2,5,10,20,50,100)) # [year] target return periods   
durations=c(1,3,6,24)*60  # [min] durations to be examined 
time_resolution=60 # [min] time resolution of the input data (vals)
thr_leftcens=0.9 # percentile threhsold for Weibull Left-censoring
Nmcmc=4000  # total number of MCMC simulations (initial number)
start_values=c(0.7,0,5,0) # starting values for mcmc sampling.
flag_like=F  # T/F perform or not the maximum likelihood method.
flag_save=T  # compute all quantiles and comparison between models.
dir_input=paste0(dir_code,'/data_input/',name_project) # results input data
dir_results=paste0('C:/Users/MATTEO/Documents/PRIN/Results/BaySMEV/',
                   case_study,'/',name_project,'/thr_',thr_leftcens) # results directory

# ONLY if estimation method "adaptMCMC" is chosen:
acc.rate=0.234  # (for "adaptMCMC")MCMC acceptance rate
scales=c(0.02,0.002,0.1,0.005) # (for "adaptMCMC") scaling of mcmc jumps for nonstat smev.
scales_gev=c(0.1,0.1,0.1) # (for "adaptMCMC") scaling of mcmc jumps for stat gev.
Nmcmc_max=40000 # (for "adaptMCMC") total number of MCMC simulations (max number)
Nslim=1 # (for "adaptMCMC") slim mcmc for statistics (to reduce lag-autocorrelation)
nburn=0.5 # (for "adaptMCMC") burn the first portion of the simulation to reduce impact of starting value




################################################################################
# PRIORS SPECIFICATION (BAYESIAN INFERENCE):
################################################################################
# Available distrib:
# "normal", "uniform", "flatprior+", "lognormal"
# with the configuration c(par1, par2, startvalue).
# For this application we apply weakly informative priors:
# Priors for nonstat (4 param) and stat SMEV: 
# Priors for all analized rainfall durations (1h, 3h, 6h, 24h):
# create a list for each duration. 
# For each list create another list for the prior of each parameter 
# (4 in case of nonstationary SMEV with Weibull distrib.)
# from  Wilson and Tuomi, 2005, GRL
priors=list() # initialise priors list for SMEV params.
#1h
priors[[1]]=list(shape_intercept =list("lognormal",log(0.7),0.4,0.7),
                 shape_slope     =list("normal", 0,0.05,0),
                 scale_intercept =list("lognormal",1.45-0.7*log(1),0.5,exp(1.45-0.7*log(1))),    
                 scale_slope     =list("normal",0,0.1039-0.0039*(1),0))
#3h
priors[[2]]=list(shape_intercept =list("lognormal",log(0.7),0.4,0.7),
                 shape_slope     =list("normal",0,0.05,0),
                 scale_intercept =list("lognormal",1.45-0.7*log(3),0.5,exp(1.45-0.7*log(3))),
                 scale_slope     =list("normal",0,0.1039-0.0039*(3),0))
#6h
priors[[3]]=list(shape_intercept =list("lognormal",log(0.7),0.4,0.7),
                 shape_slope     =list("normal",0,0.05,0),
                 scale_intercept =list("lognormal",1.45-0.7*log(6),0.5,exp(1.45-0.7*log(6))),
                 scale_slope     =list("normal",0,0.1039-0.0039*(6),0))
# 24h
priors[[4]]=list(shape_intercept=list("lognormal",log(0.7),0.4,0.7),
                 shape_slope    =list("normal",0,0.05,0),
                 scale_intercept=list("lognormal",1.45-0.7*log(24),0.5,exp(1.45-0.7*log(24))),
                 scale_slope    =list("normal",0,0.1039-0.0039*(24),0))

# Priors for stat GEV (3 param):
priors_gev<-list(loc=list("uniform",0,100,1), scale=list("uniform",0,100,1), shape=list("normal",0,0.3,0.1))
# shape ==> martins et al. 2000, beta distr (mean=0.1, stdev=0.122)
#       ==> and renard 2008, gaussian (mean=0, sd=0.3), quite poor informative










################################################################################
#                    Reading precipitation data:
################################################################################
# Creation of the folder with results.
dir.create(dir_results,recursive=T)
ISO_dirpath='E:/save_lenovo_unipd_20250711/Documents/PRIN/Data/Rain_gauges_data/Rain_Gauges_NOQC/'
ISO_regions=read.csv2(paste0(ISO_dirpath,'/ISO_IT_REGION.csv'),
                      sep=',',header=T) 
metdata_dirpath='E:/save_lenovo_unipd_20250711/Documents/PRIN/Data/Rain_gauges_data/Rain_Gauges_NOQC/METADATA'
# get list of all metadata files:
all_metadata_files=list.files(metdata_dirpath)
stat_after2001_ALL=stat_after2010_ALL=c()
stat_period10_ALL=stat_period15_ALL=stat_period20_ALL=c()
all_stations=c()
thresh_prec=200 # upper limit of precipitation values [mm/h]
missingfract=0.2 #% of data with zeros/no data for each year. 
dir.results='E:/save_lenovo_unipd_20250711/Documents/PRIN/Data/Rain_gauges_data/Rain_Gauges_NOQC/ANALYSIS/'
flag_plot=F
# metadata:
for (reg in 1:length(all_metadata_files)){
  message('###########################################')
  message(paste0('Analizing Region ',reg))
  metadata_reg_i=read.csv2(paste0(metdata_dirpath,'/',all_metadata_files[reg]),
                           sep=',',header=T)
  Region_i=metadata_reg_i$Region[1] 
  message(Region_i)
  message('###########################################')
  num_stat_i=nrow(metadata_reg_i)
  message(paste0(num_stat_i,' stations'))
  Year_Start_i=metadata_reg_i$Year_Start
  timezone_i=metadata_reg_i$TimeZone[1]
  Year_Start_i_datetime=as.POSIXct(Year_Start_i,
                                   format='%Y-%m-%d %H:%M:%S',tz=timezone_i)
  Year_End_i=metadata_reg_i$Year_End
  Year_End_i_datetime=as.POSIXct(Year_End_i,
                                 format='%Y-%m-%d %H:%M:%S',tz=timezone_i)
  # concatenate ALL STATIONS METADATA:
  all_stations=rbind(all_stations,metadata_reg_i)
}





################################################################################
# Cycle over all the stations:
################################################################################
# initialization
rm(s,Sn)
Sn<-vector(mode="list",length=nrow(all_stations))
stat_series_list=list()
for (st in 1:nrow(all_stations)){
  message('\n###########################################')
  message(paste0('Analizing Station ',st))
  Region_s=all_stations$Region[st]
  Code_s=all_stations$Code[st]
  stat_file=all_stations$File[st]
  message(paste0(Code_s,' (', Region_s,')'))
  message('###########################################')
  ISO_s=ISO_regions$ISO[ISO_regions$REGION==Region_s]
  stat_s_dirpath=paste0('E:/save_lenovo_unipd_20250711/Documents/PRIN/Data/Rain_gauges_data/Rain_Gauges_NOQC/DATA/',
                        Region_s)
  stat_series_list[[st]]=list()
  stat_series_list[[st]]$Region=Region_s
  stat_series_list[[st]]$File=stat_file
  
  #load series:
  rm(series)
  series=read.csv2(paste0(stat_s_dirpath,'/',stat_file),sep=',',header=T)
  # filter out the anomalies:
  stat_series_list[[st]]$series_filter=series[which(as.numeric(series$Prec_h)<thresh_prec),]
  Year_Start_i=all_stations$Year_Start[st]
  Year_End_i=all_stations$Year_End[st]
  timezone=all_stations$TimeZone[st]
  Year_Start_i_datetime=as.POSIXct(Year_Start_i,
                                   format='%Y-%m-%d %H:%M:%S',tz=timezone)
  Year_End_i_datetime=as.POSIXct(Year_End_i,
                                 format='%Y-%m-%d %H:%M:%S',tz=timezone)
  stat_series_list[[st]]$yrs=c()
  stat_series_list[[st]]$yrs=unique(year(stat_series_list[[st]]$series_filter$Datetime_h))
  message(paste0(length(stat_series_list[[st]]$yrs),' yrs'))

  #^****************************************************************************
  # check years with missing data (at daily scale):
  #^****************************************************************************
  # aggregate to daily scale to reduce computational time:
  precip_boulder=data.frame(date=stat_series_list[[st]]$series_filter$Datetime_h, 
                            values=as.numeric(stat_series_list[[st]]$series_filter$Prec_h))
  daily_sum_precip=c()
  daily_sum_precip<-precip_boulder %>%
    mutate(day=as.POSIXct(date,format="%Y-%m-%d")) %>%
    group_by(day) %>% # group by the day column
    summarise(daily_precip=sum(values)) %>%  # compute SUM of all precip on each day
    na.omit()
  df_daily_p=as.data.frame(daily_sum_precip)
  
  yrs_rm=c()
  for (yr in 1:length(stat_series_list[[st]]$yrs)){
    yrs_idx=yrs_Prec=Prec_nonnull=Prec_nonnan=c()
    yrs_idx=which(year(df_daily_p$day)==stat_series_list[[st]]$yrs[yr])
    yrs_Prec=as.numeric(df_daily_p$daily_precip[yrs_idx])
    Prec_nonnull=length(which(yrs_Prec>0))
    Prec_nonnan=length(which(!is.na(yrs_Prec)))
    
    # check if at least 80 % are available (at daily scale):
    if ((Prec_nonnull<= missingfract*365) | 
        (Prec_nonnan <= missingfract*365) |
        (length(yrs_Prec) <= missingfract*365)){
      yrs_rm[yr]=T
    } else {
      yrs_rm[yr]=F
    }
  }
  
  stat_series_list[[st]]$yrs_rm=c()
  stat_series_list[[st]]$yrs_rm=yrs_rm
  
  # remove yrs with less 10% data:
  stat_series_list[[st]]$yrs_ok=c()
  if (any(stat_series_list[[st]]$yrs_rm==T)){
    stat_series_list[[st]]$yrs_ok=stat_series_list[[st]]$yrs[-c(
      which(stat_series_list[[st]]$yrs_rm==T))]
  } else {
    stat_series_list[[st]]$yrs_ok=stat_series_list[[st]]$yrs
  }
  # stat_series_list[[st]]$totyrs=round((Year_End_i_datetime-Year_Start_i_datetime)/365)
  stat_series_list[[st]]$totyrs=length(stat_series_list[[st]]$yrs_ok)
  message(paste0(stat_series_list[[st]]$totyrs,
                 ' yrs ok (<',missingfract*100,'% missing data)'))
  
  
  # plot:
  flag_plot=F
  if (flag_plot){ 
    plot_timeseries_s=ggplot()+
      geom_line(data=stat_series_list[[st]]$series_filter,
                aes(x=as.POSIXct(Datetime_h,format="%Y-%m-%d %H:%M:%S",tz='UTC'),
                    y=as.numeric(Prec_h)),
                color="blue",linewidth=0.2)+
      xlab("Year")+
      theme_bw(base_size=15)+
      # annotate("text",
      #          x= as.POSIXct(stat_series_list[[st]]$series_filter$Datetime_h, 
      #                              format="%Y-%m-%d %H:%M:%S",tz='UTC')+days(1000),
      #          y= max(as.numeric(stat_series_list[[st]]$series_filter$Prec_h))-10,
      #          label=paste0("=",round(as.numeric(dist)/1000,0)," km \n",
      #                       "DEM=",DEM1," [m a.s.l.]\n",
      #                       "MAP=",round(as.numeric(MAP1),0)," mm"),
      #          color="blue",size=5)+
      scale_y_continuous(name='Precipitation [mm]')+
      ggtitle(paste0(all_stations$Code[st]," (",all_stations$Region[s],
                     ") ---- coord: ", 
                     round(as.numeric(all_stations$Lat[st]),3),"°N, ",
                     round(as.numeric(all_stations$Lon[st]),3),"°E; ",
                     round(as.numeric(all_stations$Elevation[st]),0),
                     " m a.s.l; hourly, ",
                     "from ",as.Date(Year_Start_i_datetime),
                     " to ",as.Date(Year_End_i_datetime),
                     " (",stat_series_list[[st]]$totyrs," yrs <",
                     missingfract*100,"% blanks in ", 
                     length(stat_series_list[[st]]$yrs)," yrs)"))
    print(plot_timeseries_s)
    # ggsave(plot_timeseries_s, 
    #        filename=paste0(dir.results,all_stations$Region[st],'_',
    #                        all_stations$Code[st],'.png'),
    #        width=20,height=10,dpi=100)
  }
  
  dir.res_p=paste0(dir_results,"/s_",st,"_",
                   all_stations$Code[st],"_",all_stations$Region[st], "/")
  dir.create(dir.res_p,recursive=T)
  

  
  # Begin SMEV analysis
  ##########################################################################
  print("- Extract rainfall timeseries...")
  ##########################################################################
  start_t=Sys.time()
  # Extract the time series data for the selected station:
  prec_year=array(F,dim=length(stat_series_list[[st]]$series_filter$Datetime_h))
  for (ii in 1:length(stat_series_list[[st]]$series_filter$Datetime_h)){
    if (any(stat_series_list[[st]]$yrs_ok==year(stat_series_list[[st]]$series_filter$Datetime_h[ii]))){
      prec_year[ii]=T
    }  
  }
  Prec_ts_pixel.df=data.frame(y=as.numeric(stat_series_list[[st]]$series_filter$Prec_h[which(prec_year)]),
                              year=year(stat_series_list[[st]]$series_filter$Datetime_h[which(prec_year)]),
                              timestamps=stat_series_list[[st]]$series_filter$Datetime_h[which(prec_year)]) 
  
  # Remove neg values or below a certain minimum threshold (e.g., 0.01 mm):
  Prec_ts_pixel.df$y[Prec_ts_pixel.df$y<min_rain]=0
  # saveRDS(Prec_ts_pixel.df, '~/PRIN/Scripts/smev_bayesian/example_rainfall_data.rds')
  
  if (max(Prec_ts_pixel.df$y)>0){
    s_all <- smev_analysis_main(rain.df           = Prec_ts_pixel.df,
                                name_project      = "BiasCorrection",
                                dir_results       = dir.res_p,
                                min_rain          = 0.1,
                                time_resolution   = 60,
                                min_ev_duration   = 30,
                                separation_in_min = 1440,
                                duration          = c(1,3,6,24)*60,
                                target_T          = target_T,
                                thr_leftcens      = 0.9,                   
                                estim_algorithm   = "cmdstan",
                                priors            = priors,                               
                                priors_gev        = priors_gev,                           
                                Nmcmc             = 4000,                        
                                Nmcmc_max         = 40000,                     
                                acc.rate          = 0.234,                       
                                Nslim             = 1,                           
                                nburn             = 0.5,                          
                                scales            = c(0.02,0.002,0.1,0.005),  
                                scales_gev        = c(0.1,0.1,0.1),
                                pix_row           = st, 
                                pix_col           = st,
                                flag_like         = F,                         
                                flag_saveALL      = T,                        
                                flag_gev_stat     = T,
                                flag_smev_stat    = T,
                                flag_smev_ns      = T,
                                flag_smev_partns  = F,
                                flag_model_compar = T,
                                flag_save2rds     = T)
    # remove from memory object s
    rm(s_all)
    end_t=Sys.time()
    print(paste0("comp.time=",end_t-start_t))
  } else {
    print("Station series with null values !!!")
  }
  ###############################
  # END OF ANALYSIS ON STATION s
  ###############################
}


