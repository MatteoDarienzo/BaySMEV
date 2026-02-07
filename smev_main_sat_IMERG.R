################################################################################
#                                  SMEV                                        # 
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
#             PREAMBLE (install packages and load required R modules)          #
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
sapply(paste0(dir.modules,"/",list.files(paste0(dir_code,"/modules"))),source,chdir=F)






################################################################################
#                       inputs and settings
################################################################################
case_study='Italy' # case study name (the same name of the case study folder !!!): 
sat_prod='IMERG_2001_2024' # name of results subfolder
time_resolut_prod='Hourly' # time resolution of input rainfall series
date_final="2024-12-31 23:00:00" # final date for the analysis
date_initial="2001-01-01 00:00:00" # initial date for the analysis
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
filter_fact=1
filter_fact_2=1
dir_input='C:/Users/MATTEO/Documents/PRIN/Data/Satellite_data/1_hour/' # path of satellite data
filename_input='IMERG_italy_hourly_2001_2024_aggreg_0.20deg.rds' # filename of satellite data
# dir_results=paste0('C:/Users/MATTEO/Documents/PRIN/Results/BaySMEV/',
#                    case_study,'/',sat_prod)  # results directory
dir_results=paste0('E:/save_lenovo_unipd_20250711/Documents/PRIN/Results/BaySMEV/',
                   case_study,'/',sat_prod)  # results directory





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
priors_gev<-list(loc=list("uniform",0,100,1),
                 scale=list("uniform",0,100,1),
                 shape=list("normal",0,0.3,0.1))
# shape ==> martins et al. 2000, beta distr (mean=0.1, stdev=0.122)
#       ==> and renard 2008, gaussian (mean=0, sd=0.3), quite poor informative







################################################################################
#                Reading satellite precipitation data:
################################################################################
# Creation of the folder with results.
dir.create(dir_results,recursive=T)
rm(rain_data)
rain_data<-readRDS(paste0(dir_input,'/',filename_input))





################################################################################
# Cycle over all the pixels:
################################################################################
# initialization
Prec_ts_pixel=lat_tmp=lon_tmp=c()
# define time limits for plots:
limits.X<-format(strptime(as.character(c(rain_data$timestamps_ALL_hourly[1], 
                                        rain_data$timestamps_ALL_hourly[length(
                                          rain_data$timestamps_ALL_hourly)])),
                         "%Y-%m-%d"),"%d/%m/%Y")
start_col=1
start_row=1
end_col=ncol(rain_data$tmp_array_aggr_ALL_hourly)
end_row=nrow(rain_data$tmp_array_aggr_ALL_hourly)

# cycle over pixels (IMERG: 53 x 49 x 210384)
################################
for (row in start_row:end_row) {
  for (col in start_col:end_col) {
##################################
    print("#################################################")
    print(c(row,col))
    print("#################################################")
    dir.res_p=paste0(dir_results,"/thr_",thr_leftcens,"/pix_",row,"_",col)
    dir.create(dir.res_p,recursive =T)
    rm(s_all) # remove previous results object, if any.
    
    # Choose if to analyse all pixels or just one every 'filter_fact' pixels 
    # (for computational reasons).
    # For each selected pixel then choose if to compute the full analysis,
    # of basic analysis only, by specifying filter_fact_2.
    # if not specified these filters are set to=1, meaning all pixels are analysed.
    if (!exists('filter_fact')){filter_fact=1}
    if (!exists('filter_fact_2')){filter_fact_2=1}
    
    if ((row==1 & col==1)|((row %% filter_fact==0) & (col %% filter_fact==0))){ 
      # Decide if for this pixel you want to compute all calculations or the main only:
      if ((row==1 & col==1)|((row %% filter_fact_2==0) & (col %% filter_fact_2==0))){ 
        # Performs all calculations:
        # - max likelihood mcmc
        # - max posterior mcmc (with priors)
        # - compute all quantiles and comparison
        flag_save=T
        flag_like=F  # for now, let's skip the comparison likelihood vs posterior !!!
      } else {
        # Perform only a few calculations:
        # - max posterior mcmc (with priors)
        # - a few quantiles (2, 5, 10, 20, 50, 100)
        flag_save=F
        flag_like=F
      }
      
      
      ##########################################################################
      message("- Extract rainfall timeseries...")
      ##########################################################################
      # Extract the time series data for the selected cell from 3D array:
      time_window =which(rain_data$timestamps_ALL_hourly==date_initial):which(
        rain_data$timestamps_ALL_hourly==date_final)
      Prec_ts_pixel=rain_data$tmp_array_aggr_ALL_hourly[row,col,time_window] %>% squeeze
      Prec_ts_pixel[Prec_ts_pixel==as.character(rain_data$fillvalue$value)]=NA
      Prec_ts_pixel.df=data.frame(y=Prec_ts_pixel,
                                  year=year(rain_data$timestamps_ALL_hourly[which(
                                     rain_data$timestamps_ALL_hourly==date_initial):which(
                                     rain_data$timestamps_ALL_hourly==date_final)]))
      Prec_ts_pixel.df$y=Prec_ts_pixel.df$y   # orig unit is in meters. 
      # Remove neg values or below a certain minimum threshold (e.g., 0.01 mm):
      Prec_ts_pixel.df$y[Prec_ts_pixel.df$y<min_rain]=0
      Prec_ts_pixel.df$timestamps_ALL_hourly=rain_data$timestamps_ALL_hourly[time_window]
      # plot(Prec_ts_pixel.df$year, Prec_ts_pixel.df$y)
      # saveRDS(Prec_ts_pixel.df, '~/PRIN/Scripts/smev_bayesian/example_rainfall_data.rds')
      
      if (max(Prec_ts_pixel.df$y)>0){
        # Perform SMEV (both nonstat and stat) and GEV (stat) analysis:
        # (save results into .rds file for current pixel).
        s_all <- smev_analysis_main(rain.df           = Prec_ts_pixel.df,
                                    name_project      = sat_prod,
                                    dir_results       = dir.res_p,
                                    date_initial      = date_initial, 
                                    date_final        = date_final,
                                    min_rain          = min_rain,
                                    time_resolution   = time_resolution,
                                    min_ev_duration   = min_ev_duration,
                                    separation_in_min = separation_in_min,
                                    duration          = durations,
                                    target_T          = target_T,
                                    thr_leftcens      = thr_leftcens,                   
                                    estim_algorithm   = "cmdstan",
                                    priors            = priors,                               
                                    priors_gev        = priors_gev,                           
                                    Nmcmc             = Nmcmc,
                                    pix_row           = row, 
                                    pix_col           = col,
                                    sat_prod          = sat_prod,
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
      } else {
        print("Pixel with null values !!! Please, check.")
      }
      ##########################
      # END OF ANALYSIS ON PIXEL
      ##########################
    } else {
      # keep empty element for this current pixel idx
      # Sn[[idx]] = list()
      
    }
  }
}








