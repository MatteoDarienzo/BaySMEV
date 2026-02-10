################################################################################
#                         PLOTS/MAPS OF BaySMEV                                # 
################################################################################
#                   Bayesian-SMEV (extreme analysis)  	                       #
################################################################################
# Main author:   Matteo Darienzo (matteo.darienzo@unipd.it)		           	     #
# Institure: University of Padova (Italy), Project INTENSE                     #
# Year:     2025                                                               #
# version: v. 1.0                                                              #
# last update: 06/06/2025		                                                   #
#------------------------------------------------------------------------------#
#                        Main objective of the program                         #
#------------------------------------------------------------------------------#
# This is a post-processing module of BaySMEV:                                 #
# It analizes results from BaySMEV (SMEV method, Marra et al.,2020)for rainfall#
# extreme analysis, which is a simplified version of the MEV approach (Marani  #
# and Ignaccolo,2015). This method makes use of all ordinary events, not only  #
# for a small sample of extremes (e.g., annual maxima).                        #
# More specifically, it uses a 2-parameter Weibull distribution with the       #
# left-Censoring approach.                                                     #
# Different estimation approaches are available:                               #
# - Maximum Likelihood estimation method with the fminsearch() optimisation    #
#   function as in the original matlab version of the method (please, see      #
#   https://zenodo.org/records/11934843).                                      #
# - Maximum likelihood or Posterior sampling within a Bayesian framework using #
#   an adaptive Metropolis algorithm.                                          #
# - Maximum likelihood or Posterior sampling within a Bayesian framework using #
#   an Hamiltonian MCMC algorithm with Rstan/cmdstanr R packages.              #
# Moreover, both stationary and non-stationary frameworks are available.       #
#                                                                              #
# This specific module uploads all previous results from a gridded dataset     #
# (e.g., satellite product) where for each pixel a .rds data file is available #
# Then it computes maps for different variables and parameters, and of their   #
# changes, such as: SMEV/GEV parameters estimates (maxpost, MLE) and CI,       #
# quantiles estimates for different return periods, and related decadal changes#
# in %. All of these for different durations (e.g., 1h,3h,6h,24h).             # 
#                                                                              #
#------------------------------------------------------------------------------#
#                               Acknowledgements                               #
#------------------------------------------------------------------------------#
# For the Bayesian and MCMC algoryhtm Uses ideas from methods:                 #
# - RsTooDS of Benjamin Renard (INRAE), https://zenodo.org/records/5075760     #
# - BaM of Benjamin Renard (INRAE), https://github.com/BaM-tools/RBaM          #
# - BayDERS of Matteo Darienzo(INRAE),https://github.com/MatteoDarienzo/BayDERS#
# - package rstan and cmdstanr for stan Bayesian inference and diagnostics     #
# - package adaptMCMC for adaptive metropolis algorithm.                       # 
#                                                                              #
# For the extreme analysis Uses methods:                                       #
# - MEV (Marani and Ignaccolo, 2015)                                           #
# - SMEV (Marra et al., 2019. 2020)                                            # 
# - Uses open source codes of F. Marra (University of Padua) for SMEV method   #
#   and left-censoring, https://zenodo.org/records/11934843 (SMEV methodology) #
#   and https://zenodo.org/records/15047817 (non-stationary framework)         #
#                                                                              #
# Funds:                                                                       #
# Supported by the INTENSE project (raINfall exTremEs and their impacts:       #
# from the local to the National ScalE) funded by the European Union -NextGenEU#
# within PRIN (Progetti di ricerca di Rilevante Interesse Nazionale) programme #
# (grant 2022ZC2522).                                                          #
#------------------------------------------------------------------------------#
#                                Disclaimer                                    #
# -----------------------------------------------------------------------------#
# Please, notice that is an experimental code. Further analysis is required for#
# validating the proposed tools. We are not responsible for any loss (of data, #
# profits, business, customers, orders, other costs or disturbances) derived by#
# their use in the operational practice. The authorized user accepts the risks #
# associated with using the software given its nature of free software. It is  #
# reserved to expert Users (developers or professionals) having knowledge in   #
# computer science, hydrology, extreme analysis and statistics.                #
################################################################################







################################################################################
#             PREAMBLE (install packages and load required R modules)          #
################################################################################
# Get the main directory folder automatically:         
if(!is.null(dev.list())) dev.off() # Clear plots
cat("\014") # Clear console
rm(list=ls())# Clean workspace
#Libraries used  (it automatically detects and installs missing packages):
pack=c("rstudioapi",
       "methods","lattice","gridExtra","reshape","reshape2","ggplot2","extrafont",
       "grid","gtable","chron","coda","RColorBrewer","cowplot","viridis",
       "psych","mosaicData","tidyr","lubridate","Kendall","trend","tidyverse", 
       "latex2exp","dplyr","png","pracma","CFtime","ggpubr",
       "raster","maptools","sf","terra","exactextractr","ncdf4","rnaturalearth")
install_pack<-function(x){
  for(i in x){
    #  require returns TRUE invisibly if it was able to load package
    if(!require(i,character.only=T)){
      #  If package was not able to be loaded then re-install
      install.packages(i,dependencies=T)
      #  Load package after installing
      require(i,character.only=T)
    }
  }
}
#  Then try/install packages:
install_pack(pack)
# set the directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dir_code<-getwd()
setwd(dir_code)
setwd('../')
# define the directory with the modules files:
dir.modules=paste0(getwd(),"/modules")
# Include all modules for computation:
sapply(paste0(dir.modules,"/",list.files(paste0(getwd(),"/modules"))),source,chdir=F)
setwd(dir_code)






################################################################################
#                       inputs and settings
################################################################################
case_study='Italy' # case study name (the same name of the case study folder !!!): 
sat_prod='IMERG_2001_2024' # name of results subfolder for the dataset
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
dir_input='C:/Users/MATTEO/Documents/PRIN/Data/Satellite_data/1_hour/' # folder path of satellite data
filename_input='IMERG_italy_hourly_2001_2024_aggreg_0.20deg.rds' # filename of satellite data
# dir_results=paste0('C:/Users/MATTEO/Documents/PRIN/Results/BaySMEV/',case_study,'/',sat_prod)  # results directory
dir_results=paste0('E:/save_lenovo_unipd_20250711/Documents/PRIN/Results/BaySMEV/',case_study,'/',sat_prod)  # SMEV results directory (where .rds data files are located)

sat_acron='IMERG' # acronyme of satellite products for plot titles
dir_input_RData='E:/save_lenovo_unipd_20250711/Documents/PRIN/Scripts/smev_bayesian/Case_studies/Italy/Results/IMERG_2001_2024/thr_0.9/Results_IMERG_italy_1hr_2001_2024_smevns_NEW.RData'
dir_output_RData='E:/save_lenovo_unipd_20250711/Documents/PRIN/Scripts/smev_bayesian/Case_studies/Italy/Results/IMERG_2001_2024/thr_0.9/Results_IMERG_italy_1hr_2001_2024_smevns_NEW.RData'
dir.res_maps=paste0(dir_results,"/thr_",thr_leftcens,"/maps/")
dir_shp_mask='C:/Users/MATTEO/Documents/PRIN/Qgis/georef-italy-regione/georef-italy-regione-millesime.shp' # shp with admin limits
dir_shp_world='C:/Users/MATTEO/Documents/PRIN/Qgis/WorldGeog/WorldGeog.shp'
dir_shp_user='E:/save_lenovo_unipd_20250711/Documents/PRIN/Qgis/geog_areas_paper2.shp' # some specific polygons for statistics
flag_read_only=T











################################################################################
#                     LOADING Rainfall product and set CRS:
################################################################################
# load processed satellite data from .rds file just to get the grid info, then it
# will be removed from memory in a few lines.
rain_data<-readRDS(paste0(dir_input,'/',filename_input))

# Get rain product grid (check for proper field name):
if ("tmp_array_aggr_ALL_hourly" %in% names(rain_data)){
  tmp_array=rain_data$tmp_array_aggr_ALL_hourly[,,1]
} else {
  tmp_array=rain_data$tmp_array_ALL_hourly[,,1]
}
# Get lat and lon (check for proper field name):
if ("lon_aggr" %in% names(rain_data)){
  lon=rain_data$lon_aggr
  lat=rain_data$lat_aggr
} else {
  lon=rain_data$lon
  lat=rain_data$lat
} 

# fix the crs (wgs84)
crswgs84=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
# read shp file with regions limits of Italy:
test<-read_sf(dsn=dir_shp_mask) 
# world map as sp class
gg_world=ne_countries(scale="medium",returnclass="sf")
# read shp file with World country limits:
land_shp<-read_sf(dsn=dir_shp_world)
# create lat/lon grid for map plot:
lat_tmp=lon_tmp=c()
for (row in 1:nrow(tmp_array)){
  for (col in 1:ncol(tmp_array)){
    lat_tmp=c(lat_tmp,lat[col])
    lon_tmp=c(lon_tmp,lon[row]) 
  }
}
# plot the indeces of the grid:
plot_map_indexes(tmp_array=tmp_array,
                 lon_tmp=lon_tmp,
                 lat_tmp=lat_tmp,
                 mask_shp=test,
                 mask_world=gg_world,
                 path_output=paste0(dir_results,"/map_",sat_prod,"_indexes.png"),
                 crs_map=crswgs84,
                 title_x="Longitude (°E)",
                 title_y="Latitude (°N)",
                 breaks_x=c(6.00001,10,14,18),
                 breaks_y=c(36.00001,40,44,47.9999),
                 labels_x=c("6°E","10°E","14°E","18°E"),
                 labels_y=c("36°N","40°N","44°N","48°N"),
                 limits_x=c(6,19),
                 limits_y=c(36,48))
# Remove ancillary memory
rm(rain_data)
gc() # free the memory






################################################################################
#               LOADING RESULTS FOR EACH PIXEL FROM .RDS FILES:
################################################################################
# load initialization module:
# this module inizializes all lists for the maps of different variables and  
# parameters that are computed or loaded hereafter:
source(paste0(dir_code,"/module_initialize_maps.R"))

if (flag_read_only==F){
  # Load extra info from the pixel 1,1:
  ss<-readRDS(paste0(dir_results,"/",sat_prod,"/thr_",thr_leftcens, 
               "/pix_",1,"_",1,"/s_",sat_prod,"_pix_",1,"_",1,".rds"))
  years_smev=ss$years
  halfperiod=ss$halfperiod
  Tgrid=ss$quantiles_smev[[1]]$Tgrid
  year_halfperiod=ss$year_halfperiod
  rm(ss)
  gc() 
  
  
  # run loop on the pixels:
  #^****************************************************************************
  for (row in 1:nrow(tmp_array)){
    #^**************************************************************************
    for (col in 1:ncol(tmp_array)){
      #^************************************************************************
      print(c(row,col))
      pix=pix+1
      if ((row==1&col==1)|((row %% filter_fact==0)&(col %% filter_fact==0))){
        dir.res_pix=paste0(dir_results,"/",sat_prod,"/thr_",thr_leftcens,"/pix_",row,"_",col)
        
        if (file.exists(paste0(dir.res_pix,"/s_",sat_prod,"_pix_",row,"_",col,".rds"))){
          # load rds data file for selected pixel:
          ss<-readRDS(paste0(dir.res_pix,"/s_",sat_prod,"_pix_",row,"_",col,".rds"))
          
          #^********************************************************************
          # NonStationary SMEV:
          #^********************************************************************
          a_w_quant5_Post_1h[pix]=ss$s_ns_tmp[[1]]$a_w_quant_Post[[1]]
          a_w_quant50_Post_1h[pix]=ss$s_ns_tmp[[1]]$a_w_quant_Post[[2]]
          a_w_quant95_Post_1h[pix]=ss$s_ns_tmp[[1]]$a_w_quant_Post[[3]]
          a_w_quant5_Post_3h[pix]=ss$s_ns_tmp[[2]]$a_w_quant_Post[[1]]
          a_w_quant50_Post_3h[pix]=ss$s_ns_tmp[[2]]$a_w_quant_Post[[2]]
          a_w_quant95_Post_3h[pix]=ss$s_ns_tmp[[2]]$a_w_quant_Post[[3]]
          a_w_quant5_Post_6h[pix]=ss$s_ns_tmp[[3]]$a_w_quant_Post[[1]]
          a_w_quant50_Post_6h[pix]=ss$s_ns_tmp[[3]]$a_w_quant_Post[[2]]
          a_w_quant95_Post_6h[pix]=ss$s_ns_tmp[[3]]$a_w_quant_Post[[3]]
          a_w_quant5_Post_24h[pix]=ss$s_ns_tmp[[4]]$a_w_quant_Post[[1]]
          a_w_quant50_Post_24h[pix]=ss$s_ns_tmp[[4]]$a_w_quant_Post[[2]]
          a_w_quant95_Post_24h[pix]=ss$s_ns_tmp[[4]]$a_w_quant_Post[[3]]
          
          b_w_quant5_Post_1h[pix]=ss$s_ns_tmp[[1]]$b_w_quant_Post[[1]]
          b_w_quant50_Post_1h[pix]=ss$s_ns_tmp[[1]]$b_w_quant_Post[[2]]
          b_w_quant95_Post_1h[pix]=ss$s_ns_tmp[[1]]$b_w_quant_Post[[3]]
          b_w_quant5_Post_3h[pix]=ss$s_ns_tmp[[2]]$b_w_quant_Post[[1]]
          b_w_quant50_Post_3h[pix]=ss$s_ns_tmp[[2]]$b_w_quant_Post[[2]]
          b_w_quant95_Post_3h[pix]=ss$s_ns_tmp[[2]]$b_w_quant_Post[[3]]
          b_w_quant5_Post_6h[pix]=ss$s_ns_tmp[[3]]$b_w_quant_Post[[1]]
          b_w_quant50_Post_6h[pix]=ss$s_ns_tmp[[3]]$b_w_quant_Post[[2]]
          b_w_quant95_Post_6h[pix]=ss$s_ns_tmp[[3]]$b_w_quant_Post[[3]]
          b_w_quant5_Post_24h[pix]=ss$s_ns_tmp[[4]]$b_w_quant_Post[[1]]
          b_w_quant50_Post_24h[pix]=ss$s_ns_tmp[[4]]$b_w_quant_Post[[2]]
          b_w_quant95_Post_24h[pix]=ss$s_ns_tmp[[4]]$b_w_quant_Post[[3]]
          
          a_C_quant5_Post_1h[pix]=ss$s_ns_tmp[[1]]$a_C_quant_Post[[1]]
          a_C_quant50_Post_1h[pix]=ss$s_ns_tmp[[1]]$a_C_quant_Post[[2]]
          a_C_quant95_Post_1h[pix]=ss$s_ns_tmp[[1]]$a_C_quant_Post[[3]]
          a_C_quant5_Post_3h[pix]=ss$s_ns_tmp[[2]]$a_C_quant_Post[[1]]
          a_C_quant50_Post_3h[pix]=ss$s_ns_tmp[[2]]$a_C_quant_Post[[2]]
          a_C_quant95_Post_3h[pix]=ss$s_ns_tmp[[2]]$a_C_quant_Post[[3]]
          a_C_quant5_Post_6h[pix]=ss$s_ns_tmp[[3]]$a_C_quant_Post[[1]]
          a_C_quant50_Post_6h[pix]=ss$s_ns_tmp[[3]]$a_C_quant_Post[[2]]
          a_C_quant95_Post_6h[pix]=ss$s_ns_tmp[[3]]$a_C_quant_Post[[3]]
          a_C_quant5_Post_24h[pix]=ss$s_ns_tmp[[4]]$a_C_quant_Post[[1]]
          a_C_quant50_Post_24h[pix]=ss$s_ns_tmp[[4]]$a_C_quant_Post[[2]]
          a_C_quant95_Post_24h[pix]=ss$s_ns_tmp[[4]]$a_C_quant_Post[[3]]
          
          b_C_quant5_Post_1h[pix]=ss$s_ns_tmp[[1]]$b_C_quant_Post[[1]]
          b_C_quant50_Post_1h[pix]=ss$s_ns_tmp[[1]]$b_C_quant_Post[[2]]
          b_C_quant95_Post_1h[pix]=ss$s_ns_tmp[[1]]$b_C_quant_Post[[3]]
          b_C_quant5_Post_3h[pix]=ss$s_ns_tmp[[2]]$b_C_quant_Post[[1]]
          b_C_quant50_Post_3h[pix]=ss$s_ns_tmp[[2]]$b_C_quant_Post[[2]]
          b_C_quant95_Post_3h[pix]=ss$s_ns_tmp[[2]]$b_C_quant_Post[[3]]
          b_C_quant5_Post_6h[pix]=ss$s_ns_tmp[[3]]$b_C_quant_Post[[1]]
          b_C_quant50_Post_6h[pix]=ss$s_ns_tmp[[3]]$b_C_quant_Post[[2]]
          b_C_quant95_Post_6h[pix]=ss$s_ns_tmp[[3]]$b_C_quant_Post[[3]]
          b_C_quant5_Post_24h[pix]=ss$s_ns_tmp[[4]]$b_C_quant_Post[[1]]
          b_C_quant50_Post_24h[pix]=ss$s_ns_tmp[[4]]$b_C_quant_Post[[2]]
          b_C_quant95_Post_24h[pix]=ss$s_ns_tmp[[4]]$b_C_quant_Post[[3]]
          
          a_w_MAP_ns_1h[pix]=ss$s_ns_tmp[[1]]$theta_maxpost[1] # maxpost shape intercept
          a_w_MAP_ns_3h[pix]=ss$s_ns_tmp[[2]]$theta_maxpost[1]
          a_w_MAP_ns_6h[pix]=ss$s_ns_tmp[[3]]$theta_maxpost[1]
          a_w_MAP_ns_24h[pix]=ss$s_ns_tmp[[4]]$theta_maxpost[1] 
          a_w_fmins_ns_1h[pix]=ss$s_ns_tmp[[1]]$theta_maxlik_fmins[1] # fminsearch method
          a_w_fmins_ns_3h[pix]=ss$s_ns_tmp[[2]]$theta_maxlik_fmins[1] 
          a_w_fmins_ns_6h[pix]=ss$s_ns_tmp[[3]]$theta_maxlik_fmins[1] 
          a_w_fmins_ns_24h[pix]=ss$s_ns_tmp[[4]]$theta_maxlik_fmins[1] 
          
          b_w_MAP_ns_1h[pix]=ss$s_ns_tmp[[1]]$theta_maxpost[2] # maxpost shape slope
          b_w_MAP_ns_3h[pix]=ss$s_ns_tmp[[2]]$theta_maxpost[2] 
          b_w_MAP_ns_6h[pix]=ss$s_ns_tmp[[3]]$theta_maxpost[2] 
          b_w_MAP_ns_24h[pix]=ss$s_ns_tmp[[4]]$theta_maxpost[2] 
          b_w_fmins_ns_1h[pix]=ss$s_ns_tmp[[1]]$theta_maxlik_fmins[2] # fminsearch method
          b_w_fmins_ns_3h[pix]=ss$s_ns_tmp[[2]]$theta_maxlik_fmins[2] 
          b_w_fmins_ns_6h[pix]=ss$s_ns_tmp[[3]]$theta_maxlik_fmins[2] 
          b_w_fmins_ns_24h[pix]=ss$s_ns_tmp[[4]]$theta_maxlik_fmins[2] 
          
          a_C_MAP_ns_1h[pix]=ss$s_ns_tmp[[1]]$theta_maxpost[3] # maxpost scale intercept
          a_C_MAP_ns_3h[pix]=ss$s_ns_tmp[[2]]$theta_maxpost[3]
          a_C_MAP_ns_6h[pix]=ss$s_ns_tmp[[3]]$theta_maxpost[3] 
          a_C_MAP_ns_24h[pix]=ss$s_ns_tmp[[4]]$theta_maxpost[3]
          a_C_fmins_ns_1h[pix]=ss$s_ns_tmp[[1]]$theta_maxlik_fmins[3] # fminsearch method
          a_C_fmins_ns_3h[pix]=ss$s_ns_tmp[[2]]$theta_maxlik_fmins[3] 
          a_C_fmins_ns_6h[pix]=ss$s_ns_tmp[[3]]$theta_maxlik_fmins[3] 
          a_C_fmins_ns_24h[pix]=ss$s_ns_tmp[[4]]$theta_maxlik_fmins[3] 
          
          b_C_MAP_ns_1h[pix]=ss$s_ns_tmp[[1]]$theta_maxpost[4] # maxpost scale slope
          b_C_MAP_ns_3h[pix]=ss$s_ns_tmp[[2]]$theta_maxpost[4]
          b_C_MAP_ns_6h[pix]=ss$s_ns_tmp[[3]]$theta_maxpost[4]
          b_C_MAP_ns_24h[pix]=ss$s_ns_tmp[[4]]$theta_maxpost[4] 
          b_C_fmins_ns_1h[pix]=ss$s_ns_tmp[[1]]$theta_maxlik_fmins[4] # fminsearch method
          b_C_fmins_ns_3h[pix]=ss$s_ns_tmp[[2]]$theta_maxlik_fmins[4] 
          b_C_fmins_ns_6h[pix]=ss$s_ns_tmp[[3]]$theta_maxlik_fmins[4] 
          b_C_fmins_ns_24h[pix]=ss$s_ns_tmp[[4]]$theta_maxlik_fmins[4] 
          
          #^********************************************************************
          # Stationary SMEV:
          #^********************************************************************
          a_w_stat_quant5_Post_1h[pix]=ss$s_stat_tmp[[1]]$a_w_quant_Post[[1]]
          a_w_stat_quant50_Post_1h[pix]=ss$s_stat_tmp[[1]]$a_w_quant_Post[[2]]
          a_w_stat_quant95_Post_1h[pix]=ss$s_stat_tmp[[1]]$a_w_quant_Post[[3]]
          a_w_stat_quant5_Post_3h[pix]=ss$s_stat_tmp[[2]]$a_w_quant_Post[[1]]
          a_w_stat_quant50_Post_3h[pix]=ss$s_stat_tmp[[2]]$a_w_quant_Post[[2]]
          a_w_stat_quant95_Post_3h[pix]=ss$s_stat_tmp[[2]]$a_w_quant_Post[[3]]
          a_w_stat_quant5_Post_6h[pix]=ss$s_stat_tmp[[3]]$a_w_quant_Post[[1]]
          a_w_stat_quant50_Post_6h[pix]=ss$s_stat_tmp[[3]]$a_w_quant_Post[[2]]
          a_w_stat_quant95_Post_6h[pix]=ss$s_stat_tmp[[3]]$a_w_quant_Post[[3]]
          a_w_stat_quant5_Post_24h[pix]=ss$s_stat_tmp[[4]]$a_w_quant_Post[[1]]
          a_w_stat_quant50_Post_24h[pix]=ss$s_stat_tmp[[4]]$a_w_quant_Post[[2]]
          a_w_stat_quant95_Post_24h[pix]=ss$s_stat_tmp[[4]]$a_w_quant_Post[[3]]
          
          a_C_stat_quant5_Post_1h[pix]=ss$s_stat_tmp[[1]]$a_C_quant_Post[[1]]
          a_C_stat_quant50_Post_1h[pix]=ss$s_stat_tmp[[1]]$a_C_quant_Post[[2]]
          a_C_stat_quant95_Post_1h[pix]=ss$s_stat_tmp[[1]]$a_C_quant_Post[[3]]
          a_C_stat_quant5_Post_3h[pix]=ss$s_stat_tmp[[2]]$a_C_quant_Post[[1]]
          a_C_stat_quant50_Post_3h[pix]=ss$s_stat_tmp[[2]]$a_C_quant_Post[[2]]
          a_C_stat_quant95_Post_3h[pix]=ss$s_stat_tmp[[2]]$a_C_quant_Post[[3]]
          a_C_stat_quant5_Post_6h[pix]=ss$s_stat_tmp[[3]]$a_C_quant_Post[[1]]
          a_C_stat_quant50_Post_6h[pix]=ss$s_stat_tmp[[3]]$a_C_quant_Post[[2]]
          a_C_stat_quant95_Post_6h[pix]=ss$s_stat_tmp[[3]]$a_C_quant_Post[[3]]
          a_C_stat_quant5_Post_24h[pix]=ss$s_stat_tmp[[4]]$a_C_quant_Post[[1]]
          a_C_stat_quant50_Post_24h[pix]=ss$s_stat_tmp[[4]]$a_C_quant_Post[[2]]
          a_C_stat_quant95_Post_24h[pix]=ss$s_stat_tmp[[4]]$a_C_quant_Post[[3]]
          
          a_w_MAP_stat_1h[pix]=ss$s_stat_tmp[[1]]$theta_maxpost[1] # maxpost shape intercept
          a_w_MAP_stat_3h[pix]=ss$s_stat_tmp[[2]]$theta_maxpost[1]
          a_w_MAP_stat_6h[pix]=ss$s_stat_tmp[[3]]$theta_maxpost[1]
          a_w_MAP_stat_24h[pix]=ss$s_stat_tmp[[4]]$theta_maxpost[1] 
          a_w_fmins_stat_1h[pix]=ss$s_stat_tmp[[1]]$theta_maxlik_fmins[1] # fminsearch method
          a_w_fmins_stat_3h[pix]=ss$s_stat_tmp[[2]]$theta_maxlik_fmins[1] 
          a_w_fmins_stat_6h[pix]=ss$s_stat_tmp[[3]]$theta_maxlik_fmins[1] 
          a_w_fmins_stat_24h[pix]=ss$s_stat_tmp[[4]]$theta_maxlik_fmins[1] 
          
          a_C_MAP_stat_1h[pix]=ss$s_stat_tmp[[1]]$theta_maxpost[2] # maxpost scale intercept
          a_C_MAP_stat_3h[pix]=ss$s_stat_tmp[[2]]$theta_maxpost[2]
          a_C_MAP_stat_6h[pix]=ss$s_stat_tmp[[3]]$theta_maxpost[2] 
          a_C_MAP_stat_24h[pix]=ss$s_stat_tmp[[4]]$theta_maxpost[2]
          a_C_fmins_stat_1h[pix]=ss$s_stat_tmp[[1]]$theta_maxlik_fmins[2] # fminsearch method
          a_C_fmins_stat_3h[pix]=ss$s_stat_tmp[[2]]$theta_maxlik_fmins[2] 
          a_C_fmins_stat_6h[pix]=ss$s_stat_tmp[[3]]$theta_maxlik_fmins[2] 
          a_C_fmins_stat_24h[pix]=ss$s_stat_tmp[[4]]$theta_maxlik_fmins[2] 
          
          #^********************************************************************
          # GEV stationary:
          #^********************************************************************
          a_L_gev_quant5_Post_1h[pix]=ss$s_gev_tmp[[1]]$a_L_quant_Post[[1]]
          a_L_gev_quant50_Post_1h[pix]=ss$s_gev_tmp[[1]]$a_L_quant_Post[[2]]
          a_L_gev_quant95_Post_1h[pix]=ss$s_gev_tmp[[1]]$a_L_quant_Post[[3]]
          a_L_gev_quant5_Post_3h[pix]=ss$s_gev_tmp[[2]]$a_L_quant_Post[[1]]
          a_L_gev_quant50_Post_3h[pix]=ss$s_gev_tmp[[2]]$a_L_quant_Post[[2]]
          a_L_gev_quant95_Post_3h[pix]=ss$s_gev_tmp[[2]]$a_L_quant_Post[[3]]
          a_L_gev_quant5_Post_6h[pix]=ss$s_gev_tmp[[3]]$a_L_quant_Post[[1]]
          a_L_gev_quant50_Post_6h[pix]=ss$s_gev_tmp[[3]]$a_L_quant_Post[[2]]
          a_L_gev_quant95_Post_6h[pix]=ss$s_gev_tmp[[3]]$a_L_quant_Post[[3]]
          a_L_gev_quant5_Post_24h[pix]=ss$s_gev_tmp[[4]]$a_L_quant_Post[[1]]
          a_L_gev_quant50_Post_24h[pix]=ss$s_gev_tmp[[4]]$a_L_quant_Post[[2]]
          a_L_gev_quant95_Post_24h[pix]=ss$s_gev_tmp[[4]]$a_L_quant_Post[[3]]
          
          a_w_gev_quant5_Post_1h[pix]=ss$s_gev_tmp[[1]]$a_w_quant_Post[[1]]
          a_w_gev_quant50_Post_1h[pix]=ss$s_gev_tmp[[1]]$a_w_quant_Post[[2]]
          a_w_gev_quant95_Post_1h[pix]=ss$s_gev_tmp[[1]]$a_w_quant_Post[[3]]
          a_w_gev_quant5_Post_3h[pix]=ss$s_gev_tmp[[2]]$a_w_quant_Post[[1]]
          a_w_gev_quant50_Post_3h[pix]=ss$s_gev_tmp[[2]]$a_w_quant_Post[[2]]
          a_w_gev_quant95_Post_3h[pix]=ss$s_gev_tmp[[2]]$a_w_quant_Post[[3]]
          a_w_gev_quant5_Post_6h[pix]=ss$s_gev_tmp[[3]]$a_w_quant_Post[[1]]
          a_w_gev_quant50_Post_6h[pix]=ss$s_gev_tmp[[3]]$a_w_quant_Post[[2]]
          a_w_gev_quant95_Post_6h[pix]=ss$s_gev_tmp[[3]]$a_w_quant_Post[[3]]
          a_w_gev_quant5_Post_24h[pix]=ss$s_gev_tmp[[4]]$a_w_quant_Post[[1]]
          a_w_gev_quant50_Post_24h[pix]=ss$s_gev_tmp[[4]]$a_w_quant_Post[[2]]
          a_w_gev_quant95_Post_24h[pix]=ss$s_gev_tmp[[4]]$a_w_quant_Post[[3]]
          
          a_C_gev_quant5_Post_1h[pix]=ss$s_gev_tmp[[1]]$a_C_quant_Post[[1]]
          a_C_gev_quant50_Post_1h[pix]=ss$s_gev_tmp[[1]]$a_C_quant_Post[[2]]
          a_C_gev_quant95_Post_1h[pix]=ss$s_gev_tmp[[1]]$a_C_quant_Post[[3]]
          a_C_gev_quant5_Post_3h[pix]=ss$s_gev_tmp[[2]]$a_C_quant_Post[[1]]
          a_C_gev_quant50_Post_3h[pix]=ss$s_gev_tmp[[2]]$a_C_quant_Post[[2]]
          a_C_gev_quant95_Post_3h[pix]=ss$s_gev_tmp[[2]]$a_C_quant_Post[[3]]
          a_C_gev_quant5_Post_6h[pix]=ss$s_gev_tmp[[3]]$a_C_quant_Post[[1]]
          a_C_gev_quant50_Post_6h[pix]=ss$s_gev_tmp[[3]]$a_C_quant_Post[[2]]
          a_C_gev_quant95_Post_6h[pix]=ss$s_gev_tmp[[3]]$a_C_quant_Post[[3]]
          a_C_gev_quant5_Post_24h[pix]=ss$s_gev_tmp[[4]]$a_C_quant_Post[[1]]
          a_C_gev_quant50_Post_24h[pix]=ss$s_gev_tmp[[4]]$a_C_quant_Post[[2]]
          a_C_gev_quant95_Post_24h[pix]=ss$s_gev_tmp[[4]]$a_C_quant_Post[[3]]
          
          a_L_MAP_gev_1h[pix]=ss$s_gev_tmp[[1]]$theta_maxpost[1] # maxpost loc intercept
          a_L_MAP_gev_3h[pix]=ss$s_gev_tmp[[2]]$theta_maxpost[1]
          a_L_MAP_gev_6h[pix]=ss$s_gev_tmp[[3]]$theta_maxpost[1]
          a_L_MAP_gev_24h[pix]=ss$s_gev_tmp[[4]]$theta_maxpost[1] 
          
          a_C_MAP_gev_1h[pix]=ss$s_gev_tmp[[1]]$theta_maxpost[2] # maxpost scale intercept
          a_C_MAP_gev_3h[pix]=ss$s_gev_tmp[[2]]$theta_maxpost[2]
          a_C_MAP_gev_6h[pix]=ss$s_gev_tmp[[3]]$theta_maxpost[2] 
          a_C_MAP_gev_24h[pix]=ss$s_gev_tmp[[4]]$theta_maxpost[2]
          
          a_w_MAP_gev_1h[pix]=ss$s_gev_tmp[[1]]$theta_maxpost[3] # maxpost shape intercept
          a_w_MAP_gev_3h[pix]=ss$s_gev_tmp[[2]]$theta_maxpost[3]
          a_w_MAP_gev_6h[pix]=ss$s_gev_tmp[[3]]$theta_maxpost[3] 
          a_w_MAP_gev_24h[pix]=ss$s_gev_tmp[[4]]$theta_maxpost[3]
          
          
          #^********************************************************************
          # Results of trend analysis on "n":
          #^********************************************************************
          n_ordev[pix]       =ss$n_ordev[[1]]
          lmfit_n_inter[pix] =ss$n_lmfit[[1]]$coefficients[[1]]
          lmfit_n_slope[pix] =ss$n_lmfit[[1]]$coefficients[[2]]
          signif_trend_n[pix]=ss$signif_trend_n[[1]]
          n_senss[[pix]]     =ss$n_senss[[1]]$estimates[[1]]
          pval_MK_n[[pix]]   =ss$n_MK[[1]]$p.value
          z_MK_n[[pix]]      =ss$n_MK[[1]]$statistic[[1]]
          
          
          #^********************************************************************
          # Results of trend analysis on AMAX series:
          #^********************************************************************
          AMS[[pix]]=ss$AMS  # AMAX series for each duration at pixel "pix"
          if (!is.null(ss$signif_trend[[1]])){
            signif_trend_1h[pix]  =ss$signif_trend[[1]] # significance T/F trend on AMAX
            signif_trend_3h[pix]  =ss$signif_trend[[2]] # significance T/F trend on AMAX
            signif_trend_6h[pix]  =ss$signif_trend[[3]] # significance T/F trend on AMAX
            signif_trend_24h[pix] =ss$signif_trend[[4]] # significance T/F trend on AMAX
            
            sensslope_AMAX_1h[pix] =ss$sensslope[[1]] # sensslope on AMAX
            sensslope_AMAX_3h[pix] =ss$sensslope[[2]] # sensslope on AMAX
            sensslope_AMAX_6h[pix] =ss$sensslope[[3]] # sensslope on AMAX
            sensslope_AMAX_24h[pix]=ss$sensslope[[4]] # sensslope on AMAX
            
            pval_MK_AMAX_1h =ss$res_MK[[1]]$p.value   # p-value of the MK test
            pval_MK_AMAX_3h =ss$res_MK[[2]]$p.value   # p-value of the MK test
            pval_MK_AMAX_6h =ss$res_MK[[3]]$p.value   # p-value of the MK test
            pval_MK_AMAX_24h=ss$res_MK[[4]]$p.value   # p-value of the MK test
            
            z_MK_AMAX_1h =ss$res_MK[[1]]$statistic[[1]] # z estimate of MK test
            z_MK_AMAX_3h =ss$res_MK[[2]]$statistic[[1]] # z estimate of MK test
            z_MK_AMAX_6h =ss$res_MK[[3]]$statistic[[1]] # z estimate of MK test
            z_MK_AMAX_24h=ss$res_MK[[4]]$statistic[[1]] # z estimate of MK test
            
            # res_MK[[pix]]=ss$res_MK[[1]]  # Mann-Kandall test on AMAX
          }
          
          
          #^********************************************************************
          # Significance (model selection criteria, DIC, BIC,...):
          #^********************************************************************
          # DIC:
          DIC3_ns_1h[pix]=ss$s_ns_tmp[[1]]$DIC3
          DIC3_ns_3h[pix]=ss$s_ns_tmp[[2]]$DIC3
          DIC3_ns_6h[pix]=ss$s_ns_tmp[[3]]$DIC3
          DIC3_ns_24h[pix]=ss$s_ns_tmp[[4]]$DIC3
          DIC3_stat_1h[pix]=ss$s_stat_tmp[[1]]$DIC3
          DIC3_stat_3h[pix]=ss$s_stat_tmp[[2]]$DIC3
          DIC3_stat_6h[pix]=ss$s_stat_tmp[[3]]$DIC3
          DIC3_stat_24h[pix]=ss$s_stat_tmp[[4]]$DIC3
          # BIC:
          BIC_ns_1h[pix]=ss$s_ns_tmp[[1]]$BIC
          BIC_ns_3h[pix]=ss$s_ns_tmp[[2]]$BIC
          BIC_ns_6h[pix]=ss$s_ns_tmp[[3]]$BIC
          BIC_ns_24h[pix]=ss$s_ns_tmp[[4]]$BIC
          BIC_stat_1h[pix]=ss$s_stat_tmp[[1]]$BIC
          BIC_stat_3h[pix]=ss$s_stat_tmp[[2]]$BIC
          BIC_stat_6h[pix]=ss$s_stat_tmp[[3]]$BIC
          BIC_stat_24h[pix]=ss$s_stat_tmp[[4]]$BIC
          # loglikelihood:
          loglik_fmins_stat_1h[pix]=ss$s_stat_tmp[[1]]$maxlik_fmins
          loglik_fmins_stat_3h[pix]=ss$s_stat_tmp[[2]]$maxlik_fmins
          loglik_fmins_stat_6h[pix]=ss$s_stat_tmp[[3]]$maxlik_fmins
          loglik_fmins_stat_24h[pix]=ss$s_stat_tmp[[4]]$maxlik_fmins
          loglik_fmins_ns_1h[pix]=ss$s_ns_tmp[[1]]$maxlik_fmins
          loglik_fmins_ns_3h[pix]=ss$s_ns_tmp[[2]]$maxlik_fmins
          loglik_fmins_ns_6h[pix]=ss$s_ns_tmp[[3]]$maxlik_fmins
          loglik_fmins_ns_24h[pix]=ss$s_ns_tmp[[4]]$maxlik_fmins
          # loglik_bayes_stat_1h[pix]=ss$s_stat_tmp[[1]]$maxlik_mcmc
          # loglik_bayes_stat_3h[pix]=ss$s_stat_tmp[[2]]$maxlik_mcmc
          # loglik_bayes_stat_6h[pix]=ss$s_stat_tmp[[3]]$maxlik_mcmc
          # loglik_bayes_stat_24h[pix]=ss$s_stat_tmp[[4]]$maxlik_mcmc
          # loglik_bayes_ns_1h[pix]=ss$s_ns_tmp[[1]]$maxlik_mcmc
          # loglik_bayes_ns_3h[pix]=ss$s_ns_tmp[[2]]$maxlik_mcmc
          # loglik_bayes_ns_6h[pix]=ss$s_ns_tmp[[3]]$maxlik_mcmc
          # loglik_bayes_ns_24h[pix]=ss$s_ns_tmp[[4]]$maxlik_mcmc
          
          #^********************************************************************
          # Quantiles at half-period year, 
          # (for different return periods and different durations)
          #^********************************************************************
          # 2 yrs:
          qCurve_2yr_maxpost_1h[pix]=ss$quantiles_smev[[1]]$qCurve_maxpost_halfper[which(Tgrid==2)]
          qCurve_2yr_maxpost_3h[pix]=ss$quantiles_smev[[2]]$qCurve_maxpost_halfper[which(Tgrid==2)] 
          qCurve_2yr_maxpost_6h[pix]=ss$quantiles_smev[[3]]$qCurve_maxpost_halfper[which(Tgrid==2)]
          qCurve_2yr_maxpost_24h[pix]=ss$quantiles_smev[[4]]$qCurve_maxpost_halfper[which(Tgrid==2)] 
          # 5 yrs:
          qCurve_5yr_maxpost_1h[pix]=ss$quantiles_smev[[1]]$qCurve_maxpost_halfper[which(Tgrid==5)]
          qCurve_5yr_maxpost_3h[pix]=ss$quantiles_smev[[2]]$qCurve_maxpost_halfper[which(Tgrid==5)] 
          qCurve_5yr_maxpost_6h[pix]=ss$quantiles_smev[[3]]$qCurve_maxpost_halfper[which(Tgrid==5)]
          qCurve_5yr_maxpost_24h[pix]=ss$quantiles_smev[[4]]$qCurve_maxpost_halfper[which(Tgrid==5)] 
          # 10 yrs:
          qCurve_10yr_maxpost_1h[pix]=ss$quantiles_smev[[1]]$qCurve_maxpost_halfper[which(Tgrid==10)]
          qCurve_10yr_maxpost_3h[pix]=ss$quantiles_smev[[2]]$qCurve_maxpost_halfper[which(Tgrid==10)] 
          qCurve_10yr_maxpost_6h[pix]=ss$quantiles_smev[[3]]$qCurve_maxpost_halfper[which(Tgrid==10)]
          qCurve_10yr_maxpost_24h[pix]=ss$quantiles_smev[[4]]$qCurve_maxpost_halfper[which(Tgrid==10)]  
          # 20 yrs:
          qCurve_20yr_maxpost_1h[pix]=ss$quantiles_smev[[1]]$qCurve_maxpost_halfper[which(Tgrid==20)]
          qCurve_20yr_maxpost_3h[pix]=ss$quantiles_smev[[2]]$qCurve_maxpost_halfper[which(Tgrid==20)] 
          qCurve_20yr_maxpost_6h[pix]=ss$quantiles_smev[[3]]$qCurve_maxpost_halfper[which(Tgrid==20)]
          qCurve_20yr_maxpost_24h[pix]=ss$quantiles_smev[[4]]$qCurve_maxpost_halfper[which(Tgrid==20)]  
          # 50 yrs:
          qCurve_50yr_maxpost_1h[pix]=ss$quantiles_smev[[1]]$qCurve_maxpost_halfper[which(Tgrid==50)]
          qCurve_50yr_maxpost_3h[pix]=ss$quantiles_smev[[2]]$qCurve_maxpost_halfper[which(Tgrid==50)] 
          qCurve_50yr_maxpost_6h[pix]=ss$quantiles_smev[[3]]$qCurve_maxpost_halfper[which(Tgrid==50)]
          qCurve_50yr_maxpost_24h[pix]=ss$quantiles_smev[[4]]$qCurve_maxpost_halfper[which(Tgrid==50)] 
          # 100 yrs:
          qCurve_100yr_maxpost_1h[pix]=ss$quantiles_smev[[1]]$qCurve_maxpost_halfper[which(Tgrid==100)]
          qCurve_100yr_maxpost_3h[pix]=ss$quantiles_smev[[2]]$qCurve_maxpost_halfper[which(Tgrid==100)] 
          qCurve_100yr_maxpost_6h[pix]=ss$quantiles_smev[[3]]$qCurve_maxpost_halfper[which(Tgrid==100)]
          qCurve_100yr_maxpost_24h[pix]=ss$quantiles_smev[[4]]$qCurve_maxpost_halfper[which(Tgrid==100)]
          
          #^********************************************************************
          # Median quantiles at half-period year,
          # (for 2yrs return period and different durations)
          #^********************************************************************
          qCurve_2yr_median_1h[pix]=ss$quantiles_smev[[1]]$spaghetti_qnt_df_ns_post$med_ns_post[
            which(ss$quantiles_smev[[1]]$spaghetti_qnt_df_ns_post$x==2)]
          qCurve_2yr_median_3h[pix]=ss$quantiles_smev[[2]]$spaghetti_qnt_df_ns_post$med_ns_post[
            which(ss$quantiles_smev[[2]]$spaghetti_qnt_df_ns_post$x==2)]
          qCurve_2yr_median_6h[pix]=ss$quantiles_smev[[3]]$spaghetti_qnt_df_ns_post$med_ns_post[
            which(ss$quantiles_smev[[3]]$spaghetti_qnt_df_ns_post$x==2)]
          qCurve_2yr_median_24h[pix]=ss$quantiles_smev[[4]]$spaghetti_qnt_df_ns_post$med_ns_post[
            which(ss$quantiles_smev[[4]]$spaghetti_qnt_df_ns_post$x==2)]
          
          
          #^********************************************************************
          # Uncertainty (all spsghettis) at 90%:
          #^********************************************************************
          #^********************************************************************
          #                            1h:
          #^********************************************************************
          dur=1
          all_post_shape0=ss$s_ns_tmp[[dur]]$zPost[,1]
          all_post_shape1=ss$s_ns_tmp[[dur]]$zPost[,2]
          all_post_scale0=ss$s_ns_tmp[[dur]]$zPost[,3]
          all_post_scale1=ss$s_ns_tmp[[dur]]$zPost[,4]
          # shape (SMEV):
          all_change_shape=all_post_shape1*10/(all_post_shape0+all_post_shape1*halfperiod)*100
          unc_shape_1h[pix]=quantile(all_change_shape,c(0.95))[[1]]-quantile(all_change_shape,c(0.05))[[1]]
          # scale (SMEV):
          all_change_scale=all_post_scale1*10/(all_post_scale0+all_post_scale1*halfperiod)*100
          unc_scale_1h[pix]=quantile(all_change_scale,c(0.95))[[1]]-quantile(all_change_scale,c(0.05))[[1]]
          # n (SMEV):
          unc_n_1h[pix]=(ss$n_senss[[dur]]$conf.int[[2]]*10-ss$n_senss[[dur]]$conf.int[[1]]*10)/ 
            (ss$n_lmfit[[dur]]$coefficients[[1]]+ss$n_lmfit[[dur]]$coefficients[[2]]*halfperiod)*100
          # quantile 2 years:
          # compute smev parameters at first, last and 10 year (all mcmc samaple):
          shape_t_10   = all_post_shape0 + all_post_shape1*10   # decade
          shape_t_half = all_post_shape0 + all_post_shape1*halfperiod # last year
          shape_t_last = all_post_shape0 + all_post_shape1*tail(years_smev,1) # last year
          scale_t_10   = all_post_scale0 + all_post_scale1*10   # decade
          scale_t_half = all_post_scale0 + all_post_scale1*halfperiod
          scale_t_last = all_post_scale0 + all_post_scale1*tail(years_smev,1) # last year
          n_0    = ss$n_lmfit[[dur]]$coefficients[[1]]
          n_10   = ss$n_lmfit[[dur]]$coefficients[[1]] + ss$n_lmfit[[dur]]$coefficients[[2]]*10
          n_half = ss$n_lmfit[[dur]]$coefficients[[1]] + ss$n_lmfit[[dur]]$coefficients[[2]]*halfperiod
          n_last = ss$n_lmfit[[dur]]$coefficients[[1]] + ss$n_lmfit[[dur]]$coefficients[[2]]* tail(years_smev,1)
          pr=1-1/2 # TR= 2 years
          all_qnt_2yr_0 = sapply(scale  = all_post_scale0,
                                 shape  = all_post_shape0,
                                 n_ordev= n_0,       # n_ordev,
                                 FUN    = weibull_inv,
                                 X      = pr)
          all_qnt_2yr_half = sapply(scale  = scale_t_half,
                                    shape  = shape_t_half,
                                    n_ordev= n_half,      # n_ordev,
                                    FUN    = weibull_inv,
                                    X      = pr)
          all_qnt_2yr_10 = sapply(scale  = scale_t_10,
                                  shape  = shape_t_10,
                                  n_ordev= n_10,      # n_ordev,
                                  FUN    = weibull_inv,
                                  X      = pr)
          all_qnt_2yr_last = sapply(scale  = scale_t_last,
                                    shape  = shape_t_last,
                                    n_ordev= n_last,  # n_ordev,
                                    FUN    = weibull_inv,
                                    X      = pr)
          n_decades=tail(years_smev,1)/10 # number of decades in the period
          # qnt changes for 2yrs return period:          
          all_change_qnt2yr=(all_qnt_2yr_last-all_qnt_2yr_0)/
            (all_qnt_2yr_10)*100/ n_decades
          # uncertainty at 90%:
          unc_qnt2yrs_1h[pix]=quantile(all_change_qnt2yr,c(0.95))[[1]]-
            quantile(all_change_qnt2yr,c(0.05))[[1]]
          
          #^********************************************************************
          #                            3h:
          #^********************************************************************
          dur=2 # 3h duration
          all_post_shape0=ss$s_ns_tmp[[dur]]$zPost[,1]
          all_post_shape1=ss$s_ns_tmp[[dur]]$zPost[,2]
          all_post_scale0=ss$s_ns_tmp[[dur]]$zPost[,3]
          all_post_scale1=ss$s_ns_tmp[[dur]]$zPost[,4]
          # shape (SMEV):
          all_change_shape=all_post_shape1*10/(all_post_shape0+all_post_shape1*halfperiod)*100
          unc_shape_3h[pix]=quantile(all_change_shape,c(0.95))[[1]]-quantile(all_change_shape,c(0.05))[[1]]
          # scale (SMEV):
          all_change_scale=all_post_scale1*10/(all_post_scale0+all_post_scale1*halfperiod)*100
          unc_scale_3h[pix]=quantile(all_change_scale,c(0.95))[[1]]-quantile(all_change_scale,c(0.05))[[1]]
          # n (SMEV):
          unc_n_3h[pix]=(ss$n_senss[[dur]]$conf.int[[2]]*10-ss$n_senss[[dur]]$conf.int[[1]]*10)/ 
            (ss$n_lmfit[[dur]]$coefficients[[1]]+ss$n_lmfit[[dur]]$coefficients[[2]]*halfperiod)*100
          # Quantile 2 years:
          # compute smev parameters at first, last and 10 year (all mcmc samaple):
          shape_t_10  =all_post_shape0+all_post_shape1*10   # decade
          shape_t_half=all_post_shape0+all_post_shape1*halfperiod # last year
          shape_t_last=all_post_shape0+all_post_shape1*tail(years_smev,1) # last year
          scale_t_10  =all_post_scale0+all_post_scale1*10   # decade
          scale_t_half=all_post_scale0+all_post_scale1*halfperiod
          scale_t_last=all_post_scale0+all_post_scale1*tail(years_smev,1) # last year
          n_0   =ss$n_lmfit[[dur]]$coefficients[[1]]
          n_10  =ss$n_lmfit[[dur]]$coefficients[[1]]+ss$n_lmfit[[dur]]$coefficients[[2]]*10
          n_half=ss$n_lmfit[[dur]]$coefficients[[1]]+ss$n_lmfit[[dur]]$coefficients[[2]]*halfperiod
          n_last=ss$n_lmfit[[dur]]$coefficients[[1]]+ss$n_lmfit[[dur]]$coefficients[[2]]*tail(years_smev,1)
          pr=1-1/2 # TR= 2 years
          all_qnt_2yr_0 = sapply(scale  = all_post_scale0,
                                 shape  = all_post_shape0,
                                 n_ordev= n_0,       # n_ordev,
                                 FUN    = weibull_inv,
                                 X      = pr)
          all_qnt_2yr_half = sapply(scale  = scale_t_half,
                                    shape  = shape_t_half,
                                    n_ordev= n_half,      # n_ordev,
                                    FUN    = weibull_inv,
                                    X      = pr)
          all_qnt_2yr_10 = sapply(scale  = scale_t_10,
                                  shape  = shape_t_10,
                                  n_ordev= n_10,      # n_ordev,
                                  FUN    = weibull_inv,
                                  X      = pr)
          all_qnt_2yr_last = sapply(scale  = scale_t_last,
                                    shape  = shape_t_last,
                                    n_ordev= n_last,  # n_ordev,
                                    FUN    = weibull_inv,
                                    X      = pr)
          # qnt changes for 2yrs return period:          
          all_change_qnt2yr=(all_qnt_2yr_last-all_qnt_2yr_0)/
            (all_qnt_2yr_10)*100/n_decades
          # uncertainty at 90%:
          unc_qnt2yrs_3h[pix]=quantile(all_change_qnt2yr,c(0.95))[[1]]-
            quantile(all_change_qnt2yr,c(0.05))[[1]]
          
          #^********************************************************************
          #                             6h:
          #^********************************************************************
          dur=3
          all_post_shape0=ss$s_ns_tmp[[dur]]$zPost[,1]
          all_post_shape1=ss$s_ns_tmp[[dur]]$zPost[,2]
          all_post_scale0=ss$s_ns_tmp[[dur]]$zPost[,3]
          all_post_scale1=ss$s_ns_tmp[[dur]]$zPost[,4]
          # shape:
          all_change_shape=all_post_shape1*10/(all_post_shape0+all_post_shape1*halfperiod)*100
          unc_shape_6h[pix]=quantile(all_change_shape,c(0.95))[[1]] - quantile(all_change_shape,c(0.05))[[1]]
          # scale:
          all_change_scale=all_post_scale1*10/(all_post_scale0+all_post_scale1*halfperiod)*100
          unc_scale_6h[pix]=quantile(all_change_scale,c(0.95))[[1]] - quantile(all_change_scale,c(0.05))[[1]]
          # n:
          unc_n_6h[pix]=(ss$n_senss[[dur]]$conf.int[[2]]*10 - ss$n_senss[[dur]]$conf.int[[1]]*10)/ 
            (ss$n_lmfit[[dur]]$coefficients[[1]] + ss$n_lmfit[[dur]]$coefficients[[2]]*halfperiod)*100
          # quantile 2 years:
          # compute smev parameters at first, last and 10 year (all mcmc samaple):
          shape_t_10=all_post_shape0+all_post_shape1*10   # decade
          shape_t_half=all_post_shape0+all_post_shape1*halfperiod # last year
          shape_t_last=all_post_shape0+all_post_shape1*tail(years_smev,1) # last year
          scale_t_10=all_post_scale0+all_post_scale1*10   # decade
          scale_t_half=all_post_scale0+all_post_scale1*halfperiod
          scale_t_last=all_post_scale0+all_post_scale1*tail(years_smev,1) # last year
          n_0=ss$n_lmfit[[dur]]$coefficients[[1]]
          n_10=ss$n_lmfit[[dur]]$coefficients[[1]]+ss$n_lmfit[[dur]]$coefficients[[2]]*10
          n_half=ss$n_lmfit[[dur]]$coefficients[[1]]+ss$n_lmfit[[dur]]$coefficients[[2]]*halfperiod
          n_last=ss$n_lmfit[[dur]]$coefficients[[1]]+ss$n_lmfit[[dur]]$coefficients[[2]]*tail(years_smev,1)
          pr=1-1/2   # TR= 2 years
          all_qnt_2yr_0 = sapply(scale  = all_post_scale0,
                                 shape  = all_post_shape0,
                                 n_ordev= n_0,       # n_ordev,
                                 FUN    = weibull_inv,
                                 X      = pr)
          all_qnt_2yr_half = sapply(scale  = scale_t_half,
                                    shape  = shape_t_half,
                                    n_ordev= n_half,      # n_ordev,
                                    FUN    = weibull_inv,
                                    X      = pr)
          all_qnt_2yr_10 = sapply(scale  = scale_t_10,
                                  shape  = shape_t_10,
                                  n_ordev= n_10,      # n_ordev,
                                  FUN    = weibull_inv,
                                  X      = pr)
          all_qnt_2yr_last = sapply(scale  = scale_t_last,
                                    shape  = shape_t_last,
                                    n_ordev= n_last,  # n_ordev,
                                    FUN    = weibull_inv,
                                    X      = pr)
          # qnt changes for 2yrs return period:
          all_change_qnt2yr=(all_qnt_2yr_last-all_qnt_2yr_0)/
            (all_qnt_2yr_10)*100/n_decades
          # uncertainty at 90%:
          unc_qnt2yrs_6h[pix]=quantile(all_change_qnt2yr,c(0.95))[[1]]-
            quantile(all_change_qnt2yr,c(0.05))[[1]]

          #^********************************************************************
          #                             24h:
          #^********************************************************************
          dur=4
          all_post_shape0=ss$s_ns_tmp[[dur]]$zPost[,1]
          all_post_shape1=ss$s_ns_tmp[[dur]]$zPost[,2]
          all_post_scale0=ss$s_ns_tmp[[dur]]$zPost[,3]
          all_post_scale1=ss$s_ns_tmp[[dur]]$zPost[,4]
          # shape:
          all_change_shape=all_post_shape1*10/
            (all_post_shape0+all_post_shape1*halfperiod)*100
          unc_shape_24h[pix]=quantile(all_change_shape,c(0.95))[[1]]-
            quantile(all_change_shape,c(0.05))[[1]]
          # scale:
          all_change_scale=all_post_scale1*10/ 
            (all_post_scale0+all_post_scale1*halfperiod)*100
          unc_scale_24h[pix]=quantile(all_change_scale,c(0.95))[[1]]-
            quantile(all_change_scale,c(0.05))[[1]]
          # n:
          unc_n_24h[pix]=(ss$n_senss[[dur]]$conf.int[[2]]*10-ss$n_senss[[dur]]$conf.int[[1]]*10)/ 
            (ss$n_lmfit[[dur]]$coefficients[[1]]+ss$n_lmfit[[dur]]$coefficients[[2]]*halfperiod)*100
          # quantile 2 years:
          # compute smev parameters at first, last and 10 year (all mcmc samaple):
          shape_t_10  =all_post_shape0+all_post_shape1*10   # decade
          shape_t_half=all_post_shape0+all_post_shape1*halfperiod # last year
          shape_t_last=all_post_shape0+all_post_shape1*tail(years_smev,1) # last year
          scale_t_10  =all_post_scale0+all_post_scale1*10   # decade
          scale_t_half=all_post_scale0+all_post_scale1*halfperiod
          scale_t_last=all_post_scale0+all_post_scale1*tail(years_smev,1) # last year
          n_0   =ss$n_lmfit[[dur]]$coefficients[[1]]
          n_10  =ss$n_lmfit[[dur]]$coefficients[[1]]+ss$n_lmfit[[dur]]$coefficients[[2]]*10
          n_half=ss$n_lmfit[[dur]]$coefficients[[1]]+ss$n_lmfit[[dur]]$coefficients[[2]]*halfperiod
          n_last=ss$n_lmfit[[dur]]$coefficients[[1]]+ss$n_lmfit[[dur]]$coefficients[[2]]*tail(years_smev,1)
          pr=1-1/2 # TR=2 years
          all_qnt_2yr_0=sapply(scale=all_post_scale0,
                               shape=all_post_shape0,
                               n_ordev=n_0,       # n_ordev,
                               FUN=weibull_inv,
                               X=pr)
          all_qnt_2yr_half = sapply(scale  =scale_t_half,
                                    shape  =shape_t_half,
                                    n_ordev=n_half,      # n_ordev,
                                    FUN    =weibull_inv,
                                    X      =pr)
          all_qnt_2yr_10 = sapply(scale  =scale_t_10,
                                  shape  =shape_t_10,
                                  n_ordev=n_10,      # n_ordev,
                                  FUN    =weibull_inv,
                                  X      =pr)
          all_qnt_2yr_last = sapply(scale  =scale_t_last,
                                    shape  =shape_t_last,
                                    n_ordev=n_last,  # n_ordev,
                                    FUN    =weibull_inv,
                                    X      =pr)
          # changes in 2yrs qnt per decade:
          all_change_qnt2yr=(all_qnt_2yr_last-all_qnt_2yr_0)/
            (all_qnt_2yr_10)*100/n_decades 
          # uncertainty at 90% of these changes:
          unc_qnt2yrs_24h[pix]=quantile(all_change_qnt2yr,c(0.95))[[1]]-
            quantile(all_change_qnt2yr,c(0.05))[[1]]
          
          rm(ss) # remove temporary ss object
        }
      }
    }
  }
  # save all results into new .RData file:
  base::save.image(dir_output_RData)
  
} else {
  #^****************************************************************************
  # load previous results from existing .RData file:
  #^****************************************************************************
  #^ Be careful with this: the upload of RData file includes the loading also 
  # of all settings, directories and R functions saved in the rdata file, 
  # which might be different from the current settings. 
  # In case, recompile the settings and the current module functions.
  base::load(dir_input_RData)
}



#^******************************************************************************
# Refine the directories and recompile here the R functions.
#^******************************************************************************
dir_code<-getwd()
setwd(dir_code)
setwd('../')
# define the directory with the modules files:
dir.modules=paste0(getwd(),"/modules")
# Include all modules for computation:
sapply(paste0(dir.modules,"/",list.files(paste0(getwd(),"/modules"))),source,chdir=F)
setwd(dir_code)


#^******************************************************************************
# reload the settings
#^******************************************************************************
case_study='Italy' # case study name (the same name of the case study folder !!!): 
sat_prod='IMERG_2001_2024' # name of results subfolder for the dataset
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
dir_input='C:/Users/MATTEO/Documents/PRIN/Data/Satellite_data/1_hour/' # folder path of satellite data
filename_input='IMERG_italy_hourly_2001_2024_aggreg_0.20deg.rds' # filename of satellite data
# dir_results=paste0('C:/Users/MATTEO/Documents/PRIN/Results/BaySMEV/',case_study,'/',sat_prod)  # results directory
dir_results=paste0('E:/save_lenovo_unipd_20250711/Documents/PRIN/Results/BaySMEV/',case_study,'/',sat_prod)  # SMEV results directory (where .rds data files are located)

sat_acron='IMERG' # acronyme of satellite products for plot titles
dir_input_RData='E:/save_lenovo_unipd_20250711/Documents/PRIN/Scripts/smev_bayesian/Case_studies/Italy/Results/IMERG_2001_2024/thr_0.9/Results_IMERG_italy_1hr_2001_2024_smevns_NEW.RData'
dir_output_RData='E:/save_lenovo_unipd_20250711/Documents/PRIN/Scripts/smev_bayesian/Case_studies/Italy/Results/IMERG_2001_2024/thr_0.9/Results_IMERG_italy_1hr_2001_2024_smevns_NEW.RData'
dir.res_maps=paste0(dir_results,"/thr_",thr_leftcens,"/maps/")
dir_shp_mask='C:/Users/MATTEO/Documents/PRIN/Qgis/georef-italy-regione/georef-italy-regione-millesime.shp' # shp with admin limits
dir_shp_world='C:/Users/MATTEO/Documents/PRIN/Qgis/WorldGeog/WorldGeog.shp'
dir_shp_user='E:/save_lenovo_unipd_20250711/Documents/PRIN/Qgis/geog_areas_paper2.shp' # some specific polygons for statistics
flag_read_only=T








################################################################################
#                     ORGANIZE THE RESULTS INTO GRID ARRAYS:
################################################################################
# loop on rainfall durations (in this case 4 durations: 1, 3, 6, 24 h):
# not tested in case of different durations, probably not working right now,
# to be improved!

for (dur in 1:length(durations)){
  
  #^****************************************************************************
  # Initialization for each duration:
  #^****************************************************************************
  # CI posterior:
  a_w_CI_ns[[dur]]=b_w_CI_ns[[dur]]=a_C_CI_ns[[dur]]=
    b_C_CI_ns[[dur]]=rep(NA, nrow(tmp_array)*ncol(tmp_array))
  a_w_CI_partns[[dur]]=b_w_CI_partns[[dur]]=a_C_CI_partns[[dur]]=
    b_C_CI_partns[[dur]]=rep(NA, nrow(tmp_array)*ncol(tmp_array))
  a_w_CI_stat[[dur]]=b_w_CI_stat[[dur]]=a_C_CI_stat[[dur]]=b_C_CI_stat[[dur]]=
    rep(NA,nrow(tmp_array)*ncol(tmp_array))
  a_w_CI_gev[[dur]]=a_L_CI_gev[[dur]]=a_C_CI_gev[[dur]]=
    rep(NA, nrow(tmp_array)*ncol(tmp_array))
  # % change in AMX per decade
  AMAX_change[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array))
  # % change in "n" per decade
  n_change[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array))
  # criteria
  signif_ns_DIC_diff3[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array))
  signif_ns_DIC_diff1[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array))
  signif_ns_BIC[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array))
  # likelihood ratio tests:
  pval_yy_wrt_nn[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array))
  pval_yy_wrt_nn_fmins[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array))
  signif_fmins_ns_LRT[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array)) 
  signif_bayes_ns_LRT[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array))
  # quantiles:
  quant_smev_ns_2yr_change[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array))
  quant_smev_ns_10yr_change[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array))
  quant_smev_ns_20yr_change[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array))
  quant_smev_ns_50yr_change[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array))
  quant_smev_ns_2yr_fmins_change[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array))
  quant_smev_ns_2yr_median_change[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array))
  # trend AMAX:
  signif_trend_AMAX[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array))
  sensslope_AMAX[[dur]]=rep(NA,nrow(tmp_array)*ncol(tmp_array))
  
  #^****************************************************************************
  # Organize SMEV and GEV parameters:
  #^****************************************************************************
  a_w_MAP_ns[[dur]]=eval(parse(text=paste0('a_w_MAP_ns_',durations[dur]/60,'h')))
  b_w_MAP_ns[[dur]]=eval(parse(text=paste0('b_w_MAP_ns_',durations[dur]/60,'h')))
  a_C_MAP_ns[[dur]]=eval(parse(text=paste0('a_C_MAP_ns_',durations[dur]/60,'h')))
  b_C_MAP_ns[[dur]]=eval(parse(text=paste0('b_C_MAP_ns_',durations[dur]/60,'h')))
  a_w_med_ns[[dur]]=eval(parse(text=paste0('a_w_quant50_Post_',durations[dur]/60,'h')))
  b_w_med_ns[[dur]]=eval(parse(text=paste0('b_w_quant50_Post_',durations[dur]/60,'h')))
  a_C_med_ns[[dur]]=eval(parse(text=paste0('a_C_quant50_Post_',durations[dur]/60,'h')))
  b_C_med_ns[[dur]]=eval(parse(text=paste0('b_C_quant50_Post_',durations[dur]/60,'h')))
  a_w_fmins_ns[[dur]]=eval(parse(text=paste0('a_w_fmins_ns_',durations[dur]/60,'h')))
  b_w_fmins_ns[[dur]]=eval(parse(text=paste0('b_w_fmins_ns_',durations[dur]/60,'h')))
  a_C_fmins_ns[[dur]]=eval(parse(text=paste0('a_C_fmins_ns_',durations[dur]/60,'h')))
  b_C_fmins_ns[[dur]]=eval(parse(text=paste0('b_C_fmins_ns_',durations[dur]/60,'h')))
  a_w_MAP_stat[[dur]]=eval(parse(text=paste0('a_w_MAP_stat_',durations[dur]/60,'h')))
  a_C_MAP_stat[[dur]]=eval(parse(text=paste0('a_C_MAP_stat_',durations[dur]/60, 'h')))
  a_w_med_stat[[dur]]=eval(parse(text=paste0('a_w_stat_quant50_Post_',durations[dur]/60,'h')))
  a_C_med_stat[[dur]]=eval(parse(text=paste0('a_C_stat_quant50_Post_',durations[dur]/60,'h')))
  a_w_fmins_stat[[dur]]=eval(parse(text=paste0('a_w_fmins_stat_',durations[dur]/60,'h')))
  a_C_fmins_stat[[dur]]=eval(parse(text=paste0('a_C_fmins_stat_',durations[dur]/60,'h')))
  a_L_MAP_gev[[dur]]=eval(parse(text=paste0('a_L_MAP_gev_',durations[dur]/60,'h')))
  a_w_MAP_gev[[dur]]=eval(parse(text=paste0('a_w_MAP_gev_',durations[dur]/60,'h')))
  a_C_MAP_gev[[dur]]=eval(parse(text=paste0('a_C_MAP_gev_',durations[dur]/60,'h')))
  a_L_med_gev[[dur]]=eval(parse(text=paste0('a_L_gev_quant50_Post_',durations[dur]/60,'h')))
  a_w_med_gev[[dur]]=eval(parse(text=paste0('a_w_gev_quant50_Post_',durations[dur]/60,'h')))
  a_C_med_gev[[dur]]=eval(parse(text=paste0('a_C_gev_quant50_Post_',durations[dur]/60,'h')))
  
  #^****************************************************************************
  #  Compute changes, Uncertainty, quantiles and trends for each pixel
  #^****************************************************************************
  # cycle on pixels:
  for (jj in 1:(nrow(tmp_array)*ncol(tmp_array))) {
    
    #^**************************************************************************
    #                      SIGNIFICANCE OF TRENDS IN AMAX:
    #^**************************************************************************
    signif_trend_AMAX[[dur]][jj]=eval(parse(text=paste0('signif_trend_',durations[dur]/60,'h')))[jj]
    sensslope_AMAX[[dur]][jj]=eval(parse(text=paste0('sensslope_AMAX_',durations[dur]/60,'h')))[jj]
    # Compute mean change of AMAX per decade:
    # % of change of AMAX per decade wrw halfperiod year (a sort of average change per decade):
    if (!is.na(eval(parse(text=paste0('signif_trend_',durations[dur]/60,'h')))[jj])){
      # AMAX values for duration dur:
      # b=median of "yi−m*xi" 
      b=median(AMS[[jj]][,dur]-
                 eval(parse(text=paste0('sensslope_AMAX_',durations[dur]/60,'h')))[jj]*years_smev)
      # plot(Snn[[jj]]$years,Snn[[jj]]$AMS[,dur] )
      # lines(Snn[[jj]]$years,b+sensslope_AMAX[[dur]][jj]*Snn[[jj]]$years)
      # sens.slope(Snn[[jj]]$AMS[,dur])
      AMAX_change[[dur]][jj]=eval(parse(text=paste0('sensslope_AMAX_',durations[dur]/60,'h')))[jj]*10/
        (b+eval(parse(text=paste0('sensslope_AMAX_',durations[dur]/60,'h')))[jj]*halfperiod)*100
    }
    
    #^**************************************************************************
    #                   mean change of "n" per decade:
    #^**************************************************************************
    # % of change of n per decade wrw halfperiod year (a sort of average change per decade):
    if (!is.na(signif_trend_n[jj])){  # It is the same for each duration !!!
      
      # n_change[[dur]][jj] = sensslope_n[[dur]][jj]*10 / (
      #            Snn[[jj]]$lmfit[[dur]]$coefficients[[1]] + 
      #              Snn[[jj]]$lmfit[[dur]]$coefficients[[2]]*Snn[[jj]]$halfperiod) * 100
      # 
      n_change[[dur]][jj]=lmfit_n_slope[jj]*10/
        (lmfit_n_inter[jj]+lmfit_n_slope[jj]*halfperiod)*100
    }
    
    #^**************************************************************************
    #                  SIGNIFICANCE OF NONSTATIONARY SMEV:
    #^**************************************************************************
    # Difference in DIC values ns vs stat. (min difference 3): 
    if (!is.na(eval(parse(text=paste0('DIC3_ns_',durations[dur]/60,'h')))[jj])){
      if (eval(parse(text=paste0('DIC3_ns_',durations[dur]/60,'h')))[jj]-
          eval(parse(text=paste0('DIC3_stat_',durations[dur]/60,'h')))[jj]< -3){  
        # Spiegelhalter et al. 2002 
        # (https://rss.onlinelibrary.wiley.com/doi/pdf/10.1111/1467-9868.00353)
        # Burnham and Anderson (1998) suggested models receiving AIC within 1–2 
        # of the ‘best’ deserve consideration,and 3–7 have considerably less support: 
        # these rules of thumb appear to work reasonably well for DIC. 
        signif_ns_DIC_diff3[[dur]][jj]=T # logical T/F
      } else {
        signif_ns_DIC_diff3[[dur]][jj]=F
      }
    }
    
    # Difference in DIC values ns vs stat. (min difference 1): 
    if (!is.na(eval(parse(text=paste0('DIC3_ns_',durations[dur]/60,'h')))[jj])){
      if (eval(parse(text=paste0('DIC3_ns_',durations[dur]/60,'h')))[jj]-
          eval(parse(text=paste0('DIC3_stat_',durations[dur]/60,'h')))[jj]< -1){  
        # Spiegelhalter et al. 2002 
        # (https://rss.onlinelibrary.wiley.com/doi/pdf/10.1111/1467-9868.00353)
        # Burnham and Anderson (1998) suggested models receiving AIC within 1–2 
        # of the ‘best’ deserve consideration,and 3–7 have considerably less support: 
        # these rules of thumb appear to work reasonably well for DIC. 
        signif_ns_DIC_diff1[[dur]][jj]=T # logical T/F
      } else {
        signif_ns_DIC_diff1[[dur]][jj]=F
      }
    }
    
    # Significance from the BIC values: 
    if (!is.na(eval(parse(text=paste0('BIC_ns_',durations[dur]/60,'h')))[jj])){
      if (eval(parse(text=paste0('BIC_ns_',durations[dur]/60,'h')))[jj] <
          eval(parse(text=paste0('BIC_stat_',durations[dur]/60,'h')))[jj]){
        signif_ns_BIC[[dur]][jj]=T # logical T/F
      } else {
        signif_ns_BIC[[dur]][jj]=F
      }
    }
    
    # # Significance from the Bayes Factor values: 
    # if (!is.null(Snn[[jj]]$s_ns_tmp[[dur]]$BIC) & (!is.null(Snn[[jj]]$s_stat_tmp[[dur]]$BIC))){
    #   BF_yy_wrt_nn[[dur]][jj] = exp((Snn[[jj]]$s_stat_tmp[[dur]]$BIC - Snn[[jj]]$s_ns_tmp[[dur]]$BIC)/2)
    #   
    #   if (BF_yy_wrt_nn[[dur]][jj] > 1){  # weak evidence of nonstat vs stat
    #     signif_ns_BF[[dur]][jj] = T # logical T/F
    #   } else {
    #     signif_ns_BF[[dur]][jj] = F
    #   }
    # }
    
    # Likelihood ratio test (for MLE method):
    # full ns model with respect to stationary model:
    # using stan algorithm (it is NULL if flag_like=F):
    if (!is.na(eval(parse(text=paste0(
      'loglik_bayes_stat_',durations[dur]/60,'h')))[jj])){
      pval_yy_wrt_nn[[dur]][jj]=likelihood_ratio_test(
        loglik_H0=eval(parse(text=paste0(
          'loglik_bayes_stat_',durations[dur]/60,'h')))[jj],
        loglik_H1=eval(parse(text=paste0(
          'loglik_bayes_ns_',durations[dur]/60,'h')))[jj],
        degrees_of_freedom=2)
      if (pval_yy_wrt_nn[[dur]][jj] <= 0.05){
        signif_bayes_ns_LRT[[dur]][jj]=T
      } else {
        signif_bayes_ns_LRT[[dur]][jj]=F
      }
    }
    
    if (!is.na(eval(parse(text=paste0(
      'loglik_fmins_stat_',durations[dur]/60,'h')))[jj])){
      # Using fminsearch function results:
      pval_yy_wrt_nn_fmins[[dur]][jj]=likelihood_ratio_test(
        loglik_H0= -eval(parse(text=paste0(
        'loglik_fmins_stat_',durations[dur]/60,'h')))[jj],
        loglik_H1= -eval(parse(text=paste0(
          'loglik_fmins_ns_',durations[dur]/60,'h')))[jj],
        degrees_of_freedom=2)
      if (pval_yy_wrt_nn_fmins[[dur]][jj] <= 0.05){
        signif_fmins_ns_LRT[[dur]][jj]=T
      } else {
        signif_fmins_ns_LRT[[dur]][jj]=F
      }
    }
    
    #^**************************************************************************
    #                        CREDIBILITY INTERVALS:
    #^**************************************************************************
    # Uncertainty on SMEV NonStationary parameters (4 params):
    if (!is.na(eval(parse(text=paste0('a_w_quant5_Post_',
                                      durations[dur]/60,'h')))[jj])){
      # Credibility Interval at 90%:
      a_w_CI_ns[[dur]][jj]=(eval(parse(text=paste0('a_w_quant95_Post_',
                                                   durations[dur]/60,'h')))[jj]-
                            eval(parse(text=paste0('a_w_quant5_Post_',
                                                   durations[dur]/60,'h')))[jj]) 
      # /eval(parse(text = paste0('a_w_MAP_ns_', durations[dur]/60,'h')))[jj]  # CI shape intercept
      
      b_w_CI_ns[[dur]][jj]=(eval(parse(text=paste0('b_w_quant95_Post_',
                                                   durations[dur]/60, 'h')))[jj]-
                            eval(parse(text=paste0('b_w_quant5_Post_',
                                                   durations[dur]/60, 'h')))[jj])
      # /eval(parse(text=paste0('b_w_MAP_ns_',durations[dur]/60,'h')))[jj]   # CI shape slope
      
      a_C_CI_ns[[dur]][jj]=(eval(parse(text=paste0('a_C_quant95_Post_',
                                                   durations[dur]/60,'h')))[jj] -
                            eval(parse(text=paste0('a_C_quant5_Post_',
                                                   durations[dur]/60,'h')))[jj])/
        eval(parse(text=paste0('a_C_MAP_ns_',durations[dur]/60,'h')))[jj]    # CI scale intercept
      
      b_C_CI_ns[[dur]][jj]=(eval(parse(text=paste0('b_C_quant95_Post_',
                                                   durations[dur]/60,'h')))[jj]-
                            eval(parse(text=paste0('b_C_quant5_Post_',
                                                   durations[dur]/60,'h')))[jj])
      # /eval(parse(text=paste0('b_C_MAP_ns_',durations[dur]/60,'h')))[jj]   # CI scale slope
    }
    
    
    # Uncertainty on SMEV Partial-Non-Stationary parameters (3 params):
    if (exists(paste0('a_w_partns_quant95_Post_',durations[dur]/60,'h'))){
      # Credibility Interval at 90% (posterior):
      a_w_CI_partns[[dur]][jj]=(eval(parse(text=paste0('a_w_partns_quant95_Post_',
                                                       durations[dur]/60,'h')))[jj]-
                                eval(parse(text=paste0('a_w_partns_quant5_Post_',
                                                       durations[dur]/60,'h')))[jj])  
      #/eval(parse(text=paste0('a_w_MAP_partns_',durations[dur]/60,'h')))[jj]   # CI shape intercept
      
      a_C_CI_partns[[dur]][jj]=(eval(parse(text=paste0('a_C_partns_quant95_Post_',
                                                       durations[dur]/60,'h')))[jj] -
                                eval(parse(text=paste0('a_C_partns_quant5_Post_',
                                                       durations[dur]/60,'h')))[jj])/
        eval(parse(text=paste0('a_C_MAP_partns_',durations[dur]/60,'h')))[jj]   # CI scale intercept
      
      b_C_CI_partns[[dur]][jj]=(eval(parse(text=paste0('b_C_partns_quant95_Post_',
                                                       durations[dur]/60,'h')))[jj] -
                                eval(parse(text=paste0('b_C_partns_quant5_Post_', 
                                                       durations[dur]/60,'h')))[jj])/
        eval(parse(text=paste0('b_C_MAP_partns_',durations[dur]/60,'h')))[jj]   # CI scale slope
    }
    
    
    # Uncertainty on SMEV Stationary parameters (2 params):
    if (!is.na(eval(parse(text=paste0('a_w_stat_quant95_Post_',
                                      durations[dur]/60,'h')))[jj])){
      # Credibility Interval at 90% (posterior):
      a_w_CI_stat[[dur]][jj]=(eval(parse(text=paste0('a_w_stat_quant95_Post_',
                                                     durations[dur]/60,'h')))[jj]-
                              eval(parse(text=paste0('a_w_stat_quant5_Post_',
                                                     durations[dur]/60,'h')))[jj])  
      #/eval(parse(text=paste0('a_w_MAP_stat_',durations[dur]/60,'h')))[jj]     # CI shape intercept
      
      a_C_CI_stat[[dur]][jj]=(eval(parse(text=paste0('a_C_stat_quant95_Post_',
                                                     durations[dur]/60,'h')))[jj]-
                              eval(parse(text=paste0('a_C_stat_quant5_Post_',
                                                     durations[dur]/60,'h')))[jj])/
        eval(parse(text=paste0('a_C_MAP_stat_',durations[dur]/60,'h')))[jj]     # CI scale intercept
    }
    
    
    # Uncertainty on GEV Stationary parameters (2 params):
    if (!is.na(eval(parse(text=paste0('a_w_gev_quant95_Post_',
                                      durations[dur]/60,'h')))[jj])){
      # Credibility Interval at 90% (posterior):
      a_L_CI_gev[[dur]][jj]=(eval(parse(text=paste0('a_L_gev_quant95_Post_',
                                                    durations[dur]/60,'h')))[jj]-
                             eval(parse(text=paste0('a_L_gev_quant5_Post_',
                                                    durations[dur]/60,'h')))[jj])/
        eval(parse(text=paste0('a_L_MAP_gev_',durations[dur]/60,'h')))[jj]      # CI loc intercept
      
      a_w_CI_gev[[dur]][jj]=(eval(parse(text=paste0('a_w_gev_quant95_Post_',
                                                    durations[dur]/60,'h')))[jj] -
                             eval(parse(text=paste0('a_w_gev_quant5_Post_',
                                                    durations[dur]/60,'h')))[jj])
      # /eval(parse(text = paste0('a_w_MAP_gev_', durations[dur]/60, 'h')))[jj]  # CI shape intercept
      
      a_C_CI_gev[[dur]][jj]=(eval(parse(text=paste0('a_C_gev_quant95_Post_',
                                                    durations[dur]/60,'h')))[jj] -
                             eval(parse(text=paste0('a_C_gev_quant5_Post_',
                                                    durations[dur]/60,'h')))[jj])/
        eval(parse(text = paste0('a_C_MAP_gev_',durations[dur]/60, 'h')))[jj]   # CI scale intercept
    }
    
    
    #^************************************************************************** 
    # QUANTILES:
    #^**************************************************************************
    #if (!is.na(Snn[[jj]]$quantiles_smev[[dur]])){
    
    # compute smev parameters at first, last and 10 year (MAXPOST):
    shape_t_maxpost_0=a_w_MAP_ns[[dur]][jj]    # first year 0
    shape_t_maxpost_10=a_w_MAP_ns[[dur]][jj]+b_w_MAP_ns[[dur]][jj]*10  # decade
    shape_t_maxpost_last=a_w_MAP_ns[[dur]][jj]+b_w_MAP_ns[[dur]][jj]*tail(years_smev,1) # last year
    
    scale_t_maxpost_0=a_C_MAP_ns[[dur]][jj]
    scale_t_maxpost_10=a_C_MAP_ns[[dur]][jj]+b_C_MAP_ns[[dur]][jj]*10 #decade
    scale_t_maxpost_last=a_C_MAP_ns[[dur]][jj]+b_C_MAP_ns[[dur]][jj]*tail(years_smev,1) # last year
    
    # compute smev parameters at first, last and 10 year (MAXLIKELIHOOD fmins):
    shape_t_maxlikfmins_0    =a_w_fmins_ns[[dur]][jj]       # first year 0
    shape_t_maxlikfmins_10   =a_w_fmins_ns[[dur]][jj]+b_w_fmins_ns[[dur]][jj]*10  # decade
    shape_t_maxlikfmins_half =a_w_fmins_ns[[dur]][jj]+b_w_fmins_ns[[dur]][jj]*halfperiod  # decade
    shape_t_maxlikfmins_last =a_w_fmins_ns[[dur]][jj]+b_w_fmins_ns[[dur]][jj]*tail(years_smev,1) # last year
    scale_t_maxlikfmins_0    =a_C_fmins_ns[[dur]][jj]       # first year 0
    scale_t_maxlikfmins_10   =a_C_fmins_ns[[dur]][jj]+b_C_fmins_ns[[dur]][jj]*10 #decade
    scale_t_maxlikfmins_half =a_C_fmins_ns[[dur]][jj]+b_C_fmins_ns[[dur]][jj]*halfperiod  # decade
    scale_t_maxlikfmins_last =a_C_fmins_ns[[dur]][jj]+b_C_fmins_ns[[dur]][jj]*tail(years_smev,1) # last year
    
    # compute smev parameters at first, last and 10 year (MEDIAN POSTERIOR):
    shape_t_median_0    =a_w_med_ns[[dur]][jj]       # first year 0
    shape_t_median_10   =a_w_med_ns[[dur]][jj]+b_w_med_ns[[dur]][jj]*10  # decade
    shape_t_median_half =a_w_med_ns[[dur]][jj]+b_w_med_ns[[dur]][jj]*halfperiod  # decade
    shape_t_median_last =a_w_med_ns[[dur]][jj]+b_w_med_ns[[dur]][jj]*tail(years_smev,1) # last year
    scale_t_median_0    =a_C_med_ns[[dur]][jj]       # first year 0
    scale_t_median_10   =a_C_med_ns[[dur]][jj]+b_C_med_ns[[dur]][jj]*10 #decade
    scale_t_median_half =a_C_med_ns[[dur]][jj]+b_C_med_ns[[dur]][jj]*halfperiod  # decade
    scale_t_median_last =a_C_med_ns[[dur]][jj]+b_C_med_ns[[dur]][jj]*tail(years_smev,1) # last year
    
    n_0=lmfit_n_inter[jj] # n for the first year
    n_10=lmfit_n_inter[jj]+lmfit_n_slope[jj]*10 # n after 10 years
    n_half=lmfit_n_inter[jj]+lmfit_n_slope[jj]*halfperiod # n at the halfperiod year
    n_last=lmfit_n_inter[jj]+lmfit_n_slope[jj]*tail(years_smev,1) # n for the last year
    
    
    # Get the number of decades within the period:
    n_decades=tail(years_smev,1)/10
    
    #^**************************************************************************
    # 2 YEARS RETURN PERIOD:
    #^**************************************************************************
    pr=1-1/2 # prob
    quant_smev_ns_2yr_0=sapply(scale  =scale_t_maxpost_0,
                               shape  =shape_t_maxpost_0,
                               n_ordev=n_0,     # n_ordev,
                               FUN    =weibull_inv,
                               X      =pr)
    
    quant_smev_ns_2yr_10=sapply(scale  =scale_t_maxpost_10,
                                shape  =shape_t_maxpost_10,
                                n_ordev=n_10,   # n_ordev,
                                FUN    =weibull_inv,
                                X      =pr)
    
    quant_smev_ns_2yr_last=sapply(scale  =scale_t_maxpost_last,
                                  shape  =shape_t_maxpost_last,
                                  n_ordev=n_last,  # n_ordev,
                                  FUN    =weibull_inv,
                                  X      =pr)
    
    # fminsearch:
    quant_smev_ns_2yr_fmins_0=NA
    quant_smev_ns_2yr_fmins_0=sapply(scale  =scale_t_maxlikfmins_0,
                                     shape  =shape_t_maxlikfmins_0,
                                     n_ordev=n_0,   # n_ordev,
                                     FUN    =weibull_inv,
                                     X      =pr)
    quant_smev_ns_2yr_fmins_10=sapply(scale  =scale_t_maxlikfmins_10,
                                      shape  =shape_t_maxlikfmins_10,
                                      n_ordev=n_10,  # n_ordev,
                                      FUN    =weibull_inv,
                                      X      =pr)
    # ss$quantiles_smev[[1]]$qCurve_maxlik_fmins_ns
    quant_smev_ns_2yr_fmins_half=NA
    quant_smev_ns_2yr_fmins_half=sapply(scale  =scale_t_maxlikfmins_half,
                                        shape  =shape_t_maxlikfmins_half,
                                        n_ordev=n_half,   # n_ordev,
                                        FUN    =weibull_inv,
                                        X      =pr)
    quant_smev_ns_2yr_fmins_last=NA
    quant_smev_ns_2yr_fmins_last=sapply(scale  =scale_t_maxlikfmins_last,
                                        shape  =shape_t_maxlikfmins_last,
                                        n_ordev=n_last,  # n_ordev,
                                        FUN    =weibull_inv,
                                        X      =pr)
    
    # MEDIAN post:
    quant_smev_ns_2yr_median_0=NA
    quant_smev_ns_2yr_median_0=sapply(scale  =scale_t_median_0,
                                      shape  =shape_t_median_0,
                                      n_ordev=n_0,  # n_ordev,
                                      FUN    =weibull_inv,
                                      X      =pr)
    quant_smev_ns_2yr_median_10=sapply(scale  =scale_t_median_10,
                                       shape  =shape_t_median_10,
                                       n_ordev=n_10,  # n_ordev,
                                       FUN    =weibull_inv,
                                       X      =pr)
    # ss$quantiles_smev[[1]]$qCurve_maxlik_fmins_ns
    quant_smev_ns_2yr_median_half=NA
    quant_smev_ns_2yr_median_half=sapply(scale  =scale_t_median_half,
                                         shape  =shape_t_median_half,
                                         n_ordev=n_half,  # n_ordev,
                                         FUN    =weibull_inv,
                                         X      =pr)
    quant_smev_ns_2yr_median_last=NA
    quant_smev_ns_2yr_median_last=sapply(scale  =scale_t_median_last,
                                         shape  =shape_t_median_last,
                                         n_ordev=n_last,  # n_ordev,
                                         FUN    =weibull_inv,
                                         X      =pr)
    
    # quant_smev_ns_2yr_change[[dur]][jj]=(quant_smev_ns_2yr_10-quant_smev_ns_2yr_0)/quant_smev_ns_2yr[[dur]][jj]*100
    quant_smev_ns_2yr_change[[dur]][jj]=(quant_smev_ns_2yr_last-quant_smev_ns_2yr_0)/
      eval(parse(text=paste0('qCurve_2yr_maxpost_',durations[dur]/60,'h')))[jj]*100/n_decades
    quant_smev_ns_2yr_fmins_change[[dur]][jj]=(quant_smev_ns_2yr_fmins_last-quant_smev_ns_2yr_fmins_0)/ 
      quant_smev_ns_2yr_fmins_half*100/n_decades
    quant_smev_ns_2yr_median_change[[dur]][jj]=(quant_smev_ns_2yr_median_last-quant_smev_ns_2yr_median_0)/ 
      quant_smev_ns_2yr_median_half*100/n_decades
    
    #^**************************************************************************
    # 10 YEARS:
    #^**************************************************************************
    pr=1-1/10
    quant_smev_ns_10yr_0=sapply(scale=scale_t_maxpost_0,
                                shape=shape_t_maxpost_0,
                                n_ordev=n_0, # n_ordev,
                                FUN=weibull_inv,
                                X=pr)
    quant_smev_ns_10yr_10=sapply(scale=scale_t_maxpost_10,
                                 shape=shape_t_maxpost_10,
                                 n_ordev=n_10, # n_ordev,
                                 FUN=weibull_inv,
                                 X=pr)
    quant_smev_ns_10yr_last=sapply(scale=scale_t_maxpost_last,
                                   shape=shape_t_maxpost_last,
                                   n_ordev=n_last, # n_ordev
                                   FUN=weibull_inv,
                                   X=pr)
    # quant_smev_ns_10yr_change[[dur]][jj]=(quant_smev_ns_10yr_10-quant_smev_ns_10yr_0)/quant_smev_ns_10yr[[dur]][jj]*100
    quant_smev_ns_10yr_change[[dur]][jj]=(quant_smev_ns_10yr_last-quant_smev_ns_10yr_0)/ 
      eval(parse(text=paste0('qCurve_10yr_maxpost_',durations[dur]/60,'h')))[jj]*100/n_decades
    
    #^**************************************************************************
    # 20 YEARS:
    #^**************************************************************************
    pr=1-1/20
    quant_smev_ns_20yr_0=sapply(scale=scale_t_maxpost_0,
                                shape=shape_t_maxpost_0,
                                n_ordev=n_0, # n_ordev,
                                FUN=weibull_inv,
                                X=pr)
    quant_smev_ns_20yr_10=sapply(scale=scale_t_maxpost_10,
                                 shape=shape_t_maxpost_10,
                                 n_ordev=n_10, # n_ordev,
                                 FUN=weibull_inv,
                                 X=pr)
    quant_smev_ns_20yr_last=sapply(scale=scale_t_maxpost_last,
                                   shape=shape_t_maxpost_last,
                                   n_ordev=n_last, # n_ordev
                                   FUN=weibull_inv,
                                   X=pr)
    quant_smev_ns_20yr_change[[dur]][jj]=(quant_smev_ns_20yr_last-quant_smev_ns_20yr_0)/ 
      eval(parse(text=paste0('qCurve_20yr_maxpost_',durations[dur]/60,'h')))[jj]*100/n_decades
    
    #^**************************************************************************
    # 50 YEARS:
    #^**************************************************************************
    pr=1-1/50
    quant_smev_ns_50yr_0=sapply(scale=scale_t_maxpost_0,
                                shape=shape_t_maxpost_0,
                                n_ordev=n_last, # n_ordev,
                                FUN=weibull_inv,
                                X=pr)
    quant_smev_ns_50yr_10=sapply(scale=scale_t_maxpost_10,
                                 shape=shape_t_maxpost_10,
                                 n_ordev=n_last, # n_ordev,
                                 FUN=weibull_inv,
                                 X=pr)
    quant_smev_ns_50yr_last=sapply(scale=scale_t_maxpost_last,
                                   shape=shape_t_maxpost_last,
                                   n_ordev=n_last, # n_ordev,
                                   FUN=weibull_inv,
                                   X=pr)
    # quant_smev_ns_50yr_change[[dur]][jj]=(quant_smev_ns_50yr_10-quant_smev_ns_50yr_0)/quant_smev_ns_50yr[[dur]][jj]*100
    quant_smev_ns_50yr_change[[dur]][jj]=(quant_smev_ns_50yr_last-quant_smev_ns_50yr_0)/
      eval(parse(text=paste0('qCurve_50yr_maxpost_',durations[dur]/60,'h')))[jj]*100/n_decades
    
    #^**************************************************************************
    # # CI uncertainty of quantiles:
    #^**************************************************************************
    # a_w_quant5_Post_1h + a_w_quant5_Post_1h 
    # a_w_quant95_Post_1h
    # 
    # # compute smev parameters at first, last and 10 year (MAXPOST):
    # shape_t_q5_0=eval(parse(text=paste0('a_w_quant5_Post_',durations[dur]/60,'h')))[jj]   # first year 0
    # shape_t_q95_10=eval(parse(text=paste0('a_w_quant5_Post_',durations[dur]/60,'h')))[jj]+
    # shape_t_maxpost_last=a_w_MAP_ns[[dur]][jj]+b_w_MAP_ns[[dur]][jj]*tail(years_smev,1) # last year
    # 
    # scale_t_maxpost_0=a_C_MAP_ns[[dur]][jj]
    # scale_t_maxpost_10=a_C_MAP_ns[[dur]][jj]+b_C_MAP_ns[[dur]][jj]*10 #decade
    # scale_t_maxpost_last=a_C_MAP_ns[[dur]][jj]+b_C_MAP_ns[[dur]][jj]*tail(years_smev,1) # last year
    # 
    # quant_smev_ns_2yr_CI90_change[[dur]][jj]=(quant_smev_ns_2yr_CI90_last-quant_smev_ns_2yr_CI90_0)/
    #                                           quant_smev_ns_2yr_CI90_half*100/n_decades
    
    # }
  }
}









#^******************************************************************************
# Get results of significance from trend tests:
#^******************************************************************************
# make a list of tests results 
significance=list(signif_MK_AMAX= signif_trend_AMAX,   # logical T/F
                  sensslope_AMAX= sensslope_AMAX,      # slope of the Sen's test
                  pval_MK_AMAX  = pval_MK_AMAX,        # p-value of the MK test
                  z_MK_AMAX     = z_MK_AMAX,           # z estimate of MK test
                  DIC_diff3     = signif_ns_DIC_diff3, # logical (T if nonstat mod has lower DIC, tol of 3)
                  DIC_diff1     = signif_ns_DIC_diff1, # logical (T if nonstat mod has lower DIC, tol of 1)
                  #BF           = BF                   # bayes factor (computed from BIC)
                  signif_fmins_ns_LRT = signif_fmins_ns_LRT, # logical (T if nonstat vs stat has pvalue<0.05, LRT method)
                  signif_MK_n   = signif_trend_n, # significance of trends in "n"
                  sensslope_n   = n_senss, # sen's slope of trends in "n"
                  pval_MK_n     = pval_MK_n) # MK test of trends in "n"









################################################################################
# PLOT MAPS OF ALL SMEV PARAMETERS:
################################################################################

#^******************************************************************************
# Shape intercept:
#^******************************************************************************
# NONSTATIONARY MODEL:
list_plot_shapeint_ns=plot_maps_param(
          var_list      = list(a_w_MAP_ns,a_w_CI_ns,a_w_fmins_ns),#a_w_med_ns,a_w_MaxLik_ns),
          lon_tmp       = lon_tmp,
          lat_tmp       = lat_tmp,
          durations     = durations,
          mask_shp      = land_shp,
          mask_world    = gg_world,
          flag_signif   = T,# logical
          significance  = significance,
          name_var      = "Shape intercept",
          type_estimate = c("MaxPost","CI 90%","MaxLik fmins"), #"MedianPost","MaxLik"),
          name_model    = "SMEV Nonstationary",
          acronym_model = "ns",
          path_output   = c(paste0(dir.res_maps,"/shape_int/map_shape_int_MaxPost_ns"),
                            paste0(dir.res_maps,"/shape_int/map_shape_int_CI90_ns"),
                            paste0(dir.res_maps,"/shape_int/map_shape_int_MaxLikfmins_ns")),
                            # paste0(dir.res_maps,"/shape_int/map_shape_int_medPost_ns"),
                            # paste0(dir.res_maps,"/shape_int/map_shape_int_MaxLik_ns")),
          limits        = list(c(0.5,1.2,0.1),c(0.1,0.4,0.1),c(0.5,1.2,0.1)), #c(0.4,1.4),c(0.4,1.4)),
          flag.limits   = T,
          scale_grad    = list(list(name="turbo"),
                               list(name="turbo"),
                               list(name="turbo")),
                              # list(name="turbo"),
                              # list(name="turbo")),
          crs_map       = crswgs84,
          title_x       = "Longitude (°E)",
          title_y       = "Latitude (°N)",
          breaks_x      = c(6.00001,10,14,18),
          breaks_y      = c(36.00001,40,44,47.9999),
          labels_x      = c("6°E","10°E","14°E","18°E"),
          labels_y      = c("36°N","40°N","44°N","48°N"),
          limits_x      = c(6,19),
          limits_y      = c(36,48),
          alpha_shpworld= 0.2,
          size_points   = 0.7,
          flag.saveRDS  = F)


# STATIONARY MODEL:
list_plot_shapeint_stat = plot_maps_param(
  var_list      = list(a_w_MAP_stat, a_w_CI_stat, a_w_fmins_stat),#a_w_med_stat, a_w_MaxLik_stat),
  lon_tmp       = lon_tmp,
  lat_tmp       = lat_tmp,
  durations     = durations,
  mask_shp      = land_shp,
  mask_world    = gg_world,
  flag_signif   = F, # logical
  significance  = significance,
  name_var      = "Shape intercept",
  type_estimate = c("MaxPost","CI 90%","MaxLik fmins"), #"MedianPost","MaxLik"),
  name_model    = "SMEV Stationary",
  acronym_model = "stat",
  path_output   = c(paste0(dir.res_maps,"/shape_int/map_shape_int_MaxPost_stat"),
                    paste0(dir.res_maps,"/shape_int/map_shape_int_CI90_stat"),
                    paste0(dir.res_maps,"/shape_int/map_shape_int_MaxLikfmins_stat")),
                    # paste0(dir.res_maps,"/shape_int/map_shape_int_medPost_stat"),
                    # paste0(dir.res_maps,"/shape_int/map_shape_int_MaxLik_stat")),
  limits        = list(c(0.5,1.2,0.1),c(0.1,0.4,0.1),c(0.5,1.2,0.1)), #c(0.4,1.4),c(0.4,1.4)),
  flag.limits   = T,
  scale_grad    = list(list(name="turbo"),
                       # list(name="rdbl",low="white",high="blue",mid="white",midpoint=0),
                       list(name="turbo"),
                       list(name="turbo")),
                        # list(name="turbo"),
                        # list(name="turbo")),
  crs_map       = crswgs84,
  title_x       = "Longitude (°E)",
  title_y       = "Latitude (°N)",
  breaks_x      = c(6.00001,10,14,18),
  breaks_y      = c(36.00001,40,44,47.9999),
  labels_x      = c("6°E","10°E","14°E","18°E"),
  labels_y      = c("36°N","40°N","44°N","48°N"),
  limits_x      = c(6,19),
  limits_y      = c(36,48),
  alpha_shpworld= 0.2,
  size_points   = 0.7,
  flag.saveRDS  = F)



# GEV STATIONARY MODEL:
list_plot_shapeint_gev=plot_maps_param(
  var_list      = list(a_w_MAP_gev,a_w_CI_gev), # a_w_med_gev),
  lon_tmp       = lon_tmp,
  lat_tmp       = lat_tmp,
  durations     = durations,
  mask_shp      = land_shp,
  mask_world    = gg_world,
  flag_signif   = F, # logical
  significance  = significance,
  name_var      = "Shape intercept",
  type_estimate = c("MaxPost","CI 90%"), #"MedianPost"),
  name_model    = "GEV Stationary",
  acronym_model = "gev",
  path_output   = c(paste0(dir.res_maps,"/shape_int/map_shape_int_MaxPost_gev"),
                    paste0(dir.res_maps,"/shape_int/map_shape_int_CI90_gev")),
                    # paste0(dir.res_maps, "/shape_int/map_shape_int_medPost_gev")),
  limits        = list(c(-0.2,0.5,0.1),c(0.3,0.7,0.1)), #c(0.4,1.4),c(0.4,1.4)),
  flag.limits   = T,
  scale_grad    = list(list(name="turbo"),
                       list(name="turbo")),
                      # list(name="turbo")),
  crs_map       = crswgs84,
  title_x       = "Longitude (°E)",
  title_y       = "Latitude (°N)",
  breaks_x      = c(6.00001,10,14,18),
  breaks_y      = c(36.00001,40,44,47.9999),
  labels_x      = c("6°E","10°E","14°E","18°E"),
  labels_y      = c("36°N","40°N","44°N","48°N"),
  limits_x      = c(6,19),
  limits_y      = c(36,48),
  alpha_shpworld= 0.2,
  size_points   = 0.7,
  flag.saveRDS  = F)



#^******************************************************************************
# Shape slope:
#^******************************************************************************
# NONSTATIONARY MODEL:
list_plot_shapeslope_ns=plot_maps_param(
  var_list      = list(b_w_MAP_ns,b_w_CI_ns,b_w_fmins_ns), # b_w_med_ns,b_w_MaxLik_ns),
  lon_tmp       = lon_tmp,
  lat_tmp       = lat_tmp,
  durations     = durations,
  mask_shp      = land_shp,
  mask_world    = gg_world,
  flag_signif   = T, # logical
  significance  = significance,
  name_var      = "Shape slope",
  type_estimate = c("MaxPost","CI 90%","MaxLik fmins"), #"MedianPost","MaxLik"),
  name_model    = "SMEV Nonstationary",
  acronym_model = "ns",
  path_output   = c(paste0(dir.res_maps,"/shape_slope/map_shape_slope_MaxPost_ns"),
                    paste0(dir.res_maps,"/shape_slope/map_shape_slope_CI90_ns"),
                    paste0(dir.res_maps,"/shape_slope/map_shape_slope_MaxLikfmins_ns")),
                    # paste0(dir.res_maps,"/shape_slope/map_shape_slope_medPost_ns"),
                    # paste0(dir.res_maps,"/shape_slope/map_shape_slope_MaxLik_ns")),
  limits        = list(c(-0.02, 0.02, 0.01),c(0.015,0.025, 0.001),c(-0.02, 0.02, 0.01)), #,c(-0.02,0.02),c(-0.02,0.02)),
  flag.limits   = T,
  scale_grad    = list(list(name="rdbl",low="red",high="blue",mid="white",midpoint=0),
                       list(name="turbo"),
                       list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)),
                        # list(name="rdbl",low="red",high= blue",mid="white",midpoint=0),
                        # list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)),
  crs_map       = crswgs84,
  title_x       = "Longitude (°E)",
  title_y       = "Latitude (°N)",
  breaks_x      = c(6.00001,10,14,18),
  breaks_y      = c(36.00001,40,44,47.9999),
  labels_x      = c("6°E","10°E","14°E","18°E"),
  labels_y      = c("36°N","40°N","44°N","48°N"),
  limits_x      = c(6,19),
  limits_y      = c(36,48),
  alpha_shpworld= 0.2,
  size_points   = 0.7,
  flag.saveRDS  = F)


# no limits:
list_plot_shapeslope_ns=plot_maps_param(
  var_list      = list(b_w_MAP_ns,b_w_CI_ns,b_w_fmins_ns), # b_w_med_ns,b_w_MaxLik_ns),
  lon_tmp       = lon_tmp,
  lat_tmp       = lat_tmp,
  durations     = durations,
  mask_shp      = land_shp,
  mask_world    = gg_world,
  flag_signif   = T, # logical
  significance  = significance,
  name_var      = "Shape slope",
  type_estimate = c("MaxPost","CI 90%","MaxLik fmins"), #"MedianPost","MaxLik"),
  name_model    = "SMEV Nonstationary",
  acronym_model = "ns",
  path_output   = c(paste0(dir.res_maps,"/shape_slope/map_shape_slope_MaxPost_ns_NOLIM"),
                    paste0(dir.res_maps,"/shape_slope/map_shape_slope_CI90_ns_NOLIM"),
                    paste0(dir.res_maps,"/shape_slope/map_shape_slope_MaxLikfmins_ns_NOLIM")),
                    # paste0(dir.res_maps,"/shape_slope/map_shape_slope_medPost_ns"),
                    # paste0(dir.res_maps,"/shape_slope/map_shape_slope_MaxLik_ns")),
  limits        = list(c(-0.02,0.02,0.01),c(0.015,0.025,0.001),c(-0.02,0.02,0.01)), #,c(-0.02,0.02),c(-0.02,0.02)),
  flag.limits   = F,
  scale_grad    = list(list(name="rdbl",low="red",high="blue",mid="white",midpoint=0),
                       list(name="turbo"),
                       list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)),
                        # list(name="rdbl",low="red",high="blue",mid="white",midpoint=0),
                        # list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)),
  crs_map       = crswgs84,
  title_x       = "Longitude (°E)",
  title_y       = "Latitude (°N)",
  breaks_x      = c(6.00001,10,14,18),
  breaks_y      = c(36.00001,40,44,47.9999),
  labels_x      = c("6°E","10°E","14°E","18°E"),
  labels_y      = c("36°N","40°N","44°N","48°N"),
  limits_x      = c(6,19),
  limits_y      = c(36,48),
  alpha_shpworld= 0.2,
  size_points   = 0.7,
  flag.saveRDS  = F)



#^******************************************************************************
# Scale intercept:
#^******************************************************************************
# NONSTATIONARY MODEL:
# with no limits:
list_plot_scaleint_ns = plot_maps_param(
  var_list      = list(a_C_MAP_ns, a_C_CI_ns, a_C_fmins_ns), #a_C_med_ns, , a_C_MaxLik_ns),
  lon_tmp       = lon_tmp,
  lat_tmp       = lat_tmp,
  durations     = durations,
  mask_shp      = land_shp,
  mask_world    = gg_world,
  flag_signif   = T, # logical
  significance  = significance,
  name_var      = "Scale intercept",
  type_estimate = c("MaxPost","CI 90%","MaxLik fmins"), #"MedianPost",, "MaxLik"),
  name_model    = "SMEV Nonstationary",
  acronym_model = "ns",
  path_output   = c(paste0(dir.res_maps,"/scale_int/map_scale_int_MaxPost_ns_NOLIM"),
                    paste0(dir.res_maps,"/scale_int/map_scale_int_CI90_ns_NOLIM"),
                    paste0(dir.res_maps,"/scale_int/map_scale_int_MaxLikfmins_ns_NOLIM")),
                    # paste0(dir.res_maps,"/scale_int/map_scale_int_medPost_ns"),
                    # paste0(dir.res_maps,"/scale_int/map_scale_int_MaxLik_ns")),
  limits        = list(c(0,5,1),c(0.3,0.8,0.1),c(0,5,1)), #,c(0,5),c(0,5)),
  flag.limits   = F,  # no limits for the scale, they change for different durations.
  scale_grad    = list(list(name="turbo"),
                       list(name="turbo"),
                       list(name="turbo")),
                        # list(name="turbo"),
                        # list(name="turbo")),
  crs_map       = crswgs84,
  title_x       = "Longitude (°E)",
  title_y       = "Latitude (°N)",
  breaks_x      = c(6.00001,10,14,18),
  breaks_y      = c(36.00001,40,44,47.9999),
  labels_x      = c("6°E","10°E","14°E","18°E"),
  labels_y      = c("36°N","40°N","44°N","48°N"),
  limits_x      = c(6,19),
  limits_y      = c(36,48),
  alpha_shpworld= 0.2,
  size_points   = 0.7,
  flag.saveRDS  = F)


# STATIONARY MODEL:
list_plot_scaleint_stat = plot_maps_param(
  var_list      = list(a_C_MAP_stat,a_C_CI_stat,a_C_fmins_stat), #a_C_med_stat,a_C_MaxLik_stat),
  lon_tmp       = lon_tmp,
  lat_tmp       = lat_tmp,
  durations     = durations,
  mask_shp      = land_shp,
  mask_world    = gg_world,
  flag_signif   = F, # logical
  significance  = significance,
  name_var      = "Scale intercept",
  type_estimate = c("MaxPost","CI 90%","MaxLik fmins"), #"MedianPost","MaxLik"),
  name_model    = "SMEV Stationary",
  acronym_model = "stat",
  path_output   = c(paste0(dir.res_maps,"/scale_int/map_scale_int_MaxPost_stat_NOLIM"),
                    paste0(dir.res_maps,"/scale_int/map_scale_int_CI90_stat_NOLIM"),
                    paste0(dir.res_maps,"/scale_int/map_scale_int_MaxLikfmins_stat_NOLIM")),
                    # paste0(dir.res_maps,"/scale_int/map_scale_int_medPost_stat"),
                    # paste0(dir.res_maps,"/scale_int/map_scale_int_MaxLik_stat")),
  limits        = list(c(0,5,1),c(0.01,0.03,0.01),c(0,5,1)), #c(0,5),c(0,5)),
  flag.limits   = F,  # no limits for the scale, they change for different durations.
  scale_grad    = list(list(name="turbo"),
                       list(name="turbo"),
                       list(name="turbo")),
                        # list(name="turbo"),
                        # list(name="turbo")),
  crs_map       = crswgs84,
  title_x       = "Longitude (°E)",
  title_y       = "Latitude (°N)",
  breaks_x      = c(6.00001,10,14,18),
  breaks_y      = c(36.00001,40,44,47.9999),
  labels_x      = c("6°E","10°E","14°E","18°E"),
  labels_y      = c("36°N","40°N","44°N","48°N"),
  limits_x      = c(6,19),
  limits_y      = c(36,48),
  alpha_shpworld= 0.2,
  size_points   = 0.7,
  flag.saveRDS  = F)


# GEV STATIONARY MODEL:
list_plot_scaleint_gev=plot_maps_param(
  var_list      = list(a_C_MAP_gev,a_C_CI_gev),#a_C_med_gev), #a_C_fmins_gev,a_C_MaxLik_gev),
  lon_tmp       = lon_tmp,
  lat_tmp       = lat_tmp,
  durations     = durations,
  mask_shp      = land_shp,
  mask_world    = gg_world,
  flag_signif   = F, # logical
  significance  = significance,
  name_var      = "Scale intercept",
  acronym_model = "gev",
  type_estimate = c("MaxPost","CI 90%"), # "MedianPost")
  name_model    = "GEV Stationary",
  path_output   = c(paste0(dir.res_maps,"/scale_int/map_scale_int_MaxPost_gev_NOLIM"),
                    paste0(dir.res_maps,"/scale_int/map_scale_int_CI90_gev_NOLIM"),
                    paste0(dir.res_maps,"/scale_int/map_scale_int_medPost_gev_NOLIM")),
                    # paste0(dir.res_maps,"/scale_int/map_scale_int_MaxLikfmins_gev"),
                    # paste0(dir.res_maps,"/scale_int/map_scale_int_MaxLik_gev")),
  limits        = list(c(0,5,1),c(0.01,0.03,0.01),c(0,5,1)), #c(0,5),c(0,5)),
  flag.limits   = F,  # no limits for the scale, they change for different durations.
  scale_grad    = list(list(name="turbo"),
                       list(name="turbo"),
                       list(name="turbo")),
                        # list(name="turbo"),
                        # list(name="turbo")),
  crs_map       = crswgs84,
  title_x       = "Longitude (°E)",
  title_y       = "Latitude (°N)",
  breaks_x      = c(6.00001,10,14,18),
  breaks_y      = c(36.00001,40,44,47.9999),
  labels_x      = c("6°E","10°E","14°E","18°E"),
  labels_y      = c("36°N","40°N","44°N","48°N"),
  limits_x      = c(6,19),
  limits_y      = c(36,48),
  alpha_shpworld= 0.2,
  size_points   = 0.7,
  flag.saveRDS  = F)



#^******************************************************************************
# Scale slope:
#^******************************************************************************
# NONSTATIONARY MODEL:
list_plot_scaleslope_ns = plot_maps_param(
  var_list      = list(b_C_MAP_ns,b_C_CI_ns,b_C_fmins_ns), # b_C_med_ns,b_C_MaxLik_ns),
  lon_tmp       = lon_tmp,
  lat_tmp       = lat_tmp,
  durations     = durations,
  mask_shp      = land_shp,
  mask_world    = gg_world,
  flag_signif   = T, # logical
  significance  = significance,
  name_var      = "Scale slope",
  type_estimate = c("MaxPost","CI 90%","MaxLik fmins"), # "MedianPost","MaxLik"),
  name_model    = "SMEV Nonstationary",
  acronym_model = "ns",
  path_output   = c(paste0(dir.res_maps,"/scale_slope/map_scale_slope_MaxPost_ns_NOLIM"),
                    paste0(dir.res_maps,"/scale_slope/map_scale_slope_CI90_ns_NOLIM"),
                    paste0(dir.res_maps,"/scale_slope/map_scale_slope_MaxLikfmins_ns_NOLIM")),
                    # paste0(dir.res_maps,"/scale_slope/map_scale_slope_medPost_ns"),
                    # paste0(dir.res_maps,"/scale_slope/map_scale_slope_MaxLik_ns")),
  limits        = list(c(-0.02,0.2,0.01),c(0.01,0.03,0.005),c(-0.02,0.02,0.01)), #c(-0.02, 0.02),c(-0.02, 0.02)),
  flag.limits   = F,
  scale_grad    = list(list(name="rdbl",low="red",high="blue",mid="white",midpoint=0),
                       list(name="turbo"),
                       list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)),
                        # list(name="rdbl",low="red",high="blue",mid="white",midpoint=0),
                        # list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)),
  crs_map       = crswgs84,
  title_x       = "Longitude (°E)",
  title_y       = "Latitude (°N)",
  breaks_x      = c(6.00001,10,14,18),
  breaks_y      = c(36.00001,40,44,47.9999),
  labels_x      = c("6°E","10°E","14°E","18°E"),
  labels_y      = c("36°N","40°N","44°N","48°N"),
  limits_x      = c(6,19),
  limits_y      = c(36,48),
  alpha_shpworld= 0.2,
  size_points   = 0.7,
  flag.saveRDS  = F)


# with fix limits:
list_plot_scaleslope_ns_FIX=plot_maps_param(
  var_list      = list(b_C_MAP_ns,b_C_CI_ns,b_C_fmins_ns), # b_C_med_ns,b_C_MaxLik_ns),
  lon_tmp       = lon_tmp,
  lat_tmp       = lat_tmp,
  durations     = durations,
  mask_shp      = land_shp, # test,
  mask_world    = gg_world,
  flag_signif   = T, # logical
  significance  = significance,
  name_var      = "Scale slope",
  type_estimate = c("MaxPost","CI 90%","MaxLik fmins"), #"MedianPost","MaxLik"),
  name_model    = "SMEV Nonstationary",
  acronym_model = "ns",
  path_output   = c(paste0(dir.res_maps,"/scale_slope/map_scale_slope_MaxPost_ns"),
                    paste0(dir.res_maps,"/scale_slope/map_scale_slope_CI90_ns"),
                    paste0(dir.res_maps,"/scale_slope/map_scale_slope_MaxLikfmins_ns")),
                    # paste0(dir.res_maps,"/scale_slope/map_scale_slope_medPost_ns"),
                    # paste0(dir.res_maps,"/scale_slope/map_scale_slope_MaxLik_ns")),
  limits        = list(c(-0.1,0.1,0.05),c(0.01,0.03,0.005),c(-0.1,0.1,0.05)), #c(-0.02,0.02),c(-0.02,0.02)),
  flag.limits   = T,
  scale_grad    = list(list(name="rdbl",low="red",high="blue",mid="white",midpoint=0),
                       list(name="turbo"),
                       list(name="rdbl",low="red",high="blue",mid="white", midpoint=0)),
                       #list(name="rdbl",low ="red",high="blue",mid="white",midpoint=0),
                       #list(name="rdbl",low ="red",high="blue",mid="white",midpoint=0)),
  crs_map       = crswgs84,
  title_x       = "Longitude (°E)",
  title_y       = "Latitude (°N)",
  breaks_x      = c(6.00001,10,14,18),
  breaks_y      = c(36.00001,40,44,47.9999),
  labels_x      = c("6°E","10°E","14°E","18°E"),
  labels_y      = c("36°N","40°N","44°N","48°N"),
  limits_x      = c(6,19),
  limits_y      = c(36,48),
  alpha_shpworld= 0.2,
  size_points   = 0.7,
  flag.saveRDS  = F)



#^******************************************************************************
#save all plots to .RDS Rdata file:
#^******************************************************************************
saveRDS(list(list_plot_shapeint_ns   =list_plot_shapeint_ns,
             list_plot_shapeint_stat =list_plot_shapeint_stat,
             list_plot_shapeint_gev  =list_plot_shapeint_gev,
             list_plot_shapeslope_ns =list_plot_shapeslope_ns,
             list_plot_scaleint_ns   =list_plot_scaleint_ns,
             list_plot_scaleint_stat =list_plot_scaleint_stat,
             list_plot_scaleint_gev  =list_plot_scaleint_gev,
             list_plot_scaleslope_ns =list_plot_scaleslope_ns
             ),
        file=paste0(dir.res_maps,'list_plot_smevparam_',sat_prod,'.rds'))

all.plots_smev=plot_grid(
  list_plot_shapeint_ns$map_var_ALLPIX_allvar_horiz$MaxPost,
  list_plot_shapeslope_ns$map_var_ALLPIX_allvar_horiz$MaxPost,
  list_plot_scaleint_ns$map_var_ALLPIX_allvar_horiz$MaxPost,
  list_plot_scaleslope_ns$map_var_ALLPIX_allvar_horiz$MaxPost,
  label_size=20,labels=NULL,ncol=1,nrow=4,rel_widths=c(1,1,1,1),
  #rel_heights=c(1.05,0.98,0.98,0.98),
  rel_heights = c(1,1,1,1),align='hv')
ggsave(all.plots_smev,
       filename=paste0(dir.res_maps,'all_maps_param_ALLPIX_',sat_prod,'.png'),
       width=26,height=29,dpi=200,bg="white")

all.plots_smev_LAND=plot_grid(
  list_plot_shapeint_ns$map_var_LAND_allvar_horiz$MaxPost,
  list_plot_shapeslope_ns$map_var_LAND_allvar_horiz$MaxPost,
  list_plot_scaleint_ns$map_var_LAND_allvar_horiz$MaxPost,
  list_plot_scaleslope_ns$map_var_LAND_allvar_horiz$MaxPost, # +theme(plot.title=element_text(size=20,hjust=0.5)),
  label_size=20,labels=NULL,ncol=1,nrow=4,rel_widths=c(1,1,1,1),
  # rel_heights=c(1.05,0.98,0.98,0.98),
  rel_heights=c(1,1,1,1),align='hv')
ggsave(all.plots_smev_LAND, 
       filename=paste0(dir.res_maps,'all_maps_param_LAND_',sat_prod,'.png'),
       width=26,height=29,dpi=200,bg="white")













################################################################################
# PLOT CHANGES % of SMEV PAR (MAXPOST) AMAX, quant2yrs, Shape, Scale, n 
################################################################################
# list required for the plots:
var_change_list=list(AMAX_change   = AMAX_change, # change in AMAX
                     qnt2yr_change = quant_smev_ns_2yr_change, # change in 2yrs return level
                     a_C_MAP_ns    = a_C_MAP_ns, # scale intercept maxpost NS
                     b_C_MAP_ns    = b_C_MAP_ns, # scale slope maxpost NS
                     a_w_MAP_ns    = a_w_MAP_ns, # shape intercept maxpost NS
                     b_w_MAP_ns    = b_w_MAP_ns, # shape slope maxpost NS
                     n_change      = n_change)   # change in n

plot_maps_changes_amax2yrs(var_change_list= var_change_list,
                           lon_tmp        = lon_tmp,
                           lat_tmp        = lat_tmp,
                           durations      = durations,
                           mask_shp       = land_shp, # test,
                           mask_world     = gg_world,
                           flag_signif    = T,
                           significance   = significance,
                           sat_acron      = sat_acron,
                           period         = '2001-2024',
                           year_ref       = year_halfperiod,
                           year_ref_id    = halfperiod,
                           name_file_txt  = paste0('signif_pixels_',sat_prod,'.txt'),
                           path_output    = paste0(dir.res_maps, '/amax_vs_qnt/'),
                           limits         = list(c(-20,20),c(-20,20),c(-20,20),c(-20,20),c(-20,20)),
                           flag.limits    = T,
                           crs_map        = crswgs84,
                           title_x        = "Longitude (°E)",
                           title_y        = "Latitude (°N)",
                           breaks_x       = c(6.00001,10,14,18),
                           breaks_y       = c(36.00001,40,44,47.9999),
                           labels_x       = c("6°E","10°E","14°E","18°E"),
                           labels_y       = c("36°N","40°N","44°N","48°N"),
                           limits_x       = c(6,19),
                           limits_y       = c(36,48),
                           alpha_shpworld = 0.2,
                           size_points    = 0.5,
                           flag.saveRDS   = F)





################################################################################
# PLOT UNCERTAINTY of %changes of 2yrs quantile, shape, scale (NONSTAT SMEV)
################################################################################
# list required for the plots:
var_unc_change_list=list(unc_qnt2yr_change=list(unc_qnt2yrs_1h,
                                                unc_qnt2yrs_3h,
                                                unc_qnt2yrs_6h,
                                                unc_qnt2yrs_24h),
                         unc_shape_change =list(unc_shape_1h,
                                                unc_shape_3h,
                                                unc_shape_6h,
                                                unc_shape_24h),
                         unc_scale_change =list(unc_scale_1h,
                                                unc_scale_3h,
                                                unc_scale_6h,
                                                unc_scale_24h))
plot_maps_changes_uncert(var_unc_change_list=var_unc_change_list,
                         lon_tmp        = lon_tmp,
                         lat_tmp        = lat_tmp,
                         durations      = durations,
                         mask_shp       = land_shp,
                         mask_world     = gg_world,
                         flag_signif    = T,
                         significance   = significance,
                         sat_acron      = sat_acron,
                         period         = '2001-2024',
                         year_ref       = year_halfperiod,
                         year_ref_id    = halfperiod,
                         path_output    = paste0(dir.res_maps, '/uncert/'),
                         limits         = list(c(20,50,10),c(20,30,2),c(20,50,10)),
                         flag.limits    = T,
                         crs_map        = crswgs84,
                         nbin           = 20,
                         bar.direction  = "horizontal",
                         legend.position= "bottom",
                         barwidth       = 25,
                         barheight      = 0.6,
                         title_x        = "Longitude (°E)",
                         title_y        = "Latitude (°N)",
                         breaks_x       = c(6.00001,10,14,18),
                         breaks_y       = c(36.00001,40,44,47.9999),
                         labels_x       = c("6°E","10°E","14°E","18°E"),
                         labels_y       = c("36°N","40°N","44°N","48°N"),
                         limits_x       = c(6,19),
                         limits_y       = c(36,48),
                         alpha_shpworld = 0.2,
                         size_points    = 0.5,
                         flag.saveRDS   = F)





################################################################################
# PLOT %changes of 2-20-50yrs quantiles (NONSTAT SMEV MAXPOST)
################################################################################
# list required for the plots:
var_qnt_change_list=list(quant_smev_ns_2yr_change=quant_smev_ns_2yr_change,
                         quant_smev_ns_10yr_change=quant_smev_ns_10yr_change,
                         quant_smev_ns_20yr_change=quant_smev_ns_20yr_change,
                         quant_smev_ns_50yr_change=quant_smev_ns_50yr_change)

plot_maps_changes_qnt(var_qnt_change_list=var_qnt_change_list,
                      lon_tmp        = lon_tmp,
                      lat_tmp        = lat_tmp,
                      durations      = durations,
                      mask_shp       = land_shp,
                      mask_world     = gg_world,
                      flag_signif    = T,
                      significance   = significance,
                      sat_acron      = sat_acron,
                      period         = '2001-2024',
                      year_ref       = year_halfperiod,
                      year_ref_id    = halfperiod,
                      path_output    = paste0(dir.res_maps,'/quantiles/'),
                      limits         = list(c(-20,20),c(-20,20),c(-20,20)),
                      limits_mean    = list(c(-15,15),c(-15,15),c(-15,15)),
                      flag.limits    = T,
                      crs_map        = crswgs84,
                      nbin           = 20,
                      bar.direction  = "horizontal",
                      legend.position= "bottom",
                      barwidth       = 25,
                      barheight      = 0.6,
                      size_text_plot = 6,
                      title_x        = "Longitude (°E)",
                      title_y        = "Latitude (°N)",
                      breaks_x       = c(6.00001,10,14,18),
                      breaks_y       = c(36.00001,40,44,47.9999),
                      labels_x       = c("6°E","10°E","14°E","18°E"),
                      labels_y       = c("36°N","40°N","44°N","48°N"),
                      limits_x       = c(6,19),
                      limits_y       = c(36,48),
                      alpha_shpworld = 0.2,
                      size_points    = 0.5,
                      flag.saveRDS   = F)




################################################################################
# plot of average quantiles per region with different durations
################################################################################
# list required for the plots:
var_durat_change_list= list(qnt2yr_change=quant_smev_ns_2yr_change, # change in 2yrs return level
                            a_C_MAP_ns   =a_C_MAP_ns, # scale intercept maxpost NS
                            b_C_MAP_ns   =b_C_MAP_ns, # scale slope maxpost NS
                            a_w_MAP_ns   =a_w_MAP_ns, # shape intercept maxpost NS
                            b_w_MAP_ns   =b_w_MAP_ns, # shape slope maxpost NS)
                            n_change     =n_change)   # change in n
user_polygons<-read_sf(dsn=dir_shp_user)

plot_maps_changes_durat(
              var_durat_change_list=var_durat_change_list,
              lon_tmp        = lon_tmp,
              lat_tmp        = lat_tmp,
              durations      = durations,
              mask_shp       = land_shp,
              mask_world     = gg_world,
              mask_user      = user_polygons,
              flag_signif    = T,
              significance   = significance,
              sat_acron      = sat_acron,
              period         = '2001-2024',
              year_ref       = year_halfperiod,
              year_ref_id    = halfperiod,
              path_output    = paste0(dir.res_maps, '/durat/'),
              limits         = c(-30,30),
              flag.limits    = T,
              crs_map        = crswgs84,
              title_x        = "Longitude (°E)",
              title_y        = "Latitude (°N)",
              breaks_x       = c(6.00001,10,14,18),
              breaks_y       = c(36.00001,40,44,47.9999),
              labels_x       = c("6°E","10°E","14°E","18°E"),
              labels_y       = c("36°N","40°N","44°N","48°N"),
              limits_x       = c(6,19),
              limits_y       = c(36,48),
              alpha_shpworld = 0.2,
              flag.saveRDS   = F)


