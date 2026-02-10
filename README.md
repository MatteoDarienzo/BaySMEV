# BaySMEV
BAYesian Simplified Metastatistical Extreme Value method

- v. 1.0.0

- Author: Matteo Darienzo (University of Padua, Padua, Italy)

- Contributions from: Francesco Marra (University of Padua, Padua, Italy).

- Year: 2025

- Funders: NextEU PNRR (PRIN2022 "INTENSE" project)

- Last update: 10/01/2026


# Purpose:

Software BaySMEV (BAYesian Simplified Metastatistical Extreme Value) provides the Bayesian version of SMEV method (Marra et al., 2019, Marani and Ignaccolo, 2015). 
It contains tools for rainfall extreme analysis with both stationary and non-stationary models within a Bayesian framework. MLE method is also available (as per the original SMEV method of Marra et al., 2019, https://zenodo.org/records/15047817).


# Method:

- It identifies all ordinary events within the given precipitation time series, by identifying the maximum value within the wet periods. Wet periods are defined by instants with precipitation > minimum threshold (e.g., 0.1 mm) separated by a user-defined minimum dry period (e.g., 24h).

- It computes the cdf of extremes as F^n, where F is the cdf of the ordinary events and n is the occurrence frequency of ordinary events (number of events per year).

- F is assumed as a 2-parameters Weibull distribution with a left-censoring approach (https://zenodo.org/records/11934843).

- Scale and shape parameters of the Weibull distribution are estimated by means of (choice of the user):  1) MLE with fminsearch minimization function or 2) Bayesian inference with STAN package.

- BaySMEV source codes are written in R.

- An example of main launcher ("smev_main_EXAMPLE.R") and associated ".rds" R data file with precipitation timeseries are provided in order to run the script and test its functionalities and different settings. For the use of BaySMEV to gridded datasets please refers to main launcher "smev_main_sat_IMERG.R" which performs SMEV analysis (stationary or nonstationary) for each pixel.

- For creating maps with results of the SMEV analysis computed on gridded datasets, please use postprocessing module "/BaySMEV/postprocessing/smev_results_maps_IMERG.R".



# Acknowledegments:

This research has been supported by the INTENSE project (raINfall exTremEs and their impacts: from the local to the National ScalE) funded by the European Union - Next Generation EU in the framework of PRIN (Progetti di ricerca di Rilevante Interesse Nazionale) programme (grant 2022ZC2522).

BaySMEV codes have been developed and used for paper Darienzo et al. ("Contrasting changes in extreme hourly precipitation across Italy revealed by satellite and re-analysis products, in preparation).

BaySMEV takes inspiration from matlab codes of Marra et al., 2019 downloadable from https://zenodo.org/records/11934843 (SMEV methodology, including storm separation) and https://zenodo.org/records/15047817 (for the non-stationary model implementation).

BaySMEV for the Bayesian and MCMC algoryhtm takes ideas from methods: RsTooDS of Benjamin Renard, INRAE (https://zenodo.org/records/5075760); BaM of Benjamin Renard, INRAE (https://github.com/BaM-tools/RBaM); BayDERS of Matteo Darienzo, INRAE (https://github.com/MatteoDarienzo/BayDERS); packages "rstan" and "cmdstanr" for stan Bayesian inference and diagnostics, package "adaptMCMC" for adaptive metropolis algorithm.

BaySMEV uses R package 'trend' of Thorsten Pohlert (https://cran.r-project.org/web/packages/trend/index.html) for the trend analysis (MK and Sen's slope).


BaySMEV makes use of the following R packages: "rstudioapi",
         "methods", "lattice", "gridExtra", "reshape","reshape2", "ggplot2", "extrafont",
         "grid", "gtable", "chron", "coda", "RColorBrewer", "cowplot", "viridis",
         "psych",  "mosaicData", "tidyr", "lubridate", "Kendall", "trend", "tidyverse", 
         "latex2exp", "dplyr", "png","pracma",  "CFtime", "ggpubr",
         "raster", "maptools","sf", "terra", "exactextractr", "ncdf4",  "rnaturalearth",
         "parallel", "scales", "mcmc", "evd", "stats", "R.matlab", "extRemes", "imager",
         "mmand", "adaptMCMC", "PEIP", "posterior", "loo", "bayesplot", "rstan"

# Disclaimer:

Please, notice that BaySMEV is an experimental software. Further testing and investigation are required for validating the proposed analysis and codes.
We are not responsible for any loss (of data, profits, business, customers, orders, other costs or disturbances) derived by their use in the research or operational practice. 
The authorized user accepts the risks associated with using the software given its nature of free software. 
It is reserved to expert Users (developers or professionals) having prior knowledge in computer science, hydrology and statistics.


For any question or feedback please contact us at:

- matteo.darienzo@unipd.it
- francesco.marra@unipd.it



# Example (main function, see "smev_main_EXAMPLE.R"):

Main function for SMEV analysis with most of settings with values by default (from "/module_smev_main.R").

         s_base <- smev_analysis_main(rain.df           = rain.df,
                                      name_project      = 'DEMO',
                                      dir_results       = "",
                                      date_initial      = "2001-01-01 00:00:00", 
                                      date_final        = "2024-12-31 23:00:00",
                                      time_resolution   = 60,
                                      min_ev_duration   = 30,
                                      separation_in_min = 1440,
                                      duration          = c(1,3,6,24)*60,
                                      thr_leftcens      = 0.9,
  	  	                      estim_algorithm   = "cmdstan",                  
                                      priors            = priors, # includes starting values!                               
                                      flag_smev_stat    = T,
                                      flag_smev_ns      = F,
                                      flag_smev_partns  = F,
                                      )
where:

- rain.df --> [dataframe] input rainfall timeseries dataframe. One column 'timestamps' with the dates in timestamp, and another column named 'y' with the rainfall values [mm] [real]. 

- name_project --> [character] name of the folder where the results will be saved.

- dir_results --> [character] path to the folder where to the save results.

- date_final --> [timestamp] final date for the analysis (ex: "2024-12-31 23:00:00")

- date_initial --> [timestamp] initial date for the analysis (ex: "2001-01-01 00:00:00")

- time_resolution --> [integer] Temporal resolution of the timeseries.

- min_ev_duration --> [integer] duration in minutes of an event, e.g., 30.

- separation_in_min --> [integer] min time for event separation in minutes, e.g., 1440 (one day)

- duration --> [array] defines the rainfall durations (in minutes) that you want to analyse, e.g., c(1,3,6,24)*60 for 1h, 3h,6h,24h.

- thr_leftcens --> [real] threshold percentile for the left censoring approach (e.g., 0.9 to censor 90% of data).  

- estim_algorithm --> [character] choice for the parameters estimation method is between: "cmdstan" (preferred) and "rstan" for Bayesian estimation with HMC mcmc alogrithm, "adaptMCMC" for Bayesian estimation with adaptive MCMC algorithm, and "fminsearch" for MLE with minimization algorithm.

- priors --> [list of list] list of priors for smev the parameters (shape and scale) one list for each defined duration. 
         Ex:
       
         #1h:
         priors[[1]]=list()
         priors[[1]]$shape_intercept = list("lognormal",log(0.7),0.4,0.7)
         priors[[1]]$shape_slope     = list("normal", 0,0.05,0)
         priors[[1]]$scale_intercept = list("lognormal",1.45-0.7*log(1),0.5,exp(1.45-0.7*log(1)))
         priors[[1]]$scale_slope     = list("normal",0,0.1039-0.0039*(1),0)

- start_values --> [array] array of starting values for MCMC algorithm for posterior sampling. one value for each SMEV parameter: shape int, shape slope, scale int, scale slope = c(0.7, 0, 5, 0),

- flag_smev_stat, flag_smev_ns, flag_smev_partns --> [logical] flags (T or F) fot computing different types of SMEV model: stationary, non-stationary, partial non-stationary modelds.




# Example (main with full settings):
Otherwise, the same function but with all settings is:

         s_all <- smev_analysis_main(rain.df           = rain.df,
                                     name_project      = 'DEMO',
                                     dir_results       = '',
                                     date_initial      = "2001-01-01 00:00:00", 
                                     date_final        = "2024-12-31 23:00:00",
                                     min_rain          = 0.1,
                                     time_resolution   = 60,
                                     min_ev_duration   = 30,
                                     separation_in_min = 1440,
                                     duration          = c(1,3,6,24)*60,
                                     target_T          = sort(c(exp(seq(log(1.1),log(50),.1)),2,5,10,20,50,100)),
                                     thr_leftcens      = 0.9,                   
                                     estim_algorithm   = "cmdstan",
                                     priors            = priors,  # includes starting values!                              
                                     priors_gev        = priors_gev,  # includes starting values!                           
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
                                     flag_gevstat      = T,
                                     flag_smev_stat    = T,
                                     flag_smev_ns      = T,
                                     flag_smev_partns  = F,
                                     flag_model_compar = T,
                                     flag_save2rds     = T)

where the additional inputs are:

- priors_gev --> [list of array] list of priors for GEV the parameters (location, scale and shape).

         priors_gev=NULL # initialise
         priors_gev$loc   = list("uniform",0,100,1)
         priors_gev$scale = list("uniform",0,100,1)
         priors_gev$shape = list("normal",0,0.3,0.1) # Shape ==> Martins et al.2000, Beta distrib (mean=0.1,stdev=0.122) and Renard et al.,2008, N(mean=0,sd=0.3) quite poor informative.

- Nmcmc --> [integer] Number of MCMC interations.

- Nmcmc_max --> [integer] ONLY for "adaptMCMC" estimation methods. Maximum number of MCMC iterations when not converging well.
  
- acc.rate --> [real] ONLY for "adaptMCMC" estimation methods. Acceptance rate for MCMC sampling algorithm. By default 0.234.

- Nslim --> [integer] ONLY for "adaptMCMC" estimation methods.

- nburn --> [real] ONLY for "adaptMCMC" estimation methods.

- scales --> [real] ONLY for "adaptMCMC" estimation methods.

- scales_gev --> [real] ONLY for "adaptMCMC" estimation methods.

- pix_row --> [integer] ONLY if cycling over a gridded dataset to define pixel row index.
  
- pix_col -->  [integer] ONLY if cycling over a gridded dataset to define pixel column index.

- flag_like --> [logical] T/F, T if you want to compare MaxPost vs. MLE when using Bayesian estimation method. 

- flag_gevstat --> [logical] T/F, T if you want to infer also GEV parameters for comparison with SMEV.

- flag_model_compar --> [logical] T/F, T if you want to compare the different models SMEV, MEV, GEV.

- flag_save2rds --> [logical] T/F, T if you want to save all results to a R data file format .rds

