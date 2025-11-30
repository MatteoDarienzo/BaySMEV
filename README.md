# BaySMEV
BAYesian Simplified Metastatistical Extreme Value method

- v. 1.0.0

- Author: Matteo Darienzo (University of Padua, Padua, Italy)

- Contributions from: Francesco Marra (University of Padua, Padua, Italy).

- Year: 2025

- Funders: NextEU PNRR (PRIN2022 "INTENSE" project)

- Last update: 16/10/2025


# Purpose:

Software BaySMEV (BAYesian Simplified Metastatistical Extreme Value) provides the Bayesian version of SMEV method (Marra et al., 2019, Marani and Ignaccolo, 2015). 
It contains tools for rainfall extreme analysis with both stationary and non-stationary models within a Bayesian framework. The MLE method is also here available as per the version of Marra et al., 2019 (https://zenodo.org/records/15047817).


# Method:

- It identifies all ordinary events within the given precipitation time series, by identifying the maximum value within the wet periods. Wet periods are defined by instants with precipitation > minimum threshold (e.g., 0.1 mm) separated by a user-defined minimum dry period (e.g., 24h).

- It computes the cdf of extremes as F^n, where F is the cdf of the ordinary events and n is the occurrence frequency of ordinary events (number of events per year).

- F is assumed as a 2-parameters Weibull distribution with a left-censoring approach.

- Scale and shape parameters of the Weibull distribution are estimated by means of  (choice of the user):  1) MLE with fminsearch minimization function or 2) Bayesian inference with STAN package.

- BaySMEV source codes are written in R.



# Acknowledegments:

This research has been supported by the INTENSE project (raINfall exTremEs and their impacts: from the local to the National ScalE) funded by the European Union - Next Generation EU in the framework of PRIN (Progetti di ricerca di Rilevante Interesse Nazionale) programme (grant 2022ZC2522).

BaySMEV codes have been used for paper Darienzo et al. ("Recent changes in extreme precipitation across Italy using satellite products", in preparation).

BaySMEV takes inspiration from matlab codes of Marra et al., 2019 downloadable from https://zenodo.org/records/11934843 (SMEV methodology, including storm separation) and https://zenodo.org/records/15047817 (for the non-stationary model implementation).

BaySMEV makes use of the following R packages: "rstudioapi",
         "methods", "lattice", "gridExtra", "reshape","reshape2", "ggplot2", "extrafont",
         "grid", "gtable", "chron", "coda", "RColorBrewer", "cowplot", "viridis",
         "psych",  "mosaicData", "tidyr", "lubridate", "Kendall", "trend", "tidyverse", 
         "latex2exp", "dplyr", "png","pracma",  "CFtime", "ggpubr",
         "raster", "maptools","sf", "terra", "exactextractr", "ncdf4",  "rnaturalearth".


# Disclaimer:

Please, notice that BaySMEV is an experimental software. Further testing and investigation are required for validating the proposed analysis. 
We are not responsible for any loss (of data, profits, business, customers, orders, other costs or disturbances) derived by their use in the research or operational practice. 
The authorized user accepts the risks associated with using the software given its nature of free software. 
It is reserved to expert Users (developers or professionals) having prior knowledge in computer science, hydrology and statistics.

For any question or feedback please contact us at:

- matteo.darienzo@unipd.it



# Example:

Simplified function for SMEV analysis with most of settings with values by default. Import function from "/module_smev_main.R".

         s_base <- smev_analysis_main(rain.df           = rain.df,
                                      name_project      = 'DEMO',
                                      dir_results       = dir_results,
                                      date_initial      = date_initial, 
                                      date_final        = date_final,
                                      time_resolution   = 60,
                                      min_ev_duration   = 30,
                                      separation_in_min = 1440,
                                      duration          = c(1,3,6,24)*60,
                                      thr_leftcens      = 0.9,                   
                                      priors            = priors,                               
                                      start_values      = c(0.7, 0, 5, 0),
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

- priors --> [list of array] list of priors for smev the parameters (shape and scale) one list for each defined duration. 
         ex:
       
         #1h
         priors[[1]]=list()
         priors[[1]]$shape_intercept =list("lognormal",log(0.7),0.4,0.7)
         priors[[1]]$shape_slope     =list("normal", 0,0.05,0)
         priors[[1]]$scale_intercept =list("lognormal",1.45-0.7*log(1),0.5,exp(1.45-0.7*log(1)))
         priors[[1]]$scale_slope     =list("normal",0,0.1039-0.0039*(1),0)

- start_values --> [array] array of starting values for MCMC algorithm for posterior sampling. one value for each SMEV parameter: shape int, shape slope, scale int, scale slope = c(0.7, 0, 5, 0),

- flag_smev_stat, flag_smev_ns, flag_smev_partns --> [logical] flags (T or F) fot computing different types of SMEV model: stationary, non-stationary, partial non-stationary modelds.
