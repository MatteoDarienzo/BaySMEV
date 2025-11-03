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
It contains tools for rainfall extreme analysis with both stationary and non-stationary models within a Bayesian framework, which provides results with quantitative uncertainty. 



# Method:

- It identifies all ordinary events within the given precipitation time series, by identifying the maximum value within the wet periods. Wet periods are defined by instants with precipitation > minimum threshold (e.g., 0.1 mm) separated by a user-defined minimum dry period (e.g., 24h).

- It computes the cdf of extremes as F^n, where F is the cdf of the ordinary events and n is the occurrence frequency of ordinary events (number of events per year).

- F is assumed as a 2-parameters Weibull distribution with a left-censoring approach.

- Scale and shape parameters of the Weibull distribution are estimated by means of 1) MLE with fminsearch minimization function or 2) Bayesian inference with STAN package.

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
