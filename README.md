# BaySMEV
BAYesian Simplified Metastatistical Extreme Value method

v. 1.0.0
Author: Matteo Darienzo (University of Padua, Padua, Italy)

Contributions from: Francesco Marra.

Year: 2025

Funders: NextEU PNRR (PRIN2022 "INTENSE" project)

Last update: 16/10/2025

Purpose:

Software BaySMEV (BAYesian Simplified Metastatistical Extreme Value) provides tools for rainfall extreme analysis with both stationary and nonstationary models within a Bayesian framework, which provides results with qunatitative uncertainty.
It is the Bayesian version of SMEV method (Marra et al., 2019, Marani and Ignaccolo, 2015).
Source codes are written in R.


Acknowledegments:

This research has been supported by the INTENSE project (raINfall exTremEs and their impacts: from the local to the National ScalE) funded by the European Union - Next Generation EU in the framework of PRIN (Progetti di ricerca di Rilevante Interesse Nazionale) programme (grant 2022ZC2522).
BaySMEV codes are used for paper in preparation Darienzo et al. ("Recent changes in extreme precipitation across Italy using satellite products").
BaySMEV codes make use of other referenced codes of Marra et al., 2019 downloadable from https://zenodo.org/records/11934843 (SMEV methodology) and https://zenodo.org/records/15047817 (non-stationary model) and of Darienzo et al., 2021 downloadable from https://github.com/MatteoDarienzo/BayDERS (segmentation method) and https://github.com/BaM-tools/BaM (Bayesian modelling).

It makes use of the following R packages: "rstudioapi",
         "methods", "lattice", "gridExtra", "reshape","reshape2", "ggplot2", "extrafont",
         "grid", "gtable", "chron", "coda", "RColorBrewer", "cowplot", "viridis",
         "psych",  "mosaicData", "tidyr", "lubridate", "Kendall", "trend", "tidyverse", 
         "latex2exp", "dplyr", "png","pracma",  "CFtime", "ggpubr",
         "raster", "maptools","sf", "terra", "exactextractr", "ncdf4",  "rnaturalearth".
