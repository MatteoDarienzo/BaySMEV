library(trend)


################################################################################
trend_AMAX =function(s, 
                     alpha=0.05, 
                     dir.res_trend, 
                     dur=NULL){
################################################################################
  # initialization:
  flag_save=flag_like=signif_trend=res_MK=senss=sensslope=0
  
  
  # MANN-KENDALL AND SENSS SLOPE TESTS
  ####################################
  # The null hypothesis is that the data come from a population with 
  # independent realizations and are identically distributed.
  # For the two sided test, the alternative hypothesis is that the data 
  # follow a monotonic trend.
  print('Mann-Kendall trend analysis on AMAX:')
  # If continuity = TRUE then a continuity correction will be employed:
  res_MK = trend::mk.test(na.omit(s$AMS[,dur]), 
                          alternative = c("two.sided", "greater", "less"), continuity = TRUE)
  
  # Computes Sen’s slope for linear rate of change and corresponding confidence intervalls
  # print("Using R package trend with function 'sens.slope()'... ")
  senss = sens.slope(na.omit(s$AMS[,dur]))
  
  
  
  # Significance level at alpha (e.g., 0.05)
  ##########################################
  signif_trend=F
  if (res_MK$p.value <= alpha){
    # Significant trend detected:
    print('Signif. MK trend detected on AMAX!')
    # Sen's Slope test (ss>0 positive trend, <0 negative trend):
    if (senss$estimates[["Sen's slope"]] >0){
      print("Signif. Increasing significant trend detected on AMAX!")
    } else {
      print("Signif. Decreasing significant trend detected on AMAX!")
    }
    signif_trend = T
  } else {
    print('No trend detected on AMAX')
    
  }
  
  if (signif_trend == T){
    flag_save=T
    flag_like=T
    print('Perform all calculations')
  }
  

  # plot trend:
  #############
  if ((!is.null(s$row)) | (!is.null(s$col)) | (!is.null(s$sat_prod))){
    nfplot_trend = paste0(dir.res_trend, "/trend_AMAX_", durations[dur]/60, "h_", s$sat_prod, "_", s$row,"_",s$col,".png")
    title_trend_plot = paste0('Trend on AMAX rainfall - duration ',
                        s$durations[dur]/60, 'h [mm/h] \n', 
                        '(', s$sat_prod, ', ', year(s$time[1]), ' - ', 
                        year( tail(s$time,1)), ', pix:', s$row, '-', s$col, ')')
  } else {
    # generic case
    nfplot_trend= paste0(dir.res_trend, "/trend_AMAX_", durations[dur]/60, "h_", s$name_project, ".png")
    title_trend_plot = paste0('Trend on AMAX rainfall - duration ', 
                        s$durations[dur]/60, 'h [mm/h] \n', 
                        '(', s$name_project, ", ", year(s$time[1]), ' - ',  year( tail(s$time,1)), ')')
  }
  
  print(paste0("saving figure with trend at: ", nfplot_trend))
  df_trend_AMAX = data.frame(t=unique(s$t), y=na.omit(s$AMS[,dur]), date=s$blocks)

  trendAMAXplot = trendAMAX.plot(df_trend_AMAX = df_trend_AMAX,
                                 res_MK        = res_MK,
                                 senss         = senss,
                                 title         = title_trend_plot,
                                 pathout       = nfplot_trend)
  trendAMAXplot$duration = paste0(s$durations[dur]/60, 'h')

  return(list(trendAMAXplot= trendAMAXplot,
              flag_save    = flag_save,
              flag_like    = flag_like,
              signif_trend = signif_trend,
              res_MK       = res_MK,
              senss        = senss,
              sensslope    = senss$estimates[["Sen's slope"]]))
}















################################################################################
trend_n = function(s_obj         = list(), 
                   alpha         = 0.05, 
                   dir.res_trend = '', 
                   dur           = 1,
                   nfplot        ='.png',
                   title_plot    =''){
################################################################################
  
  # Linear regression:
  ####################
  n_lmfit = lm(s_obj$ncount[[dur]] ~ s_obj$years)
  
  # plot lin regression:
  png(filename = paste0(dir.res_trend, "/lin_regr_n_", paste0(durations[dur]/60, "h"),
                        nfplot), width=700, height=500)
  plot(x    = s_obj$years_unique,
       y    = s_obj$ncount[[dur]], 
       main = paste0("Number of ordinary events per year, n, \n duration=",
                    durations[dur]/60, "h ", title_plot), 
       xlab="Year", ylab="number of ord. events, n", lwd=2,
       cex.lab=1.6, cex.axis=1.2, cex=1.5, cex.main=1.5, cex.sub=2)
  lines(x = s_obj$years_unique, 
        y = n_lmfit$coefficients[1] + n_lmfit$coefficients[2]*s_obj$years, 
        col='red', lwd=2)
  dev.off()
  
  
  
  # MK trend on n:
  ################
  # MK test & Sen's slope test applied to n:
  n_MK = trend::mk.test(na.omit(s_obj$ncount[[dur]]), 
                        alternative = c("two.sided", "greater", "less"), 
                        continuity = T)
  
  # Computes Sen’s slope for linear rate of change and corresponding confidence intervals
  n_senss = sens.slope(na.omit(s_obj$ncount[[dur]]))
  
  df_trend_n = data.frame(t    = unique(s_obj$t),
                          y    = na.omit(s_obj$ncount[[dur]]), 
                          date = s_obj$blocks)
  
  trend_n_plot = trend_n.plot(df_trend_n = df_trend_n,
                              res_MK     = n_MK,
                              senss      = n_senss,
                              title      = paste0('Trend on number of ordinary events, n - duration ', 
                                                  durations[dur]/60, 'h [mm/h] \n', title_plot),
                              pathout    = paste0(dir.res_trend, "/trend_n_", durations[dur]/60, "h", nfplot))
  trend_n_plot$duration = paste0(durations[dur]/60, 'h')
  
  
  # Significance level at 0.05 (alpha):
  signif_trend_n=F
  if (n_MK$p.value <= alpha){
    # Significant trend on n detected:
    print('Signif. MK trend detected on n')
    # Sen's Slope test (ss>0 positive trend, <0 negative trend):
    if (n_senss$estimates[["Sen's slope"]] >0){
      print("Signif. Increasing significant trend detected on n")
    } else {
      print("Signif. Decreasing significant trend detected on n!")
    }
    signif_trend_n = T
  } else {
    print('No trend detected on n')
  }

  
  return(list(n_lmfit        = n_lmfit,
              trend_n_plot   = trend_n_plot,
              signif_trend_n = signif_trend_n,
              n_MK           = n_MK,
              n_senss        = n_senss,
              n_sensslope    = n_senss$estimates[["Sen's slope"]]))
}
