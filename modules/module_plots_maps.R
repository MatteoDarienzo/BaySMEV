
################################################################################
# plot the indeces of the studied grid
################################################################################
plot_map_indexes = function(tmp_array,
                            lon_tmp,
                            lat_tmp,
                            mask_shp=NULL,
                            mask_world,
                            path_output,
                            crs_map=crswgs84){
################################################################################
  pix=0;index=c();
  for (row in 1:nrow(tmp_array)){
    for (col in 1:ncol(tmp_array)){
      pix=pix+1
      index$row[pix]=row
      index$col[pix]=col
    }
  }
  # create dataframe for the map
  indexes_df<-data.frame(x=lon_tmp,y=lat_tmp,z=index$row,h=index$col)
  # trend_df$z[is.na(trend_df$z)]=0
  # indexes_df_NEW =indexes_df
  indexes_df_NEW=data.frame(x=round(indexes_df$x,3),
                            y=round(indexes_df$y,3),
                            z=round(indexes_df$z,3),
                            h=round(indexes_df$h,3))
  # get raster from dataframe
  rast_indexes=rasterFromXYZ(indexes_df_NEW,crs=crswgs84)
  # plot
  map_indexes=ggplot()+ 
    theme_bw()+
    scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
    scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
    geom_tile(data=indexes_df_NEW,aes(x=x, y=y),fill=NA,alpha=1,color="black")+
    annotate("text",x=indexes_df_NEW$x-0.05,y=indexes_df_NEW$y,
             label=indexes_df_NEW$z,color="red",size=1.7)+
    annotate("text",x=indexes_df_NEW$x+0.05,y=indexes_df_NEW$y,
             label=indexes_df_NEW$h,color="blue",size=1.7)+
    geom_sf(data = mask_world,color="grey30",fill="black",alpha=0.2) +
    geom_sf(data = test,color="black",fill="black",alpha=0.3)+
    coord_sf(xlim=c(6, 19),ylim=c(36,48),expand=F)+
    scale_fill_gradient2(low="red",high="blue",mid="white",midpoint=0,na.value=NA)+
    theme(
       panel.grid.major=element_blank()
      ,panel.grid.minor=element_blank()
      ,axis.text.x=element_text(size=13)
      ,axis.text.y=element_text(size=13)
      ,legend.text=element_text(size=15)
      ,legend.title=element_text(size=20))

  # save map to png 
  ggsave(map_indexes,
         filename=path_output,
         width=20,height=15,dpi=200)
}












################################################################################
plot_map_quantiles <- function(plotlist, title, pathout){
################################################################################
  #^* GOAL: plot the precipitation time series
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [integer]
  #^*    2. [integer]
  #^*    3. [character]
  #^* OUT
  #^*    1. [real]
  #^****************************************************************************
  #^* REF.: 
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************
  trend_df_NEW=data.frame(x=round(trend_df$x,3),
                          y=round(trend_df$y,3), 
                          z=round(trend_df$z,3))
  rast_trend=rasterFromXYZ(trend_df_NEW, crs=crswgs84) 
  writeRaster(rast_trend,
              filename=paste0("C:/Users/MATTEO/Documents/PRIN/Qgis/trend_rain_",
                              acron_sat[sat_prod], ".tif"),
              format="GTiff",overwrite=T)
  
  ## crop and mask (only Italy extent):
  r=rast_trend
  r2<-crop(r, extent(test))
  r3<-mask(r2,test)
  test_spdf<-as(r3,"SpatialPixelsDataFrame")
  test_df<-as.data.frame(test_spdf)
  trend_df_NEW[trend_df_NEW==9999] =NA
  
  # plot map of EU with trend results 
  # (pixel level NOT AGGREGATED):  
  map_res_trend=ggplot() + 
    theme_bw()+
    scale_x_continuous(expand = c(0, 0), name = "Longitude (°E)")+
    scale_y_continuous(expand = c(0, 0), name = "Latitude (°N)")+
    geom_tile(data=trend_df_NEW, aes(x=x, y=y, fill=z), alpha=1) + 
    #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
    geom_sf(data = gg_world, fill = NA, color = "grey50") +
    geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+
    coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)+
    scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0,
                         na.value = NA, name = "Signif.\nSen's\nslope")+
    theme(
      ,legend.text     = element_text(size=15)
      ,legend.title    = element_text(size=20))
  ggsave(map_res_trend, filename =pathout, width = 15, height =15, dpi = 200)
}









################################################################################
plot_maps_gridplot <- function(plotlist, title, pathout){
################################################################################
  # create title object:
  title1 <- ggdraw() + 
    draw_label(title, fontface = 'bold', x = 0.5,   hjust = 0.5, vjust = 0.5, size=20)
  # theme(
  #   # add margin on the left of the drawing canvas,
  #   # so title is aligned with left edge of first plot
  #   plot.margin = margin(0, 0, 0, 7)
  # )
  
  # ,legend.title    = element_text(size=30),
  # ,legend.position = 'bottom'
  # ,axis.title      = element_text(size=40)
  # ,axis.text.x     = element_text(size=40)) 
  
  all.plots = plot_grid( plotlist[[1]]+theme(legend.title=element_blank()),#+theme(legend.position="none"),
                         plotlist[[2]]+theme(legend.title=element_blank()),#+theme(legend.position="none"),
                         plotlist[[3]]+theme(legend.title=element_blank()),#+theme(legend.position="none"),
                         plotlist[[4]]+theme(legend.title=element_blank()),
                         label_size=20,
                         labels=c("1h","3h","6h","24h"),
                         ncol=2,nrow=2,
                         rel_heights=c(1,1,1,1),
                         align='hv') 
  all.plots.title=plot_grid(title1,
                            all.plots,
                            ncol=1,nrow=2,
                            rel_heights=c(0.1,1))
  ggsave(all.plots.title,filename=pathout,width=14,height=16,dpi=200,bg="white")
}









# ################################################################################
# # plor maps with smev parameters 
# ################################################################################
# plot_maps_param_test= function(var_list,
#                                lon_tmp,
#                                lat_tmp,
#                                durations,
#                                mask_shp,
#                                flag_signif=T,
#                                significance,
#                                name_var="Shape intercept",
#                                type_estimate="MaxPost",
#                                name_model = "SMEV Nonstationary",
#                                path_output,
#                                limits,
#                                flag.limits=T,
#                                scale_grad,
#                                crs_map = crswgs84){
# ################################################################################
#       # create dataframe for the map plot:
#       a_w_MAP_df <- data.frame(x=lon_tmp, y=lat_tmp, z=var_list)
#       # trend_df$z[is.na(trend_df$z)]=0
#       a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
#                                   y=round(a_w_MAP_df$y,3), 
#                                   z=round(a_w_MAP_df$z,3))
#       ggplot() + 
#         theme_bw()+
#         scale_x_continuous(expand = c(0, 0), name = "Longitude (°E)")+
#         scale_y_continuous(expand = c(0, 0), name = "Latitude (°N)")+
#         geom_tile(data=a_w_MAP_df, aes(x=x, y=y, fill=z), alpha=1) + 
#         #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
#         geom_sf(data = gg_world, fill = NA, color = "grey70") +
#         geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+
#         coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)+
#         theme(
#           # legend.position = "none"
#           ,legend.text     = element_text(size=15)
#           ,legend.title    = element_text(size=20))
# }







################################################################################
plot_maps_param = function(var_list      = list(),
                           lon_tmp       = NULL,
                           lat_tmp       = NULL,
                           durations     = c(60,180,360,1440), # in minutes
                           mask_shp      = NULL,
                           mask_world    = NULL,
                           flag_signif   = F,
                           significance  = NULL,
                           name_var      = "Shape intercept",
                           type_estimate = "MaxPost",
                           name_model    = "SMEV Nonstationary",
                           acronym_model = "ns",
                           path_output   = "/map_test",
                           limits        = list(),
                           flag.limits   = F,
                           scale_grad    = list(),
                           crs_map       = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")){
################################################################################
  
  # Initialize the lists output:
  map_var_LAND_allvar<-vector(mode = "list",length=length(type_estimate))
  names(map_var_LAND_allvar)=c(type_estimate)
  map_var_ALLPIX_allvar<-vector(mode = "list",length=length(type_estimate))
  names(map_var_ALLPIX_allvar)=c(type_estimate)
  map_var_LAND_allvar_horiz<-vector(mode = "list",length=length(type_estimate))
  names(map_var_LAND_allvar_horiz)=c(type_estimate)
  map_var_ALLPIX_allvar_horiz<-vector(mode = "list",length=length(type_estimate))
  names(map_var_ALLPIX_allvar_horiz)=c(type_estimate)
  
  # a few plot settings for legend:
  nbin=20
  bar.direction="horizontal"
  legend.position="bottom"
  barwidth=23.5
  barheight=0.6
  
  
  # For each var/parameter "name_var"
  #^****************************************************************************
  for (v in 1:length(var_list)){
  #^****************************************************************************
    message(paste0("- creating maps for '",name_var,"' (",type_estimate[v],
                   "- ",name_model,") for all durations..."))
    # initialize variable plot list:
    map_var_ALLPIX=list()
    map_var_LAND=list()
    a_w_MAP_df=c()
    
    # For each rainfall duration
    #^**************************************************************************
    for (dur in 1:length(durations)){
    #^**************************************************************************
      # create dataframe for the map plot:
      a_w_MAP_df<-data.frame(x=lon_tmp, 
                             y=lat_tmp,
                             z=var_list[[v]][[dur]])
      # trend_df$z[is.na(trend_df$z)]=0
      a_w_MAP_df_NEW=a_w_MAP_df 
      # a_w_MAP_df_NEW=data.frame(x=round(a_w_MAP_df$x,3),
      #                           y=round(a_w_MAP_df$y,3),
      #                           z=round(a_w_MAP_df$z,3))
      rast_a_w_MAP=rasterFromXYZ(a_w_MAP_df_NEW,crs=crs_map)
      
      # detect pixels with signif. increasing AMAX
      signif_AMAX_pos_df<-data.frame(x=lon_tmp, 
                                     y=lat_tmp, 
                                     z=significance$signif_MK_AMAX[[dur]], 
                                     d=significance$sensslope_AMAX[[dur]])
      signif_AMAX_pos_df$z[signif_AMAX_pos_df$z==F]=NA
      signif_AMAX_pos_df$z[signif_AMAX_pos_df$d<0]=NA
      signif_AMAX_pos_df=na.omit(signif_AMAX_pos_df)
      signif_AMAX_pos_df=subset(signif_AMAX_pos_df,z!=0)  
      
      # detect pixels with signif. decreasing AMAX
      signif_AMAX_neg_df<-data.frame(x=lon_tmp, 
                                     y=lat_tmp, 
                                     z=significance$signif_MK_AMAX[[dur]], 
                                     d=significance$sensslope_AMAX[[dur]])
      signif_AMAX_neg_df$z[signif_AMAX_neg_df$z==F]=NA
      signif_AMAX_neg_df$z[signif_AMAX_neg_df$d>=0]=NA
      signif_AMAX_neg_df=na.omit(signif_AMAX_neg_df)
      signif_AMAX_neg_df=subset(signif_AMAX_neg_df,z!=0)  
      
      # detect pixels with signific nonstationarity (using DIC):
      # difference in DIC >3:
      signif_ns_DIC_df_3<-data.frame(x=lon_tmp,
                                     y=lat_tmp,
                                     z=significance$DIC_diff3[[dur]])
      # signif_ns_DIC_df_3 <- data.frame(x=round(lon_tmp,3),
      #                                  y=round(lat_tmp,3), 
      #                                  z=significance$DIC_diff3[[dur]])
      signif_ns_DIC_df_3$z[signif_ns_DIC_df_3$z==F]=NA
      signif_ns_DIC_df_3=na.omit(signif_ns_DIC_df_3)
      signif_ns_DIC_df_3=subset(signif_ns_DIC_df_3,z!=0)
      # difference in DIC >1:
      signif_ns_DIC_df_1 <-data.frame(x=lon_tmp, 
                                      y=lat_tmp, 
                                      z=significance$DIC_diff1[[dur]])
      signif_ns_DIC_df_1$z[signif_ns_DIC_df_1$z==F]=NA
      signif_ns_DIC_df_1=na.omit(signif_ns_DIC_df_1)
      signif_ns_DIC_df_1=subset(signif_ns_DIC_df_1,z!=0)
      
      
      #^************************************************************************
      #^                   ALL PIXELS (land & sea) 
      #^************************************************************************
      map_a_w_MAP = ggplot() + 
        theme_bw()+
        # scale_x_continuous(expand = c(0, 0), name = "Longitude (°E)")+
        # scale_y_continuous(expand = c(0, 0), name = "Latitude (°N)")+
        scale_x_continuous(expand = c(0, 0), name = "Longitude (°E)", 
                           breaks=c(6.00001, 10, 14, 18), 
                           labels=c("6°E", "10°E", "14°E", "18°E"))+
        scale_y_continuous(expand = c(0, 0), name = "Latitude (°N)",  
                           breaks=c(36.00001, 40, 44, 47.9999), 
                           labels=c("36°N", "40°N", "44°N", "48°N"))+
        geom_tile(data=a_w_MAP_df_NEW, aes(x=x, y=y, fill=z), alpha=1) + 
        #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
        geom_sf(data = mask_world, fill = NA, color = "grey70") +
        geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+
        coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
        
        if (flag_signif){
          # 2do: implement the significance from likelihood ratio test for fmins MLE method !!!
          map_a_w_MAP=map_a_w_MAP+
          # geom_point(data=signif_AMAX_pos_df, aes(x=x, y=y, size=d), fill="blue", color="blue", shape=24)+
          # geom_point(data=signif_AMAX_neg_df, aes(x=x, y=y, size=d), fill="red", color="red",  shape=25)+
          # geom_point(data=signif_ns_DIC_df_1, aes(x=x, y=y, size=z),   fill="gray90", color="gray90", shape=21, size = 0.7)+
          geom_point(data=signif_ns_DIC_df_3, aes(x=x, y=y, size=z),   fill="black", color="black", shape=21, size = 0.7)
          #geom_point(data=signif_ns_DIC_df, aes(x=x, y=y, size=z), shape = 1, size = 1)
        }
        map_a_w_MAP=map_a_w_MAP+
        theme(
           panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank()
          ,axis.text.x      = element_text(size=13)
          ,axis.text.y      = element_text(size=13)
          ,legend.position  = legend.position
          ,legend.text      = element_text(size=15)
          ,legend.title     = element_text(size=20))
      
        
      if (flag.limits){
        if (scale_grad[[v]]$name=="turbo"){
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_viridis(option = "turbo",
                               limits = c(limits[[v]][1],limits[[v]][2]), 
                               name   = '', # paste0(name_var,' (',type_estimate[v],'), ',name_model),
                               oob    = scales::squish,
                               breaks = seq(limits[[v]][1], limits[[v]][2], limits[[v]][3]),
                               labels = seq(limits[[v]][1], limits[[v]][2], limits[[v]][3]),
                               guide  = guide_colourbar(nbin=nbin,
                                                        raster=T,
                                                        display = "rectangles",
                                                        direction=bar.direction,
                                                        barwidth=barwidth, #grid::unit(0.6, "npc"),
                                                        barheight=barheight, # grid::unit(1, "npc"), # - grid::unit(, "line"), #40,
                                                        frame.colour=c("NA"),
                                                        frame.linewidth=1))
        } else {
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name   = '', # paste0(name_var,' (',type_estimate[v],'), ',name_model),
                                 limits = c(limits[[v]][1],limits[[v]][2]), 
                                 oob    = scales::squish,
                                 breaks = seq(limits[[v]][1], limits[[v]][2], limits[[v]][3]),
                                 labels = seq(limits[[v]][1], limits[[v]][2], limits[[v]][3]),
                                 guide  = guide_colourbar(nbin=nbin,
                                                          raster=T,
                                                          display = "rectangles",
                                                          direction=bar.direction,
                                                          barwidth=barwidth, #grid::unit(0.6, "npc"),
                                                          barheight=barheight, # grid::unit(1, "npc"), # - grid::unit(, "line"), #40,
                                                          frame.colour=c("NA"),
                                                          frame.linewidth=1))
        }
      } else {
        if (scale_grad[[v]]$name=="turbo"){
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_viridis(option="turbo",
                               name   = '', # paste0(name_var,' (',type_estimate[v],'), ',name_model),
                               guide  = guide_colourbar(nbin=nbin,
                                                        raster=T,
                                                        display = "rectangles",
                                                        direction=bar.direction,
                                                        barwidth=barwidth, #grid::unit(0.6, "npc"),
                                                        barheight=barheight, # grid::unit(1, "npc"), # - grid::unit(, "line"), #40,
                                                        frame.colour=c("NA"),
                                                        frame.linewidth=1))
        } else {
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name   = '', # paste0(name_var,' (',type_estimate[v],'), ',name_model),
                                 guide  = guide_colourbar(nbin=nbin,
                                                          raster=T,
                                                          display = "rectangles",
                                                          direction=bar.direction,
                                                          barwidth=barwidth, #grid::unit(0.6, "npc"),
                                                          barheight=barheight, # grid::unit(1, "npc"), # - grid::unit(, "line"), #40,
                                                          frame.colour=c("NA"),
                                                          frame.linewidth=1))
        }
      }
        # scale_fill_gradient2(low = scale_grad[[v]]$low, 
        #                      high =scale_grad[[v]]$high, 
        #                      mid = scale_grad[[v]]$mid, 
        #                      midpoint = scale_grad[[v]]$midpoint,
        #                      na.value = NA,
        #                      name = paste0(name_var, ' (', type_estimate[v],
        #                                    '), ', name_model),
        #                      limits = limits[[v]]) +


      # ggsave(map_a_w_MAP, filename =paste0(dir_results, "/", folder_name, 
      #                                        "/map_a_w_MAP", ".png"),
      #        width = 15, height =15, dpi = 200)
      # map_a_w_MAP
      map_var_ALLPIX=c(map_var_ALLPIX,list(map_a_w_MAP))
      
      
      
      #^************************************************************************
      #^                       ONLY LAND PIXELS
      #^************************************************************************
      # Replace nan with 1:
      a_w_MAP_df$z[a_w_MAP_df$z==9999] =NA
      a_w_MAP_df_nonan <- a_w_MAP_df
      a_w_MAP_df_nonan$z[is.na(a_w_MAP_df_nonan$z)]=0
      # colMap <- colorRampPalette(c("yellow", "red"))(3)
      a_w_MAP_df_nonan_NEW = a_w_MAP_df_nonan
      # a_w_MAP_df_nonan_NEW = data.frame(x=round(a_w_MAP_df_nonan$x,3),
      #                                   y=round(a_w_MAP_df_nonan$y,3),
      #                                   z=round(a_w_MAP_df_nonan$z,3))
      rast_a_w_MAP_nonan =rasterFromXYZ(a_w_MAP_df_nonan_NEW, crs=crs_map)
      ## crop and mask (only LAND extent):
      r=rast_a_w_MAP_nonan
      r2 <- crop(r, extent(mask_shp))
      r3 <- mask(r2, mask_shp)
      test_spdf <- as(r3, "SpatialPixelsDataFrame")
      test_df <- as.data.frame(test_spdf)


      rast_signif_DIC3_LAND = rasterFromXYZ(signif_ns_DIC_df_3, crs=crs_map) 
      ## crop and mask (only LAND extent):
      r=rast_signif_DIC3_LAND
      r2 <- crop(r, extent(mask_shp))
      r3 <- mask(r2, mask_shp)
      test_spdf <- as(r3, "SpatialPixelsDataFrame")
      test_df_signif_DIC3 <- as.data.frame(test_spdf)
      
      # rast_signif_DIC1_LAND = rasterFromXYZ(signif_ns_DIC_df_1, crs=crs_map) 
      # ## crop and mask (only LAND extent):
      # r=rast_signif_DIC1_LAND
      # r2 <- crop(r, extent(land_shp))
      # r3 <- mask(r2, land_shp)
      # test_spdf <- as(r3, "SpatialPixelsDataFrame")
      # test_df_signif_DIC1 <- as.data.frame(test_spdf)
      
      map_a_w_MAP_LAND= ggplot() +
        # geom_sf(data = mask_world, fill = NA, color = "grey50")+
        theme_bw()+
        # scale_x_continuous(expand = c(0, 0), name = "Longitude (°E)")+
        # scale_y_continuous(expand = c(0, 0), name = "Latitude (°N)")+
        scale_x_continuous(expand = c(0, 0), name = "Longitude (°E)", 
                           breaks=c(6.00001, 10, 14, 18), 
                           labels=c("6°E", "10°E", "14°E", "18°E"))+
        scale_y_continuous(expand = c(0, 0), name = "Latitude (°N)",  
                           breaks=c(36.00001, 40, 44, 47.9999), 
                           labels=c("36°N", "40°N", "44°N", "48°N"))+
        geom_tile(data=na.omit(test_df), aes(x=x, y=y, fill=z), alpha=1) +
        #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
        geom_sf(data = mask_world, fill = NA, color = "grey70")+ #, linewidth=0.5) +
        geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+
        coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)+
        theme(
           panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank()
          ,axis.text.x      = element_text(size=13)
          ,axis.text.y      = element_text(size=13)
          ,legend.position  = legend.position
          ,legend.text      = element_text(size=15)
          ,legend.title     = element_text(size=20))

      if (flag_signif){
        # 2do: implement the significance from likelihood ratio test for fmins MLE method !!!
        map_a_w_MAP_LAND=map_a_w_MAP_LAND+
        #geom_point(data=signif_ns_DIC_df_1, aes(x=x, y=y, size=z), 
        #           fill="black", color="black", shape=21, size = 0.7)
        geom_point(data=test_df_signif_DIC3, aes(x=x, y=y, size=z),
                   fill="black", color="black", shape=21, size = 0.7)
        }
      if (flag.limits){
        if (scale_grad[[v]]$name=="turbo"){
          map_a_w_MAP_LAND=map_a_w_MAP_LAND+
          scale_fill_viridis(option="turbo",
                             limits = c(limits[[v]][1],limits[[v]][2]), 
                             name   = '', # paste0(name_var,' (',type_estimate[v],'), ',name_model),
                             oob    = scales::squish,
                             breaks = seq(limits[[v]][1], limits[[v]][2], limits[[v]][3]),
                             labels = seq(limits[[v]][1], limits[[v]][2], limits[[v]][3]),
                             guide  = guide_colourbar(nbin=nbin,
                                                      raster=T,
                                                      display = "rectangles",
                                                      direction=bar.direction,
                                                      barwidth=barwidth, #grid::unit(0.6, "npc"),
                                                      barheight=barheight, # grid::unit(1, "npc"), # - grid::unit(, "line"), #40,
                                                      frame.colour=c("NA"),
                                                      frame.linewidth=1))
        } else {
          map_a_w_MAP_LAND=map_a_w_MAP_LAND+
          scale_fill_gradient2(low = scale_grad[[v]]$low,
                               high =scale_grad[[v]]$high,
                               mid = scale_grad[[v]]$mid,
                               midpoint = scale_grad[[v]]$midpoint,
                               na.value = NA,
                               name   = '', # paste0(name_var,' (',type_estimate[v],'), ',name_model),
                               limits = c(limits[[v]][1],limits[[v]][2]), 
                               oob    = scales::squish,
                               breaks = seq(limits[[v]][1], limits[[v]][2], limits[[v]][3]),
                               labels = seq(limits[[v]][1], limits[[v]][2], limits[[v]][3]),
                               guide  = guide_colourbar(nbin=nbin,
                                                        raster=T,
                                                        display = "rectangles",
                                                        direction=bar.direction,
                                                        barwidth=barwidth, #grid::unit(0.6, "npc"),
                                                        barheight=barheight, # grid::unit(1, "npc"), # - grid::unit(, "line"), #40,
                                                        frame.colour=c("NA"),
                                                        frame.linewidth=1))
        }
      } else {
        if (scale_grad[[v]]$name=="turbo"){
          map_a_w_MAP_LAND=map_a_w_MAP_LAND+
            scale_fill_viridis(option="turbo",
                               name   = '', # paste0(name_var,' (',type_estimate[v],'), ',name_model),
                               guide  = guide_colourbar(nbin=nbin,
                                                        raster=T,
                                                        display = "rectangles",
                                                        direction=bar.direction,
                                                        barwidth=barwidth, #grid::unit(0.6, "npc"),
                                                        barheight=barheight, # grid::unit(1, "npc"), # - grid::unit(, "line"), #40,
                                                        frame.colour=c("NA"),
                                                        frame.linewidth=1))

        } else {
          map_a_w_MAP_LAND=map_a_w_MAP_LAND+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name   = '', # paste0(name_var,' (',type_estimate[v],'), ',name_model),
                                 guide  = guide_colourbar(nbin=nbin,
                                                          raster=T,
                                                          display = "rectangles",
                                                          direction=bar.direction,
                                                          barwidth=barwidth, #grid::unit(0.6, "npc"),
                                                          barheight=barheight, # grid::unit(1, "npc"), # - grid::unit(, "line"), #40,
                                                          frame.colour=c("NA"),
                                                          frame.linewidth=1))
        }
      }

      # theme(legend.position ="bottom")
      # ggsave(map_a_w_MAP_LAND, filename =paste0(dir_results, "/",  folder_name,
      #                                        "/map_a_w_MAP_LAND_", ".png"),
      #        width = 15, height =15, dpi = 200)
      # map_a_w_MAP_LAND

      map_var_LAND = c(map_var_LAND, list(map_a_w_MAP_LAND))
    }
  
    
    parent_dir= dirname(path_output[1])
    if (!dir.exists(parent_dir)){
      dir.create(parent_dir, recursive = T)
    }

    
    
    
    ##############
    # ALL PIXELS:
    ##############
    # merging the plots for different durations:
    plot_maps_gridplot(plotlist=map_var_ALLPIX, 
                       title= paste0(name_var," (", type_estimate[v], ") - ", name_model), 
                       pathout=paste0(path_output[v], '_', sat_prod, '_ALLPIX.png'))
    
    # merging the plots for different durations i horizontal with legend:
    legend_ALLPIX <- cowplot::get_plot_component(map_var_ALLPIX[[1]]+
                                          theme(legend.title = element_blank(),
                                                plot.margin = margin(t=0, r=0, b=0, l=0)),
                                          'guide-box-right', return_all=T)
    name2print= toupper(scan(text=name_var, what="")) # get the name of variable for title
    title_ALLPIX = cowplot::ggdraw() + draw_label(paste0(name2print[1], '\n ', name2print[2]), 
                                          fontface='bold', 
                                          x=0.5, hjust=0.5, vjust=0.5, size=26)
    # title_ALLPIX = plot_grid(cowplot::ggdraw() +
    #                             draw_label(paste0(name2print[1], '\n ', name2print[2]), 
    #                                         fontface='bold', 
    #                                         x=0.5, hjust=0.5, vjust=3, size=26),
    #                          # legend_ALLPIX,          
    #                          labels=NULL,
    #                          ncol =1, nrow=1,
    #                          align = 'h')
    
    all.plots_vardur_ALLPIX = plot_grid(
      plotlist=list(title_ALLPIX,
                    map_var_ALLPIX[[1]]+
                      ggtitle('1h')+ # this should be changed if other durations !!!!!!!!
                      theme(legend.title = element_blank(), 
                            #legend.position="none",
                            plot.margin = margin(t=0, r=0, b=10, l=0),
                            axis.text.x = element_text(size=12), # element_blank(),
                            axis.text.y = element_text(size=12),
                            axis.title  = element_blank(), 
                            plot.title  = element_text(face="bold", hjust = 0.5, size=26)),
                    map_var_ALLPIX[[2]] +
                      ggtitle('3h')+
                      theme(legend.title = element_blank(), 
                            # legend.position="none", 
                            axis.text.x = element_text(size=12), # element_blank(), 
                            axis.text.y = element_blank(),
                            plot.margin = margin(t=0, r=0, b=10, l=0),
                            axis.title  = element_blank(),
                            plot.title  = element_text(face="bold", hjust = 0.5, size=26)),
                    map_var_ALLPIX[[3]]+
                      ggtitle('6h')+
                      theme(legend.title = element_blank(), 
                            # legend.position="none",
                            plot.margin = margin(t=0, r=0, b=10, l=0),
                            axis.text.x = element_text(size=12), # element_blank(), 
                            axis.text.y = element_blank(),
                            axis.title  = element_blank(),
                            plot.title  = element_text(face="bold", hjust = 0.5, size=26)),
                    map_var_ALLPIX[[4]]+
                      ggtitle('24h')+
                      theme(legend.title = element_blank(), 
                            # legend.position="none",
                            plot.margin = margin(t=0, r=10, b=10, l=0),
                            axis.text.x = element_text(size=12), # element_blank(), 
                            axis.text.y = element_blank(),
                            axis.title  = element_blank(),
                            plot.title  = element_text(face="bold", hjust = 0.5, size=26))),
      label_size = 23,
      labels=NULL, # c("", "1h", "3h", "6h", "24h"),
      ncol =5, nrow=1, rel_widths = c(0.5,1,1,1,1),
      align = 'hv')
    # all.plots_vardur:
    ggsave(all.plots_vardur_ALLPIX, 
           filename=paste0(path_output[v], '_', sat_prod, '_horiz_ALLPIX.png'),
           width = 24, height =8, dpi = 200, bg = "white")
    
    
    
    
    
    ##############
    # ONLY LAND:
    ##############
    # merging the plots for different durations:
    plot_maps_gridplot(plotlist=map_var_LAND, 
                       title= paste0(name_var," (", type_estimate[v], ") - ", name_model), 
                       pathout=paste0(path_output[v], '_', sat_prod, '_LAND.png'))
    
    # merging the plots for different durations i horizontal with legend:
    legend_LAND <- cowplot::get_plot_component(map_var_LAND[[1]]+
                                            theme(legend.title = element_blank(),
                                                  plot.margin = margin(t=0, r=0, b=0, l=0)),
                                          'guide-box-right', return_all=T)
    title_LAND = cowplot::ggdraw() + draw_label(paste0(name2print[1], '\n ', name2print[2]), 
                                                  fontface='bold', 
                                                  x=0.5, hjust=0.5, vjust=0.5, size=26)
    # title_LAND = plot_grid(cowplot::ggdraw() +
    #                          draw_label(paste0(name2print[1], '\n ', name2print[2]), 
    #                                     fontface='bold', 
    #                                     x=0.5, hjust=0.5, vjust=3, size=26),
    #                        legend,          
    #                        labels=NULL,
    #                        ncol =1, nrow=2,
    #                        align = 'h')
    all.plots_vardur_LAND = plot_grid(
      plotlist=list(title_LAND,
                    map_var_LAND[[1]]+
                      ggtitle('1h')+
                      theme(legend.title = element_blank(), 
                            # legend.position="none",
                            plot.margin = margin(t=0, r=0, b=10, l=0),
                            axis.text.x = element_text(size=12), # element_blank(),
                            axis.text.y = element_text(size=12),
                            axis.title = element_blank(), 
                            plot.title = element_text(face="bold", hjust = 0.5, size=26)),
                    map_var_LAND[[2]]+
                      ggtitle('3h')+
                      theme(legend.title = element_blank(), 
                            # legend.position="none", 
                            axis.text.x = element_text(size=12), #element_blank(), 
                            axis.text.y = element_blank(),
                            plot.margin = margin(t=0, r=0, b=10, l=0),
                            axis.title  = element_blank(),
                            plot.title  = element_text(face="bold", hjust = 0.5, size=26)),
                    map_var_LAND[[3]]+
                      ggtitle('6h')+
                      theme(legend.title = element_blank(), 
                            # legend.position="none",
                            plot.margin = margin(t=0, r=0, b=10, l=0),
                            axis.text.x = element_text(size=12), # element_blank(), 
                            axis.text.y = element_blank(),
                            axis.title  = element_blank(),
                            plot.title  = element_text(face="bold", hjust = 0.5, size=26)),
                    map_var_LAND[[4]]+
                      ggtitle('24h')+
                      theme(legend.title = element_blank(), 
                            # legend.position="none",
                            plot.margin = margin(t=0, r=10, b=10, l=0),
                            axis.text.x = element_text(size=12), # element_blank(), 
                            axis.text.y = element_blank(),
                            axis.title  = element_blank(),
                            plot.title  = element_text(face="bold", hjust = 0.5, size=26))),
      label_size = 23,
      labels=NULL, # c("", "1h", "3h", "6h", "24h"),
      ncol =5, nrow=1, 
      rel_widths = c(0.5,1,1,1,1),
      align = 'hv')
    
    ggsave(all.plots_vardur_LAND, 
           filename=paste0(path_output[v], '_', sat_prod, '_horiz_LAND.png'),
           width = 24, height =8, dpi = 200, bg = "white")
    ############################################################################
    
    
    # save plots to return:
    map_var_LAND_allvar[[v]] = list() 
    map_var_LAND_allvar[[v]] = map_var_LAND
    map_var_ALLPIX_allvar[[v]] = list() 
    map_var_ALLPIX_allvar[[v]]  = map_var_ALLPIX
    
    map_var_LAND_allvar_horiz[[v]] = list() 
    map_var_LAND_allvar_horiz[[v]] = all.plots_vardur_LAND
    map_var_ALLPIX_allvar_horiz[[v]] = list() 
    map_var_ALLPIX_allvar_horiz[[v]]  = all.plots_vardur_ALLPIX
    
  } 
  


  return(list(map_var_LAND_allvar         = map_var_LAND_allvar,
              map_var_ALLPIX_allvar       = map_var_ALLPIX_allvar,
              map_var_LAND_allvar_horiz   = map_var_LAND_allvar_horiz,
              map_var_ALLPIX_allvar_horiz = map_var_ALLPIX_allvar_horiz))
}
















# ################################################################################
# plot_maps_param_paper = function(var_list,
#                                  lon_tmp,
#                                  lat_tmp,
#                                  durations,
#                                  mask_shp,
#                                  flag_signif=T,
#                                  significance,
#                                  name_var="Shape intercept",
#                                  type_estimate="MaxPost",
#                                  name_model = "SMEV Nonstationary",
#                                  path_output,
#                                  limits,
#                                  flag.limits=T,
#                                  scale_grad,
#                                  crs_map = crswgs84){
#   ################################################################################
#   # var = list() of matrices for each duration (e.g., a_w_MAP_ns)
#   
#   
#   for (v in 1:length(var_list)){
#     
#     # initialize variable plot list:
#     map_var_EU = list()
#     map_var_ITA = list()
#     a_w_MAP_df=c()
#     
#     
#     # durations
#     for (dur in 1:length(durations)){
#       
#       ############
#       # MAP
#       ############
#       # create dataframe for the map plot:
#       a_w_MAP_df <- data.frame(x=lon_tmp, y=lat_tmp, z=var_list[[v]][[dur]])
#       # trend_df$z[is.na(trend_df$z)]=0
#       a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
#                                   y=round(a_w_MAP_df$y,3), 
#                                   z=round(a_w_MAP_df$z,3))
#       rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
#       
#       # increasing AMAX
#       signif_AMAX_pos_df <- data.frame(x=lon_tmp, y=lat_tmp, 
#                                        z=significance$signif_MK_AMAX[[dur]], 
#                                        d=significance$sensslope_AMAX[[dur]])
#       signif_AMAX_pos_df$z[signif_AMAX_pos_df$z==F]=NA
#       signif_AMAX_pos_df$z[signif_AMAX_pos_df$d<0] =NA
#       signif_AMAX_pos_df = na.omit(signif_AMAX_pos_df)
#       signif_AMAX_pos_df = subset(signif_AMAX_pos_df, z != 0)  
#       
#       # decreasing AMAX:
#       signif_AMAX_neg_df <- data.frame(x=lon_tmp, y=lat_tmp, 
#                                        z=significance$signif_MK_AMAX[[dur]], 
#                                        d=significance$sensslope_AMAX[[dur]])
#       signif_AMAX_neg_df$z[signif_AMAX_neg_df$z==F]=NA
#       signif_AMAX_neg_df$z[signif_AMAX_neg_df$d>=0]=NA
#       signif_AMAX_neg_df = na.omit(signif_AMAX_neg_df)
#       signif_AMAX_neg_df = subset(signif_AMAX_neg_df, z != 0)  
#       
#       # nonstationary model is more likely than stationary (DIC)
#       signif_ns_DIC_df_3 <- data.frame(x=lon_tmp, y=lat_tmp, z=significance$DIC_diff3[[dur]])
#       signif_ns_DIC_df_3$z[signif_ns_DIC_df_3$z==F]=NA
#       signif_ns_DIC_df_3 = na.omit(signif_ns_DIC_df_3)
#       signif_ns_DIC_df_3 = subset(signif_ns_DIC_df_3, z != 0)      
#       
#       signif_ns_DIC_df_1 <- data.frame(x=lon_tmp, y=lat_tmp, z=significance$DIC_diff1[[dur]])
#       signif_ns_DIC_df_1$z[signif_ns_DIC_df_1$z==F]=NA
#       signif_ns_DIC_df_1 = na.omit(signif_ns_DIC_df_1)
#       signif_ns_DIC_df_1 = subset(signif_ns_DIC_df_1, z != 0)  
#       
#       
#       
#       # ALL TILE 
#       #####################################
#       map_a_w_MAP = ggplot() + 
#         theme_bw()+
#         scale_x_continuous(expand = c(0, 0), name = "Longitude (°E)")+
#         scale_y_continuous(expand = c(0, 0), name = "Latitude (°N)")+
#         geom_tile(data=a_w_MAP_df_NEW, aes(x=x, y=y, fill=z), alpha=0.7) + 
#         #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
#         geom_sf(data = gg_world, fill = NA, color = "grey70")+ #, linewidth=0.5) +
#         geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+ #, linewidth=0.5)+
#         coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
#       
#       if (flag_signif){
#         map_a_w_MAP=map_a_w_MAP+
#           # geom_point(data=signif_AMAX_pos_df, aes(x=x, y=y, size=d), fill="blue", color="blue", shape=24)+
#           # geom_point(data=signif_AMAX_neg_df, aes(x=x, y=y, size=d), fill="red", color="red",  shape=25)+
#           #geom_point(data=signif_ns_DIC_df_1, aes(x=x, y=y, size=z),   fill="black", color="black", shape=21, size = 0.7)
#           geom_point(data=signif_ns_DIC_df_3, aes(x=x, y=y, size=z),   fill="black", color="black", shape=21, size = 0.7)
#         #geom_point(data=signif_ns_DIC_df, aes(x=x, y=y, size=z), shape = 1, size = 1)
#         
#       }
#       map_a_w_MAP=map_a_w_MAP+
#         theme(
#           # legend.position = "none"
#           ,legend.text     = element_text(size=15)
#           ,legend.title    = element_text(size=20))
#       
#       
#       if (flag.limits){
#         if (scale_grad[[v]]$name=="turbo"){
#           map_a_w_MAP=map_a_w_MAP+
#             scale_fill_viridis(option="turbo",
#                                limits = limits[[v]],
#                                name = paste0(name_var, ' (', type_estimate[v],
#                                              '), ', name_model),
#                                oob    = scales::squish)
#           
#         } else {
#           map_a_w_MAP=map_a_w_MAP+
#             scale_fill_gradient2(low = scale_grad[[v]]$low,
#                                  high =scale_grad[[v]]$high,
#                                  mid = scale_grad[[v]]$mid,
#                                  midpoint = scale_grad[[v]]$midpoint,
#                                  na.value = NA,
#                                  name = paste0(name_var, '(',
#                                                type_estimate[v], '), ', name_model),
#                                  limits = limits[[v]],
#                                  oob    = scales::squish)
#           
#         }
#         
#       } else {
#         if (scale_grad[[v]]$name=="turbo"){
#           
#           map_a_w_MAP=map_a_w_MAP+
#             scale_fill_viridis(option="turbo",
#                                name = paste0(name_var, ' (', type_estimate[v],
#                                              '), ', name_model))
#         } else {
#           map_a_w_MAP=map_a_w_MAP+
#             scale_fill_gradient2(low = scale_grad[[v]]$low,
#                                  high =scale_grad[[v]]$high,
#                                  mid = scale_grad[[v]]$mid,
#                                  midpoint = scale_grad[[v]]$midpoint,
#                                  na.value = NA,
#                                  name = paste0(name_var, '(',
#                                                type_estimate[v], '), ', name_model))
#         }
#       }
#       # scale_fill_gradient2(low = scale_grad[[v]]$low, 
#       #                      high =scale_grad[[v]]$high, 
#       #                      mid = scale_grad[[v]]$mid, 
#       #                      midpoint = scale_grad[[v]]$midpoint,
#       #                      na.value = NA,
#       #                      name = paste0(name_var, ' (', type_estimate[v],
#       #                                    '), ', name_model),
#       #                      limits = limits[[v]]) +
#       
#       
#       # ggsave(map_a_w_MAP, filename =paste0(dir_results, "/", folder_name, 
#       #                                        "/map_a_w_MAP", ".png"),
#       #        width = 15, height =15, dpi = 200)
#       # map_a_w_MAP
#       map_var_EU = c(map_var_EU, list(map_a_w_MAP))
#       
#       
#       
#       
#       
#       
#       
#       # ITALY (pixel level NOT AGGREGATED):  
#       ####################################
#       # Replace nan with 1:
#       a_w_MAP_df$z[a_w_MAP_df$z==9999] =NA
#       a_w_MAP_df_nonan <- a_w_MAP_df
#       a_w_MAP_df_nonan$z[is.na(a_w_MAP_df_nonan$z)]=0
#       # colMap <- colorRampPalette(c("yellow", "red"))(3)
#       a_w_MAP_df_nonan_NEW = data.frame(x=round(a_w_MAP_df_nonan$x,3),
#                                         y=round(a_w_MAP_df_nonan$y,3), 
#                                         z=round(a_w_MAP_df_nonan$z,3))
#       rast_a_w_MAP_nonan =rasterFromXYZ(a_w_MAP_df_nonan_NEW, crs=crs_map) 
#       ## crop and mask (only Italy extent):
#       r=rast_a_w_MAP_nonan
#       r2 <- crop(r, extent(mask_shp))
#       r3 <- mask(r2, mask_shp)
#       test_spdf <- as(r3, "SpatialPixelsDataFrame")
#       test_df <- as.data.frame(test_spdf)
#       
#       map_a_w_MAP_ITALY= ggplot() + 
#         # geom_sf(data = gg_world, fill = NA, color = "grey50")+
#         theme_bw()+
#         scale_x_continuous(expand = c(0, 0), name = "Longitude (°E)")+
#         scale_y_continuous(expand = c(0, 0), name = "Latitude (°N)")+
#         geom_tile(data=na.omit(test_df), aes(x=x, y=y, fill=z), alpha=1) + 
#         #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
#         geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+
#         coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)+
#         theme(
#           # legend.position = "none"
#           ,legend.text     = element_text(size=15)
#           ,legend.title    = element_text(size=20))
#       
#       
#       if (flag.limits){
#         if (scale_grad[[v]]$name=="turbo"){
#           map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
#             scale_fill_viridis(option="turbo",
#                                limits = limits[[v]],
#                                name = paste0(name_var, ' (', type_estimate[v],
#                                              '), ', name_model),
#                                oob    = scales::squish)
#         } else {
#           map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
#             scale_fill_gradient2(low = scale_grad[[v]]$low,
#                                  high =scale_grad[[v]]$high,
#                                  mid = scale_grad[[v]]$mid,
#                                  midpoint = scale_grad[[v]]$midpoint,
#                                  na.value = NA,
#                                  name = paste0(name_var, '(',
#                                                type_estimate[v], '), ', name_model),
#                                  limits = limits[[v]],
#                                  oob    = scales::squish)
#         }
#         
#       } else {
#         if (scale_grad[[v]]$name=="turbo"){
#           
#           map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
#             scale_fill_viridis(option="turbo",
#                                name = paste0(name_var, ' (', type_estimate[v],
#                                              '), ', name_model))
#         } else {
#           map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
#             scale_fill_gradient2(low = scale_grad[[v]]$low,
#                                  high =scale_grad[[v]]$high,
#                                  mid = scale_grad[[v]]$mid,
#                                  midpoint = scale_grad[[v]]$midpoint,
#                                  na.value = NA,
#                                  name = paste0(name_var, '(',
#                                                type_estimate[v], '), ', name_model))
#         }
#       }
#       
#       # ggsave(map_a_w_MAP_ITALY, filename =paste0(dir_results, "/",  folder_name, 
#       #                                        "/map_a_w_MAP_noEU_", ".png"),
#       #        width = 15, height =15, dpi = 200)
#       # map_a_w_MAP_ITALY
#       
#       map_var_ITA = c(map_var_ITA, list(map_a_w_MAP_ITALY))
#       
#       
#     }
#     
#     parent_dir= dirname(path_output[1])
#     if (!dir.exists(parent_dir)){
#       dir.create(parent_dir, recursive = T)
#     }
#     
#     
#     
#     
#     # merging the plots for different durations:
#     plot_maps_gridplot(plotlist=map_var_EU,
#                        title= paste0(name_var," (", type_estimate[v], ") - ", name_model),
#                        pathout=paste0(path_output[v], '_EU.png'))
#     plot_maps_gridplot(plotlist=map_var_ITA, 
#                        title= paste0(name_var," (", type_estimate[v], ") - ", name_model), 
#                        pathout=paste0(path_output[v], '_ITA.png'))
#     
#     
#     
#     # Plot only 1h and 24 hours durations for paper/presentations:
#     ##############################################################
#     # merging the plots (cutting space and with labels on the top):
#     title1 <- ggdraw() + 
#       draw_label(paste0(name_var," (", type_estimate[v], ") - ", name_model),
#                  fontface = 'bold', x = 0.5,   hjust = 0.5, vjust = 0.5, size=20)
#     # theme(
#     #   # add margin on the left of the drawing canvas,
#     #   # so title is aligned with left edge of first plot
#     #   plot.margin = margin(0, 0, 0, 7)
#     # )
#     
#     all.plots = plot_grid( map_var_EU[[1]] + theme( legend.title = element_blank()
#                                                    ,legend.position="none"
#                                                    ,legend.margin = margin(0, 0, 0, l=0) # turned off for alignment
#                                                    ,axis.title.y =element_blank() #element_text(size = 15)
#                                                    ,axis.text.y = element_text(size=13)
#                                                    ,axis.title.x = element_blank() # element_text(size = 15)
#                                                    ,axis.text.x = element_text(size=13)
#                                                    ,plot.margin = margin(t=15, r=20, b=0, l=5)),
#                            
#                            map_var_EU[[4]] + theme( legend.title = element_blank()
#                                                    ,axis.title.y = element_blank()
#                                                    ,axis.text.y = element_blank()
#                                                    ,axis.title.x = element_blank()
#                                                    ,axis.text.x = element_text(size=13)
#                                                    ,plot.margin = margin(t=15, r=-5, b=0, l=-10)),
#                            label_size = 25,
#                            labels = c("    1h", "24h"),
#                            nrow=1, 
#                            hjust = -4,
#                            rel_heights = c(1, 1),
#                            rel_widths = c(1, 1),
#                            align = 'h') 
#     all.plots
#     all.plots.title = plot_grid( title1,
#                                  all.plots,
#                                  ncol =1, nrow=2,
#                                  rel_heights = c(0.05, 1))
#     ggsave(all.plots.title, filename=paste0(path_output[v], '_EU_3.png'), width = 14, height =9, dpi = 150, bg = "white")
#   
#   }  
# }
# 









  


################################################################################
plot_maps_change = function(var_int_list,
                            var_slope_list,
                            year_ref      = 2011,
                            year_ref_id   = 13,
                            lon_tmp,
                            lat_tmp,
                            durations,
                            mask_shp,
                            flag_signif   = T, # logical
                            significance  = significance,
                            name_var="Shape intercept",
                            type_estimate="MaxPost",
                            name_model = "SMEV Nonstationary",
                            path_output,
                            limits,
                            flag.limits=T,
                            scale_grad,
                            crs_map){
################################################################################
  scale_change = list()
  
  
  for (v in 1:length(var_int_list)){
    
    # initialize variable plot list:
    map_var_EU = list()
    map_var_ITA = list()
    a_w_MAP_df=c()
    scale_change[[v]] =list()
    
    # durations
    for (dur in 1:length(durations)){

      # % of change of slope per decade wrw halfperiod year (a sort of average change per decade):
      scale_change[[v]][[dur]] = var_slope_list[[v]][[dur]]*10 / (var_int_list[[v]][[dur]] +
                                                                  var_slope_list[[v]][[dur]]*year_ref_id)*100
                    
      # create dataframe for the map plot:
      a_w_MAP_df <- data.frame(x=lon_tmp, y=lat_tmp, z=scale_change[[v]][[dur]])
      
      # trend_df$z[is.na(trend_df$z)]=0
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3), 
                                  z=round(a_w_MAP_df$z,3))
      rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
      
      
      # # Find pixels with significant nonstationarity:
      # trend_df <- data.frame(x=lon_tmp, y=lat_tmp, z=signif_trend)
      # trend_df_NEW = data.frame(x=round(trend_df$x,3),
      #                           y=round(trend_df$y,3), 
      #                           z=round(trend_df$z,3))
      # df_NEW_SIGNIF_EU = trend_df_NEW
      # df_NEW_SIGNIF_EU[df_NEW_SIGNIF_EU==9999] =NA
      # df_NEW_SIGNIF_EU = na.omit(df_NEW_SIGNIF_EU)
      # df_NEW_SIGNIF_EU = subset(df_NEW_SIGNIF_EU, z != 0)  
      
      # increasing AMAX
      signif_AMAX_pos_df <- data.frame(x=lon_tmp, y=lat_tmp,
                                       z=significance$signif_MK_AMAX[[dur]], 
                                       d=significance$sensslope_AMAX[[dur]])
      signif_AMAX_pos_df$z[signif_AMAX_pos_df$z==F]=NA
      signif_AMAX_pos_df$z[signif_AMAX_pos_df$d<0] =NA
      signif_AMAX_pos_df = na.omit(signif_AMAX_pos_df)
      signif_AMAX_pos_df = subset(signif_AMAX_pos_df, z != 0)  
      
      # decreasing AMAX:
      signif_AMAX_neg_df <- data.frame(x=lon_tmp, y=lat_tmp, 
                                       z=significance$signif_MK_AMAX[[dur]],
                                       d=significance$sensslope_AMAX[[dur]])
      signif_AMAX_neg_df$z[signif_AMAX_neg_df$z==F]=NA
      signif_AMAX_neg_df$z[signif_AMAX_neg_df$d>=0]=NA
      signif_AMAX_neg_df = na.omit(signif_AMAX_neg_df)
      signif_AMAX_neg_df = subset(signif_AMAX_neg_df, z != 0)  
      
      # nonstationary model is more likely than stationary (DIC)
      signif_ns_DIC_df_3 <- data.frame(x=lon_tmp, y=lat_tmp, z=significance$DIC_diff3[[dur]])
      signif_ns_DIC_df_3$z[signif_ns_DIC_df_3$z==F]=NA
      signif_ns_DIC_df_3 = na.omit(signif_ns_DIC_df_3)
      signif_ns_DIC_df_3 = subset(signif_ns_DIC_df_3, z != 0)      
      
      signif_ns_DIC_df_1 <- data.frame(x=lon_tmp, y=lat_tmp, z=significance$DIC_diff1[[dur]])
      signif_ns_DIC_df_1$z[signif_ns_DIC_df_1$z==F]=NA
      signif_ns_DIC_df_1 = na.omit(signif_ns_DIC_df_1)
      signif_ns_DIC_df_1 = subset(signif_ns_DIC_df_1, z != 0)  
      
      
      # WHOLE TILE 
      #####################################
      map_a_w_MAP = ggplot() + 
        theme_bw()+
        scale_x_continuous(expand = c(0, 0), name = "Longitude (°E)")+
        scale_y_continuous(expand = c(0, 0), name = "Latitude (°N)")+
        geom_tile(data=a_w_MAP_df_NEW, aes(x=x, y=y, fill=z), alpha=1) + 
        #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
        geom_sf(data = gg_world, fill = NA, color = "grey70") +
        geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+
        
        # geom_point(data =df_NEW_SIGNIF_EU, aes(x=x, y=y, size= z), 
        #            fill="black", shape = 19, size = 1) +  # pixels with signif. nonstationarity
        
        coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
      if (flag_signif){
        map_a_w_MAP=map_a_w_MAP+
          geom_point(data=signif_AMAX_pos_df, aes(x=x, y=y, size=d), fill="blue", color="blue", shape=24)+
          geom_point(data=signif_AMAX_neg_df, aes(x=x, y=y, size=d), fill="red", color="red",  shape=25)+
          geom_point(data=signif_ns_DIC_df_1, aes(x=x, y=y, size=z),   fill="gray90", color="gray90", shape=21, size = 0.7)+
          geom_point(data=signif_ns_DIC_df_3, aes(x=x, y=y, size=z),   fill="black", color="black", shape=21, size = 0.7)
        #geom_point(data=signif_ns_DIC_df, aes(x=x, y=y, size=z), shape = 1, size = 1)
        
      }
      map_a_w_MAP=map_a_w_MAP+
        # scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0,
        #                      na.value = NA, name = paste0(name_var, ' (', type_estimate[v], '), ', name_model))+
        theme(
          # legend.position = "none"
          ,legend.text     = element_text(size=15)
          ,legend.title    = element_text(size=20))
      if (flag.limits){
        if (scale_grad[[v]]$name=="turbo"){
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_viridis(option="turbo",
                               limits = limits[[v]],
                               name = paste0(name_var, ' (', type_estimate[v],
                                             '), ', name_model),
                               oob    = scales::squish)
          
        } else {
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name = paste0(name_var, '(',
                                               type_estimate[v], '), ', name_model),
                                 limits = limits[[v]],
                                 oob    = scales::squish)
          
        }
        
      } else {
        if (scale_grad[[v]]$name=="turbo"){
          
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_viridis(option="turbo",
                               name = paste0(name_var, ' (', type_estimate[v],
                                             '), ', name_model))
        } else {
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name = paste0(name_var, '(',
                                               type_estimate[v], '), ', name_model))
        }
      }
      # ggsave(map_a_w_MAP, filename =paste0(dir_results, "/", folder_name, 
      #                                        "/map_a_w_MAP", ".png"),
      #        width = 15, height =15, dpi = 200)
      # map_a_w_MAP
      map_var_EU = c(map_var_EU, list(map_a_w_MAP))
      
      
      
      
      
      
      # ITALY (pixel level NOT AGGREGATED):  
      #####################################
      # Replace nan with 1:
      a_w_MAP_df$z[a_w_MAP_df$z==9999] =NA
      a_w_MAP_df_nonan <- a_w_MAP_df
      a_w_MAP_df_nonan$z[is.na(a_w_MAP_df_nonan$z)]=0
      # colMap <- colorRampPalette(c("yellow", "red"))(3)
      a_w_MAP_df_nonan_NEW = data.frame(x=round(a_w_MAP_df_nonan$x,3),
                                        y=round(a_w_MAP_df_nonan$y,3), 
                                        z=round(a_w_MAP_df_nonan$z,3))
      rast_a_w_MAP_nonan =rasterFromXYZ(a_w_MAP_df_nonan_NEW, crs=crs_map) 
      ## crop and mask (only Italy extent):
      r=rast_a_w_MAP_nonan
      r2 <- crop(r, extent(mask_shp))
      r3 <- mask(r2, mask_shp)
      test_spdf <- as(r3, "SpatialPixelsDataFrame")
      test_df <- as.data.frame(test_spdf)
      
      map_a_w_MAP_ITALY= ggplot() + 
        # geom_sf(data = gg_world, fill = NA, color = "grey50")+
        theme_bw()+
        scale_x_continuous(expand = c(0, 0), name = "Longitude (°E)")+
        scale_y_continuous(expand = c(0, 0), name = "Latitude (°N)")+
        geom_tile(data=na.omit(test_df), aes(x=x, y=y, fill=z), alpha=1) + 
        #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
        geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+
        coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)+
        # scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0,
        #                      na.value = NA, name = paste0(name_var, '(', type_estimate[v], '), ', name_model))+
        theme(
          # legend.position = "none"
          ,legend.text     = element_text(size=15)
          ,legend.title    = element_text(size=20))
      if (flag.limits){
        if (scale_grad[[v]]$name=="turbo"){
          map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
            scale_fill_viridis(option="turbo",
                               limits = limits[[v]],
                               name = paste0(name_var, ' (', type_estimate[v],
                                             '), ', name_model),
                               oob    = scales::squish)
        } else {
          map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name = paste0(name_var, '(',
                                               type_estimate[v], '), ', name_model),
                                 limits = limits[[v]],
                                 oob    = scales::squish)
        }
        
      } else {
        if (scale_grad[[v]]$name=="turbo"){
          
          map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
            scale_fill_viridis(option="turbo",
                               name = paste0(name_var, ' (', type_estimate[v],
                                             '), ', name_model))
        } else {
          map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name = paste0(name_var, '(',
                                               type_estimate[v], '), ', name_model))
        }
      }
      
      # ggsave(map_a_w_MAP_ITALY, filename =paste0(dir_results, "/",  folder_name, 
      #                                        "/map_a_w_MAP_noEU_", ".png"),
      #        width = 15, height =15, dpi = 200)
      # map_a_w_MAP_ITALY
      
      map_var_ITA = c(map_var_ITA, list(map_a_w_MAP_ITALY))
      
    }
    
    parent_dir= dirname(path_output[1])
    if (!dir.exists(parent_dir)){
      dir.create(parent_dir, recursive = T)
    }
    
    # merging the plots for different durations:
    plot_maps_gridplot(plotlist=map_var_EU, 
                       title= paste0(name_var," (", type_estimate[v], ") - ", name_model), 
                       pathout=paste0(path_output[v], '_EU.png'))
    
    plot_maps_gridplot(plotlist=map_var_ITA, 
                       title= paste0(name_var," (", type_estimate[v], ") - ", name_model), 
                       pathout=paste0(path_output[v], '_ITA.png'))
  }  
  
}














################################################################################
plot_maps_change_paper = function(var_int_list,
                                  var_slope_list,
                                  year_ref      = 2011,
                                  year_ref_id   = 13,
                                  lon_tmp,
                                  lat_tmp,
                                  durations,
                                  mask_shp,
                                  flag_signif   = T, # logical
                                  significance  = significance,
                                  name_var="Shape intercept",
                                  type_estimate="MaxPost",
                                  name_model = "SMEV Nonstationary",
                                  path_output,
                                  limits,
                                  flag.limits=T,
                                  scale_grad,
                                  crs_map){
  ################################################################################
  scale_change = list()
  
  
  for (v in 1:length(var_int_list)){
    
    # initialize variable plot list:
    map_var_EU = list()
    map_var_ITA = list()
    a_w_MAP_df=c()
    scale_change[[v]] =list()
    
    # durations
    for (dur in 1:length(durations)){
      
      # % of change of slope per decade wrw halfperiod year (a sort of average change per decade):
      scale_change[[v]][[dur]] = var_slope_list[[v]][[dur]]*10 / (var_int_list[[v]][[dur]] +
                                                                    var_slope_list[[v]][[dur]]*year_ref_id)*100
      
      # create dataframe for the map plot:
      a_w_MAP_df <- data.frame(x=lon_tmp, y=lat_tmp, z=scale_change[[v]][[dur]])
      
      # trend_df$z[is.na(trend_df$z)]=0
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3), 
                                  z=round(a_w_MAP_df$z,3))
      rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
      
      
      # # Find pixels with significant nonstationarity:
      # trend_df <- data.frame(x=lon_tmp, y=lat_tmp, z=signif_trend)
      # trend_df_NEW = data.frame(x=round(trend_df$x,3),
      #                           y=round(trend_df$y,3), 
      #                           z=round(trend_df$z,3))
      # df_NEW_SIGNIF_EU = trend_df_NEW
      # df_NEW_SIGNIF_EU[df_NEW_SIGNIF_EU==9999] =NA
      # df_NEW_SIGNIF_EU = na.omit(df_NEW_SIGNIF_EU)
      # df_NEW_SIGNIF_EU = subset(df_NEW_SIGNIF_EU, z != 0)  
      
      # increasing AMAX
      signif_AMAX_pos_df <- data.frame(x=lon_tmp, y=lat_tmp,
                                       z=significance$signif_MK_AMAX[[dur]], 
                                       d=significance$sensslope_AMAX[[dur]])
      signif_AMAX_pos_df$z[signif_AMAX_pos_df$z==F]=NA
      signif_AMAX_pos_df$z[signif_AMAX_pos_df$d<0] =NA
      signif_AMAX_pos_df = na.omit(signif_AMAX_pos_df)
      signif_AMAX_pos_df = subset(signif_AMAX_pos_df, z != 0)  
      
      # decreasing AMAX:
      signif_AMAX_neg_df <- data.frame(x=lon_tmp, y=lat_tmp, 
                                       z=significance$signif_MK_AMAX[[dur]],
                                       d=significance$sensslope_AMAX[[dur]])
      signif_AMAX_neg_df$z[signif_AMAX_neg_df$z==F]=NA
      signif_AMAX_neg_df$z[signif_AMAX_neg_df$d>=0]=NA
      signif_AMAX_neg_df = na.omit(signif_AMAX_neg_df)
      signif_AMAX_neg_df = subset(signif_AMAX_neg_df, z != 0)  
      
      # nonstationary model is more likely than stationary (DIC)
      signif_ns_DIC_df_3 <- data.frame(x=lon_tmp, y=lat_tmp, z=significance$DIC_diff3[[dur]])
      signif_ns_DIC_df_3$z[signif_ns_DIC_df_3$z==F]=NA
      signif_ns_DIC_df_3 = na.omit(signif_ns_DIC_df_3)
      signif_ns_DIC_df_3 = subset(signif_ns_DIC_df_3, z != 0)      
      
      signif_ns_DIC_df_1 <- data.frame(x=lon_tmp, y=lat_tmp, z=significance$DIC_diff1[[dur]])
      signif_ns_DIC_df_1$z[signif_ns_DIC_df_1$z==F]=NA
      signif_ns_DIC_df_1 = na.omit(signif_ns_DIC_df_1)
      signif_ns_DIC_df_1 = subset(signif_ns_DIC_df_1, z != 0)  
      
      
      # WHOLE TILE 
      #####################################
      map_a_w_MAP = ggplot() + 
        theme_bw()+
        scale_x_continuous(expand = c(0, 0), name = "Longitude (°E)")+
        scale_y_continuous(expand = c(0, 0), name = "Latitude (°N)")+
        geom_tile(data=a_w_MAP_df_NEW, aes(x=x, y=y, fill=z), alpha=0.7) + 
        #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
        geom_sf(data = gg_world, fill = NA, color = "grey70",  linewidth=0.5) +
        geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2, linewidth=0.5)+
        
        # geom_point(data =df_NEW_SIGNIF_EU, aes(x=x, y=y, size= z), 
        #            fill="black", shape = 19, size = 1) +  # pixels with signif. nonstationarity
        
        coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
      if (flag_signif){
        map_a_w_MAP=map_a_w_MAP+
          # geom_point(data=signif_AMAX_pos_df, aes(x=x, y=y, size=d), fill="blue", color="blue", shape=24)+
          # geom_point(data=signif_AMAX_neg_df, aes(x=x, y=y, size=d), fill="red", color="red",  shape=25)+
          #geom_point(data=signif_ns_DIC_df_1, aes(x=x, y=y, size=z),   fill="gray90", color="gray90", shape=21, size = 0.7)+
          geom_point(data=signif_ns_DIC_df_3, aes(x=x, y=y, size=z),   fill="black", color="black", shape=21, size = 0.7)
          # geom_point(data=signif_ns_DIC_df, aes(x=x, y=y, size=z), shape = 1, size = 1)
        
      }
      map_a_w_MAP=map_a_w_MAP+
        # scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0,
        #                      na.value = NA, name = paste0(name_var, ' (', type_estimate[v], '), ', name_model))+
        theme(
          # legend.position = "none"
          ,legend.text     = element_text(size=15)
          ,legend.title    = element_text(size=20))
      if (flag.limits){
        if (scale_grad[[v]]$name=="turbo"){
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_viridis(option="turbo",
                               limits = limits[[v]],
                               name = paste0(name_var, ' (', type_estimate[v],
                                             '), ', name_model),
                               oob    = scales::squish)
          
        } else {
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name = paste0(name_var, '(',
                                               type_estimate[v], '), ', name_model),
                                 limits = limits[[v]],
                                 oob    = scales::squish)
          
        }
        
      } else {
        if (scale_grad[[v]]$name=="turbo"){
          
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_viridis(option="turbo",
                               name = paste0(name_var, ' (', type_estimate[v],
                                             '), ', name_model))
        } else {
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name = paste0(name_var, '(',
                                               type_estimate[v], '), ', name_model))
        }
      }
      # ggsave(map_a_w_MAP, filename =paste0(dir_results, "/", folder_name, 
      #                                        "/map_a_w_MAP", ".png"),
      #        width = 15, height =15, dpi = 200)
      # map_a_w_MAP
      map_var_EU = c(map_var_EU, list(map_a_w_MAP))
      
      
      
      
      
      
      # ITALY (pixel level NOT AGGREGATED):  
      #####################################
      # Replace nan with 1:
      a_w_MAP_df$z[a_w_MAP_df$z==9999] =NA
      a_w_MAP_df_nonan <- a_w_MAP_df
      a_w_MAP_df_nonan$z[is.na(a_w_MAP_df_nonan$z)]=0
      # colMap <- colorRampPalette(c("yellow", "red"))(3)
      a_w_MAP_df_nonan_NEW = data.frame(x=round(a_w_MAP_df_nonan$x,3),
                                        y=round(a_w_MAP_df_nonan$y,3), 
                                        z=round(a_w_MAP_df_nonan$z,3))
      rast_a_w_MAP_nonan =rasterFromXYZ(a_w_MAP_df_nonan_NEW, crs=crs_map) 
      ## crop and mask (only Italy extent):
      r=rast_a_w_MAP_nonan
      r2 <- crop(r, extent(mask_shp))
      r3 <- mask(r2, mask_shp)
      test_spdf <- as(r3, "SpatialPixelsDataFrame")
      test_df <- as.data.frame(test_spdf)
      
      map_a_w_MAP_ITALY= ggplot() + 
        # geom_sf(data = gg_world, fill = NA, color = "grey50")+
        theme_bw()+
        scale_x_continuous(expand = c(0, 0), name = "Longitude (°E)")+
        scale_y_continuous(expand = c(0, 0), name = "Latitude (°N)")+
        geom_tile(data=na.omit(test_df), aes(x=x, y=y, fill=z), alpha=1) + 
        #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
        geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+
        coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)+
        # scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0,
        #                      na.value = NA, name = paste0(name_var, '(', type_estimate[v], '), ', name_model))+
        theme(
          # legend.position = "none"
          ,legend.text     = element_text(size=15)
          ,legend.title    = element_text(size=20))
      if (flag.limits){
        if (scale_grad[[v]]$name=="turbo"){
          map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
            scale_fill_viridis(option="turbo",
                               limits = limits[[v]],
                               name = paste0(name_var, ' (', type_estimate[v],
                                             '), ', name_model),
                               oob    = scales::squish)
        } else {
          map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name = paste0(name_var, '(',
                                               type_estimate[v], '), ', name_model),
                                 limits = limits[[v]],
                                 oob    = scales::squish)
        }
        
      } else {
        if (scale_grad[[v]]$name=="turbo"){
          
          map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
            scale_fill_viridis(option="turbo",
                               name = paste0(name_var, ' (', type_estimate[v],
                                             '), ', name_model))
        } else {
          map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name = paste0(name_var, '(',
                                               type_estimate[v], '), ', name_model))
        }
      }
      
      # ggsave(map_a_w_MAP_ITALY, filename =paste0(dir_results, "/",  folder_name, 
      #                                        "/map_a_w_MAP_noEU_", ".png"),
      #        width = 15, height =15, dpi = 200)
      # map_a_w_MAP_ITALY
      
      map_var_ITA = c(map_var_ITA, list(map_a_w_MAP_ITALY))
      
    }
    
    parent_dir= dirname(path_output[1])
    if (!dir.exists(parent_dir)){
      dir.create(parent_dir, recursive = T)
    }
    
    # merging the plots for different durations:
    plot_maps_gridplot(plotlist=map_var_EU, 
                       title= paste0(name_var," (", type_estimate[v], ") - ", name_model), 
                       pathout=paste0(path_output[v], '_EU.png'))
    
    plot_maps_gridplot(plotlist=map_var_ITA, 
                       title= paste0(name_var," (", type_estimate[v], ") - ", name_model), 
                       pathout=paste0(path_output[v], '_ITA.png'))
    
    
    
    # Plot only 1h and 24 hours durations for paper/presentations:
    ##############################################################
    # merging the plots (cutting space and with labels on the top):
    title1 <- ggdraw() + 
      draw_label(paste0(name_var," (", type_estimate[v], ") - ", name_model),
                 fontface = 'bold', x = 0.5,   hjust = 0.5, vjust = 0.5, size=20)
    # theme(
    #   plot.margin = margin(0, 0, 0, 7)
    # )
    
    all.plots = plot_grid( map_var_EU[[1]] + theme( legend.title = element_blank()
                                                    ,legend.position="none"
                                                    ,legend.margin = margin(0, 0, 0, l=0) # turned off for alignment
                                                    ,axis.title.y =element_blank() #element_text(size = 15)
                                                    ,axis.text.y = element_text(size=13)
                                                    ,axis.title.x = element_blank() # element_text(size = 15)
                                                    ,axis.text.x = element_text(size=13)
                                                    ,plot.margin = margin(t=15, r=20, b=0, l=5)),
                           
                           map_var_EU[[4]] + theme( legend.title = element_blank()
                                                    ,axis.title.y = element_blank()
                                                    ,axis.text.y = element_blank()
                                                    ,axis.title.x = element_blank()
                                                    ,axis.text.x = element_text(size=13)
                                                    ,plot.margin = margin(t=15, r=-5, b=0, l=-10)),
                           label_size = 25,
                           labels = c("    1h", "24h"),
                           nrow=1, 
                           hjust = -4,
                           rel_heights = c(1, 1),
                           rel_widths = c(1, 1),
                           align = 'h') 
    all.plots
    all.plots.title = plot_grid( title1,
                                 all.plots,
                                 ncol =1, nrow=2,
                                 rel_heights = c(0.05, 1))
    ggsave(all.plots.title, filename=paste0(path_output[v], '_EU_3.png'), width = 14, height =9, dpi = 150, bg = "white")
    
    
  } 
}

















################################################################################
plot_maps_nchange = function(var_list,
                             lon_tmp,
                             lat_tmp,
                             durations,
                             mask_shp,
                             flag_signif=T,
                             significance,
                             name_var="Shape intercept",
                             type_estimate="MaxPost",
                             name_model = "SMEV Nonstationary",
                             path_output,
                             limits,
                             flag.limits=T,
                             scale_grad,
                             crs_map = crswgs84){
  ################################################################################
  # var = list() of matrices for each duration (e.g., a_w_MAP_ns)
  
  
  for (v in 1:length(var_list)){
    
    # initialize variable plot list:
    map_var_EU = list()
    map_var_ITA = list()
    a_w_MAP_df=c()
    
    # durations
    for (dur in 1:length(durations)){
      
      ############
      # MAP
      ############
      # create dataframe for the map plot:
      a_w_MAP_df <- data.frame(x=lon_tmp, y=lat_tmp, z=var_list[[v]][[dur]])
      # trend_df$z[is.na(trend_df$z)]=0
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3), 
                                  z=round(a_w_MAP_df$z,3))
      rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
      
      # increasing AMAX
      signif_AMAX_pos_df <- data.frame(x=lon_tmp, y=lat_tmp, 
                                       z=significance$signif_MK_AMAX[[dur]], 
                                       d=significance$sensslope_AMAX[[dur]])
      signif_AMAX_pos_df$z[signif_AMAX_pos_df$z==F]=NA
      signif_AMAX_pos_df$z[signif_AMAX_pos_df$d<0] =NA
      signif_AMAX_pos_df = na.omit(signif_AMAX_pos_df)
      signif_AMAX_pos_df = subset(signif_AMAX_pos_df, z != 0)  
      
      # decreasing AMAX:
      signif_AMAX_neg_df <- data.frame(x=lon_tmp, y=lat_tmp,
                                       z=significance$signif_MK_AMAX[[dur]],
                                       d=significance$sensslope_AMAX[[dur]])
      signif_AMAX_neg_df$z[signif_AMAX_neg_df$z==F]=NA
      signif_AMAX_neg_df$z[signif_AMAX_neg_df$d>=0]=NA
      signif_AMAX_neg_df = na.omit(signif_AMAX_neg_df)
      signif_AMAX_neg_df = subset(signif_AMAX_neg_df, z != 0)  
      
      
      # nonstationary model is more likely than stationary (DIC)
      signif_ns_DIC_df_3 <- data.frame(x=lon_tmp, y=lat_tmp, z=significance$DIC_diff3[[dur]])
      signif_ns_DIC_df_3$z[signif_ns_DIC_df_3$z==F]=NA
      signif_ns_DIC_df_3 = na.omit(signif_ns_DIC_df_3)
      signif_ns_DIC_df_3 = subset(signif_ns_DIC_df_3, z != 0)      
      
      signif_ns_DIC_df_1 <- data.frame(x=lon_tmp, y=lat_tmp, z=significance$DIC_diff1[[dur]])
      signif_ns_DIC_df_1$z[signif_ns_DIC_df_1$z==F]=NA
      signif_ns_DIC_df_1 = na.omit(signif_ns_DIC_df_1)
      signif_ns_DIC_df_1 = subset(signif_ns_DIC_df_1, z != 0)  
      

      # signif "n"
      signif_n_df  = Snn[[jj]]$signif_trend_n[[dur]]
      signif_n_df <- data.frame(x=lon_tmp, y=lat_tmp,
                                    z=significance$signif_MK_n[[dur]],
                                    d=significance$sensslope_n[[dur]])
      signif_n_df$z[signif_n_df$z==F]=NA
      signif_n_df = na.omit(signif_n_df)
      signif_n_df = subset(signif_n_df, z != 0)
      
      # # increasing "n"
      # signif_n_pos_df  = Snn[[jj]]$signif_trend_n[[dur]]
      # signif_n_pos_df <- data.frame(x=lon_tmp, y=lat_tmp, 
      #                               z=significance$signif_MK_n[[dur]], 
      #                               d=significance$sensslope_n[[dur]])
      # signif_n_pos_df$z[signif_n_pos_df$z==F]=NA
      # signif_n_pos_df$z[signif_n_pos_df$d<0] =NA
      # signif_n_pos_df = na.omit(signif_n_pos_df)
      # signif_n_pos_df = subset(signif_n_pos_df, z != 0)  
      # 
      # # decreasing "n":
      # signif_n_neg_df <- data.frame(x=lon_tmp, y=lat_tmp,
      #                                  z=significance$signif_MK_n[[dur]],
      #                                  d=significance$sensslope_n[[dur]])
      # signif_n_neg_df$z[signif_n_neg_df$z==F]=NA
      # signif_n_neg_df$z[signif_n_neg_df$d>=0]=NA
      # signif_n_neg_df = na.omit(signif_n_neg_df)
      # signif_n_neg_df = subset(signif_n_neg_df, z != 0)  
      # 

      
      # ALL TILE 
      #####################################
      map_a_w_MAP = ggplot() + 
        theme_bw()+
        scale_x_continuous(expand = c(0, 0), name = "Longitude (°E)")+
        scale_y_continuous(expand = c(0, 0), name = "Latitude (°N)")+
        geom_tile(data=a_w_MAP_df_NEW, aes(x=x, y=y, fill=z), alpha=1) + 
        #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
        geom_sf(data = gg_world, fill = NA, color = "grey70") +
        geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+
        coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
      
      if (flag_signif){
        map_a_w_MAP=map_a_w_MAP+
          #geom_point(data=signif_n_pos_df, aes(x=x, y=y, size=d), fill="blue", color="blue", shape=24)+
          # geom_point(data=signif_n_neg_df, aes(x=x, y=y, size=d), fill="red", color="red",  shape=25)+
          geom_point(data=signif_n_df, aes(x=x, y=y, size=z),   fill="black", color="black", shape=21, size = 0.5)
        #geom_point(data=signif_ns_DIC_df, aes(x=x, y=y, size=z), shape = 1, size = 1)
        
      }
      map_a_w_MAP=map_a_w_MAP+
        theme(
          # legend.position = "none"
          ,legend.text     = element_text(size=15)
          ,legend.title    = element_text(size=20))
      
      
      if (flag.limits){
        if (scale_grad[[v]]$name=="turbo"){
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_viridis(option="turbo",
                               limits = limits[[v]],
                               name = paste0(name_var, ' (', type_estimate[v],
                                             '), ', name_model),
                               oob    = scales::squish)
          
        } else {
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name = paste0(name_var, '(',
                                               type_estimate[v], '), ', name_model),
                                 limits = limits[[v]],
                                 oob    = scales::squish)
          
        }
        
      } else {
        if (scale_grad[[v]]$name=="turbo"){
          
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_viridis(option="turbo",
                               name = paste0(name_var, ' (', type_estimate[v],
                                             '), ', name_model))
        } else {
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name = paste0(name_var, '(',
                                               type_estimate[v], '), ', name_model))
        }
      }
      # scale_fill_gradient2(low = scale_grad[[v]]$low, 
      #                      high =scale_grad[[v]]$high, 
      #                      mid = scale_grad[[v]]$mid, 
      #                      midpoint = scale_grad[[v]]$midpoint,
      #                      na.value = NA,
      #                      name = paste0(name_var, ' (', type_estimate[v],
      #                                    '), ', name_model),
      #                      limits = limits[[v]]) +
      
      
      # ggsave(map_a_w_MAP, filename =paste0(dir_results, "/", folder_name, 
      #                                        "/map_a_w_MAP", ".png"),
      #        width = 15, height =15, dpi = 200)
      # map_a_w_MAP
      map_var_EU = c(map_var_EU, list(map_a_w_MAP))
      
      
      
      
      
      
      
      # ITALY (pixel level NOT AGGREGATED):  
      ####################################
      # Replace nan with 1:
      a_w_MAP_df$z[a_w_MAP_df$z==9999] =NA
      a_w_MAP_df_nonan <- a_w_MAP_df
      a_w_MAP_df_nonan$z[is.na(a_w_MAP_df_nonan$z)]=0
      # colMap <- colorRampPalette(c("yellow", "red"))(3)
      a_w_MAP_df_nonan_NEW = data.frame(x=round(a_w_MAP_df_nonan$x,3),
                                        y=round(a_w_MAP_df_nonan$y,3), 
                                        z=round(a_w_MAP_df_nonan$z,3))
      rast_a_w_MAP_nonan =rasterFromXYZ(a_w_MAP_df_nonan_NEW, crs=crs_map) 
      ## crop and mask (only Italy extent):
      r=rast_a_w_MAP_nonan
      r2 <- crop(r, extent(mask_shp))
      r3 <- mask(r2, mask_shp)
      test_spdf <- as(r3, "SpatialPixelsDataFrame")
      test_df <- as.data.frame(test_spdf)
      
      map_a_w_MAP_ITALY= ggplot() + 
        # geom_sf(data = gg_world, fill = NA, color = "grey50")+
        theme_bw()+
        scale_x_continuous(expand = c(0, 0), name = "Longitude (°E)")+
        scale_y_continuous(expand = c(0, 0), name = "Latitude (°N)")+
        geom_tile(data=na.omit(test_df), aes(x=x, y=y, fill=z), alpha=1) + 
        #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
        geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+
        coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)+
        theme(
          # legend.position = "none"
          ,legend.text     = element_text(size=15)
          ,legend.title    = element_text(size=20))
      
      
      if (flag.limits){
        if (scale_grad[[v]]$name=="turbo"){
          map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
            scale_fill_viridis(option="turbo",
                               limits = limits[[v]],
                               name = paste0(name_var, ' (', type_estimate[v],
                                             '), ', name_model),
                               oob    = scales::squish)
        } else {
          map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name = paste0(name_var, '(',
                                               type_estimate[v], '), ', name_model),
                                 limits = limits[[v]],
                                 oob    = scales::squish)
        }
        
      } else {
        if (scale_grad[[v]]$name=="turbo"){
          
          map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
            scale_fill_viridis(option="turbo",
                               name = paste0(name_var, ' (', type_estimate[v],
                                             '), ', name_model))
        } else {
          map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name = paste0(name_var, '(',
                                               type_estimate[v], '), ', name_model))
        }
      }
      
      map_var_ITA = c(map_var_ITA, list(map_a_w_MAP_ITALY))
      
      
    }
    
    parent_dir= dirname(path_output[1])
    if (!dir.exists(parent_dir)){
      dir.create(parent_dir, recursive = T)
    }
    
    # merging the plots for different durations:
    plot_maps_gridplot(plotlist=map_var_EU, 
                       title= paste0(name_var," (", type_estimate[v], ") - ", name_model), 
                       pathout=paste0(path_output[v], '_EU.png'))
    
    plot_maps_gridplot(plotlist=map_var_ITA, 
                       title= paste0(name_var," (", type_estimate[v], ") - ", name_model), 
                       pathout=paste0(path_output[v], '_ITA.png'))
  }  
  
}
























################################################################################
plot_maps_AMAXchange = function(var_list,
                             lon_tmp,
                             lat_tmp,
                             durations,
                             mask_shp,
                             flag_signif=T,
                             significance,
                             name_var="Shape intercept",
                             type_estimate="MaxPost",
                             name_model = "SMEV Nonstationary",
                             path_output,
                             limits,
                             flag.limits=T,
                             scale_grad,
                             crs_map = crswgs84){
  ################################################################################
  # var = list() of matrices for each duration (e.g., a_w_MAP_ns)
  
  
  for (v in 1:length(var_list)){
    
    # initialize variable plot list:
    map_var_EU = list()
    map_var_ITA = list()
    a_w_MAP_df=c()
    
    # durations
    for (dur in 1:length(durations)){
      
      ############
      # MAP
      ############
      # create dataframe for the map plot:
      a_w_MAP_df <- data.frame(x=lon_tmp, y=lat_tmp, z=var_list[[v]][[dur]])
      # trend_df$z[is.na(trend_df$z)]=0
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3), 
                                  z=round(a_w_MAP_df$z,3))
      rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
      
      # increasing AMAX
      signif_AMAX_pos_df <- data.frame(x=lon_tmp, y=lat_tmp, 
                                       z=significance$signif_MK_AMAX[[dur]], 
                                       d=significance$sensslope_AMAX[[dur]])
      signif_AMAX_pos_df$z[signif_AMAX_pos_df$z==F]=NA
      signif_AMAX_pos_df$z[signif_AMAX_pos_df$d<0] =NA
      signif_AMAX_pos_df = na.omit(signif_AMAX_pos_df)
      signif_AMAX_pos_df = subset(signif_AMAX_pos_df, z != 0)  
      
      # decreasing AMAX:
      signif_AMAX_neg_df <- data.frame(x=lon_tmp, y=lat_tmp,
                                       z=significance$signif_MK_AMAX[[dur]],
                                       d=significance$sensslope_AMAX[[dur]])
      signif_AMAX_neg_df$z[signif_AMAX_neg_df$z==F]=NA
      signif_AMAX_neg_df$z[signif_AMAX_neg_df$d>=0]=NA
      signif_AMAX_neg_df = na.omit(signif_AMAX_neg_df)
      signif_AMAX_neg_df = subset(signif_AMAX_neg_df, z != 0)  
      
      
      # nonstationary model is more likely than stationary (DIC)
      signif_ns_DIC_df_3 <- data.frame(x=lon_tmp, y=lat_tmp, z=significance$DIC_diff3[[dur]])
      signif_ns_DIC_df_3$z[signif_ns_DIC_df_3$z==F]=NA
      signif_ns_DIC_df_3 = na.omit(signif_ns_DIC_df_3)
      signif_ns_DIC_df_3 = subset(signif_ns_DIC_df_3, z != 0)      
      
      signif_ns_DIC_df_1 <- data.frame(x=lon_tmp, y=lat_tmp, z=significance$DIC_diff1[[dur]])
      signif_ns_DIC_df_1$z[signif_ns_DIC_df_1$z==F]=NA
      signif_ns_DIC_df_1 = na.omit(signif_ns_DIC_df_1)
      signif_ns_DIC_df_1 = subset(signif_ns_DIC_df_1, z != 0)  
      
      
      # signif AMAX"
      signif_amax_df <- data.frame(x=lon_tmp, y=lat_tmp,
                                   z=significance$signif_MK_AMAX[[dur]],
                                   d=significance$sensslope_AMAX[[dur]])
      signif_amax_df$z[signif_amax_df$z==F]=NA
      signif_amax_df = na.omit(signif_amax_df)
      signif_amax_df = subset(signif_amax_df, z != 0)
      
      
      
      
      
      # ALL TILE 
      #####################################
      map_a_w_MAP = ggplot() + 
        theme_bw()+
        scale_x_continuous(expand = c(0, 0), name = "Longitude (°E)")+
        scale_y_continuous(expand = c(0, 0), name = "Latitude (°N)")+
        geom_tile(data=a_w_MAP_df_NEW, aes(x=x, y=y, fill=z), alpha=0.7) + 
        #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
        geom_sf(data = gg_world, fill = NA, color = "grey70") +
        geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+
        coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
        
      if (flag_signif){
        map_a_w_MAP=map_a_w_MAP+
          #geom_point(data=signif_n_pos_df, aes(x=x, y=y, size=d), fill="blue", color="blue", shape=24)+
          # geom_point(data=signif_n_neg_df, aes(x=x, y=y, size=d), fill="red", color="red",  shape=25)+
          geom_point(data=signif_amax_df, aes(x=x, y=y, size=z),   fill="black", color="black", shape=21, size = 0.7)
        #geom_point(data=signif_ns_DIC_df, aes(x=x, y=y, size=z), shape = 1, size = 1)
        
      }
      map_a_w_MAP=map_a_w_MAP+
        theme(
          # legend.position = "none"
          ,legend.text     = element_text(size=15)
          ,legend.title    = element_text(size=20))
      
      
      if (flag.limits){
        if (scale_grad[[v]]$name=="turbo"){
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_viridis(option="turbo",
                               limits = limits[[v]],
                               name = paste0(name_var, ' (', type_estimate[v],
                                             '), ', name_model),
                               oob    = scales::squish)
          
        } else {
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name = paste0(name_var, '(',
                                               type_estimate[v], '), ', name_model),
                                 limits = limits[[v]],
                                 oob    = scales::squish)
          
        }
        
      } else {
        if (scale_grad[[v]]$name=="turbo"){
          
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_viridis(option="turbo",
                               name = paste0(name_var, ' (', type_estimate[v],
                                             '), ', name_model))
        } else {
          map_a_w_MAP=map_a_w_MAP+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name = paste0(name_var, '(',
                                               type_estimate[v], '), ', name_model))
        }
      }
      # scale_fill_gradient2(low = scale_grad[[v]]$low, 
      #                      high =scale_grad[[v]]$high, 
      #                      mid = scale_grad[[v]]$mid, 
      #                      midpoint = scale_grad[[v]]$midpoint,
      #                      na.value = NA,
      #                      name = paste0(name_var, ' (', type_estimate[v],
      #                                    '), ', name_model),
      #                      limits = limits[[v]]) +
      
      
      # ggsave(map_a_w_MAP, filename =paste0(dir_results, "/", folder_name, 
      #                                        "/map_a_w_MAP", ".png"),
      #        width = 15, height =15, dpi = 200)
      # map_a_w_MAP
      map_var_EU = c(map_var_EU, list(map_a_w_MAP))
      
      
      
      
      
      
      
      # ITALY (pixel level NOT AGGREGATED):  
      ####################################
      # Replace nan with 1:
      a_w_MAP_df$z[a_w_MAP_df$z==9999] =NA
      a_w_MAP_df_nonan <- a_w_MAP_df
      a_w_MAP_df_nonan$z[is.na(a_w_MAP_df_nonan$z)]=0
      # colMap <- colorRampPalette(c("yellow", "red"))(3)
      a_w_MAP_df_nonan_NEW = data.frame(x=round(a_w_MAP_df_nonan$x,3),
                                        y=round(a_w_MAP_df_nonan$y,3), 
                                        z=round(a_w_MAP_df_nonan$z,3))
      rast_a_w_MAP_nonan =rasterFromXYZ(a_w_MAP_df_nonan_NEW, crs=crs_map) 
      ## crop and mask (only Italy extent):
      r=rast_a_w_MAP_nonan
      r2 <- crop(r, extent(mask_shp))
      r3 <- mask(r2, mask_shp)
      test_spdf <- as(r3, "SpatialPixelsDataFrame")
      test_df <- as.data.frame(test_spdf)
      
      map_a_w_MAP_ITALY= ggplot() + 
        # geom_sf(data = gg_world, fill = NA, color = "grey50")+
        theme_bw()+
        scale_x_continuous(expand = c(0, 0), name = "Longitude (°E)")+
        scale_y_continuous(expand = c(0, 0), name = "Latitude (°N)")+
        geom_tile(data=na.omit(test_df), aes(x=x, y=y, fill=z), alpha=1) + 
        #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
        geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+
        coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)+
        theme(
          # legend.position = "none"
          ,legend.text     = element_text(size=15)
          ,legend.title    = element_text(size=20))
      
      
      if (flag.limits){
        if (scale_grad[[v]]$name=="turbo"){
          map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
            scale_fill_viridis(option="turbo",
                               limits = limits[[v]],
                               name = paste0(name_var, ' (', type_estimate[v],
                                             '), ', name_model),
                               oob    = scales::squish)
        } else {
          map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name = paste0(name_var, '(',
                                               type_estimate[v], '), ', name_model),
                                 limits = limits[[v]],
                                 oob    = scales::squish)
        }
        
      } else {
        if (scale_grad[[v]]$name=="turbo"){
          
          map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
            scale_fill_viridis(option="turbo",
                               name = paste0(name_var, ' (', type_estimate[v],
                                             '), ', name_model))
        } else {
          map_a_w_MAP_ITALY=map_a_w_MAP_ITALY+
            scale_fill_gradient2(low = scale_grad[[v]]$low,
                                 high =scale_grad[[v]]$high,
                                 mid = scale_grad[[v]]$mid,
                                 midpoint = scale_grad[[v]]$midpoint,
                                 na.value = NA,
                                 name = paste0(name_var, '(',
                                               type_estimate[v], '), ', name_model))
        }
      }
      
      map_var_ITA = c(map_var_ITA, list(map_a_w_MAP_ITALY))
      
      
    }
    
    parent_dir= dirname(path_output[1])
    if (!dir.exists(parent_dir)){
      dir.create(parent_dir, recursive = T)
    }
    
    # merging the plots for different durations:
    plot_maps_gridplot(plotlist=map_var_EU, 
                       title= paste0(name_var," (", type_estimate[v], ") - ", name_model), 
                       pathout=paste0(path_output[v], '_EU.png'))
    
    plot_maps_gridplot(plotlist=map_var_ITA, 
                       title= paste0(name_var," (", type_estimate[v], ") - ", name_model), 
                       pathout=paste0(path_output[v], '_ITA.png'))
  }  
  
}























################################################################################
# PLOT % CHANGES of AMAX, quant2yrs, Shape, Scale, n (SMEV MAXPOST)  
################################################################################
plot_maps_changes_amax2yrs=function(
                             var_change_list = list(),   
                             lon_tmp       = NULL,
                             lat_tmp       = NULL,
                             durations     = c(60,180,360,1440), # in minutes
                             mask_shp      = NULL,
                             mask_world    = NULL,
                             flag_signif   = F,
                             significance  = NULL,
                             sat_acron     = 'TEST',
                             period        = '',
                             year_ref      = NA,
                             year_ref_id   = NA,
                             name_file_txt = paste0('signif_pixels.txt'),
                             path_output   = "/map_test",
                             limits        = list(c(-20, 20)),
                             flag.limits   = T,
                             scale_grad    = list(),
                             crs_map       = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")){
################################################################################
  # initialize paths
  dir.res_maps_amax_qnt =path_output
  dir.create(dir.res_maps_amax_qnt)
  cat(paste0(sat_acron,' - ', period), file = paste0(dir.res_maps_amax_qnt, name_file_txt), sep="\n")
  
  # initialize output lists
  pixels_signif_AMAX=pixels_signif_AMAX_LAND=list()
  percent_signif_AMAX=percent_signif_AMAX_LAND=list()
  pixels_signif_DIC3=pixels_signif_DIC3_LAND=list()
  percent_signif_DIC3=percent_signif_DIC3_LAND=list()
  pixels_signif_DIC1=pixels_signif_DIC1_LAND=list()
  percent_signif_DIC1=percent_signif_DIC1_LAND=list()
  pixels_signif_n=pixels_signif_n_LAND=list()
  percent_signif_n=percent_signif_n_LAND=list()
  
  # cycle on rainfall durations (4 durations implemented for instance!)
  for (dur in 1:4){
    
    message(paste0('DURATION = ', durations[dur]/60, 'h'))
    cat(paste0("##########################################\n DURATION=", 
               durations[dur]/60, "h \n##########################################"), 
        file = paste0(dir.res_maps_amax_qnt, name_file_txt), sep="\n", append=T)
    
    
    ############################################################################
    #                               AMAX:
    ############################################################################
    cat(paste0("\nAMAX:"), file = paste0(dir.res_maps_amax_qnt, name_file_txt),
        sep="\n", append=T)
    # Initialize plot:
    flag_signif   = T
    name_var      = "AMAX change per decade"
    type_estimate = c("sensslope")
    name_model    = "AMAX"
    flag.limits   = T
    scale_grad    = list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)
    crs_map       = crswgs84
    v             = 1 # only MaxPost is considered for instance!
    
    # initialize variable plot list:
    AMAXchange    = AMAXchange_LAND = c()
    a_w_MAP_df    = c() # reinitialize variable!
    rast_a_w_map_LAND_nonan = c()
    
    # create dataframe for the map plot:
    a_w_MAP_df <- data.frame(x=lon_tmp,
                             y=lat_tmp, 
                             z=var_change_list$AMAX_change[[dur]])
    # trend_df$z[is.na(trend_df$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_NEW=a_w_MAP_df
    } else {
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3),
                                  z=round(a_w_MAP_df$z,3))
    }
    rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
    
    # increasing AMAX
    signif_AMAX_pos_df <- data.frame(x=lon_tmp, y=lat_tmp, 
                                     z=significance$signif_MK_AMAX[[dur]], 
                                     d=significance$sensslope_AMAX[[dur]])
    signif_AMAX_pos_df$z[signif_AMAX_pos_df$z==F]=NA
    signif_AMAX_pos_df$z[signif_AMAX_pos_df$d<0] =NA
    signif_AMAX_pos_df = na.omit(signif_AMAX_pos_df)
    signif_AMAX_pos_df = subset(signif_AMAX_pos_df, z != 0)  
    
    # decreasing AMAX:
    signif_AMAX_neg_df <- data.frame(x=lon_tmp, y=lat_tmp,
                                     z=significance$signif_MK_AMAX[[dur]],
                                     d=significance$sensslope_AMAX[[dur]])
    signif_AMAX_neg_df$z[signif_AMAX_neg_df$z==F]=NA
    signif_AMAX_neg_df$z[signif_AMAX_neg_df$d>=0]=NA
    signif_AMAX_neg_df = na.omit(signif_AMAX_neg_df)
    signif_AMAX_neg_df = subset(signif_AMAX_neg_df, z != 0)  
    
    
    # nonstationary model is more likely than stationary (DIC)
    signif_ns_DIC_df_3 <- data.frame(x=lon_tmp, 
                                     y=lat_tmp, 
                                     z=significance$DIC_diff3[[dur]])
    signif_ns_DIC_df_3$z[signif_ns_DIC_df_3$z==F]=NA
    signif_ns_DIC_df_3 = na.omit(signif_ns_DIC_df_3)
    signif_ns_DIC_df_3 = subset(signif_ns_DIC_df_3, z != 0)      
    
    signif_ns_DIC_df_1 <- data.frame(x=lon_tmp, y=lat_tmp, 
                                     z=significance$DIC_diff1[[dur]])
    signif_ns_DIC_df_1$z[signif_ns_DIC_df_1$z==F]=NA
    signif_ns_DIC_df_1 = na.omit(signif_ns_DIC_df_1)
    signif_ns_DIC_df_1 = subset(signif_ns_DIC_df_1, z != 0)  
    
    
    # signif AMAX"
    if (sat_acron=='GSMAP'){
      signif_amax_df <- data.frame(x=lon_tmp,
                                   y=lat_tmp,
                                   z=significance$signif_MK_AMAX[[dur]],
                                   d=significance$sensslope_AMAX[[dur]])
    } else {
      signif_amax_df <- data.frame(x=round(lon_tmp,3),
                                   y=round(lat_tmp,3),
                                   z=significance$signif_MK_AMAX[[dur]],
                                   d=significance$sensslope_AMAX[[dur]])
    }
    signif_amax_df$z[signif_amax_df$z==F]=NA
    signif_amax_df = na.omit(signif_amax_df)
    signif_amax_df = subset(signif_amax_df, z != 0)
    
    # compute percentage respect total amount of pixels
    pixels_signif_AMAX[[dur]] = length(signif_amax_df$z)
    total_pixels = length(a_w_MAP_df$z)  # segm_df$z
    percent_signif_AMAX[[dur]] = pixels_signif_AMAX[[dur]]/total_pixels*100
    print(paste0("pixels with signif. trend in AMAX (WHOLE) = ", pixels_signif_AMAX[[dur]], 
                 " / ", total_pixels, " (", round(percent_signif_AMAX[[dur]],1), "%)"))
    cat(paste0("pixels with signif. trend in AMAX (WHOLE) = ", pixels_signif_AMAX[[dur]], 
               " / ", total_pixels, " (", round(percent_signif_AMAX[[dur]],1), "%)"), 
        file = paste0(dir.res_maps_amax_qnt,name_file_txt), sep="\n", append=T)
    
    
    AMAXchange = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=a_w_MAP_df_NEW, aes(x=x, y=y, fill=z), alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
      geom_sf(data = gg_world, fill = NA, color = "grey70") +
      geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+
      coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
    
    if (flag_signif){
      AMAXchange=AMAXchange+
        #geom_point(data=signif_n_pos_df, aes(x=x, y=y, size=d), fill="blue", color="blue", shape=24)+
        # geom_point(data=signif_n_neg_df, aes(x=x, y=y, size=d), fill="red", color="red",  shape=25)+
        geom_point(data=signif_amax_df, aes(x=x, y=y, size=z),   fill="black", color="black", shape=21, size = 0.4)
      #geom_point(data=signif_ns_DIC_df, aes(x=x, y=y, size=z), shape = 1, size = 1)
    }
    AMAXchange=AMAXchange+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y=element_text(size=13)
        # legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    # only land:
    ############
    # Replace nan with 1:
    a_w_MAP_df_LAND=a_w_MAP_df_NEW
    a_w_MAP_df_LAND$z[a_w_MAP_df_LAND$z==9999] =NA
    a_w_MAP_df_LAND_nonan <- a_w_MAP_df_LAND
    a_w_MAP_df_LAND_nonan$z[is.na(a_w_MAP_df_LAND_nonan$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan
    } else {
      a_w_MAP_df_LAND_nonan_NEW=data.frame(x=round(a_w_MAP_df_LAND_nonan$x,3),
                                           y=round(a_w_MAP_df_LAND_nonan$y,3), 
                                           z=round(a_w_MAP_df_LAND_nonan$z,3))
    }
    rast_a_w_map_LAND_nonan=rasterFromXYZ(a_w_MAP_df_LAND_nonan_NEW, crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2 <- crop(r, extent(mask_shp), snap='in')
    r3 <- mask(r2, mask_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df <- as.data.frame(test_spdf)
    
    rast_signif_amax_LAND =rasterFromXYZ(signif_amax_df, crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_signif_amax_LAND
    r2 <- crop(r, as(extent(mask_shp), 'SpatialPolygons')) #extent(mask_shp))
    r3 <- mask(r2, mask_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df_signif_amax <- as.data.frame(test_spdf)
    
    
    pixels_signif_AMAX_LAND[[dur]] = length(test_df_signif_amax$d)
    total_pixels_LAND = length(test_df$z)
    percent_signif_AMAX_LAND[[dur]] = pixels_signif_AMAX_LAND[[dur]]/total_pixels_LAND*100
    print(paste0("pixels with signif. trend in AMAX (LAND) = ", pixels_signif_AMAX_LAND[[dur]], 
                 " / ", total_pixels_LAND, " (", round(percent_signif_AMAX_LAND[[dur]],1), "%)"))
    cat(paste0("pixels with signif. trend in AMAX (LAND) = ", pixels_signif_AMAX_LAND[[dur]], 
               " / ", total_pixels_LAND, " (", round(percent_signif_AMAX_LAND[[dur]],1), "%)"), 
        file = paste0(dir.res_maps_amax_qnt, name_file_txt), sep="\n", append=T)
    
    
    AMAXchange_LAND = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=na.omit(test_df), aes(x=x, y=y, fill=z), alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
      geom_sf(data = gg_world, fill = NA, color = "grey70") +
      geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+
      coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
    if (flag_signif){
      AMAXchange_LAND=AMAXchange_LAND+
        geom_point(data=test_df_signif_amax, aes(x=x, y=y, size=z),  
                   fill="black", color="black", shape=21, size = 0.4)
    }
    AMAXchange_LAND=AMAXchange_LAND+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y=element_text(size=13)
        # legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    

    ############################################################################
    #                        2yrs quantile:
    ############################################################################
    cat(paste0("\n2yrs qnt:"), 
        file=paste0(dir.res_maps_amax_qnt, name_file_txt),
        sep="\n",append=T)
    
    flag_signif   = T
    flag.limits   = T
    scale_grad    = list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)
    crs_map       = crswgs84
    v             = 2
    name_var      = "Change in 2yrs quantiles per decade"
    type_estimate = c("MaxPost")
    name_model    = "SMEV Nonstationary"
    
    # initialize variable plot list:
    twoyrchange   = twoyrchange_LAND = c()
    a_w_MAP_df    = c()
    rast_a_w_map_LAND_nonan=c()
    
    # create dataframe for the map plot:
    a_w_MAP_df <- data.frame(x=lon_tmp, 
                             y=lat_tmp, 
                             z=var_change_list$qnt2yr_change[[dur]])
    # trend_df$z[is.na(trend_df$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_NEW=a_w_MAP_df
    } else {
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3), 
                                  z=round(a_w_MAP_df$z,3))
    }
    rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
    
    # nonstationary model is more likely than stationary (DIC)
    if (sat_acron=='GSMAP'){
      signif_ns_DIC_df_3 <- data.frame(x=lon_tmp,
                                       y=lat_tmp, 
                                       z=significance$DIC_diff3[[dur]])   
    } else {
      signif_ns_DIC_df_3 <- data.frame(x=round(lon_tmp,3),
                                       y=round(lat_tmp,3), 
                                       z=significance$DIC_diff3[[dur]])
    }
    signif_ns_DIC_df_3$z[signif_ns_DIC_df_3$z==F]=NA
    signif_ns_DIC_df_3 = na.omit(signif_ns_DIC_df_3)
    signif_ns_DIC_df_3 = subset(signif_ns_DIC_df_3, z != 0)      
    
    if (sat_acron=='GSMAP'){
      signif_ns_DIC_df_1 <- data.frame(x=lon_tmp,
                                       y=lat_tmp,
                                       z=significance$DIC_diff1[[dur]])
    } else {
      signif_ns_DIC_df_1 <- data.frame(x=round(lon_tmp,3),
                                       y=round(lat_tmp,3),
                                       z=significance$DIC_diff1[[dur]])
    }
    signif_ns_DIC_df_1$z[signif_ns_DIC_df_1$z==F]=NA
    signif_ns_DIC_df_1 = na.omit(signif_ns_DIC_df_1)
    signif_ns_DIC_df_1 = subset(signif_ns_DIC_df_1, z != 0)  
    
    # compute percentage respect total amount of pixels:
    pixels_signif_DIC3[[dur]] = length(signif_ns_DIC_df_3$z)
    total_pixels = length(a_w_MAP_df$z)  # segm_df$z
    percent_signif_DIC3[[dur]] = pixels_signif_DIC3[[dur]]/total_pixels*100
    print(paste0("pixels with signif. NS deltaDIC>=3 (WHOLE) = ", pixels_signif_DIC3[[dur]], 
                 " / ", total_pixels, " (", round(percent_signif_DIC3[[dur]],1), "%)"))
    cat(paste0("pixels with signif. NS deltaDIC>=3 (WHOLE) = ", pixels_signif_DIC3[[dur]], 
               " / ", total_pixels, " (", round(percent_signif_DIC3[[dur]],1), "%)"), 
        file = paste0(dir.res_maps_amax_qnt, name_file_txt), sep="\n", append=T)
    # compute percentage respect total amount of pixels:
    pixels_signif_DIC1[[dur]] = length(signif_ns_DIC_df_1$z)
    total_pixels = length(a_w_MAP_df$z)  # segm_df$z
    percent_signif_DIC1[[dur]] = pixels_signif_DIC1[[dur]]/total_pixels*100
    print(paste0("pixels with signif. NS deltaDIC>=1 (WHOLE) = ", pixels_signif_DIC1[[dur]], 
                 " / ", total_pixels, " (", round(percent_signif_DIC1[[dur]],1), "%)"))
    cat(paste0("pixels with signif. NS deltaDIC>=1 (WHOLE) = ", pixels_signif_DIC1[[dur]], 
               " / ", total_pixels, " (", round(percent_signif_DIC1[[dur]],1), "%)"), 
        file = paste0(dir.res_maps_amax_qnt, name_file_txt), sep="\n", append=T)
    
    
    twoyrchange = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=a_w_MAP_df_NEW, aes(x=x, y=y, fill=z), alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
      geom_sf(data = gg_world, fill = NA, color = "grey70")+ #,linewidth=0.5) +
      geom_sf(data = test, colour = "black",fill=NA,alpha=0.2)+ #,linewidth=0.5)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)
    if (flag_signif){
      twoyrchange=twoyrchange+
        # geom_point(data=signif_AMAX_pos_df, aes(x=x, y=y, size=d), fill="blue", color="blue", shape=24)+
        # geom_point(data=signif_AMAX_neg_df, aes(x=x, y=y, size=d), fill="red", color="red",  shape=25)+
        #geom_point(data=signif_ns_DIC_df_1, aes(x=x, y=y, size=z),   fill="black", color="black", shape=21, size = 0.7)
        geom_point(data=signif_ns_DIC_df_3, aes(x=x, y=y,size=z),   
                   fill="black", color="black",shape=21,size=0.4)
      #geom_point(data=signif_ns_DIC_df, aes(x=x, y=y, size=z), shape = 1, size = 1)
    }
    twoyrchange=twoyrchange+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y=element_text(size=13)
        # legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    # only land:
    ############
    # Replace nan with 1:
    a_w_MAP_df_LAND=a_w_MAP_df_NEW
    a_w_MAP_df_LAND$z[a_w_MAP_df_LAND$z==9999] =NA
    a_w_MAP_df_LAND_nonan <- a_w_MAP_df_LAND
    a_w_MAP_df_LAND_nonan$z[is.na(a_w_MAP_df_LAND_nonan$z)]=0
    # colMap <- colorRampPalette(c("yellow", "red"))(3)
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan
    } else {
      a_w_MAP_df_LAND_nonan_NEW =data.frame(x=round(a_w_MAP_df_LAND_nonan$x,3),
                                            y=round(a_w_MAP_df_LAND_nonan$y,3), 
                                            z=round(a_w_MAP_df_LAND_nonan$z,3))
    }
    rast_a_w_map_LAND_nonan =rasterFromXYZ(a_w_MAP_df_LAND_nonan_NEW, 
                                          crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2 <- crop(r, extent(mask_shp))
    r3 <- mask(r2, mask_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df <- as.data.frame(test_spdf)
    
    rast_signif_DIC3_LAND = rasterFromXYZ(signif_ns_DIC_df_3, crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_signif_DIC3_LAND
    r2 <- crop(r, extent(mask_shp))
    r3 <- mask(r2, mask_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df_signif_DIC3 <- as.data.frame(test_spdf)
    
    rast_signif_DIC1_LAND = rasterFromXYZ(signif_ns_DIC_df_1, crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_signif_DIC1_LAND
    r2 <- crop(r, extent(mask_shp))
    r3 <- mask(r2, mask_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df_signif_DIC1 <- as.data.frame(test_spdf)
    
    # compute percentage respect total amount of pixels:
    pixels_signif_DIC3_LAND[[dur]] = length(test_df_signif_DIC3$z)
    total_pixels = length(test_df$z)  # segm_df$z
    percent_signif_DIC3_LAND[[dur]] = pixels_signif_DIC3_LAND[[dur]]/total_pixels*100
    print(paste0("pixels with signif. NS deltaDIC>=3 (LAND) = ", pixels_signif_DIC3_LAND[[dur]], 
                 " / ", total_pixels, " (", round(percent_signif_DIC3_LAND[[dur]],1), "%)"))
    cat(paste0("pixels with signif. NS deltaDIC>=3 (LAND) = ", pixels_signif_DIC3_LAND[[dur]], 
               " / ", total_pixels, " (", round(percent_signif_DIC3_LAND[[dur]],1), "%)"), 
        file = paste0(dir.res_maps_amax_qnt, name_file_txt), sep="\n", append=T)
    pixels_signif_DIC1_LAND[[dur]] = length(test_df_signif_DIC1$z)
    percent_signif_DIC1_LAND[[dur]] = pixels_signif_DIC1_LAND[[dur]]/total_pixels*100
    print(paste0("pixels with signif. NS deltaDIC>=1 (LAND) = ", pixels_signif_DIC1_LAND[[dur]], 
                 " / ", total_pixels, " (", round(percent_signif_DIC1_LAND[[dur]],1), "%)"))
    cat(paste0("pixels with signif. NS deltaDIC>=1 (LAND) = ", pixels_signif_DIC1_LAND[[dur]], 
               " / ", total_pixels, " (", round(percent_signif_DIC1_LAND[[dur]],1), "%)"), 
        file = paste0(dir.res_maps_amax_qnt, name_file_txt), sep="\n", append=T)
    
    
    twoyrchange_LAND = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=na.omit(test_df), aes(x=x, y=y, fill=z), alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
      geom_sf(data=gg_world, fill=NA, color="grey70")+ #, linewidth=0.5) +
      geom_sf(data=test, colour="black", fill=NA, alpha=0.2)+ #, linewidth=0.5)+
      coord_sf(xlim=c(6,19), ylim=c(36,48), expand=F)
    if (flag_signif){
      twoyrchange_LAND=twoyrchange_LAND+
        # geom_point(data=signif_ns_DIC_df_1, aes(x=x, y=y, size=z),   
        #            fill="black", color="black", shape=21, size = 0.7)
        geom_point(data=test_df_signif_DIC3, aes(x=x, y=y, size=z),  
                   fill="black", color="black", shape=21, size=0.4)
    }
    twoyrchange_LAND=twoyrchange_LAND+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y=element_text(size=13)
        # legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    
    
    
    ############################################################################
    #                                 scale:
    ############################################################################
    cat(paste0("\nScale:"), 
        file=paste0(dir.res_maps_amax_qnt, name_file_txt),
        sep="\n",append=T)
    
    # limits      = list(c(-30, 30))
    flag_signif   = T
    flag.limits   = T
    scale_grad    = list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)
    crs_map       = crswgs84
    v             = 3
    name_var      = "Change in scale per decade"
    type_estimate = c("MaxPost")
    name_model    = "SMEV Nonstationary"
    year_ref      = year_halfperiod
    year_ref_id   = halfperiod
    
    # initialize variable plot list:
    scalechange       = scalechange_LAND = c()
    a_w_MAP_df        = c()
    scale_change      = c()
    rast_a_w_map_LAND_nonan=c()
    
    # % of change of slope per decade wrw halfperiod year (a sort of average change per decade):
    scale_change = var_change_list$b_C_MAP_ns[[dur]]*10 /
           (var_change_list$a_C_MAP_ns[[dur]]+var_change_list$b_C_MAP_ns[[dur]]*year_ref_id)*100
    # create dataframe for the map plot:
    a_w_MAP_df <- data.frame(x=lon_tmp, 
                             y=lat_tmp, 
                             z=scale_change)
    # trend_df$z[is.na(trend_df$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_NEW=a_w_MAP_df
    } else {
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3), 
                                  z=round(a_w_MAP_df$z,3))
    }
    rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
    
    # map:
    scalechange = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=a_w_MAP_df_NEW, aes(x=x, y=y, fill=z), alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
      geom_sf(data = gg_world, fill = NA, color="grey70")+ #,  linewidth=0.5) +
      geom_sf(data = test, colour = "black", fill=NA,alpha=0.2)+ #, linewidth=0.5)+ 
      # geom_point(data =df_NEW_SIGNIF_EU, aes(x=x, y=y, size= z), 
      #            fill="black", shape = 19, size=1)+  #pixels with signif. nonstationarity
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)
    if (flag_signif){
      scalechange=scalechange+
        geom_point(data=signif_ns_DIC_df_3, aes(x=x, y=y, size=z),
                   fill="black",color="black",shape=21,size=0.4)
    }
    scalechange=scalechange+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y=element_text(size=13)
        # legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    # only land:
    ############
    # Replace nan with 1:
    a_w_MAP_df_LAND=a_w_MAP_df_NEW
    a_w_MAP_df_LAND$z[a_w_MAP_df_LAND$z==9999] =NA
    a_w_MAP_df_LAND_nonan <- a_w_MAP_df_LAND
    a_w_MAP_df_LAND_nonan$z[is.na(a_w_MAP_df_LAND_nonan$z)]=0
    # colMap <- colorRampPalette(c("yellow", "red"))(3)
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan
    } else {
      a_w_MAP_df_LAND_nonan_NEW =data.frame(x=round(a_w_MAP_df_LAND_nonan$x,3),
                                            y=round(a_w_MAP_df_LAND_nonan$y,3), 
                                            z=round(a_w_MAP_df_LAND_nonan$z,3))
    }
    rast_a_w_map_LAND_nonan =rasterFromXYZ(a_w_MAP_df_LAND_nonan_NEW, 
                                          crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2 <- crop(r, extent(mask_shp))
    r3 <- mask(r2, mask_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df <- as.data.frame(test_spdf)
    
    
    scalechange_LAND = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=na.omit(test_df), aes(x=x, y=y, fill=z), alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
      geom_sf(data = gg_world, fill = NA, color = "grey70")+ #,  linewidth=0.5) +
      geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+ #, linewidth=0.5)+ 
      coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
    if (flag_signif){
      scalechange_LAND=scalechange_LAND+
        geom_point(data=test_df_signif_DIC3, aes(x=x, y=y, size=z),   
                   fill="black", color="black", shape=21, size = 0.4)
    }
    scalechange_LAND=scalechange_LAND+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y=element_text(size=13)
        # legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    
    
    ############################################################################
    #                             SHAPE:
    ############################################################################
    cat(paste0("\nShape:"), 
        file = paste0(dir.res_maps_amax_qnt,name_file_txt), 
        sep="\n", append=T)
    
    flag_signif   = T
    flag.limits   = T
    scale_grad    = list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)
    crs_map       = crswgs84
    v             = 4
    name_var      = "Change in shape per decade"
    type_estimate = c("MaxPost")
    name_model    = "SMEV Nonstationary"

    # initialize variable plot list:
    shapechange       = shapechange_LAND = c()
    a_w_MAP_df        = c()
    shape_change      = c()
    rast_a_w_map_LAND_nonan = c()
    
    # % of change of slope per decade wrw halfperiod year (a sort of average change per decade):
    shape_change = var_change_list$b_w_MAP_ns[[dur]]*10 /
                    (var_change_list$a_w_MAP_ns[[dur]]+var_change_list$b_w_MAP_ns[[dur]]*year_ref_id)*100
    # create dataframe for the map plot:
    a_w_MAP_df <- data.frame(x=lon_tmp, y=lat_tmp, z=shape_change)
    # trend_df$z[is.na(trend_df$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_NEW=a_w_MAP_df
    } else {
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3), 
                                  z=round(a_w_MAP_df$z,3))
    }
    rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
    
    # map:
    shapechange = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=a_w_MAP_df_NEW,aes(x=x,y=y,fill=z),alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z),alpha=1)+
      geom_sf(data=gg_world,fill=NA,color="grey70")+ #,linewidth=0.5) +
      geom_sf(data=test,colour="black",fill=NA,alpha=0.2)+ #,linewidth=0.5)+ 
      # geom_point(data =df_NEW_SIGNIF_EU,aes(x=x,y=y,size=z), 
      #            fill="black",shape=19,size=1)+  #pixels with signif. nonstationarity
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)
    
    if (flag_signif){
      shapechange=shapechange+
        geom_point(data=signif_ns_DIC_df_3, aes(x=x,y=y,size=z),  
                   fill="black",color="black",shape=21,size=0.4)
    }
    shapechange=shapechange+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y=element_text(size=13)
        # legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    # only land:
    ############
    # Replace nan with 1:
    a_w_MAP_df_LAND=a_w_MAP_df_NEW
    a_w_MAP_df_LAND$z[a_w_MAP_df_LAND$z==9999] =NA
    a_w_MAP_df_LAND_nonan <- a_w_MAP_df_LAND
    a_w_MAP_df_LAND_nonan$z[is.na(a_w_MAP_df_LAND_nonan$z)]=0
    # colMap <- colorRampPalette(c("yellow", "red"))(3)
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan
    } else {
      a_w_MAP_df_LAND_nonan_NEW =data.frame(x=round(a_w_MAP_df_LAND_nonan$x,3),
                                            y=round(a_w_MAP_df_LAND_nonan$y,3), 
                                            z=round(a_w_MAP_df_LAND_nonan$z,3))
    }
    rast_a_w_map_LAND_nonan =rasterFromXYZ(a_w_MAP_df_LAND_nonan_NEW, 
                                           crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2 <- crop(r,extent(mask_shp))
    r3 <- mask(r2,mask_shp)
    test_spdf <- as(r3,"SpatialPixelsDataFrame")
    test_df <- as.data.frame(test_spdf)
    
    shapechange_LAND = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=na.omit(test_df),aes(x=x,y=y,fill=z),alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
      geom_sf(data=gg_world,fill=NA,color="grey70")+  #,linewidth=0.5) +
      geom_sf(data=test,colour="black",fill=NA,alpha=0.2)+ #,linewidth=0.5)+ 
      # geom_point(data =df_NEW_SIGNIF_EU,aes(x=x,y=y,size=z), 
      #            fill="black",shape=19,size=1)+ # pixels with signif. nonstationarity
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)
    
    if (flag_signif){
      shapechange_LAND=shapechange_LAND+
        geom_point(data=test_df_signif_DIC3,aes(x=x,y=y,size=z),  
                   fill="black",color="black",shape=21,size=0.4)
    }
    shapechange_LAND=shapechange_LAND+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y=element_text(size=13)
        # legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    
    
    
    
    ############################################################################
    #                          N. ORD. EVENTS, n:
    ############################################################################
    cat(paste0("\nNumb.ordin.events,n:"),file=paste0(dir.res_maps_amax_qnt,name_file_txt), 
        sep="\n", append=T)
    
    flag_signif   = T
    significance  = significance
    name_var      = "n"
    type_estimate = c("sensslope")
    name_model    = "n, number of ordinary events per year"
    flag.limits   = T
    scale_grad    = list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)
    crs_map       = crswgs84
    v             = 5
    
    # initialize variable plot list:
    nchange       = nchange_LAND = c()
    a_w_MAP_df    = c()
    rast_a_w_map_LAND_nonan = c()
    
    # create dataframe for the map plot:
    a_w_MAP_df <- data.frame(x=lon_tmp,
                             y=lat_tmp, 
                             z=var_change_list$n_change[[dur]])
    # trend_df$z[is.na(trend_df$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_NEW=a_w_MAP_df
    } else {
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3), 
                                  z=round(a_w_MAP_df$z,3))
    }
    rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
    
    # signif "n":
    # signif_n_df  = signif_trend_n # same for each duration!
    if (sat_acron=='GSMAP'){
      signif_n_df <- data.frame(x=lon_tmp,
                                y=lat_tmp,
                                z=signif_trend_n,
                                d=significance$sensslope_n)
    } else {
      signif_n_df <- data.frame(x=round(lon_tmp,3),
                                y=round(lat_tmp,3),
                                z=signif_trend_n,
                                d=significance$sensslope_n)
    }
    signif_n_df$z[signif_n_df$z==F]=NA
    signif_n_df = na.omit(signif_n_df)
    signif_n_df = subset(signif_n_df, z != 0)
    
    # compute percentage respect total amount of pixels
    pixels_signif_n[[dur]] = length(signif_n_df$z)
    total_pixels = length(a_w_MAP_df$z)  # segm_df$z
    percent_signif_n[[dur]] = pixels_signif_n[[dur]]/total_pixels*100
    print(paste0("pixels with signif. trend in n (WHOLE) = ", pixels_signif_n[[dur]], 
                 " / ", total_pixels, " (", round(percent_signif_n[[dur]],1), "%)"))
    cat(paste0("pixels with signif. trend in n (WHOLE) = ", pixels_signif_n[[dur]], 
               " / ", total_pixels, " (", round(percent_signif_n[[dur]],1), "%)"), 
        file = paste0(dir.res_maps_amax_qnt, name_file_txt), sep="\n", append=T)
    
    nchange = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=a_w_MAP_df_NEW, aes(x=x, y=y, fill=z), alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
      geom_sf(data = gg_world, fill = NA, color = "grey70")+ #, linewidth=0.5) +
      geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+ #, linewidth=0.5)+
      coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
    if (flag_signif){
      nchange=nchange+
        geom_point(data=signif_n_df, aes(x=x, y=y, size=z),   fill="black", color="black", shape=21, size = 0.5)
      # geom_point(data=signif_ns_DIC_df_3, aes(x=x, y=y, size=z),   fill="black", color="black", shape=21, size = 0.4)
    }
    nchange=nchange+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y=element_text(size=13)
        # legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    
    # only land:
    ############
    # Replace nan with 1:
    a_w_MAP_df_LAND=a_w_MAP_df_NEW
    a_w_MAP_df_LAND$z[a_w_MAP_df_LAND$z==9999] =NA
    a_w_MAP_df_LAND_nonan <- a_w_MAP_df_LAND
    a_w_MAP_df_LAND_nonan$z[is.na(a_w_MAP_df_LAND_nonan$z)]=0
    # colMap <- colorRampPalette(c("yellow", "red"))(3)
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan
    } else {
      a_w_MAP_df_LAND_nonan_NEW =data.frame(x=round(a_w_MAP_df_LAND_nonan$x,3),
                                            y=round(a_w_MAP_df_LAND_nonan$y,3), 
                                            z=round(a_w_MAP_df_LAND_nonan$z,3))
    }
    rast_a_w_map_LAND_nonan =rasterFromXYZ(a_w_MAP_df_LAND_nonan_NEW, 
                                          crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2 <- crop(r, extent(mask_shp))
    r3 <- mask(r2, mask_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df <- as.data.frame(test_spdf)
    
    rast_signif_n_LAND = rasterFromXYZ(signif_n_df, crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_signif_n_LAND
    r2 <- crop(r, extent(mask_shp))
    r3 <- mask(r2, mask_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df_signif_n <- as.data.frame(test_spdf)
    
    # compute percentage respect total amount of pixels:
    pixels_signif_n_LAND[[dur]] = length(test_df_signif_n$z)
    total_pixels = length(test_df$z)  # segm_df$z
    percent_signif_n_LAND[[dur]] = pixels_signif_n_LAND[[dur]]/total_pixels*100
    print(paste0("pixels with signif. trend in n (LAND) = ", pixels_signif_n_LAND[[dur]], 
                 " / ", total_pixels, " (", round(percent_signif_n_LAND[[dur]],1), "%)"))
    cat(paste0("pixels with signif. trend in n (LAND) = ", pixels_signif_n_LAND[[dur]], 
               " / ", total_pixels, " (", round(percent_signif_n_LAND[[dur]],1), "%)"), 
        file = paste0(dir.res_maps_amax_qnt,name_file_txt), sep="\n", append=T)
    
    
    nchange_LAND = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=na.omit(test_df), aes(x=x, y=y, fill=z), alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
      geom_sf(data = gg_world, fill = NA, color = "grey70")+ #, linewidth=0.5) +
      geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+ #, linewidth=0.5)+
      coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
    if (flag_signif){
      nchange_LAND=nchange_LAND+
        geom_point(data=test_df_signif_n, aes(x=x, y=y, size=z), 
                   fill="black", color="black", shape=21, size = 0.5)
    }
    nchange_LAND=nchange_LAND+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y=element_text(size=13)
        # legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    
    ############################################################################
    # final plot titles:
    title0 <- ggdraw() +
      draw_label('',     fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25)+
      theme(plot.margin=margin(t=0,r=0,b=3,l=5))
    title1 <- ggdraw() +
      draw_label('AMAX', fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25)+
      theme(plot.margin=margin(t=0,r=0,b=3,l=5))
    title2 <- ggdraw() +
      draw_label('2yrs', fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25)+
      theme(plot.margin=margin(t=0,r=0,b=3,l=5))
    title3 <- ggdraw() +
      draw_label('scale',fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25)+
      theme(plot.margin=margin(t=0,r=0,b=3,l=5))
    title4 <- ggdraw() +
      draw_label('shape',fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25)+
      theme(plot.margin=margin(t=0,r=0,b=3,l=5))
    title5 <- ggdraw() +
      draw_label('n',    fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25)+
      theme(plot.margin=margin(t=0,r=0,b=3,l=5))
    all.titles = plot_grid(title0,title1,title2,title3,title4,title5,
                           labels=NULL,
                           ncol=1,nrow=6,
                           rel_heights=c(0.2,1,1,1,1,1),
                           align='hv')
    # merge final plot:
    all.plots = plot_grid( ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25),
                           AMAXchange  + theme(legend.title=element_blank(), 
                                               legend.position="none", 
                                               plot.margin=margin(t=0,r=0,b=3,l=0),
                                               axis.title=element_blank()),                       
                           twoyrchange + theme(legend.title=element_blank(), 
                                               legend.position="none", 
                                               plot.margin=margin(t=0,r=0,b=3,l=0),
                                               axis.title=element_blank()),
                           scalechange + theme(legend.title=element_blank(),
                                               legend.position="none", 
                                               plot.margin=margin(t=0,r=0,b=3,l=0),
                                               axis.title=element_blank()),
                           shapechange + theme(legend.title=element_blank(),
                                               legend.position="none", 
                                               plot.margin=margin(t=0,r=0,b=3,l=0),
                                               axis.title=element_blank()),
                           nchange     + theme(legend.title=element_blank(), 
                                               legend.position="none", 
                                               plot.margin=margin(t=0,r=0,b=3,l=0),
                                               axis.title=element_blank()),
                           # label_size=20,
                           labels=NULL, # c("","2yrs","scale","shape","n"),
                           ncol=1,nrow=6,
                           rel_heights=c(0.2,1,1,1,1,1),
                           align='hv')
    all.all.plots = plot_grid(  # create title object:
      all.titles,
      all.plots, # + theme(plot.title=element_text(size=20,hjust=0.5))+ggtitle(paste0(sat_acron,'\n')),
      label_size=20,
      labels=NULL,
      ncol=2,nrow=1,
      rel_widths=c(0.15,1),
      align='hv',
      hjust=-0.5,label_x=0)
    ggsave(all.all.plots, 
           filename=paste0(dir.res_maps_amax_qnt,'all_maps_',sat_prod,'_',durations[dur]/60,'h.png'),
           width=8,height=24,dpi=200,bg="white")
    
    
    # ONLY LAND:
    ############
    all.plots_LAND = plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25),
                               AMAXchange_LAND  + theme(legend.title=element_blank(), 
                                                        legend.position="none", 
                                                        plot.margin=margin(t=0,r=0,b=3,l=0),
                                                        axis.title=element_blank()),                       
                               twoyrchange_LAND + theme(legend.title=element_blank(), 
                                                        legend.position="none", 
                                                        plot.margin=margin(t=0,r=0,b=3,l=0),
                                                       axis.title=element_blank()),
                               scalechange_LAND + theme(legend.title=element_blank(), 
                                                        legend.position="none", 
                                                        plot.margin=margin(t=0,r=0,b=3,l=0),
                                                       axis.title=element_blank()),
                               shapechange_LAND + theme(legend.title=element_blank(), 
                                                        legend.position="none",
                                                        plot.margin=margin(t=0,r=0,b=3,l=0),
                                                       axis.title=element_blank()),
                               nchange_LAND     + theme(legend.title=element_blank(), 
                                                        legend.position="none", 
                                                        plot.margin=margin(t=0,r=0,b=3,l=0),
                                                        axis.title=element_blank()),
                               labels = NULL, # c("","2yrs","scale","shape","n"),
                               ncol=1,nrow=6,
                               rel_heights=c(0.2,1,1,1,1,1),
                               align='hv')
    all.all.plots_LAND=plot_grid(  # create title object:
      all.titles,
      all.plots_LAND, # + theme(plot.title = element_text(size= 20, hjust = 0.5)) + ggtitle("CMORPH\n"),
      label_size=20,
      labels=NULL,
      ncol=2,nrow=1,
      rel_widths=c(0.15,1),
      align='hv',
      hjust=-0.5,label_x=0)
    ggsave(all.all.plots_LAND,
           filename=paste0(dir.res_maps_amax_qnt,'all_maps_LAND_',sat_prod,'_', 
                           durations[dur]/60,'h.png'),
           width=8,height=24,dpi=200,bg="white")
    
    
    # Save lists of plots into .rds R file:
    saveRDS(list( AMAXchange           = AMAX_change[[dur]],
                  AMAXchange_plot      = AMAXchange,
                  AMAXchange_plot_LAND = AMAXchange_LAND,
                  twoyrchange          = quant_smev_ns_2yr_median_change[[dur]],
                  twoyrchange_plot     = twoyrchange,
                  twoyrchange_plot_LAND= twoyrchange_LAND,
                  scalechange          = scale_change[[1]][[dur]],
                  scalechange_plot     = scalechange,
                  scalechange_plot_LAND= scalechange_LAND,
                  shapechange          = shape_change[[1]][[dur]],
                  shapechange_plot     = shapechange,
                  shapechange_plot_LAND= shapechange_LAND,
                  nchange              = n_change[[dur]],
                  nchange_plot         = nchange,
                  nchange_plot_LAND    = nchange_LAND,
                  lon_tmp              = lon_tmp,
                  lat_tmp              = lat_tmp,
                  durat                = paste0(durations[dur]/60, 'h')), 
            file=paste0(dir.res_maps_amax_qnt,'all_maps_',sat_prod,'_',durations[dur]/60,'h.rds'))
 } # end loop on durations
  
  
  # Save lists of SIGNIF info into .rds R file:
  pixels_signif_list=list(percent_signif_AMAX     = percent_signif_AMAX,
                          pixels_signif_AMAX      = pixels_signif_AMAX,
                          pixels_signif_AMAX_LAND  = pixels_signif_AMAX_LAND,
                          percent_signif_AMAX     = percent_signif_AMAX,
                          percent_signif_AMAX_LAND = percent_signif_AMAX_LAND,
                          pixels_signif_DIC3      = pixels_signif_DIC3,
                          pixels_signif_DIC3_LAND  = pixels_signif_DIC3_LAND,
                          percent_signif_DIC3     = percent_signif_DIC3,
                          percent_signif_DIC3_LAND = percent_signif_DIC3_LAND,
                          pixels_signif_DIC1      = pixels_signif_DIC1,
                          pixels_signif_DIC1_LAND  = pixels_signif_DIC1_LAND,
                          percent_signif_DIC1     = percent_signif_DIC1,
                          percent_signif_DIC1_LAND = percent_signif_DIC1_LAND,
                          pixels_signif_n         = pixels_signif_n,
                          pixels_signif_n_LAND     = pixels_signif_n_LAND,
                          percent_signif_n        = percent_signif_n,
                          percent_signif_n_LAND    = percent_signif_n_LAND)
  saveRDS(pixels_signif_list,
          file=paste0(dir.res_maps_amax_qnt, 'signif_pixels_', sat_prod, '.rds'))
}


























################################################################################
# PLOT UNCERTAINTY of %changes of 2yrs quantile, shape, scale (NONSTAT SMEV)
################################################################################
plot_maps_changes_uncert=function(
                              var_unc_change_list = list(),   
                              lon_tmp        = NULL,
                              lat_tmp        = NULL,
                              durations      = c(60,180,360,1440), # in minutes
                              mask_shp       = NULL,
                              mask_world     = NULL,
                              flag_signif    = F,
                              significance   = NULL,
                              sat_acron      = 'TEST',
                              period         = '',
                              year_ref       = NA,
                              year_ref_id    = NA,
                              name_file_txt  = paste0('signif_pixels.txt'),
                              path_output    = "/map_test",
                              limits         = list(c(-20, 20)),
                              flag.limits    = T,
                              scale_grad     = list(),
                              crs_map        = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),
                              nbin           = 20,
                              bar.direction  = "horizontal",
                              legend.position= "bottom",
                              barwidth       = 25,
                              barheight      = 0.6){
################################################################################
  # initialize paths
  dir.res_maps_uncert =path_output
  dir.create(dir.res_maps_uncert)
  

  #  cycle on 4 rainfall durations
  for (dur in 1:4){
    message(paste0('DURATION = ', durations[dur]/60, 'h'))
    
    # stat.significant nonstat. pixels:
    if (sat_acron=='GSMAP'){
      signif_ns_DIC_df_3 <- data.frame(x=lon_tmp,
                                       y=lat_tmp,
                                       z=significance$DIC_diff3[[dur]])   
    } else {
      signif_ns_DIC_df_3 <- data.frame(x=round(lon_tmp,3),
                                       y=round(lat_tmp,3),
                                       z=significance$DIC_diff3[[dur]])
    }
    signif_ns_DIC_df_3$z[signif_ns_DIC_df_3$z==F]=NA
    signif_ns_DIC_df_3 = na.omit(signif_ns_DIC_df_3)
    signif_ns_DIC_df_3 = subset(signif_ns_DIC_df_3, z != 0)
    
    if (sat_acron=='GSMAP'){
      signif_ns_DIC_df_1 <- data.frame(x=lon_tmp,
                                       y=lat_tmp,
                                       z=significance$DIC_diff1[[dur]])   
    } else {
      signif_ns_DIC_df_1 <- data.frame(x=round(lon_tmp,3),
                                       y=round(lat_tmp,3),
                                       z=significance$DIC_diff1[[dur]])
    }
    signif_ns_DIC_df_1$z[signif_ns_DIC_df_1$z==F]=NA
    signif_ns_DIC_df_1 = na.omit(signif_ns_DIC_df_1)
    signif_ns_DIC_df_1 = subset(signif_ns_DIC_df_1, z != 0)
    
    
    ##############################################################################
    #                   Uncert in changes of 2yrs quantile:
    ##############################################################################
    flag_signif   = T
    flag.limits   = T
    scale_grad    = list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)
    crs_map       = crswgs84
    v             = 1
    name_var      = "Change in 2yrs quantiles per decade"
    type_estimate = c("MaxPost")
    name_model    = "SMEV Nonstationary"
    
    # initialize:
    a_w_MAP_df=c()
    rast_a_w_map_LAND_nonan=c()
    twoyrchange_UNCERT=twoyrchange_LAND_UNCERT=c()
    
    # create dataframe for the map plot:
    a_w_MAP_df <- data.frame(x=lon_tmp, 
                             y=lat_tmp, 
                             z=var_unc_change_list$unc_qnt2yr_change[[dur]])
    # eval(parse(text=(paste0('unc_qnt2yrs_', durations[dur]/60, 'h')))))
    
    # trend_df$z[is.na(trend_df$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_NEW=a_w_MAP_df
    } else {
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3), 
                                  z=round(a_w_MAP_df$z,3))
    }
    rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
    
    twoyrchange_UNCERT=ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=a_w_MAP_df_NEW,aes(x=x,y=y,fill=z),alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z),alpha=1)+
      geom_sf(data=gg_world, fill=NA,color="grey70")+ #,linewidth=0.5) +
      geom_sf(data=test,colour="black",fill=NA, alpha=0.2)+ #,linewidth=0.5)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)
    if (flag_signif){
      twoyrchange_UNCERT=twoyrchange_UNCERT+
        geom_point(data=signif_ns_DIC_df_3,aes(x=x,y=y,size=z),   
                   fill="black",color="black",shape=21,size=0.5)
    }
    twoyrchange_UNCERT=twoyrchange_UNCERT+
      theme(
         panel.grid.major=element_blank()
        ,axis.title.x=element_blank()
        ,axis.text.x=element_blank() # element_text(size=13)
        ,axis.text.y=element_text(size=13)
        ,legend.position="bottom"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_viridis(option="turbo",
                         limits=c(limits[[v]][1],limits[[v]][2]),
                         name='', # paste0(name_var,' (',type_estimate[v],'), ',name_model),
                         oob=scales::squish,
                         breaks=seq(limits[[v]][1],limits[[v]][2],limits[[v]][3]),
                         labels=seq(limits[[v]][1],limits[[v]][2],limits[[v]][3]),
                         guide=guide_colourbar(nbin=nbin,
                                                  raster=F,
                                                  display="rectangles",
                                                  direction=bar.direction,
                                                  barwidth=barwidth, #grid::unit(0.6, "npc"),
                                                  barheight=barheight, # grid::unit(1, "npc"), # - grid::unit(, "line"), #40,
                                                  frame.colour=c("NA"),
                                                  frame.linewidth=1))
    # guide  = guide_colourbar(nbin=20,
    #                          raster=T,
    #                          display = "rectangles",
    #                          direction=bar.direction,
    #                          barwidth=barwidth, #grid::unit(0.6, "npc"),
    #                          barheight=barheight, # grid::unit(1, "npc"), # - grid::unit(, "line"), #40,
    #                          frame.colour=c("NA"),
    #                          frame.linewidth=1))
    # scale_fill_gradient2(low = scale_grad[[v]]$low,
    #              high =scale_grad[[v]]$high,
    #              mid = scale_grad[[v]]$mid,
    #              midpoint = scale_grad[[v]]$midpoint,
    #              na.value = NA,
    #              name = paste0(name_var, '(',
    #                            type_estimate[v], '), ', name_model),
    #              limits = limits[[v]],
    #              oob    = scales::squish)
    
    # only land:
    ############
    # Replace nan with 1:
    a_w_MAP_df_LAND=a_w_MAP_df_NEW
    a_w_MAP_df_LAND$z[a_w_MAP_df_LAND$z==9999] =NA
    a_w_MAP_df_LAND_nonan <- a_w_MAP_df_LAND
    a_w_MAP_df_LAND_nonan$z[is.na(a_w_MAP_df_LAND_nonan$z)]=0
    # colMap <- colorRampPalette(c("yellow", "red"))(3)
    a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan
    # a_w_MAP_df_LAND_nonan_NEW = data.frame(x=round(a_w_MAP_df_LAND_nonan$x,3),
    #                                       y=round(a_w_MAP_df_LAND_nonan$y,3), 
    #                                       z=round(a_w_MAP_df_LAND_nonan$z,3))
    rast_a_w_map_LAND_nonan =rasterFromXYZ(a_w_MAP_df_LAND_nonan_NEW, 
                                           crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2 <- crop(r, extent(land_shp))
    r3 <- mask(r2, land_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df <- as.data.frame(test_spdf)
    
    rast_signif_DIC3_LAND = rasterFromXYZ(signif_ns_DIC_df_3, crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_signif_DIC3_LAND
    r2 <- crop(r, extent(land_shp))
    r3 <- mask(r2, land_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df_signif_DIC3 <- as.data.frame(test_spdf)
    
    rast_signif_DIC1_LAND = rasterFromXYZ(signif_ns_DIC_df_1, crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_signif_DIC1_LAND
    r2 <- crop(r, extent(land_shp))
    r3 <- mask(r2, land_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df_signif_DIC1 <- as.data.frame(test_spdf)
    
    
    twoyrchange_LAND_UNCERT = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=na.omit(test_df), aes(x=x, y=y, fill=z), alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
      geom_sf(data = gg_world, fill = NA, color = "grey70")+ #, linewidth=0.5) +
      geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+ #, linewidth=0.5)+
      coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
    if (flag_signif){
      twoyrchange_LAND_UNCERT=twoyrchange_LAND_UNCERT+
        geom_point(data=test_df_signif_DIC3, aes(x=x, y=y, size=z),  
                   fill="black", color="black", shape=21, size = 0.5) #0.4)
    }
    twoyrchange_LAND_UNCERT=twoyrchange_LAND_UNCERT+
      theme(
        panel.grid.major=element_blank()
        ,axis.title.x=element_blank()
        ,axis.text.x=element_blank() # element_text(size=13)
        ,axis.text.y=element_text(size=13)
        ,legend.position="bottom"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_viridis(option="turbo",
                         limits=c(limits[[v]][1],limits[[v]][2]),
                         name='', # paste0(name_var,' (',type_estimate[v],'), ',name_model),
                         oob=scales::squish,
                         breaks=seq(limits[[v]][1],limits[[v]][2],limits[[v]][3]),
                         labels=seq(limits[[v]][1],limits[[v]][2],limits[[v]][3]),
                         guide=guide_colourbar(nbin=nbin,
                                               raster=F,
                                               display="rectangles",
                                               direction=bar.direction,
                                               barwidth=barwidth, #grid::unit(0.6, "npc"),
                                               barheight=barheight, # grid::unit(1, "npc"), # - grid::unit(, "line"), #40,
                                               frame.colour=c("NA"),
                                               frame.linewidth=1))
    
    ############################################################################
    #                   uncert in changes of shape param SMEV
    ############################################################################
    flag_signif   = T
    flag.limits   = T
    scale_grad    = list(list(name="rdbl",low="red",high="blue",mid="white",midpoint=0))
    crs_map       = crswgs84
    v             = 2
    name_var      = "Change in shape per decade"
    type_estimate = c("MaxPost")
    name_model    = "SMEV Nonstationary"
    
    # initialize:
    a_w_MAP_df=c()
    rast_a_w_map_LAND_nonan=c()
    shapechange_UNCERT=shapechange_LAND_UNCERT=c()
    
    # create dataframe for the map plot:
    a_w_MAP_df <- data.frame(x=lon_tmp, 
                             y=lat_tmp, 
                             z=eval(parse(text=(paste0('unc_shape_', durations[dur]/60, 'h')))))
    # trend_df$z[is.na(trend_df$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_NEW=a_w_MAP_df
    } else {
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3), 
                                  z=round(a_w_MAP_df$z,3))
    }
    rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
    
    shapechange_UNCERT = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=a_w_MAP_df_NEW, aes(x=x,y=y,fill=z),alpha=1)+ 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z),alpha=1)+
      geom_sf(data=gg_world,fill=NA,color="grey70")+ #, linewidth=0.5)+
      geom_sf(data=test,colour="black",fill=NA,alpha=0.2)+ #, linewidth=0.5)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)
    if (flag_signif){
      shapechange_UNCERT=shapechange_UNCERT+
        geom_point(data=signif_ns_DIC_df_3,aes(x=x,y=y,size=z),   
                   fill="black",color="black",shape=21,size = 0.5)
    }
    shapechange_UNCERT=shapechange_UNCERT+
      theme(
        panel.grid.major=element_blank()
        ,axis.title.x=element_blank()
        ,axis.text.x=element_blank() # element_text(size=13)
        ,axis.text.y=element_text(size=13)
        ,legend.position="bottom"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_viridis(option="turbo",
                         limits=c(limits[[v]][1],limits[[v]][2]),
                         name='', # paste0(name_var,' (',type_estimate[v],'), ',name_model),
                         oob=scales::squish,
                         breaks=seq(limits[[v]][1],limits[[v]][2],limits[[v]][3]),
                         labels=seq(limits[[v]][1],limits[[v]][2],limits[[v]][3]),
                         guide=guide_colourbar(nbin=nbin,
                                               raster=F,
                                               display="rectangles",
                                               direction=bar.direction,
                                               barwidth=barwidth, #grid::unit(0.6, "npc"),
                                               barheight=barheight, # grid::unit(1, "npc"), # - grid::unit(, "line"), #40,
                                               frame.colour=c("NA"),
                                               frame.linewidth=1))
    # only land:
    ############
    # Replace nan with 1:
    a_w_MAP_df_LAND=a_w_MAP_df_NEW
    a_w_MAP_df_LAND$z[a_w_MAP_df_LAND$z==9999] =NA
    a_w_MAP_df_LAND_nonan <- a_w_MAP_df_LAND
    a_w_MAP_df_LAND_nonan$z[is.na(a_w_MAP_df_LAND_nonan$z)]=0
    # colMap <- colorRampPalette(c("yellow", "red"))(3)
    a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan
    # a_w_MAP_df_LAND_nonan_NEW = data.frame(x=round(a_w_MAP_df_LAND_nonan$x,3),
    #                                       y=round(a_w_MAP_df_LAND_nonan$y,3), 
    #                                       z=round(a_w_MAP_df_LAND_nonan$z,3))
    rast_a_w_map_LAND_nonan =rasterFromXYZ(a_w_MAP_df_LAND_nonan_NEW, 
                                           crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2 <- crop(r, extent(land_shp))
    r3 <- mask(r2, land_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df <- as.data.frame(test_spdf)
    
    rast_signif_DIC3_LAND = rasterFromXYZ(signif_ns_DIC_df_3, crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_signif_DIC3_LAND
    r2 <- crop(r, extent(land_shp))
    r3 <- mask(r2, land_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df_signif_DIC3 <- as.data.frame(test_spdf)
    
    rast_signif_DIC1_LAND = rasterFromXYZ(signif_ns_DIC_df_1, crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_signif_DIC1_LAND
    r2 <- crop(r, extent(land_shp))
    r3 <- mask(r2, land_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df_signif_DIC1 <- as.data.frame(test_spdf)
    
    
    shapechange_LAND_UNCERT = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=na.omit(test_df), aes(x=x, y=y, fill=z), alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
      geom_sf(data = gg_world, fill = NA, color = "grey70")+ #, linewidth=0.5) +
      geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+ #, linewidth=0.5)+
      coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
    if (flag_signif){
      shapechange_LAND_UNCERT=shapechange_LAND_UNCERT+
        geom_point(data=test_df_signif_DIC3, aes(x=x, y=y, size=z),  
                   fill="black", color="black", shape=21, size = 0.5) #0.4)
    }
    shapechange_LAND_UNCERT=shapechange_LAND_UNCERT+
      theme(
        panel.grid.major=element_blank()
        ,axis.title.x=element_blank()
        ,axis.text.x=element_blank() # element_text(size=13)
        ,axis.text.y=element_text(size=13)
        ,legend.position="bottom"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_viridis(option="turbo",
                         limits=c(limits[[v]][1],limits[[v]][2]),
                         name='', # paste0(name_var,' (',type_estimate[v],'), ',name_model),
                         oob=scales::squish,
                         breaks=seq(limits[[v]][1],limits[[v]][2],limits[[v]][3]),
                         labels=seq(limits[[v]][1],limits[[v]][2],limits[[v]][3]),
                         guide=guide_colourbar(nbin=nbin,
                                               raster=F,
                                               display="rectangles",
                                               direction=bar.direction,
                                               barwidth=barwidth, #grid::unit(0.6, "npc"),
                                               barheight=barheight, # grid::unit(1, "npc"), # - grid::unit(, "line"), #40,
                                               frame.colour=c("NA"),
                                               frame.linewidth=1))
    

    ########################################################################
    #                   uncert in changes of scale param 
    ########################################################################
    flag_signif   = T
    flag.limits   = T
    scale_grad    = list(list(name="rdbl",low="red",high="blue",mid="white",midpoint=0))
    crs_map       = crswgs84
    v             = 3
    name_var      = "Change in Scale per decade"
    type_estimate = c("MaxPost")
    name_model    = "SMEV Nonstationary"
    
    # initialize:
    a_w_MAP_df=c()
    rast_a_w_map_LAND_nonan=c()
    scalechange_UNCERT=scalechange_LAND_UNCERT=c()
    
    # create dataframe for the map plot:
    a_w_MAP_df <- data.frame(x=lon_tmp,
                             y=lat_tmp,
                             z=var_unc_change_list$unc_scale_change[[dur]])
    # trend_df$z[is.na(trend_df$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_NEW=a_w_MAP_df
    } else {
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3), 
                                  z=round(a_w_MAP_df$z,3))
    }
    rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
    
    scalechange_UNCERT = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=a_w_MAP_df_NEW, aes(x=x, y=y, fill=z), alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
      geom_sf(data = gg_world, fill = NA, color = "grey70")+ #, linewidth=0.5) +
      geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+ #, linewidth=0.5)+
      coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
    if (flag_signif){
      scalechange_UNCERT=scalechange_UNCERT+
        geom_point(data=signif_ns_DIC_df_3, aes(x=x, y=y, size=z),   
                   fill="black", color="black", shape=21, size = 0.5)
    }
    scalechange_UNCERT=scalechange_UNCERT+
      theme(
        panel.grid.major=element_blank()
        #,axis.title.x=element_blank()
        #,axis.text.x=element_blank() # element_text(size=13)
        ,axis.text.y=element_text(size=13)
        ,legend.position="bottom"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_viridis(option="turbo",
                         limits=c(limits[[v]][1],limits[[v]][2]),
                         name='', # paste0(name_var,' (',type_estimate[v],'), ',name_model),
                         oob=scales::squish,
                         breaks=seq(limits[[v]][1],limits[[v]][2],limits[[v]][3]),
                         labels=seq(limits[[v]][1],limits[[v]][2],limits[[v]][3]),
                         guide=guide_colourbar(nbin=nbin,
                                               raster=F,
                                               display="rectangles",
                                               direction=bar.direction,
                                               barwidth=barwidth, #grid::unit(0.6, "npc"),
                                               barheight=barheight, # grid::unit(1, "npc"), # - grid::unit(, "line"), #40,
                                               frame.colour=c("NA"),
                                               frame.linewidth=1))
    # only land:
    ############
    # Replace nan with 1:
    a_w_MAP_df_LAND=a_w_MAP_df_NEW
    a_w_MAP_df_LAND$z[a_w_MAP_df_LAND$z==9999] =NA
    a_w_MAP_df_LAND_nonan <- a_w_MAP_df_LAND
    a_w_MAP_df_LAND_nonan$z[is.na(a_w_MAP_df_LAND_nonan$z)]=0
    # colMap <- colorRampPalette(c("yellow", "red"))(3)
    a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan
    # a_w_MAP_df_LAND_nonan_NEW = data.frame(x=round(a_w_MAP_df_LAND_nonan$x,3),
    #                                       y=round(a_w_MAP_df_LAND_nonan$y,3), 
    #                                       z=round(a_w_MAP_df_LAND_nonan$z,3))
    rast_a_w_map_LAND_nonan =rasterFromXYZ(a_w_MAP_df_LAND_nonan_NEW, 
                                           crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2 <- crop(r, extent(land_shp))
    r3 <- mask(r2, land_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df <- as.data.frame(test_spdf)
    
    rast_signif_DIC3_LAND = rasterFromXYZ(signif_ns_DIC_df_3, crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_signif_DIC3_LAND
    r2 <- crop(r, extent(land_shp))
    r3 <- mask(r2, land_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df_signif_DIC3 <- as.data.frame(test_spdf)
    
    rast_signif_DIC1_LAND = rasterFromXYZ(signif_ns_DIC_df_1, crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_signif_DIC1_LAND
    r2 <- crop(r, extent(land_shp))
    r3 <- mask(r2, land_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df_signif_DIC1 <- as.data.frame(test_spdf)
    
    
    scalechange_LAND_UNCERT=ggplot()+ 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=na.omit(test_df), aes(x=x, y=y, fill=z), alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z),alpha=1)+
      geom_sf(data=gg_world,fill=NA,color="grey70")+ #,linewidth=0.5) +
      geom_sf(data=test,colour="black",fill=NA,alpha=0.2)+ #,linewidth=0.5)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)
    if (flag_signif){
      scalechange_LAND_UNCERT=scalechange_LAND_UNCERT+
        geom_point(data=test_df_signif_DIC3, aes(x=x,y=y,size=z),  
                   fill="black",color="black",shape=21,size=0.5) #0.4)
    }
    scalechange_LAND_UNCERT=scalechange_LAND_UNCERT+
      theme(
        panel.grid.major=element_blank()
        # ,axis.title.x=element_blank()
        # ,axis.text.x=element_blank() # element_text(size=13)
        ,axis.text.y=element_text(size=13)
        ,legend.position="bottom"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_viridis(option="turbo",
                         limits=c(limits[[v]][1],limits[[v]][2]),
                         name='', # paste0(name_var,' (',type_estimate[v],'), ',name_model),
                         oob=scales::squish,
                         breaks=seq(limits[[v]][1],limits[[v]][2],limits[[v]][3]),
                         labels=seq(limits[[v]][1],limits[[v]][2],limits[[v]][3]),
                         guide=guide_colourbar(nbin=nbin,
                                               raster=F,
                                               display="rectangles",
                                               direction=bar.direction,
                                               barwidth=barwidth, #g rid::unit(0.6, "npc"),
                                               barheight=barheight, # grid::unit(1, "npc"), # - grid::unit(, "line"), #40,
                                               frame.colour=c("NA"),
                                               frame.linewidth=1))
    
    # ##########################################################################
    title0 <- ggdraw()+
      draw_label('',                   fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25)+
      theme(plot.margin=margin(t=0,r=0,b=3,l=5))
    title1 <- ggdraw()+
      draw_label('2yrs\n return level',fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25)+
      theme(plot.margin=margin(t=0,r=0,b=3,l=5))
    title2 <- ggdraw()+
      draw_label('Shape',              fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25)+
      theme(plot.margin=margin(t=0,r=0,b=3,l=5))
    title3 <- ggdraw()+
      draw_label('Scale',              fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25)+
      theme(plot.margin = margin(t=0,r=0,b=3,l=5))
    all.titles = plot_grid(title0,title1,title2,title3,
                           labels=NULL,
                           ncol=1,nrow=4,
                           rel_heights=c(0.25,1,1,1),
                           align='hv')
    # merge final plot:
    all.plots_ALLPIX=plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25),
                               twoyrchange_UNCERT + theme(legend.title=element_blank(),
                                                          # legend.position="none", 
                                                          plot.margin=margin(t=0,r=0,b=20,l=0),
                                                          axis.title=element_blank()),
                               shapechange_UNCERT + theme(legend.title=element_blank(), 
                                                          # legend.position="none",
                                                          plot.margin=margin(t=5,r=0,b=20,l=0),
                                                          axis.title=element_blank()),
                               scalechange_UNCERT + theme(legend.title=element_blank(), 
                                                          # legend.position="none", 
                                                          plot.margin=margin(t=5,r=0,b=20,l=0),
                                                          axis.title=element_blank()),
                               labels=NULL,
                               ncol=1,nrow=4,
                               rel_heights = c(0.2,1,1,1),
                               align='hv')
    all.all.uncert.plots_ALLPIX=plot_grid(
      all.titles,
      all.plots_ALLPIX, # +theme(plot.title=element_text(size=20,hjust=0.5))+ggtitle("CMORPH\n"),
      label_size=20,
      labels=NULL,
      ncol=2,nrow=1,
      rel_widths=c(0.25,1),
      align='hv',
      hjust=-0.5,label_x=0)
    ggsave(all.all.uncert.plots_ALLPIX, 
           filename=paste0(dir.res_maps_uncert,'/all_maps_uncert_',
                           sat_prod,'_ALLPIX_',durations[dur]/60,'h.png'), 
           width=8,height=24,dpi=200,bg="white")
    
    
    # Only land pixels:
    ###################
    all.plots_LAND=plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25),
                             twoyrchange_LAND_UNCERT + theme(legend.title= element_blank(),
                                                             # legend.position="none", 
                                                               plot.margin=margin(t=0,r=0,b=20,l=0),
                                                               axis.title=element_blank()),
                             shapechange_LAND_UNCERT + theme(legend.title=element_blank(), 
                                                               # legend.position="none",
                                                               plot.margin=margin(t=5,r=0,b=20,l=0),
                                                               axis.title=element_blank()),
                             scalechange_LAND_UNCERT + theme(legend.title=element_blank(), 
                                                               #legend.position="none", 
                                                               plot.margin=margin(t=5,r=0,b=20,l=0),
                                                               axis.title=element_blank()),
                             labels=NULL,
                             ncol=1,nrow=4,
                             rel_heights = c(0.2,1,1,1),
                             align='hv')
    all.all.uncert.plots_LAND= plot_grid(  # create title object:
      all.titles,
      all.plots_LAND, # +theme(plot.title=element_text(size=20,hjust=0.5))+ggtitle("CMORPH\n"),
      label_size=22,
      labels=NULL,
      ncol=2,nrow=1,
      rel_widths=c(0.2,1),
      align='hv',
      hjust=-0.5,label_x=0)
    ggsave(all.all.uncert.plots_LAND, 
           filename=paste0(dir.res_maps_uncert,'/all_maps_uncert_',sat_prod,
                           '_LAND_',durations[dur]/60,'h.png'), 
           width=8,height=24,dpi=200,bg="white")
    
    
    # Save results and plots to .rds R file:
    saveRDS(list( twoyrchange_UNCERT      = twoyrchange_UNCERT,
                  twoyrchange_LAND_UNCERT = twoyrchange_LAND_UNCERT,
                  shapechange_UNCERT      = shapechange_UNCERT,
                  shapechange_LAND_UNCERT = shapechange_LAND_UNCERT,
                  scalechange_UNCERT      = scalechange_UNCERT,
                  scalechange_LAND_UNCERT = scalechange_LAND_UNCERT,
                  lon_tmp                 = lon_tmp,
                  lat_tmp                 = lat_tmp,
                  durat                   = paste0(durations[dur]/60, 'h')), 
            file=paste0(dir.res_maps_uncert,'/all_maps_uncert_',sat_prod,'_',
                        durations[dur]/60,'h.rds'))
  } # end loop on durations
  ##############################################################################
}

























################################################################################
# PLOT %changes of 2-20-50yrs quantiles (NONSTAT SMEV)
################################################################################
plot_maps_changes_qnt=function(
                        var_qnt_change_list = list(),   
                        lon_tmp        = NULL,
                        lat_tmp        = NULL,
                        durations      = c(60,180,360,1440), # in minutes
                        mask_shp       = NULL,
                        mask_world     = NULL,
                        flag_signif    = F,
                        significance   = NULL,
                        sat_acron      = 'TEST',
                        period         = '',
                        year_ref       = NA,
                        year_ref_id    = NA,
                        name_file_txt  = paste0('signif_pixels.txt'),
                        path_output    = "/map_test",
                        limits         = list(c(-20, 20)),
                        limits_mean    = list(c(-15, 15)),
                        flag.limits    = T,
                        scale_grad     = list(),
                        crs_map        = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),
                        nbin           = 20,
                        bar.direction  = "horizontal",
                        legend.position= "bottom",
                        barwidth       = 25,
                        barheight      = 0.6,
                        size_text_plot = 6){
################################################################################
  # initialize paths
  dir.res_maps_qnt =path_output
  dir.create(dir.res_maps_qnt)

  
  for (dur in 1:4){
    message(paste0('DURATION = ', durations[dur]/60, 'h'))
    
    # nonstationary model is more likely than stationary (DIC)
    if (sat_acron=='GSMAP'){
      signif_ns_DIC_df_3 <- data.frame(x=lon_tmp,
                                       y=lat_tmp,
                                       z=significance$DIC_diff3[[dur]])   
    } else {
      signif_ns_DIC_df_3 <- data.frame(x=round(lon_tmp,3),
                                       y=round(lat_tmp,3),
                                       z=round(significance$DIC_diff3[[dur]],3))
    }
    signif_ns_DIC_df_3$z[signif_ns_DIC_df_3$z==F]=NA
    signif_ns_DIC_df_3 = na.omit(signif_ns_DIC_df_3)
    signif_ns_DIC_df_3 = subset(signif_ns_DIC_df_3, z != 0)      
    
    if (sat_acron=='GSMAP'){
      signif_ns_DIC_df_1 <- data.frame(x=lon_tmp,
                                       y=lat_tmp,
                                       z=significance$DIC_diff1[[dur]])   
    } else {
      signif_ns_DIC_df_1 <- data.frame(x=round(lon_tmp,3),
                                       y=round(lat_tmp,3),
                                       z=round(significance$DIC_diff1[[dur]],3))
    }
    signif_ns_DIC_df_1$z[signif_ns_DIC_df_1$z==F]=NA
    signif_ns_DIC_df_1 = na.omit(signif_ns_DIC_df_1)
    signif_ns_DIC_df_1 = subset(signif_ns_DIC_df_1, z != 0)  
    
  
    
    ############################################################################
    #                        2yrs quantile:
    ############################################################################
    flag_signif   = T
    flag.limits   = T
    scale_grad    = list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)
    crs_map       = crswgs84
    v             = 1
    name_var      = "Change in 2yrs quantiles per decade"
    type_estimate = c("MaxPost")
    name_model    = "SMEV Nonstationary"
    
    # initialize:
    a_w_MAP_df=c()
    rast_a_w_map_LAND_nonan=c()
    twoyrchange=twoyrchange_LAND=c()
    twoyrchange_LAND_mean=twoyrchange_LAND_mean_signif=c()
    
    # create dataframe for the map plot:
    a_w_MAP_df <- data.frame(x=lon_tmp, 
                             y=lat_tmp, 
                             z=var_qnt_change_list$quant_smev_ns_2yr_change[[dur]])
    # trend_df$z[is.na(trend_df$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_NEW=a_w_MAP_df   
    } else {
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3),
                                  z=round(a_w_MAP_df$z,3))
    }
    # rast_a_w_MAP=rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
    
    twoyrchange=ggplot()+
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=a_w_MAP_df_NEW,aes(x=x,y=y,fill=z),alpha=1)+ 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
      geom_sf(data=gg_world,fill=NA,color="grey70")+ #,linewidth=0.5) +
      geom_sf(data=test,colour="black",fill=NA,alpha=0.2)+ #,linewidth=0.5)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)
    if (flag_signif){
      twoyrchange=twoyrchange+
        # geom_point(data=signif_ns_DIC_df_1, 
        #            aes(x=x, y=y, size=z),fill="black",color="black",shape=21,size=0.7)
        geom_point(data=signif_ns_DIC_df_3,aes(x=x,y=y,size=z),  
                   fill="black",color="black",shape=21,size=0.4)
    }
    twoyrchange=twoyrchange+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        #legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    # only land:
    ####################
    # Replace nan with 1:
    a_w_MAP_df_LAND=a_w_MAP_df_NEW
    a_w_MAP_df_LAND$z[a_w_MAP_df_LAND$z==9999] =NA
    a_w_MAP_df_LAND_nonan <- a_w_MAP_df_LAND
    a_w_MAP_df_LAND_nonan$z[is.na(a_w_MAP_df_LAND_nonan$z)]=0
    # colMap <- colorRampPalette(c("yellow", "red"))(3)
    # a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan   
    } else {
      a_w_MAP_df_LAND_nonan_NEW = data.frame(x=round(a_w_MAP_df_LAND_nonan$x,3),
                                            y=round(a_w_MAP_df_LAND_nonan$y,3),
                                            z=round(a_w_MAP_df_LAND_nonan$z,3))
    }
    rast_a_w_map_LAND_nonan =rasterFromXYZ(a_w_MAP_df_LAND_nonan_NEW, 
                                          crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2 <- crop(r, extent(mask_shp))
    r3 <- mask(r2, mask_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df <- as.data.frame(test_spdf)
    
    rast_signif_DIC3_LAND = rasterFromXYZ(signif_ns_DIC_df_3, crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_signif_DIC3_LAND
    r2 <- crop(r, extent(mask_shp))
    r3 <- mask(r2, land_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df_signif_DIC3 <- as.data.frame(test_spdf)
    
    rast_signif_DIC1_LAND = rasterFromXYZ(signif_ns_DIC_df_1, crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_signif_DIC1_LAND
    r2 <- crop(r, extent(mask_shp))
    r3 <- mask(r2, mask_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df_signif_DIC1 <- as.data.frame(test_spdf)
    
    
    twoyrchange_LAND = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=na.omit(test_df), aes(x=x, y=y, fill=z), alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
      geom_sf(data = gg_world, fill = NA, color = "grey70")+ #, linewidth=0.5) +
      geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+ #, linewidth=0.5)+
      coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
    if (flag_signif){
      twoyrchange_LAND=twoyrchange_LAND+
        # geom_point(data=test_df_signif_DIC1, aes(x=x, y=y, size=z),   fill="black", color="black", shape=21, size = 0.7)
        geom_point(data=test_df_signif_DIC3, aes(x=x, y=y, size=z),  
                   fill="black", color="black", shape=21, size = 0.4)
    }
    twoyrchange_LAND=twoyrchange_LAND+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        #legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    
    # Aggregate to administrative regions:
    ######################################
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2<-crop(r,extent(land_shp))
    r3<-mask(r2,land_shp)
    
    mavg<-terra::extract(r3,test,weights=T,fun=mean,na.rm=T)
    mavg_1=unlist(lapply(mavg, function(x) if (!is.null(x)) mean(x,na.rm=T) else NA))
    output=data.frame(mavg)
    # print(output)
    test$z=mavg
    
    mavg<-terra::extract(r3,land_shp,weights=T,fun=mean,na.rm=T)
    mavg_1=unlist(lapply(mavg, function(x) if (!is.null(x)) mean(x,na.rm=T) else NA))
    output=data.frame(mavg)
    # print(output)
    land_shp_NEW=land_shp
    land_shp_NEW$z= mavg
    
    twoyrchange_LAND_mean=ggplot()+
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_sf(data=land_shp_NEW,colour="black",alpha=1)+
      geom_sf(data=test,aes(fill=z),colour="black",alpha=1)+
      geom_sf_text(
        data=test,aes(label=round(z,0)),
        stat="sf_coordinates", position = "nudge",
        parse=F,check_overlap=F,
        na.rm=F,show.legend=NA,inherit.aes=T,
        fun.geometry = NULL,size=size_text_plot)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y = element_text(size=13)
        #legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name="Mean changes per decade of \n2yrs return levels [%]",
                           limits=limits_mean[[v]],
                           oob=scales::squish)
    
    # Get % of signif pixels for each region:
    ##############################################
    if (sat_acron=='GSMAP'){
      signif_ns_DIC_df_3 <- data.frame(x=lon_tmp,
                                       y=lat_tmp, 
                                       z=significance$DIC_diff3[[dur]])  
    } else {
      signif_ns_DIC_df_3 <- data.frame(x=round(lon_tmp,3),
                                       y=round(lat_tmp,3), 
                                       z=significance$DIC_diff3[[dur]])
    }
    r=r2=r3=c()
    r=rasterFromXYZ(signif_ns_DIC_df_3,crs=crs_map) 
    r2<-crop(r,extent(test))
    r3<-mask(r2,test)
    signif_reg=terra::extract(r,test,weights=F,fun=sum,na.rm=T)
    signif_reg_df=data.frame(reg_name=test$reg_name,
                             signif_pixels_1=signif_reg)
    signif_reg_df$pix_reg_1=unlist(lapply(terra::extract(r,test,weights=F,na.rm=T,cells=T),length))
    # percent of pixels:
    signif_reg_df$percent_signif_1=round(signif_reg_df$signif_pixels_1/signif_reg_df$pix_reg_1*100,2)
    # number of signif pixels per region (accounting for coverage fraction)
    signif_reg_df$signif_pixels_2=ceil(exactextractr::exact_extract(r3,test,
                                      function(values,coverage_fraction) 
                                       sum(values*coverage_fraction,
                                       na.rm=T),progress=F))
    # Number of pixels per region:
    signif_reg_df$pix_reg_2=ceil(exactextractr::exact_extract(r3,test, 
                                  function(values, coverage_fraction) 
                                    sum(coverage_fraction),progress=F))
    # Percent:
    signif_reg_df$percent_signif_2=round(signif_reg_df$signif_pixels_2/signif_reg_df$pix_reg_2*100,2)
    test$signif=signif_reg_df$percent_signif_2
    # # Number of distinct raster values within the polygon (coverage fractions are ignored)
    # signif_reg_df$dist_val = exact_extract(r, test, function(values, coverage_fraction) length(unique(values)))
    # signif_reg_df
    # Number of distinct raster values in cells more than 10% covered by the polygon
    # signif_reg = exact_extract(r, test, function(values, coverage_fraction) length(unique(values[coverage_fraction > 0.1])))
    
    twoyrchange_LAND_mean_signif = ggplot()+
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_sf(data=land_shp_NEW, colour="black",alpha=1)+
      geom_sf(data=test,aes(fill=z),colour="black",alpha=1)+
      # geom_point(data=test_df_signif_DIC3,aes(x=x,y=y,size=z),  
      #            fill="black", color="black",shape=21,size=1)+
      geom_sf_text(
        data=test,aes(label=paste0(round(signif,0),'%')), 
        stat="sf_coordinates",position="nudge",
        parse=F,check_overlap=F,
        na.rm=F,show.legend=NA,inherit.aes=T,
        fun.geometry = NULL,size=size_text_plot)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y = element_text(size=13)
        #legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name="Mean changes per decade of \n2yrs return levels [%]",
                           limits=limits_mean[[v]],
                           oob=scales::squish)

    
    
    
    ############################################################################
    #                 10yrs quantile (computed but not plotted):
    ############################################################################
    flag_signif   = T
    flag.limits   = T
    scale_grad    = list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)
    crs_map       = crswgs84
    v             = 2
    name_var      = "Change in 10yrs quantiles per decade"
    type_estimate = c("MaxPost")
    name_model    = "SMEV Nonstationary"

    # initialize:
    a_w_MAP_df=c()
    rast_a_w_map_LAND_nonan=c()
    tenyrchange=tenyrchange_LAND=c()
    tenyrchange_LAND_mean=tenyrchange_LAND_mean_signif=c()
    
    # create dataframe for the map plot:
    a_w_MAP_df <- data.frame(x=lon_tmp, 
                             y=lat_tmp, 
                             z=var_qnt_change_list$quant_smev_ns_10yr_change[[dur]])
    # trend_df$z[is.na(trend_df$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_NEW=a_w_MAP_df 
    } else {
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3),
                                  z=round(a_w_MAP_df$z,3))
    }
    rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)

    signif_ns_DIC_df_3 <- data.frame(x=lon_tmp,
                                     y=lat_tmp, 
                                     z=significance$DIC_diff3[[dur]])
    signif_ns_DIC_df_3$z[signif_ns_DIC_df_3$z==F]=NA
    signif_ns_DIC_df_3 = na.omit(signif_ns_DIC_df_3)
    signif_ns_DIC_df_3 = subset(signif_ns_DIC_df_3, z != 0)     
    
    tenyrchange = ggplot() +
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)",
                         breaks=c(6.00001,10,14,18),
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",
                         breaks=c(36.00001,40,44,47.9999),
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=a_w_MAP_df_NEW,aes(x=x, y=y,fill=z),alpha=1)+
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z),alpha=1)+
      geom_sf(data=gg_world,fill=NA,color="grey70")+ #,linewidth=0.5) +
      geom_sf(data=test,colour="black",fill=NA,alpha=0.2)+ #,linewidth=0.5)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)
    if (flag_signif){
      tenyrchange=tenyrchange+
        #geom_point(data=signif_ns_DIC_df_1,aes(x=x,y=y,size=z),
        #           fill="black",color="black",shape=21,size=0.7)
        geom_point(data=signif_ns_DIC_df_3,aes(x=x,y=y,size=z),
                   fill="black",color="black",shape=21,size=0.4)
    }
    tenyrchange=tenyrchange+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        #legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    # only land:
    ############
    # Replace nan with 1:
    a_w_MAP_df_LAND=a_w_MAP_df_NEW
    a_w_MAP_df_LAND$z[a_w_MAP_df_LAND$z==9999] =NA
    a_w_MAP_df_LAND_nonan <- a_w_MAP_df_LAND
    a_w_MAP_df_LAND_nonan$z[is.na(a_w_MAP_df_LAND_nonan$z)]=0
    # colMap <- colorRampPalette(c("yellow", "red"))(3)
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan 
    } else {
      a_w_MAP_df_LAND_nonan_NEW = data.frame(x=round(a_w_MAP_df_LAND_nonan$x,3),
                                            y=round(a_w_MAP_df_LAND_nonan$y,3),
                                            z=round(a_w_MAP_df_LAND_nonan$z,3))
    }
    rast_a_w_map_LAND_nonan =rasterFromXYZ(a_w_MAP_df_LAND_nonan_NEW,
                                          crs=crs_map)
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2<-crop(r,extent(land_shp))
    r3<-mask(r2,land_shp)
    test_spdf<-as(r3,"SpatialPixelsDataFrame")
    test_df<-as.data.frame(test_spdf)

    tenyrchange_LAND=ggplot()+
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)",
                         breaks=c(6.00001,10,14,18),
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",
                         breaks=c(36.00001,40,44,47.9999),
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=na.omit(test_df),aes(x=x,y=y,fill=z),alpha=1) +
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z),alpha=1)+
      geom_sf(data=gg_world,fill=NA,color="grey70")+ #,linewidth=0.5) +
      geom_sf(data=test,colour="black",fill=NA,alpha=0.2)+ #,linewidth=0.5)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)
    if (flag_signif){
      tenyrchange_LAND=tenyrchange_LAND+
        # geom_point(data=test_df_signif_DIC1,aes(x=x,y=y,size=z),
        #            fill="black",color="black",shape=21,size=0.7)
        geom_point(data=test_df_signif_DIC3, aes(x=x,y=y,size=z),
                   fill="black",color="black",shape=21,size=0.4)
    }
    tenyrchange_LAND=tenyrchange_LAND+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        #legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)

    # Aggregate to administrative regions:
    ######################################
    r=r2=r3=c()
    signif_reg_df=c()
    r=rast_a_w_map_LAND_nonan
    r2<-crop(r,extent(land_shp))
    r3<-mask(r2,land_shp)
    mavg<-terra::extract(r3,test,weights=T,fun=mean,na.rm=T)
    mavg_1=unlist(lapply(mavg, function(x) if (!is.null(x)) mean(x,na.rm=T) else NA))
    output=data.frame(mavg)
    # print(output)
    test$z=mavg
    
    mavg<-terra::extract(r3,land_shp,weights=T,fun=mean,na.rm=T)
    mavg_1=unlist(lapply(mavg, function(x) if (!is.null(x)) mean(x,na.rm=T) else NA))
    output=data.frame(mavg)
    # print(output)
    land_shp_NEW=land_shp
    land_shp_NEW$z= mavg
    
    tenyrchange_LAND_mean=ggplot()+
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_sf(data=land_shp_NEW,colour="black",alpha=1)+
      geom_sf(data=test,aes(fill=z),colour="black",alpha=1)+
      geom_sf_text(
        data=test,aes(label=round(z,0)),
        stat="sf_coordinates", position = "nudge",
        parse=F,check_overlap=F,
        na.rm=F,show.legend=NA,inherit.aes=T,
        fun.geometry = NULL,size=size_text_plot)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y = element_text(size=13)
        #legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name="Mean changes per decade of \n10yrs return levels [%]",
                           limits=limits_mean[[v]],
                           oob=scales::squish)
    
    # Get % of signif pixels for each region:
    ##############################################
    if (sat_acron=='GSMAP'){
      signif_ns_DIC_df_3 <- data.frame(x=lon_tmp,
                                       y=lat_tmp, 
                                       z=significance$DIC_diff3[[dur]])  
    } else {
      signif_ns_DIC_df_3 <- data.frame(x=round(lon_tmp,3),
                                       y=round(lat_tmp,3), 
                                       z=significance$DIC_diff3[[dur]])
    }
    r=r2=r3=c()
    signif_reg_df=c()
    r=rasterFromXYZ(signif_ns_DIC_df_3,crs=crs_map) 
    r2<-crop(r,extent(test))
    r3<-mask(r2,test)
    signif_reg=terra::extract(r,test,weights=F,fun=sum,na.rm=T)
    signif_reg_df=data.frame(reg_name=test$reg_name,
                             signif_pixels_1=signif_reg)
    signif_reg_df$pix_reg_1=unlist(lapply(terra::extract(r,test,weights=F,na.rm=T,cells=T),length))
    # percent of pixels:
    signif_reg_df$percent_signif_1=round(signif_reg_df$signif_pixels_1/signif_reg_df$pix_reg_1*100,2)
    # number of signif pixels per region (accounting for coverage fraction)
    signif_reg_df$signif_pixels_2=ceil(exactextractr::exact_extract(r3,test,
                                        function(values,coverage_fraction) 
                                          sum(values*coverage_fraction,
                                              na.rm=T),progress=F))
    # Number of pixels per region:
    signif_reg_df$pix_reg_2=ceil(exactextractr::exact_extract(r3,test, 
                                      function(values, coverage_fraction) 
                                        sum(coverage_fraction),progress=F))
    # Percent:
    signif_reg_df$percent_signif_2=round(signif_reg_df$signif_pixels_2/signif_reg_df$pix_reg_2*100,2)
    test$signif=signif_reg_df$percent_signif_2
    # # Number of distinct raster values within the polygon (coverage fractions are ignored)
    # signif_reg_df$dist_val = exact_extract(r, test, function(values, coverage_fraction) length(unique(values)))
    # signif_reg_df
    # Number of distinct raster values in cells more than 10% covered by the polygon
    # signif_reg = exact_extract(r, test, function(values, coverage_fraction) length(unique(values[coverage_fraction > 0.1])))
    
    tenyrchange_LAND_mean_signif = ggplot()+
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_sf(data=land_shp_NEW, colour="black",alpha=1)+
      geom_sf(data=test,aes(fill=z),colour="black",alpha=1)+
      # geom_point(data=test_df_signif_DIC3,aes(x=x,y=y,size=z),  
      #            fill="black", color="black",shape=21,size=1)+
      geom_sf_text(
        data=test,aes(label=paste0(round(signif,0),'%')), 
        stat="sf_coordinates",position="nudge",
        parse=F,check_overlap=F,
        na.rm=F,show.legend=NA,inherit.aes=T,
        fun.geometry = NULL,size=size_text_plot)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y = element_text(size=13)
        #legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name="Mean changes per decade of \n10yrs return levels [%]",
                           limits=limits_mean[[v]],
                           oob=scales::squish)
    
    
    
    
    
    
    
    ############################################################################
    #                        20yrs quantile:
    ############################################################################
    flag_signif   = T
    flag.limits   = T
    scale_grad    = list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)
    crs_map       = crswgs84
    v             = 2
    name_var      = "Change in 20yrs quantiles per decade"
    type_estimate = c("MaxPost")
    name_model    = "SMEV Nonstationary"
    
    # initialize:
    a_w_MAP_df=c()
    rast_a_w_map_LAND_nonan=c()
    twentyyrchange=twentyyrchange_LAND=c()
    twentyyrchange_LAND_mean=twentyyrchange_LAND_mean_signif=c()
    
    # create dataframe for the map plot:
    a_w_MAP_df <- data.frame(x=lon_tmp,
                             y=lat_tmp, 
                             z=var_qnt_change_list$quant_smev_ns_20yr_change[[dur]])
    # trend_df$z[is.na(trend_df$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_NEW=a_w_MAP_df 
    } else {
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3), 
                                  z=round(a_w_MAP_df$z,3))
    }
    rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
    
    signif_ns_DIC_df_3 <- data.frame(x=lon_tmp,
                                     y=lat_tmp, 
                                     z=significance$DIC_diff3[[dur]])
    signif_ns_DIC_df_3$z[signif_ns_DIC_df_3$z==F]=NA
    signif_ns_DIC_df_3 = na.omit(signif_ns_DIC_df_3)
    signif_ns_DIC_df_3 = subset(signif_ns_DIC_df_3, z != 0)     
    
    twentyyrchange = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=a_w_MAP_df_NEW, aes(x=x, y=y, fill=z), alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
      geom_sf(data = gg_world, fill = NA, color = "grey70")+ #, linewidth=0.5) +
      geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+ #, linewidth=0.5)+
      coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
    if (flag_signif){
      twentyyrchange=twentyyrchange+
        #geom_point(data=signif_ns_DIC_df_1, aes(x=x, y=y, size=z),   fill="black", color="black", shape=21, size = 0.7)
        geom_point(data=signif_ns_DIC_df_3, aes(x=x, y=y, size=z),   fill="black", color="black", shape=21, size = 0.4)
    }
    twentyyrchange=twentyyrchange+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        #legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    # only land:
    ############
    # Replace nan with 1:
    a_w_MAP_df_LAND=a_w_MAP_df_NEW
    a_w_MAP_df_LAND$z[a_w_MAP_df_LAND$z==9999] =NA
    a_w_MAP_df_LAND_nonan <- a_w_MAP_df_LAND
    a_w_MAP_df_LAND_nonan$z[is.na(a_w_MAP_df_LAND_nonan$z)]=0
    # colMap <- colorRampPalette(c("yellow", "red"))(3)
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan
    } else {
      a_w_MAP_df_LAND_nonan_NEW=data.frame(x=round(a_w_MAP_df_LAND_nonan$x,3),
                                           y=round(a_w_MAP_df_LAND_nonan$y,3), 
                                           z=round(a_w_MAP_df_LAND_nonan$z,3))
    }
    rast_a_w_map_LAND_nonan =rasterFromXYZ(a_w_MAP_df_LAND_nonan_NEW, 
                                          crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2 <- crop(r, extent(land_shp))
    r3 <- mask(r2, land_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df <- as.data.frame(test_spdf)
    
    
    twentyyrchange_LAND = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=na.omit(test_df), aes(x=x, y=y, fill=z), alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z), alpha=1)+
      geom_sf(data = gg_world, fill = NA, color = "grey70")+ #, linewidth=0.5) +
      geom_sf(data = test, colour = "black", fill = NA, alpha = 0.2)+ #, linewidth=0.5)+
      coord_sf(xlim=c(6, 19), ylim=c(36, 48), expand = F)
    if (flag_signif){
      twentyyrchange_LAND=twentyyrchange_LAND+
        # geom_point(data=test_df_signif_DIC1, aes(x=x, y=y, size=z),   fill="black", color="black", shape=21, size = 0.7)
        geom_point(data=test_df_signif_DIC3, aes(x=x, y=y, size=z),  
                   fill="black", color="black", shape=21, size = 0.4)
    }
    twentyyrchange_LAND=twentyyrchange_LAND+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        #legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    
    # Aggregate to administrative regions:
    ######################################
    r=r2=r3=c()
    signif_reg_df=c()
    r=rast_a_w_map_LAND_nonan
    r2<-crop(r,extent(land_shp))
    r3<-mask(r2,land_shp)
    mavg<-terra::extract(r3,test,weights=T,fun=mean,na.rm=T)
    mavg_1=unlist(lapply(mavg, function(x) if (!is.null(x)) mean(x,na.rm=T) else NA))
    output=data.frame(mavg)
    # print(output)
    test$z=mavg
    
    mavg<-terra::extract(r3,land_shp,weights=T,fun=mean,na.rm=T)
    mavg_1=unlist(lapply(mavg, function(x) if (!is.null(x)) mean(x,na.rm=T) else NA))
    output=data.frame(mavg)
    # print(output)
    land_shp_NEW=land_shp
    land_shp_NEW$z= mavg
    
    twentyyrchange_LAND_mean=ggplot()+
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_sf(data=land_shp_NEW,colour="black",alpha=1)+
      geom_sf(data=test,aes(fill=z),colour="black",alpha=1)+
      geom_sf_text(
        data=test,aes(label=round(z,0)),
        stat="sf_coordinates", position = "nudge",
        parse=F,check_overlap=F,
        na.rm=F,show.legend=NA,inherit.aes=T,
        fun.geometry = NULL,size=size_text_plot)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y = element_text(size=13)
        #legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name="Mean changes per decade of \n20yrs return levels [%]",
                           limits=limits_mean[[v]],
                           oob=scales::squish)
    
    # Get % of signif pixels for each region:
    ##############################################
    if (sat_acron=='GSMAP'){
      signif_ns_DIC_df_3 <- data.frame(x=lon_tmp,
                                       y=lat_tmp, 
                                       z=significance$DIC_diff3[[dur]])  
    } else {
      signif_ns_DIC_df_3 <- data.frame(x=round(lon_tmp,3),
                                       y=round(lat_tmp,3), 
                                       z=significance$DIC_diff3[[dur]])
    }
    r=r2=r3=c()
    signif_reg_df=c()
    r=rasterFromXYZ(signif_ns_DIC_df_3,crs=crs_map) 
    r2<-crop(r,extent(test))
    r3<-mask(r2,test)
    signif_reg=terra::extract(r,test,weights=F,fun=sum,na.rm=T)
    signif_reg_df=data.frame(reg_name=test$reg_name,
                             signif_pixels_1=signif_reg)
    signif_reg_df$pix_reg_1=unlist(lapply(terra::extract(r,test,weights=F,na.rm=T,cells=T),length))
    # percent of pixels:
    signif_reg_df$percent_signif_1=round(signif_reg_df$signif_pixels_1/signif_reg_df$pix_reg_1*100,2)
    # number of signif pixels per region (accounting for coverage fraction)
    signif_reg_df$signif_pixels_2=ceil(exactextractr::exact_extract(r3,test,
                                        function(values,coverage_fraction) 
                                          sum(values*coverage_fraction,
                                              na.rm=T),progress=F))
    # Number of pixels per region:
    signif_reg_df$pix_reg_2=ceil(exactextractr::exact_extract(r3,test, 
                                function(values, coverage_fraction) 
                                  sum(coverage_fraction),progress=F))
    # Percent:
    signif_reg_df$percent_signif_2=round(signif_reg_df$signif_pixels_2/signif_reg_df$pix_reg_2*100,2)
    test$signif=signif_reg_df$percent_signif_2
    # # Number of distinct raster values within the polygon (coverage fractions are ignored)
    # signif_reg_df$dist_val = exact_extract(r, test, function(values, coverage_fraction) length(unique(values)))
    # signif_reg_df
    # Number of distinct raster values in cells more than 10% covered by the polygon
    # signif_reg = exact_extract(r, test, function(values, coverage_fraction) length(unique(values[coverage_fraction > 0.1])))
    
    twentyyrchange_LAND_mean_signif = ggplot()+
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_sf(data=land_shp_NEW, colour="black",alpha=1)+
      geom_sf(data=test,aes(fill=z),colour="black",alpha=1)+
      # geom_point(data=test_df_signif_DIC3,aes(x=x,y=y,size=z),  
      #            fill="black", color="black",shape=21,size=1)+
      geom_sf_text(
        data=test,aes(label=paste0(round(signif,0),'%')), 
        stat="sf_coordinates",position="nudge",
        parse=F,check_overlap=F,
        na.rm=F,show.legend=NA,inherit.aes=T,
        fun.geometry = NULL,size=size_text_plot)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y = element_text(size=13)
        #legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name="Mean changes per decade of \n20yrs return levels [%]",
                           limits=limits_mean[[v]],
                           oob=scales::squish)
    
    
    
    

    ############################################################################
    #                          50yrs quantile:
    ############################################################################
    flag_signif   = T
    flag.limits   = T
    scale_grad    = list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)
    crs_map       = crswgs84
    v             = 3
    name_var      = "Change in 50yrs quantiles per decade"
    type_estimate = c("MaxPost")
    name_model    = "SMEV Nonstationary"
    
    # initialize:
    a_w_MAP_df=c()
    rast_a_w_map_LAND_nonan=c()
    fiftyyrchange=fiftyyrchange_LAND=c()
    fiftyyrchange_LAND_mean=fiftyyrchange_LAND_mean_signif=c()
    
    # create dataframe for the map plot:
    a_w_MAP_df <- data.frame(x=lon_tmp, 
                             y=lat_tmp, 
                             z=var_qnt_change_list$quant_smev_ns_50yr_change[[dur]])
    # trend_df$z[is.na(trend_df$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_NEW=a_w_MAP_df 
    } else {
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3), 
                                  z=round(a_w_MAP_df$z,3))
    }
    rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
    
    signif_ns_DIC_df_3 <- data.frame(x=lon_tmp,
                                     y=lat_tmp, 
                                     z=significance$DIC_diff3[[dur]])
    signif_ns_DIC_df_3$z[signif_ns_DIC_df_3$z==F]=NA
    signif_ns_DIC_df_3 = na.omit(signif_ns_DIC_df_3)
    signif_ns_DIC_df_3 = subset(signif_ns_DIC_df_3, z != 0)     
    
    fiftyyrchange = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=a_w_MAP_df_NEW,aes(x=x,y=y,fill=z),alpha=1) + 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z),alpha=1)+
      geom_sf(data=gg_world,fill=NA,color="grey70")+ #,linewidth=0.5) +
      geom_sf(data=test,colour="black",fill=NA,alpha=0.2)+ #,linewidth=0.5)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)
    if (flag_signif){
      fiftyyrchange=fiftyyrchange+
        #geom_point(data=signif_ns_DIC_df_1, aes(x=x, y=y, size=z),   fill="black", color="black", shape=21, size = 0.7)
        geom_point(data=signif_ns_DIC_df_3,aes(x=x,y=y,size=z),
                   fill="black",color="black",shape=21,size=0.4)
    }
    fiftyyrchange=fiftyyrchange+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        #legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    ####################à
    # Replace nan with 1:
    a_w_MAP_df_LAND=a_w_MAP_df_NEW
    a_w_MAP_df_LAND$z[a_w_MAP_df_LAND$z==9999] =NA
    a_w_MAP_df_LAND_nonan <- a_w_MAP_df_LAND
    a_w_MAP_df_LAND_nonan$z[is.na(a_w_MAP_df_LAND_nonan$z)]=0
    # colMap <- colorRampPalette(c("yellow", "red"))(3)
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan
    } else {
      a_w_MAP_df_LAND_nonan_NEW=data.frame(x=round(a_w_MAP_df_LAND_nonan$x,3),
                                           y=round(a_w_MAP_df_LAND_nonan$y,3), 
                                           z=round(a_w_MAP_df_LAND_nonan$z,3))
    }
    rast_a_w_map_LAND_nonan=rasterFromXYZ(a_w_MAP_df_LAND_nonan_NEW, 
                                          crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2<-crop(r, extent(land_shp))
    r3<-mask(r2,land_shp)
    test_spdf<-as(r3, "SpatialPixelsDataFrame")
    test_df<-as.data.frame(test_spdf)
    
    fiftyyrchange_LAND = ggplot() + 
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_tile(data=na.omit(test_df),aes(x=x,y=y,fill=z),alpha=1)+ 
      #geom_raster(data=trend_df,aes(x=x,y=y,fill=z),alpha=1)+
      geom_sf(data=gg_world,fill=NA,color="grey70")+ #,linewidth=0.5) +
      geom_sf(data=test,colour="black",fill=NA,alpha=0.2)+ #,linewidth=0.5)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)
    if (flag_signif){
      fiftyyrchange_LAND=fiftyyrchange_LAND+
        # geom_point(data=test_df_signif_DIC1,aes(x=x,y=y,size=z),   
        #            fill="black",color="black",shape=21,size=0.7)
        geom_point(data=test_df_signif_DIC3,aes(x=x,y=y,size=z),  
                   fill="black",color="black",shape=21,size=0.4)
    }
    fiftyyrchange_LAND=fiftyyrchange_LAND+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        #legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name=paste0(name_var,'(',type_estimate,'), ',name_model),
                           limits=limits[[v]],
                           oob=scales::squish)
    
    # Aggregate to administrative regions:
    ######################################
    r=r2=r3=c()
    signif_reg_df=c()
    r=rast_a_w_map_LAND_nonan
    r2<-crop(r,extent(land_shp))
    r3<-mask(r2,land_shp)
    mavg<-terra::extract(r3,test,weights=T,fun=mean,na.rm=T)
    mavg_1=unlist(lapply(mavg, function(x) if (!is.null(x)) mean(x,na.rm=T) else NA))
    output=data.frame(mavg)
    # print(output)
    test$z=mavg
    
    mavg<-terra::extract(r3,land_shp,weights=T,fun=mean,na.rm=T)
    mavg_1=unlist(lapply(mavg, function(x) if (!is.null(x)) mean(x,na.rm=T) else NA))
    output=data.frame(mavg)
    # print(output)
    land_shp_NEW=land_shp
    land_shp_NEW$z= mavg
    
    fiftyyrchange_LAND_mean=ggplot()+
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_sf(data=land_shp_NEW,colour="black",alpha=1)+
      geom_sf(data=test,aes(fill=z),colour="black",alpha=1)+
      geom_sf_text(
        data=test,aes(label=round(z,0)),
        stat="sf_coordinates", position = "nudge",
        parse=F,check_overlap=F,
        na.rm=F,show.legend=NA,inherit.aes=T,
        fun.geometry = NULL,size=size_text_plot)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y = element_text(size=13)
        #legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name="Mean changes per decade of \n50yrs return levels [%]",
                           limits=limits_mean[[v]],
                           oob=scales::squish)
    
    # Get % of signif pixels for each region:
    ##############################################
    if (sat_acron=='GSMAP'){
      signif_ns_DIC_df_3 <- data.frame(x=lon_tmp,
                                       y=lat_tmp, 
                                       z=significance$DIC_diff3[[dur]])  
    } else {
      signif_ns_DIC_df_3 <- data.frame(x=round(lon_tmp,3),
                                       y=round(lat_tmp,3), 
                                       z=significance$DIC_diff3[[dur]])
    }
    r=r2=r3=c()
    signif_reg_df=c()
    r=rasterFromXYZ(signif_ns_DIC_df_3,crs=crs_map) 
    r2<-crop(r,extent(test))
    r3<-mask(r2,test)
    signif_reg=terra::extract(r,test,weights=F,fun=sum,na.rm=T)
    signif_reg_df=data.frame(reg_name=test$reg_name,
                             signif_pixels_1=signif_reg)
    signif_reg_df$pix_reg_1=unlist(lapply(terra::extract(r,test,weights=F,na.rm=T,cells=T),length))
    # percent of pixels:
    signif_reg_df$percent_signif_1=round(signif_reg_df$signif_pixels_1/signif_reg_df$pix_reg_1*100,2)
    # number of signif pixels per region (accounting for coverage fraction)
    signif_reg_df$signif_pixels_2=ceil(exactextractr::exact_extract(r3,test,
                                        function(values,coverage_fraction) 
                                          sum(values*coverage_fraction,
                                              na.rm=T),progress=F))
    # Number of pixels per region:
    signif_reg_df$pix_reg_2=ceil(exactextractr::exact_extract(r3,test, 
                                  function(values, coverage_fraction) 
                                    sum(coverage_fraction),progress=F))
    # Percent:
    signif_reg_df$percent_signif_2=round(signif_reg_df$signif_pixels_2/signif_reg_df$pix_reg_2*100,2)
    test$signif=signif_reg_df$percent_signif_2
    # # Number of distinct raster values within the polygon (coverage fractions are ignored)
    # signif_reg_df$dist_val = exact_extract(r, test, function(values, coverage_fraction) length(unique(values)))
    # signif_reg_df
    # Number of distinct raster values in cells more than 10% covered by the polygon
    # signif_reg = exact_extract(r, test, function(values, coverage_fraction) length(unique(values[coverage_fraction > 0.1])))
    
    fiftyyrchange_LAND_mean_signif = ggplot()+
      theme_bw()+
      # scale_x_continuous(expand=c(0,0),name="Longitude (°E)")+
      # scale_y_continuous(expand=c(0,0),name="Latitude (°N)")+
      scale_x_continuous(expand=c(0,0),name="Longitude (°E)", 
                         breaks=c(6.00001,10,14,18), 
                         labels=c("6°E","10°E","14°E","18°E"))+
      scale_y_continuous(expand=c(0,0),name="Latitude (°N)",  
                         breaks=c(36.00001,40,44,47.9999), 
                         labels=c("36°N","40°N","44°N","48°N"))+
      geom_sf(data=land_shp_NEW, colour="black",alpha=1)+
      geom_sf(data=test,aes(fill=z),colour="black",alpha=1)+
      # geom_point(data=test_df_signif_DIC3,aes(x=x,y=y,size=z),  
      #            fill="black", color="black",shape=21,size=1)+
      geom_sf_text(
        data=test,aes(label=paste0(round(signif,0),'%')), 
        stat="sf_coordinates",position="nudge",
        parse=F,check_overlap=F,
        na.rm=F,show.legend=NA,inherit.aes=T,
        fun.geometry = NULL,size=size_text_plot)+
      coord_sf(xlim=c(6,19),ylim=c(36,48),expand=F)+
      theme(
        panel.grid.major=element_blank()
        ,axis.text.x=element_text(size=13)
        ,axis.text.y = element_text(size=13)
        #legend.position="none"
        ,legend.text=element_text(size=15)
        ,legend.title=element_text(size=20))+
      scale_fill_gradient2(low=scale_grad$low,
                           high=scale_grad$high,
                           mid=scale_grad$mid,
                           midpoint=scale_grad$midpoint,
                           na.value=NA,
                           name="Mean changes per decade of \n50yrs return levels [%]",
                           limits=limits_mean[[v]],
                           oob=scales::squish)
    ############################################################################
    title0 <- ggdraw() +
      draw_label('',     fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25)+
      theme(plot.margin=margin(t=0,r=0,b=3,l=5))
    title1 <- ggdraw() +
      draw_label('2yrs', fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25)+
      theme(plot.margin=margin(t=0,r=0,b=3,l=5))
    title2 <- ggdraw() +
      draw_label('10yrs',fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25)+
      theme(plot.margin=margin(t=0,r=0,b=3,l=5))
    title3 <- ggdraw() +
      draw_label('20yrs',fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25)+
      theme(plot.margin=margin(t=0,r=0,b=3,l=5))
    title4 <- ggdraw() +
      draw_label('50yrs',fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25)+
      theme(plot.margin=margin(t=0,r=0,b=3,l=5))
    all.titles = plot_grid(title0,title1,title3,title4,
                           labels=NULL,
                           ncol=1,nrow=4,
                           rel_heights=c(0.25,1,1,1),
                           align='hv')
    # ALL PIXELS:
    all.plots=plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25),
                           twoyrchange+theme(legend.title=element_blank(), 
                                             legend.position="none", 
                                             plot.margin=margin(t=0,r=0,b=3,l=0),
                                             axis.title=element_blank()),
                           # tenyrchange+theme(legend.title=element_blank(), 
                           #                   legend.position="none", 
                           #                   plot.margin=margin(t=0,r=0,b=3,l=0),
                           #                   axis.title=element_blank()),
                           twentyyrchange+theme(legend.title=element_blank(), 
                                                legend.position="none", 
                                                plot.margin=margin(t=0,r=0,b=3,l=0),
                                                axis.title=element_blank()),
                           fiftyyrchange+theme(legend.title=element_blank(), 
                                               legend.position="none", 
                                               plot.margin=margin(t=0,r=0,b=3,l=0),
                                               axis.title=element_blank()),
                           labels=NULL,
                           ncol=1,nrow=5,
                           rel_heights=c(0.2, 1, 1, 1),
                           align = 'hv')
    
    
    all.all.qnt.plots=plot_grid(  # create title object:
        all.titles,
        all.plots, # + theme(plot.title = element_text(size= 20, hjust = 0.5)) + ggtitle("CMORPH\n"),
        label_size=25,
        labels=NULL,
        ncol=2,nrow=1,
        rel_widths=c(0.2,1),
        align='hv',
        hjust=-0.5,label_x=0)
    ggsave(all.all.qnt.plots,
           filename=paste0(dir.res_maps_qnt,'all_maps_qnt_',sat_prod,
                           '_',durations[dur]/60,'h.png'), 
           width=8,height=24,dpi=200,bg="white")
    
    
    # ONLY LAND:
    all.plots_LAND = plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25),
                               twoyrchange_LAND+theme(legend.title=element_blank(), 
                                                      legend.position="none", 
                                                      plot.margin=margin(t=0,r=0,b=3,l=0),
                                                      axis.title=element_blank()),
                               # tenyrchange_LAND+theme(legend.title=element_blank(), 
                               #                        legend.position="none", 
                               #                        plot.margin=margin(t=0,r=0,b=3,l=0),
                               #                        axis.title=element_blank()),
                               twentyyrchange_LAND+theme(legend.title=element_blank(), 
                                                         legend.position="none", 
                                                         plot.margin=margin(t=0,r=0,b=3,l=0),
                                                         axis.title=element_blank()),
                               fiftyyrchange_LAND+theme(legend.title=element_blank(), 
                                                        legend.position="none", 
                                                        plot.margin=margin(t=0,r=0,b=3,l=0),
                                                        axis.title=element_blank()),
                               labels=NULL,
                               ncol=1,nrow=4,
                               rel_heights=c(0.2,1,1,1),
                               align='hv')
    all.all.qnt.plots_LAND = plot_grid(
      all.titles,
      all.plots_LAND, # + theme(plot.title = element_text(size= 20, hjust = 0.5)) + ggtitle("CMORPH\n"),
      label_size=20,
      labels=NULL,
      ncol=2,nrow=1,
      rel_widths=c(0.2,1),
      align='hv',
      hjust=-0.5,label_x=0)
    
    ggsave(all.all.qnt.plots_LAND, 
           filename=paste0(dir.res_maps_qnt,'all_maps_qnt_LAND_',sat_prod,
                           '_',durations[dur]/60,'h.png'), 
           width=8,height=24,dpi=200,bg="white")
    
    
    # mean region:
    all.plots_LAND_mean = plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=25),
                                    twoyrchange_LAND_mean+theme(legend.title=element_blank(), 
                                                                legend.position="none", 
                                                                plot.margin=margin(t=0,r=0,b=3,l=0),
                                                                axis.title=element_blank()),
                                    # tenyrchange_LAND_mean+theme(legend.title=element_blank(), 
                                    #                        legend.position="none", 
                                    #                        plot.margin=margin(t=0,r=0,b=3,l=0),
                                    #                        axis.title=element_blank()),
                                    twentyyrchange_LAND_mean+theme(legend.title=element_blank(), 
                                                                   legend.position="none", 
                                                                   plot.margin=margin(t=0,r=0,b=3,l=0),
                                                                   axis.title=element_blank()),
                                    fiftyyrchange_LAND_mean+theme(legend.title=element_blank(), 
                                                                  legend.position="none", 
                                                                  plot.margin=margin(t=0,r=0,b=3,l=0),
                                                                  axis.title=element_blank()),
                                    labels=NULL,
                                    ncol=1,nrow=4,
                                    rel_heights = c(0.2,1,1,1),
                                    align='hv')
  
    all.all.qnt.plots_LAND_mean=plot_grid(
        all.titles,
        all.plots_LAND_mean, # + theme(plot.title = element_text(size= 20, hjust = 0.5)) + ggtitle("CMORPH\n"),
        label_size=20,
        labels=NULL,
        ncol=2,nrow=1,
        rel_widths=c(0.2,1),
        align='hv',
        hjust=-0.5,label_x=0)
    ggsave(all.all.qnt.plots_LAND_mean, 
           filename=paste0(dir.res_maps_qnt,'/all_maps_qnt_',sat_prod,
                           '_LAND_mean_',durations[dur]/60,'h.png'), 
           width=8,height=24,dpi=200,bg="white")
    
    
    # mean region + significant pixels:
    all.plots_LAND_mean_signif=plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=20),
                                         twoyrchange_LAND_mean_signif+theme(legend.title=element_blank(), 
                                                                            legend.position="none", 
                                                                            plot.margin=margin(t=0,r=0,b=3,l=0),
                                                                            axis.title=element_blank()),
                                         # tenyrchange_LAND_mean_signif +theme(legend.title   = element_blank(),
                                         #                                     legend.position= "none", 
                                         #                                     plot.margin    = margin(t=0, r=0, b=3, l=0),
                                         #                                     axis.title     = element_blank()),
                                         twentyyrchange_LAND_mean_signif+theme(legend.title=element_blank(), 
                                                                               legend.position="none", 
                                                                               plot.margin=margin(t=0,r=0,b=3,l=0),
                                                                               axis.title=element_blank()),
                                         fiftyyrchange_LAND_mean_signif+theme(legend.title=element_blank(), 
                                                                              legend.position="none", 
                                                                              plot.margin=margin(t=0,r=0,b=3,l=0),
                                                                              axis.title=element_blank()),
                                         labels=NULL,
                                         ncol=1,nrow=4,
                                         rel_heights = c(0.2,1,1,1),
                                         align='hv')
    
    all.all.qnt.plots_LAND_mean_signif = plot_grid(
      all.titles,
      all.plots_LAND_mean_signif,
      label_size=20,
      labels=NULL,
      ncol=2,nrow=1,
      rel_widths=c(0.2,1),
      align='hv',
      hjust=-0.5,label_x=0)
    ggsave(all.all.qnt.plots_LAND_mean_signif, 
           filename=paste0(dir.res_maps_qnt,'/all_maps_qnt_',sat_prod,
                           '_LAND_mean_signif_',durations[dur]/60,'h.png'), 
           width=8,height=24,dpi=200,bg="white")
  
    
    
    # save results and plots to .rds R file:
    saveRDS(list( twoyrchange_plot                  = twoyrchange,
                  twoyrchange_plot_LAND             = twoyrchange_LAND,
                  twoyrchange_plot_LAND_mean        = twoyrchange_LAND_mean,
                  twoyrchange_plot_LAND_mean_signif = twoyrchange_LAND_mean_signif,
                  # tenyrchange_plot                  = tenyrchange,
                  # tenyrchange_plot_LAND             = tenyrchange_LAND,
                  # tenyrchange_plot_LAND_mean        = tenyrchange_LAND_mean,
                  # tenyrchange_plot_LAND_mean_signif = tenyrchange_LAND_mean_signif,
                  twentyyrchange_plot               = twentyyrchange,
                  twentyyrchange_plot_LAND          = twentyyrchange_LAND,
                  twentyyrchange_plot_LAND_mean     = twentyyrchange_LAND_mean,
                  twentyyrchange_plot_LAND_mean_signif = twentyyrchange_LAND_mean_signif,
                  fiftyyrchange_plot                = fiftyyrchange,
                  fiftyyrchange_plot_LAND           = fiftyyrchange_LAND,
                  fiftyyrchange_plot_LAND_mean      = fiftyyrchange_LAND_mean,
                  fiftyyrchange_plot_LAND_mean_signif = fiftyyrchange_LAND_mean_signif,
                  lon_tmp                           = lon_tmp,
                  lat_tmp                           = lat_tmp,
                  durat                             = paste0(durations[dur]/60, 'h')), 
            file=paste0(dir.res_maps_qnt,'/all_maps_qnt_',
                        sat_prod,'_',durations[dur]/60,'h.rds'))
    
  } # end loop on durations
  ##############################################################################
}
























################################################################################
# plot of average quantiles per region with different durations
################################################################################
plot_maps_changes_durat=function(
                              var_durat_change_list = list(),   
                              lon_tmp        = NULL,
                              lat_tmp        = NULL,
                              durations      = c(60,180,360,1440), # in minutes
                              mask_shp       = NULL,
                              mask_world     = NULL,
                              mask_user      = NULL,
                              flag_signif    = F,
                              significance   = NULL,
                              sat_acron      = 'TEST',
                              period         = '',
                              year_ref       = NA,
                              year_ref_id    = NA,
                              path_output    = "/map_test",
                              limits         = c(-20, 20),
                              flag.limits    = T,
                              crs_map        = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
                              ){
################################################################################
  
  # initialize paths
  dir.res_maps_durat =path_output
  dir.create(dir.res_maps_durat)
  
  # initiliaze the lists with results for each duration:
  change_scale=change_shape=change_n=change_qnt2yrs=list()
  change_scale_stdev=change_shape_stdev=change_n_stdev=change_qnt2yrs_stdev=change_qnt2yrs_stdev_zones=list()
  change_qnt2yrs_zones=change_qnt2yrs_stdev_zones=list()
  change_shape_zones=change_shape_stdev_zones=list()
  change_scale_zones=change_scale_stdev_zones=list()
  change_n_zones=change_n_stdev_zones=list()  
  
  for (dur in 1:4){
    message(paste0('DURATION = ', durations[dur]/60, 'h'))
    
    
    ############################################################################
    #                        2yrs quantile:
    ############################################################################
    flag_signif   = T
    flag.limits   = T
    scale_grad    = list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)
    crs_map       = crswgs84
    v             = 1
    name_var      = "Change in 2yrs RL per decade"
    type_estimate = c("MaxPost")
    name_model    = "SMEV Nonstationary"
    
    # initialize:
    a_w_MAP_df=c()
    
    # create dataframe for the map plot:
    a_w_MAP_df <- data.frame(x=lon_tmp,
                             y=lat_tmp, 
                             z=quant_smev_ns_2yr_change[[dur]])
    # trend_df$z[is.na(trend_df$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_NEW=a_w_MAP_df 
    } else {
      a_w_MAP_df_NEW=data.frame(x=round(a_w_MAP_df$x,3),
                                y=round(a_w_MAP_df$y,3), 
                                z=round(a_w_MAP_df$z,3))
    }
    rast_a_w_MAP=rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
    
    # nonstationary model is more likely than stationary (DIC)
    if (sat_acron=='GSMAP'){
      signif_ns_DIC_df_3 <- data.frame(x=lon_tmp,
                                       y=lat_tmp, 
                                       z=significance$DIC_diff3[[dur]])
    } else {
      signif_ns_DIC_df_3 <- data.frame(x=round(lon_tmp,3),
                                       y=round(lat_tmp,3), 
                                       z=significance$DIC_diff3[[dur]])
    }
    signif_ns_DIC_df_3$z[signif_ns_DIC_df_3$z==F]=NA
    signif_ns_DIC_df_3 = na.omit(signif_ns_DIC_df_3)
    signif_ns_DIC_df_3 = subset(signif_ns_DIC_df_3, z != 0)    
    
    if (sat_acron=='GSMAP'){
      signif_ns_DIC_df_1 <- data.frame(x=lon_tmp,
                                       y=lat_tmp, 
                                       z=significance$DIC_diff1[[dur]])
    } else {
      signif_ns_DIC_df_1 <- data.frame(x=round(lon_tmp,3),
                                       y=round(lat_tmp,3),
                                       z=significance$DIC_diff1[[dur]])
    }
    signif_ns_DIC_df_1$z[signif_ns_DIC_df_1$z==F]=NA
    signif_ns_DIC_df_1 = na.omit(signif_ns_DIC_df_1)
    signif_ns_DIC_df_1 = subset(signif_ns_DIC_df_1, z != 0)  
    
    # only land:
    ############
    # Replace nan with 1:
    a_w_MAP_df_LAND=a_w_MAP_df_NEW
    a_w_MAP_df_LAND$z[a_w_MAP_df_LAND$z==9999] =NA
    a_w_MAP_df_LAND_nonan <- a_w_MAP_df_LAND
    a_w_MAP_df_LAND_nonan$z[is.na(a_w_MAP_df_LAND_nonan$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan
    } else {
      a_w_MAP_df_LAND_nonan_NEW =data.frame(x=round(a_w_MAP_df_LAND_nonan$x,3),
                                            y=round(a_w_MAP_df_LAND_nonan$y,3), 
                                            z=round(a_w_MAP_df_LAND_nonan$z,3))
    }
    rast_a_w_map_LAND_nonan=rasterFromXYZ(a_w_MAP_df_LAND_nonan_NEW, crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2 <- crop(r, extent(land_shp))
    r3 <- mask(r2, land_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df <- as.data.frame(test_spdf)
    
    # rast_signif_DIC3_LAND = rasterFromXYZ(signif_ns_DIC_df_3, crs=crs_map) 
    # r=rast_signif_DIC3_LAND
    # r2 <- crop(r, extent(land_shp))
    # r3 <- mask(r2, land_shp)
    # test_spdf <- as(r3, "SpatialPixelsDataFrame")
    # test_df_signif_DIC3 <- as.data.frame(test_spdf)
    # 
    # rast_signif_DIC1_LAND = rasterFromXYZ(signif_ns_DIC_df_1, crs=crs_map) 
    # r=rast_signif_DIC1_LAND
    # r2 <- crop(r, extent(land_shp))
    # r3 <- mask(r2, land_shp)
    # test_spdf <- as(r3, "SpatialPixelsDataFrame")
    # test_df_signif_DIC1 <- as.data.frame(test_spdf)
    
    # Aggregate to administrative regions:
    mavg=c()
    mavg<-terra::extract(r3,test,weights=T,fun=mean,na.rm=T)
    mavg_1=unlist(lapply(mavg, function(x) if (!is.null(x)) mean(x,na.rm=T) else NA))
    output=data.frame(mavg)
    change_qnt2yrs[[dur]]=mavg
    mavg_sd<-terra::extract(r3,test,weights=F,fun=sd,na.rm=T)
    change_qnt2yrs_stdev[[dur]]=mavg_sd
    
    # Aggregate to SPECIFIC POLYGONS:
    if (!is.null(mask_user)){
      mavg_zones=c()
      mavg_zones<-terra::extract(r3,user_polygons,weights=T,fun=mean,na.rm=T)
      mavg_1_zones=unlist(lapply(mavg_zones, function(x) if (!is.null(x)) mean(x,na.rm=T) else NA))
      output_zones=data.frame(mavg_zones)
      change_qnt2yrs_zones[[dur]]=mavg_zones
      mavg_sd_zones<-terra::extract(r3,user_polygons,weights=F,fun=sd,na.rm=T)
      change_qnt2yrs_stdev_zones[[dur]]=mavg_sd_zones
    }    
    
    
    
    
  
    
    ############################################################################
    #                                 scale:
    ############################################################################
    flag_signif   = T
    flag.limits   = T
    scale_grad    = list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)
    crs_map       = crswgs84
    v             = 2
    name_var      = "Change in scale per decade"
    type_estimate = c("MaxPost")
    name_model    = "SMEV Nonstationary"
    year_ref      = year_halfperiod
    year_ref_id   = halfperiod
    
    # initialize:
    a_w_MAP_df=c()
    scale_change = c()
    
    # initialize variable plot list:
    map_var_EU = list()
    map_var_LAND = list()
    
    # % of change of slope per decade wrw halfperiod year (a sort of average change per decade):
    scale_change = var_durat_change_list$b_C_MAP_ns[[dur]]*10 /
      (var_durat_change_list$a_C_MAP_ns[[dur]]+var_durat_change_list$b_C_MAP_ns[[dur]]*year_ref_id)*100
    
    # create dataframe for the map plot:
    a_w_MAP_df <- data.frame(x=lon_tmp, 
                             y=lat_tmp, 
                             z=scale_change)
    # trend_df$z[is.na(trend_df$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_NEW = a_w_MAP_df
    } else {
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3), 
                                  z=round(a_w_MAP_df$z,3))
    }
    rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
    
    # only land:
    ############
    # Replace nan with 1:
    a_w_MAP_df_LAND=a_w_MAP_df_NEW
    a_w_MAP_df_LAND$z[a_w_MAP_df_LAND$z==9999] =NA
    a_w_MAP_df_LAND_nonan <- a_w_MAP_df_LAND
    a_w_MAP_df_LAND_nonan$z[is.na(a_w_MAP_df_LAND_nonan$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan
    } else {
      a_w_MAP_df_LAND_nonan_NEW = data.frame(x=round(a_w_MAP_df_LAND_nonan$x,3),
                                             y=round(a_w_MAP_df_LAND_nonan$y,3), 
                                             z=round(a_w_MAP_df_LAND_nonan$z,3))
    }
    rast_a_w_map_LAND_nonan =rasterFromXYZ(a_w_MAP_df_LAND_nonan_NEW, 
                                           crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2 <- crop(r, extent(land_shp))
    r3 <- mask(r2, land_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df <- as.data.frame(test_spdf)
    
    # Aggregate to administrative regions:
    mavg<-terra::extract(r3,test,weights=T,fun=mean,na.rm=T)
    mavg_1=unlist(lapply(mavg,function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA))
    output=data.frame(mavg)
    change_scale[[dur]]=mavg
    mavg_sd<-terra::extract(r3,test,weights=F,fun=sd,na.rm=T)
    change_scale_stdev[[dur]]=mavg_sd
    
    # Aggregate to SPECIFIC POLYGONS:
    if (!is.null(mask_user)){
      mavg_zones=c()
      mavg_zones<-terra::extract(r3,user_polygons,weights=T,fun=mean,na.rm=T)
      mavg_1_zones=unlist(lapply(mavg_zones, function(x) if (!is.null(x)) mean(x,na.rm=T) else NA))
      output_zones=data.frame(mavg_zones)
      change_scale_zones[[dur]]=mavg_zones
      mavg_sd_zones<-terra::extract(r3,user_polygons,weights=F,fun=sd,na.rm=T)
      change_scale_stdev_zones[[dur]]=mavg_sd_zones
    }    
    
    
    
    
    
    
    
    ############################################################################
    #                             SHAPE:
    ############################################################################
    flag_signif   = T
    flag.limits   = T
    scale_grad    = list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)
    crs_map       = crswgs84
    v             = 3
    name_var      = "Change in shape per decade"
    type_estimate = c("MaxPost")
    name_model    = "SMEV Nonstationary"
    year_ref      = year_halfperiod
    year_ref_id   = halfperiod
    
    # initialize:
    a_w_MAP_df=c()
    shape_change = c()
    
    # initialize variable plot list:
    map_var_EU = list()
    map_var_LAND = list()
    
    # % of change of slope per decade wrw halfperiod year (a sort of average change per decade):
    shape_change = var_durat_change_list$b_w_MAP_ns[[dur]]*10 /
      (var_durat_change_list$a_w_MAP_ns[[dur]]+var_durat_change_list$b_w_MAP_ns[[dur]]*year_ref_id)*100
    
    # create dataframe for the map plot:
    a_w_MAP_df <- data.frame(x=lon_tmp, 
                             y=lat_tmp, 
                             z=shape_change)
    # trend_df$z[is.na(trend_df$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_NEW=a_w_MAP_df
    } else {
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3), 
                                  z=round(a_w_MAP_df$z,3))
    }
    rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
    
    # only land:
    ############
    # Replace nan with 1:
    a_w_MAP_df_LAND=a_w_MAP_df_NEW
    a_w_MAP_df_LAND$z[a_w_MAP_df_LAND$z==9999] =NA
    a_w_MAP_df_LAND_nonan <- a_w_MAP_df_LAND
    a_w_MAP_df_LAND_nonan$z[is.na(a_w_MAP_df_LAND_nonan$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan
    } else {
      a_w_MAP_df_LAND_nonan_NEW = data.frame(x=round(a_w_MAP_df_LAND_nonan$x,3),
                                             y=round(a_w_MAP_df_LAND_nonan$y,3), 
                                             z=round(a_w_MAP_df_LAND_nonan$z,3))
    }
    rast_a_w_map_LAND_nonan =rasterFromXYZ(a_w_MAP_df_LAND_nonan_NEW, 
                                           crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2 <- crop(r, extent(land_shp))
    r3 <- mask(r2, land_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df <- as.data.frame(test_spdf)
    
    # Aggregate to amministrative regions:
    mavg<-terra::extract(r3,test,weights=T,fun=mean,na.rm=T)
    mavg_1=unlist(lapply(mavg,function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA))
    output=data.frame(mavg)
    change_shape[[dur]]=mavg
    mavg_sd<-terra::extract(r3,test,weights=F,fun=sd,na.rm=T)
    change_shape_stdev[[dur]]=mavg_sd
    
    # Aggregate to SPECIFIC POLYGONS:
    if (!is.null(mask_user)){
      mavg_zones=c()
      mavg_zones<-terra::extract(r3,user_polygons,weights=T,fun=mean,na.rm=T)
      mavg_1_zones=unlist(lapply(mavg_zones, function(x) if (!is.null(x)) mean(x,na.rm=T) else NA))
      output_zones=data.frame(mavg_zones)
      change_shape_zones[[dur]]=mavg_zones
      mavg_sd_zones<-terra::extract(r3,user_polygons,weights=F,fun=sd,na.rm=T)
      change_shape_stdev_zones[[dur]]=mavg_sd_zones
    }   
    
    
    
    
    
    
    
    
    
    ############################################################################
    #                               n
    ############################################################################
    flag_signif   = T
    flag.limits   = T
    scale_grad    = list(name="rdbl",low="red",high="blue",mid="white",midpoint=0)
    crs_map       = crswgs84
    v             = 3
    name_var      = "Change in n per decade"
    type_estimate = c("MaxPost")
    name_model    = "SMEV Nonstationary"
    year_ref      = year_halfperiod
    year_ref_id   = halfperiod
    
    # initialize:
    a_w_MAP_df=c()
    change_n=c()
    # initialize variable plot list:
    map_var_EU = list()
    map_var_LAND = list()
    
    # create dataframe for the map plot:
    a_w_MAP_df <- data.frame(x=lon_tmp, 
                             y=lat_tmp, 
                             z=n_change[[dur]])
    # trend_df$z[is.na(trend_df$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_NEW=a_w_MAP_df
    } else {
      a_w_MAP_df_NEW = data.frame(x=round(a_w_MAP_df$x,3),
                                  y=round(a_w_MAP_df$y,3), 
                                  z=round(a_w_MAP_df$z,3))
    }
    rast_a_w_MAP =rasterFromXYZ(a_w_MAP_df_NEW, crs=crs_map)
    
    # only land:
    ############
    # Replace nan with 1:
    a_w_MAP_df_LAND=a_w_MAP_df_NEW
    a_w_MAP_df_LAND$z[a_w_MAP_df_LAND$z==9999] =NA
    a_w_MAP_df_LAND_nonan <- a_w_MAP_df_LAND
    a_w_MAP_df_LAND_nonan$z[is.na(a_w_MAP_df_LAND_nonan$z)]=0
    if (sat_acron=='GSMAP'){
      a_w_MAP_df_LAND_nonan_NEW=a_w_MAP_df_LAND_nonan
    } else {
      a_w_MAP_df_LAND_nonan_NEW = data.frame(x=round(a_w_MAP_df_LAND_nonan$x,3),
                                             y=round(a_w_MAP_df_LAND_nonan$y,3), 
                                             z=round(a_w_MAP_df_LAND_nonan$z,3))
    }
    rast_a_w_map_LAND_nonan =rasterFromXYZ(a_w_MAP_df_LAND_nonan_NEW, 
                                           crs=crs_map) 
    ## crop and mask (only Italy extent):
    r=r2=r3=c()
    r=rast_a_w_map_LAND_nonan
    r2 <- crop(r, extent(land_shp))
    r3 <- mask(r2, land_shp)
    test_spdf <- as(r3, "SpatialPixelsDataFrame")
    test_df <- as.data.frame(test_spdf)
    
    # Aggregate to amministrative regions:
    mavg<-terra::extract(r3,test,weights=T,fun=mean,na.rm=T)
    mavg_1=unlist(lapply(mavg,function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA))
    output=data.frame(mavg)
    change_n[[dur]]=mavg
    mavg_sd<-terra::extract(r3,test,weights=F,fun=sd,na.rm=T)
    change_n_stdev[[dur]]=mavg_sd
    
    # Aggregate to SPECIFIC POLYGONS:
    if (!is.null(mask_user)){
      mavg_zones=c()
      mavg_zones<-terra::extract(r3,user_polygons,weights=T,fun=mean,na.rm=T)
      mavg_1_zones=unlist(lapply(mavg_zones, function(x) if (!is.null(x)) mean(x,na.rm=T) else NA))
      output_zones=data.frame(mavg_zones)
      change_n_zones[[dur]]=mavg_zones
      mavg_sd_zones<-terra::extract(r3,user_polygons,weights=F,fun=sd,na.rm=T)
      change_n_stdev_zones[[dur]]=mavg_sd_zones
    }   
    ##########################  
  } # end loop on durations
  
  
  
  
  # save results and plots to .rds R file:
  message('Saving results to .rds file...')
  saveRDS(list( change_shape_stdev   = change_shape_stdev,
                change_shape         = change_shape,
                change_shape_stdev_zones=change_shape_stdev_zones,
                change_shape_zones   = change_shape_zones,
                change_scale_stdev   = change_scale_stdev,
                change_scale         = change_scale,
                change_scale_stdev_zones=change_scale_stdev_zones,
                change_scale_zones   = change_scale_zones,
                change_qnt2yrs_stdev = change_qnt2yrs_stdev,
                change_qnt2yrs       = change_qnt2yrs,
                change_qnt2yrs_stdev_zones=change_qnt2yrs_stdev_zones,
                change_qnt2yrs_zones = change_qnt2yrs_zones,
                change_n             = change_n,
                change_n_stdev       = change_n_stdev,
                change_n_zones       = change_n_zones,
                change_n_stdev_zones = change_n_stdev_zones,
                lon_tmp              = lon_tmp,
                lat_tmp              = lat_tmp,
                durat                = paste0(durations[dur]/60, 'h')), 
          file=paste0(dir.res_maps_durat,'/all_maps_unc_durat_',sat_prod,'.rds'))
  
  
  
  
  
  
  # one plot for each ZONE:
  ##############################
  changes_with_duration_LAND_mean=list()
  df_changes_durat_zones =list()
  df_changes_durat_zones_merged=c()
  message('Plotting merged results aggregating at zone level...')
  for (z in 1:nrow(mask_user)){
  ############################
    change_shape_unlist = unlist(lapply(change_shape_zones, function(l) l[[z]]))  
    change_scale_unlist = unlist(lapply(change_scale_zones, function(l) l[[z]]))  
    change_n_unlist = unlist(lapply(change_n_zones, function(l) l[[z]]))  
    change_qnt2yrs_unlist = unlist(lapply(change_qnt2yrs_zones, function(l) l[[z]]))
    
    change_shape_stdev_unlist = unlist(lapply(change_shape_stdev_zones, function(l) l[[z]]))  
    change_scale_stdev_unlist = unlist(lapply(change_scale_stdev_zones, function(l) l[[z]]))  
    change_n_stdev_unlist = unlist(lapply(change_n_stdev_zones, function(l) l[[z]]))
    change_qnt2yrs_stdev_unlist = unlist(lapply(change_qnt2yrs_stdev_zones, function(l) l[[z]]))
    
    df_changes_durat_zones[[z]] =data.frame(dur         =c(1, 2, 3, 4), 
                                            shape_mean  =change_shape_unlist,
                                            shape_sd    =change_shape_stdev_unlist,
                                            scale_mean  =change_scale_unlist,
                                            scale_sd    =change_scale_stdev_unlist,
                                            n_mean      =change_n_unlist,
                                            n_sd        =change_n_stdev_unlist,
                                            qnt2yrs_mean=change_qnt2yrs_unlist,
                                            qnt2yrs_sd  =change_qnt2yrs_stdev_unlist)
    df_changes_durat_zones_merged=rbind(df_changes_durat_zones_merged,
                                   data.frame(dur=c(1, 2, 3, 4), 
                                              var_name= c(rep("Shape",4),rep("Scale",4),rep("n",4), rep("qnt2yrs",4)),
                                              stdev = c(change_shape_stdev_unlist,
                                                        change_scale_stdev_unlist,
                                                        change_n_stdev_unlist,
                                                        change_qnt2yrs_stdev_unlist),
                                              mean = c(change_shape_unlist,
                                                       change_scale_unlist,
                                                       change_n_unlist,
                                                       change_qnt2yrs_unlist),
                                              reg_i    =z,
                                              id       =mask_user$id[z],
                                              reg_name =mask_user$name[z]))
  }
  
  
  df_changes_durat_zones_merged_sort=df_changes_durat_zones_merged %>%arrange(id)

  ##############################################################################
  shape_col="shape"; scale_col="scale"; n_col="n"; qnt2yrs_col="2yrs RL"
  color.curves = c("shape"  ="black",
                   "scale"  ="red",
                   "n"      ="black",
                   "2yrs rt"="blue")
  
  plot_facets=ggplot(data=df_changes_durat_zones_merged_sort, aes(x=dur,y=mean)) +
    # geom_jitter() +
    theme_bw()+
    scale_x_continuous(expand=c(0,0),name="Durations [h]",
                       labels=c("1","3","6","24"),
                       breaks=c(1,2,3,4),
                       limits=c(0.8,4.2))+
    scale_y_continuous(expand=c(0,0),name="Changes [%/decade]",limits=limits,oob=scales::squish)+
    #facet_wrap(~reg_name, nrow = 1, ncol=10)+
    facet_wrap(~factor(reg_name,levels=unique(df_changes_durat_zones_merged_sort$reg_name)),
               nrow = 1, ncol=10)+ #,scale="free_y")+
    geom_hline(aes(yintercept=0), color="gray50", linetype="dashed", linewidth=0.5)+
    geom_line(aes(x=dur, y=mean, color=var_name, linetype=var_name), linewidth=3)+
    #geom_ribbon(aes(x=dur, ymin=mean-2*stdev, ymax=mean+2*stdev, fill=var_name), alpha=0.2)+
    #geom_ribbon(aes(x=dur, ymin=mean-stdev, ymax=mean+stdev, fill=var_name), alpha=0.2)+
    geom_point(aes(x=dur, y=mean, fill=var_name), shape=21, size=7)+
    coord_cartesian(clip = 'off')+
    # scale_fill_manual(name= "Percent changes in:\n",
    #                   labels =c(shape_col, scale_col, qnt2yrs_col),
    #                   breaks =c(shape_col, scale_col, qnt2yrs_col),  
    #                   values =color.curves,
    #                   guide    = guide_legend(override.aes = list(
    #                     linetype = c("solid", "solid", "solid"),
    #                     shape    = c(21, 21, 21))))+
    theme(
      panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
      ,axis.title       = element_text(size=35)
      ,axis.text.x      = element_text(size=30)
      ,axis.text.y      = element_text(size=30)
      ,legend.text      = element_text(size=40)
      ,legend.title     = element_blank() #element_text(size=20)
      ,legend.position  = "right"
      ,legend.direction = "vertical"
      ,strip.text       = element_text(size=40)) #,angle=90))
  ggsave(plot_facets, 
         filename=paste0(dir.res_maps_durat, '/all_maps_changes_durat_',sat_prod,'_ZONES.png'),
         width = 46, height =7, dpi = 150, bg = "white")
  
  plot_facets_all=plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',
                                                x=0.5,hjust=0.5,vjust=0.5,size=30),
                            plot_facets,
                            labels=NULL, # c("", "2yrs", "scale", "shape", "n"),
                            ncol=2,nrow=1,
                            rel_widths=c(0.05,1),
                            align='hv')
  ggsave(plot_facets_all, 
         filename=paste0(dir.res_maps_durat, '/all_maps_changes_durat_',sat_prod,'_ZONES_bis.png'),
         width = 46, height =7, dpi = 150, bg = "white")
  
  
  plot_facets_all_zones=plot_grid(plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=40),
                                          plot_facets+theme(
                                            #legend.position="none"
                                            axis.title.y=element_blank()
                                            ,axis.text.x=element_blank()
                                            ,axis.title.x=element_blank()
                                            ,strip.text.x=element_text(size=40)),#,angle=90)), 
                                          #,strip.background=element_blank()),
                                          labels=NULL,
                                          ncol=2,nrow=1,
                                          rel_widths=c(0.07,1),
                                          align='hv'),
                                plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=40),
                                          plot_facets+theme(
                                            #legend.position="none"
                                            axis.title.y=element_blank()
                                            ,axis.text.x=element_blank()
                                            ,axis.title.x=element_blank(),
                                            ,strip.background=element_blank()
                                            ,strip.text.x=element_blank()),
                                          labels=NULL, # c("", "2yrs", "scale", "shape", "n"),
                                          ncol=2,nrow=1,
                                          rel_widths=c(0.07,1),
                                          align='hv'),
                                plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=40),
                                          plot_facets+theme(
                                            #legend.position="none"
                                            axis.title.y=element_blank()
                                            ,axis.text.x=element_blank()
                                            ,axis.title.x=element_blank(),
                                            ,strip.background=element_blank()
                                            ,strip.text.x = element_blank()),
                                          labels=NULL, # c("", "2yrs", "scale", "shape", "n"),
                                          ncol=2,nrow=1,
                                          rel_widths=c(0.07,1),
                                          align='hv'),
                                plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=40),
                                          plot_facets+theme(
                                            #legend.position="none"
                                            axis.title.y=element_blank()
                                            # ,axis.text.x=element_blank()
                                            # ,axis.title.x=element_blank(),
                                            ,strip.background=element_blank()
                                            ,strip.text.x = element_blank()),
                                          labels=NULL, # c("", "2yrs", "scale", "shape", "n"),
                                          ncol=2,nrow=1,
                                          rel_widths=c(0.07,1),
                                          align='hv'),
                                labels=NULL,
                                ncol=1,nrow=4,
                                rel_heights = c(1.1,1,1,1.1),
                                align='hv')
  ggsave(plot_facets_all_zones, 
         filename=paste0(dir.res_maps_durat, '/all_maps_changes_durat_',sat_prod,'_bis_ZONES.png'),
         width = 47, height =30, dpi = 150, bg = "white")
  
  
  
  # # vertical alignment:
  # #####################
  # plot_facets_vert=ggplot(data=df_changes_durat_zones_merged, aes(x=dur,y=mean)) +
  #   # geom_jitter() +
  #   theme_bw()+
  #   scale_x_continuous(expand=c(0,0),name="Durations [h]",
  #                      labels=c("1","3","6","24"),
  #                      breaks=c(1,2,3,4),
  #                      limits=c(0.8,4.2))+
  #   scale_y_continuous(expand=c(0,0),name="Changes [%/decade]",limits=limits,oob=scales::squish)+
  #   facet_wrap(~reg_name, ncol = 1, nrow=10)+
  #   geom_hline(aes(yintercept=0), color="gray50", linetype="dashed")+
  #   geom_line(aes(x=dur, y=mean, color=var_name, linetype=var_name), linewidth=2)+
  #   #geom_ribbon(aes(x=dur, ymin=mean-stdev, ymax=mean+stdev, fill=var_name), alpha=0.2)+
  #   #geom_ribbon(aes(x=dur, ymin=mean-2*stdev, ymax=mean+2*stdev, fill=var_name), alpha=0.2)+
  #   geom_point(aes(x=dur, y=mean, fill=var_name), shape=21, size=5)+
  #   coord_cartesian(clip = 'off')+
  #   # scale_fill_manual(name= "Percent changes in:\n",
  #   #                   labels =c(shape_col, scale_col, qnt2yrs_col),
  #   #                   breaks =c(shape_col, scale_col, qnt2yrs_col),  
  #   #                   values =color.curves,
  #   #                   guide    = guide_legend(override.aes = list(
  #   #                     linetype = c("solid", "solid", "solid"),
  #   #                     shape    = c(21, 21, 21))))+
  #   theme(
  #     panel.grid.major = element_blank()
  #     ,panel.grid.minor = element_blank()
  #     ,axis.title       = element_text(size=30)
  #     ,axis.text.x      = element_text(size=25)
  #     ,axis.text.y      = element_text(size=25)
  #     ,legend.text      = element_text(size=40)
  #     ,legend.title     = element_blank() #element_text(size=20)
  #     ,legend.position  = "right"
  #     ,legend.direction = "vertical"
  #     ,strip.text       = element_text(size=40))# ,angle=90))
  # ggsave(plot_facets_vert, 
  #        filename=paste0(dir.res_maps_durat, '/all_maps_changes_durat_',sat_prod,'_ZONES_vert.png'),
  #        width = 10, height =40, dpi = 200, bg = "white")
  # 
  # plot_facets_all_vert=plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',
  #                                               x=0.5,hjust=0.5,vjust=0.5,size=30),
  #                           plot_facets_vert,
  #                           labels=NULL, # c("", "2yrs", "scale", "shape", "n"),
  #                           ncol=2,nrow=1,
  #                           rel_widths=c(0.04,1),
  #                           align='hv')
  # ggsave(plot_facets_all_vert, 
  #        filename=paste0(dir.res_maps_durat, '/all_maps_changes_durat_',sat_prod,'_ZONES_vert_bis.png'),
  #        width = 10, height =40, dpi = 200, bg = "white")
  # 
  # 
  # 
  # plot_facets_all_zones_vert=plot_grid(plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=40),
  #                                                plot_facets_vert+theme(
  #                                             #legend.position="none"
  #                                             axis.title.y=element_blank()
  #                                             ,axis.text.x=element_blank()
  #                                             ,axis.title.x=element_blank()
  #                                             ,strip.text.x=element_text(size=40,angle=90)), 
  #                                           #,strip.background=element_blank()),
  #                                           labels=NULL,
  #                                           ncol=2,nrow=1,
  #                                           rel_widths=c(0.07,1),
  #                                           align='hv'),
  #                                 plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=40),
  #                                           plot_facets_vert+theme(
  #                                             #legend.position="none"
  #                                             axis.title.y=element_blank()
  #                                             ,axis.text.x=element_blank()
  #                                             ,axis.title.x=element_blank(),
  #                                             ,strip.background=element_blank()
  #                                             ,strip.text.x=element_blank()),
  #                                           labels=NULL, # c("", "2yrs", "scale", "shape", "n"),
  #                                           ncol=2,nrow=1,
  #                                           rel_widths=c(0.07,1),
  #                                           align='hv'),
  #                                 plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=40),
  #                                           plot_facets_vert+theme(
  #                                             #legend.position="none"
  #                                             axis.title.y=element_blank()
  #                                             ,axis.text.x=element_blank()
  #                                             ,axis.title.x=element_blank(),
  #                                             ,strip.background=element_blank()
  #                                             ,strip.text.x = element_blank()),
  #                                           labels=NULL, # c("", "2yrs", "scale", "shape", "n"),
  #                                           ncol=2,nrow=1,
  #                                           rel_widths=c(0.07,1),
  #                                           align='hv'),
  #                                 plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=40),
  #                                           plot_facets_vert+theme(
  #                                             #legend.position="none"
  #                                             axis.title.y=element_blank()
  #                                             # ,axis.text.x=element_blank()
  #                                             # ,axis.title.x=element_blank(),
  #                                             ,strip.background=element_blank()
  #                                             ,strip.text.x = element_blank()),
  #                                           labels=NULL, # c("", "2yrs", "scale", "shape", "n"),
  #                                           ncol=2,nrow=1,
  #                                           rel_widths=c(0.07,1),
  #                                           align='hv'),
  #                                 labels=NULL,
  #                                 ncol=4,nrow=1,
  #                                 rel_heights = c(1,1,1,1),
  #                                 align='hv')
  # ggsave(plot_facets_all_zones_vert, 
  #        filename=paste0(dir.res_maps_durat, '/all_maps_changes_durat_',sat_prod,'_vert_bis_ZONES.png'),
  #        width = 45, height =30, dpi = 150, bg = "white")
  # 
  # 
  

  
  
  
  
  
  
  # one plot for each region:
  ###########################
  changes_with_duration_LAND_mean=list()
  df_changes_durat =list()
  df_changes_durat_merged=c()
  numb_regions = 20
  # some modifications to df names:
  test$reg_name[9]="E.R."
  test$reg_name[18]="F.V.G."
  test$reg_name[19]="Trentino A.A."
  message('Plotting merged results aggregating at administrative level (Regions)...')
  
  
  for (reg in 1:numb_regions){
  ############################
    change_shape_unlist = unlist(lapply(change_shape, function(l) l[[reg]]))  
    change_scale_unlist = unlist(lapply(change_scale, function(l) l[[reg]]))  
    change_qnt2yrs_unlist = unlist(lapply(change_qnt2yrs, function(l) l[[reg]]))
    change_shape_stdev_unlist = unlist(lapply(change_shape_stdev, function(l) l[[reg]]))  
    change_scale_stdev_unlist = unlist(lapply(change_scale_stdev, function(l) l[[reg]]))  
    change_qnt2yrs_stdev_unlist = unlist(lapply(change_qnt2yrs_stdev, function(l) l[[reg]]))
    
    df_changes_durat[[reg]] =data.frame(dur=c(1, 2, 3, 4), 
                                        shape_mean  =change_shape_unlist,
                                        shape_sd    =change_shape_stdev_unlist,
                                        scale_mean  =change_scale_unlist,
                                        scale_sd    =change_scale_stdev_unlist,
                                        qnt2yrs_mean=change_qnt2yrs_unlist,
                                        qnt2yrs_sd  =change_qnt2yrs_stdev_unlist)
    
    # df_changes_durat_merged =rbind(df_changes_durat_merged,
    #                                data.frame(dur=c(1, 2, 3, 4), 
    #                                     shape_mean  =change_shape_unlist,
    #                                     shape_sd    =change_shape_stdev_unlist,
    #                                     scale_mean  =change_scale_unlist,
    #                                     scale_sd    =change_scale_stdev_unlist,
    #                                     qnt2yrs_mean=change_qnt2yrs_unlist,
    #                                     qnt2yrs_sd  =change_qnt2yrs_stdev_unlist,
    #                                     reg_id      =reg,
    #                                     reg         =test$reg_name[reg]))
    df_changes_durat_merged =rbind(df_changes_durat_merged,
                                   data.frame(dur=c(1, 2, 3, 4), 
                                              var_name= c(rep("Shape",4),rep("Scale",4), rep("qnt2yrs",4)),
                                              stdev = c(change_shape_stdev_unlist,
                                                        change_scale_stdev_unlist,
                                                        change_qnt2yrs_stdev_unlist),
                                              mean = c(change_shape_unlist,
                                                       change_scale_unlist,
                                                       change_qnt2yrs_unlist),
                                              reg_i    =reg,
                                              reg_name =test$reg_name[reg]))
    
    shape_col="shape"; scale_col="scale"; qnt2yrs_col="2yrs RL"
    color.curves = c("shape"  ="black",
                     "scale"  ="red",
                     "2yrs rt"="blue")
    # changes_with_duration_LAND_mean[[reg]] = ggplot(data=df_changes_durat[[reg]])+
    #   theme_bw()+
    #   scale_x_continuous(expand = c(0, 0), name = "Durations", labels=c("1h", "3h", "6h", "24h"))+
    #   scale_y_continuous(expand=c(0,0),name="Changes [%/decade]",limits=limits,oob=scales::squish)+
    #   geom_hline(aes(yintercept=0), color="gray20", linetype="dashed")+
    #   geom_line(aes(x=dur, y=shape_mean), color="black", linewidth=1.1)+
    #   geom_ribbon(aes(x=dur, ymin=shape_mean-2*shape_sd, ymax=shape_mean+2*shape_sd), fill="black", alpha=0.1)+
    #   geom_line(aes(x=dur, y=scale_mean), color="red", linewidth=1.1)+
    #   geom_ribbon(aes(x=dur, ymin=scale_mean-2*scale_sd, ymax=scale_mean+2*scale_sd), fill="red", alpha=0.1)+
    #   geom_line(aes(x=dur, y=qnt2yrs_mean), color="blue", linewidth=1.1)+
    #   geom_ribbon(aes(x=dur, ymin=qnt2yrs_mean-2*qnt2yrs_sd, ymax=qnt2yrs_mean+2*qnt2yrs_sd), fill="blue", alpha=0.1)+
    #   geom_point(aes(x=dur, y=shape_mean, fill=shape_col), shape=21, size=3)+
    #   geom_point(aes(x=dur, y=scale_mean, fill=scale_col), shape=21, size=3)+
    #   geom_point(aes(x=dur, y=qnt2yrs_mean, fill=qnt2yrs_col), shape=21,size=3)+
    #   coord_cartesian(clip = 'off')+
    #   scale_fill_manual(name= "Percent changes in:\n",
    #                     labels =c(shape_col, scale_col, qnt2yrs_col),
    #                     breaks =c(shape_col, scale_col, qnt2yrs_col),  
    #                     values =color.curves,
    #                     guide    = guide_legend(override.aes = list(
    #                       linetype = c("solid", "solid", "solid"),
    #                       shape    = c(21, 21, 21))))+
    #   theme(
    #     panel.grid.major = element_blank()
    #     ,panel.grid.minor = element_blank()
    #     ,axis.title       = element_text(size=17)
    #     ,axis.text.x      = element_text(size=17)
    #     ,axis.text.y      = element_text(size=17)
    #     ,legend.text      = element_text(size=20)
    #     ,legend.title     = element_text(size=20)
    #     ,legend.position  = "none" # "right"
    #     ,legend.direction = "vertical")
  }
  # changes_with_duration_grid_plot= plot_grid(plotlist=changes_with_duration_LAND_mean,
  #                                            label_size=15,
  #                                            labels=as.list(test$reg_name),
  #                                            ncol =4,nrow=5,
  #                                            rel_heights=1,
  #                                            align ='hv')
  # ggsave(changes_with_duration_grid_plot,
  #        filename=paste0(dir.res_maps_durat, '/all_maps_changes_durat_',sat_prod,'.png'),
  #        width = 10, height =24, dpi = 200, bg = "white")
  
  

  
  
  
  
  ################################################################################
  # plot only for one product and for all variables (scale, shape, qnt2yrs) + uncert:
  plot_facets_reg =ggplot(data=df_changes_durat_merged, aes(x=dur,y=mean)) +
    # geom_jitter() +
    theme_bw()+
    scale_x_continuous(expand=c(0,0),name="Durations [h]",
                       labels=c("1","3","6","24"),
                       breaks=c(1,2,3,4),
                       limits=c(0.8,4.2))+
    scale_y_continuous(expand=c(0,0),name="Changes [%/decade]",limits=limits,oob=scales::squish)+
    facet_wrap(~reg_name, nrow = 1, ncol=20)+
    geom_hline(aes(yintercept=0), color="gray20", linetype="dashed")+
    geom_line(aes(x=dur, y=mean, color=var_name, linetype=var_name), linewidth=2)+
    geom_ribbon(aes(x=dur, ymin=mean-2*stdev, ymax=mean+2*stdev, fill=var_name), alpha=0.2)+
    geom_point(aes(x=dur, y=mean, fill=var_name), shape=21, size=5)+
    coord_cartesian(clip = 'off')+
    # scale_fill_manual(name= "Percent changes in:\n",
    #                   labels =c(shape_col, scale_col, qnt2yrs_col),
    #                   breaks =c(shape_col, scale_col, qnt2yrs_col),  
    #                   values =color.curves,
    #                   guide    = guide_legend(override.aes = list(
    #                     linetype = c("solid", "solid", "solid"),
    #                     shape    = c(21, 21, 21))))+
    theme(
      panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
      ,axis.title       = element_text(size=25)
      ,axis.text.x      = element_text(size=25)
      ,axis.text.y      = element_text(size=25)
      ,legend.text      = element_text(size=25)
      ,legend.title     = element_blank() #element_text(size=20)
      ,legend.position  = "right"
      ,legend.direction = "vertical"
      ,strip.text       = element_text(size=30,angle=90))
  ggsave(plot_facets_reg, 
         filename=paste0(dir.res_maps_durat, '/all_maps_changes_durat_',sat_prod,'_reg.png'),
         width = 42, height =7, dpi = 200, bg = "white")
  
  plot_facets_reg_all=plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',
                                                x=0.5,hjust=0.5,vjust=0.5,size=30),
                            plot_facets_reg,
                            labels=NULL, # c("", "2yrs", "scale", "shape", "n"),
                            ncol=2,nrow=1,
                            rel_widths=c(0.04,1),
                            align='hv')
  ggsave(plot_facets_all, 
         filename=paste0(dir.res_maps_durat, '/all_maps_changes_durat_',sat_prod,'_reg_ALL.png'),
         width = 42, height =7, dpi = 150, bg = "white")
  

  
  plot_facets_reg_all_all=plot_grid(plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=40),
                                              plot_facets_reg+theme(
                                            #legend.position="none"
                                            axis.title.y=element_blank()
                                            ,axis.text.x=element_blank()
                                            ,axis.title.x=element_blank()
                                            ,strip.text.x=element_text(size=40,angle=90)), 
                                          #,strip.background=element_blank()),
                                          labels=NULL,
                                          ncol=2,nrow=1,
                                          rel_widths=c(0.07,1),
                                          align='hv'),
                                plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=40),
                                          plot_facets_reg+theme(
                                            #legend.position="none"
                                            axis.title.y=element_blank()
                                            ,axis.text.x=element_blank()
                                            ,axis.title.x=element_blank(),
                                            ,strip.background=element_blank()
                                            ,strip.text.x=element_blank()),
                                          labels=NULL, # c("", "2yrs", "scale", "shape", "n"),
                                          ncol=2,nrow=1,
                                          rel_widths=c(0.07,1),
                                          align='hv'),
                                plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=40),
                                          plot_facets_reg+theme(
                                            #legend.position="none"
                                            axis.title.y=element_blank()
                                            ,axis.text.x=element_blank()
                                            ,axis.title.x=element_blank(),
                                            ,strip.background=element_blank()
                                            ,strip.text.x = element_blank()),
                                          labels=NULL, # c("", "2yrs", "scale", "shape", "n"),
                                          ncol=2,nrow=1,
                                          rel_widths=c(0.07,1),
                                          align='hv'),
                                plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=40),
                                          plot_facets_reg+theme(
                                            #legend.position="none"
                                            axis.title.y=element_blank()
                                            # ,axis.text.x=element_blank()
                                            # ,axis.title.x=element_blank(),
                                            ,strip.background=element_blank()
                                            ,strip.text.x = element_blank()),
                                          labels=NULL, # c("", "2yrs", "scale", "shape", "n"),
                                          ncol=2,nrow=1,
                                          rel_widths=c(0.07,1),
                                          align='hv'),
                                labels=NULL,
                                ncol=1,nrow=4,
                                rel_heights = c(1,1,1,1),
                                align='hv')
  ggsave(plot_facets_reg_all_all, 
         filename=paste0(dir.res_maps_durat, '/all_maps_changes_durat_',sat_prod,'_reg_ALLALL.png'),
         width = 45, height =30, dpi = 150, bg = "white")
  
  
  
  
  
  # plot_facets_all_all_v =plot_grid(plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=40),
  #                                            plot_facets+
  #                                              facet_wrap(~reg_name,nrow=20,ncol=1)+
  #                                              theme(
  #                                                legend.position="none"
  #                                                ,axis.text.y=element_blank()
  #                                                ,axis.title.y=element_blank()),
  #                                            # ,strip.text.x=element_text(size=30)), 
  #                                            #,strip.background=element_blank()),
  #                                            labels=NULL,
  #                                            ncol=1,nrow=2,
  #                                            rel_heights=c(0.05,1),
  #                                            align='hv'),
  #                                  plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=40),
  #                                            plot_facets+
  #                                              facet_wrap(~reg_name,nrow=20,ncol=1)+
  #                                              theme(
  #                                                legend.position="none"
  #                                                ,axis.text.y=element_blank()
  #                                                ,axis.title.y=element_blank()),
  #                                            # ,strip.background=element_blank()),
  #                                            #,strip.text.x=element_blank()),
  #                                            labels=NULL, # c("", "2yrs", "scale", "shape", "n"),
  #                                            ncol=1,nrow=2,
  #                                            rel_heights=c(0.05,1),
  #                                            align='hv'),
  #                                  plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=40),
  #                                            plot_facets+
  #                                              facet_wrap(~reg_name,nrow=20,ncol=1)+
  #                                              theme(
  #                                                legend.position="none"
  #                                                ,axis.text.y=element_blank()
  #                                                ,axis.title.y=element_blank()),
  #                                            # ,strip.background=element_blank()),
  #                                            # ,strip.text.x = element_blank()),
  #                                            labels=NULL, # c("", "2yrs", "scale", "shape", "n"),
  #                                            ncol=1,nrow=2,
  #                                            rel_heights=c(0.05,1),
  #                                            align='hv'),
  #                                  plot_grid(ggdraw()+draw_label(sat_acron,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=40),
  #                                            plot_facets+
  #                                              facet_wrap(~reg_name,nrow=20,ncol=1)+
  #                                              theme(
  #                                                # legend.position="none"
  #                                                ,axis.text.y=element_blank()
  #                                                ,axis.title.y=element_blank()),
  #                                            # ,strip.background=element_blank()
  #                                            # ,strip.text.x = element_blank()),
  #                                            labels=NULL, # c("", "2yrs", "scale", "shape", "n"),
  #                                            ncol=1,nrow=2,
  #                                            rel_heights=c(0.05,1),
  #                                            align='hv'),
  #                                  labels=NULL,
  #                                  ncol=4,nrow=1,
  #                                  rel_widths=c(1,1,1,1.4),
  #                                  align='hv')
  # ggsave(plot_facets_all_all_v, 
  #        filename=paste0(dir.res_maps_durat, '/all_maps_changes_durat_',
  #                        sat_prod,'_bis_ALLALL_vert.png'),
  #        width=20,height=40,dpi=200,bg ="white")
  
  
  
  # save R data file of changes per duration and per region:
  saveRDS(list(df_changes_durat        = df_changes_durat,
               df_changes_durat_merged = df_changes_durat_merged, 
               plot_facets_reg         = plot_facets_reg,
               plot_facets_reg_all     = plot_facets_reg_all,
               df_changes_durat_zones  = df_changes_durat_zones,
               df_changes_durat_zones_merged = df_changes_durat_zones_merged,
               df_changes_durat_zones_merged_sort = df_changes_durat_zones_merged_sort, 
               plot_facets             = plot_facets,
               plot_facets_all         = plot_facets_all,
               plot_facets_all_zones   = plot_facets_all_zones,
               product_name            = sat_prod),
          paste0(dir.res_maps_durat, 'all_changes_durat_',sat_prod,'_zones.rds'))
  
}










