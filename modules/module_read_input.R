library(png)
library(ncdf4)
library(CFtime)
library(rnaturalearth)
library(raster)


################################################################################
load_sat_data_netcdf=function(sPathMeteoIN, 
                              namfnc, 
                              txtname,
                              sTimeName= "time",
                              sVarName= "PRE",
                              sat_prod,
                              dir_results,
                              flag.check=F){
################################################################################
  #^* GOAL: load satellite data from netcdf file
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
  # inizialize objects and variables:
  lon=lat=time=tunits=tmp_array=fillvalue=c();
  cf=timestamps=time_cf=dates_numeric=tmp_slice=c()
  
  # READ NETCDF:
  ncfname<-paste0(sPathMeteoIN, namfnc)
  print(paste0("Reading Satellite Product ",ncfname))
  # open file:
  ncin<-nc_open(ncfname)
  # print(ncin)
  # Save metadata to a text file
  {
    sink(paste0(sPathMeteoIN,txtname))
    print(ncin)
    sink()
  }
  # get longitude and latitude
  lon<-ncvar_get(ncin,"lon")
  nlon<-dim(lon)
  lat<-ncvar_get(ncin,"lat")
  nlat<-dim(lat)
  print(paste0("Size = ",nlon," lon; ",nlat," lat"))
  # get time
  time<-ncvar_get(ncin,sTimeName)
  tunits<-ncatt_get(ncin,sTimeName,"units")
  nt<-dim(time)
  # get rainfall variable
  tmp_array<-ncvar_get(ncin,sVarName)
  dlname<-ncatt_get(ncin,sVarName,"long_name")
  dunits<-ncatt_get(ncin,sVarName,"units")
  fillvalue<-ncatt_get(ncin,sVarName,"_FillValue")
  dim(tmp_array)
  # get global attributes
  title<-ncatt_get(ncin,0,"title")
  institution<-ncatt_get(ncin,0,"institution")
  datasource<-ncatt_get(ncin,0,"source")
  references<-ncatt_get(ncin,0,"references")
  history<-ncatt_get(ncin,0,"history")
  Conventions<-ncatt_get(ncin,0,"Conventions")
  # close netcdf file
  nc_close(ncin)
  
  # Time conversion to CFtime class:
  cf<-CFtime(tunits$value,calendar="proleptic_gregorian",time) 
  timestamps<-as_timestamp(cf) # get character-string times
  time_cf<-CFparse(cf,timestamps) # parse the string into date components
  # tmp_array[tmp_array==fillvalue$value,] <- NA # replace netCDF fill values with NA's
  # # head(as.vector(tmp_array[,,1]))
  # length(na.omit(as.vector(tmp_array[,,1])))
  dates_numeric=(as.numeric(as.POSIXct(strptime(timestamps, "%Y-%m-%d"))) - 
                  as.numeric(as.POSIXct(strptime(timestamps[1], "%Y-%m-%d"))))/86400
  
  # # save slice to geotiff for check:
  # ##################################
  # slice=1
  # r <- raster(tmp_array[, , slice])
  # proj4string(r)=CRS("+proj=longlat +datum=WGS84 +no_defs")
  # extent(r)=c(6, 19, 36, 48)
  # # aggregating by factor 2 (ompute the mean value):
  # writeRaster(r, filename=paste0(sPathMeteoIN, "/CMORPH_orig_0.25deg.tif"), format="GTiff", overwrite=T)
  
  
  # check one slice:
  ##################
  slice=50
  fnameplot=paste0(dir_results,"/",sat_prod,"/slice_",slice,"_",sat_prod,".png") 
  dir.create(paste0(dir_results,"/",sat_prod))
  print(paste0("saving figure with one frame of sat. data to file: ",fnameplot))
  tmp_slice<-tmp_array[,,slice]
  timestamps[slice]
  tmp_slice[tmp_slice=="NaN"]=NA
  
  # if (sat_prod==5){
  #   tmp_array<-tmp_array[,ncol(tmp_slice):1,]
  #   # tmp_array<-tmp_array[nrow(tmp_slice):1,,]
  #   lat=sort(lat)
  #   tmp_slice<-tmp_array[,,slice]
  #   tmp_slice[tmp_slice=="NaN"]=NA
  # } else if (sat_prod==3){
  #   lat=sort(lat)
  # }
  
  # add world map as sp class:
  gg_world=ne_countries(scale="medium",returnclass="sf") 
  data(wrld_simpl)
  
  # plot map just as an example:
  png(filename=fnameplot,width=600,height=600)
  par(mfrow=c(1,1))
  image(lon,lat,tmp_slice,col=rev(brewer.pal(10,"RdBu")))
  # plot(wrld_simpl,add=TRUE)
  plot(gg_world,fill=NA,color="grey70",add=T)
  dev.off()
  
  
  # image(lon[17:32],lat[35:44], tmp_slice[17:32, 35:44], col=rev(brewer.pal(10,"RdBu")))
  # data(wrld_simpl)
  # plot(wrld_simpl,add=TRUE)
   
  # tmp_slice[tmp_slice=="NaN"]=NA
  # tmp_slice[tmp_slice==0]=NA
  # image(lon,lat, tmp_slice, col=rev(brewer.pal(10,"RdBu")))
  # data(wrld_simpl)
  # plot(wrld_simpl,add=TRUE)
  
  # points(rep(lon[210], 19),lat[48:66], pch = 19, col='red', cex=1)
  # points(rep(lon[213], 16),lat[57:72], pch = 19, col='red', cex=1)
  # points(rep(lon[222], 12),lat[59:70], pch = 19, col='red', cex=1)
  # dev.off()
  
  
  if (flag.check){
    # CHECK THE PRODUCT by computing the cumulative values:
    ########################################################
    # tmp_sum =  colSums(tmp_array, dims=2, na.rm=T)
    tmp_sum= tmp_array[,,1]
    for (row in 1:nrow(tmp_array)) {
      for (col in 1:ncol(tmp_array)) {
        tmp_sum[row,col] = mean(tmp_array[row,col,], na.rm=T)
  
      }
    }
  
    fnameplot_mean= paste0(dir_results, "/", sat_prod,"/mean_", sat_prod, ".png")
    png(filename=fnameplot_mean, width=600, height=600)
    par(mfrow=c(1,1))
    image(lon,lat, tmp_sum, col=rev(brewer.pal(10,"RdBu")))
    data(wrld_simpl)
    # plot(wrld_simpl,add=TRUE)
    plot(gg_world, fill = NA, color = "grey70", add=T)
    dev.off()
  }   
  
  return(list(lon=lon,
              lat=lat, 
              timestamps=timestamps,
              dates_numeric=dates_numeric,
              fillvalue=fillvalue,
              tmp_array=tmp_array))
}
 