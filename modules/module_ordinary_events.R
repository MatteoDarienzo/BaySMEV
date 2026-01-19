library(imager)
library(lubridate)


################################################################################
event_separation_dry_spell <- function(s, separation, min_ev_duration){
################################################################################
  #^* GOAL: SEPARATION pf ordinary events, defines DRY time intervals 
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [list] s --> list containing variable values "vals", instants "time",
  #^*                    and time step in minutes "time_resolution". 
  #^*    2. [integer] separation --> min time (minutes) for event separation,ex:1440 
  #^*    3. [integer] min_ev_duration --> min duration (minutes) of events, ex:30
  #^* OUT
  #^*    1. [list] s list updated with new var "is_event", "yr", "to" and "from"
  #^****************************************************************************
  #^* REF.: translated function in R from open source matlab codes of Francesco
  #^*       Marra (University of Padua, Italy), SMEV method:
  #^*       from https://zenodo.org/records/11934843 (SMEV methodology) and
  #^*       from https://zenodo.org/records/15047817 (non-stationary framework). 
  #^*       procedure described in Marra et al., 2019. 2020
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************
  # SEPARATION defines by 'separation' DRY time intervals 
  # defines start and end of all events separated by at least
  ##############################################################################
  
  # Identify the storms:
  ######################
  # s.is_event = erode(dilate(c(rep(0,separation),s.vals>0, 
  #                    rep(0,separation)),rep(1,separation-1)),rep(1,separation-1))
  mm=as.cimg(c(rep(0,separation), s$vals>0, rep(0,separation)),
             length(c(rep(0,separation), s$vals>0, rep(0,separation))),1)
  nn = as.cimg(rep(1,separation),length(rep(1,separation)),1)
  # use dilate function to dilate the packed bin image using structuring element of ones
  m_dil=imager::dilate(mm,nn)
  # use erode function to erode the packed bin image using structuring element of ones
  s$is_event=imager::erode(m_dil,nn)
  
  # s$is_event = erode(dilate(c(rep(0,separation+1),s$vals>0,rep(0,separation)),
  #                           rep(1,separation)),rep(1,separation+1))
  # erosion(dilation(c(rep(0,separation+1), s$vals>0, rep(0,separation)),  
  #                  rep(1,separation)),rep(1,separation))

  # Remove fist and last events:
  s$is_event=s$is_event[(1+separation):(length(s$is_event)-separation)]
  # Define function to detect data sequences of (0,1) or (1,0):
  "%seq_in%" = function(b,a) which(
    sapply(1:(length(a)-length(b)+1),function(i) all(a[i:(i+length(b)-1)]==b))) 
  # where there is the sequence 0 1, take the second index
  s$from=c(0,1) %seq_in% s$is_event+1  
  # where there is the sequence 1 0, take the first index
  s$to  =c(1,0) %seq_in% s$is_event    
  
  # removes events that end after the timeseries ends (that have 1 at the last index)
  if (tail(s$is_event,1)==1){ s$from = s$from[-length(s$from)]}  
  # removes events that start before the timeseries starts (that have 1 at the first index)
  if (s$is_event[1]==1){ s$to = s$to[-1] }                       
  
  # Remove ordinary events across non-subsequent years:
  toremove = which(lubridate::year(s$time[s$to])-lubridate::year(s$time[s$from])>=2)
  if (any(toremove)){
    s$is_event[s$from[toremove]:s$to[toremove]]=0
    s$from=s$from[-toremove]
    s$to=s$to[-toremove]
  }
  
  # get the year of each data
  s$yr=lubridate::year(s$time[s$to])
  
  # Removes too short events (< min_ev_duration):
  is_short = rep(F,length(s$from))
  for (ev in 1:length(s$from)){
    is_short[ev]=(s$to[ev]-s$from[ev]+1)*s$time_resolution < min_ev_duration;
  }
  
  if (length(which(is_short==T)) >0){
    s$from=s$from[-which(is_short==T)]
    s$to=s$to[-which(is_short==T)]
    s$yr=s$yr[-which(is_short==T)]
  }
  return(s)
}




