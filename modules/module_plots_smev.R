


################################################################################
plot_rain_series<-function(time,
                           vals, 
                           title_x, title_y,
                           main_title, 
                           path2save){
################################################################################
  #^* GOAL: plot the precipitation time series
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [real] time -->
  #^*    2. [real] vals -->
  #^*    3. [character] title_x -->
  #^*    4. [character] title_y -->
  #^*    5. [character] main_title -->
  #^*    6. [character] path2save -->
  #^* OUT
  #^*    1.
  #^****************************************************************************
  #^* REF.: 
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************
  rain_series=ggplot()+
    geom_line(aes(x=as.Date(time),y=vals),color="blue",size=0.2)+
    coord_cartesian(clip="off")+
    scale_y_continuous(name=title_y,expand=c(0,0))+
    scale_x_date(expand=c(0,0))+
    labs(x=title_x,y=title_y)+
    theme_bw(base_size=11)+
    theme(
      panel.grid.minor=element_blank()
      ,panel.grid.major=element_blank()
      ,plot.title=element_text(hjust=0.5,face="bold",size=12)
    )+ggtitle(main_title)
  # save plot:
  ggsave(rain_series,filename=path2save,width=10,height=4,dpi=150)
}









################################################################################
compar_quantiles_lik_post= function(q2Lik, q2Post, dir_results){
################################################################################
  #^* GOAL: plot the precipitation time series
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN 
  #^*    1. [] q2Lik -->
  #^*    2. [] q2Post -->
  #^*    3. [character] dir_results -->
  #^* OUT
  #^*    1. [real] -->
  #^****************************************************************************
  #^* REF.: 
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************
  # PLOTS:
  png(filename=paste0(dir_results, "/hist_lik_post_q2_yy_1_a.png"),
      width=1000,height=1000)
  par(mfrow=c(2,1))
  hist(q2Lik,50,xlab='q2',main='Without Priors')
  hist(q2Post,50,main='With Priors')
  dev.off()
  png(filename=paste0(dir_results,"/hist_lik_post_q2_yy_1_b.png"),
      width=1000,height=1000)
  c1<-"blue" # rgb(173,216,230,max=255,alpha=80,names="lt.blue")
  c2<-"red"  # rgb(255,192,203,max=255,alpha=80,names="lt.pink")
  # Plot 1st histogram using a transparent color
  plot(hist(q2Lik,50,xlab='q2',main='Without Priors',plot=F),col=c1) 
  # Add 2nd histogram using different color
  plot(hist(q2Post,50,main='With Priors',plot=F),col=c2,add=T) 
  dev.off()
  
  pdf_ALL=ggplot()+
   geom_density(aes(x=q2Lik),color="blue",fill=NA,lwd=1.2)+
   geom_density(aes(x=q2Post),color="red",fill=NA,lwd=1.2)+
   theme_bw(base_size=20)
  ggsave(pdf_ALL,filename=paste0(dir_results,"/hist_lik_post_q2_yy_1_c.png"))
  pdf_ALL
}







################################################################################
plot_quantiles_gridplot <- function(plotlist,
                                    title,
                                    pathout,
                                    flag.plot=T){
################################################################################
  #^* GOAL: plot quantiles in gridplot (one plot for eahc duration)
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [integer] plotlist --> 
  #^*    2. [integer] title --> 
  #^*    3. [character] pathout --> 
  #^*    4. [character] flag.plot -->
  #^* OUT
  #^*    1.
  #^****************************************************************************
  #^* REF.: 
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************
  if (flag.plot){
    # create title object:
    title1<-ggdraw()+ 
      draw_label(title,fontface='bold',x=0.5,hjust=0.5,vjust=0.5,size=20)
    # Final plot with all plots combined:
    # glegend<-get_plot_component(spaghettiplot,'guide-box',return_all=T)
    # glegend<-get_legend2(plotlist[[1]])
    # ,legend.title=element_text(size=30),
    # ,legend.position='bottom'
    # ,axis.title=element_text(size=40)
    # ,axis.text.x=element_text(size=40)) 
    all.plots=NULL
    
    if (length(plotlist)==1){
      all.plots=plot_grid(plotlist[[1]]+theme(plot.title=element_blank()),
                          label_size=15,
                          labels=c(plotlist[[1]]$duration),
                          ncol=2,nrow=2,rel_heights=1, 
                          align='hv')
    } else if (length(plotlist)==2){
      all.plots = plot_grid(plotlist[[1]]+theme(legend.position="none",
                                                plot.title=element_blank()),
                            plotlist[[2]]+theme(plot.title=element_blank()),
                            label_size=15,
                            labels=c(plotlist[[1]]$duration,
                                     plotlist[[2]]$duration),
                            ncol=2,nrow=1,
                            rel_heights=c(1,1),
                            align='hv')
    } else if (length(plotlist)==3){
      all.plots = plot_grid(plotlist[[1]]+theme(legend.position="none", 
                                                plot.title=element_blank()),
                            plotlist[[2]]+theme(legend.position="none", 
                                                 plot.title=element_blank()),
                            plotlist[[3]]+theme(plot.title=element_blank()),
                            label_size=15,
                            labels=c(plotlist[[1]]$duration,
                                     plotlist[[2]]$duration,
                                     plotlist[[3]]$duration),
                            ncol=2,nrow=2,
                            rel_heights=c(1,1,1),
                            align='hv')
    } else if (length(plotlist)==4){
      all.plots = plot_grid( plotlist[[1]]+theme(legend.position="none",
                                                 plot.title=element_blank()),
                             plotlist[[2]]+theme(legend.position="none",
                                                 plot.title=element_blank()),
                             plotlist[[3]]+theme(legend.position="none",
                                                 plot.title=element_blank()),
                             plotlist[[4]]+theme(plot.title=element_blank()),
                             label_size=15,
                             labels=c(plotlist[[1]]$duration, 
                                      plotlist[[2]]$duration,
                                      plotlist[[3]]$duration,
                                      plotlist[[4]]$duration),
                             ncol=2,nrow=2,rel_heights=c(1,1,1,1),align='hv')
    }    
    if (!is.null(all.plots)){
      all.plots.title = plot_grid(title1,
                                  all.plots,
                                  ncol=1,
                                  nrow=2,
                                  rel_heights=c(0.08,1),
                                  align='hv')
      # all.plots.title
      ggsave(all.plots.title,filename=pathout,width=12,height=12,dpi=150,bg="white")
    }
  }
}












################################################################################
model.selection.plot <- function(criteria.df, 
                                 title,
                                 name_output){
################################################################################
  #^* GOAL: plot the information criteria with different models
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [dataframe] criteria.df  -->
  #^*    2. [character] title  -->
  #^*    3. [character] name_output  -->
  #^* OUT
  #^*    1. [real] bic.plot -->
  #^****************************************************************************
  #^* REF.: 
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************
  if (!is.null(criteria.df$WAIC)){
    # define legend items andf colors:
    col =c("AIC"="blue","BIC"="red","HQC"="purple","DIC1 (mean)"="green",
           "DIC1 (MAP)"="forestgreen","DIC2 (mean)"="orange","DIC3"="black",
           "WAIC"="cyan","LOOIC"="brown") 
    AICcol="AIC";BICcol="BIC";HQCcol="HQC";DIC1col="DIC1 (mean)";
    DIC1col_b="DIC1 (MAP)";DIC2col="DIC2 (mean)";DIC3col="DIC3";WAICcol="WAIC"
    LOOICcol="LOOIC"
    
    # plot all criteria:
    bic.plot<-ggplot(data=criteria.df) + 
      theme_light(base_size=20)+
      geom_line(aes(x=x,y=BIC,color=BICcol,group=BICcol),size=0.2,na.rm=T)+
      geom_line(aes(x=x,y=AIC,color=AICcol,group=AICcol),size=0.2,na.rm=T)+
      geom_line(aes(x=x,y=HQC,color=HQCcol,group=HQCcol),size=0.2,na.rm=T)+
      geom_line(aes(x=x,y=DIC1,color=DIC1col,group=DIC1col),size=0.2,na.rm=T)+
      geom_line(aes(x=x,y=DIC1_b,color=DIC1col_b,group=DIC1col_b),size=0.2,na.rm=T)+
      geom_line(aes(x=x,y=DIC2,color=DIC2col,group=DIC2col),size=0.2,na.rm=T)+
      geom_line(aes(x=x,y=DIC3,color=DIC3col,group=DIC3col),size=0.2,na.rm=T)+
      geom_line(aes(x=x,y=WAIC,color=WAICcol,group=WAICcol),size=0.2,na.rm=T)+
      geom_line(aes(x=x,y=LOOIC,color=LOOICcol,group=LOOICcol),size=0.2,na.rm=T)+
      geom_point(aes(x=x,y=BIC,color=BICcol,group=BICcol),shape=0,size=5,na.rm=T)+
      geom_point(aes(x=x,y=AIC,color=AICcol,group=AICcol),shape=2,size=5,na.rm=T)+
      geom_point(aes(x=x,y=HQC,color=HQCcol,group=HQCcol),shape=6,size=5,na.rm=T)+
      geom_point(aes(x=x,y=DIC1,color=DIC1col,group=DIC1col),shape=5,size=5, na.rm=T)+
      geom_point(aes(x=x,y=DIC1_b,color=DIC1col_b,group=DIC1col_b),shape=5,size=5,na.rm=T) +
      geom_point(aes(x=x,y=DIC2,color=DIC2col,group=DIC2col),shape=5,size=5,na.rm=T) +
      geom_point(aes(x=x,y=DIC3,color=DIC3col,group=DIC3col),shape=5,size=5,na.rm=T) +
      geom_point(aes(x=x,y=WAIC,color=WAICcol,group=WAICcol),shape=4,size=5,na.rm=T) +
      geom_point(aes(x=x,y=LOOIC,color=LOOICcol,group=LOOICcol),shape=3,size=5,na.rm=T) +
      # points with min value:
      geom_point(aes(x=BICmin,y=BIC[BICmin]),shape=15,size=6,
                 color="red",na.rm=T)+
      geom_point(aes(x=AICmin,y=AIC[AICmin]),shape=17,size=6,
                 color="blue", na.rm=T) +
      geom_point(aes(x=HQCmin,y=HQC[HQCmin]),shape=25,size=6,
                 color="purple",fill="purple",na.rm=T)+
      geom_point(aes(x=DIC1min,y=DIC1[DIC1min]),shape=23,size=6,
                 color="green",fill="green",na.rm=T)+
      geom_point(aes(x=DIC1min_b,y=DIC1_b[DIC1min_b]),shape=23,size=6,
                 color="forestgreen",fill="forestgreen",na.rm=T)+
      geom_point(aes(x=DIC2min,y=DIC2[DIC2min]),shape=23,size=6, 
                 color="orange",fill="orange",na.rm=T)+
      geom_point(aes(x=DIC3min,y=DIC3[DIC3min]),shape=23,size=6,
                 color="black",fill="black",na.rm=T)+
      geom_point(aes(x=WAICmin,y=WAIC[WAICmin]),shape=4,size=6, 
                 color="cyan",fill="cyan",na.rm=T)+
      geom_point(aes(x=LOOICmin,y=WAIC[LOOICmin]),shape=3,size=6,
                 color="brown",fill="brown",na.rm=T)+
      scale_x_discrete(name="SMEV models",expand=c(0,0))+
      #scale_y_continuous(name="BIC,AIC, QC,DIC",expand=c(0,0), 
      #                   limits=c(min_grid,max_grid),
      #                   breaks=seq(min_grid,max_grid,50))+ 
      #                   # breaks=scales::pretty_breaks(n=2))+
      scale_y_continuous(name="Criterion value",expand=c(0,0)) +
      xlab("Number of segments") + 
      ylab("Model selection criteria")+
      coord_cartesian(clip='off')+
      scale_colour_manual(name=element_blank(),
                          values=col, 
                          breaks=c(AICcol,BICcol,HQCcol,DIC1col,DIC1col_b,
                                   DIC2col,DIC3col,WAICcol,LOOICcol),
                          guide = guide_legend(override.aes = list(
                            linetype=c("solid","solid","solid","solid","solid", 
                                       "solid","solid","solid","solid"),
                            shape=c(2,0,6,5,5,5,5,4,3)),
                            size=c(5,5,5,5,5,5,5,5,5)))+
      ggtitle(title)+
      theme(   plot.title=element_text(hjust=0.5,size=12) 
              ,axis.text=element_text(size=13) 
              ,axis.title=element_text(size=13)
              ,plot.margin=margin(c(1,1,1,1),unit="cm")
              ,panel.grid.major=element_line(size=0.3,linetype="dashed")
              ,panel.grid.minor=element_blank()
              #,panel.grid.major=element_blank()
              #,panel.grid.minor=element_blank()
              ,panel.background=element_rect(fill="transparent") 
              #,axis.line=element_line(colour="black")
              ,axis.ticks=element_line(colour="black")
              ,text=element_text(size=13)
              ,legend.key=element_rect(colour="transparent",fill="transparent")
              ,legend.background=element_rect(colour="transparent",fill="transparent")
              #,legend.title=element_text(colour="blue",size=16,face="bold"),
              ,legend.position="bottom"
              ,legend.text=element_text(size=7,margin=margin(r=0,l=0,unit="pt"))
              ,legend.box="horizontal")+
      guides(colour=guide_legend(nrow=1))
  } else {
    # define legend items:
    col=c("AIC"="blue","BIC"="red","HQC"="purple","DIC1 (mean)"="green",
          "DIC1 (MAP)"="forestgreen","DIC2 (mean)"="orange","DIC3"="black")  
    AICcol="AIC";BICcol="BIC";HQCcol="HQC";DIC1col="DIC1 (mean)";
    DIC1col_b="DIC1 (MAP)";DIC2col="DIC2 (mean)";DIC3col="DIC3"
    
    # plot all criteria:
    bic.plot<-ggplot(data=criteria.df)+ 
      theme_light(base_size=20)+
      # all lines between points:
      geom_line(aes(x=x,y=BIC,color=BICcol,group=BICcol),size=0.2,na.rm=T) +
      geom_line(aes(x=x,y=AIC,color=AICcol,group=AICcol),size=0.2,na.rm=T) +
      geom_line(aes(x=x,y=HQC,color=HQCcol,group=HQCcol),size=0.2,na.rm=T) +
      geom_line(aes(x=x,y=DIC1,color=DIC1col,group=DIC1col),size=0.2,na.rm=T) +
      geom_line(aes(x=x,y=DIC1_b,color=DIC1col_b,group=DIC1col_b),size=0.2,na.rm=T) +
      geom_line(aes(x=x,y=DIC2,color=DIC2col,group=DIC2col),size=0.2,na.rm=T) +
      geom_line(aes(x=x,y=DIC3,color=DIC3col,group=DIC3col),size=0.2,na.rm=T) +
      geom_point(aes(x=x,y=BIC,color=BICcol,group=BICcol),shape=0,size=5,na.rm=T) +
      geom_point(aes(x=x,y=AIC,color=AICcol,group=AICcol),shape=2,size=5,na.rm=T) +
      geom_point(aes(x=x,y=HQC,color=HQCcol,group=HQCcol),shape=6,size=5,na.rm=T) +
      geom_point(aes(x=x,y=DIC1,color=DIC1col,group=DIC1col),shape=5,size=5,na.rm=T) +
      geom_point(aes(x=x,y=DIC1_b,color=DIC1col_b,group=DIC1col_b),shape=5,size=5,na.rm=T) +
      geom_point(aes(x=x,y=DIC2,color=DIC2col,group=DIC2col),shape=5,size=5,na.rm=T) +
      geom_point(aes(x=x,y=DIC3,color=DIC3col,group=DIC3col),shape=5,size=5,na.rm=T) +
      # points with min value:
      geom_point(aes(x=BICmin,y=BIC[BICmin]),shape=15,size=6,
                 color="red",na.rm=T)+
      geom_point(aes(x=AICmin,y=AIC[AICmin]),shape=17,size=6,
                 color="blue",na.rm=T)+
      geom_point(aes(x=HQCmin,y=HQC[HQCmin]),shape=25,size=6,
                 color="purple",fill="purple",na.rm=T)+
      geom_point(aes(x=DIC1min,y=DIC1[DIC1min]),shape=23,size=6,
                 color="green",fill="green",na.rm=T)+
      geom_point(aes(x=DIC1min_b,y=DIC1_b[DIC1min_b]),shape=23,size=6,
                 color="forestgreen",fill="forestgreen",na.rm=T)+
      geom_point(aes(x=DIC2min,y=DIC2[DIC2min]),shape=23,size=6,
                 color="orange",fill="orange",na.rm=T)+
      geom_point(aes(x=DIC3min,y=DIC3[DIC3min]),shape=23,size=6,
                 color="black",fill="black",na.rm=T)+
      scale_x_discrete(name="SMEV models", expand=c(0,0)) +
      #scale_y_continuous(name="BIC,AIC,HQC,DIC",expand=c(0,0),
      #                   limits=c(min_grid,max_grid), 
      #                   breaks=seq(min_grid,max_grid,50))+
      #                   # breaks = scales::pretty_breaks(n = 2)) + 
      scale_y_continuous(name="Criterion value",expand=c(0,0))+
      xlab("Number of segments")+
      ylab("Model selection criteria")+
      coord_cartesian(clip='off')+
      scale_colour_manual(name=element_blank(),
                          values=col, 
                          breaks=c(AICcol,BICcol,HQCcol,DIC1col,DIC1col_b,
                                   DIC2col,DIC3col),
                          guide=guide_legend(override.aes=list(
                            linetype=c("solid","solid","solid","solid", 
                                       "solid","solid","solid"),
                            shape=c(2,0,6,5,5,5,5)),
                            size=c(5,5,5,5,5,5,5)))+
      ggtitle(title)+
      theme(plot.title=element_text(hjust=0.5,size=12) 
            ,axis.text=element_text(size=13) 
            ,axis.title=element_text(size=13)
            ,plot.margin=margin(c(1,1,1,1),unit="cm")
            ,panel.grid.major=element_line(size=0.3,linetype="dashed")
            ,panel.grid.minor=element_blank()
            #,panel.grid.major=element_blank(), 
            #,panel.grid.minor=element_blank(),
            ,panel.background=element_rect(fill="transparent") 
            ,axis.ticks=element_line(colour="black")
            ,text=element_text(size=13)
            ,legend.key=element_rect(colour="transparent",fill="transparent")
            ,legend.background=element_rect(colour="transparent",fill="transparent")
            #,axis.line= element_line(colour = "black"),
            #,legend.title= element_text(colour="blue", size=16, face="bold"),
            ,legend.position="bottom"
            ,legend.text=element_text(size=7,margin=margin(r=0,l=0,unit="pt"))
            ,legend.box="horizontal")+
      guides(colour=guide_legend(nrow=1))
  }
  ggsave(bic.plot,filename=name_output,width=8,height=5,dpi=110)
  # pdf(paste0(dir.res,"/Criteria_",name_output,".pdf"),useDingbats=F,
  #     width=6,height=4,dpi=100)
  # print(bic.plot)
  # dev.off()
  return(bic.plot)
}









################################################################################
bayesfactor.plot <- function(bayes_factor.df, # values of BF
                             title,
                             name_output){
################################################################################
  #^* GOAL: plot the Bayes Factor between two models
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [dataframe] bayes_factor.df --> 
  #^*    2. [character] title --> 
  #^*    3. [character] name_output --> 
  #^* OUT
  #^*    1. [plot] BF.plot --> 
  #^****************************************************************************
  #^* REF.: 
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  # how to read BF values (Jeffreys scale of evidence):
  # Bij < 1/100       ==> extreme evidence for Mj
  # 1/100< Bij < 1/30 ==> very strong evidence for Mj
  # 1/30< Bij < 1/10  ==> strong evidence for Mj
  # 1/10< Bij < 1/3   ==> moderate evidence for Mj
  # 1/3 < Bij < 1     ==> weak evidence for Mj
  # 1   < Bij < 3     ==> weak evidence for Mi
  # 3   < Bij < 10    ==> moderate evidence for Mi
  # 10  < Bij < 30    ==> strong evidence for Mi
  # 30  < Bij < 100   ==> very strong evidence for Mi
  # Bij > 100         ==> extreme evidence for Mi
  # or other classification here:
  # https://journals.sagepub.com/doi/epub/10.1177/0149206314560412
  #^****************************************************************************
  BF.plot <- ggplot(data=bayes_factor.df)+ 
    geom_hline(yintercept=c(1/30,1/10,1/3,1,3,10),color="red",linetype="dashed")+
    annotate("text",x=1,y=11,label="Strong evidence of 1 vs 2",color="pink1",size=5)+
    annotate("text",x=1,y=7,label="Moderate evidence of 1 vs 2",color="pink1",size=5)+
    annotate("text",x=1,y=2,label="Weak evidence of 1 vs 2",color="pink1",size=5)+
    annotate("text",x=1,y=0.5,label="Weak evidence of 2 vs 1",color="pink1",size=5)+
    annotate("text",x=1,y=0.2,label="Moderate evidence of 2 vs 1",color="pink1",size=5)+
    annotate("text",x=1,y=1/20,label="Strong evidence of 2 vs 1",color="pink1",size=5)+
    annotate("text",x=1,y=1/100,label="Extreme evidence of 2 vs 1",color="pink1",size=5)+
    geom_point(aes(x=x,y=BF),shape=18,size=5,na.rm=T)+
    scale_x_discrete(name="Compared SMEV models") +
    scale_y_log10(name="Bayes Factor value")+ # in log scale !!!
    annotation_logticks(base=10,sides="l",scaled=T,colour="black",size=0.3,linetype=1)+
    coord_cartesian(clip='off')+
    theme_light(base_size=20)+
    theme(plot.title=element_text(hjust=0.5,size=12)
          ,axis.text=element_text(size=12) 
          ,axis.title=element_text(size=12)
          ,plot.margin=margin(c(1,1,1,1),unit="cm")
          ,panel.grid.major=element_blank()
          ,panel.grid.minor=element_blank()
          ,panel.background=element_rect(fill="transparent") 
          ,axis.ticks=element_line(colour="black")
          ,text=element_text(size=13)
          ,legend.key=element_rect(colour="transparent",fill="transparent")
          ,legend.background=element_rect(colour="transparent",fill="transparent")
          #,axis.line=element_line(colour="black"),
          #,legend.title=element_text(colour="blue",size=16,face="bold"),
          ,legend.position="bottom"
          ,legend.text=element_text(size=9)
          ,legend.box="horizontal")+ggtitle(title)
  ggsave(BF.plot,filename=name_output,width=8,height=5,dpi=110)
  return(list(BF.plot))
}










################################################################################
trendAMAX.plot <- function(df_trend_AMAX, # df with the AMAX series
                           res_MK,
                           senss,
                           title,
                           pathout){
################################################################################
  #^* GOAL:  Analyse and plot trends in rainfall time series (AMAX in this case)
  #^*        MK trend analysis + Sen's Slope test
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
  mksen_AMAX=ggplot(
    data=df_trend_AMAX,
    mapping=aes(x=year(as.Date(as.character(date),format="%Y")),y=y))+
    geom_point()+
    geom_line(linewidth=0.1)+
    geom_smooth(method='lm',se=T,alpha=0.1,fill="blue")+ # sen's estimator
    annotate("text",
             x=year(as.Date(as.character(tail(df_trend_AMAX$date,1)-3),format="%Y")),  
             y=max(df_trend_AMAX$y)-(max(df_trend_AMAX$y)-min(df_trend_AMAX$y))*1/10, 
             label=paste0("z: ",round(res_MK$statistic[["z"]],3),
                          "\ntau: ",round(res_MK$estimates[["tau"]],3),
                          "\np-value: ",round(res_MK$p.value,3),
                          "\nSen's slope: ",round(senss$estimates[["Sen's slope"]],3)),  
             color="red",size=4)+
    theme_bw(base_size=13)+
    xlab("Year")+ylab("Rainfall [mm/h]")+
    coord_cartesian(clip='off')+
    ggtitle(title)+
    theme(
       panel.grid.major=element_blank()
      ,panel.grid.minor=element_blank()
      ,axis.text=element_text(size=12) 
      ,plot.title=element_text(hjust=0.5,size=12)
      ,plot.margin=margin(c(1,1,1,1),unit="cm"))
  ggsave(mksen_AMAX,filename=pathout,width=8,height=4,dpi=110)
  return(mksen_AMAX)
}







################################################################################
trend_n.plot <- function(df_trend_n, # df with the AMAX series
                           res_MK,
                           senss,
                           title,
                           pathout){
################################################################################
  #^* GOAL:  Analyse and plot trends in rainfall time series ("n" in this case)
  #^*        MK trend analysis + Sen's Slope test
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [dataframe] df_trend_n -->
  #^*    2. [] res_MK --> MK estimates
  #^*    3. [] senss --> Sen's slope estimates
  #^*    4. [character] title --> 
  #^*    5. [character] pathout -->
  #^* OUT
  #^*    1. [real] mksen_n
  #^****************************************************************************
  #^* REF.:
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************
  mksen_n=ggplot(
    data=df_trend_n, 
    mapping=aes(x=year(as.Date(as.character(date),format="%Y")),y=y))+
    geom_point()+
    geom_line(linewidth=0.1)+
    geom_smooth(method='lm',se=T,alpha=0.1,fill="blue")+ # sen's estimator
    annotate("text",
             x=year(as.Date(as.character(tail(df_trend_n$date,1)-3),format="%Y")),  
             y=max(df_trend_n$y)-(max(df_trend_n$y)-min(df_trend_n$y))*1/10, 
             label=paste0("z: ",round(res_MK$statistic[["z"]],3),
                          "\ntau: ",round(res_MK$estimates[["tau"]],3),
                          "\np-value: ",round(res_MK$p.value,3),
                          "\nSen's slope: ",round(senss$estimates[["Sen's slope"]],3)),  
             color="red",size=4)+
    theme_bw(base_size=13)+
    xlab("Year")+ylab("Number of ordinary events, n")+
    coord_cartesian(clip='off')+
    ggtitle(title)+
    theme(
      panel.grid.major=element_blank()
      ,panel.grid.minor=element_blank()
      ,axis.text=element_text(size=12) 
      ,plot.title=element_text(hjust=0.5,size=12)
      ,plot.margin=margin(c(1,1,1,1),unit="cm"))
  ggsave(mksen_n,filename=pathout,width=8,height =4,dpi=110)
  return(mksen_n)
}
















################################################################################
plot_spaghetti_smev<- function(flg_save=T,
                               s,     # object containing all SMEV computations
                               quant, # object containing quantiles
                               id=1,
                               dir.res_pix){
################################################################################
  #^* GOAL: plot the SPAGHETTI from SMEV Bayesian
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [logical] flg_save --> flag to activate saving plots, otherwise just
  #^*                              save the duration
  #^*    2. [list of lists] s --> main object containing all SMEV results
  #^*    3. [list] quant --> object containing quantiles
  #^*    4. [integer] id --> duration id
  #^*    5. [character] dir.res_pix --> path of folder where to save plots
  #^* OUT
  #^*    1. [list] spaghettiplot_smev -->list object containing quantiles 
  #^*                                    (spaghetti) plots
  #^****************************************************************************
  #^* REF.: 
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************
  if (flg_save==T){
    # define text and colors for legend items:
    spaghetti_ns="SpaghettiPost (SMEV nonstat)"
    mediancurve_ns="MedianPost (SMEV nonstat)"
    Maxlikecurve_ns="MaxLik (SMEV nonstat)"
    MAPcurve_ns="MaxPost (SMEV nonstat)"
    uncert90_ns="90% unc.Post (SMEV nonstat)"
    uncert90_stat="90% unc.Post (SMEV stat)"
    uncert90_gev= "90% unc.Post (GEV stat)"
    MAPcurve_partns="MaxPost (SMEV part.nonstat)"
    MAPcurve_stat="MaxPost (SMEV stat)"
    Maxlikecurve_ns_fmins="MaxLik fmins (SMEV nonstat)"
    Maxlikecurve_partns_fmins="MaxLik fmins (SMEV part.nonstat)"
    Maxlikecurve_stat_fmins="MaxLik fmins (SMEV stat)"
    MAPcurvegev_stat="MaxPost (GEV stat)"
    empir="Empir. return levels"
    color.curves=color.breaks=color.linetype=color.shape=color.Utot=c()
    
    # plot:
    spaghettiplot_smev= ggplot()
    if (!is.null(quant$spaghetti_df_melt_ns_post)){
      # plot only spaghetti of NS smev model:
      spaghettiplot_smev=spaghettiplot_smev+
      geom_line(data=quant$spaghetti_df_melt_ns_post,
                aes(x=grid,y=value,color=spaghetti_ns,group=series),
                linetype="solid",size=0.05)
      color.curves=c(color.curves,"SpaghettiPost (SMEV nonstat)"="gray60")
      color.breaks=c(color.breaks,spaghetti_ns)
      color.linetype=c(color.linetype,"solid")
      color.shape=c(color.shape, NA)
    }        
    # non-stationary SMEV:
    if (!is.null(quant$spaghetti_qnt_df_ns_post)){
      spaghettiplot_smev=spaghettiplot_smev+
        geom_line(data=quant$spaghetti_qnt_df_ns_post, 
                  aes(x=x,y=med_ns_post,color=mediancurve_ns),
                  linetype="solid",size=1)
      color.curves=c(color.curves,"MedianPost (SMEV nonstat)"="forestgreen")
      color.breaks=c(color.breaks,mediancurve_ns)
      color.linetype=c(color.linetype,"solid")
      color.shape=c(color.shape, NA)
    }
    if (!is.null(quant$spaghetti_qnt_df_ns_post)){
      spaghettiplot_smev=spaghettiplot_smev+
        geom_line(data=quant$spaghetti_qnt_df_ns_post, 
                  aes(x=x,y=low_ns_post,color=uncert90_ns),
                  linetype="dashed",size=0.5)+
        geom_line(data=quant$spaghetti_qnt_df_ns_post, 
                  aes(x=x,y=high_ns_post),
                  color='blue',linetype="dashed",size=0.5)
       # geom_ribbon(data=quant$spaghetti_qnt_df,
       #              aes(x=x,ymin=low,ymax=high),fill="blue",alpha=0.1)+
      color.curves=c(color.curves,"90% unc.Post (SMEV nonstat)"="blue")
      color.breaks=c(color.breaks,uncert90_ns)
      color.linetype=c(color.linetype, "dashed")
      color.shape=c(color.shape, NA)
    }
    if (!is.null(quant$spaghetti_qnt_df_ns_lik)){
      spaghettiplot_smev=spaghettiplot_smev+
        geom_line(data=quant$spaghetti_qnt_df_ns_lik,
                  aes(x=x,y=MAP_ns_lik,color=Maxlikecurve_ns),size=0.7)
      color.curves=c(color.curves,"MaxLik (SMEV nonstat)"="orange")
      color.breaks=c(color.breaks,Maxlikecurve_ns)
      color.linetype=c(color.linetype,"solid")
      color.shape=c(color.shape,NA)
    }
    if (!is.null(quant$spaghetti_qnt_df_ns_post)){
      spaghettiplot_smev=spaghettiplot_smev+
        geom_line(data=quant$spaghetti_qnt_df_ns_post, 
                  aes(x=x,y=MAP_ns_post,color=MAPcurve_ns),size=0.7)
      color.curves=c(color.curves,"MaxPost (SMEV nonstat)"="blue")
      color.breaks=c(color.breaks,MAPcurve_ns)
      color.linetype=c(color.linetype,"solid")
      color.shape=c(color.shape,NA)
    }
    # stationary SMEV:
    if (!is.null(quant$spaghetti_qnt_df_stat_post)){
      spaghettiplot_smev=spaghettiplot_smev+
        geom_line(data=quant$spaghetti_qnt_df_stat_post,
                  aes(x=x,y=low_stat_post,color=uncert90_stat),
                  linetype="dashed",size=0.5)+
        geom_line(data=quant$spaghetti_qnt_df_stat_post, 
                  aes(x=x,y=high_stat_post),
                  color='red',linetype="dashed",size=0.5)
      color.curves=c(color.curves,"90% unc.Post (SMEV stat)"="red")
      color.breaks=c(color.breaks,uncert90_stat)
      color.linetype=c(color.linetype,"dashed")
      color.shape=c(color.shape,NA)
    }
    if (!is.null(quant$spaghetti_qnt_df_stat_post)){
      spaghettiplot_smev=spaghettiplot_smev+
        geom_line(data=quant$spaghetti_qnt_df_stat_post,
                  aes(x=x,y=MAP_stat_post,color=MAPcurve_stat), 
                  linetype="solid",size=0.7)
      color.curves=c(color.curves,"MaxPost (SMEV stat)"="red")
      color.breaks=c(color.breaks,MAPcurve_stat)
      color.linetype=c(color.linetype,"solid")
      color.shape=c(color.shape,NA)
    }
    # Partial non-stationary SMEV:
    if (!is.null(quant$spaghetti_qnt_df_partns_post)){
      spaghettiplot_smev=spaghettiplot_smev+
        geom_line(data=quant$spaghetti_qnt_df_partns_post,
                  aes(x=x,y=MAP_partns_post,color=MAPcurve_partns),
                  linetype="solid",size=0.7)
      color.curves=c(color.curves,"MaxPost (SMEV part.nonstat)"="purple")
      color.breaks=c(color.breaks,MAPcurve_partns)
      color.linetype=c(color.linetype,"solid")
      color.shape=c(color.shape,NA)
    }
    # stationary GEV:
    if (!is.null(quant$spaghetti_qnt_df_gev_post)){
      spaghettiplot_smev=spaghettiplot_smev+
        geom_line(data=quant$spaghetti_qnt_df_gev_post, 
                  aes(x=x,y=low_gev_post,color=uncert90_gev),
                  linetype="dashed",size=0.5)+
        geom_line(data=quant$spaghetti_qnt_df_gev_post, 
                  aes(x=x,y=high_gev_post),color='green',
                  linetype="dashed",size=0.5)
      color.curves=c(color.curves,"90% unc.Post (GEV stat)"="green")
      color.breaks=c(color.breaks,uncert90_gev)
      color.linetype=c(color.linetype,"dashed")
      color.shape=c(color.shape,NA)
    }
    if (!is.null(quant$spaghetti_qnt_df_gev_post)){
      spaghettiplot_smev=spaghettiplot_smev+
        geom_line(data=quant$spaghetti_qnt_df_gev_post, 
                  aes(x=x,y=MAP_gev_post,color=MAPcurvegev_stat),
                  linetype="solid",size=0.7)
      color.curves=c(color.curves,"MaxPost (GEV stat)"="green")
      color.breaks=c(color.breaks,MAPcurvegev_stat)
      color.linetype=c(color.linetype,"solid")
      color.shape =c(color.shape, NA)
    }
    # maxlik with optimisation "fminsearch function": 
    if (!is.null(quant$qCurve_maxlik_fmins_stat)){
      spaghettiplot_smev=spaghettiplot_smev+
        geom_line(aes(x=quant$Tgrid,y=quant$qCurve_maxlik_fmins_stat,
                      color=Maxlikecurve_stat_fmins),
                  linetype="solid",size=0.7)
      color.curves=c(color.curves,"MaxLik fmins (SMEV stat)"="black")
      color.breaks=c(color.breaks,Maxlikecurve_stat_fmins)
      color.linetype=c(color.linetype,"solid")
      color.shape =c(color.shape, NA)
    }
    if (!is.null(quant$qCurve_maxlik_fmins_partns)){
      spaghettiplot_smev=spaghettiplot_smev+
        geom_line(aes(x=quant$Tgrid,y=quant$qCurve_maxlik_fmins_partns,
                      color=Maxlikecurve_partns_fmins), 
                  linetype="solid",size=0.7)
      color.curves=c(color.curves,"MaxLik fmins (SMEV part.nonstat)"="brown")
      color.breaks=c(color.breaks,Maxlikecurve_partns_fmins)
      color.linetype=c(color.linetype,"solid")
      color.shape=c(color.shape,NA)
    }
    if (!is.null(quant$qCurve_maxlik_fmins_ns)){
      spaghettiplot_smev=spaghettiplot_smev+
        geom_line(aes(x=quant$Tgrid,y=quant$qCurve_maxlik_fmins_ns,
                      color=Maxlikecurve_ns_fmins),
                  linetype="solid",size=0.7)
      color.curves=c(color.curves,"MaxLik fmins (SMEV nonstat)"="pink")
      color.breaks=c(color.breaks,Maxlikecurve_ns_fmins)
      color.linetype=c(color.linetype,"solid")
      color.shape=c(color.shape,NA)
    }
    # emprical return levels:
    if (!is.null(s$empir_qnt_df[[id]])){
      spaghettiplot_smev=spaghettiplot_smev+
        geom_point(data=s$empir_qnt_df[[id]],
                   aes(x=x,y=y,color=empir),shape=21,fill="black",size=2)
      color.curves=c(color.curves,"Empir. return levels"="black")
      color.breaks=c(color.breaks,empir)
      color.linetype=c(color.linetype,"blank")
      color.shape=c(color.shape,21)
    }
    
    
    if ((!is.null(s$row))|(!is.null(s$col))|(!is.null(s$sat_prod))){
      nfplot=paste0('/spaghetti_post_qnt_',s$sat_prod,'_',s$row,'_',s$col,
                    '_', s$durations[id]/s$time_resolution, 'h.png')
      title_plot=paste0('Quantiles spaghetti -- duration=',
             s$durations[id]/s$time_resolution, 'h \n (', s$sat_prod, ', ', 
             year(s$time[1]),' - ',year( tail(s$time,1)),
             ', pix:', s$row,'-', s$col,')')
    } else {
      nfplot=paste0('/spaghetti_post_qnt_',s$name_project,
                    '_', s$durations[id]/s$time_resolution, 'h.png')
      title_plot=paste0('Quantiles spaghetti -- duration=',
                        s$durations[id]/s$time_resolution, 'h \n (',
                        s$name_project, ', ', 
                        year(s$time[1]),' - ',year(tail(s$time,1)),')')
    }
    
    spaghettiplot_smev=spaghettiplot_smev+
      xlab("Return period T [years]")+
      ylab("Rainfall [mm/h]")+
      theme_bw(base_size=12)+
      geom_line()+
      scale_y_continuous(expand=c(0,0))+
      # scale_x_log10(breaks=c(1,2,5,10,20,50,100),expand=c(0,0),limits=c(1,100))+
      # scale_x_log10(breaks=c(1,2,5,10,50),expand=c(0,0),limits=c(1,50))+
      scale_x_log10(limits= c(1,100),expand=c(0,0))+
      annotation_logticks(base=10,sides="b",scaled=T,colour="black",
                          size=0.3,linetype=1)+
      coord_cartesian(clip='off')+
      scale_colour_manual(name=element_blank(),
                          values=color.curves,
                          breaks=color.breaks,
                          labels=color.breaks, 
                          guide=guide_legend(override.aes=list(
                            linetype=color.linetype,
                            shape=color.shape))) +
      theme(panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.title=element_blank()
            ,legend.position=c(0.28, 0.7)
            ,legend.text=element_text(size=10,
                                      margin=margin(r=0,t=1,b=0,l=1,unit="pt"))
            ,plot.title=element_text(hjust=0.5)
            ,plot.margin=margin(c(0.5,1,0.5,1),unit="cm")
            ,legend.key.size=grid::unit(0.9,'lines'))+
      ggtitle(title_plot)
    spaghettiplot_smev$duration=paste0(s$durations[id]/60,'h')
    # save plot:
    ggsave(spaghettiplot_smev,
           filename=paste0(dir.res_pix,nfplot),
           width=8,height=6,dpi=130)
  } else {
    # no plot:
    spaghettiplot_smev=NULL
    spaghettiplot_smev$duration=paste0(s$durations[id]/60,'h')
  }
  return(spaghettiplot_smev)
}







################################################################################
plot_quantiles_smev<- function(flg_save=T,
                               s,           # object containing all SMEV results
                               quant,       # object containing quantiles
                               id,          # duration id
                               dir.res_pix){# directory for saving plot
################################################################################
  #^* GOAL: plot the precipitation time series
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [logical] flg_save --> flag to activate saving plots, otherwise just
  #^*                              save the duration
  #^*    2. [list of lists] s --> main object containing all SMEV results
  #^*    3. [list] quant --> object containing quantiles
  #^*    4. [integer] id --> duration id
  #^*    5. [character] dir.res_pix --> path of folder where to save plots
  #^* OUT
  #^*    1. [list] quantileplot_smev --> list object containing quantiles
  #^*                                    (uncertainty ribbons) plots
  #^****************************************************************************
  #^* REF.: 
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************
  # QUANTILES PLOT
  if (flg_save==T){
    # define text and colors for legend items:
    mediancurve_ns="MedianPost (SMEV nonstat)"
    MAPcurve_ns="MaxPost (SMEV nonstat)"
    MAPcurve_stat="MaxPost (SMEV stat)"
    MAPcurvegev_stat="MaxPost (GEV stat)"
    utot90_ns="90% unc.Post (SMEV nonstat)"
    utot95_ns="95% unc.Post (SMEV nonstat)"
    utot90_stat="90% unc.Post (SMEV stat)"
    utot90_gev="90% unc.Post (GEV stat)"
    empir="Empir. return levels"
    color.curves=color.breaks=color.linetype=color.shape=color.Utot=c()
    # PLOT:
    quantileplot_smev= ggplot()
    # Full non-stationary SMEV:
    if (!is.null(quant$spaghetti_qnt_df_ns_post)){
      quantileplot_smev=quantileplot_smev+
        geom_line(data=quant$spaghetti_qnt_df_ns_post, 
                  aes(x=x,y=med_ns_post,color=mediancurve_ns),
                  linetype="solid", size=1)+
        geom_ribbon(data=quant$spaghetti_qnt_df_ns_post, 
                    aes(x=x,ymin=low2.5_ns_post,ymax=high97.5_ns_post,
                        fill=utot95_ns),alpha=0.1) +
        geom_ribbon(data=quant$spaghetti_qnt_df_ns_post,
                    aes(x=x,ymin=low_ns_post, ymax=high_ns_post, 
                        fill=utot90_ns), alpha=0.2) +
        geom_line(data=quant$spaghetti_qnt_df_ns_post, 
                  aes(x=x,y=low_ns_post, color=utot90_ns),
                  linetype="dashed", size=0.5)+
        geom_line(data=quant$spaghetti_qnt_df_ns_post, 
                  aes(x=x,y=high_ns_post),color='blue', 
                  linetype="dashed", size=0.5)+
        geom_line(data=quant$spaghetti_qnt_df_ns_post,
                  aes(x=x,y=MAP_ns_post,color=MAPcurve_ns),size=0.7)
        # geom_ribbon(data=quant$spaghetti_qnt_df, 
        #             aes(x=x,ymin=low,ymax=high),fill="blue",alpha=0.1)
      color.curves=c(color.curves,"MedianPost (SMEV nonstat)"  ="forestgreen",
                                  "MaxPost (SMEV nonstat)"     ="blue",
                                  "90% unc.Post (SMEV nonstat)"="blue")
      color.breaks= c(color.breaks, mediancurve_ns, MAPcurve_ns, utot90_ns)
      color.linetype =c(color.linetype, "solid", "solid", "dashed")
      color.shape =c(color.shape,NA,NA,NA)
      color.Utot=c(color.Utot,
                   "95% unc.Post (SMEV nonstat)"="blue",
                   "90% unc.Post (SMEV nonstat)"="blue")
    }
    # Stationary SMEV:
    if (!is.null(quant$spaghetti_qnt_df_stat_post)){
      quantileplot_smev=quantileplot_smev+
        # geom_ribbon(data=quant$spaghetti_qnt_df_stat_post,
        #             aes(x=x, ymin=low_stat_post, ymax=high_stat_post, 
        #                 fill=utot90_stat), alpha=0.2)+
        geom_line(data=quant$spaghetti_qnt_df_stat_post, 
                  aes(x=x, y=low_stat_post, color=utot90_stat),
                  linetype="dashed", size=0.5)+
        geom_line(data=quant$spaghetti_qnt_df_stat_post,
                  aes(x=x, y=high_stat_post),color='red',
                  linetype="dashed", size=0.5)+
        geom_line(data=quant$spaghetti_qnt_df_stat_post,
                  aes(x=x, y=MAP_stat_post, color=MAPcurve_stat), linetype="solid", size=0.7)
      color.curves=c(color.curves, "MaxPost (SMEV stat)"         = "red",
                                   "90% unc.Post (SMEV stat)"    = "red")
      color.breaks= c(color.breaks,MAPcurve_stat,utot90_stat)
      color.linetype =c(color.linetype,"solid","dashed")
      color.shape =c(color.shape,NA,NA)
      color.Utot=c(color.Utot,"90% unc.Post (SMEV stat)"="red")
    }
    # Stationary GEV:
    if (!is.null(quant$spaghetti_qnt_df_gev_post)){
      quantileplot_smev=quantileplot_smev+
        # geom_ribbon(data=quant$spaghetti_qnt_df_gev_post,   
        #             aes(x=x, ymin=low_gev_post, ymax=high_gev_post, 
        #              fill=utot90_gev), alpha=0.2) +
        geom_line(data=quant$spaghetti_qnt_df_gev_post,
                  aes(x=x,y=low_gev_post,color=utot90_gev), 
                  linetype="dashed",size=0.5)+
        geom_line(data=quant$spaghetti_qnt_df_gev_post, 
                  aes(x=x,y=high_gev_post),color='green',
                  linetype="dashed",size=0.5)+
        geom_line(data=quant$spaghetti_qnt_df_gev_post,
                  aes(x=x,y=MAP_gev_post,color=MAPcurvegev_stat),
                  linetype="solid",size=0.7)
      color.curves=c(color.curves,"MaxPost (GEV stat)"="green",
                                  "90% unc.Post (GEV stat)"="green")
      color.breaks= c(color.breaks, MAPcurvegev_stat, utot90_gev)
      color.linetype =c(color.linetype, "solid", "dashed")
      color.shape =c(color.shape, NA, NA)
      # color.Utot   = c(color.Utot,
      #                  "90% unc.Post (GEVstat)"    = "green")
    }
    # Empirical return levels:
    if (!is.null(s$empir_qnt_df[[id]])){
      quantileplot_smev=quantileplot_smev+
        geom_point(data=s$empir_qnt_df[[id]],  
                   aes(x=x,y=y,color=empir),shape=21,fill="black",size=2)
      color.curves=c(color.curves,"Empir. return levels"="black")
      color.breaks=c(color.breaks,empir)
      color.linetype=c(color.linetype,"blank")
      color.shape=c(color.shape,21)
    }
    
    if ((!is.null(s$row))|(!is.null(s$col))|(!is.null(s$sat_prod))){
      nfplot=paste0('/unc90_post_qnt_SMEV_',s$sat_prod,'_',s$row,'_',s$col,
                    '_', s$durations[id]/s$time_resolution, 'h.png')
      title_plot=paste0('Quantiles spaghetti SMEV -- duration=',
                        s$durations[id]/s$time_resolution, 'h \n (', s$sat_prod, ', ', 
                        year(s$time[1]),' - ',year( tail(s$time,1)),
                        ', pix:', s$row,'-', s$col,')')
    } else {
      nfplot=paste0('/unc90_post_qnt_SMEV_',s$name_project,
                    '_', s$durations[id]/s$time_resolution, 'h.png')
      title_plot=paste0('Quantiles SMEV -- duration=',
                        s$durations[id]/s$time_resolution, 'h \n (',
                        s$name_project, ', ', 
                        year(s$time[1]),' - ',year(tail(s$time,1)),')')
    }
    
    quantileplot_smev=quantileplot_smev+
      xlab("Return period T [years]") +
      ylab("Rainfall [mm/h]") +
      theme_bw(base_size=12)+
      geom_line()+
      # scale_y_continuous(expand=c(0,0),breaks=seq(5,40,5))+
      # scale_x_log10(breaks=c(1,2,5,10,20,50,100),expand=c(0,0),limits=c(1,100))+
      # scale_x_log10(breaks=c(1,2,5,10,50),expand=c(0,0),limits=c(1,50))+
      scale_x_log10(limits=c(1,100),expand=c(0,0)) +
      annotation_logticks(base=10, sides="b",scaled=T,
                          colour="black",size=0.3,linetype=1)+
      coord_cartesian(clip='off')+
      scale_fill_manual(name=element_blank(),
                        values=color.Utot)+
      scale_colour_manual(name=element_blank(),
                          values=color.curves,
                          breaks=color.breaks,
                          labels=color.breaks,
                          guide=guide_legend(override.aes=list(
                            linetype=color.linetype,
                            shape=color.shape)))+
      theme(panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,plot.margin=margin(c(0.5,1,0.5,1),unit="cm")
            ,legend.title=element_blank()
            ,legend.position=c(0.28,0.7)
            ,legend.text=element_text(size=10,
                                      margin=margin(r=0,l=2,t=0,b=0,unit="pt")) 
            ,plot.title=element_text(hjust=0.5)
             ,legend.key.size=grid::unit(0.8,'lines'))+
      ggtitle(title_plot)
    # save plot:
    ggsave(quantileplot_smev,
           filename=paste0(dir.res_pix,nfplot),                
           width=8,height=6,dpi=130)
    quantileplot_smev$duration=paste0(s$durations[id]/60, 'h')
  } else {
    # no plot:
    quantileplot_smev=NULL
    quantileplot_smev$duration=paste0(s$durations[id]/60,'h')
  }
  return(quantileplot_smev) 
}














################################################################################
plot_spaghetti_mev<- function(flg_save=T,
                              s,      # object containing all SMEV computations
                              quant,  # object containting MEV quantiles
                              id=1,
                              dir.res_pix){
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
  # SPAGHETTI PLOT
  if ( flg_save== T){
    # define text and colors for legend items:
    ##########################################
    Utot = "90% unc MEV (nonstat SMEV)"
    Utot95 = "95% unc MEV (nonstat SMEV)"
    
    spaghetti_MEV     = "Spaghetti MEV (nonstat SMEV)" 
    uncert90_MEV      = "90% unc MEV (nonstat SMEV)"
    Mediancurve_MEV   = "Median MEV (nonstat SMEV)"
    MAPcurve_nsSMEV   = "MaxPost nonstat SMEV"
    MAPcurve_statSMEV = "MaxPost stat SMEV"
    MAPcurve_gev      = "MaxPost stat GEV"
    empir             = "Empirical return levels"
    color.Utot        = c("90% unc MEV (nonstat SMEV)"= "blue")
                          #"95% unc MEV (nonstat SMEV)"= "blue")
    color.curves = c("Spaghetti MEV (nonstat SMEV)"  = "gray40",
                     "Median MEV (nonstat SMEV)"     = "red",
                     # "90% unc MEV (nonstat SMEV)"    = "blue",
                     # "MaxPost nonstat SMEV"              = "red",
                     # "MaxPost stat SMEV"                 = "orange",
                     # "MaxPost Stat GEV"                  = "green",
                     "Empirical return levels"           = "black")
    
    if ((!is.null(s$row))|(!is.null(s$col))|(!is.null(s$sat_prod))){
      nfplot=paste0('/spaghetti_post_qnt_MEV_',s$sat_prod,'_',s$row,'_',s$col,
                    '_', s$durations[id]/s$time_resolution, 'h.png')
      title_plot=paste0('Quantiles spaghetti MEV (from non-stationary SMEV) -- duration=',
                        s$durations[id]/s$time_resolution, 'h \n (', s$sat_prod, ', ', 
                        year(s$time[1]),' - ',year( tail(s$time,1)),
                        ', pix:', s$row,'-', s$col,')')
    } else {
      nfplot=paste0('/spaghetti_post_qnt_MEV_',s$name_project,
                    '_', s$durations[id]/s$time_resolution, 'h.png')
      title_plot=paste0('Quantiles spaghetti MEV (from non-stationary SMEV) -- duration=',
                        s$durations[id]/s$time_resolution, 'h \n (',
                        s$name_project, ', ', 
                        year(s$time[1]),' - ',year(tail(s$time,1)),')')
    }
    
    
    spaghettiplot_mev = ggplot() +
      geom_line(data=quant$spaghetti_mev_df_melt,  aes(x=grid,y=value, group = series, color=spaghetti_MEV), size=0.1)+
      geom_ribbon(data=quant$spaghetti_mev_qnt_df, aes(x=x, ymin=low, ymax=high, fill=Utot), alpha=0.1) +
      #geom_line(data=quant$spaghetti_mev_qnt_df,   aes(x=x, y =low,   color=uncert90_MEV), linetype="dashed", size=0.5)+
      #geom_line(data=quant$spaghetti_mev_qnt_df,   aes(x=x, y =high), color='blue', linetype="dashed", size=0.5)+
      geom_line(data=quant$spaghetti_mev_qnt_df,   aes(x=x, y =med,   color=Mediancurve_MEV), size=1.5)+
      # geom_line(data=quant$spaghetti_mev_qnt_df  aes(x=x, y=MAP,    color=MAPcurve_statSMEV), size=1.5)+
      # geom_line(data=quant$spaghetti_mev_qnt_df, aes(x=x, y=gev,    color=MAPcurve_gev), size=1.5)+
      geom_point(data=s$empir_qnt_df[[id]],        aes(x=x,y=y,       color=empir), shape=21, fill= "black", size=2)+
      xlab("Return period T [years]") +
      ylab("Rainfall [mm/h]") +
      theme_bw(base_size = 12)+
      geom_line()+
      scale_y_continuous(expand = c(0,0))+ #, breaks = seq(5,40,5))+
      #scale_x_log10(breaks = c(1,2,5,10,20,50,100), expand = c(0,0), limits= c(1,100)) +
      scale_x_log10(limits= c(1,100),  expand = c(0,0)) +
      # scale_x_log10(breaks = c(1,2,5,10,50), expand = c(0,0), limits= c(1,50)) +
      # scale_x_log10(limits= c(1,50),  expand = c(0,0)) +
      annotation_logticks(base = 10, sides = "b", scaled = TRUE,  colour = "black", size = 0.3, linetype = 1) +
      coord_cartesian(clip = 'off') +
      scale_fill_manual(name = element_blank(),
                        values = color.Utot) +
      scale_colour_manual(name     = element_blank(),
                          values   = color.curves,
                          breaks   = c(spaghetti_MEV, Mediancurve_MEV, empir),
                          labels   =c(spaghetti_MEV, Mediancurve_MEV, empir),
                          guide    = guide_legend(override.aes = list(
                            linetype = c("solid", "solid", "blank"),
                            shape    = c(NA, NA, 21)))) +
      theme(,panel.grid.major= element_blank()
            ,panel.grid.minor= element_blank()
            ,legend.title     = element_blank()
            ,legend.position  = c(0.3, 0.8)
            ,legend.text      = element_text(size=10, margin=margin(r=0, t=0, b=0, l=2, unit="pt"))
            ,plot.title       = element_text(hjust = 0.5)
            ,plot.margin      = margin(c(0.5, 1, 0.5, 1), unit = "cm")
            ,legend.key.size  = grid::unit(1.2, 'lines'))+
      ggtitle(title_plot)
    
    # save plot:
    ggsave(spaghettiplot_mev, 
           filename=paste0(dir.res_pix,nfplot),                
           width=8,height=6,dpi=130)
    
    spaghettiplot_mev$duration =  paste0(s$durations[id]/60, 'h')
    
  } else {
    
    # no smev inversion and no MEV plot
    spaghettiplot_mev = NULL
    spaghettiplot_mev$duration =  paste0(s$durations[id]/60, 'h')
    
  }
  
  return(spaghettiplot_mev)
}




    
    



################################################################################
plot_quantiles_mev<- function(flg_save=T,
                              s,      # object containing all SMEV computations
                              quant,  # object containting MEV quantiles
                              id=1,
                              dir.res_pix){
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
  
  # QUANTILES PLOT
  
  if ( flg_save== T){  
    
    # define text and colors for legend items:
    ##########################################
    Utot = "90% unc MEV (nonstat SMEV)"
    Utot95 = "95% unc MEV (nonstat SMEV)"
    
    Mediancurve_MEV = "Median MEV (nonstat SMEV)"
    MAPcurve_nsSMEV = "MaxPost nonstat SMEV"
    MAPcurve_statSMEV = "MaxPost stat SMEV"
    MAPcurve_gev = "MaxPost stat GEV"
    empir = "Empirical return levels"
    color.Utot   = c("90% unc MEV (nonstat SMEV)"= "blue",
                     "95% unc MEV (nonstat SMEV)"= "blue")
    color.curves = c("Median MEV (nonstat SMEV)"     = "red",
                     # "MaxPost nonstat SMEV"              = "red",
                     # "MaxPost stat SMEV"                 = "orange",
                     # "MaxPost Stat GEV"                  = "green",
                     "Empirical return levels"           = "black")
    
    if ((!is.null(s$row))|(!is.null(s$col))|(!is.null(s$sat_prod))){
      nfplot=paste0('/unc90_post_qnt_MEV_',s$sat_prod,'_',s$row,'_',s$col,
                    '_', s$durations[id]/s$time_resolution, 'h.png')
      title_plot=paste0('Quantiles MEV (from non-stationary SMEV) -- duration=',
                        s$durations[id]/s$time_resolution, 'h \n (', s$sat_prod, ', ', 
                        year(s$time[1]),' - ',year( tail(s$time,1)),
                        ', pix:', s$row,'-', s$col,')')
    } else {
      nfplot=paste0('/unc90_post_qnt_MEV_',s$name_project,
                    '_', s$durations[id]/s$time_resolution, 'h.png')
      title_plot=paste0('Quantiles MEV (from non-stationary SMEV) -- duration=',
                        s$durations[id]/s$time_resolution, 'h \n (',
                        s$name_project, ', ', 
                        year(s$time[1]),' - ',year(tail(s$time,1)),')')
    }
    
    quantileplot_mev = ggplot() +
      geom_ribbon(data=quant$spaghetti_mev_qnt_df, aes(x=x, ymin=low2.5, ymax=high97.5, fill=Utot95), alpha=0.1) +
      geom_ribbon(data=quant$spaghetti_mev_qnt_df, aes(x=x, ymin=low, ymax=high, fill=Utot), alpha=0.2) +
      geom_line(data=quant$spaghetti_mev_qnt_df,   aes(x=x, y=low),  color='blue', linetype="dashed", size=0.5)+
      geom_line(data=quant$spaghetti_mev_qnt_df,   aes(x=x, y=high), color='blue', linetype="dashed", size=0.5)+
      geom_line(data=quant$spaghetti_mev_qnt_df,   aes(x=x, y=med,   color=Mediancurve_MEV), size=1.5)+
      # geom_line(data=quant$spaghetti_mev_qnt_df  aes(x=x, y=MAP,   color=MAPcurve_statSMEV), size=1.5)+
      # geom_line(data=quant$spaghetti_mev_qnt_df, aes(x=x, y=gev,   color=MAPcurve_gev), size=1.5)+
      geom_point(data=s$empir_qnt_df[[id]],        aes(x=x,y=y,      color=empir), shape=21, fill= "black", size=2)+
      xlab("Return period T [years]") +
      ylab("Rainfall [mm/h]") +
      theme_bw(base_size=12)+
      geom_line()+
      scale_y_continuous(expand=c(0,0))+
      # scale_x_log10(breaks = c(1,2,5,10,20,50,100), expand = c(0,0), limits= c(1,100)) +
      # scale_x_log10(breaks = c(1,2,5,10, 50, 100), expand = c(0,0), limits= c(1,100)) +
      # scale_x_log10(limits= c(1,50),  expand = c(0,0)) +
      scale_x_log10(limits=c(1,100), expand=c(0,0)) +
      annotation_logticks(base=10, sides="b", scaled=T, colour="black", size=0.3, linetype=1) +
      coord_cartesian(clip='off')+
      scale_fill_manual(name   = element_blank(),
                        values = color.Utot) +
      scale_colour_manual(name     = element_blank(),
                          values   = color.curves,
                          breaks   = c(Mediancurve_MEV, empir),
                          labels   = c(Mediancurve_MEV, empir),
                          guide    = guide_legend(override.aes = list(
                            linetype=c("solid", "blank"),
                            shape=c(NA,21)))) +
      theme(    panel.grid.major=element_blank()
                ,panel.grid.minor=element_blank()
                ,plot.margin=margin(c(0.5,1,0.5,1),unit="cm")
                ,legend.title=element_blank()
                ,legend.position=c(0.3,0.8)
                ,legend.text=element_text(size=10,margin=margin(r=0,l=2,t=0,b=0,unit="pt")) 
                ,plot.title=element_text(hjust=0.5)
                ,legend.key.size=grid::unit(1.2,'lines'))+
      ggtitle(title_plot)
    # save plot:
    ggsave(quantileplot_mev,
           filename=paste0(dir.res_pix,nfplot),                
           width=8,height=6,dpi=130)
    quantileplot_mev$duration=paste0(s$durations[id]/60,'h')
  } else {
    # no smev inversion and no MEV plot
    quantileplot_mev=NULL
    quantileplot_mev$duration=paste0(s$durations[id]/60,'h')
  }
  return(quantileplot_mev)
}















