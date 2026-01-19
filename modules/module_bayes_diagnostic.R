library(cmdstanr) # STAN
library(posterior) # read STAN results
library(bayesplot) # for mcmc traceplots
library(loo) # for ELPD computation 
library(ggplot2) # for plots
library(PEIP) # for chi-square test
library(coda) # for mcmc traceplots
library(reshape2) # with function melt() for dataframes
library(psych) # with pairs.panels()  



################################################################################
likelihood_ratio_test = function(loglik_H0, loglik_H1, degrees_of_freedom){
################################################################################
  #^* GOAL: Computes the statistical significance of model H1 with respect to
  #^*       model H0 depending on the number of additional parameters
  #^*       (degrees_of_freedom)
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^*             Francesco Marra, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [real] loglik_H0 --> log-likelihood of model H0
  #^*    2. [real] loglik_H1 --> log-likelihood of model H1
  #^*    3. [integer] degrees_of_freedom --> number of additional params
  #^* OUT
  #^*    1. [real] pval
  #^****************************************************************************
  #^* REF.: uses PEIP package for chi square 
  #^*       https://www.rdocumentation.org/packages/PEIP/versions/2.2-5
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  #^****************************************************************************  

  # likelihood ratio = -2*log(Lik_H0/Lik_H1)
  pval = 1 - PEIP::chi2cdf(-2*( loglik_H0 - loglik_H1 ), 
                           degrees_of_freedom)
  
  return(pval)
}





################################################################################
mcmc.plots= function(zPost, dir_results){
################################################################################
  #^* GOAL: plot the mcmc pdfs
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [integer] zPost --> 
  #^*    2. [character] dir_results --> 
  #^* OUT
  #^*    1.
  #^****************************************************************************
  #^* REF.: 
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS: 
  #^****************************************************************************
  for (par in 1:ncol(zPost)){
    pdf_param=ggplot()+
      geom_density(aes(x=zPost[,par]),color="blue",fill="blue",alpha=0.3,lwd=1.2)+
      geom_density(aes(x=zPrior[,par]),color="red",fill="red",alpha=0.3,lwd=1.2)+
      theme_bw(base_size=20)
    ggsave(pdf_param,
           filename=paste0(dir_results,"/pdf_param_",par,"_post_vs_prior_yy_1.png"))
  }
}





################################################################################
calc.dic <- function(zPost, lik, post, lik.fun,...){
################################################################################
  #^* GOAL: Estimate Deviance Information Criterion (DIC)
  #^*       DIC is based on the deviance D(theta) = (-2*log(p(x|theta))) and a
  #^*       penalization using the effective number of parameters (pD).
  #^*       ==> DIC = Dthetabar + 2*pD
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [array of reals] zPost --> matrix of posterior samples
  #^*    2. [array of reals] lik --> vector of likelihood of posterior samples
  #^*    2. [array of reals] post --> vector of the log-posterior
  #^*    3. [function] lik.fun --> function that compute the likelihood
  #^*    4. ...   --> other parameters or variables that are passed to 'lik.fun'
  #^* OUT
  #^*    list()
  #^*    1. DIC --> Deviance Information Criterion
  #^*    2. IC --> Bayesian Predictive Information Criterion
  #^*    3. pD --> Effective number of parameters (pD = Dbar - Dhat)
  #^*    4. pV --> Effective number of parameters (pV = var(D)/2)
  #^*    5. Dbar --> Expected value of the deviance over the posterior
  #^*    6. Dhat --> Deviance at the mean posterior estimate
  #^****************************************************************************
  #^* REF.: Bayesian Data Analysis.
  #^*       Gelman, A., Carlin, J., Stern, H., and Rubin D. Second Edition, 2003.
  #^*       Bayesian predictive information criterion for the evaluation of 
  #^*       hierarchical Bayesian and empirical Bayes models. 
  #^*       Biometrika, 2007
  #^*       Benjamin Renard (BaM, https://github.com/BaM-tools/BaM/issues/5)
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  #^*       thetabar could be the maxpost or the mean/median, it is a choice! 
  #^*       (Spiegelhalter et al. 2002).
  #^*       Three different versions of DIC found in literature
  #^*       (https://github.com/BaM-tools/BaM/issues/5):
  #^*
  #^*       1) DIC1 =D(thetahat)+ 2*[E(D(theta)) - D(thetahat))] 
  #^*               =2*E(D(theta))- D(thetahat)   
  #^*              A) -4*mean(loglike)+ 2*loglike(mean(thetaposterior))
  #^*              B) -4*mean(loglike)+ 2*loglike(thetamaxpost)
  #^*       2) DIC2 =D(thetahat) + 2*(1/2* var(D(theta)))
  #^*              A) -2*loglike(mean(thetaposterior))+2*(1/2*var(-2*loglike))
  #^*              B) -2*loglike(thetamaxpost) + 2 *(1/2*var(-2*loglike))
  #^*       3) DIC3 =E(D(theta) +  1/2*var(D(theta)) 
  #^*                =-2*mean(log(p(x|theta))) + 1/2*var(-2*log(p(x|theta))) 
  #^*               (it's the average between DIC1 and DIC2)
  #^*              A) -2*mean(loglike) + 1/2 *var(-2*loglike))
  #^****************************************************************************
  # Considering the posterior LIKELIHOOD:
  # using posterior mean:
  D.bar <- -2*mean(lik)  # mean likelihood
  theta.bar <- apply(zPost,2,mean) 
  # deviance of likelihood at mean set of parameters
  D.hat <- -2*lik.fun(theta.bar,...)  
  pD <- D.bar-D.hat
  # half variance of likelihood
  pV <- var(-2*lik)/2  
  DIC1=pD+D.bar
  DIC2=D.hat+2*pV
  DIC3=D.bar+pV
  D.bar <- -2*mean(lik)  # mean likelihood
  # using MAP (maximum a posteriori):
  theta.bar <- zPost[which.max(post),]
  # deviance of likelihood at MAP set of parameters
  D.hat <- -2*lik.fun(theta.bar,...)  
  pD <- D.bar-D.hat
  # half variance of likelihood
  pV <- var(-2*lik)/2  
  DIC1_b=pD+D.bar
  DIC2_b=D.hat+2*pV
  return(list( DIC1   = DIC1,
               DIC1_b = DIC1_b,
               DIC2   = DIC2,
               DIC2_b = DIC2_b,
               DIC3   = DIC3,
               IC     = 2*pD+D.bar,
               pD     = pD,
               pV     = pV,
               Dbar   = D.bar,
               Dhat   = D.hat))
}







################################################################################
calc.bic <- function(zPost, lik, post, y) {
################################################################################
  #^* GOAL: Bayesian Information Criterion (BIC)
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [array of reals] zPost --> matrix of posterior samples (k=num.param)
  #^*    2. [array of reals] lik --> vector of likelihood of posterior samples
  #^*    3. [array of reals] post --> vector of the log-posterior
  #^*    4. [array of reals] y --> vector of observations (and n=length)
  #^* OUT
  #^*    list()
  #^*    1. BIC --> Bayesian Information Criterion = k*ln(n) - 2*ln(max(L))
  #^*    2. BIC_b --> Bayesian Information Criterion = k*ln(n) - 2*ln(max(Post))
  #^****************************************************************************
  #^* REF.: Schwarz, G. (1978) Estimating the Dimension of a Model, 
  #^*       Annals of Statistics 6, 461--464.
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  #^****************************************************************************
  penaltyterm=log(length(y))*ncol(zPost)
  likeliterm=-2*max(lik)
  posteriorterm=-2*max(post)
  return(list(BIC  =penaltyterm+likeliterm,
              BIC_b=penaltyterm+posteriorterm))
}







################################################################################
calc.aic <- function(zPost, lik, post){
################################################################################
  #^* GOAL: Akaike Information Criterion (AIC)
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [array of reals] zPost --> matrix of posterior samples (k=num.param)
  #^*    2. [array of reals] lik --> vector of likelihood of posterior samples
  #^*    3. [array of reals] post --> vector of the log-posterior
  #^* OUT
  #^*    list()
  #^*    1. AIC --> Akaike Information Criterion = 2*k - 2*ln(max(L))
  #^*    2. AIC_b --> Akaike Information Criterion = 2*k - 2*ln(max(Post))
  #^****************************************************************************
  #^* REF.: Akaike, H. (1973),"Information theory and an extension of the maximum 
  #^*       likelihood principle", 
  #^*       in Petrov,B.N.; Csáki,F., 2nd International Symposium on Information 
  #^*       Theory,1971, Budapest: Akadémiai Kiadó, pp. 267–281.
  #^*       Springer, DOI: 10.1007/978-1-4612-0919-5_38
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  #^****************************************************************************
  penaltyterm=2*ncol(zPost)
  likeliterm=-2*max(lik)
  posteriorterm=-2*max(post)
  return(list(AIC  =penaltyterm+likeliterm,
              AIC_b=penaltyterm+posteriorterm))
}





################################################################################
calc.hqc <- function(zPost, lik, post, y) {
################################################################################
  #^* GOAL: Hannan-Quinn information Criterion (HQC)
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [array of reals] zPost --> matrix of posterior samples (k=num.param)
  #^*    2. [array of reals] lik --> vector of likelihood of posterior samples
  #^*    3. [array of reals] post --> vector of the log-posterior
  #^*    4. [array of reals] y --> vector of observations (and n=length)
  #^* OUT
  #^*    list()
  #^*    1. HQC -->  = -2*L_{max} + 2*k*ln(ln(n))
  #^*    2. HQC_b -->  = -2*L_{Post} + 2*k*ln(ln(n))
  #^****************************************************************************
  #^* REF.: Hannan, E. J.; Quinn, B. G. (1979-01-01). 
  #^*       "The Determination of the Order of an Autoregression". 
  #^*       Journal of the Royal Statistical Society: Series B (Methodological). 
  #^*       41(2):190–195. doi:10.1111/j.2517-6161.1979.tb01072.x. ISSN0035-9246.
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  #^****************************************************************************
  penaltyterm=log(log(length(y)))*2*ncol(zPost)
  likeliterm=-2*max(lik)
  posteriorterm=-2*max(post)
  return(list(HQC  =penaltyterm+likeliterm,
              HQC_b=penaltyterm+posteriorterm))
}







   
################################################################################
cmdstan.diagn = function(stanfit, dir_res, name_output){
################################################################################
  #^* GOAL: diagnostics of mcmc samples from cmdstan (STAN) alogorithm
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [stan object] stanfit --> object of STAN results
  #^*    2. [character] dir_res --> path of folder where to save results
  #^*    3. [character] name_output --> strng to differentiate results files
  #^* OUT
  #^*    list()
  #^****************************************************************************
  #^* REF.: cmdstan, bayesplot, loo packages
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  #^****************************************************************************
  # summary statistics of mcmc chains from stan object:
  # stanfit$summary(variables=c("alpha_0","alpha_1","sigma_0","sigma_1"),"mean","sd")
  # stanfit$diagnostic_summary()
  
  # plot traceplots of mcmc chains:
  draws<-stanfit$draws()
  list_of_draws=as_draws_df(draws)
  color_scheme_set("mix-blue-red")
  pars_trace=bayesplot::mcmc_trace(
                draws,
                pars=c("alpha_0","alpha_1","sigma_0","sigma_1"),
                n_warmup=0,
                facet_args=list(nrow=4,labeller=label_parsed,ncol=1,
                                strip.position="left"),
                size=1)+
                facet_text(size=15) 
  ggsave(pars_trace,
         filename=paste0(dir_res,"/traceplots_mcmc_",name_output,".png"),
         width=6,height=4,dpi=120)
  
  # Criteria for expected predictive accuracy:
  # PSIS-LOO computation, estimates of:
  # --> "elpd_loo" (expected log predictive density)
  # --> p_loo (effective number of parameters)
  # --> looic=−2*elpd_loo (LOO criterion, conventional scale of "deviance")
  ELPD=loo::loo(stanfit$draws("log_lik"),  
                r_eff=relative_eff(stanfit$draws("log_lik")),
                merge_chains=T,
                save_psis=T)
  # WAIC information criterion:
  WAIC=loo::waic(stanfit$draws("log_lik"))
}
  









################################################################################
mcmc.diagn = function(zPost, 
                      priors_id,
                      dir_res, 
                      name_output,
                      save=T){
################################################################################
  #^* GOAL: diagnostics of mcmc samples (study of convergence)
  #^*       NON_STATIONARY SMEV (shape int, shape slope, scale int, scale slope)
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [array of reals] zPost --> matrix of posterior samples
  #^*    2. [list of lists] priors_id --> priors list
  #^*    3. [character] dir_res --> path of folder where to save results
  #^*    4. [character] name_output --> string to differentiate results files
  #^*    5. [logical] save --> flag for saving plots
  #^* OUT
  #^*    1. [logical] convergence --> TRUE if good MCMC convergence
  #^****************************************************************************
  #^* REF.: uses coda, psych, stats, reshape2 packages
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  ##############################################################################

  #^****************************************************************************
  #^                        study of convergence
  #^****************************************************************************
  len=length(zPost[,1])
  gel=NULL
  convergence=T
  nam_param=c("shape intercept","shape slope","scale intercept","scale slope")
  for (i in 1:4){
    # Gelman test (divide mcmc in 4 chains):
    mcmc.list=list(zPost[1:(len/4),i], 
                   zPost[(len/4+1):(len/2),i], 
                   zPost[(len/2+1):(len*3/4),i], 
                   zPost[(len*3/4+1):len,i])
    mcmc.list=as.mcmc.list(lapply(mcmc.list, coda::mcmc))
    gel[[i]]=coda::gelman.diag(x=mcmc.list,confidence=0.95,transform=F)
    par(mar=c(6,6,2,5))  # Set the margin on all sides 
    coda::gelman.plot(x=mcmc.list,auto.layout=F) 
    title(main=names(zPost)[i],cex.main=2)
    #acfplot(x = mcmc.list, lag.max=200)
    if (gel[[i]]$psrf[1]>1.2){ # Gelman factor should stay > 1.2 
      # GELMAN FACTOR: limit is suggested around 1.2
      convergence = F
    }
    dev.off()
  }
  
  #^****************************************************************************
  #^                             Plots (if save=T)
  #^****************************************************************************
  if (save){
  ##########
    # Plot traceplots:
    png(filename=paste0(dir_res,"/traceplots_mcmc_",name_output,".png"),
        width=1000,height=1000)
    par(mfrow=c(4,1))
    coda::traceplot(as.mcmc(zPost[,1]),main=nam_param[1],smooth=T,cex.lab=2,
                    cex.axis=1,cex.main=2,cex.sub=1.5)
    coda::traceplot(as.mcmc(zPost[,2]),main=nam_param[2],smooth=T,cex.lab=2,
                    cex.axis=1,cex.main=2,cex.sub=1.5)
    coda::traceplot(as.mcmc(zPost[,3]),main=nam_param[3],smooth=T,cex.lab=2,
                    cex.axis=1,cex.main=2,cex.sub=1.5)
    coda::traceplot(as.mcmc(zPost[,4]),main=nam_param[4],smooth=T,cex.lab=2,
                    cex.axis=1,cex.main=2,cex.sub=1.5)
    dev.off()
  
    # mcmc autocorrelation:
    png(paste0(dir_res,"/autocorr_mcmc_",name_output,".png"),
        width=1000,height=1000)
    par(mfrow=c(4,1))
    for (i in 1:(ncol(zPost))){
      par(mar=c(6,6,2,5))  # Set the margin on all sides 
      coda::autocorr.plot(x=zPost[,i],auto.layout=F,cex.lab=3,
                    cex.axis=3,cex.main=3,cex.sub=2)
      title(main=names(zPost)[i],cex.main=2)
    }
    dev.off()
  
    # Gelman test (divide mcmc in 4 chains):
    png(paste0(dir_res,"/gelman_mcmc_",name_output,".png"),
        width=1000,height=1000)
    par(mfrow=c(4,1))
    dev.off()
    

    # Compute and plot Gelman Reduction Factor:
    Rc="Reduction factor Rc"
    Ru="Upper confidence limit Ru" 
    limit="threshold"
    col.gelman=c("Reduction factor Rc"      ="blue", 
                 "Upper confidence limit Ru"="red", 
                 "threshold"                ="red")
    Rc.Ru.plot<-ggplot()+
      coord_cartesian(xlim=c(1,length(gel)))+
      #scale_y_continuous(breaks=seq(0,Qmax,Qmax/5))+
      scale_x_continuous(breaks=seq(1,length(gel),1)) +
      geom_point(mapping=aes(x=seq(1,length(gel),1),
                             unlist(gel,use.names=F)[seq(1,length(gel)*2,2)], 
                             col=Rc),shape=0,size=4,stroke=0.5)+
      geom_point(mapping=aes(x=seq(1,length(gel),1),
                             unlist(gel,use.names=F)[seq(2,length(gel)*2,2)], 
                             col=Ru),shape=4,size=4,stroke=0.5)+
      theme_bw(base_size=12)+
      geom_hline(aes(yintercept=1.2,col=limit),lwd=1,linetype="dashed")+
      scale_colour_manual(name="Legend", 
                          values=col.gelman ,
                          breaks=c(Rc,Ru,limit),
                          guide=guide_legend(override.aes=list(
                            linetype=c("blank","blank","dashed"),
                            shape=c(0,4,NA)), 
                            title.position="top"))+
      xlab("Parameter index")+ylab("Gelman factor")+
      theme(text=element_text(size=12)
            ,panel.grid.major=element_line(size=0.4,linetype="dashed")
            ,panel.grid.minor=element_blank()
            ,legend.position="bottom"
            ,legend.text=element_text(size=10)
            ,legend.justification=c(0,1)
            ,legend.box="vertical"
            ,legend.box.just="left"
            ,legend.text.align=0)
    ggsave(Rc.Ru.plot,
           filename=paste0(dir_res,"/Rc_Ru_gelman_mcmc_",name_output,".png"),
           width=6,height=4,dpi=130)
  
    # Density plots mcmc:
    df=as.data.frame(zPost)
    names(df)=nam_param
    if (!is.null(priors_id)){
      # priors pdf:
      dPrior=matrix(NA,nrow=10000,ncol=4)
      grid=matrix(NA,10000,ncol(dPrior)) 
      for(i in 1:ncol(dPrior)){
        grid[,i]=seq(min(df[,i]),max(df[,i]),length.out=10000)
      }    
      # compute pdf of specified priors:
      # shape intercept:
      if (priors_id$shape_intercept[[1]]=='uniform'){
        dPrior[,1]=dunif(x=grid[,1],
                         min=priors_id$shape_intercept[[2]], 
                         max=priors_id$shape_intercept[[3]])
      } else if (priors_id$shape_intercept[[1]]=='normal'){
        dPrior[,1]=dnorm(x=grid[,1],
                         mean=priors_id$shape_intercept[[2]],
                         sd=priors_id$shape_intercept[[3]])
      } else if (priors_id$shape_intercept[[1]]=='lognormal'){
        dPrior[,1]=dlnorm(x=grid[,1], 
                          meanlog=priors_id$shape_intercept[[2]], 
                          sdlog=priors_id$shape_intercept[[3]])
      } else if (priors_id$shape_intercept[[1]]=='flatprior+'){
        dPrior[,1] = rep(0, grid[,1])
      } else {
        stop("prior distribution not available for shape intercept !!!")
      }
      # shape slope:
      if (priors_id$shape_slope[[1]]=='uniform'){
        dPrior[,2]=dunif(x=grid[,2],
                         min=priors_id$shape_slope[[2]],
                         max=priors_id$shape_slope[[3]])
      } else if (priors_id$shape_slope[[1]]=='normal'){
        dPrior[,2]=dnorm(x=grid[,2],
                         mean=priors_id$shape_slope[[2]],
                         sd=priors_id$shape_slope[[3]])
      } else if (priors_id$shape_slope[[1]]=='lognormal'){
        dPrior[,2]=dlnorm(x=grid[,2], 
                          meanlog=round(Transf_Gauss_lognorm(
                            E=priors_id$shape_slope[[2]],
                            stdev=priors_id$shape_slope[[3]])$mu,digits=4),
                          sdlog=round(Transf_Gauss_lognorm(
                            E=priors_id$shape_slope[[2]],
                            stdev=priors_id$shape_slope[[3]])$sd,digits=4))
      } else if (priors_id$shape_slope[[1]]=='flatprior+'){
        dPrior[,2]=rep(0,grid[,2])
      } else {
        stop("prior distribution not available for shape slope !!!")
      }
      # scale intercept:
      if (priors_id$scale_intercept[[1]]=='uniform'){
        dPrior[,3]=dunif(x=grid[,3],
                         min=priors_id$scale_intercept[[2]],
                         max=priors_id$scale_intercept[[3]])
      } else if (priors_id$scale_intercept[[1]]=='normal'){
        dPrior[,3]=dnorm(x=grid[,3],
                         mean=priors_id$scale_intercept[[2]], 
                         sd=priors_id$scale_intercept[[3]])
      } else if (priors_id$scale_intercept[[1]]=='lognormal'){
        dPrior[,3]=dlnorm(x=grid[,3], 
                          meanlog=priors_id$scale_intercept[[2]], 
                          sdlog=priors_id$scale_intercept[[3]])
      } else if (priors_id$scale_intercept[[1]]=='flatprior+'){
        dPrior[,3]=rep(0,grid[,3])
      } else {
        stop("prior distribution not available for scale intercept !!!")
      }
      # scale slope:
      if (priors_id$scale_slope[[1]]=='uniform'){
        dPrior[,4]=dunif(x=grid[,4],
                         min=priors_id$scale_slope[[2]],
                         max=priors_id$scale_slope[[3]])
      } else if (priors_id$scale_slope[[1]]=='normal'){
        dPrior[,4]=dnorm(x=grid[,4], 
                         mean=priors_id$scale_slope[[2]], 
                         sd=priors_id$scale_slope[[3]])
      } else if (priors_id$scale_slope[[1]]=='lognormal'){
        dPrior[,4]=dlnorm(x=grid[,4], 
                          meanlog=round(Transf_Gauss_lognorm(
                            E=priors_id$scale_slope[[2]], 
                            stdev=priors_id$scale_slope[[3]])$mu,digits=4),
                          sdlog=round(Transf_Gauss_lognorm(
                            E=priors_id$scale_slope[[2]], 
                            stdev=priors_id$scale_slope[[3]])$sd,digits=4))
      } else if (priors_id$scale_slope[[1]]=='flatprior+'){
        dPrior[,4]=rep(0,grid[,4])
      } else {
        stop("prior distribution not available for scale slope !!!")
      }
      # create unique dataframe of priors:
      dfprior=as.data.frame(dPrior)
      names(dfprior)=nam_param
      priorDF=dfprior
      priorDF.grid=data.frame(grid)
      names(priorDF.grid)=nam_param
      priorDF.grid.bis=cbind(priorDF.grid, indx=1:nrow(priorDF.grid))
      priorDF.grid.bis.m=melt(priorDF.grid.bis, id.vars='indx')
      priorDF.bis=cbind(priorDF, indx=1:nrow(priorDF))
      priorDF.bis.m=cbind(melt(priorDF.bis, id.vars='indx'),
                          grid=priorDF.grid.bis.m$value)
    }
    X=cbind(df, indx=1:nrow(df))
    Xm=melt(X, id.vars='indx')
    
    # plot pdf:
    pdf_param=ggplot(Xm)+
      # geom_density(aes(x= zPost),color="blue",fill="blue",alpha=0.3,lwd=1.2)+
      geom_density(aes(value,fill=variable),colour=NA,alpha=0.8)
    if (!is.null(priors_id)){
      pdf_param=pdf_param+
      geom_line(data=priorDF.bis.m,aes(x=grid,y=value),colour="black")
    }
    pdf_param=pdf_param+
      theme_bw(base_size=15)+
      facet_wrap(~variable, scales='free',ncol=2)+
      theme(legend.position="none"
            ,axis.text=element_text(size=10)
            ,axis.title.x=element_blank()
            ,panel.grid.minor=element_blank())
    ggsave(pdf_param,
           filename=paste0(dir_res,"/pdf_param_mcmc_",name_output,".png"),
           width=6,height=6,dpi=100)
  
    # Plot scatterplots:
    # m=GGally::ggpairs(as.data.frame(zPost[[iter]]),
    #                   diag=list(continuous='density'),
    #                   lower=list(continuous='points'),
    #                   upper=list(continuous='cor'),
    #                   axisLabels='show')
    # png(filename=paste0(dir_res,"/scatterplots_mcmc_",name_output,".png"),
    #     width=1000,height=1000)
    # pairs(zPost)
    # dev.off()

    # Plot correlations between parameters
    png(filename=paste0(dir_res, "/corr_matrix_mcmc_",name_output,".png"),
        width=10,height=7,units='in',res=100)
    plt=psych::pairs.panels(df,
                            method="pearson", # correlation method
                            hist.col="#00AFBB", # color 
                            density=T,  # show density plots
                            ellipses=T) # show correlation ellipses
    dev.off()
  }
  return(convergence)
}











################################################################################
mcmc.diagn_stat=function(zPost,
                         priors_id,
                         dir_res,
                         name_output,
                         save=T){
################################################################################
  #^* GOAL: diagnostics of mcmc samples (study of convergence)
  #^*       STATIONARY SMEV (shape int, scale int)
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [array of reals] zPost --> matrix of posterior samples
  #^*    2. [list of lists] priors_id --> priors list
  #^*    3. [character] dir_res --> path of folder where to save results
  #^*    4. [character] name_output --> string to differentiate results files
  #^*    5. [logical] save --> flag for saving plots
  #^* OUT
  #^*    1. [logical] convergence --> TRUE if good MCMC convergence
  #^****************************************************************************
  #^* REF.: uses coda, psych, stats, reshape2 packages
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  ##############################################################################
  
  #^****************************************************************************
  #^                        study of convergence
  #^****************************************************************************
  len=length(zPost[,1])
  gel=NULL
  convergence=T
  nam_param=c("shape","scale")
  # Gelman test (divide mcmc in 4 chains):
  for (i in 1:2){
    mcmc.list=list(zPost[1:(len/4),i], 
                   zPost[(len/4+1):(len/2),i], 
                   zPost[(len/2+1):(len*3/4),i], 
                   zPost[(len*3/4+1):len,i])
    mcmc.list=as.mcmc.list(lapply(mcmc.list, coda::mcmc))
    gel[[i]]=coda::gelman.diag(x=mcmc.list,confidence=0.95,transform=F)
    par(mar=c(6,6,2,5))   # Set the margin on all sides 
    coda::gelman.plot(x=mcmc.list,auto.layout=F,cex.lab=2,cex.axis=2,
                      cex.main=2,cex.sub=2) 
    title(main=names(zPost)[i],cex.main=2)
    #acfplot(x=mcmc.list,lag.max=200)
    if (gel[[i]]$psrf[1] > 1.2) { 
      # GELMAN FACTOR: limit is suggested around 1.2
      convergence=F  # ===> bad convergence !!!
    }
  }  
  dev.off()
  
  #^****************************************************************************
  #^                         Plots (if save=T)
  #^****************************************************************************
  if (save){
    # Plot traceplots:
    png(filename=paste0(dir_res, "/traceplots_mcmc_",name_output,".png"), 
        width=1000,height=500)
    par(mfrow=c(2,1))
    coda::traceplot(as.mcmc(zPost[,1]),main="shape",smooth=T,cex.lab=2,
                    cex.axis=1,cex.main=2,cex.sub=1.5)
    coda::traceplot(as.mcmc(zPost[,2]),main="scale",smooth=T,cex.lab=2,
                    cex.axis=1,cex.main=2,cex.sub=1.5)
    dev.off()
    
    # mcmc autocorrelation:
    png(paste0(dir_res,"/autocorr_mcmc_",name_output,".png"),
        width=1000,height=500)
    par(mfrow=c(2,1))
    for (i in 1:(ncol(zPost))){
      par(mar=c(6,6,2,5))  # Set the margin on all sides 
      coda::autocorr.plot(x=zPost[,i],auto.layout=F,cex.lab=2,cex.axis=2,
                          cex.main=2,cex.sub=2)
      # title(main=names(zPost)[i],cex.main=2)
    }
    dev.off()
  
    # Gelman test (divide mcmc in 4 chains):
    png(paste0(dir_res,"/gelman_mcmc_",name_output,".png"),
        width=1000,height=500)
    par(mfrow=c(2,1))
    dev.off()

    # Compute and plot Gelman Reduction Factor:
    Rc= "Reduction factor Rc"
    Ru="Upper confidence limit Ru" 
    limit= "threshold"
    col.gelman = c("Reduction factor Rc"       = "blue", 
                   "Upper confidence limit Ru" = "red", 
                   "threshold"                 = "red")
    Rc.Ru.plot <-ggplot()+
      coord_cartesian(xlim=c(1,length(gel)))+
      scale_x_continuous(breaks=seq(1,length(gel),1)) +
      geom_point(mapping=aes(x=seq(1,length(gel),1),
                             unlist(gel,use.names=F)[seq(1,length(gel)*2,2)], 
                             col=Rc),shape=0,size=4,stroke=0.5) +
      geom_point(mapping=aes(x=seq(1,length(gel),1),
                             unlist(gel,use.names=F)[seq(2,length(gel)*2,2)], 
                             col=Ru),shape=4,size=4,stroke=0.5) +
      theme_bw(base_size=15)+
      geom_hline(aes(yintercept=1.2,col=limit),lwd=1,linetype="dashed")+
      scale_colour_manual(name="Legend",
                          values=col.gelman,
                          breaks=c(Rc,Ru,limit),
                          guide=guide_legend(override.aes=list(
                            linetype=c("blank","blank","dashed"),
                            shape=c(0,4,NA)), 
                            title.position="top"))+
      xlab("Parameter index")+ 
      ylab("Gelman factor")+
      theme(text=element_text(size=15)
            ,panel.grid.major=element_line(size=0.4,linetype="dashed")
            ,panel.grid.minor=element_blank()
            ,legend.position="bottom"
            ,legend.text=element_text(size=10)
            ,legend.justification=c(0,1)
            ,legend.box="vertical"
            ,legend.box.just="left"
            ,legend.text.align=0)
    ggsave(Rc.Ru.plot, 
           filename=paste0(dir_res,"/Rc_Ru_gelman_mcmc_",name_output,".png"),
           width=6,height=4,dpi=130)
    
    # Density plots mcmc:
    df=as.data.frame(zPost)
    names(df)=nam_param
    if (!is.null(priors_id)){
      # priors pdf:
      dPrior=matrix(NA,nrow=10000,ncol=2)
      grid=matrix(NA,10000,ncol(dPrior)) 
      for(i in 1:ncol(dPrior)){
        grid[,i]=seq(min(df[,i]),max(df[,i]),length.out=10000)
      }
      # compute pdf of specified priors:
      # shape intercept:
      if (priors_id$shape_intercept[[1]]=='uniform'){
        dPrior[,1]=dunif(x=grid[,1], 
                         min=priors_id$shape_intercept[[2]], 
                         max=priors_id$shape_intercept[[3]])
      } else if (priors_id$shape_intercept[[1]]=='normal'){
        dPrior[,1]=dnorm(x=grid[,1],
                         mean=priors_id$shape_intercept[[2]],
                         sd=priors_id$shape_intercept[[3]])
      } else if (priors_id$shape_intercept[[1]]=='lognormal'){
        dPrior[,1]=dlnorm(x=grid[,1], 
                          meanlog=round(Transf_Gauss_lognorm(
                            E=priors_id$shape_intercept[[2]],
                            stdev=priors_id$shape_intercept[[3]])$mu,digits=4),
                          sdlog=round(Transf_Gauss_lognorm(
                            E=priors_id$shape_intercept[[2]], 
                            stdev=priors_id$shape_intercept[[3]])$sd,digits=4))
      } else if (priors_id$shape_intercept[[1]]=='flatprior+'){
        dPrior[,1]=rep(0,grid[,1])
      } else {
        stop("prior distribution not available for shape intercept !!!")
      }
      # scale intercept:
      if (priors_id$scale_intercept[[1]]=='uniform'){
        dPrior[,2]=dunif(x=grid[,2],
                         min=priors_id$scale_intercept[[2]], 
                         max=priors_id$scale_intercept[[3]])
      } else if (priors_id$scale_intercept[[1]]=='normal'){
        dPrior[,2]=dnorm(x=grid[,2], 
                           mean=priors_id$scale_intercept[[2]],
                           sd=priors_id$scale_intercept[[3]])
      } else if (priors_id$scale_intercept[[1]]=='lognormal'){
        dPrior[,2]=dlnorm(x=grid[,2], 
                          meanlog=round(Transf_Gauss_lognorm(
                            E=priors_id$scale_intercept[[2]], 
                            stdev=priors_id$scale_intercept[[3]])$mu,digits=4),
                          sdlog=round(Transf_Gauss_lognorm(
                            E=priors_id$scale_intercept[[2]],
                            stdev=priors_id$scale_intercept[[3]])$sd,digits=4))
      } else if (priors_id$scale_intercept[[1]]=='flatprior+'){
        dPrior[,2]=rep(0, grid[,2])
      } else {
        stop("prior distribution not available for scale intercept !!!")
      }
      # create unique dataframe of priors:
      dfprior=as.data.frame(dPrior)
      names(dfprior)=nam_param
      priorDF=dfprior
      priorDF.grid=data.frame(grid)
      names(priorDF.grid)=nam_param
      priorDF.grid.bis=cbind(priorDF.grid,indx=1:nrow(priorDF.grid))
      priorDF.grid.bis.m=melt(priorDF.grid.bis, id.vars='indx')
      priorDF.bis=cbind(priorDF, indx=1:nrow(priorDF))
      priorDF.bis.m=cbind(melt(priorDF.bis, id.vars='indx'), 
                          grid=priorDF.grid.bis.m$value)
    }
    X=cbind(df,indx=1:nrow(df))
    Xm=melt(X,id.vars='indx')
    # plot pdf:
    pdf_param=ggplot(Xm) +
      geom_density(aes(value,fill=variable),colour=NA,alpha=0.8)
    if (!is.null(priors_id)){
      pdf_param=pdf_param +
        geom_line(data=priorDF.bis.m,aes(x=grid,y=value),colour="black")
    }
    pdf_param=pdf_param+
      theme_bw(base_size=15)+
      facet_wrap(~variable,scales='free',ncol=2)+
      theme(legend.position="none"
            ,axis.text=element_text(size=10)
            ,axis.title.x=element_blank()
            ,panel.grid.minor=element_blank())
    ggsave(pdf_param,
           filename=paste0(dir_res,"/pdf_param_mcmc_",name_output,".png"),
           width=6,height=6,dpi=100)

    # Correlations between parameters
    png(filename=paste0(dir_res,"/corr_matrix_mcmc_",name_output,".png"),
        width=10,height=7,units='in',res=100)
    plt=psych::pairs.panels(df,
                     method="pearson", # correlation method
                     hist.col="#00AFBB", # color
                     density=T,  # show density plots
                     ellipses=T) # show correlation ellipses
    dev.off()
    # Plot scatterplots:
    # png(filename=paste0(dir_res,"/scatterplots_mcmc_",name_output,".png"), 
    #     width=1000, height=1000)
    # pairs(zPost)
    # dev.off()
  }
  return(convergence)
}









################################################################################
mcmc.diagn_partns = function(zPost, 
                             priors_id,
                             dir_res, 
                             name_output,
                             save=T){
################################################################################
  #^* GOAL: diagnostics of mcmc samples (study of convergence)
  #^*       PARTIAL NON_STATIONARY SMEV (shape int, scale int, scale slope)
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [array of reals] zPost --> matrix of posterior samples
  #^*    2. [list of lists] priors_id --> priors list
  #^*    3. [character] dir_res --> path of folder where to save results
  #^*    4. [character] name_output --> string to differentiate results files
  #^*    5. [logical] save --> flag for saving plots
  #^* OUT
  #^*    1. [logical] convergence --> TRUE if good MCMC convergence
  #^****************************************************************************
  #^* REF.: uses coda, psych, stats, reshape2 packages
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  ##############################################################################
  
  #^****************************************************************************
  #^                        study of convergence
  #^****************************************************************************
  len=length(zPost[,1])
  gel=NULL
  convergence=T
  nam_param=c("shape","scale intercept","scale slope")
  # Gelman test (divide mcmc in 4 chains):
  for (i in 1:3){
    mcmc.list  = list(zPost[1:(len/4),i], 
                      zPost[(len/4+1):(len/2),i], 
                      zPost[(len/2+1):(len*3/4),i], 
                      zPost[(len*3/4+1):len,i])
    mcmc.list=as.mcmc.list(lapply(mcmc.list, coda::mcmc))
    gel[[i]]=coda::gelman.diag(x=mcmc.list,confidence=0.95,transform=F)
    par(mar=c(6,6,2,5))  # Set the margin on all sides 
    coda::gelman.plot(x=mcmc.list,auto.layout=F) 
    title(main=nam_param[i],cex.main=2)
    #acfplot(x = mcmc.list, lag.max=200)
    if (gel[[i]]$psrf[1]>1.2) { 
      # GELMAN FACTOR: limit is suggested around 1.2
      convergence=F  # ===> bad convergence !!!
    }
    dev.off()
  }
  
  #^****************************************************************************
  #^                             Plots (if save=T)
  #^****************************************************************************
  if (save){
    # Plot traceplots:
    png(filename=paste0(dir_res,"/traceplots_mcmc_",name_output,".png"),
        width=1000,height=700)
    par(mfrow=c(3,1))
    coda::traceplot(as.mcmc(zPost[,1]),main="shape",smooth=T,
                    cex.lab=2,cex.axis=1,cex.main=2,cex.sub=1.5)
    coda::traceplot(as.mcmc(zPost[,2]),main="scale intercept",smooth=T,
                    cex.lab=2,cex.axis=1,cex.main=2,cex.sub=1.5)
    coda::traceplot(as.mcmc(zPost[,3]),main="scale slope",smooth=T, 
                    cex.lab=2,cex.axis=1,cex.main=2,cex.sub=1.5)
    dev.off()
    
    # mcmc autocorrelation:
    png(paste0(dir_res,"/autocorr_mcmc_",name_output,".png"),
        width=1000,height=700)
    par(mfrow=c(3,1))
    for (i in 1:(ncol(zPost))) {
      par(mar=c(6,6,2,5))  # Set the margin on all sides 
      autocorr.plot(x=zPost[,i],auto.layout=F,
                    cex.lab=3,cex.axis=3,cex.main=3,cex.sub=2)
      title(main=nam_param[i],cex.main=2)
    }
    dev.off()
 
    png(paste0(dir_res,"/gelman_mcmc_",name_output,".png"),
        width=1000,height=700)
    par(mfrow=c(3,1))
    dev.off()
    
    # Compute and plot Gelman Reduction Factor:
    Rc="Reduction factor Rc"
    Ru="Upper confidence limit Ru" 
    limit="threshold"
    col.gelman=c("Reduction factor Rc"      ="blue", 
                 "Upper confidence limit Ru"="red", 
                 "threshold"                ="red")
    Rc.Ru.plot<-ggplot()+
      coord_cartesian(xlim=c(1,length(gel)))+
      scale_x_continuous(breaks=seq(1,length(gel),1))+
      geom_point(mapping=aes(x=seq(1,length(gel),1),
                             unlist(gel,use.names=F)[seq(1,length(gel)*2,2)], 
                             col=Rc),shape=0,size=4,stroke=0.5)+
      geom_point(mapping=aes(x=seq(1,length(gel),1),
                             unlist(gel,use.names=F)[seq(2,length(gel)*2,2)], 
                             col=Ru),shape=4,size=4,stroke=0.5)+
      theme_bw(base_size=15)+
      geom_hline(aes(yintercept=1.2,col=limit),lwd=1,linetype="dashed")+
      scale_colour_manual(name="Legend", 
                          values=col.gelman ,
                          breaks=c(Rc,Ru,limit),
                          guide=guide_legend(override.aes=list(
                            linetype=c("blank","blank","dashed"),
                            shape=c(0,4,NA)), 
                            title.position="top")) +
      xlab("Parameter index")+ylab("Gelman factor")+
      theme(text=element_text(size=15)
            ,panel.grid.major=element_line(size=0.4,linetype="dashed")
            ,panel.grid.minor=element_blank()
            ,legend.position="bottom"
            ,legend.text=element_text(size=10)
            ,legend.justification=c(0,1)
            ,legend.box="vertical"
            ,legend.box.just="left"
            ,legend.text.align=0)
    ggsave(Rc.Ru.plot,
           filename=paste0(dir_res,"/Rc_Ru_gelman_mcmc_",name_output,".png"),
           width=6,height=4,dpi=130)

    # Density plots mcmc:
    df=as.data.frame(zPost)
    names(df)=nam_param
    if (!is.null(priors_id)){
      # priors pdf:
      dPrior=matrix(NA,nrow=10000,ncol=3)
      grid=matrix(NA,10000,ncol(dPrior)) 
      for(i in 1:ncol(dPrior)){
        grid[,i]=seq(min(df[,i]),max(df[,i]),length.out=10000)
      }    
      # compute pdf of specified priors:
      # shape intercept:
      if (priors_id$shape_intercept[[1]]=='uniform'){
        dPrior[,1]=dunif(x=grid[,1], 
                         min=priors_id$shape_intercept[[2]], 
                         max=priors_id$shape_intercept[[3]])
      } else if (priors_id$shape_intercept[[1]]=='normal'){
        dPrior[,1]=dnorm(x=grid[,1], 
                         mean=priors_id$shape_intercept[[2]],
                         sd=priors_id$shape_intercept[[3]])
      } else if (priors_id$shape_intercept[[1]]=='lognormal'){
        dPrior[,1]=dlnorm(x=grid[,1], 
                          meanlog=round(Transf_Gauss_lognorm(
                            E=priors_id$shape_intercept[[2]], 
                            stdev=priors_id$shape_intercept[[3]])$mu,digits=4),
                          sdlog=round(Transf_Gauss_lognorm(
                            E=priors_id$shape_intercept[[2]],
                            stdev=priors_id$shape_intercept[[3]])$sd,digits =4))
      } else if (priors_id$shape_intercept[[1]]=='flatprior+'){
        dPrior[,1]=rep(0,grid[,1])
      } else {
        stop("prior distribution not available for shape intercept !!!")
      }
      # scale intercept:
      if (priors_id$scale_intercept[[1]]=='uniform'){
        dPrior[,2]=dunif(x=grid[,2], 
                         min=priors_id$scale_intercept[[2]], 
                         max=priors_id$scale_intercept[[3]])
      } else if (priors_id$scale_intercept[[1]]=='normal'){
        dPrior[,2]=dnorm(x=grid[,2], 
                         mean=priors_id$scale_intercept[[2]],
                         sd=priors_id$scale_intercept[[3]])
      } else if (priors_id$scale_intercept[[1]]=='lognormal'){
        dPrior[,2]=dlnorm(x=grid[,2], 
                          meanlog=round(Transf_Gauss_lognorm(
                            E=priors_id$scale_intercept[[2]],
                            stdev=priors_id$scale_intercept[[3]])$mu,digits=4),
                          sdlog=round(Transf_Gauss_lognorm(
                            E=priors_id$scale_intercept[[2]], 
                            stdev=priors_id$scale_intercept[[3]])$sd,digits=4))
      } else if (priors_id$scale_intercept[[1]]=='flatprior+'){
        dPrior[,2]=rep(0,grid[,2])
      } else {
        stop("prior distribution not available for scale intercept !!!")
      }
      # scale slope:
      if (priors_id$scale_slope[[1]]=='uniform'){
        dPrior[,3]=dunif(x=grid[,3], 
                         min=priors_id$scale_slope[[2]], 
                         max=priors_id$scale_slope[[3]])
      } else if (priors_id$scale_slope[[1]]=='normal'){
        dPrior[,3]=dnorm(x=grid[,3], 
                         mean=priors_id$scale_slope[[2]], 
                         sd=priors_id$scale_slope[[3]])
      } else if (priors_id$scale_slope[[1]]=='lognormal'){
        dPrior[,3]=dlnorm(x=grid[,3], 
                            meanlog=round(Transf_Gauss_lognorm(
                              E=priors_id$scale_slope[[2]], 
                              stdev=priors_id$scale_slope[[3]])$mu,digits=4),
                            sdlog=round(Transf_Gauss_lognorm(
                              E=priors_id$scale_slope[[2]], 
                              stdev=priors_id$scale_slope[[3]])$sd,digits=4))
      } else if (priors_id$scale_slope[[1]]=='flatprior+'){
        dPrior[,3]=rep(0,grid[,3])
      } else {
        stop("prior distribution not available for scale slope !!!")
      }
      # create unique dataframe of priors:
      dfprior=as.data.frame(dPrior)
      names(dfprior)=nam_param
      priorDF=dfprior
      priorDF.grid=data.frame(grid)
      names(priorDF.grid)=nam_param
      priorDF.grid.bis=cbind(priorDF.grid, indx=1:nrow(priorDF.grid))
      priorDF.grid.bis.m=melt(priorDF.grid.bis, id.vars='indx')
      priorDF.bis=cbind(priorDF, indx=1:nrow(priorDF))
      priorDF.bis.m=cbind(melt(priorDF.bis, id.vars='indx'), 
                          grid=priorDF.grid.bis.m$value)
    }
    X=cbind(df,indx=1:nrow(df))
    Xm=melt(X,id.vars='indx')
    # plot pdf:
    pdf_param = ggplot(Xm) +
      # geom_density(aes(x=zPost),color="blue",fill="blue",alpha=0.3,lwd=1.2)+
      geom_density(aes(value,fill=variable),colour=NA,alpha=0.8)
    if (!is.null(priors_id)){
      pdf_param=pdf_param+
        geom_line(data=priorDF.bis.m,aes(x=grid,y=value),colour="black")
    }
    pdf_param=pdf_param+
      theme_bw(base_size=15)+
      facet_wrap(~variable,scales='free',ncol=3)+
      theme(legend.position="none"
            ,axis.text=element_text(size=10)
            ,axis.title.x=element_blank()
            ,panel.grid.minor=element_blank())
    ggsave(pdf_param, 
           filename=paste0(dir_res,"/pdf_param_mcmc_",name_output,".png"),
           width=6,height=6,dpi=100)
    
    # Correlations between parameters
    png(filename=paste0(dir_res,"/corr_matrix_mcmc_",name_output,".png"),
        width=10,height=7,units='in',res=100)
    plt=psych::pairs.panels(df,
                            method="pearson", # correlation method
                            hist.col="#00AFBB", # color
                            density=T,  # show density plots
                            ellipses=T # show correlation ellipses
    )
    dev.off()
    # Plot scatterplots:
    # png(filename=paste0(dir_res,"/scatterplots_mcmc_",name_output,".png"), 
    #     width=1000,height=1000)
    # pairs(zPost)
    # dev.off()
  }
  return(convergence)
}










################################################################################
mcmc.diagn_stat_gev = function(zPost,
                               priors_id,
                               dir_res,
                               name_output,
                               save=T){
################################################################################
  #^* GOAL: diagnostics of mcmc samples (study of convergence)
  #^*       STATIONARY GEV (loc, scale, shape)
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [array of reals] zPost --> matrix of posterior samples
  #^*    2. [list of lists] priors_id --> priors list
  #^*    3. [character] dir_res --> path of folder where to save results
  #^*    4. [character] name_output --> string to differentiate results files
  #^*    5. [logical] save --> flag for saving plots
  #^* OUT
  #^*    1. [logical] convergence --> TRUE if good MCMC convergence
  #^****************************************************************************
  #^* REF.: uses coda, psych, stats, reshape2 packages
  #^****************************************************************************
  #^* to do: 
  #^****************************************************************************
  #^* COMMENTS:
  ##############################################################################
  
  #^****************************************************************************
  #^                        Study of convergence
  #^****************************************************************************
  len=length(zPost[,1])
  gel=NULL
  convergence=T
  nam_param=c("loc","scale","shape")
  # Gelman test (divide mcmc in 4 chains):
  for (i in 1:3){
    mcmc.list=list(zPost[1:(len/4),i], 
                   zPost[(len/4+1):(len/2),i], 
                   zPost[(len/2+1):(len*3/4),i], 
                   zPost[(len*3/4+1):len,i])
    mcmc.list=as.mcmc.list(lapply(mcmc.list, coda::mcmc))
    gel[[i]]=coda::gelman.diag(x=mcmc.list,confidence=0.95,transform=F)
    par(mar=c(6,6,2,5))  # Set the margin on all sides 
    coda::gelman.plot(x=mcmc.list,auto.layout=F,
                      cex.lab=2,cex.axis=2,cex.main=2,cex.sub=2) 
    title(main=names(zPost)[i], cex.main = 2)
    #acfplot(x=mcmc.list,lag.max=200)
    if (gel[[i]]$psrf[1] > 1.2) {  
      # GELMAN FACTOR: limit is suggested around 1.2
      convergence = F  # ===> bad convergence !!!
    }
  }
  dev.off()
  
  
  #^****************************************************************************
  #^                             Plots (if save=T)
  #^****************************************************************************
  if (save){
    # Plot traceplots:
    png(filename=paste0(dir_res,"/traceplots_mcmc_",name_output,".png"), 
        width=1000,height=700)
    par(mfrow=c(3,1))
    coda::traceplot(as.mcmc(zPost[,1]),main="loc",smooth=T, 
                    cex.lab=2,cex.axis=1,cex.main=2,cex.sub=1.5)
    coda::traceplot(as.mcmc(zPost[,2]),main="scale",smooth=T,
                    cex.lab=2,cex.axis=1,cex.main=2,cex.sub=1.5)
    coda::traceplot(as.mcmc(zPost[,3]),main="shape",smooth=T,
                    cex.lab=2,cex.axis=1,cex.main=2,cex.sub=1.5)
    dev.off()
  
    # mcmc autocorrelation:
    png(paste0(dir_res,"/autocorr_mcmc_",name_output,".png"), 
        width=1000,height=700)
    par(mfrow=c(3,1))
    for (i in 1:(ncol(zPost))){
      par(mar=c(6,6,2,5))  # Set the margin on all sides 
      coda::autocorr.plot(x=zPost[,i],auto.layout=F,
                          cex.lab=2,cex.axis=2,cex.main=2,cex.sub=2)
      # title(main=names(zPost)[i],cex.main=2)
    }
    dev.off()
  
    png(paste0(dir_res,"/gelman_mcmc_",name_output,".png"),
        width=1000,height=700)
    par(mfrow=c(3,1))
    dev.off()
    
    # Compute and plot Gelman Reduction Factor:
    Rc= "Reduction factor Rc"
    Ru="Upper confidence limit Ru" 
    limit="threshold"
    col.gelman=c("Reduction factor Rc"      ="blue", 
                 "Upper confidence limit Ru"="red", 
                 "threshold"                ="red")
    Rc.Ru.plot<-ggplot()+
      coord_cartesian(xlim=c(1,length(gel)))+
      scale_x_continuous(breaks=seq(1, length(gel),1))+
      geom_point(mapping=aes(x=seq(1,length(gel),1),
                             unlist(gel,use.names=F)[seq(1,length(gel)*2,2)], 
                             col=Rc),shape=0, size=4,stroke=0.5)+
      geom_point(mapping=aes(x=seq(1,length(gel),1), 
                             unlist(gel,use.names=F)[seq(2,length(gel)*2,2)], 
                             col=Ru),shape=4,size=4,stroke=0.5)+
      theme_bw(base_size=15)+
      geom_hline(aes(yintercept=1.2,col=limit),lwd=1,linetype="dashed")+
      scale_colour_manual(name="Legend", 
                          values=col.gelman,
                          breaks=c(Rc,Ru,limit),
                          guide=guide_legend(override.aes=list(
                            linetype=c("blank","blank","dashed"),
                            shape=c(0,4,NA)), 
                            title.position="top"))+
      xlab("Parameter index")+ylab("Gelman factor")+
      theme(text=element_text(size=15)
            ,panel.grid.major=element_line(size=0.4,linetype="dashed")
            ,panel.grid.minor=element_blank()
            ,legend.position="bottom"
            ,legend.text=element_text(size=10) 
            ,legend.justification=c(0,1)
            ,legend.box="vertical"
            ,legend.box.just="left"
            ,legend.text.align=0)
    ggsave(Rc.Ru.plot, 
           filename=paste0(dir_res,"/Rc_Ru_gelman_mcmc_",name_output,".png"),
           width=6,height=4,dpi=130)
    Rc.Ru.plot

    # Density plots mcmc:
    df=as.data.frame(zPost)
    names(df)=nam_param
    if (!is.null(priors_id)){
      # priors pdf:
      dPrior=matrix(NA, nrow= 10000, ncol=3)
      grid=matrix(NA, 10000, ncol(dPrior)) 
      for(i in 1:ncol(dPrior)){
        grid[,i]=seq(min(df[,i]), max(df[,i]),length.out=10000)
      }    
      # compute pdf of specified priors:
      # loc:
      if (priors_id$loc[[1]]=='uniform'){
        dPrior[,1]=dunif(x=grid[,1], 
                         min=priors_id$loc[[2]], 
                         max=priors_id$loc[[3]])
      } else if (priors_id$loc[[1]]=='normal'){
        dPrior[,1]=dnorm(x=grid[,1], 
                         mean=priors_id$loc[[2]], 
                         sd=priors_id$loc[[3]])
      } else if (priors_id$loc[[1]]=='lognormal'){
        dPrior[,1]=dlnorm(x=grid[,1], 
                          meanlog=round(Transf_Gauss_lognorm(
                            E=priors_id$loc[[2]], 
                            stdev=priors_id$loc[[3]])$mu,digits=4),
                          sdlog=round(Transf_Gauss_lognorm(
                            E=priors_id$loc[[2]],
                            stdev=priors_id$loc[[3]])$sd,digits=4))
      } else if (priors_id$loc[[1]]=='flatprior+'){
        dPrior[,1]=rep(0,length(grid[,1]))
      } else {
        stop("prior distribution not available for location !!!")
      }
      # scale:
      if (priors_id$scale[[1]]=='uniform'){
        dPrior[,2]=dunif(x=grid[,2],
                         min=priors_id$scale[[2]],
                         max=priors_id$scale[[3]])
      } else if (priors_id$scale[[1]]=='normal'){
        dPrior[,2]=dnorm(x=grid[,2],
                         mean=priors_id$scale[[2]], 
                         sd=priors_id$scale[[3]])
      } else if (priors_id$scale[[1]]=='lognormal'){
        dPrior[,2]=dlnorm(x=grid[,2], 
                          meanlog=round(Transf_Gauss_lognorm(
                            E=priors_id$scale[[2]], 
                            stdev=priors_id$scale[[3]])$mu,digits=4),
                          sdlog=round(Transf_Gauss_lognorm(
                            E=priors_id$scale[[2]],
                            stdev=priors_id$scale[[3]])$sd,digits=4))
      } else if (priors_id$scale[[1]]=='flatprior+'){
        dPrior[,2]=rep(0,length(grid[,2]))
      } else {
        stop("prior distribution not available for scale !!!")
      }
      # shape:
      if (priors_id$shape[[1]]=='uniform'){
        dPrior[,3]=dunif(x=grid[,3], 
                         min=priors_id$shape[[2]], 
                         max=priors_id$shape[[3]])
      } else if (priors_id$shape[[1]]=='normal'){
        dPrior[,3]=dnorm(x=grid[,3],
                         mean=priors_id$shape[[2]], 
                         sd=priors_id$shape[[3]])
      } else if (priors_id$shape[[1]]=='lognormal'){
        dPrior[,3]=dlnorm(x=grid[,3], 
                          meanlog=round(Transf_Gauss_lognorm(
                            E=priors_id$shape[[2]],
                            stdev=priors_id$shape[[3]])$mu,digits=4),
                          sdlog=round(Transf_Gauss_lognorm(
                            E=priors_id$shape[[2]], 
                            stdev=priors_id$shape[[3]])$sd,digits=4))
      } else if (priors_id$shape[[1]]=='flatprior+'){
        dPrior[,3]=rep(0,length(grid[,3]))
      } else {
        stop("prior distribution not available for shape !!!")
      }
      dfprior=as.data.frame(dPrior)
      names(dfprior)=nam_param
      priorDF=dfprior
      priorDF.grid=data.frame(grid)
      names(priorDF.grid)=nam_param
      priorDF.grid.bis=cbind(priorDF.grid, indx=1:nrow(priorDF.grid))
      priorDF.grid.bis.m=melt(priorDF.grid.bis, id.vars='indx')
      priorDF.bis=cbind(priorDF,indx=1:nrow(priorDF))
      priorDF.bis.m=cbind(melt(priorDF.bis, id.vars='indx'), 
                          grid=priorDF.grid.bis.m$value)
    }
    X=cbind(df,indx=1:nrow(df))
    Xm=melt(X,id.vars='indx')
    # plot pdf:
    pdf_param = ggplot(Xm)+
      geom_density(aes(value,fill=variable),colour=NA,alpha=0.8)
    if (!is.null(priors_id)){
      pdf_param=pdf_param+
        geom_line(data=priorDF.bis.m,aes(x=grid,y=value),colour="black")
    }
    pdf_param = pdf_param+
      theme_bw(base_size=15)+
      facet_wrap(~variable,scales='free',ncol=3)+
      theme(legend.position="none"
            ,axis.text=element_text(size=10)
            ,axis.title.x=element_blank()
            ,panel.grid.minor=element_blank())
    ggsave(pdf_param,
           filename=paste0(dir_res,"/pdf_param_mcmc_",name_output,".png"),
           width=6,height=6,dpi=100)
    
    # Plot correlations between parameters
    png(filename=paste0(dir_res, "/corr_matrix_mcmc_", name_output,".png"),
        width=10,height=7,units='in',res=100)
    plt=psych::pairs.panels(df,
                            method="pearson", # correlation method
                            hist.col="#00AFBB", # color
                            density=T,  # show density plots
                            ellipses=T) # show correlation ellipses
    dev.off()
    # Plot scatterplots:
    # png(filename=paste0(dir_res,"/scatterplots_mcmc_",name_output,".png"),
    #     width=1000,height=1000)
    # pairs(zPost)
    # dev.off()
  }
  return(convergence)
}









################################################################################
model_comparison = function(s                = s, 
                            flag_smev_ns     = T, 
                            flag_smev_stat   = T,
                            flag_smev_partns = F,                           
                            flg_save         = T,
                            dur              = 1, 
                            sat_prod         = NULL,
                            pix_row          = NULL,
                            pix_col          = NULL,
                            dir.res_pix      ='',
                            nfplot           ='.png',
                            title_plot       = ''){
################################################################################  
  #^* GOAL: models comparison (plot and save information criteria)
  #^****************************************************************************
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy
  #^****************************************************************************
  #^* CREATED/MODIFIED: 2025/11/28
  #^****************************************************************************
  #^* IN
  #^*    1. [list object] s --> main object with all results of the analysis 
  #^*    2. [logical] flag_smev_ns --> flag if smev full NS is activated
  #^*    3. [logical] flag_smev_stat --> flag if smev stationary is activated
  #^*    3. [logical] flag_smev_partns --> flag if smev partial NS is activated
  #^*    5. [logical] flg_save --> flag to save the plots and criteria to "s"
  #^*    6. [integer] dur --> duration index
  #^*    7. [character] sat_prod --> if satellite product, name of product
  #^*    8. [integer] pix_row --> if gridded dataset, row index
  #^*    9. [integer] pix_col --> if gridded dataset, column index
  #^*    10. [character] dir.res_pix --> path of folder where to save plots
  #^*    11. [character] nfplot --> filename for plots
  #^*    12. [character] title_plot --> title for plots
  #^* OUT
  #^*    1. [list of plots] modelselectionplot -->
  #^*    2. [list of plots] bayesfactorplot --> 
  #^*    3. [list object] s --> updated s with criteria and BF plots
  #^****************************************************************************
  #^* REF.: uses coda, psych, stats, reshape2 packages
  #^****************************************************************************
  #^* to do: BAYES FACTOR is for instance computed using the BIC, next steps are
  #^*        to compute marginal likelihood.
  #^****************************************************************************
  #^* COMMENTS:
  ##############################################################################
  
  
  #^****************************************************************************
  #                      LIKELIHOOD RATIO TEST 
  #^****************************************************************************
  # in case of Maximum Likelihood Estimation approach:
  if (flag_smev_ns & flag_smev_stat){
    # full ns model with respect to stationary model:
    # using stan algorithm (it is NULL if flag_like=F):
    s$pval_yy_wrt_nn[[dur]]=likelihood_ratio_test(
      loglik_H0=s$s_stat_tmp[[dur]]$maxlik_mcmc,
      loglik_H1=s$s_ns_tmp[[dur]]$maxlik_mcmc,
      degrees_of_freedom=2)
    # using fminsearch function (MLE method):
    s$pval_yy_wrt_nn_fmins[[dur]]=likelihood_ratio_test(
      loglik_H0=-s$s_stat_tmp[[dur]]$maxlik_fmins,
      loglik_H1=-s$s_ns_tmp[[dur]]$maxlik_fmins,
      degrees_of_freedom=2)
  }
  if (flag_smev_stat & flag_smev_partns){
    # partial-nonstat. model with respect to stationary model:
    # using stan algorithm (it is NULL if flag_like=F):
    s$pval_ny_wrt_nn[[dur]]=likelihood_ratio_test(
      loglik_H0=s$s_stat_tmp[[dur]]$maxlik_mcmc,
      loglik_H1=s$s_partns_tmp[[dur]]$maxlik_mcmc,
      degrees_of_freedom=1)
    # using fminsearch function (MLE method):
    s$pval_ny_wrt_nn_fmins[[dur]]=likelihood_ratio_test(
      loglik_H0=-s$s_stat_tmp[[dur]]$maxlik_fmins,
      loglik_H1=-s$s_partns_tmp[[dur]]$maxlik_fmins,
      degrees_of_freedom=1)
  }
  if (flag_smev_ns & flag_smev_partns){
    # full nonstat. model with respect to partial-nonstat model:
    # using stan algorithm (it is NULL if flag_like=F):
    s$pval_yy_wrt_ny[[dur]]=likelihood_ratio_test(
      loglik_H0=s$s_partns_tmp[[dur]]$maxlik_mcmc,
      loglik_H1=s$s_ns_tmp[[dur]]$maxlik_mcmc,
      degrees_of_freedom=1)
    # using fminsearch function (MLE method):
    s$pval_yy_wrt_ny_fmins[[dur]]=likelihood_ratio_test(
      loglik_H0=-s$s_partns_tmp[[dur]]$maxlik_fmins,
      loglik_H1=-s$s_ns_tmp[[dur]]$maxlik_fmins,
      degrees_of_freedom=1)
  }
  
  
  #^****************************************************************************
  #              compare criteria DIC/BIC/AIC/WAIC/LOO:
  #^****************************************************************************
  # If one of the Bayesian approaches is used (analyze the posterior sample), 
  # compute also the model selection criteria (e.g., BIC, DIC, ...):
  modelselectionplot=NULL
  if (any(!is.null(s$s_ns_tmp[[dur]]$BIC),
          !is.null(s$s_partns_tmp[[dur]]$BIC), 
          !is.null(s$s_stat_tmp[[dur]]$BIC))){
    # Get criteria and their min values:
    BIC=c(s$s_ns_tmp[[dur]]$BIC,   
          s$s_partns_tmp[[dur]]$BIC,   
          s$s_stat_tmp[[dur]]$BIC)
    BIC_b=c(s$s_ns_tmp[[dur]]$BIC_b,  
            s$s_partns_tmp[[dur]]$BIC_b,  
            s$s_stat_tmp[[dur]]$BIC_b)
    AIC=c(s$s_ns_tmp[[dur]]$AIC,   
          s$s_partns_tmp[[dur]]$AIC,   
          s$s_stat_tmp[[dur]]$AIC)
    AIC_b=c(s$s_ns_tmp[[dur]]$AIC_b,  
            s$s_partns_tmp[[dur]]$AIC_b,
            s$s_stat_tmp[[dur]]$AIC_b)
    HQC=c(s$s_ns_tmp[[dur]]$HQC,    
          s$s_partns_tmp[[dur]]$HQC,  
          s$s_stat_tmp[[dur]]$HQC)
    HQC_b=c(s$s_ns_tmp[[dur]]$HQC_b,  
            s$s_partns_tmp[[dur]]$HQC_b,
            s$s_stat_tmp[[dur]]$HQC_b)
    DIC1=c(s$s_ns_tmp[[dur]]$DIC1,  
           s$s_partns_tmp[[dur]]$DIC1,  
           s$s_stat_tmp[[dur]]$DIC1)
    DIC1_b=c(s$s_ns_tmp[[dur]]$DIC1_b, 
             s$s_partns_tmp[[dur]]$DIC1_b,
             s$s_stat_tmp[[dur]]$DIC1_b)
    DIC2=c(s$s_ns_tmp[[dur]]$DIC2,  
           s$s_partns_tmp[[dur]]$DIC2,  
           s$s_stat_tmp[[dur]]$DIC2)
    DIC2_b=c(s$s_ns_tmp[[dur]]$DIC2_b,
             s$s_partns_tmp[[dur]]$DIC2_b,
             s$s_stat_tmp[[dur]]$DIC2_b)
    DIC3=c(s$s_ns_tmp[[dur]]$DIC3,
           s$s_partns_tmp[[dur]]$DIC3,
           s$s_stat_tmp[[dur]]$DIC3)
    WAIC=c(s$s_ns_tmp[[dur]]$WAIC[3],
           s$s_partns_tmp[[dur]]$WAIC[3],
           s$s_stat_tmp[[dur]]$WAIC[3])
    pWAIC=c(s$s_ns_tmp[[dur]]$WAIC[2],
            s$s_partns_tmp[[dur]]$WAIC[2],
            s$s_stat_tmp[[dur]]$WAIC[2])
    elpdWAIC=c(s$s_ns_tmp[[dur]]$WAIC[1],
               s$s_partns_tmp[[dur]]$WAIC[1],
               s$s_stat_tmp[[dur]]$WAIC[1])
    LOOIC=c(s$s_ns_tmp[[dur]]$all_ELPD[3],
            s$s_partns_tmp[[dur]]$all_ELPD[3],
            s$s_stat_tmp[[dur]]$all_ELPD[3])
    ploo=c(s$s_ns_tmp[[dur]]$all_ELPD[2],
           s$s_partns_tmp[[dur]]$all_ELPD[2],
           s$s_stat_tmp[[dur]]$all_ELPD[2])
    epldLOO=c(s$s_ns_tmp[[dur]]$all_ELPD[1],
              s$s_partns_tmp[[dur]]$all_ELPD[1],
              s$s_stat_tmp[[dur]]$all_ELPD[1])
    
    # save to results object "s":
    if (flg_save){
    ##############
      # create dataframe with the different criteria:
      x=c()
      if (flag_smev_ns){x=c(x,"nonstat")}
      if (flag_smev_partns){x=c(x,"partial nonstat")}
      if (flag_smev_stat){x=c(x,"stat")}
      
      s$criteria.df[[dur]]=data.frame(x=x)
      if(!is.null(BIC)){ 
        s$criteria.df[[dur]]$BIC=BIC
        s$criteria.df[[dur]]$BICmin=which.min(BIC) 
      }
      if(!is.null(BIC_b)){
        s$criteria.df[[dur]]$BIC_b=BIC_b
        s$criteria.df[[dur]]$BICmin_b=which.min(BIC_b) 
      }
      if(!is.null(AIC)){
        s$criteria.df[[dur]]$AIC=AIC
        s$criteria.df[[dur]]$AICmin=which.min(AIC)
      }
      if(!is.null(AIC_b)){
        s$criteria.df[[dur]]$AIC_b=AIC_b
        s$criteria.df[[dur]]$AICmin_b=which.min(AIC_b)
      }
      if(!is.null(HQC)){ 
        s$criteria.df[[dur]]$HQC=HQC
        s$criteria.df[[dur]]$HQCmin=which.min(HQC) 
      }
      if(!is.null(HQC_b)){
        s$criteria.df[[dur]]$HQC_b=HQC_b
        s$criteria.df[[dur]]$HQCmin_b=which.min(HQC_b)
      }
      if(!is.null(DIC1)){
        s$criteria.df[[dur]]$DIC1=DIC1
        s$criteria.df[[dur]]$DIC1min=which.min(DIC1) 
      } 
      if(!is.null(DIC1_b)){ 
        s$criteria.df[[dur]]$DIC1_b=DIC1_b
        s$criteria.df[[dur]]$DIC1min_b= which.min(DIC1_b) 
      }
      if(!is.null(DIC2)){ 
        s$criteria.df[[dur]]$DIC2=DIC2
        s$criteria.df[[dur]]$DIC2min= which.min(DIC2) 
      } 
      if(!is.null(DIC2_b)){ 
        s$criteria.df[[dur]]$DIC2_b=DIC2_b
        s$criteria.df[[dur]]$DIC2min_b= which.min(DIC2_b) 
      }
      if(!is.null(DIC3)){
        s$criteria.df[[dur]]$DIC3=DIC3
        s$criteria.df[[dur]]$DIC3min=which.min(DIC3) 
      }
      if(!is.null(elpdWAIC)){ 
        s$criteria.df[[dur]]$elpdWAIC=elpdWAIC
        s$criteria.df[[dur]]$elpdWAICmin=which.min(elpdWAIC) 
      }
      if(!is.null(pWAIC)){ 
        s$criteria.df[[dur]]$pWAIC=pWAIC
        s$criteria.df[[dur]]$pWAICmin=which.min(pWAIC) 
      } 
      if(!is.null(WAIC)){ 
        s$criteria.df[[dur]]$WAIC=WAIC
        s$criteria.df[[dur]]$WAICmin=which.min(WAIC) 
      }
      if(!is.null(LOOIC)){ 
        s$criteria.df[[dur]]$LOOIC=LOOIC
        s$criteria.df[[dur]]$LOOICmin=which.min(LOOIC)
      }
      
      # plot criteria values:
      if (ncol(s$criteria.df[[dur]])>1){
        title1=paste0('Model selection criteria on SMEV models -- ',
                      s$durations[dur]/60, 'h [mm/h] \n',title_plot)
        name_output1=paste0(dir.res_pix,"/Criteria_",s$durations[dur]/60,"h",nfplot)
        # plot:
        modelselectionplot=model.selection.plot(criteria.df=s$criteria.df[[dur]],
                                                title=title1,
                                                name_output=name_output1)
        # add "duration" to the criteria plot object:
        modelselectionplot$duration=paste0(durations[dur]/60,'h')
      }
    }
  }
      
      
      
  #^****************************************************************************
  # BAYES FACTOR:
  #^****************************************************************************
  # not implemented the original formulation with the marginal Likelihood, yet!
  # A simplified version of Bayes Factor from ERIC-JAN WAGENMAKERS (2007),
  # "A practical solution to the pervasive problems of p values".
  # where BF is computed from the BIC values of the two comparing models.
  s$BF_yy_wrt_nn=exp((BIC[3]-BIC[1])/2)
  s$BF_yy_wrt_ny=exp((BIC[2]-BIC[1])/2)
  s$BF_ny_wrt_nn=exp((BIC[3]-BIC[2])/2)
  s$BF_nn_wrt_yy=exp((BIC[1]-BIC[3])/2)
  s$BF_ny_wrt_yy=exp((BIC[1]-BIC[2])/2)
  s$BF_nn_wrt_ny=exp((BIC[2]-BIC[3])/2)
  bayesfactorplot=NULL
  
  # plot BF values:
  if (flg_save){
  ##############
      x=c()
      BF=c()
      if (flag_smev_ns & flag_smev_stat){
        x=c(x, "nonstat vs. stat")
        BF=c(BF, s$BF_yy_wrt_nn)
      }
      if (flag_smev_ns & flag_smev_partns){
        x=c(x,"nonstat vs. partial nonstat")
        BF=c(BF, s$BF_yy_wrt_ny)
      }
      if (flag_smev_partns & flag_smev_stat){
        x=c(x,"partial nonstat vs. stat")
        BF=c(BF, s$BF_ny_wrt_nn)
      }

      if (length(x)>0){
        s$bayes_factor.df[[dur]]=data.frame(x=x, BF=BF)
        
      if (ncol(s$bayes_factor.df[[dur]])>1){
        s$bayes_factor.df[[dur]]$x=x
        # if (nrow(s$bayes_factor.df[[dur]])){
        title1=paste0('Bayes Factor on SMEV models -- duration=',
                      s$durations[dur]/60,'h [mm/h] \n',title_plot)
        name_output1=paste0(dir.res_pix,"/BF_",s$durations[dur]/60,"h",nfplot)
        # plot
        bayesfactorplot=bayesfactor.plot(bayes_factor.df=s$bayes_factor.df[[dur]],
                                         title=title1,
                                         name_output=name_output1)
        # add "duration" to the BayesFactor plot object:
        bayesfactorplot$duration=paste0(durations[dur]/60,'h')
        }
      }
  }

  return(list(modelselectionplot=modelselectionplot,
              bayesfactorplot=bayesfactorplot))
}







