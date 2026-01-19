# needed packages:
library(adaptMCMC) # adaptive MCMC metropolis algorithm
library(mcmc) # metropolis algorithm
library(evd) # with the GEV pdf, cdf, ...
library(ggplot2) # for plots
library(bayesplot) # for mcmc diagnostic plots
library(rstan) # STAN
library(cmdstanr) # CMD STAN 
library(posterior) # posterior mcmc diagnostics
library(bayesplot) # posterior mcmc diagnostics
library(loo) # posterior mcmc diagnostics
library(PEIP)
library(coda)



################################################################################
#^*                               READ ME: 
################################################################################
#^* STAN package for Bayesian inference (HMC algorithm) requires prior 
#^* installation "install_cmdstan(cores=2)" and creation of a local directory 
#^* where the STAN models will be saved and compiled:
#^*    install.packages("cmdstanr",repos=c('https://stan-dev.r-universe.dev',
#^*                      getOption("repos")))
#^*    check_cmdstan_toolchain(fix=T,quiet=T)
#^*    install_cmdstan(cores=2)
#^*    cmdstan_path()
#^*    cmdstan_version()
################################################################################





################################################################################
smev_stat= function(algorithm="fminsearch", #"fminsearch","rstan","cmdstan","adaptMCMC"
                    Nmcmc=4000,
                    Nmcmc_max=40000,
                    nburn=0.5,
                    Nslim=1, 
                    scales=c(0.02,0.002,0.1,0.005),
                    priors_id=list(shape_intercept=list("lognormal",log(0.7),0.4,0.7),
                                   scale_intercept=list("lognormal",1.45-0.7*log(1),0.5,exp(1.45-0.7*log(1)))),
                    acc.rate=0.234,
                    y=c(), 
                    t=c(), 
                    thr=NA,
                    dir.res_id="", 
                    dur=1, 
                    flag_save=T,
                    flag_like=F,
                    flag_bayes=T){
################################################################################
  #^* GOAL: estimation of SMEV parameters STATIONARY MODEL                     #
  #^*       both scale and shape constant with t.                              #
  #^***************************************************************************#
  #^* PROGRAMMER: Matteo Darienzo, University of Padua, Italy                  #
  #^***************************************************************************#
  #^* CREATED/MODIFIED: 2025/11/28                                             #
  #^***************************************************************************#
  #^*                               INPUTS                                     #
  #^***************************************************************************#
  #^*   1. [character] algorithm --> algorithms for the parameters estimation: #
  #^*                               "cmdstan", "fminsearch","rstan",adaptMCMC" #
  #^*   2. [integer] Nmcmc --> tot number of MCMC simulations, e.g., 4000      # 
  #^*   3. [integer] Nmcmc_max --> (only if "adaptMCMC" algorithm is selected) #
  #^*                               max total number of MCMC sim, e.g., 30000  # 
  #^*   4 . [real] nburn --> burn first fraction of the MCMC iterations to     #
  #^*                        reduce impact of start value on mcmc convergence. #    
  #^*   5. [integer] Nslim --> (only if "adaptMCMC" algorithm is selected)     #
  #^*                           slim mcmc for posterior statistics (to reduce  #
  #^*                           lag-autocorrelation)                           # 
  #^*   6. [array of reals] scales -->(only if "adaptMCMC" algorithm is chosen #
  #^*                                    scaling setting of mcmc jumps to      #
  #^*                                    explore SMEV parameters distribution: #
  #^*                                    =c(0.02, 0.002, 0.1, 0.005)           #
  #^*   7. [list of lists] priors_id --> list of priors for SMEV parameters,   #
  #^*                                  independently if stat or nonstat model  # 
  #^*                                  (shape intercept, shape slope,          #
  #^*                                   scale intercpet, scale slope).         #
  #^*                                  one list for each selected duration, and#
  #^*                                  one list for each parameter, including  #
  #^*                                  distrib.name,mean,stdev,init value, ex: #
  #^*                                  priors_id=list(); priors_id[[1]]=list() #
  #^*                                  priors_id[[1]]$shape_intercept=list(    #
  #^*                                    "lognormal",log(0.7),0.4,0.7)         #
  #^*                                  priors_id[[1]]$shape_slope=list(        #
  #^*                                     "normal",0,0.05,0)                   #
  #^*                                  priors_id[[1]]$scale_intercept=list(...)#    
  #^*                                  priors_id[[1]]$scale_slope=list(...)    #  
  #^*   8.  [real] acc.rate --> (only if "adaptMCMC" algorithm is selected)    #
  #^*                           targeted acceptance rate for mcmc iter. (0.234)#
  #^*   9.  [array of reals] y --> array or ordinary events values             #
  #^*   10. [array of integers] t --> array of ordinary events years (from 0)  #
  #^*   11. [character] thr --> threshold percentile for left censoring (0.9)  #
  #^*   12. [character] dir.res_id --> path of folder with results             #
  #^*   13. [array of integers] dur --> durations (in minutes) of event window #
  #^*                                        =c(1,3,6,24)*60                   #
  #^*   14. [logical] flag_save --> flag for saving results                    #
  #^*   15. [logical] flag_like --> flag to activate maximum likelihood method #
  #^*   16. [logical] flag_bayes --> flag to activate Bayesain maxpost method  #
  #^*   17. [extra] ...                                                        #
  #^***************************************************************************#
  #^*                               OUTPUTS                                    #
  #^***************************************************************************#
  #^*   1. [list] s_tmp -->list object containing results of estimation:       #
  #^*                     "theta_maxlik_fmins" (set of param maximizing Lik    #
  #^*                     "maxlik_fmins" (Maximum Likelihood from fminsearch)  #
  #^*                                                                          #
  #^*                     "acc.rate" (input targeted acceptance rate)          #
  #^*                                                                          #
  #^*                     "mcmc_Lik" (mcmc likelihood sample)                  #
  #^*                     "zLik" (sampled mcmc likelihood parameters)          #
  #^*                     "zLik_acceptance"(Likelihood MCMC acceptance rate)   #
  #^*                     "theta_maxlik_mcmc" (set of param maximizing Lik mcmc#
  #^*                     "maxlik_mcmc" (Maximum Likelihood from MCMC)         # 
  #^*                     "a_w_quant_Like" (shape inter 5,50,95% qnt from like)#
  #^*                     "b_w_quant_Like" (shape slope 5,50,95% qnt from like)#
  #^*                     "a_C_quant_Like" (scale inter 5,50,95% qnt from like)#
  #^*                     "b_C_quant_Like" (sclae slope 5,50,95% qnt from like)#
  #^*                                                                          #
  #^*                     "mcmc_Post" (parameters posterior sample)            #
  #^*                     "zPost" (sampled posterior parameters)               #
  #^*                     "zPost_p" (sampled log-posterior)                    #
  #^*                     "zPost_acceptance"(posterior MCMC acceptance rate)   #
  #^*                     "init_values_post (initial values of MCMC chains)    #
  #^*                     "converg_Post" (convergence of posterior mcmc)       #
  #^*                     "a_w_quant_Post" (shape inter 5,50,95% qnt from post)#
  #^*                     "b_w_quant_Post" (shape slope 5,50,95% qnt from post)#
  #^*                     "a_C_quant_Post" (scale inter 5,50,95% qnt from post)#
  #^*                     "b_C_quant_Post" (scale slope 5,50,95% qnt from post)#
  #^*                                                                          #
  #^*                     "DIC1", "DIC1_b" (two versions of DIC from           # 
  #^*                     "DIC2", "DIC2_b" (two versions of DIC from           #
  #^*                     "DIC3" (DIC crit. from                               #
  #^*                     "AIC", "AIC_b" (two versions of AIC, Aikake,1974)    # 
  #^*                     "BIC", "BIC_b" (two versions of BIC, Schwarz,1978)   #
  #^*                     "HQC", "HQC_b" (Hannan-Quinn 1979)                   #
  #^***************************************************************************#
  #^*                          ACKNOWLEDGEMENTS:                               #
  #^***************************************************************************#
  #^* For the extreme analysis it uses methods: MEV (Marani and Ignaccolo,2015)#
  #^* and its simplified version SMEV (Marra et al., 2019. 2020).              #
  #^* It uses open source codes of Francesco Marra (University of Padua) for   #
  #^* storm separation, identification of ordinary events, SMEV method and     #
  #^* left-censoring approach, and non-stationary model freely available at:   #
  #^* https://zenodo.org/records/11934843 (SMEV methodology)                   #
  #^* https://zenodo.org/records/15047817 (non-stationary SMEV)                #
  #^*                                                                          #
  #^* For the Bayesian and MCMC algorithms it uses some ideas from:            #
  #^* - RsTooDS of Benjamin Renard (INRAE): https://zenodo.org/records/5075760 #
  #^* - BaM of Benjamin Renard (INRAE): https://github.com/BaM-tools/RBaM      #
  #^* - BayDERS of Matteo Darienzo(INRAE):                                     #
  #^*   https://github.com/MatteoDarienzo/BayDERS                              #
  #^* - packages "rstan" and "cmdstanr" for stan Bayesian inference + MCMC alg.#
  #^* - package "adaptMCMC" for adaptive metropolis algorithm from R package of#
  #^*   A. Scheidegger (doi:10.1007/s11222-011-9269-5)                         #
  #^***************************************************************************#
  #^*                               READ ME:                                   #
  #^***************************************************************************#
  #^* STAN package for Bayesian inference (HMC algorithm) requires prior       #
  #^* installation "install_cmdstan(cores=2)" and creation of a local directory#
  #^* where the STAN models will be saved and compiled.                        #
  #^* install.packages("cmdstanr",repos=c('https://stan-dev.r-universe.dev',   #
  #^*                   getOption("repos")))                                   #
  #^* check_cmdstan_toolchain(fix=T,quiet=T)                                   #
  #^* install_cmdstan(cores=2)                                                 #
  #^* cmdstan_path()                                                           #
  #^* cmdstan_version()                                                        #
  #^***************************************************************************#
  #^*                               TO DO                                      #
  #^***************************************************************************#
  #^* marginal likelihood for Bayes Factor                                     #
  #^* PSIS-LOO leave-one-out cross-validation                                  #
  #^* WAIC criterion                                                           #
  #^***************************************************************************#
  #^*                              DISCLAIMER:                                 #
  #^***************************************************************************#
  #^* Please, notice that is an experimental code. Further analysis is required#
  #^* for validating the proposed tools and scripts. We are not responsible for#
  #^* any loss (e.g.,of data, profits, business, customers, orders, other costs#
  #^* or disturbances) derived by their use in the operational practice.       #
  #^* The authorized user accepts the risks associated with using the software #
  #^* given its nature of free software. It is reserved to expert Users        #
  #^* (developers or professionals) having knowledge in computer science,      #
  #^* hydrology, extreme analysis and statistics.                              #
  ##############################################################################
  rm(s_tmp)
  s_tmp=NULL
  dir.create(dir.res_id,recursive=T)
  
  
  #^* Selected algorithm for parameters estimation:
  #^* choice between "cmdstan" and "rstan" for Bayesian estimation with STAN
  #^* package and Hamiltonian MCMC algorithm.
  #^* "cmdstan" is preferred, since computationally more efficient.
  #^* "fminsearch" for MLE method with a minimization optimization function.
  #^* "adaptMCMC" for Bayesian estimation with adaptive metropolis algorithm.
  
  
  #^****************************************************************************
  if (algorithm=="fminsearch"){
  #^****************************************************************************
    print("Max Likelihood with minimisation function")
    # Use fminsearch function (as in matlab SMEV version of Marra 2019) to 
    # compute maximum likelihood method with minimisation optimization function:
    fmins=NULL
    initval=c(priors_id$shape_intercept[[4]], # shape start value
              priors_id$scale_intercept[[4]]) # scale start value
    fmins=pracma::fminsearch(f=logLike1_stat,
                             x0=initval, 
                             y=y,
                             t=t,
                             thr=thr,
                             maxiter=5000) # HARDCODED maxiter!!
    s_tmp$theta_maxlik_fmins=fmins$xmin
    s_tmp$maxlik_fmins=fmins$fmin
  
    
    
    
    
    
  #^****************************************************************************
  } else if (algorithm =="adaptMCMC"){
  #^****************************************************************************
    
    ################
    if (flag_like){
    ################
      # MCMC metropolis algorithm (maximum LIKELIHOOD method):
      print("Max Likelihood")
      converg=F
      Nmcmc_tmp=Nmcmc/nburn*5 #accounting for thinning each 5 mcmc.
      zLik_ALL=zLik_p_ALL=c()
      it=0
      while ((converg ==F) & (Nmcmc_tmp<Nmcmc_max)){
        it=it+1
        set.seed(1)
        initval= c(priors_id$shape_intercept[[4]],priors_id$scale_intercept[[4]])
        mcmc_Lik=adaptMCMC::MCMC(logLike_stat, 
                                 n               = Nmcmc_tmp, 
                                 init            = initval,
                                 scale           = c(scales[1],scales[3]), 
                                 adapt           = T, 
                                 showProgressBar = F,
                                 acc.rate        = acc.rate, 
                                 y               = y,
                                 t               = t,
                                 thr             = thr)
        
        # zLik_ALL=rbind(zLik_ALL, mcmc_Lik$batch)
        zLik_ALL=rbind(zLik_ALL, mcmc_Lik$samples)
        zLik_p_ALL=rbind(zLik_p_ALL, mcmc_Lik$log.p)
        
        zLik= zLik_ALL[(nburn*nrow(zLik_ALL)+1):nrow(zLik_ALL),]
        zLik_p= zLik_p_ALL[(nburn*length(zLik_p_ALL)+1):length(zLik_p_ALL)]
        
        zLik= apply(zLik, 2, function(x) x[seq(1, length(x), Nslim)])
        zLik_p= zLik_p[seq(1, length(zLik_p), Nslim)]
        
        # print(paste0("------- accept.rate = ", floor(mcmc_Lik$accept*100), " %"))
        converg=mcmc.diagn_stat(zPost=zLik, 
                                priors_id=NULL,
                                dir_res=dir.res_id,
                                name_output=paste0('lik_',it), 
                                save=T)
        s_tmp$converg_Lik= converg
        converg=T
        if (converg==F){
          print(paste0("-------------------- bad converg.: increasing mcmc"))
          Nmcmc_tmp=30000
        }
      }
      s_tmp$mcmc_Lik=mcmc_Lik
      s_tmp$zLik=zLik
      s_tmp$zLik_acceptance=floor(mcmc_Lik$accept*100)
      
      # Compute likelihood to find maximum likelihood:
      # like=c()
      # for (tt in 1:nrow(zLik)){
      #   like[tt]=logLike(theta=zLik[tt,],y=y,t=t,thr=thr)
      # }
      
      s_tmp$theta_maxlik_mcmc=zLik[which.max(zLik_p),]
      s_tmp$maxlik_mcmc=max(zLik_p)
      
      # Credibility Intervals:
      # shape intercept
      s_tmp$a_w_quant_Like=quantile(zLik[,1],probs=c(0.05,0.5,0.95))
      # scale intercept
      s_tmp$a_C_quant_Like=quantile(zLik[,3],probs=c(0.05,0.5,0.95))
    }

  
  
    ################
    if (flag_bayes){
    ################
      # MCMC metropolis algorithm (maximum posterior method, BAYES THEOREM):
      print("Max Posterior")
      zPost_ALL=zPost_p_ALL=c()
      
      # ADAPTIVE METROPOLIS ALGORITHM
        it=0
        converg=F
        Nmcmc_tmp =Nmcmc/nburn*5 #accounting for thinning each 5 mcmc
        while ((converg ==F) & (Nmcmc_tmp < Nmcmc_max)){
          it=it+1
          set.seed(1)
          # metropolis algortyhm MCMC
          initval=c(priors_id$shape_intercept[[4]],priors_id$scale_intercept[[4]])
          mcmc_Post=adaptMCMC::MCMC(logPost_stat,
                                    n               = Nmcmc_tmp,
                                    init            = initval,
                                    scale           = c(scales[1], scales[3]),
                                    adapt           = T,
                                    showProgressBar = F,
                                    acc.rate        = acc.rate,
                                    priors          = priors_id,
                                    y               = y,
                                    t               = t,
                                    thr             = thr)
          # zPost_ALL=rbind(zPost_ALL, mcmc_Post$batch)
          zPost_ALL=rbind(zPost_ALL, mcmc_Post$samples)
          zPost_p_ALL=rbind(zPost_p_ALL, mcmc_Post$log.p)
          zPost=zPost_ALL[(nburn*nrow(zPost_ALL)+1):nrow(zPost_ALL),]
          zPost_p=zPost_p_ALL[(nburn*length(zPost_p_ALL)+1):length(zPost_p_ALL)]
          zPost=apply(zPost,2,function(x) x[seq(1,length(x),Nslim)])
          zPost_p=zPost_p[seq(1,length(zPost_p),Nslim)]
          # print(paste0("--------- accept.rate = ",floor(mcmc_Post$accept*100)," %"))
          converg = mcmc.diagn_stat(zPost       = zPost,
                                    priors_id   = priors_id,
                                    dir_res     = dir.res_id,
                                    name_output = paste0('post_',it),
                                    save        = T)
          s_tmp$converg_Post=converg
          converg=T
          if (converg==F){
            print(paste0("------------------ bad converg.: increasing mcmc"))
            Nmcmc_tmp=30000
          }
        }
        #s_tmp <- list(acc.rate=acc.rate)
        s_tmp$acc.rate = acc.rate
        s_tmp$mcmc_Post = mcmc_Post
        s_tmp$zPost= zPost
        s_tmp$zPost_acceptance = floor(mcmc_Post$accept*100)
    }
      
    
    #fminsearch (MLE)
    #################
    # Just for comparison, compute also fminsearch() as Marra et al., 2019:
    fmins=NULL
    initval=c(priors_id$shape_intercept[[4]],
              priors_id$scale_intercept[[4]])
    fmins=pracma::fminsearch(f=logLike1_stat, 
                             x0=initval,
                             y=y, 
                             t=t, 
                             thr=thr, 
                             maxiter=5000)
    s_tmp$theta_maxlik_fmins=fmins$xmin
    s_tmp$maxlik_fmins=fmins$fmin
    
    
    
  #^****************************************************************************
  } else if (algorithm== "rstan"){
  #^****************************************************************************
      #fminsearch (MLE)
      #################
      print("Max Likelihood with minimisation function")
      # Just for comparison, compute also fminsearch() as Marra et al., 2019:
      fmins=NULL
      initval=c(priors_id$shape_intercept[[4]],
                priors_id$scale_intercept[[4]])
      fmins=pracma::fminsearch(f=logLike1_stat, 
                               x0=initval,
                               y=y, 
                               t=t, 
                               thr=thr, 
                               maxiter=5000)
      s_tmp$theta_maxlik_fmins=fmins$xmin
      s_tmp$maxlik_fmins=fmins$fmin
    
      
      # Bayesian
      print("Apply Rstan package, Hamiltonian mcmc...")
      smodel = paste0(
        "data {
        int<lower=0> N_uncen;
        vector[N_uncen] y_uncen; 
        vector[N_uncen] t_uncen;
        int<lower=0> N_cen;
        vector[N_cen] t_cen;
        real<upper=min(y_uncen)> thr;
      }
      parameters {
        real<lower=0> alpha;
        real<lower=0> sigma;
      }
      model {
        // priors on component parameters
        alpha ~ ",  priors_id$shape_intercept[[1]], "(", priors_id$shape_intercept[[2]] , ",", priors_id$shape_intercept[[3]], ");  
        sigma ~ ",  priors_id$scale_intercept[[1]], "(", priors_id$scale_intercept[[2]] , ",", priors_id$scale_intercept[[3]], "); 
        
        target += weibull_lpdf(y_uncen | alpha, sigma);
        target += N_cen * weibull_lcdf(thr | alpha, sigma);
      }
      ")
      #       generated quantities{
      #         int<lower=0> k;
      #         vector[N] log_lik;
      #         k=0;
      #         for(i in 1:N_uncen){
      #           k=k+1;
      #           log_lik[k] = weibull_lpdf(y_uncen[i] | alpha, sigma);
      #         }
      #         for(i in 1:N_cen){
      #           k=k+1;
      #           log_lik[k] = weibull_lcdf(thr | alpha, sigma);
      #         }
      # 		  }"
      #      )
      file<-file.path(cmdstan_path(),"examples","smev")
      f<-write_stan_file(smodel,
                         dir=file,
                         basename=paste0("stat_smev2_dur",durations[dur]/60,"h.stan"))
      mod=0
      mod<-cmdstan_model(f)
      # mod$compile(force_recompile = T)
      # mod$print()
      # mod$exe_file()

      ord_data <- list(
        N_cen=length(y[y<thr]),
        t_cen=t[y<thr],
        N_uncen=length(y[y>=thr]),
        y_uncen=y[y>=thr],
        t_uncen=t[y>=thr],
        thr=thr,
        N=length(y)
      )
      init_f_stat<-function () list(alpha=rlnorm(1,priors_id$shape_intercept[[2]],0.1), 
                                    sigma=rlnorm(1,priors_id$scale_intercept[[2]],0.1))
      rstan_options(auto_write=T)
      # options(mc.cores=parallel::detectCores())
      Nmcmc_tmp =Nmcmc 
      fit2=0
      fit2<-rstan::stan(
        file=paste0(file, "/stat_smev2_dur",durations[dur]/60,"h.stan"), # Stan program
        data=ord_data, # named list of data
        chains=4, # num of Markov chains
        warmup=500, # num of warmup iterat per chain
        iter=1500, # tot num of iterat per chain
        cores=8, # number of cores (could use one per chain)
        refresh=0, # no progress shown
        seed=1,
        thin=1, # thinning to reduce mcmc autocorrelation
        init=init_f_stat,
        algorithm="NUTS") # "Fixed_param", "NUTS", "HMC"

      # alpha --> shape
      # sigma --> scale
      rstan::traceplot(fit2,pars=c("alpha","sigma"),inc_warmup=F,nrow=2)
      # scatterplots with correlations
      pairs(fit2,pars=c("alpha","sigma"),las=2)
      list_of_draws<-rstan::extract(fit2)
    
      
      # organise posterior sample:
      zPost=cbind(list_of_draws$alpha,list_of_draws$sigma)
      zPost_p=list_of_draws$lp__
      zPost=apply(zPost,2,function(x) x[seq(1,length(x),1)]) # slimming
      zPost_p=zPost_p[seq(1,length(zPost_p),1)]
      s_tmp$zPost=zPost
      s_tmp$zPost_p=zPost_p
      converg=mcmc.diagn_stat(zPost=zPost,
                              priors_id=priors_id,
                              dir_res=dir.res_id,
                              name_output=paste0('post_1'),
                              save=T)
      s_tmp$converg_Post=converg
      
      # compute posterior likelihood for the performance criteria:
      marg.like=c()
      for (tt in 1:nrow(zPost)){
        marg.like[tt]=logLike_stat(theta=zPost[tt,],y=y,t=t,thr=thr)
      }
      
      # compute maxpost:
      s_tmp$theta_maxpost=zPost[which.max(zPost_p),]
      s_tmp$maxpost=max(zPost_p)
      
      # Criterion using posterior likelihood or posterior (different versions):
      # DIC: Gelman et al, 2004 "Bayesian data analysis"
      all_dic=calc.dic(zPos=zPost,
                       lik=marg.like,
                       post=zPost_p,
                       lik.fun=logLike_stat,
                       t=t,
                       y=y,
                       thr=thr)   
      s_tmp$DIC1  =all_dic$DIC1
      s_tmp$DIC1_b=all_dic$DIC1_b
      s_tmp$DIC2  =all_dic$DIC2
      s_tmp$DIC2_b=all_dic$DIC2_b
      s_tmp$DIC3  =all_dic$DIC3
      
      # AIC (Aikake, 1974):
      all_aic=calc.aic(zPost=zPost,lik=marg.like,post=zPost_p)
      s_tmp$AIC=all_aic$AIC
      s_tmp$AIC_b=all_aic$AIC_b
      
      # BIC (Schwarz, 1978):
      all_bic=calc.bic(zPost=zPost,lik=marg.like,post=zPost_p, y=y)
      s_tmp$BIC=all_bic$BIC
      s_tmp$BIC_b=all_bic$BIC_b
      
      # HQC (Hannan-Quinn, 1979):
      all_hqc=calc.hqc(zPost=zPost,lik=marg.like,post=zPost_p, y=y)
      s_tmp$HQC=all_hqc$HQC
      s_tmp$HQC_b=all_hqc$HQC_b
      
      # Credibility Intervals of posterior sample:
      # shape intercept
      s_tmp$a_w_quant_Post=quantile(zPost[,1],probs=c(0.05,0.5,0.95))
      # scale intercept
      s_tmp$a_C_quant_Post=quantile(zPost[,2],probs=c(0.05,0.5,0.95))
      
      
      
      # TO DO:
      #^********
      # PSIS-LOO computation:  Leave-one-out cross-validation (LOO) 
      # - "elpd_loo" (expected log-predictive density)
      # - p_loo (effective number of parameters)
      # - looic = − 2*elpd_loo  (LOO information criterion as "deviance").
      # the last line of output gives k, the reliability of the LOO approx,
      # (see Vehtari, Simpson, Gelman, Yao, and Gabry (2019)). 
      
      # if log_lik saved in cmdstanfit object through "generated quantities":
      # all_ELPD1=loo::loo(fit$draws("log_lik"),
      #                   r_eff=relative_eff(fit$draws("log_lik")),
      #                   merge_chains=T,
      #                   save_psis = T,
      #                   cores = getOption("mc.cores", 4))
      
      # # otherwise we compute manually the log-likelihood for each observation:
      # log_lik_mat=0
      # log_lik_mat <- sapply(1:nrow(zPost), function(i) logLike_stat_loo(thetat=zPost[i,], y, t, thr))
      # # s_tmp$log_lik_mat=log_lik_mat
      # 
      # all_ELPD = loo::loo(t(log_lik_mat),
      #                     save_psis = T,
      #                     cores = getOption("mc.cores", 4))
      # s_tmp$all_ELPD = all_ELPD$estimates
      # 
      # # WAIC information criterion:
      # # WAIC=loo::waic(fit$draws("log_lik"))
      # WAIC=loo::waic(t(log_lik_mat))
      # s_tmp$WAIC= WAIC$estimates
      # 
      # # marginal likelihood for Bayes Factor:
      #########################################
      # .... to be implemented yet !!!
      
      
  #^****************************************************************************
  } else if (algorithm =="cmdstan"){
  #^****************************************************************************
    
    
    #fminsearch (MLE)
    #################
    print("Max Likelihood with minimisation function")
    # Just for comparison, compute also fminsearch() as Marra et al., 2019:
    fmins=NULL
    initval=c(priors_id$shape_intercept[[4]],
              priors_id$scale_intercept[[4]])
    fmins=pracma::fminsearch(f=logLike1_stat, 
                             x0=initval,
                             y=y, 
                             t=t, 
                             thr=thr, 
                             maxiter=5000)
    s_tmp$theta_maxlik_fmins=fmins$xmin
    s_tmp$maxlik_fmins=fmins$fmin
    
    
    print("Apply cmdstanr package, running hamiltonian mcmc...")
    
    if (flag_like){
      # perform also maximum LIKELIHOOD to be compared with posterior sample:
      print("Max Likelihood")
      converg=F
      # NOT DONE YET!!!
      #
    }
    
    if (flag_bayes){
      
      smodel = paste0(
      "data {
        int<lower=0> N_uncen;
        vector[N_uncen] y_uncen; 
        vector[N_uncen] t_uncen;
        int<lower=0> N_cen;
        vector[N_cen] t_cen;
        real <upper=min(y_uncen)> thr;
      }
      parameters {
        real<lower=0> alpha;
        real<lower=0> sigma;
      }
      model {
        // priors on component parameters
        alpha ~ ",priors_id$shape_intercept[[1]],"(",priors_id$shape_intercept[[2]],",",priors_id$shape_intercept[[3]],");  
        sigma ~ ",priors_id$scale_intercept[[1]],"(",priors_id$scale_intercept[[2]],",",priors_id$scale_intercept[[3]],"); 
        
        target += weibull_lpdf(y_uncen | alpha, sigma);
        target += N_cen * weibull_lcdf(thr | alpha, sigma);
      }
      ")
      #       generated quantities{
      #         int<lower=0> k;
      #         vector[N] log_lik;
      #         k=0;
      #         for(i in 1:N_uncen){
      #           k=k+1;
      #           log_lik[k] = weibull_lpdf(y_uncen[i] | alpha, sigma);
      #         }
      #         for(i in 1:N_cen){
      #           k=k+1;
      #           log_lik[k] = weibull_lcdf(thr | alpha, sigma);
      #         }
      # 		  }"
      #      )
      
      # y_uncen ~ weibull(alpha, sigma);
      # real <upper=min(y_uncen)> thr;
      file<-file.path(cmdstan_path(),"examples","smev")
      f<-write_stan_file(smodel,
                         dir=file,
                         basename=paste0("stat_smev2_dur",durations[dur]/60,"h.stan"))
      mod=0
      mod<-cmdstan_model(f)
      # mod$compile(force_recompile = T)
      # mod$print()
      # mod$exe_file()
      ord_data<-list(
        N_cen=length(y[y<thr]),
        # y_cen=y[y<thr],
        t_cen=t[y<thr],
        N_uncen=length(y[y>=thr]),
        y_uncen=y[y>=thr],
        t_uncen=t[y>=thr],
        thr=thr,
        N=length(y)
      )
      init_f_stat<-function () list(alpha=rlnorm(1,priors_id$shape_intercept[[2]],0.1), 
                                    sigma=rlnorm(1,priors_id$scale_intercept[[2]],0.1))
      # init_f_stat<-list(list(alpha=0.7,sigma=1)
      fit=0
      nchain=0
      Nmcmc_tmp=Nmcmc
      while (nchain<4){
        print(paste0("generate ",Nmcmc_tmp/4," sim. per chain (4 chains), burned first half, tot= ",Nmcmc_tmp))
        fit <- mod$sample(
          data=ord_data,
          seed=123,
          chains=4,
          parallel_chains=4, # getOption("mc.cores",8),
          iter_warmup=500, # Nmcmc_tmp/4,
          iter_sampling=Nmcmc_tmp/4,
          show_messages=F,
          show_exceptions=F,
          init=init_f_stat,
          save_warmup=F,
          refresh=0 # print update every 500 iter.
        )
        # fit$summary(variables = c("alpha", "sigma"), "mean", "sd")
        # fit$diagnostic_summary()
        draws<-fit$draws()
        list_of_draws=as_draws_df(draws)
        nchain=length(unique(list_of_draws$.chain))
        if (nchain<4){
          message("one or more chain did not converge !")
        }      
      }

      s_tmp$init_values_post = fit$init()
      # s_tmp$draws= draws
      
      # traceplots:
      color_scheme_set("mix-blue-red")
      gg_traceplots=bayesplot::mcmc_trace(
        draws,
        pars=c("alpha","sigma"),
        n_warmup=0,
        facet_args=list(nrow=2,labeller=label_parsed,ncol=1,strip.position="left"))+
        facet_text(size=15)+theme_bw(base_size = 12)+
        theme(panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank())
      ggsave(gg_traceplots,filename=paste0(dir.res_id,"/trace_chains_mcmc_post_1.png"),
             width=10,height=4,dpi=80)
      
      # intervals param:
      gg_areas=mcmc_areas(draws,
                          point_est=c("mean"),
                          border_size=0.5,
                          area_method="equal height",
                          pars=c("alpha","sigma"),
                          prob=0.95)+theme_bw(base_size = 20)
      ggsave(gg_areas,filename=paste0(dir.res_id,"/areas_param_mcmc_post_1.png"),
             width=6,height=4,dpi=80)
      
      # pdf param:
      gg_pdf=mcmc_dens_overlay(draws, 
                               pars=c("alpha","sigma"))+
        facet_text(size=15)+theme_bw(base_size=20)
      ggsave(gg_pdf, filename=paste0(dir.res_id,"/pdf_param_mcmc_chains_post_1.png"),
             width=10,height=4,dpi=80)
      

      #save stan fit object:
      # s_tmp$fit = fit
    }    
    
    # organise posterior sample:
    zPost=cbind(list_of_draws$alpha,list_of_draws$sigma)
    zPost_p=list_of_draws$lp__
    zPost=apply(zPost,2,function(x) x[seq(1,length(x),1)]) # slimming
    zPost_p=zPost_p[seq(1,length(zPost_p),1)]
    s_tmp$zPost=zPost
    s_tmp$zPost_p=zPost_p
    converg=mcmc.diagn_stat(zPost=zPost,
                            priors_id=priors_id,
                            dir_res=dir.res_id,
                            name_output=paste0('post_1'),
                            save=T)
    s_tmp$converg_Post=converg
    
    # compute posterior likelihood for the performance criteria:
    marg.like=c()
    for (tt in 1:nrow(zPost)){
      marg.like[tt]=logLike_stat(theta=zPost[tt,],y=y,t=t,thr=thr)
    }
    
    # compute maxpost:
    s_tmp$theta_maxpost=zPost[which.max(zPost_p),]
    s_tmp$maxpost=max(zPost_p)
    
    # Criterion using posterior likelihood or posterior (different versions):
    # DIC: Gelman et al, 2004 "Bayesian data analysis"
    all_dic=calc.dic(zPos=zPost,
                     lik=marg.like,
                     post=zPost_p,
                     lik.fun=logLike_stat,
                     t=t,
                     y=y,
                     thr=thr)   
    s_tmp$DIC1  =all_dic$DIC1
    s_tmp$DIC1_b=all_dic$DIC1_b
    s_tmp$DIC2  =all_dic$DIC2
    s_tmp$DIC2_b=all_dic$DIC2_b
    s_tmp$DIC3  =all_dic$DIC3
    
    # AIC (Aikake, 1974):
    all_aic=calc.aic(zPost=zPost,lik=marg.like,post=zPost_p)
    s_tmp$AIC=all_aic$AIC
    s_tmp$AIC_b=all_aic$AIC_b
    
    # BIC (Schwarz, 1978):
    all_bic=calc.bic(zPost=zPost,lik=marg.like,post=zPost_p, y=y)
    s_tmp$BIC=all_bic$BIC
    s_tmp$BIC_b=all_bic$BIC_b
    
    # HQC (Hannan-Quinn, 1979):
    all_hqc=calc.hqc(zPost=zPost,lik=marg.like,post=zPost_p, y=y)
    s_tmp$HQC=all_hqc$HQC
    s_tmp$HQC_b=all_hqc$HQC_b
    
    # Credibility Intervals of posterior sample:
    # shape intercept
    s_tmp$a_w_quant_Post=quantile(zPost[,1],probs=c(0.05,0.5,0.95))
    # scale intercept
    s_tmp$a_C_quant_Post=quantile(zPost[,2],probs=c(0.05,0.5,0.95))
    
    
    
    # TO DO:
    #^********
    # PSIS-LOO computation:  Leave-one-out cross-validation (LOO) 
    # - "elpd_loo" (expected log-predictive density)
    # - p_loo (effective number of parameters)
    # - looic = − 2*elpd_loo  (LOO information criterion as "deviance").
    # the last line of output gives k, the reliability of the LOO approx,
    # (see Vehtari, Simpson, Gelman, Yao, and Gabry (2019)). 
    
    # if log_lik saved in cmdstanfit object through "generated quantities":
    # all_ELPD1=loo::loo(fit$draws("log_lik"),
    #                   r_eff=relative_eff(fit$draws("log_lik")),
    #                   merge_chains=T,
    #                   save_psis = T,
    #                   cores = getOption("mc.cores", 4))
    
    # # otherwise we compute manually the log-likelihood for each observation:
    # log_lik_mat=0
    # log_lik_mat <- sapply(1:nrow(zPost), function(i) logLike_stat_loo(thetat=zPost[i,], y, t, thr))
    # # s_tmp$log_lik_mat=log_lik_mat
    # 
    # all_ELPD = loo::loo(t(log_lik_mat),
    #                     save_psis = T,
    #                     cores = getOption("mc.cores", 4))
    # s_tmp$all_ELPD = all_ELPD$estimates
    # 
    # # WAIC information criterion:
    # # WAIC=loo::waic(fit$draws("log_lik"))
    # WAIC=loo::waic(t(log_lik_mat))
    # s_tmp$WAIC= WAIC$estimates
    # 
    # # marginal likelihood for Bayes Factor:
    #########################################
    # .... to be implemented yet !!!
  }
  
  return(s_tmp)
}

    
    
