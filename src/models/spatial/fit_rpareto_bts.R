rm(list=ls())
library(tidyverse)
library(doParallel)
library(foreach)
setwd("~/Extreme-Irish-Summer-Temperatures/")


# code here is adapted from R package mvPot: Multivariate Peaks-over-Threshold Modelling for Spatial Extreme Events, author: de Fondeville
# https://cran.r-project.org/web/packages/mvPot/index.html

fit_bts_mod = function(marg_mod, bts_seq, num_quantiles, change_to_90 = FALSE){
  init_vals = read_csv(paste0("output/rpareto_model_fits/true_rpareto_fits_model_",marg_mod, ".csv"),
                       col_names = F) %>% 
    .[nrow(.),] %>% # get last row in this file
    unlist %>% as.numeric()
  
  HRD_ll = function(dt, lcs, vr, conditioned.site = 1, zeta, gamma, eta){
    nlocs = nrow(lcs) # number of sites
    loc.id.pairs = expand.grid(seq(nlocs) ,seq(nlocs))
    loc.pairs = cbind(lcs[loc.id.pairs[,1],], lcs[loc.id.pairs[,2],]) # all pairs of locations
    
    dstncs = loc.pairs %>%
      apply(MARGIN = 1, FUN = function(x) sqrt((x[3] - x[1])^2 + (x[4] - x[2])^2))

    Lambda = (vr(dstncs)) %>% matrix(nrow = nlocs, byrow = T)
    conditioned.site.lambda = Lambda[conditioned.site,]
    Lambda = Lambda[-conditioned.site, -conditioned.site]
    
    psi = outer(conditioned.site.lambda[-1], conditioned.site.lambda[-1], "+") - Lambda # "covar" matrix

    tryCatch({
      inversePsi = psi %>% solve()  # inverse of "covarianve" matrix
      
      detPsi = determinant(psi)$modulus[1] # calculates log det by default
      
      omegas = log(t(sapply(dt, function(x) x[-conditioned.site]))/sapply(dt, "[[", conditioned.site)) %>% sweep(2, conditioned.site.lambda[-1], "+")

      if(nlocs == 2){
        omegas = t(omegas)
      }
      
      summ = apply(omegas, MARGIN = 1, function(x){t(x) %*% inversePsi %*% (x)}) %>% sum
      length(dt)*detPsi+ summ # LL
      
    }, error=function(e){
      print("omg caugth error")
      return (2^30)
    })
    
  }
  
  ngll = function(pars){

    if( pars[2] < 0) return (2^100)
    if( pars[1] < 0) return (2^100)

    LL <- foreach(i= 1:length(exceedances), .combine = 'c') %dopar% {
      my.vario <- function(h){
        if(stat){
          alpha = pars[1]
          beta = pars[2] 
        }else{
          alpha = pars[1] + pars[3]*temporal_covar[i]
          beta = pars[2] 
        }
        
        if(alpha <= 0){
          return (2^100)
        } 
        
        if( beta <= 0){
          return (2^100)
        } 
        
        nu = 0.2
        res = rep(NA, length(h))
        res[h == 0] = 0
        res[h != 0] = alpha*(1 - ((((h[h != 0] /beta)^nu) * besselK(x = (h[h != 0] / beta), nu = nu))/(gamma(nu)*2^(nu - 1))))
        res
      }
      HRD_ll(lcs = exceedances_locs[i][[1]], dt = exceedances[i], vr = my.vario, zeta = 1, gamma = 1, eta = 1)
    }
    LL = sum(LL)
    print(LL)
    LL
  }
  
  
  for(bs in bts_seq){
    print(bs)
    
    if(file.exists(paste0("data/processed/data_for_rpareto/bootstraps/data_for_rpareto_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_", bs))){
      
      data_for_rpareto = readRDS(paste0("data/processed/data_for_rpareto/bootstraps/data_for_rpareto_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_", bs))
      obs_to_rem = data_for_rpareto$exceedances %>% map(~{sum(.x<=0) <= 0}) %>% unlist  # remove obs with 0
      .GlobalEnv$exceedances <- data_for_rpareto$exceedances[obs_to_rem]
      .GlobalEnv$exceedances_locs <- data_for_rpareto$exceedances_locs[obs_to_rem]
      
      .GlobalEnv$HRD_ll <- HRD_ll
      .GlobalEnv$stat <- TRUE

      tryCatch({
        set.seed(12341234)
        cl <- parallel::makeForkCluster(12)
        doParallel::registerDoParallel(cl)
        
        fit <- optimx::optimx(par = init_vals, fn = ngll ,hessian = F, method = 'L-BFGS-B') # Variance + range constant
        parallel::stopCluster(cl)
        

        tibble(bs, fit$p1, fit$p2) %>%
          write_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_", marg_mod, "_L-BFGS-B.csv"), append = T)

      }, error=function(e){
        tibble(bs) %>%
          write_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_failed_model_", marg_mod,"_frechet.csv"), append = T)
      })
    }
  }
}

# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(1,65), num_quantiles = 30)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(66,120), num_quantiles = 30)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(121,190), num_quantiles = 30)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(191,300), num_quantiles = 30)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(301,280), num_quantiles = 30)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(281,300), num_quantiles = 30)})
