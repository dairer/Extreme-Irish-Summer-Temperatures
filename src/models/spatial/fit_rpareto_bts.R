rm(list=ls())
library(tidyverse)
library(doParallel)
library(foreach)
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")

# read_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_", "mod_1", ".csv-Mead.csv"),
#          col_names = c('bts', 'alpha', 'beta')) 

fit_bts_mod = function(marg_mod, bts_seq, num_quantiles, change_to_90 = FALSE){
  
  # bts_fitted = read_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_", marg_mod, "_L-BFGS-B.csv"),
  #                       col_names = c('bts', 'alpha', 'beta')) %>%
  #   pull(bts)
  # 
  # bts_seq = bts_seq[!(bts_seq %in% bts_fitted)] # dont refit!
  
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
      # for some reason when we have 2 sites the matrix is transposed...
      
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
    
    #print(pars)
    if( pars[2] < 0) return (2^100)
    if( pars[1] < 0) return (2^100)
    # if( pars[1] <= 10) return (2^100)
    
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

      # if(change_to_90){
      #   # find 50% of data
      #   # remove obs below this
      #   
      #   new_thres = quantile(data_for_rpareto$extreme_dates$cost,0.25)
      #   cost_of_each_event = data_for_rpareto$exceedances %>% map(mean) %>% unlist
      #   .GlobalEnv$exceedances <- data_for_rpareto$exceedances[cost_of_each_event > new_thres]
      #   .GlobalEnv$exceedances_locs <- data_for_rpareto$exceedances_locs[cost_of_each_event > new_thres]
      # }else{
        .GlobalEnv$exceedances <- data_for_rpareto$exceedances[obs_to_rem]
        .GlobalEnv$exceedances_locs <- data_for_rpareto$exceedances_locs[obs_to_rem]
      #}
      
      .GlobalEnv$HRD_ll <- HRD_ll
      .GlobalEnv$stat <- TRUE

      tryCatch({
        set.seed(12341234)
        cl <- parallel::makeForkCluster(12)
        doParallel::registerDoParallel(cl)
        
        fit <- optimx::optimx(par = init_vals, fn = ngll ,hessian = F, method = 'L-BFGS-B') # Variance + range constant
        #fit <- optimx::optimx(par = init_vals, fn = ngll ,hessian = F, method = 'Nelder-Mead')#, method = 'L-BFGS-B') # Variance + range constant
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

# job::job({fit_bts_mod(marg_mod = "mod_1", bts_seq = seq(1,100), num_quantiles = 15)})
# job::job({fit_bts_mod(marg_mod = "mod_1", bts_seq = seq(101,200), num_quantiles = 15)})
# job::job({fit_bts_mod(marg_mod = "mod_1", bts_seq = seq(201,300), num_quantiles = 15)})
job::job({fit_bts_mod(marg_mod = "mod_1", bts_seq = seq(301,400), num_quantiles = 15)})
job::job({fit_bts_mod(marg_mod = "mod_1", bts_seq = seq(401,500), num_quantiles = 15)})

# job::job({fit_bts_mod(marg_mod = "mod_2", bts_seq = seq(1,100), num_quantiles = 15)})
# job::job({fit_bts_mod(marg_mod = "mod_2", bts_seq = seq(101,200), num_quantiles = 15)})
# job::job({fit_bts_mod(marg_mod = "mod_2", bts_seq = seq(201,300), num_quantiles = 15)})
job::job({fit_bts_mod(marg_mod = "mod_2", bts_seq = seq(301,400), num_quantiles = 15)})
job::job({fit_bts_mod(marg_mod = "mod_2", bts_seq = seq(401,500), num_quantiles = 15)})





# ----- 11 october
job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(1,100), num_quantiles = 25)})
job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(101,200), num_quantiles = 25)})


# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(1,65), num_quantiles = 25)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(66,120), num_quantiles = 25)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(121,190), num_quantiles = 25)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(191,250), num_quantiles = 25)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(251,280), num_quantiles = 25)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(281,300), num_quantiles = 25)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(301,330), num_quantiles = 25)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(331,360), num_quantiles = 25)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(361,390), num_quantiles = 25)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(391,420), num_quantiles = 25)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(421,450), num_quantiles = 25)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(451,455), num_quantiles = 25)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(456,460), num_quantiles = 25)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(461,470), num_quantiles = 25)})
# job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(471,480), num_quantiles = 25)})




# job::job({fit_bts_mod(marg_mod = "mod_2", bts_seq = seq(1,40), num_quantiles = 15)})
# job::job({fit_bts_mod(marg_mod = "mod_2", bts_seq = seq(41,80), num_quantiles = 15)})
# job::job({fit_bts_mod(marg_mod = "mod_2", bts_seq = seq(81,120), num_quantiles = 15)})
# job::job({fit_bts_mod(marg_mod = "mod_2", bts_seq = seq(121,160), num_quantiles = 15)})
# job::job({fit_bts_mod(marg_mod = "mod_2", bts_seq = seq(161,200), num_quantiles = 15)})
# job::job({fit_bts_mod(marg_mod = "mod_2", bts_seq = seq(201,250), num_quantiles = 15)})
# job::job({fit_bts_mod(marg_mod = "mod_2", bts_seq = seq(251,300), num_quantiles = 15)})
# # job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(161,180), num_quantiles = 15)})





#job::job({fit_bts_mod(marg_mod = "mod_4", bts_seq = seq(181,200), num_quantiles = 15)})
# 
variogram_model_score <- function(h, vario_params){
  alpha <-  (vario_params[1])
  beta <- (vario_params[2])
  nu = 0.2
  res = rep(NA, length(h))
  res[h == 0] = 0
  res[h != 0] = alpha*(1 - ((((h[h != 0] /beta)^nu) * besselK(x = (h[h != 0] / beta), nu = nu))/(gamma(nu)*2^(nu - 1))))
  res
}
# 
# 
# 
# 
# res = read_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_", 'mod_4', ".csv-Mead.csv"),
#          col_names = c('bts', 'alpha', 'beta'))
# 
# 
# res$alpha %>% hist
# res$beta %>% hist
distances = seq(0.005, 3.5, length.out = 100)
# 
# my_mat = c()
# for(i in seq(nrow(res))){
#   my_mat = rbind(my_mat, 2 * (1 - pnorm(sqrt(variogram_model_score(distances, c(res[i,]$alpha, res[i,]$beta)) / 2))))
# }
# 
# tibble(distances,
#        true =  2 * (1 - pnorm(sqrt(variogram_model_score(distances, c(9.773973, 220.926720)) / 2))),
#        lower = my_mat %>% apply(MARGIN = 2, FUN = quantile, 0.025),
#        upper = my_mat %>% apply(MARGIN = 2, FUN = quantile, 0.975)) %>%
#   ggplot()+
#   geom_line(aes(distances, true))+
#   geom_ribbon(aes(distances, ymin = lower, ymax = upper), alpha = 0.25)+
#   geom_line(aes(distances, upper), alpha = 0.25)+
#   geom_line(aes(distances, lower), alpha = 0.25)
# 
# 
# 
# 
# res$bts %>% sort
# 
# 
# # ----- TRUE estimate from model
tibble(distances,
       mod_4a = 2 * (1 - pnorm(sqrt(variogram_model_score(distances, c( 12,450)) / 2))),
       mod_4b = 2 * (1 - pnorm(sqrt(variogram_model_score(distances, c( 12,450)) / 2)))) %>%
  ggplot()+
  geom_line(aes(distances, mod_4a, col = 'mod_4a'))+
  geom_line(aes(distances, mod_4b, col = 'mod_4b'))
# 
# 
# 
# 
# 
# 
# init_vals = read_csv(paste0("output/rpareto_model_fits/true_rpareto_fits_model_",'mod_4', ".csv"),
#                      col_names = F) %>% 
#   .[nrow(.),] %>% # get last row in this file
#   unlist %>% as.numeric()
# 
# read_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_", 'mod_4', ".csv-Mead.csv"),
#          col_names = c('bts', 'alpha', 'beta')) %>%
#   ggplot()+
#   geom_histogram(aes(alpha))+
#   geom_vline(aes(xintercept = init_vals[1]))
# 
# 
# read_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_", 'mod_4', ".csv-Mead.csv"),
#          col_names = c('bts', 'alpha', 'beta')) %>%
#   filter(beta>100) %>%
#   ggplot()+
#   geom_density(aes(beta))+
#   geom_vline(aes(xintercept = init_vals[2]))+
#   xlim(c(300, 1600))
# 
# init_vals = read_csv(paste0("output/rpareto_model_fits/true_rpareto_fits_model_",'mod_1', ".csv"),
#                      col_names = F) %>% 
#   .[nrow(.),] %>% # get last row in this file
#   unlist %>% as.numeric()
# 
# 
# 
# 
# read_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_mod_1.csv-Mead_thresh_90.csv"),
#          col_names = c('bts', 'alpha', 'beta')) %>%
#   ggplot()+
#   geom_histogram(aes(alpha))+
#   geom_vline(aes(xintercept = init_vals[1]))
# 
# read_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_mod_1.csv-Mead_thresh_90.csv"),
#          col_names = c('bts', 'alpha', 'beta')) %>%
#   ggplot()+
#   geom_histogram(aes(beta))+
#   geom_vline(aes(xintercept = init_vals[2]))
# 
# 
# 
# read_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_mod_1.csv-Mead_thresh_90.csv"),
#          col_names = c('bts', 'alpha', 'beta')) %>%
#   pull(alpha) %>%
#   quantile(c(0.025, 0.975))
# 
# 
# 
# 
# 
# 
# 
# 
# init_vals = read_csv(paste0("output/rpareto_model_fits/true_rpareto_fits_model_",'mod_1', ".csv"),
#                      col_names = F) %>% 
#   .[nrow(.),] %>% # get last row in this file
#   unlist %>% as.numeric()
# 
# 
# 
# 
# read_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_", 'mod_1', ".csv-Mead.csv"),
#          col_names = c('bts', 'alpha', 'beta')) %>%
#   pull(alpha) %>%
#   quantile(c(0.025, 0.975))
# 
# read_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_", 'mod_1', ".csv-Mead.csv"),
#          col_names = c('bts', 'alpha', 'beta')) %>%
#   ggplot()+
#   geom_histogram(aes(alpha))+
#   geom_vline(aes(xintercept = init_vals[1]))
# 
# 
# 
# read_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_", 'mod_1', ".csv-Mead.csv"),
#          col_names = c('bts', 'alpha', 'beta')) %>%
#   ggplot()+
#   geom_histogram(aes(beta))+
#   geom_vline(aes(xintercept = init_vals[2]))
