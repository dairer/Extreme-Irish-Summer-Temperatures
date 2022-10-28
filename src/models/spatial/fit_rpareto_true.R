rm(list=ls())
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland//")
library(tidyverse)
library(doParallel)
library(foreach)


marg_mod = "mod_1"
#marg_mod = "mod_2"
#marg_mod = "mod_4"




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
  inversePsi = psi %>% solve()  # inverse of "covarianve" matrix
  
  detPsi = determinant(psi)$modulus[1] # calculates log det by default
  
  omegas = log(t(sapply(dt, function(x) x[-conditioned.site]))/sapply(dt, "[[", conditioned.site)) %>% sweep(2, conditioned.site.lambda[-1], "+")
  # for some reason when we have 2 sites the matrix is transposed...
  
  if(nlocs == 2){
    omegas = t(omegas)
  }
  
  summ = apply(omegas, MARGIN = 1, function(x){t(x) %*% inversePsi %*% (x)}) %>% sum
  length(dt)*detPsi+ summ # LL
}

# negative likelihood function
ngll = function(pars){
  
  print(pars)
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
data_for_rpareto = readRDS(paste0("data/processed/data_for_rpareto/true/data_for_rpareto_", marg_mod))
.GlobalEnv$exceedances <- data_for_rpareto$exceedances
.GlobalEnv$exceedances_locs <- data_for_rpareto$exceedances_locs
.GlobalEnv$HRD_ll <- HRD_ll


# new_thres = quantile(data_for_rpareto$extreme_dates$cost,0.25)
# cost_of_each_event = data_for_rpareto$exceedances %>% map(mean) %>% unlist
# .GlobalEnv$exceedances <- data_for_rpareto$exceedances[cost_of_each_event > new_thres]
# .GlobalEnv$exceedances_locs <- data_for_rpareto$exceedances_locs[cost_of_each_event > new_thres]
# 

.GlobalEnv$stat <- TRUE
cl <<- parallel::makeForkCluster(25)
doParallel::registerDoParallel(cl)
fit <<- optimx::optimx(par = c(1,1), fn = ngll ,hessian = F, method = 'Nelder-Mead') # Variance + range constant
#14.49639, 778.83897
parallel::stopCluster(cl)

# method = 'L-BFGS-B'
tibble(fit$p1, fit$p2) %>%
  write_csv(paste0("output/rpareto_model_fits/true_rpareto_fits_model_", marg_mod, ".csv"), append = T)

# tibble(fit$value) %>%
#   write_csv(paste0("corrections/data/processed_data/data_for_r_pareto/true/stationary_LL_value_",marg_mod,"_Nelder-Mead_frechet.csv"), append = T)
# 
# 
# # 
variogram_model_score <- function(h, vario_params){
  alpha <-  (vario_params[1])
  beta <- (vario_params[2])
  nu = 0.2
  res = rep(NA, length(h))
  res[h == 0] = 0
  res[h != 0] = alpha*(1 - ((((h[h != 0] /beta)^nu) * besselK(x = (h[h != 0] / beta), nu = nu))/(gamma(nu)*2^(nu - 1))))
  res
}


#        mod_4b = 2 * (1 - pnorm(sqrt(variogram_model_score(distances, c( 6.777471273861808,84.70033942987607)) / 2)))) %>%

distances = seq(0.05,3.75, length.out = 100)
# ----- TRUE estimate from model
tibble(distances,
       #mod_4a = 2 * (1 - pnorm(sqrt(variogram_model_score(distances, c( 15.23729 ,705.03887)) / 2))),
       mod_4b = 2 * (1 - pnorm(sqrt(variogram_model_score(distances, c( 9.773973, 220.923958)) / 2)))) %>%
  ggplot()+
 # geom_line(aes(distances, mod_4a, col = 'mod_4a'))+
  geom_line(aes(distances, mod_4b, col = 'mod_4b'))
# 
# # 
# 3.089766425860744,11.584922449102129
# 
# # distances = seq(0.05,4.5, length.out = 100)
# # # ----- TRUE estimate from model
# # tibble(distances,
# #        mod_4 = 2 * (1 - pnorm(sqrt(variogram_model_score(distances, c( 3.089766, 11.58477)) / 2))),
# #        mod_1 = 2 * (1 - pnorm(sqrt(variogram_model_score(distances, c( 17.45752, 991.03576)) / 2))),
# #        mod_1_bts = 2 * (1 - pnorm(sqrt(variogram_model_score(distances, c( 16.86171, 991.02767)) / 2)))) %>%
# #   ggplot()+
# #   geom_line(aes(distances, mod_4, col = 'mod_4'))+
# #   geom_line(aes(distances, mod_1_bts, col = 'mod_1_bts'))+
# #   geom_line(aes(distances, mod_1, col = 'mod_1'))
