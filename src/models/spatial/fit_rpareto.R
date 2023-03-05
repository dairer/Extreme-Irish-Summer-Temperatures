# code here is adapted from R package mvPot: Multivariate Peaks-over-Threshold Modelling for Spatial Extreme Events, author: de Fondeville
# https://cran.r-project.org/web/packages/mvPot/index.html


rm(list=ls())
setwd("~/Extreme-Irish-Summer-Temperatures/")
library(tidyverse)
library(doParallel)
library(foreach)

marg_mod = "mod_0"
#marg_mod = "mod_1"
#marg_mod = "mod_2"

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

.GlobalEnv$stat <- TRUE
cl <<- parallel::makeForkCluster(25)
doParallel::registerDoParallel(cl)
fit <<- optimx::optimx(par = c(1,1), fn = ngll ,hessian = F, method = 'Nelder-Mead') # Variance + range constant
parallel::stopCluster(cl)

# method = 'L-BFGS-B'
tibble(fit$p1, fit$p2) %>%
  write_csv(paste0("output/rpareto_model_fits/true_rpareto_fits_model_", marg_mod, ".csv"), append = T)
