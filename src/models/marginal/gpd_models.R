
#  - - - - - - - - - MODEL 0
ngll_0 = function(par){
  scale_est = exp(par[1] + par[2]*clim_scale)
  
  shape_est =  par[3]  
  
  if(any(scale_est <= 0)) return(2^30)
  sum(evd::dgpd(x = excess_dat, loc=0, scale = scale_est, shape=shape_est, log=T)) #+ dnorm(shape_est, mean = 0, sd = 0.01, log = T))
}

fit_mod_0 = function(this_dat, this_clim_scale, initial_pars = c(0.35 , 0.6, -0.1520767)){
  excess_dat <<- this_dat
  clim_scale <<- log(this_clim_scale) 
  optim(par = initial_pars, fn = ngll_0, control = list(fnscale = -1))$par
}

ngll_0_fix_shape = function(par){
  scale_est = exp(par[1] + par[2]*clim_scale)
  if(any(scale_est <= 0)) return(-2^30)
  if(any(scale_est > -1/shape_est)) return(-2^30)
  if(any((1+shape_est*excess_dat/scale_est)< 0)) return(-2^30)
  sum(evd::dgpd(x = excess_dat, loc=0, scale = scale_est, shape=shape_est, log=T)) 
}

fit_mod_0_fix_shape  = function(this_dat, this_clim_scale, this_shape_est, initial_pars = c(0.35 , 0.6)){
  print(initial_pars)
  excess_dat <<- this_dat
  clim_scale <<- log(this_clim_scale) 
  shape_est <<- this_shape_est
  optim(par = c(0.35 , 0.6), fn = ngll_0_fix_shape, control = list(fnscale = -1))$par
}

my_predict_0 = function(estimates_pars, this_clim_scale){
  tibble(scale = exp(estimates_pars[1] + estimates_pars[2]*log(this_clim_scale)),
         shape = estimates_pars[3])
}

rl_mod_0 = function(estimates_pars, rl_quantile, thresh, this_clim_scale){
  estimated_scale = exp(estimates_pars[1] + estimates_pars[2]*log(this_clim_scale))
  estimated_shape = estimates_pars[3]
  return(thresh + estimated_scale * ((1-rl_quantile)^(-estimated_shape) - 1)/estimated_shape)
}

# - - - - - - - - - - MODEL 1
ngll_1 = function(par){
  if(par[2] < 0) return(2^30)
  scale_est = exp(par[1] + par[2]*clim_scale+ par[3]*loess_temp_anom)
  shape_est = par[4]
  
  if(any(scale_est <= 0)) return(2^30)
  if(any(scale_est > -1/shape_est)) return(2^30)
  if(any((1+shape_est*excess_dat/scale_est)< 0)) return(2^30)
  -sum(evd::dgpd(x = excess_dat, loc=0, scale = scale_est, shape=shape_est, log=T))
}

fit_mod_1 = function(this_dat, this_clim_scale, this_loess_temp_anom, initial_pars = c(0.3510713,  0.7598344, 0.3735851, -0.1429355)){
  excess_dat <<- this_dat
  clim_scale <<- log(this_clim_scale)
  loess_temp_anom <<- this_loess_temp_anom
  optim(par =initial_pars, fn = ngll_1)$par
}

ngll_1_fix_shape = function(par){
  if(par[2] < 0) return(2^30)
  scale_est = exp(par[1] + par[2]*clim_scale+ par[3]*loess_temp_anom)
  
  if(any(scale_est <= 0)) return(2^30)
  if(any(scale_est > -1/shape_est)) return(2^30)
  if(any((1+shape_est*excess_dat/scale_est)< 0)) return(2^30)
  -sum(evd::dgpd(x = excess_dat, loc=0, scale = scale_est, shape=shape_est, log=T))
}

fit_mod_1_fix_shape  = function(this_dat, this_clim_scale, this_loess_temp_anom,this_shape_est, initial_pars = c(0.3510713,  0.7598344, 0.3735851)){
  excess_dat <<- this_dat
  clim_scale <<- log(this_clim_scale) 
  loess_temp_anom <<- this_loess_temp_anom
  
  shape_est <<- this_shape_est
  optim(par = initial_pars, fn = ngll_1_fix_shape)$par
}

my_predict_1 = function(estimates_pars, this_clim_scale, this_loess_temp_anom){
  tibble(scale = exp(estimates_pars[1] + estimates_pars[2]*log(this_clim_scale)+ estimates_pars[3]*this_loess_temp_anom),
         shape = estimates_pars[4])
}

rl_mod_1 = function(estimates_pars, rl_quantile, thresh, this_clim_scale, this_loess_temp_anom){
  estimated_scale = exp(estimates_pars[1] + estimates_pars[2]*log(this_clim_scale)+ estimates_pars[3]*this_loess_temp_anom)
  estimated_shape = estimates_pars[4]
  return(thresh + estimated_scale * ((1-rl_quantile)^(-estimated_shape) - 1)/estimated_shape)
}


# - - - - - - - - - - MODEL 2
ngll_2 = function(par){
  scale_est = exp(par[1] + par[2]*clim_scale + par[3]*dist_sea + par[4]*loess_temp_anom + par[5]*dist_sea*loess_temp_anom) 
  shape_est = par[6]
  if(any(scale_est <= 0)) return(2^30)
  if(any(scale_est > -1/shape_est)) return(2^30)
  if(any((1+shape_est*excess_dat/scale_est)< 0)) return(2^30)
  -sum(evd::dgpd(x = excess_dat, loc=0, scale = scale_est, shape=shape_est, log=T))
}

fit_mod_2 = function(this_dat, this_clim_scale, this_loess_temp_anom, this_dist_sea, initial_pars = c( 0.0713, 0.700, 0.0513, 0.123, 0.00117, -0.15)){
  excess_dat <<- this_dat
  clim_scale <<- log(this_clim_scale)
  dist_sea <<- log(1+this_dist_sea)
  loess_temp_anom <<- this_loess_temp_anom
  
 optim(par = initial_pars, fn = ngll_2)$par  
}

my_predict_2 = function(estimates_pars, this_clim_scale, this_loess_temp_anom, this_dist_sea){
  tibble(scale = exp(estimates_pars[1] + estimates_pars[2]*log(this_clim_scale)+
                       estimates_pars[3]*log(1+this_dist_sea)+
                       estimates_pars[4]*this_loess_temp_anom+
                       estimates_pars[5]*log(1+this_dist_sea)*this_loess_temp_anom),
         shape = estimates_pars[6])
}

ngll_2_fix_shape = function(par){
  scale_est = exp(par[1] + par[2]*clim_scale + par[3]*dist_sea + par[4]*loess_temp_anom + par[5]*dist_sea*loess_temp_anom) 
  if(any(scale_est <= 0)) return(2^30)
  if(any(scale_est > -1/shape_est)) return(2^30)
  if(any((1+shape_est*excess_dat/scale_est)< 0)) return(2^30)
  -sum(evd::dgpd(x = excess_dat, loc=0, scale = scale_est, shape=shape_est, log=T))
}

fit_mod_2_fix_shape  = function(this_dat, this_clim_scale, this_loess_temp_anom, this_dist_sea, this_shape_est, initial_pars = c( 0.0713, 0.700, 0.0513, 0.123, 0.00117)){
  excess_dat <<- this_dat
  clim_scale <<- log(this_clim_scale)
  dist_sea <<- log(1+this_dist_sea)
  loess_temp_anom <<- this_loess_temp_anom
  shape_est <<- this_shape_est
  optim(par = initial_pars, fn = ngll_2_fix_shape)$par
}

rl_mod_2= function(estimates_pars, rl_quantile, thresh, this_clim_scale, this_loess_temp_anom, this_dist_sea){
  estimated_scale =exp(estimates_pars[1] + estimates_pars[2]*log(this_clim_scale)+
                         estimates_pars[3]*log(1+this_dist_sea)+
                         estimates_pars[4]*this_loess_temp_anom+
                         estimates_pars[5]*log(1+this_dist_sea)*this_loess_temp_anom)
  estimated_shape = estimates_pars[6]
  return(thresh + estimated_scale * ((1-rl_quantile)^(-estimated_shape) - 1)/estimated_shape)
  
}
