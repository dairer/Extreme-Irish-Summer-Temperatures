# # Author: Dáire Healy
# # Date: 17 Feb 2022
# 
# Fits and saves quantile regression model and lambda estimates
gc()
rm(list = ls())

library(tidyverse)
library(evgam)
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")

num_quantiles = 100
obs_data = readRDS(paste0("data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))





# ---- get covariates for prediction
temporal_covariates = obs_data %>%
  dplyr::select(year, loess_temp_anom, loess_glob_temp_anom) %>%
  unique() %>%
  arrange(year)

quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
fit_clim_quants = T # bool, re-estimate clim quantiles?

if(fit_clim_quants){
  file.remove(paste0("data/processed/qunatile_model_fit_pars_num_quantiles_",num_quantiles,".csv"))
  for(q in seq_along(quantiles_to_estimate_bulk)){
    #obs_data_for_quant_reg = obs_data
    obs_data_for_quant_reg = rlang::duplicate(obs_data, shallow = FALSE)
    zeta = obs_data_for_quant_reg$quantile[[1]][q] # quantile to estimate
    print(paste0("Fitting  tau = ", zeta))

    obs_data_for_quant_reg$value = obs_data_for_quant_reg$value %>% lapply(`[[`, q) %>% unlist
    qunatile_model_fit <- evgam(maxtp ~ value + loess_temp_anom, obs_data_for_quant_reg,
                                  family = "ald", ald.args = list(tau = zeta))
    
    # save parameter estimates
    tibble(tau = zeta,
           beta_0 = qunatile_model_fit$location$coefficients[1],
           beta_1 = qunatile_model_fit$location$coefficients[2],
           beta_2 = qunatile_model_fit$location$coefficients[3]) %>%
      write_csv(paste0("data/processed/qunatile_model_fit_pars_num_quantiles_",num_quantiles,".csv"), append = T)
  }
}

# --- read in fitted quantile regression coefficients
quant_reg_model_pars = read_csv(paste0("data/processed/qunatile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
         col_names = c('tau', 'beta_0', 'beta_1', 'beta_2'))

# # --- creates a tibble with each station and its quantule model
obs_smoothed_quantiles = obs_data %>%
  group_by(Station) %>%
  group_map(~{

    # --- get the climate qunatile estimates closest to current station
    clim_vals <<- obs_data %>%
      filter(Station == .x$Station[1]) %>%
      dplyr::select(quantile, value) %>%
      unique() %>%
      pull(value) %>%
      unlist()

    # predict quantile for each year and site
    quant_reg_pars = quant_reg_model_pars %>%
      arrange(tau)

    res = c()
    for(q in seq_along(quantiles_to_estimate_bulk)){
      qpars = quant_reg_pars[q,]
      res = rbind(res,
                  tibble(quantile =  qpars$tau,
                         year = temporal_covariates$year,
                         loess_temp_anom = temporal_covariates$loess_temp_anom,
                         quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2)*(temporal_covariates$loess_temp_anom)))
    }

    print(paste0("Interpolating quantile estimates for ", .x$Station[1]))
    # interpolate quantiles over tau for each year
    res %>%
      group_by(year) %>%
      group_map(~{
        tibble(year = .x$year[1],
               tau_to_temp = list(splinefun(.x$quantile,.x$quant_value,  method = 'monoH.FC')),
               temp_to_tau = list(splinefun(.x$quant_value,.x$quantile,  method = 'monoH.FC')))
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble() %>%
      mutate(Station = .x$Station[1])

  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()



# save quantile models
obs_smoothed_quantiles %>% saveRDS(paste0("output/quant_models_num_quantiles_",num_quantiles,".csv"))
#obs_smoothed_quantiles = readRDS("output/quant_models")



# Calculate lambda
lambda_thresh_ex = obs_data %>%
  group_by(Station) %>%
  group_map(~{

    print(.x$Station[1])

    print(obs_smoothed_quantiles%>%
            filter(Station == .x$Station[1]))

    thresh_exceedance_9 = obs_smoothed_quantiles%>%
      filter(Station == .x$Station[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$threshold_9[1], x))


    tibble(Station = .x$Station[1],
           year = temporal_covariates$year,
           thresh_exceedance_9 = 1-thresh_exceedance_9)
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()




lambda_thresh_ex %>%
  write_csv(paste0("data/processed/thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv"))


# ------------ get splines on clim scale

clim_grid = read_csv("data/processed/clim_scale_grid.csv")
clim_dat_full = read_csv("~/JRSS_organised_code/corrections/data/processed_data/full_clim_data.csv")

# estimate empiracle quantules for climate data
clim_quantiles = clim_dat_full %>%
  group_by(id) %>%
  group_map(~{
    res = c()
    for(q in quantiles_to_estimate_bulk){
      res = rbind(res, tibble(id = .x$id[1],
                              Long = .x$Long[1],
                              Lat = .x$Lat[1],
                              Long.projected = .x$Long.projected[1],
                              Lat.projected = .x$Lat.projected[1],
                              quantile = q,
                              value = as.numeric(quantile(.x$maxtp, q))))
    }
    res
  }, .keep = T)%>%
  plyr::rbind.fill() %>%
  as.tibble()

clim_quantiles_subset = clim_quantiles %>%
  # filter(id %in% obs_data$id) %>%
  group_by(id) %>%
  group_map(~{
    
    tibble(id = .x$id[1], quantile = list(.x$quantile), value = list(.x$value))
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()


clim_grid = clim_grid %>%
  left_join(clim_quantiles_subset)

clim_date_w_quantile_mod = c()


# here
for(i in seq(nrow(clim_grid))){
  print(clim_grid[i,])

  this_data_for_pred = tibble(year = c(1931, 1942, 2020, 2022)) %>%
    left_join(temporal_covariates) %>%
    mutate(Long = clim_grid[i,]$Long,
           Lat = clim_grid[i,]$Lat,
           Long.projected = clim_grid[i,]$Long.projected,
           Lat.projected = clim_grid[i,]$Lat.projected,
           id = clim_grid[i,]$id,
           # threshold = clim_grid[i,]$threshold,
           # clim_scale = clim_grid[i,]$clim_scale,
           quantile = clim_grid[i,]$quantile,
           value = clim_grid[i,]$value) 
  
  # predict quantile for each year and site
  quant_reg_pars = quant_reg_model_pars %>%
    arrange(tau)

  clim_vals = clim_grid[i,]$value[[1]]
  res = c()
  for(q in seq_along(quantiles_to_estimate_bulk)){
    qpars = quant_reg_pars[q,]
    res = rbind(res,
                tibble(quantile =  qpars$tau,
                       year = this_data_for_pred$year,
                       loess_temp_anom = this_data_for_pred$loess_temp_anom,
                       quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2)*(this_data_for_pred$loess_temp_anom)))
  }
  
  res = res %>%
    group_by(year) %>%
    group_map(~{
      tibble(year = .x$year[1],
             tau_to_temp = list(splinefun(.x$quantile,.x$quant_value,  method = 'monoH.FC')),
             temp_to_tau = list(splinefun(.x$quant_value,.x$quantile,  method = 'monoH.FC')))
    }, .keep = T) %>%
    plyr::rbind.fill() %>%
    as_tibble() 
  
  clim_date_w_quantile_mod = rbind(clim_date_w_quantile_mod, 
                                   this_data_for_pred %>% left_join(res))

}

clim_date_w_quantile_mod %>%
  saveRDS(paste0("output/quant_models_clim_num_quantiles_",num_quantiles,".csv"))

# thresold at these points 

qunatile_model_fit_9 = readRDS("output/threshold_model_9")
clim_date_w_quantile_mod = readRDS(paste0("output/quant_models_clim_num_quantiles_",num_quantiles,".csv"))

# get clim threshold values
clim_thresh = clim_dat_full %>%
  group_by(id) %>%
  summarise(clim_thresh_value_9 = quantile(maxtp, 0.9))

clim_date_w_quantile_mod = clim_date_w_quantile_mod %>% left_join(clim_thresh)
clim_date_w_quantile_mod$threshold_9= predict(qunatile_model_fit_9, clim_date_w_quantile_mod)$location

#---- get lambda for climate model
# Calculate lambda
lambda_thresh_ex = clim_date_w_quantile_mod %>%
  group_by(id) %>%
  group_map(~{
    
    thresh_exceedance_9 = clim_date_w_quantile_mod%>%
      filter(id == .x$id[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$threshold_9[1], x))
    
    
    tibble(id = .x$id[1],
           year = c(1931, 1942, 2020, 2022),
           thresh_exceedance_9 = 1-thresh_exceedance_9)  
    
    
    
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()


lambda_thresh_ex %>% 
  write_csv(paste0("data/processed/climate_thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv"))

clim_date_w_quantile_mod %>% 
  left_join(lambda_thresh_ex) %>%
  saveRDS(paste0("output/quant_models_clim_num_quantiles_",num_quantiles,".csv"))


# read_csv("data/processed/thresh_exceedance_lambda_num_quantiles_100.csv")
# 


# 
# # ----- plot thresh exceedance with uncertainty and effect of temporal covariate
# 
# 
# 
# beta_2_plot=bts_quant_reg %>%
#   filter(bts != 96) %>%
#   group_by(tau) %>%
#   summarise(beta_a_2_u = quantile(beta_a_2, 0.975),
#             beta_a_2_l = quantile(beta_a_2, 0.025)) %>%
#   ggplot()+
#   #geom_segment(aes(x = tau, xend = tau,  y = beta_a_2_l, yend = beta_a_2_u), alpha = 0.25)+
#   geom_ribbon(aes(x = tau, ymin = beta_a_2_l, ymax = beta_a_2_u), alpha = 0.25)+
#   geom_line(data = quant_reg_pars %>% filter(tau<0.9, tau> 0.05), aes(tau,beta_a_2))+
#   ylim(-2, 2.7)+
#   theme_minimal(12)+
#   labs(y = expression(beta[2]),
#        x = expression(tau))+
#   theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+
#   scale_x_continuous(breaks=c(0.05, 0.25, 0.5, 0.75, 0.95))
# 
# 
# lambda_plot = true_lambda %>%
#   ggplot()+
#   geom_ribbon(data = bts_lam, aes(x = year, ymin = thresh_exceedance_c_l, ymax = thresh_exceedance_c_u), alpha = 0.25)+
#   geom_line(data = true_lambda %>% group_by(year) %>% summarise(thresh_exceedance = mean(thresh_exceedance)), aes(year, thresh_exceedance))+
#   theme_minimal(12)+
#   labs(colour = expression(s),
#        linetype =  expression(s),
#        y = expression(lambda),
#        x = "Year")+
#   theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
# 
# 
# 
# plot_fig_3 = gridExtra::grid.arrange(beta_2_plot, lambda_plot, nrow =1)
# # # Author: Dáire Healy
# # # Date: 17 Feb 2022
# # 
# # # Description: This script prepares observational temperature data for marginal 
# # # modelling. This involves incorporating all covariates as descriped in the 
# # # paper associated to this work.
# # # Note: bulk modelling in preformed in this script. 
# # gc()
# # rm(list = ls())
# # library(tidyverse)
# # library(vroom)
# # library(evgam)
# # library(mgcv)
# # library(raster) # package for netcdf manipulation
# # 
# # setwd("~/JRSS_organised_code/")
# # 
# # # i. ========  ========  Global parameters ========  ========  ========
# # 
# # # marginal threshold 
# # marg_thresh_quantile = 0.8 # equiv to 95th quantile of highest 1/4 of data
# # 
# # 
# # # 1. ========  ========  Read in all data  ========  ========  ======== 
# # # 1a. ---- for plotting 
# # # colour palatte 
# # my_pal = c( 
# #   '#062c30', # extra dark 
# #   '#003f5c',
# #   '#2f4b7c',
# #   '#665191',
# #   '#a05195',
# #   '#d45087',
# #   '#f95d6a',
# #   '#ff7c43',
# #   '#ffa600')
# # 
# # # map of ireland sf
# # ireland_sf <-rnaturalearth::ne_countries(type = "map_units",
# #                                          scale='large', # resolution of map
# #                                          returnclass = 'sf', # spatial object
# #                                          continent = "europe") %>%
# #   .[.$name %in% c('Ireland', 'N. Ireland'),]
# # 
# # # 1b. ---- observational data
# # obs_data = read_csv("~/JRSS_Paper/processed_data/obs_all_data.csv")
# # 
# # # removing data that is not fit for modeking (has been descretised)
# # obs_data = obs_data %>%
# #   filter(!(Station == "aldergrove" &  lubridate::year(date) <=1960),
# #          !(Station == "altnahinch-filters" &  lubridate::year(date) <=1970),
# #          !(Station == "armagh" &  lubridate::year(date) <=1970),
# #          !(Station == "ballypatrick-forest" &  lubridate::year(date) <=1968),
# #          !(Station == "hillsborough" &  lubridate::year(date) <=1970),
# #          !(Station == "lough-navar-forest" &  lubridate::year(date) <=1968),
# #          !(Station == "murlough" &  lubridate::year(date) <=1970),
# #          !(Station == "stormont-castle" &  lubridate::year(date) <=1970),
# #          !(Station == "waterford (tycor)" &  lubridate::year(date) <=1962),
# #          !(Station == "littleton ii b. na m." &  lubridate::year(date) <=1962),
# #          !(Station == "helens-bay" &  lubridate::year(date) <=1970)) %>%
# #   filter(!Station %in% c("cavan (drumconnick)", 
# #                          "recess (cloonacartan)",
# #                          "tuam (airglooney)",
# #                          "dooks",
# #                          "athy (chanterlands)")) %>%
# #   filter(!(Station == "helens-bay" &  maxtp >= 32),
# #          !(Station == "portglenone" &  maxtp >= 30))  %>%
# #   filter(!(Station == "dublin (phoenix park)" &  as.Date(date) > 
# #              (obs_data %>%
# #                 filter(Station == "phoenixpark") %>%
# #                 pull(date) %>%
# #                 as.Date() %>%
# #                 min)))%>%
# #   mutate(Station = as.character(Station),
# #          date = as.Date(date),
# #          year = lubridate::year(date)) %>%
# #   filter(lubridate::month(date) %in% months_to_study) # just take summer months
# # 
# # # project coordinates
# # my_coords = obs_data %>% dplyr::select(Long, Lat)
# # sp::coordinates(my_coords)<- c("Long", "Lat")
# # 
# # # add Coordinate Reference System, which is WGS84
# # sp::proj4string(my_coords) <- sp::CRS("+proj=longlat +datum=WGS84") 
# # 
# # # project to utm zone 29
# # proj_cords <- sp::spTransform(my_coords, 
# #                               sp::CRS(paste0("+proj=utm +zone=29 ellps=WGS84"))) 
# # proj_cords = proj_cords %>% sp::coordinates()/100000 # scale coordinates
# # obs_data$Long.projected = proj_cords[,1]
# # obs_data$Lat.projected = proj_cords[,2]
# # 
# # # 1b. ---- climate data
# # historic = vroom("data/raw_data/clim/r12i1p1-ICHEC-EC-EARTH-(Ireland)CLMcom-CLM-CCLM4-8-17 (EU).csv")
# # clim_data = historic %>%
# #   filter(lubridate::month(date) %in% months_to_study) # only take summer months 
# # rm(list = c('historic')) # clear memory
# # 
# # # project coordinates
# # my_coords = clim_data %>% dplyr::select(Long, Lat)
# # sp::coordinates(my_coords)<- c("Long", "Lat")
# # sp::proj4string(my_coords) <- sp::CRS("+proj=longlat +datum=WGS84") 
# # proj_cords <- sp::spTransform(my_coords, 
# #                               sp::CRS(paste0("+proj=utm +zone=29 ellps=WGS84"))) 
# # proj_cords = proj_cords %>% sp::coordinates()/100000 
# # clim_data$Long.projected = proj_cords[,1]
# # clim_data$Lat.projected = proj_cords[,2]
# # clim_data$id = clim_data %>% group_indices(Long, Lat)
# # 
# # 
# # # 1c. Add temporal covariates (HADCRUT5)
# # hadcrut_data = raster::brick("corrections/HadCRUT.5.0.1.0.anomalies.ensemble_mean.nc")
# # 
# # hadcrut_loc_dat = raster::coordinates(hadcrut_data) %>% as_tibble()
# # names(hadcrut_loc_dat) = c('Long', 'Lat')
# # 
# # hadcrut_data = hadcrut_data %>% values() %>% as_tibble()
# # hadcrut_data$Long = hadcrut_loc_dat$Long
# # hadcrut_data$Lat = hadcrut_loc_dat$Lat
# # names(hadcrut_data) = names(hadcrut_data) %>% str_remove_all("X")
# # 
# # hadcrut_data = hadcrut_data %>%
# #   tidyr::pivot_longer(c(-Long, -Lat), names_to = "date", values_to = "temp_anom") %>% # Making each row an observation
# #   mutate(date = date %>% lubridate::ymd())
# # 
# # 
# # irealnd_hadcrut = hadcrut_data %>%
# #   filter(Long > -11, Long < -6,   Lat > 50, Lat < 56) %>%
# #   mutate(year = lubridate::year(date)) %>%
# #   mutate(month = lubridate::month(date)) %>%
# #   filter(month %in% c(6,7,8)) 
# # 
# # global_hadcrut = hadcrut_data %>% 
# #   drop_na() %>% 
# #   mutate(year = lubridate::year(date))%>%
# #   mutate(month = lubridate::month(date)) %>%
# #   filter(month %in% c(6,7,8)) %>% 
# #   group_by(year, month) %>%
# #   summarise(temp_anom = mean(temp_anom)) %>% ungroup() %>% unique()
# # 
# # 
# # loess_fit_irel = predict(loess(temp_anom ~ year, data = irealnd_hadcrut), se=T)
# # loess_fit_glob = loess(temp_anom ~ year, data = global_hadcrut)
# # 
# # irealnd_hadcrut$loess_temp_anom = loess_fit_irel$fit
# # irealnd_hadcrut$loess_temp_anom_u = loess_fit_irel$fit+1.96*(loess_fit_irel$se.fit)
# # irealnd_hadcrut$loess_temp_anom_l = loess_fit_irel$fit-1.96*(loess_fit_irel$se.fit)
# # 
# # loess_glob_preds = predict(loess_fit_glob, irealnd_hadcrut, se = T)
# # irealnd_hadcrut$loess_glob_temp_anom = loess_glob_preds$fit
# # irealnd_hadcrut$loess_glob_temp_anom_u = loess_glob_preds$fit+1.96*(loess_glob_preds$se.fit)
# # irealnd_hadcrut$loess_glob_temp_anom_l = loess_glob_preds$fit-1.96*(loess_glob_preds$se.fit)
# # 
# # 
# # irealnd_hadcrut %>%
# #   ggplot()+
# #   geom_point(aes(year, temp_anom))+
# #   geom_ribbon(aes(x = year,ymin = loess_temp_anom_l, ymax = loess_temp_anom_u, fill = 'Ireland'), alpha = 0.25)+
# #   geom_ribbon(aes(x = year,ymin = loess_glob_temp_anom_l, ymax = loess_glob_temp_anom_u, fill = 'Global'), alpha = 0.25)+
# #   geom_line(aes(year, loess_temp_anom, col = 'Ireland'), size = 1)+
# #   geom_line(aes(year, loess_glob_temp_anom, col = 'Global'), linetype = 'dashed', size = 1)+
# #   theme_minimal()+
# #   theme_minimal(12)+
# #   theme(legend.position = 'none')+
# #   labs(x = 'year',
# #        y = "Temperature anomaly")
# # 
# # 
# # obs_data = obs_data %>%
# #   left_join((irealnd_hadcrut %>%
# #                dplyr::select(year, loess_temp_anom, loess_glob_temp_anom) %>%
# #                unique()), by = 'year')
# # 
# # 
# # 
# # # 1d. get climate subset that matches up with observational stations 
# # #  ------- calculate the distance between two points
# # dist_h <- function(long1, lat1, long2, lat2) {
# #   sqrt((long1 - long2)^2 + (lat1 - lat2)^2)
# # }
# # 
# # obs_sites = obs_data %>%
# #   dplyr::select(Station, Long, Lat, Long.projected, Lat.projected) %>%
# #   unique()
# # 
# # clim_sites = clim_data %>%
# #   dplyr::select(id, Long, Lat, Long.projected, Lat.projected) %>%
# #   unique()
# # 
# # # itterate through all stations and find climate grid point
# # closest_irel = c()
# # for(i in seq(nrow(obs_sites))){
# #   smallest_dist = 9999999999
# #   id_of_smallest_dist_irel = 9999999999
# #   
# #   for(j in seq(nrow(clim_sites))){
# #     this_dist = dist_h(obs_sites[i,]$Long.projected,
# #                        obs_sites[i,]$Lat.projected,
# #                        clim_sites[j,]$Long.projected,
# #                        clim_sites[j,]$Lat.projected)
# #     
# #     if(this_dist < smallest_dist){
# #       smallest_dist = this_dist
# #       id_of_smallest_dist_irel = clim_sites$id[j]
# #     }
# #   }
# #   closest_irel = c(closest_irel, id_of_smallest_dist_irel)
# # }
# # 
# # # add column "id" to obs data which is the id of the closest clim grid point
# # obs_sites$id = closest_irel
# # obs_data = obs_data %>%
# #   left_join(obs_sites)
# # 
# # 
# # # procedure
# # # 1. define \tau = 0.1, 0.15, ... , 0.99
# # # 2. caluate \tau-th qunautle at each site for climate model data q_\tau(X_c(s))
# # # 3. Preform spatio temporal qunatile regression on staiton data for each \tau
# # # ----- X_o(s,t; \tau) ~ ALD(q_\tau(X_c(s)), hadcrut(t); \tau) 
# # 
# # 
# # 
# # 
# # 
# # 
# # # ---- only use climate model data as spatial covariate
# # # ---- take empirical quantiles of climate model
# # quantiles_to_estimate_bulk = seq(0.1,0.999,length.out = 30)
# # 
# # clim_quantiles = clim_data %>%
# #   group_by(id) %>%
# #   group_map(~{
# #     res = c()
# #     for(q in quantiles_to_estimate_bulk){
# #       res = rbind(res, tibble(id = .x$id[1],
# #                               Long = .x$Long[1],
# #                               Lat = .x$Lat[1],
# #                               Long.projected = .x$Long.projected[1],
# #                               Lat.projected = .x$Lat.projected[1],
# #                               quantile = q,
# #                               value = as.numeric(quantile(.x$maxtp, q))))
# #     }
# #     res
# #   }, .keep = T)%>%
# #   plyr::rbind.fill() %>%
# #   as.tibble()
# # 
# # 
# # 
# # 
# # clim_quantiles_subset = clim_quantiles %>%
# #   filter(id %in% obs_data$id) %>%
# #   group_by(id) %>%
# #   group_map(~{
# #     
# #     tibble(id = .x$id[1], quantile = list(.x$quantile), value = list(.x$value))
# #   }, .keep = T) %>%
# #   plyr::rbind.fill() %>%
# #   as_tibble()
# # 
# # 
# # obs_data = obs_data %>%
# #   dplyr::select(-value, -quantile) %>%
# #   left_join(clim_quantiles_subset)
# # 
# # 
# # # for each tau, preform constatn, spatial, temporal + spatio-temporal regression 
# # obs_q_mod = c()
# # for(q in seq_along(quantiles_to_estimate_bulk)){
# #   print(q)
# #   obs_data_for_quant_reg = obs_data
# #   
# #   zeta = obs_data$quantile[[1]][q]
# #   obs_data_for_quant_reg$value = obs_data$value %>% lapply(`[[`, q) %>% unlist
# #   
# #   qunatile_model <-  maxtp ~ value + s(loess_temp_anom)
# #   qunatile_model_fit <- evgam(qunatile_model,
# #                               obs_data_for_quant_reg,
# #                               family = "ald",
# #                               ald.args = list(tau = zeta))
# #   
# #   obs_q_mod = rbind(obs_q_mod, tibble(q, quant_mod = list(qunatile_model_fit)))
# # }
# # 
# # 
# # 
# # 
# # loess_temp_anom_vals = obs_data %>%
# #   dplyr::select(year, loess_temp_anom) %>%
# #   unique() %>%
# #   pull(loess_temp_anom) %>%
# #   as.numeric()
# # 
# # 
# # # 
# # # 
# # # 
# # # get_obs_quantiles = function(my_stations, my_years){
# # #   
# # #   irealnd_hadcrut %>%
# # #     dplyr::select(year, loess_temp_anom) %>%
# # #     filter(year %in% my_years) %>%
# # #     unique() %>%
# # #     pull(loess_temp_anom) %>%
# # #     as.numeric()
# # #   
# # #   
# # #   obs_data %>%
# # #     filter(Station %in% my_stations) %>%
# # #     dplyr::select(quantile, value)
# # #   
# # # }
# # 
# # 
# # # temporal covariate 
# # loess_temp_anom_vals = irealnd_hadcrut %>%
# #   dplyr::select(year, loess_temp_anom) %>%
# #   filter(year %in% seq(1942, 2020)) %>%
# #   unique() %>%
# #   pull(loess_temp_anom) %>%
# #   as.numeric()
# # 
# # 
# # 
# # obs_smoothed_quantiles = obs_data %>%
# #   group_by(Station) %>%
# #   group_map(~{
# #     
# #     print(.x$Station[1])
# #     # get the climate qunatile estimates closest to current station
# #     clim_vals = obs_data %>%
# #       filter(Station == .x$Station[1]) %>%
# #       dplyr::select(quantile, value) %>%
# #       unique() %>%
# #       pull(value) %>%
# #       unlist()
# #     
# #     
# #     # for each tau, for each year predict obs quantile
# #     res = c()
# #     for(q in seq_along(quantiles_to_estimate_bulk)){
# #       res = rbind(res,tibble(year = seq(1942, 2020), 
# #                              quantile =  quantiles_to_estimate_bulk[q], 
# #                              value = obs_q_mod$quant_mod %>% lapply(predict, tibble(value = clim_vals[q], loess_temp_anom = loess_temp_anom_vals)) %>% .[[1]] %>% .$location))
# #     }
# #     
# #     
# #     # for each year interpolate over quantiles
# #     res %>% group_by(year) %>%
# #       group_map(~{
# #         
# #         tibble(year = .x$year[1], quant_spline = list(splinefun(.x$value,.x$quantile,  method = 'monoH.FC')))
# #         
# #       }, .keep = T) %>%
# #       plyr::rbind.fill() %>%
# #       as_tibble() %>%
# #       mutate(Station = .x$Station[1])
# #     
# #   }, .keep = T) %>%
# #   plyr::rbind.fill() %>%
# #   as_tibble()
# # 
# # 
# # # --- get thresold at each site
# # clim_thresh_values = clim_data %>%
# #   filter(id %in% obs_data$id) %>%
# #   group_by(id) %>%
# #   summarise(clim_thresh_value = quantile(maxtp, 0.8))
# # 
# # obs_data = obs_data %>%
# #   left_join(clim_thresh_values)
# # 
# # qunatile_model <-  maxtp ~ clim_thresh_value
# # qunatile_model_fit <- evgam(qunatile_model,
# #                             obs_data,
# #                             family = "ald",
# #                             ald.args = list(tau = 0.8))
# # obs_data$threshold = qunatile_model_fit$location$fitted
# # 
# # 
# # 
# # thresh_exceedance_prob = obs_data %>%
# #   group_by(Station) %>%
# #   group_map(~{
# #     
# #     thresh_exceedance = obs_smoothed_quantiles%>%
# #       filter(Station == .x$Station[1]) %>%
# #       pull(quant_spline) %>%
# #       sapply(function(x) sapply(.x$threshold[1], x))
# #     
# #     tibble(Station = .x$Station[1],year = seq(1942, 2020), thresh_exceedance = 1-thresh_exceedance) 
# #     
# #   }, .keep = T) %>%
# #   plyr::rbind.fill() %>%
# #   as_tibble()
# # 
# # 
# # thresh_exceedance_prob %>%
# #   ggplot()+
# #   geom_line(aes(year, thresh_exceedance, group = Station), alpha = 0.25) +
# #   theme_minimal()+
# #   labs(x = "Year",
# #        y = expression(lambda[t]))
# # 
# # 
# # 
# # # 
# # # thresh_exceedance_prob %>%
# # #   ggplot()+
# # #   geom_line(aes(quantile, value, group = year, col = year))+
# # #   geom_hline(yintercept = 20.1)+
# # #   theme_minimal()+
# # #   xlim(0.75, 0.9)+
# # #   ylim(18, 22)+
# # #   labs(x = "Quantile",
# # #        y = "Temperature",
# # #        col = "Year")
# # # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # # clim_data$maxtp %>% var
# # # 
# # # clim_data = clim_data %>% mutate(year = lubridate::year(date))
# # # lm(maxtp~year, data = clim_data) %>% residuals() %>% var
# # # 
# # # gam(maxtp~s(year), data = clim_data)%>% residuals() %>% var
# # # 
# # # 1 - (8.129504/8.45141)
# # 
# # # # look at one climate grid point
# # # slope_varying_quant_low = clim_data %>%
# # #   group_by(id) %>%
# # #   group_map(~{
# # #     print(.x$id[1])
# # #     for_qnt_reg = .x %>%
# # #       mutate(year = lubridate::year(date)) %>%
# # #       left_join(irealnd_hadcrut %>% dplyr::select(year, loess_temp_anom) %>% unique)
# # #     
# # #     qunatile_model_fit <- evgam(maxtp ~ loess_temp_anom,
# # #                                 for_qnt_reg,
# # #                                 family = "ald",
# # #                                 ald.args = list(tau = 0.1))
# # #     summary(qunatile_model_fit) 
# # #     
# # #     
# # #     tibble(.x$id[1], Long = .x$Long[1],   Lat = .x$Lat[1], slope = qunatile_model_fit$location$coefficients[2])
# # #   }, .keep = T) %>%
# # #   plyr::rbind.fill() %>%
# # #   as.tibble
# # # 
# # # gridExtra::grid.arrange(slope_varying_quant_low %>%
# # #                           ggplot()+
# # #                           geom_point(aes(Long, Lat, col = slope))+
# # #                           labs(title = "0.1 qnt")+
# # #                           theme_minimal()+
# # #                           coord_map(),
# # #                         slope_varying_quant %>%
# # #                           ggplot()+
# # #                           geom_point(aes(Long, Lat, col = slope))+
# # #                           labs(title = "0.8 qnt")+
# # #                           theme_minimal()+
# # #                           coord_map(), nrow = 1)
# # # temporally changing quantiles
# # 
# # obs_data_for_pred = obs_data %>%
# #   dplyr::select(Station, id) %>% unique()
# # 
# # obs_data_for_pred = tibble(Station = rep(obs_data_for_pred$Station,length(seq(1940, 2020))) %>% sort,
# #                   year = rep(seq(1940, 2020), nrow(obs_data_for_pred))) %>%
# #   left_join(obs_data_for_pred) %>%
# #   left_join(irealnd_hadcrut %>% dplyr::select(year, loess_temp_anom) %>% unique) 
# # 
# # obs_sites_quantile_estimtes = clim_quantiles %>%
# #   group_by(quantile) %>%
# #   group_map(~{
# #     obs_data_for_quant_reg = .x %>%
# #       right_join(obs_sites %>% 
# #                    dplyr::select(Station, id)) %>%
# #       dplyr::select(Station, quantile, value) %>%
# #       right_join(obs_data) %>%
# #       dplyr::select(Station, maxtp,
# #                     quantile, value,
# #                     Long, Lat,
# #                     Long.projected, Lat.projected, 
# #                     year,
# #                     loess_temp_anom)
# #     
# #     paste0("quantile_model_",signif(.x$quantile[1], 2)) %>% print
# #     zeta = obs_data_for_quant_reg$quantile[1]
# # 
# #     
# #     qunatile_model <-  maxtp ~ value + loess_temp_anom
# # 
# #     
# #     qunatile_model_fit <- evgam(qunatile_model,
# #                                 obs_data_for_quant_reg,
# #                                 family = "ald",
# #                                 ald.args = list(tau = zeta))
# # 
# #     # save model
# #     saveRDS(qunatile_model_fit,
# #             file = paste0("corrections/bulk_model/models/qunatile_mod_",signif(.x$quantile[1], 2)))
# # 
# # 
# #     
# #     this_obs_data_for_pred = obs_data_for_pred %>% left_join(.x %>%
# #                                         dplyr::select(id, value) %>%
# #                                         unique())
# #     
# #     this_obs_data_for_pred$quantile = zeta
# #     this_obs_data_for_pred$quantile_value = predict(qunatile_model_fit, this_obs_data_for_pred)$location    
# #     
# #     rm(ls = "qunatile_model_fit")
# #     this_obs_data_for_pred
# #     #  
# #   }, .keep = T)%>%
# #   plyr::rbind.fill() %>%
# #   as.tibble()
# # 
# # 
# # # interpolate and save each model
# # obs_sites_quantile_estimtes %>%
# #   group_by(Station, year) %>%
# #   group_map(~{
# #     
# #     this_fitted_mod = stats::splinefun(x = .x$quantile, 
# #                                   y = .x$quantile_value, 
# #                                   method = "monoH.FC",
# #                                   ties = mean)
# #     
# #     this_fitted_mod %>% saveRDS(paste0("corrections/quantile_regression_models/obs_site_year_mods/",.x$Station[1], "_", .x$year[1]))
# #   }, .keep = T) 
# # 
# # 
# # 
# # obs_sites_quantile_estimtes %>%
# #   filter(Station %in% c('claremorris', "mullingar", "malin_head", "armagh")) %>%
# #   ggplot()+
# #   geom_line(aes(quantile, quantile_value, col = year, group = year))+
# #   theme_minimal()+
# #   facet_wrap(~Station, nrow = 1)+
# #   ggplot2::scale_colour_gradientn(colours = my_pal)
# # 
# # 
# # for_interpolation = obs_sites_quantile_estimtes %>%
# #   filter(Station == 'claremorris') %>%
# #   filter(year == 2020) %>%
# #   dplyr::select(quantile, quantile_value)
# #   
# # 
# # 
# # 
# # 
# # 
# # 
# # ?crps
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # obs_bulk_models = obs_sites_quantile_estimtes %>%
# #   group_by(Station) %>%
# #   group_map(~{
# #     
# #     print(.x$Station[1])
# #     this_site_thresh = obs_data %>%
# #       filter(Station == .x$Station[1]) %>%
# #       pull(threshold) %>%
# #       .[1] %>%
# #       print()
# #     
# #     this_site_data = obs_data %>%
# #       filter(Station == .x$Station[1]) %>%
# #       pull(maxtp)
# #     
# #     #this_thresh_quantile = ecdf(this_site_data)(this_site_thresh) %>% print
# #     this_thresh_quantile = 0.8
# #     
# #     # -- make sure nothing greater than threshold
# #     .x = .x %>%
# #       filter(qunantile_value < this_site_thresh,
# #              quantile < this_thresh_quantile)
# #     
# #     data_for_quant_spline = tibble(Station = .x$Station[1],
# #                                    quantile = c(.x$quantile, this_thresh_quantile),
# #                                    qunantile_value = c(.x$qunantile_value, this_site_thresh))
# #     
# #     k <- 10 #num knots
# #     knots <- data.frame(x = seq(0, this_thresh_quantile,length=k)) 
# #     sm <- smoothCon(s(quantile,k=k,bs="cr"), data_for_quant_spline, knots=knots)[[1]]
# #     
# #     # value to force spline through
# #     fixed_y_val = this_site_thresh
# #     
# #     X <- sm$X[, -k]
# #     S <- sm$S[[1]][-k, -k]  
# #     off <-  data_for_quant_spline$qunantile_value*0 + (fixed_y_val)
# #     
# #     # fit spline
# #     gam_mod <- gam(data_for_quant_spline$qunantile_value ~ X -1 + offset(off), paraPen = list(X = list(S)))
# #     
# #     tibble(Station = .x$Station[1], sm = list(sm), k, gam_mod = list(gam_mod), fixed_y_val, this_thresh_quantile)
# #     
# #   }, .keep = T) %>%
# #   plyr::rbind.fill() %>%
# #   as_tibble()
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # quantile_reg_mod = readRDS("corrections/bulk_model/models/qunatile_mod_0.79")
# # 
# # data_for_testing = clim_quantiles %>%
# #   filter(quantile >= 0.78, quantile <= 0.795) 
# # 
# # 
# # data_for_testing = rbind(data_for_testing%>%
# #          mutate(year = 1942,loess_temp_anom = irealnd_hadcrut %>% filter(year == 1942) %>% pull(loess_temp_anom) %>% unique),
# #        data_for_testing%>%
# #            mutate(year = 2020,loess_temp_anom = irealnd_hadcrut %>% filter(year == 2020) %>% pull(loess_temp_anom) %>% unique))
# # 
# # 
# #     
# #     
# # data_for_testing$predicted_quantile = predict(quantile_reg_mod, data_for_testing)$location
# # 
# # data_for_testing %>% ggplot()+
# #   geom_point(aes(Long, Lat, col = predicted_quantile)) +
# #   coord_map()+
# #   ggplot2::scale_colour_gradientn(colours = my_pal)+
# #   facet_wrap(~year)+
# #   theme_minimal()+
# #   labs(title = "0.8th quantile")
# # 
# # 
# # 
# # gridExtra::grid.arrange(plt_1942,plt_2020, nrow=1)
# # 
# # 
# # 
# # 
# # # print("fitting threshold")
# # thresh_mod <-  maxtp ~ s(Long.projected, Lat.projected) 
# # thresh_mod_fit <- evgam(thresh_mod, obs_data, family = "ald", ald.args = list(tau = marg_thresh_quantile))
# # thresh_mod_fit %>% saveRDS("corrections/thresh_mod_fit")
# # #thresh_mod_fit = readRDS("corrections/thresh_mod_fit")
# # obs_data$threshold = thresh_mod_fit$location$fitted 
# # 
# # obs_data %>%
# #   dplyr::select(Long, Lat, threshold) %>%
# #   unique() %>% ggplot()+
# #   geom_point(aes(Long, Lat, col = threshold))
# # 
# # 
# # # 2. ========  ========  Bulk Model        ========  ========  ======== 
# # # 2a. ---- Climate model
# # # calculate quantiles for climate data
# # print("declustering climate data")
# # 
# # # 
# # # extreme_clim_data = clim_data %>%
# # #   filter(id %in% obs_sites$id) %>%
# # #   left_join(all_dates) %>%
# # #   filter(maxtp>threshold)
# # # 
# # # all_dates = tibble(date = seq(extreme_clim_data$date %>% min,
# # #                               extreme_clim_data$date %>% max, by = 'day'))
# # # 
# # # date_groups = map(seq(ceiling(nrow(all_dates)/5)), function(x) rep(x,5)) %>% unlist
# # # 
# # # date_groups = date_groups[1:nrow(all_dates)]
# # # all_dates$date_groups = date_groups
# # # 
# # # extreme_clim_data = extreme_clim_data %>%
# # #   left_join(all_dates) %>%
# # #   group_by(date_groups, id) %>%
# # #   filter(maxtp == max(maxtp)) %>%
# # #   ungroup() %>%
# # #   select(-date_groups)
# # # 
# # 
# # # 
# # # clim_data_subset_declustered = rbind(clim_data %>%
# # #                                        filter(id %in% obs_sites$id) %>% 
# # #                                        filter(maxtp<=threshold), extreme_data)
# # # 
# # # 
# # 
# # 
# # 
# # # ---- last run 5/mar 13:31
# # # clim_data_subset_declustered = clim_data %>%
# # #   filter(id %in% obs_sites$id) %>%
# # #   left_join(all_dates) %>%
# # #   group_by(date_groups, id) %>%
# # #   filter(maxtp == max(maxtp)) %>%
# # #   ungroup()
# # 
# # # extreme_clim_data %>%
# # #    write_csv("corrections/clim_data_subset_declustered.csv")
# # #clim_data_subset_declustered = read_csv("corrections/clim_data_subset_declustered.csv")
# # 
# # 
# # 
# # # figure out where in the data (what quantile) the threshold is in the
# # # declustered data and do quantule regression upto this 
# # 
# # # quantile for each site
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # # Quantiles to at which to preform quantile regression for bulk model
# # quantiles_to_estimate_bulk = seq(0.01, 0.8, length.out = 30)
# # 
# # print("calculating quantiles 0-marg_thresh_quantile for climate data ... ")
# # 
# # 
# # # taking empirical quantiles from climate model data
# # clim_quantiles = clim_data %>%
# #   group_by(id) %>%
# #   group_map(~{
# #     res = c()
# #     for(q in quantiles_to_estimate_bulk){
# #       res = rbind(res, tibble(id = .x$id[1],
# #                               Long = .x$Long[1],
# #                               Lat = .x$Lat[1],
# #                               Long.projected = .x$Long.projected[1],
# #                               Lat.projected = .x$Lat.projected[1],
# #                               quantile = q,
# #                               value = as.numeric(quantile(.x$maxtp, q))))
# #     }
# #     res
# #   }, .keep = T)%>%
# #   plyr::rbind.fill() %>%
# #   as.tibble()
# # 
# # 
# # 
# # # 2b. ---- Observational data
# # # # --- quantile regression with climate quantile as covariate (~2 hours)
# # # # --- only need to run once, results saved to "quantile_regression_models"
# # # # # last ran 7/3/22 0900
# # # 
# # # 
# # # # ----------
# # clim_quantiles %>%
# #   group_by(quantile) %>%
# #   group_map(~{
# #     
# #     obs_data_for_quant_reg = .x %>%
# #       right_join(obs_sites %>% 
# #                    dplyr::select(Station, id)) %>%
# #       dplyr::select(Station, quantile, value) %>%
# #       right_join(obs_data) %>%
# #       dplyr::select(Station, maxtp,
# #                     quantile, value,
# #                     Long, Lat,
# #                     Long.projected, Lat.projected)
# #     
# #     paste0("quantile_model_",signif(.x$quantile[1], 2)) %>% print
# #     zeta = obs_data_for_quant_reg$quantile[1]
# #     
# #     print(obs_data_for_quant_reg)
# #     #qunatile_model <-  maxtp ~ s(Long.projected, Lat.projected) + value
# #     qunatile_model <-  maxtp ~ value
# #     
# #     qunatile_model_fit <- evgam(qunatile_model,
# #                                 obs_data_for_quant_reg,
# #                                 family = "ald",
# #                                 ald.args = list(tau = zeta))
# #     
# #     saveRDS(qunatile_model_fit,
# #             file = paste0("corrections/quantile_regression_models/quantile_model_",
# #                           signif(.x$quantile[1], 2)))
# #     rm(ls = "qunatile_model_fit")
# #   }, .keep = T)
# # 
# # 
# # 
# # 
# # 
# # # Interpolate quantiles
# # #----- predict obs threshold with fitted ALD models
# # # read in each quantile regression model, and get its predicted value for each 
# # # obsevational station
# # # using the fitted ALD quantile regression models, get estimated quantiles for 
# # # each site and each quantile.
# # 
# # new_mod = readRDS("corrections/quantile_regression_models/quantile_model_")
# # summary(new_mod)
# # 
# # print("interpolating qualtiles")
# # obs_sites_quantile_estimtes = clim_quantiles %>%
# #   group_by(quantile) %>%
# #   group_map(~{
# #     obs_data_for_quant_reg = .x %>%
# #       right_join(obs_sites %>% dplyr::select(Station, id)) %>%
# #       dplyr::select(Station, quantile, value) %>%
# #       right_join(obs_data) %>%
# #       dplyr::select(Station, maxtp,
# #                     quantile, value,
# #                     Long, Lat,
# #                     Long.projected, Lat.projected)
# #     
# #     qunatile_model_fit = readRDS(paste0("corrections/quantile_regression_models/quantile_model_",
# #                                         signif(.x$quantile[1], 2)))
# #     
# #     obs_data_for_quant_reg$qunantile_value = qunatile_model_fit$location$fitted
# #     
# #     obs_data_for_quant_reg %>%
# #       dplyr::select(Station, quantile, qunantile_value,
# #                     Long, Lat,
# #                     Long.projected, Lat.projected) %>%
# #       unique()
# #   }, .keep = T) %>%
# #   plyr::rbind.fill() %>%
# #   as_tibble()
# # 
# # 
# # # take threshold informed by climate data
# # obs_data = obs_data %>%
# #   left_join(obs_sites_quantile_estimtes %>% filter(quantile == 0.8)  %>%
# #               dplyr::select(Station, threshold_w_clim = qunantile_value))
# # 
# # 
# # 
# # # do monotonic interpoalation 
# # # splinefun {stats} - monoH.FC
# # 
# # # interpolate quantuiles (without clim threshold)
# # obs_bulk_models = obs_sites_quantile_estimtes %>%
# #   group_by(Station) %>%
# #   group_map(~{
# #     
# #     print(.x$Station[1])
# #     this_site_thresh = obs_data %>%
# #       filter(Station == .x$Station[1]) %>%
# #       pull(threshold) %>%
# #       .[1] %>%
# #       print()
# #     
# #     this_site_data = obs_data %>%
# #       filter(Station == .x$Station[1]) %>%
# #       pull(maxtp)
# #     
# #     #this_thresh_quantile = ecdf(this_site_data)(this_site_thresh) %>% print
# #     this_thresh_quantile = 0.8
# #     
# #     # -- make sure nothing greater than threshold
# #     .x = .x %>%
# #       filter(qunantile_value < this_site_thresh,
# #              quantile < this_thresh_quantile)
# #     
# #     data_for_quant_spline = tibble(Station = .x$Station[1],
# #                                    quantile = c(.x$quantile, this_thresh_quantile),
# #                                    qunantile_value = c(.x$qunantile_value, this_site_thresh))
# #     
# #     k <- 10 #num knots
# #     knots <- data.frame(x = seq(0, this_thresh_quantile,length=k)) 
# #     sm <- smoothCon(s(quantile,k=k,bs="cr"), data_for_quant_spline, knots=knots)[[1]]
# #     
# #     # value to force spline through
# #     fixed_y_val = this_site_thresh
# #     
# #     X <- sm$X[, -k]
# #     S <- sm$S[[1]][-k, -k]  
# #     off <-  data_for_quant_spline$qunantile_value*0 + (fixed_y_val)
# #     
# #     # fit spline
# #     gam_mod <- gam(data_for_quant_spline$qunantile_value ~ X -1 + offset(off), paraPen = list(X = list(S)))
# #     
# #     tibble(Station = .x$Station[1], sm = list(sm), k, gam_mod = list(gam_mod), fixed_y_val, this_thresh_quantile)
# #     
# #   }, .keep = T) %>%
# #   plyr::rbind.fill() %>%
# #   as_tibble()
# # 
# # 
# # 
# # 
# # 
# # 
# # # interpolate quantuiles (with clim threshold)
# # obs_bulk_models_with_thresh = obs_sites_quantile_estimtes %>%
# #   group_by(Station) %>%
# #   group_map(~{
# #     
# #     print(.x$Station[1])
# #     this_site_thresh = obs_data %>%
# #       filter(Station == .x$Station[1]) %>%
# #       pull(threshold_w_clim) %>%
# #       .[1] %>%
# #       print()
# #     
# #     this_site_data = obs_data %>%
# #       filter(Station == .x$Station[1]) %>%
# #       pull(maxtp)
# #     
# #     #this_thresh_quantile = ecdf(this_site_data)(this_site_thresh) %>% print
# #     this_thresh_quantile = 0.8
# #     
# #     # -- make sure nothing greater than threshold
# #     
# #     data_for_quant_spline = tibble(Station = .x$Station[1],
# #                                    quantile = .x$quantile,
# #                                    qunantile_value = .x$qunantile_value)
# #     
# #     k <- 10 #num knots
# #     knots <- data.frame(x = seq(0, this_thresh_quantile,length=k)) 
# #     sm <- smoothCon(s(quantile,k=k,bs="cr"), data_for_quant_spline, knots=knots)[[1]]
# #     
# #     # value to force spline through
# #     fixed_y_val = this_site_thresh
# #     
# #     X <- sm$X[, -k]
# #     S <- sm$S[[1]][-k, -k]  
# #     off <-  data_for_quant_spline$qunantile_value*0 + (fixed_y_val)
# #     
# #     # fit spline
# #     gam_mod <- gam(data_for_quant_spline$qunantile_value ~ X -1 + offset(off), paraPen = list(X = list(S)))
# #     tibble(Station = .x$Station[1], sm = list(sm), k, gam_mod = list(gam_mod), fixed_y_val, this_thresh_quantile)
# #     
# #   }, .keep = T) %>%
# #   plyr::rbind.fill() %>%
# #   as_tibble()
# # 
# # 
# # 
# # 
# # 
# # # ---- test back transformation
# # # this_site_thresh = obs_data %>%
# # #   filter(Station == "claremorris") %>%
# # #   pull(threshold) %>%
# # #   .[1] 
# # 
# # # this_site_data = obs_data %>%
# # #   filter(Station == "claremorris") %>%
# # #   pull(maxtp)
# # # 
# # # this_thresh_quantile = ecdf(this_site_data)(this_site_thresh)
# # # 
# # # 
# # # obs_sites_quantile_estimtes_test = obs_sites_quantile_estimtes %>%
# # #   filter(Station == "claremorris") %>%
# # #   filter(qunantile_value < this_site_thresh,
# # #          quantile < 0.8)
# # # 
# # # 
# # # data_for_quant_spline = tibble(Station = obs_sites_quantile_estimtes_test$Station[1],
# # #                                quantile = c(obs_sites_quantile_estimtes_test$quantile, 0.8),
# # #                                qunantile_value = c(obs_sites_quantile_estimtes_test$qunantile_value, this_site_thresh))
# # # 
# # # 
# # # data_for_quant_spline %>% 
# # #   ggplot() + 
# # #   geom_point(aes(quantile, qunantile_value))
# # # 
# # # 
# # # 
# # # 
# # 
# # 
# # 
# # 
# # 
# # 
# # obs_bulk_models = readRDS("~/JRSS_organised_code/corrections/bulk_splines_at_obs_sites")
# # 
# # 
# # new_mod = readRDS("corrections/quantile_regression_models/quantile_model_0.01")
# # 
# # 
# # new_mod %>%
# #   filter(Station == "claremorris")
# # 
# # 
# # testing_mod = obs_bulk_models %>%
# #   filter(Station == "strabane")
# # 
# # 
# # testing_mod$Station %>% unique
# # 
# # 
# # quantile_to_estimate_up_to_test = testing_mod$this_thresh_quantile
# # sm_test = testing_mod$sm %>% .[[1]]
# # k_test = testing_mod$k
# # gam_mod_test = testing_mod$gam_mod %>% .[[1]]
# # fixed_y_val_test = testing_mod$fixed_y_val
# # 
# # newdata <- data.frame(x = seq(0, 0.8, length.out = 100))
# # newdata_tmp <- Predict.matrix(sm_test, data.frame(quantile = newdata$x))
# # newdata_tmp <- newdata_tmp[, -k_test]
# # newdata$y_pred_fit1 <- (newdata_tmp %*% coef(gam_mod_test))[, 1] + (fixed_y_val_test)
# # 
# # obs_data %>%
# #   filter(Station == "strabane") %>% pull(threshold)
# # tibble(newdata) %>%
# #   ggplot()+
# #   geom_line(aes(x, y_pred_fit1))+
# #   geom_hline(yintercept = 20.40378)+
# #   geom_point(data = tibble(x = seq(0,0.8, length.out = 20),
# #                            y = obs_data %>%
# #                              filter(Station =="strabane") %>%
# #                              pull(maxtp) %>%
# #                              quantile(seq(0,0.8, length.out = 20)) %>%
# #                              as.numeric), aes(x,y))
# # # 
# # # 
# # # 
# # # 
# # # 
# # # 
# # # 
# # # data_for_back_transformation = obs_data %>%
# # #   filter(Station == 'claremorris') %>%
# # #   filter(maxtp < threshold) %>%
# # #   pull(maxtp)
# # # 
# # # newdata <- data.frame(x = seq(0, 0.8, length.out = length(data_for_back_transformation)))
# # # newdata_tmp <- Predict.matrix(sm_test, data.frame(quantile = newdata$x))
# # # newdata_tmp <- newdata_tmp[, -k_test]
# # # newdata$y_pred_fit1 <- (newdata_tmp %*% coef(gam_mod_test))[, 1] + (fixed_y_val_test)
# # # 
# # # qunant_preds = (newdata_tmp %*% coef(gam_mod_test))[, 1] + (fixed_y_val_test)
# # # 
# # # res = rep(NA, length(data_for_back_transformation))
# # # 
# # # for(i in seq(length(data_for_back_transformation))){
# # #   print(data_for_back_transformation[i])
# # #   
# # #   res[i] = newdata$x[which.min(abs(qunant_preds - data_for_back_transformation[i]))]
# # #   
# # # }
# # # summary(res)
# # # hist(res)
# # #
# # # 
# # # obs_data %>%
# # #   filter(Station == 'katesbridge') %>%
# # #   pull(maxtp) %>% quantile(seq(0,1, length.out = 25))
# # # 
# # # 
# # # newdata %>% 
# # #   ggplot()+
# # #   geom_point(aes(x, y_pred_fit1))+
# # #   geom_point(data = obs_sites_quantile_estimtes %>% 
# # #                filter(Station == "katesbridge") %>%
# # #                filter(quantile < quantile_to_estimate_up_to_test),aes(quantile, qunantile_value), col = 'red')
# # # 
# # # 
# # # # ----
# # # obs_sites_quantile_estimtes %>%
# # #   filter(Station == "katesbridge") %>%
# # #   ggplot()+
# # #   geom_point(aes(quantile, qunantile_value), col = 'red')+
# # #   geom_point(data = obs_data %>%
# # #                filter(Station == "katesbridge") %>%
# # #                summarise(quantiles_to_estimate_bulk, q = quantile(maxtp, quantiles_to_estimate_bulk)),
# # #              aes(quantiles_to_estimate_bulk, q))
# # # 
# # # -------- END OF TEST
# # 
# # 
# # 
# # # 
# # # print("replace 95th quantile with empirical one (threshold for obs. sites)")
# # # # replace 95th quantile with empirical one (threshold for obs. sites)
# # # # this reduces the jump between the bulk and tail model.
# # # # add min observation as 0th quantile
# # # obs_sites_quantile_estimtes = rbind(obs_data %>%
# # #         group_by(Station) %>%
# # #         mutate(qunantile_value = min(maxtp)) %>%
# # #         dplyr::select(Station, Long, Lat, Long.projected, Lat.projected,qunantile_value) %>%
# # #         mutate(quantile = 0) %>% 
# # #         ungroup() %>%
# # #         unique(),
# # #       obs_sites_quantile_estimtes %>%
# # #         filter(quantile != quantile_to_estimate_up_to))
# # # 
# # # 
# # #       # obs_data %>%
# # #       #   dplyr::select(Station, qunantile_value = threshold, Long, Lat, Long.projected, Lat.projected) %>%
# # #       #   mutate(quantile = quantile_to_estimate_up_to) %>%
# # #       #   unique())
# # # 
# # 
# # # 
# # # obs_bulk_models = obs_sites_quantile_estimtes %>%
# # #   group_by(Station) %>%
# # #   group_map(~{
# # #     
# # #     
# # #     this_site_thresh = obs_data %>%
# # #       filter(Station == .x$Station[1]) %>%
# # #       pull(threshold) %>% .[1] 
# # #     
# # #     this_site_maxtp = obs_data %>%
# # #       filter(Station == .x$Station[1]) %>%
# # #       pull(maxtp) 
# # # 
# # #     quantile_to_estimate_up_to = signif(ecdf(this_site_maxtp)(this_site_thresh), 2)
# # # 
# # #     
# # #     k <- 10 #num knots
# # #     knots <- data.frame(x = seq(0, quantile_to_estimate_up_to,length=k)) 
# # #     sm <- smoothCon(s(quantile,k=k,bs="cr"), .x, knots=knots)[[1]]
# # #     
# # #     # value to force spline through
# # #     fixed_y_val = this_site_thresh
# # #     
# # #     X <- sm$X[, -k]
# # #     S <- sm$S[[1]][-k, -k]  
# # #     off <-  .x$qunantile_value*0 + (fixed_y_val)
# # #     
# # #     # fit spline
# # #     gam_mod <- gam(.x$qunantile_value ~ X -1 + offset(off), paraPen = list(X = list(S)))
# # #     
# # #     tibble(Station = .x$Station[1], sm = list(sm), k, gam_mod = list(gam_mod), fixed_y_val, quantile_to_estimate_up_to)
# # #   }, .keep = T) %>%
# # #   plyr::rbind.fill() %>%
# # #   as_tibble()
# # 
# # 
# # # # ========HOW TO making new predictions!
# # # newdata <- data.frame(x = seq(0, quantile_to_estimate_up_to, length.out = 100))
# # # newdata_tmp <- Predict.matrix(sm, data.frame(quantile = newdata$x))
# # # newdata_tmp <- newdata_tmp[, -k]
# # # newdata$y_pred_fit1 <- (newdata_tmp %*% coef(gam_mod))[, 1] + (fixed_y_val)
# # 
# # 
# # #===== old method - doesnt go through threshold
# # # obs_bulk_models = obs_sites_quantile_estimtes %>%
# # #   group_by(Station) %>%
# # #   group_map(~{
# # #     quantile_model = mgcv::gam(qunantile_value~s(quantile), data = .x) 
# # # 
# # #     .x %>% 
# # #       select(Station, Long, Lat, Long.projected, Lat.projected) %>%
# # #       unique() %>%
# # #       mutate(quantile_model = list(quantile_model))
# # #   }, .keep = T) %>%
# # #   plyr::rbind.fill() %>%
# # #   as_tibble()
# # # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # print("save bulk models")
# # # save bulk models
# # obs_bulk_models %>% saveRDS("corrections/bulk_splines_at_obs_sites")
# # obs_bulk_models_with_thresh %>% saveRDS("corrections/bulk_splines_at_obs_sites_w_clim_thresh")
# # 
# # # # plot to check looks ok
# # # tibble(quantile = seq(0.01, marg_thresh_quantile, length.out=100),
# # #        est = obs_bulk_models %>% filter(Station == "mullingar") %>%
# # #          pull(quantile_model) %>%
# # #          .[[1]] %>%
# # #          predict(tibble(quantile = seq(0.01, marg_thresh_quantile, length.out=100)))) %>%
# # #   ggplot()+
# # #   geom_line(aes(quantile, est))+
# # #   geom_point(data = obs_sites_quantile_estimtes %>% filter(Station == "mullingar"), 
# # #              aes(quantile,qunantile_value ))
# # 
# # 
# # 
# # 
# # # 3. ========  ========  Tail Model        ========  ========  ======== 
# # print("model climate tail")
# # # 3a. ---- Climate model
# # clim_data_extreme = clim_data %>%
# #   group_by(id) %>%
# #   mutate(threshold = quantile(maxtp, marg_thresh_quantile),
# #          excess = maxtp - threshold) %>%
# #   filter(excess > 0) %>%
# #   ungroup()
# # 
# # ngll = function(par){
# #   if(par <= 0) return(2^30)
# #   if(par > -1/shape_param) return(2^30)
# #   if(any((1+shape_param*this.dat/par)< 0)) return(2^30)
# #   
# #   -sum(evd::dgpd(x = this.dat, loc=0, scale = par, shape=shape_param, log=T))
# # }
# # 
# # estiamte_scale_fixed_shape = function(x,shape_c){
# #   this.dat <<- x
# #   shape_param <<- shape_c
# #   optim(par = c(0), fn = ngll, method = 'Brent', lower=0, upper = 5)
# # }
# # 
# # # fit scale parameter to each location with constant shap
# # num_sites = clim_data_extreme$id %>% unique()
# # scales = c()
# # loglik_sum = c()
# # 
# # 
# # 
# # 
# # 
# # for(potential_shape in potential_shape_values_climate){
# #   print(potential_shape)
# #   
# #   loglik = c()
# #   for(i in (clim_data_extreme$id %>% unique())){
# #     this_clim_extm_irel = clim_data_extreme %>%
# #       filter(id == i) %>% pull(excess)
# #     
# #     
# #     model_fit = estiamte_scale_fixed_shape(this_clim_extm_irel, potential_shape)
# #     estimated_scale = model_fit$par
# #     
# #     scales = c(scales, estimated_scale)
# #     loglik = c(loglik, model_fit$value)
# #   }
# #   loglik_sum = c(loglik_sum, sum(loglik))
# # }
# # 
# # optimal_shape = potential_shape_values_climate[which.min(loglik_sum)]
# # optimal_shape
# # 
# # 
# # #optimal_shape = -0.2038076
# # print(optimal_shape)
# # 
# # scales = c()
# # for(i in (clim_data_extreme$id %>% unique())){
# #   this_clim_extm_irel = clim_data_extreme %>%
# #     filter(id == i) %>% pull(excess)
# #   
# #   model_fit = estiamte_scale_fixed_shape(this_clim_extm_irel, optimal_shape)
# #   estimated_scale = model_fit$par
# #   scales = c(scales, estimated_scale)
# # }
# # 
# # clim_data_extreme %>% 
# #   dplyr::select(Long, Lat, id) %>%
# #   unique() %>% 
# #   mutate(scale = scales) %>%
# #   ggplot()+
# #   geom_point(aes(Long, Lat, col = scale), size = 2.5)+
# #   geom_sf(data = ireland_sf, alpha = 0, col = 'black')+
# #   ggplot2::scale_color_gradientn(colors = my_pal)+
# #   theme_minimal(12)+
# #   labs(col = expression(sigma))
# # 
# # 
# # clim_data_extreme %>% 
# #   dplyr::select(Long, Lat, Long.projected, Lat.projected, id, threshold) %>%
# #   unique() %>% 
# #   mutate(clim_scale = scales) %>%
# #   write_csv("corrections/clim_scale_grid.csv")
# # 
# # 
# # obs_sites = obs_sites %>%
# #   left_join(clim_data_extreme %>% 
# #               dplyr::select(id) %>%
# #               unique() %>% 
# #               mutate(clim_scale = scales))
# # 
# # obs_sites %>%
# #   ggplot()+
# #   geom_point(aes(Long, Lat, col = clim_scale), size = 2.5)+
# #   geom_sf(data = ireland_sf, alpha = 0, col = 'black')+
# #   ggplot2::scale_color_gradientn(colors = my_pal)+
# #   theme_minimal(12)+
# #   labs(col = expression(sigma))
# # 
# # obs_data = obs_data %>%
# #   left_join(obs_sites)
# # 
# # 
# # obs_data %>% write_csv("corrections/obs_data.csv")
# # 
# # 
# # # 
# # # # 4. ========  ========  Temporally decluster ========  ========
# # # # do i need to do this?
# # # 
# # # 
# # # # 5. ========  ========  Marginal Model ========  ========
# # # # here fit the finalised marginal model
# # # 
# # # 
# # # # 6. ========  ========  Threshold exceedance model ========  ========
# # # 
# # # system.time({
# # #   binom_model <-  mgcv::gam((maxtp>threshold_wo_clim) ~ 1, # same formula as scale 
# # #                             data=all_data,
# # #                             family=binomial,
# # #                             method="REML",
# # #                             control=list(nthreads=4))
# # # })
# # # 
# # # 
# # # 
# # # 
# # # 
# # # 
# # # 
# # # 
# # # #
# # # 
# # # 
# # # 
# # # 
# # # # generate list of dates for each bootstrap dataset
# # # 
# # # 
# # # 
# # # 
# # # 
# # # 
# # # obs_extreme = read_csv("corrections/obs_data.csv") %>%
# # #   mutate(excess = maxtp - threshold) %>%
# # #   filter(excess > 0)
# # # 
# # # 
# # # 
# # # 
# # # evgam::evgam(formula = list(excess ~1,~1), family = 'gpd', data = obs_extreme) %>%
# # #   BIC
# # # 
# # # evgam::evgam(formula = list(excess ~clim_scale,~1), family = 'gpd', data = obs_extreme) %>%
# # #   BIC
# # # 
# # # evgam::evgam(formula = list(excess ~clim_scale+s(co2)+s(year)+s(gmt),~1), family = 'gpd', data = obs_extreme) %>%
# # #   BIC
# # # 
# # # # 3a. ---- Climate model
# # # # 2b. Observational data 
# # # 
# # # 
# # # 
# # # 
# # # 
# # # # 2. Estiamte quantiles from
# # # 
# # # # 3.
# # # 
# # # # 4.
# # # 
# # # # 5.
# # # 
# # # # 6.
# # # 
# # # # 7.
# # # 
# # # #read in all climate data
# # # 
# # # 
# # # # read in all observational data 