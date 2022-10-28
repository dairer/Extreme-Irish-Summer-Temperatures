
library(tidyverse)
marg_mod = 'mod_4'

calc_chi_true = function(marg_mod, yr, tmp){
  #num_samples = 8000 --- model 4
  set.seed(123456)
  num_samples = 4000 #--- model 1+2
  
  sites = read_csv("data/processed/obs_pairs_with_dist.csv") %>% sample_n(num_samples)

  obs_sites = read_csv("data/processed/obs_data.csv") %>%
    dplyr::select(Station, Long.projected, Lat.projected) %>%
    unique() %>%
    left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% dplyr::select(Station, dist_sea))
  
  grid_simulated = as.data.frame(read_csv("data/processed/obs_grid_simulated_on.csv")) %>%
    left_join(obs_sites) %>% as.tibble()
  
  frechet_val = grid_simulated %>%
    left_join(read_csv(paste0("output/obs_sites_extreme_temps_frechet_scale_",marg_mod,".csv")) %>% filter(temp == tmp, year == yr)) %>% 
    pull(frechet_value)
  
  
  if(file.exists(paste0("output/simulations/rescaled_simulations_on_obs_grid/true/",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp))){
    my_simulations = readRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/true/",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp))
    
    for(s in seq(nrow(sites))){
      
      if((s %% 100) == 0){
        print(s)
      }
      
      id1 = sites[s,]$V1
      id2 = sites[s,]$V2
      
      # ---- to get all sims at loc with id v1 ---> unlist(lapply(my_simulations, "[[", sites[s,]$V1))
      id1_exceeds = (unlist(lapply(my_simulations, "[[", which(grid_simulated$Station == id1))) > frechet_val[which(grid_simulated$Station == id1)])
      id2_exceeds = (unlist(lapply(my_simulations, "[[", which(grid_simulated$Station == id2))) > frechet_val[which(grid_simulated$Station == id1)])
      
      tibble(sites[s,] %>% mutate(chi =sum(id1_exceeds & id2_exceeds)/sum(id1_exceeds) )) %>%
        write_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_model_",marg_mod,"_yr_",yr, "_conditioned_on_",tmp,".csv"),append = T)
    }
  }
}


job::job({calc_chi_true(marg_mod = "mod_4", yr = 2020, tmp = 30)})
job::job({calc_chi_true(marg_mod = "mod_4", yr = 2020, tmp = 29)})
job::job({calc_chi_true(marg_mod = "mod_4", yr = 2020, tmp = 28)})
job::job({calc_chi_true(marg_mod = "mod_4", yr = 1942, tmp = 30)})
job::job({calc_chi_true(marg_mod = "mod_4", yr = 1942, tmp = 29)})
job::job({calc_chi_true(marg_mod = "mod_4", yr = 1942, tmp = 28)})

job::job({calc_chi_true(marg_mod = "mod_1", yr = 2020, tmp = 30)})
job::job({calc_chi_true(marg_mod = "mod_1", yr = 2020, tmp = 29)})
job::job({calc_chi_true(marg_mod = "mod_1", yr = 2020, tmp = 28)})
job::job({calc_chi_true(marg_mod = "mod_1", yr = 1942, tmp = 30)})
job::job({calc_chi_true(marg_mod = "mod_1", yr = 1942, tmp = 29)})
job::job({calc_chi_true(marg_mod = "mod_1", yr = 1942, tmp = 28)})

job::job({calc_chi_true(marg_mod = "mod_2", yr = 2020, tmp = 30)})
job::job({calc_chi_true(marg_mod = "mod_2", yr = 2020, tmp = 29)})
job::job({calc_chi_true(marg_mod = "mod_2", yr = 2020, tmp = 28)})
job::job({calc_chi_true(marg_mod = "mod_2", yr = 1942, tmp = 30)})
job::job({calc_chi_true(marg_mod = "mod_2", yr = 1942, tmp = 29)})
job::job({calc_chi_true(marg_mod = "mod_2", yr = 1942, tmp = 28)})



# read_csv("output/simulations/simulation_summary/chi_data_scale_clim_grid_model_mod_4_yr_1942_min_temp_28_conditioned_on_27.csv",
#          col_names = c('s1', 's2', 'dist', 'chi')) %>%
#   drop_na()
# 
# 
# # job::job({calc_chi_true(marg_mod = "mod_4", yr = 2020, tmp = 30)})
# # job::job({calc_chi_true(marg_mod = "mod_4", yr = 2020, tmp = 29)})
# # job::job({calc_chi_true(marg_mod = "mod_4", yr = 2020, tmp = 28)})
# # job::job({calc_chi_true(marg_mod = "mod_4", yr = 1942, tmp = 30)})
# # job::job({calc_chi_true(marg_mod = "mod_4", yr = 1942, tmp = 29)})
# # job::job({calc_chi_true(marg_mod = "mod_4", yr = 1942, tmp = 28)})
# 



calc_chi_bts = function(marg_mod, yr, tmp, bts_seq){
  #num_samples = 8000 -- model 4
  num_samples = 500 #-- model 1+2
  set.seed(123456)
  
  
  # bts_seq = bts_seq[!(bts_seq %in% (read_csv(paste0("output/simulations/simulation_summary/bootstrap_chi_data_scale_model_",marg_mod,"_bts_yr_",yr, "_conditioned_on_",tmp, ".csv"),
  #                                            col_names = c('bts', 's1', 's2', 'h', 'chi')) %>%
  #                                     pull(bts) %>% unique))]


  sites = read_csv("data/processed/obs_pairs_with_dist.csv") %>% sample_n(num_samples)

  obs_sites = read_csv("data/processed/obs_data.csv") %>%
    dplyr::select(Station, Long.projected, Lat.projected) %>%
    unique() %>%
    left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% dplyr::select(Station, dist_sea))

  grid_simulated = as.data.frame(read_csv("data/processed/obs_grid_simulated_on.csv")) %>%
    left_join(obs_sites) %>% as.tibble()


  for(bts in bts_seq){
    print("new bootstrap ... ")
    print(bts)
    
    if(file.exists(paste0("output/simulations/rescaled_simulations_on_obs_grid/bootstraps/bootstrap_",bts ,"_model_", marg_mod,"_yr_", yr,"_min_temp_conditioned_on_", tmp))){
      
      print("yay bts")
      frechet_val = grid_simulated %>%
        left_join(read_csv(paste0("output/extreme_frechet_values_bts/obs_sites_extreme_temps_bootstraps_frechet_scale_", marg_mod,"_bts_", bts ,".csv"),
                           col_names = c('year', 'Station', 'bts', 'temp', 'frechet_value')) %>% filter(temp == tmp, year == yr)) %>%
        pull(frechet_value)
  
      my_simulations = readRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/bootstraps/bootstrap_",bts ,"_model_", marg_mod,"_yr_", yr,"_min_temp_conditioned_on_", tmp))
      
      for(s in seq(nrow(sites))){
        
        id1 = sites[s,]$V1
        id2 = sites[s,]$V2
        
        # ---- to get all sims at loc with id v1 ---> unlist(lapply(my_simulations, "[[", sites[s,]$V1))
        id1_exceeds = (unlist(lapply(my_simulations, "[[", which(grid_simulated$Station == id1))) > frechet_val[which(grid_simulated$Station == id1)])
        id2_exceeds = (unlist(lapply(my_simulations, "[[", which(grid_simulated$Station == id2))) > frechet_val[which(grid_simulated$Station == id1)])
        
        
        tibble(bts = bts, sites[s,] %>% mutate(chi =sum(id1_exceeds & id2_exceeds)/sum(id1_exceeds) )) %>%
          write_csv(paste0("output/simulations/simulation_summary/bootstrap_chi_data_scale_model_",marg_mod,"_bts_yr_",yr, "_conditioned_on_",tmp, ".csv"),append = T)
      }
    }
  }
}



job::job({calc_chi_bts(marg_mod = "mod_1", yr = 2020, tmp = 29, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_1", yr = 2020, tmp = 29, bts_seq = seq(401,500))})
job::job({calc_chi_bts(marg_mod = "mod_1", yr = 2020, tmp = 28, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_1", yr = 2020, tmp = 28, bts_seq = seq(401,500))})
job::job({calc_chi_bts(marg_mod = "mod_1", yr = 2020, tmp = 30, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_1", yr = 2020, tmp = 30, bts_seq = seq(401,500))})
job::job({calc_chi_bts(marg_mod = "mod_1", yr = 1942, tmp = 29, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_1", yr = 1942, tmp = 29, bts_seq = seq(401,500))})
job::job({calc_chi_bts(marg_mod = "mod_1", yr = 1942, tmp = 28, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_1", yr = 1942, tmp = 28, bts_seq = seq(401,500))})
job::job({calc_chi_bts(marg_mod = "mod_1", yr = 1942, tmp = 30, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_1", yr = 1942, tmp = 30, bts_seq = seq(401,500))})

job::job({calc_chi_bts(marg_mod = "mod_2", yr = 2020, tmp = 29, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_2", yr = 2020, tmp = 29, bts_seq = seq(401,500))})
job::job({calc_chi_bts(marg_mod = "mod_2", yr = 2020, tmp = 28, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_2", yr = 2020, tmp = 28, bts_seq = seq(401,500))})
job::job({calc_chi_bts(marg_mod = "mod_2", yr = 2020, tmp = 30, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_2", yr = 2020, tmp = 30, bts_seq = seq(401,500))})



job::job({calc_chi_bts(marg_mod = "mod_2", yr = 1942, tmp = 29, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_2", yr = 1942, tmp = 29, bts_seq = seq(401,500))})
job::job({calc_chi_bts(marg_mod = "mod_2", yr = 1942, tmp = 28, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_2", yr = 1942, tmp = 28, bts_seq = seq(401,500))})
job::job({calc_chi_bts(marg_mod = "mod_2", yr = 1942, tmp = 30, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_2", yr = 1942, tmp = 30, bts_seq = seq(401,500))})

# # job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 30, bts_seq = seq(1,100))})
# # job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 30, bts_seq = seq(101,200))})
# # job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 30, bts_seq = seq(201,300))})
# # job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 30, bts_seq = seq(301,400))})
# # job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 30, bts_seq = seq(401,500))})
# # 
# # job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 29, bts_seq = seq(1,100))})
# # job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 29, bts_seq = seq(101,200))})
# # job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 29, bts_seq = seq(201,300))})
# # job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 29, bts_seq = seq(301,400))})
# # job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 29, bts_seq = seq(401,500))})
# # 
# # job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 28, bts_seq = seq(1,100))})
# # job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 28, bts_seq = seq(101,200))})
# # job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 28, bts_seq = seq(201,300))})
# # job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 28, bts_seq = seq(301,400))})
# # job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 28, bts_seq = seq(401,500))})


job::job({calc_chi_bts(marg_mod = "mod_4", yr = 1942, tmp = 30, bts_seq = seq(1,100))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 1942, tmp = 30, bts_seq = seq(101,200))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 1942, tmp = 30, bts_seq = seq(201,300))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 1942, tmp = 30, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 1942, tmp = 30, bts_seq = seq(401,500))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 1942, tmp = 29, bts_seq = seq(1,100))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 1942, tmp = 29, bts_seq = seq(101,200))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 1942, tmp = 29, bts_seq = seq(201,300))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 1942, tmp = 29, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 1942, tmp = 29, bts_seq = seq(401,500))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 1942, tmp = 28, bts_seq = seq(1,100))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 1942, tmp = 28, bts_seq = seq(101,200))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 1942, tmp = 28, bts_seq = seq(201,300))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 1942, tmp = 28, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 1942, tmp = 28, bts_seq = seq(401,500))})

job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 30, bts_seq = seq(1,100))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 30, bts_seq = seq(101,200))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 30, bts_seq = seq(201,300))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 30, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 30, bts_seq = seq(401,500))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 29, bts_seq = seq(1,100))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 29, bts_seq = seq(101,200))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 29, bts_seq = seq(201,300))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 29, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 29, bts_seq = seq(401,500))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 28, bts_seq = seq(1,100))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 28, bts_seq = seq(101,200))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 28, bts_seq = seq(201,300))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 28, bts_seq = seq(301,400))})
job::job({calc_chi_bts(marg_mod = "mod_4", yr = 2020, tmp = 28, bts_seq = seq(401,500))})







# # ------ PLOT MODELS


marg_mod = 'mod_4'


chi_bts = rbind(read_csv(paste0("output/simulations/simulation_summary/bootstrap_chi_data_scale_model_", marg_mod,"_bts_yr_1942_conditioned_on_28.csv"),
                         col_names = c('bts', 's1', 's2', 'distance', 'chi')) %>% mutate(temp = "28°C", year = '1942'),
                read_csv(paste0("output/simulations/simulation_summary/bootstrap_chi_data_scale_model_", marg_mod,"_bts_yr_1942_conditioned_on_29.csv"),
                         col_names = c('bts', 's1', 's2', 'distance', 'chi')) %>% mutate(temp = "29°C", year = '1942'),
                read_csv(paste0("output/simulations/simulation_summary/bootstrap_chi_data_scale_model_", marg_mod,"_bts_yr_1942_conditioned_on_30.csv"),
                         col_names = c('bts', 's1', 's2', 'distance', 'chi')) %>% mutate(temp = "30°C", year = '1942'),
                read_csv(paste0("output/simulations/simulation_summary/bootstrap_chi_data_scale_model_", marg_mod,"_bts_yr_2020_conditioned_on_28.csv"),
                         col_names = c('bts', 's1', 's2', 'distance', 'chi')) %>% mutate(temp = "28°C", year = '2020'),
                read_csv(paste0("output/simulations/simulation_summary/bootstrap_chi_data_scale_model_", marg_mod,"_bts_yr_2020_conditioned_on_29.csv"),
                         col_names = c('bts', 's1', 's2', 'distance', 'chi')) %>% mutate(temp = "29°C", year = '2020'),
                read_csv(paste0("output/simulations/simulation_summary/bootstrap_chi_data_scale_model_", marg_mod,"_bts_yr_2020_conditioned_on_30.csv"),
                         col_names = c('bts', 's1', 's2', 'distance', 'chi')) %>% mutate(temp = "30°C", year = '2020'))

chi_bts_summety = chi_bts %>%
  mutate(dist_bin = cut(distance, breaks=seq(0, 4, length.out = 20))) %>%
  group_by(dist_bin) %>%
  summarise(bts, year, temp, chi, distance = mean(distance))  %>%
  group_by(distance, year, temp, bts) %>%
  drop_na() %>%
  summarise(mn = mean(chi)) %>%
  group_by(distance, year, temp)  %>%
  summarise(upper = quantile(mn, 0.975),
            lower = quantile(mn, 0.025))

chi_true = rbind(read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_model_", marg_mod,"_yr_1942_conditioned_on_28.csv"),
                          col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "28°C", year = '1942'),
                 read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_model_", marg_mod,"_yr_1942_conditioned_on_29.csv"),
                          col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "29°C", year = '1942'),
                 read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_model_", marg_mod,"_yr_1942_conditioned_on_30.csv"),
                          col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "30°C", year = '1942'),
                 read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_model_", marg_mod,"_yr_2020_conditioned_on_28.csv"),
                          col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "28°C", year = '2020'),
                 read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_model_", marg_mod,"_yr_2020_conditioned_on_29.csv"),
                          col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "29°C", year = '2020'),
                 read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_model_", marg_mod,"_yr_2020_conditioned_on_30.csv"),
                          col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "30°C", year = '2020'))
# 
chi_true_summary = chi_true %>%
  mutate(dist_bin = cut(distance, breaks=seq(0, 4, length.out = 20))) %>%
  group_by(dist_bin) %>%
  summarise(year, temp, chi, distance = mean(distance))  %>%
  group_by(distance, year, temp) %>%
  drop_na() %>%
  summarise(mn = mean(chi))
# 
# 
# # ----- Uncondition
prob_2020_28 = read_csv(paste0("output/prob_observing_T_anywhere_clim_", marg_mod,".csv"))  %>%
  filter(temps == 28) %>% pull(prob_2020)

prob_2020_29 =read_csv(paste0("output/prob_observing_T_anywhere_clim_", marg_mod,".csv")) %>%
  filter(temps == 29) %>% pull(prob_2020)

prob_2020_30 = read_csv(paste0("output/prob_observing_T_anywhere_clim_", marg_mod,".csv"))  %>%
  filter(temps == 30) %>% pull(prob_2020)


prob_1942_28 = read_csv(paste0("output/prob_observing_T_anywhere_clim_", marg_mod,".csv"))  %>%
  filter(temps == 28) %>% pull(prob_1942)

prob_1942_29 = read_csv(paste0("output/prob_observing_T_anywhere_clim_", marg_mod,".csv"))  %>%
  filter(temps == 29) %>% pull(prob_1942)

prob_1942_30 = read_csv(paste0("output/prob_observing_T_anywhere_clim_", marg_mod,".csv"))  %>%
  filter(temps == 30) %>% pull(prob_1942)

bts_unconditioned = chi_bts
true_unconditioned = chi_true

bts_unconditioned = rbind(bts_unconditioned %>% filter(year == 1942, temp == '28°C') %>% mutate(chi = prob_1942_28*chi),
                          bts_unconditioned %>% filter(year == 1942, temp == '29°C') %>% mutate(chi = prob_1942_29*chi),
                          bts_unconditioned %>% filter(year == 1942, temp == '30°C') %>% mutate(chi = prob_1942_30*chi),
                          bts_unconditioned %>% filter(year == 2020, temp == '28°C') %>% mutate(chi = prob_2020_28*chi),
                          bts_unconditioned %>% filter(year == 2020, temp == '29°C') %>% mutate(chi = prob_2020_29*chi),
                          bts_unconditioned %>% filter(year == 2020, temp == '30°C') %>% mutate(chi = prob_2020_30*chi))


true_unconditioned = rbind(true_unconditioned %>% filter(year == 1942, temp == '28°C') %>% mutate(chi = prob_1942_28*chi),
                           true_unconditioned %>% filter(year == 1942, temp == '29°C') %>% mutate(chi = prob_1942_29*chi),
                           true_unconditioned %>% filter(year == 1942, temp == '30°C') %>% mutate(chi = prob_1942_30*chi),
                           true_unconditioned %>% filter(year == 2020, temp == '28°C') %>% mutate(chi = prob_2020_28*chi),
                           true_unconditioned %>% filter(year == 2020, temp == '29°C') %>% mutate(chi = prob_2020_29*chi),
                           true_unconditioned %>% filter(year == 2020, temp == '30°C') %>% mutate(chi = prob_2020_30*chi))




bts_dat = rbind(chi_bts %>%
                  mutate(lab = "conditioned"),
                bts_unconditioned %>%
                  mutate(lab = "unconditioned"))


true_data = rbind(chi_true %>%
                    mutate(lab = "conditioned"),
                  true_unconditioned %>%
                    mutate(lab = "unconditioned"))


chi_true_summary = true_data %>%
  mutate(dist_bin = cut(distance, breaks=seq(0, 4, length.out = 20))) %>%
  group_by(dist_bin) %>%
  summarise(year, temp, chi,lab, distance = mean(distance))  %>%
  group_by(distance, year, temp,lab) %>%
  drop_na() %>%
  summarise(mn = mean(chi))


chi_bts_summet = bts_dat %>%
  mutate(dist_bin = cut(distance, breaks=seq(0, 4, length.out = 20))) %>%
  group_by(dist_bin) %>%
  summarise(bts, year, temp, chi, lab,distance = mean(distance))  %>%
  group_by(distance, year, temp, bts, lab) %>%
  drop_na() %>%
  summarise(mn = mean(chi)) %>%
  group_by(distance, year, temp, lab)  %>%
  summarise(upper = quantile(mn, 0.975),
            lower = quantile(mn, 0.025))



plt = chi_bts_summet %>%
  ggplot()+
  geom_ribbon(aes(x = distance*100, ymin = lower, ymax = upper, fill = year), alpha = 0.25)+
  geom_smooth(data = chi_true_summary, aes(x = distance*100, y = mn,col = year, linetype = year), se=F)+
  facet_grid(lab~temp,scale = 'free')+
  labs(x = "Distance (Km)",
       y = expression(chi[o]^D),
       col = "Year",
       shape = "Year",
       linetype = "Year")+
  xlim(0, 375)+
  theme_minimal(12)+
  theme(axis.text.x = element_text(angle = 30),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = 'none',
        strip.text.y = element_text(size=0))


ggsave(paste0("output/figs/chi_data_scale_mod_",marg_mod,".pdf"), height = 4, width = 7)


library(ggthemes)
library(extrafont)
library(remotes)
remotes::install_version("Rttf2pt1", version = "1.3.8")
extrafont::font_import()
# Load fonts
loadfonts(quiet = T)
fonts()

plt = chi_bts_summet %>%
  filter(lab == 'unconditioned') %>%
  ggplot()+
  geom_ribbon(aes(x = distance*100, ymin = lower, ymax = upper, fill = year), alpha = 0.25)+
  geom_smooth(data = chi_true_summary %>% filter(lab == 'unconditioned'), aes(x = distance*100, y = mn,col = year, linetype = year), se=F)+
  facet_wrap(~temp)+
  labs(x = "Distance (Km)",
       y = expression(chi[o]^D),
       col = "Year",
       shape = "Year",
       linetype = "Year")+
  xlim(0, 375)+
  theme_minimal(12)+
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text(angle = 30),
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    strip.text.y = element_text(size=0),
    legend.position = 'none',
    plot.background = element_rect(fill = "transparent", color = NA), 
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    axis.text = element_text(color = 'white'),
    axis.title = element_text(color = 'white'),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    axis.line = element_line(color = 'lightgrey'))

ggsave(plot = plt, filename = "chi_data_scale.png", bg = 'transparent', width = 8, height = 3)
#   



