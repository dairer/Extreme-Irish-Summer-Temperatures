
library(tidyverse)
marg_mod = 'mod_2'

calc_chi_true = function(marg_mod, yr, tmp){
  set.seed(123456)
  num_samples = 8000
  
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
      id2_exceeds = (unlist(lapply(my_simulations, "[[", which(grid_simulated$Station == id2))) > frechet_val[which(grid_simulated$Station == id2)])
      
      tibble(sites[s,] %>% mutate(chi =sum(id1_exceeds & id2_exceeds)/sum(id1_exceeds) )) %>%
        write_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_model_",marg_mod,"_yr_",yr, "_conditioned_on_",tmp,".csv"),append = T)
    }
  }
}

# job::job({calc_chi_true(marg_mod = "mod_2", yr = 2020, tmp = 30)})
# job::job({calc_chi_true(marg_mod = "mod_2", yr = 2020, tmp = 29)})
# job::job({calc_chi_true(marg_mod = "mod_2", yr = 2020, tmp = 28)})
# job::job({calc_chi_true(marg_mod = "mod_2", yr = 1942, tmp = 30)})
# job::job({calc_chi_true(marg_mod = "mod_2", yr = 1942, tmp = 29)})
# job::job({calc_chi_true(marg_mod = "mod_2", yr = 1942, tmp = 28)})



calc_chi_bts = function(marg_mod, yr, tmp, bts_seq){
  num_samples = 8000 --
  set.seed(123456)

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
        id2_exceeds = (unlist(lapply(my_simulations, "[[", which(grid_simulated$Station == id2))) > frechet_val[which(grid_simulated$Station == id2)])
        
        
        tibble(bts = bts, sites[s,] %>% mutate(chi =sum(id1_exceeds & id2_exceeds)/sum(id1_exceeds) )) %>%
          write_csv(paste0("output/simulations/simulation_summary/bootstrap_chi_data_scale_model_",marg_mod,"_bts_yr_",yr, "_conditioned_on_",tmp, ".csv"),append = T)
      }
    }
  }
}


# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 2020, tmp = 30, bts_seq = seq(1,100))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 2020, tmp = 30, bts_seq = seq(101,200))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 2020, tmp = 30, bts_seq = seq(201,300))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 2020, tmp = 29, bts_seq = seq(1,100))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 2020, tmp = 29, bts_seq = seq(101,200))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 2020, tmp = 29, bts_seq = seq(201,300))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 2020, tmp = 28, bts_seq = seq(1,100))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 2020, tmp = 28, bts_seq = seq(101,200))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 2020, tmp = 28, bts_seq = seq(201,300))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 1942, tmp = 30, bts_seq = seq(1,100))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 1942, tmp = 30, bts_seq = seq(101,200))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 1942, tmp = 30, bts_seq = seq(201,300))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 1942, tmp = 29, bts_seq = seq(1,100))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 1942, tmp = 29, bts_seq = seq(101,200))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 1942, tmp = 29, bts_seq = seq(201,300))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 1942, tmp = 28, bts_seq = seq(1,100))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 1942, tmp = 28, bts_seq = seq(101,200))})
# job::job({calc_chi_bts(marg_mod = "mod_2", yr = 1942, tmp = 28, bts_seq = seq(201,300))})


# # ------ PLOT MODELS
marg_mod = 'mod_2'
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

chi_true_summary = chi_true %>%
  mutate(dist_bin = cut(distance, breaks=seq(0, 4, length.out = 20))) %>%
  group_by(dist_bin) %>%
  summarise(year, temp, chi, distance = mean(distance))  %>%
  group_by(distance, year, temp) %>%
  drop_na() %>%
  summarise(mn = mean(chi))


# # ----- Uncondition
prob_2020_28 = read_csv(paste0("output/prob_extreme_temp_imp_samp_",marg_mod,".csv"),
                        col_names = c('temp', 'p_1942', 'p_2020'))  %>%
  filter(temp == 28) %>% pull(p_2020)

prob_2020_29  = read_csv(paste0("output/prob_extreme_temp_imp_samp_",marg_mod,".csv"),
                         col_names = c('temp', 'p_1942', 'p_2020'))  %>%
  filter(temp == 29) %>% pull(p_2020)

prob_2020_30 = read_csv(paste0("output/prob_extreme_temp_imp_samp_",marg_mod,".csv"),
                           col_names = c('temp', 'p_1942', 'p_2020'))  %>%
  filter(temp == 30) %>% pull(p_2020)

prob_1942_28 = read_csv(paste0("output/prob_extreme_temp_imp_samp_",marg_mod,".csv"),
                        col_names = c('temp', 'p_1942', 'p_2020'))  %>%
  filter(temp == 28) %>% pull(p_1942)

prob_1942_29 = read_csv(paste0("output/prob_extreme_temp_imp_samp_",marg_mod,".csv"),
                        col_names = c('temp', 'p_1942', 'p_2020'))  %>%
  filter(temp == 29) %>% pull(p_1942)

prob_1942_30 = read_csv(paste0("output/prob_extreme_temp_imp_samp_",marg_mod,".csv"),
                        col_names = c('temp', 'p_1942', 'p_2020'))  %>%
  filter(temp == 30) %>% pull(p_1942)

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

dat_for_rat = chi_true_summary %>%
  ungroup() %>%
  filter(distance < 1.16, distance > 1) %>%
  filter(lab == 'unconditioned')  %>%
  dplyr::select(year, mn, temp) 

dat_for_rat[dat_for_rat$year == 2020 & dat_for_rat$temp == "28°C",]$mn/dat_for_rat[dat_for_rat$year == 1942 & dat_for_rat$temp == "28°C",]$mn
dat_for_rat[dat_for_rat$year == 2020 & dat_for_rat$temp == "29°C",]$mn/dat_for_rat[dat_for_rat$year == 1942 & dat_for_rat$temp == "29°C",]$mn
dat_for_rat[dat_for_rat$year == 2020 & dat_for_rat$temp == "30°C",]$mn/dat_for_rat[dat_for_rat$year == 1942 & dat_for_rat$temp == "30°C",]$mn

# 2.8, 3.5 and 4.7
plt = chi_bts_summet %>%
  ggplot()+
  geom_ribbon(aes(x = distance*100, ymin = lower, ymax = upper, fill = year), alpha = 0.25)+
  geom_smooth(data = chi_true_summary, aes(x = distance*100, y = mn,col = year, linetype = year), se=F)+
  facet_grid(lab~temp,scale = 'free')+
  labs(x = "Distance (km)",
       y = expression(chi[o]),
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
