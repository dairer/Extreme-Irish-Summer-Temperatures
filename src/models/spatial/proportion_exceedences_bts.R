rm(list=ls())
library(tidyverse)
library(scales)
library(latex2exp)

calc_prop_exceedence_bts = function(marg_mod, tmp_conditioned_on, bts_range){

  setwd("~/Extreme-Irish-Summer-Temperatures/")
  
  obs_sites = read_csv("data/processed/obs_data.csv") %>%
    dplyr::select(Station, Long.projected, Lat.projected) %>% 
    unique
  
  for(bts in bts_range){
    print(paste0("Fitting bootstrap number: ", bts))
    yr = 2020

    if(file.exists(paste0("output/simulations/rescaled_simulations_on_obs_grid/bootstraps/bootstrap_", bts, "_model_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_", tmp_conditioned_on))){
      
      my_simulations  = readRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/bootstraps/bootstrap_", bts, "_model_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_", tmp_conditioned_on))
      prop_ex_2020 = c()
      
      for(tmp in c(tmp_conditioned_on)){
      
        frechet_val = read_csv("data/processed/obs_grid_simulated_on.csv") %>%
          left_join(obs_sites) %>%
          left_join(read_csv(paste0("output/extreme_frechet_values_bts/obs_sites_extreme_temps_bootstraps_frechet_scale_",marg_mod,"_bts_",bts, ".csv"),
                             col_names = c('year', 'Station', 'bts', 'temp', 'frechet_value')) %>% filter(temp == tmp, year == yr)) %>%
          pull(frechet_value)
        
        frechet_val[frechet_val == -Inf] = Inf
        
        num_exceed_tmp = my_simulations %>%
          map(~{ sum(.x > frechet_val)}) %>%
          unlist()
        
        prop_ex_2020 = c(prop_ex_2020, mean(num_exceed_tmp[num_exceed_tmp>0]/182))
      }
      
      
      yr = 1942
      my_simulations  = readRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/bootstraps/bootstrap_", bts, "_model_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_", tmp_conditioned_on))
      prop_ex_1942 = c()
      
      for(tmp in c(tmp_conditioned_on)){
        frechet_val = read_csv("data/processed/obs_grid_simulated_on.csv") %>%
          left_join(obs_sites) %>%
          left_join(read_csv(paste0("output/extreme_frechet_values_bts/obs_sites_extreme_temps_bootstraps_frechet_scale_",marg_mod,"_bts_",bts, ".csv"),
                             col_names = c('year', 'Station', 'bts', 'temp', 'frechet_value')) %>% filter(temp == tmp, year == yr)) %>%
          pull(frechet_value)
        frechet_val[frechet_val == -Inf] = Inf
        
        num_exceed_tmp = my_simulations %>% 
          map(~{ sum(.x > frechet_val)}) %>%
          unlist()
        
        prop_ex_1942 = c(prop_ex_1942, mean(num_exceed_tmp[num_exceed_tmp>0]/182))
      }
      
      tibble(temp = c(tmp_conditioned_on), prop_ex_1942, prop_ex_2020, bts = bts) %>%
        write_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_bootstrap_",marg_mod,"_conditioned_on_", tmp_conditioned_on, ".csv"), append = T)
    }
  }
  
}


tmp_conditioned_on = 35

job::job({calc_prop_exceedence_bts('mod_2',tmp_conditioned_on, seq(1,300))})



# ===== Plot results

marg_mod = 'mod_2'
bts_dat = rbind(read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_bootstrap_",marg_mod ,"_conditioned_on_27.csv"),
                         col_names = c('temp', 'prop_ex_1942', 'prop_ex_2020', 'bts')),
                read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_bootstrap_",marg_mod ,"_conditioned_on_28.csv"),
                       col_names = c('temp', 'prop_ex_1942', 'prop_ex_2020', 'bts')),
                read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_bootstrap_",marg_mod ,"_conditioned_on_29.csv"),
                       col_names = c('temp', 'prop_ex_1942', 'prop_ex_2020', 'bts')),
                read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_bootstrap_",marg_mod ,"_conditioned_on_30.csv"),
                      col_names = c('temp', 'prop_ex_1942', 'prop_ex_2020', 'bts')),
                read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_bootstrap_",marg_mod ,"_conditioned_on_31.csv"),
                      col_names = c('temp', 'prop_ex_1942', 'prop_ex_2020', 'bts')),
                read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_bootstrap_",marg_mod ,"_conditioned_on_32.csv"),
                      col_names = c('temp', 'prop_ex_1942', 'prop_ex_2020', 'bts')),
                read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_bootstrap_",marg_mod ,"_conditioned_on_33.csv"),
                      col_names = c('temp', 'prop_ex_1942', 'prop_ex_2020', 'bts')),
                read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_bootstrap_",marg_mod ,"_conditioned_on_34.csv"),
                      col_names = c('temp', 'prop_ex_1942', 'prop_ex_2020', 'bts')),
                read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_bootstrap_",marg_mod ,"_conditioned_on_35.csv"),
                      col_names = c('temp', 'prop_ex_1942', 'prop_ex_2020', 'bts')))


bts_dat = bts_dat %>%
  group_by(temp) %>%
  summarise(prop_ex_2020_u = quantile(prop_ex_2020, 0.975, na.rm = T),
            prop_ex_2020_l = quantile(prop_ex_2020, 0.025, na.rm = T),
            prop_ex_1942_u = quantile(prop_ex_1942, 0.975, na.rm = T),
            prop_ex_1942_l = quantile(prop_ex_1942, 0.025, na.rm = T))


true_dat = rbind(read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_true_",marg_mod ,"_temp_condtioned_on_27.csv")),
                read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_true_",marg_mod ,"_temp_condtioned_on_28.csv")),
                read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_true_",marg_mod ,"_temp_condtioned_on_29.csv")),
                read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_true_",marg_mod ,"_temp_condtioned_on_30.csv")),
                read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_true_",marg_mod ,"_temp_condtioned_on_31.csv")),
                read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_true_",marg_mod ,"_temp_condtioned_on_32.csv")),
                read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_true_",marg_mod ,"_temp_condtioned_on_33.csv")),
                read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_true_",marg_mod ,"_temp_condtioned_on_34.csv")),
                read_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_true_",marg_mod ,"_temp_condtioned_on_35.csv")))


prob_observing_T = read_csv(paste0("output/prob_extreme_temp_imp_samp_",marg_mod,".csv"),
         col_names = c('temp', 'p_1942', 'p_2020'))%>%
  pivot_longer(-temp) %>%
  mutate(year = str_remove(name,"p_")) 

prob_observing_T$value = prob_observing_T$value

prob_observing_T = rbind(prob_observing_T %>% filter(name == 'p_1942') %>% 
                           mutate(typ = 'actual') %>% dplyr::select(-name),
                         prob_observing_T %>% filter(name == 'p_2020') %>% 
                           mutate(typ = 'actual') %>% dplyr::select(-name)) %>%
  mutate(year = as.numeric(year)) %>%
  rename(unconditional_factor = value)


bts_dat = bts_dat %>%
  pivot_longer(-temp) %>%
  mutate(year = str_remove(name,"prop_ex_")) %>%
  mutate(year = str_remove(year,"_u")) %>%
  mutate(year = str_remove(year,"_l"))


bts_dat = rbind(bts_dat %>% filter(name == 'prop_ex_2020_u') %>% mutate(typ = 'upper') %>% dplyr::select(-name),
                bts_dat %>% filter(name == 'prop_ex_2020_l') %>% mutate(typ = 'lower') %>% dplyr::select(-name),
                bts_dat %>% filter(name == 'prop_ex_1942_u') %>% mutate(typ = 'upper') %>% dplyr::select(-name),
                bts_dat %>% filter(name == 'prop_ex_1942_l') %>% mutate(typ = 'lower') %>% dplyr::select(-name)) %>%
  mutate(year = as.numeric(year))


bts_dat = bts_dat %>%
  left_join(prob_observing_T %>% dplyr::select(-typ) %>% unique, by= c('year', 'temp'))

true_dat = true_dat %>%
  pivot_longer(-temp) %>%
  mutate(year = str_remove(name,"prop_ex_")) 
true_dat$year = as.numeric(true_dat$year)
true_dat = true_dat %>%
  left_join(prob_observing_T)


bts_dat$unconditioned = bts_dat$value*bts_dat$unconditional_factor
true_dat$unconditioned = true_dat$value*true_dat$unconditional_factor


bts_dat = rbind(bts_dat %>%
  dplyr::select(temp, value, year, typ) %>%
  pivot_wider(names_from = c(year, typ), values_from = value) %>%
  mutate(lab = 'conditioned'),
  bts_dat %>%
    dplyr::select(temp, unconditioned, year, typ) %>%
    pivot_wider(names_from = c(year, typ), values_from = unconditioned) %>%
    mutate(lab = 'unconditioned'))

true_dat = rbind(true_dat %>%
                  dplyr::select(temp, value, year, typ) %>%
                  pivot_wider(names_from = c(year, typ), values_from = value) %>%
                  mutate(lab = 'conditioned'),
                 true_dat %>%
                  dplyr::select(temp, unconditioned, year, typ) %>%
                  pivot_wider(names_from = c(year, typ), values_from = unconditioned) %>%
                  mutate(lab = 'unconditioned'))


plt = gridExtra::grid.arrange(bts_dat %>% 
  filter(lab == 'conditioned') %>%
  ggplot()+
  geom_ribbon(aes(temp, ymin = `2020_lower`, ymax = `2020_upper`, fill = '2020'), alpha = 0.25)+
  geom_ribbon(aes(temp, ymin = `1942_lower`, ymax = `1942_upper`, fill = '1942'), alpha = 0.25)+
  geom_smooth(data = true_dat %>% filter(lab == 'conditioned'),
              aes(temp, `1942_actual`, col = '1942', linetype = '1942'), se = F)+
  geom_smooth(data = true_dat %>% filter(lab == 'conditioned'),
              aes(temp, `2020_actual`, col = '2020', linetype = '2020'), se = F)+
  scale_x_continuous(limits = c(27, 34.5),
                     breaks = c(28,  30,  32, 34),
                     label = paste0(c(28,  30,  32,34),"°C"))+
  labs(x = "Temperature",
       y = expression(E[o]),
       col = "Year",
       fill = "Year",
       linetype = "Year")+
  theme_minimal(12)+
  theme(axis.text.x = element_text(angle = 0),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = 'none',
        strip.text.x = element_text(size=0)),
bts_dat %>% 
  filter(lab == 'unconditioned') %>%
  ggplot()+
  geom_ribbon(aes(temp, ymin = `2020_lower`, ymax = `2020_upper`, fill = '2020'), alpha = 0.25)+
  geom_ribbon(aes(temp, ymin = `1942_lower`, ymax = `1942_upper`, fill = '1942'), alpha = 0.25)+
  geom_smooth(data = true_dat %>% filter(lab == 'unconditioned'),
              aes(temp, `1942_actual`, col = '1942', linetype = '1942'), se = F)+
  geom_smooth(data = true_dat %>% filter(lab == 'unconditioned'),
              aes(temp, `2020_actual`, col = '2020', linetype = '2020'), se = F)+
  scale_x_continuous(limits = c(27, 34.5),
                     breaks = c(28,  30,  32, 34),
                     label = paste0(c(28,  30,  32,34),"°C"))+
  labs(x = "Temperature",
      y = "",
       col = "Year",
       fill = "Year",
       linetype = "Year")+
  theme_minimal(12)+
  theme(axis.text.x = element_text(angle = 0),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = 'none',
        strip.text.x = element_text(size=0))+
  facet_wrap(~lab, scales = 'free')+
   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                 labels = trans_format("log10", math_format(10^.x))),
widths = c(1, 1), nrow = 1)

ggsave(plot =  plt, paste0("output/figs/prop_exceedance_",marg_mod,".pdf"), height = 3, width =8)
