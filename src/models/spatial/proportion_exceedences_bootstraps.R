rm(list=ls())
library(tidyverse)


calc_prop_exceedence_bts = function(marg_mod, tmp_conditioned_on, bts_range){
  
  setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")
  
  obs_sites = read_csv("data/processed/obs_data.csv") %>%
    dplyr::select(Station, Long.projected, Lat.projected) %>% 
    unique
  
  for(bts in bts_range){
    print(paste0("Fitting bootstrap number: ", bts))
    yr = 2020
    
   # file.exists(paste0("output/simulations/rescaled_simulations_on_obs_grid/bootstraps/bootstrap_", bts, "_model_",marg_mod,"_yr_",1942, "_min_temp_conditioned_on_27"))
      
    if(file.exists(paste0("output/simulations/rescaled_simulations_on_obs_grid/bootstraps/bootstrap_", bts, "_model_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_", tmp_conditioned_on))){
      
      my_simulations  = readRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/bootstraps/bootstrap_", bts, "_model_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_", tmp_conditioned_on))
      prop_ex_2020 = c()
      
      #for(tmp in seq(27, 36, by = 0.25)){
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
      
      #tibble(temp = seq(27, 36, by = 0.25), prop_ex_1942, prop_ex_2020, bts = bts) %>%
      tibble(temp = c(tmp_conditioned_on), prop_ex_1942, prop_ex_2020, bts = bts) %>%
        write_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_bootstrap_",marg_mod,"_conditioned_on_", tmp_conditioned_on, ".csv"), append = T)
    }
  }
  
}


tmp_conditioned_on = 35

job::job({calc_prop_exceedence_bts('mod_1',tmp_conditioned_on, seq(301,325))})
job::job({calc_prop_exceedence_bts('mod_1',tmp_conditioned_on, seq(326,350))})
job::job({calc_prop_exceedence_bts('mod_1',tmp_conditioned_on, seq(351,375))})
job::job({calc_prop_exceedence_bts('mod_1',tmp_conditioned_on, seq(376,400))})
job::job({calc_prop_exceedence_bts('mod_1',tmp_conditioned_on, seq(401,425))})
job::job({calc_prop_exceedence_bts('mod_1',tmp_conditioned_on, seq(426,450))})
job::job({calc_prop_exceedence_bts('mod_1',tmp_conditioned_on, seq(451,475))})
job::job({calc_prop_exceedence_bts('mod_1',tmp_conditioned_on, seq(479,500))})

job::job({calc_prop_exceedence_bts('mod_2',tmp_conditioned_on, seq(301,325))})
job::job({calc_prop_exceedence_bts('mod_2',tmp_conditioned_on, seq(326,350))})
job::job({calc_prop_exceedence_bts('mod_2',tmp_conditioned_on, seq(351,375))})
job::job({calc_prop_exceedence_bts('mod_2',tmp_conditioned_on, seq(376,400))})
job::job({calc_prop_exceedence_bts('mod_2',tmp_conditioned_on, seq(401,425))})
job::job({calc_prop_exceedence_bts('mod_2',tmp_conditioned_on, seq(426,450))})
job::job({calc_prop_exceedence_bts('mod_2',tmp_conditioned_on, seq(451,475))})
job::job({calc_prop_exceedence_bts('mod_2',tmp_conditioned_on, seq(476,500))})









job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,  seq(1,25))})
job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,  seq(26,50))})
job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,  seq(51,75))})
job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,  seq(76,100))})

job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,   seq(101,125))})
job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,  seq(126,150))})
job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,   seq(151,175))})
job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,   seq(176,200))})

job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,   seq(201,225))})
job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,   seq(226,250))})
job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,   seq(251,275))})
job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,  seq(276,300))})

job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,  seq(301,325))})
job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,   seq(326,350))})
job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,   seq(351,375))})
job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,  seq(376,400))})

job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,  seq(401,425))})
job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,  seq(426,450))})
job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,  seq(451,475))})
job::job({calc_prop_exceedence_bts('mod_4', tmp_conditioned_on,   seq(476,500))})


# job::job({calc_prop_exceedence_bts('mod_4', seq(1,25))})
# job::job({calc_prop_exceedence_bts('mod_4', seq(26,50))})
# job::job({calc_prop_exceedence_bts('mod_4', seq(51,75))})
# job::job({calc_prop_exceedence_bts('mod_4', seq(76,100))})
# 
# job::job({calc_prop_exceedence_bts('mod_4', seq(101,125))})
# job::job({calc_prop_exceedence_bts('mod_4', seq(126,150))})
# job::job({calc_prop_exceedence_bts('mod_4', seq(151,175))})
# job::job({calc_prop_exceedence_bts('mod_4', seq(176,200))})
# 
# job::job({calc_prop_exceedence_bts('mod_4', seq(201,225))})
# job::job({calc_prop_exceedence_bts('mod_4', seq(226,250))})
# job::job({calc_prop_exceedence_bts('mod_4', seq(251,275))})
# job::job({calc_prop_exceedence_bts('mod_4', seq(276,300))})
# 
# job::job({calc_prop_exceedence_bts('mod_4', seq(301,325))})
# job::job({calc_prop_exceedence_bts('mod_4', seq(326,350))})
# job::job({calc_prop_exceedence_bts('mod_4', seq(351,375))})
# job::job({calc_prop_exceedence_bts('mod_4', seq(376,400))})
# 
# job::job({calc_prop_exceedence_bts('mod_4', seq(401,425))})
# job::job({calc_prop_exceedence_bts('mod_4', seq(426,450))})
# job::job({calc_prop_exceedence_bts('mod_4', seq(451,475))})
# job::job({calc_prop_exceedence_bts('mod_4', seq(476,500))})

# job::job({calc_prop_exceedence_bts('mod_1', seq(301,325))})
# job::job({calc_prop_exceedence_bts('mod_1', seq(326,350))})
# job::job({calc_prop_exceedence_bts('mod_1', seq(351,375))})
# job::job({calc_prop_exceedence_bts('mod_1', seq(376,400))})
# 
# job::job({calc_prop_exceedence_bts('mod_1', seq(401,425))})
# job::job({calc_prop_exceedence_bts('mod_1', seq(426,450))})
# job::job({calc_prop_exceedence_bts('mod_1', seq(451,475))})
#job::job({calc_prop_exceedence_bts('mod_1', seq(479,500))})
 
# job::job({calc_prop_exceedence_bts('mod_2', seq(301,325))})
# job::job({calc_prop_exceedence_bts('mod_2', seq(326,350))})
# job::job({calc_prop_exceedence_bts('mod_2', seq(351,375))})
# job::job({calc_prop_exceedence_bts('mod_2', seq(376,400))})
# 
# job::job({calc_prop_exceedence_bts('mod_2', seq(401,425))})
# job::job({calc_prop_exceedence_bts('mod_2', seq(426,450))})
# job::job({calc_prop_exceedence_bts('mod_2', seq(451,475))})
# job::job({calc_prop_exceedence_bts('mod_2', seq(476,500))})
#

marg_mod = 'mod_4'
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





prob_observing_T = read_csv(paste0("output/prob_observing_T_anywhere_clim_",marg_mod,".csv")) %>%
  rename(temp = temps) %>%
  pivot_longer(-temp) %>%
  mutate(year = str_remove(name,"prob_")) %>%
  mutate(year = str_remove(year,"_u")) %>%
  mutate(year = str_remove(year,"_l"))


prob_observing_T = rbind(prob_observing_T %>% filter(name == 'prob_1942') %>% mutate(typ = 'actual') %>% dplyr::select(-name),
      prob_observing_T %>% filter(name == 'prob_2020') %>% mutate(typ = 'actual') %>% dplyr::select(-name),
      prob_observing_T %>% filter(name == 'prob_1942_u') %>% mutate(typ = 'upper') %>% dplyr::select(-name),
      prob_observing_T %>% filter(name == 'prob_2020_u') %>% mutate(typ = 'upper') %>% dplyr::select(-name),
      prob_observing_T %>% filter(name == 'prob_1942_l') %>% mutate(typ = 'lower') %>% dplyr::select(-name),
      prob_observing_T %>% filter(name == 'prob_2020_l') %>% mutate(typ = 'lower') %>% dplyr::select(-name)) %>%
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
  left_join(prob_observing_T)

true_dat = true_dat %>%
  pivot_longer(-temp) %>%
  mutate(year = str_remove(name,"prop_ex_")) 
true_dat$typ = 'actual'
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



plt = bts_dat %>%
  ggplot()+
  geom_ribbon(aes(temp, ymin = `2020_lower`, ymax = `2020_upper`, fill = '2020'), alpha = 0.25)+
  geom_ribbon(aes(temp, ymin = `1942_lower`, ymax = `1942_upper`, fill = '1942'), alpha = 0.25)+
  geom_smooth(data = true_dat,
              aes(temp, `1942_actual`, col = '1942', linetype = '1942'), se = F)+
  geom_smooth(data = true_dat,
              aes(temp, `2020_actual`, col = '2020', linetype = '2020'), se = F)+
  scale_x_continuous(limits = c(27, 34.5),
                     breaks = c(28,  30,  32, 34),
                     label = paste0(c(28,  30,  32,34),"°C"))+
  labs(x = "Temperature",
       y = expression(E[o]^D),
       col = "Year",
       fill = "Year",
       linetype = "Year")+
  theme_minimal(12)+
  theme(axis.text.x = element_text(angle = 0),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = 'none',
        strip.text.x = element_text(size=0))+
  facet_wrap(~lab, scales = 'free')


library(scales)


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
       y = expression(E[o]^`D|T`),
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
       y = expression(E[o]^D),
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
                labels = trans_format("log10", math_format(10^.x))), nrow = 1)

ggsave(plot =  plt, paste0("output/figs/prop_exceedance_",marg_mod,".pdf"), height = 3, width =8)




true_dat %>% filter(lab == 'unconditioned') %>%
  filter(temp == 28)





0.00407/0.00140



true_dat %>% filter(lab == 'unconditioned') %>%
  mutate(dif = `2020_actual`/`1942_actual`)



# bts_dat %>% 
#   filter(lab == 'unconditioned') %>%
#   mutate(`2020_upper` = 1/(92*`2020_upper`) ,
#          `2020_lower` = 1/(92*`2020_lower`),
#          `1942_upper` = 1/(92*`1942_upper`),
#          `1942_lower` = 1/(92*`1942_lower`))
#          
plt = bts_dat %>% 
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
       y = expression(E[o]^D),
       col = "Year",
       fill = "Year",
       linetype = "Year")+
  theme_minimal(12)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  theme(
    axis.text.x = element_text(angle = 0),
    axis.title.y = element_text(angle = 0, vjust = 0.5),
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

ggsave(plot = plt, filename = "prop_exceedence.png", bg = 'transparent', width = 4.5, height = 3.5)
#   
# bts_dat %>% 
#   filter(lab == 'unconditioned') %>%
#   ggplot()+
#   geom_segment(aes(x = (temp), xend = (temp), y = 1/(92*`2020_lower`), yend = 1/(92*`2020_upper`), col = '2020'))+
#   geom_segment(aes(x = temp, xend = temp, y = 1/(92*`1942_lower`), yend = 1/(92*`1942_upper`), col = '1942'))+
#   geom_point(data = true_dat %>% filter(lab == 'unconditioned'), aes(x = temp, y = 1/(92*`1942_actual`), col = '1942'))+
#   geom_point(data = true_dat %>% filter(lab == 'unconditioned'), aes(x = (temp), y = 1/(92*`2020_actual`), col = '2020'))+
#   theme_minimal()+
#   labs(x = "Temperature",
#        y = "Return period (Years)",
#        fill = 'Year', col = "Year")+
#   theme_minimal(12)+
#   theme(legend.position = 'none')+
#   scale_x_continuous(limits = c(27, 34),
#                      breaks = c(25, 28,  30,  32, 34),
#                      label = paste0(c(26, 28,  30,  32,34),"°C"))+
#   scale_y_log10(limits = c(1, 23000),breaks = c(1, 10,  100, 1000, 10000),labels = c(1, 10,  100, 1000, 10000))+ coord_flip()
# 
# 
# 
#   
# 
# 637084
# 
# # 
# # ggsave("output/figs/prop_exceedance_mod_4.pdf", height = 3, width =6)
# 
# 
# 
# 
# bts_dat %>%
#   ggplot()+
#   geom_ribbon(aes(temp, ymin = prop_ex_2020_l, ymax = prop_ex_2020_u, fill = '2020'), alpha = 0.25)+
#   geom_ribbon(aes(temp, ymin = prop_ex_1942_l, ymax = prop_ex_1942_u, fill = '1942'), alpha = 0.25)+
#   geom_smooth(data =true_data,
#             aes(temp, prop_ex_1942, col = '1942', linetype = '1942'), se = F)+
#   geom_smooth(data = true_data,
#             aes(temp, prop_ex_2020, col = '2020', linetype = '2020'), se = F)+
#   scale_x_continuous(limits = c(27, 34.5),
#                      breaks = c(28,  30,  32, 34),
#                      label = paste0(c(28,  30,  32,34),"°C"))+
#   labs(x = "Temperature",
#        y = expression(E[o]^D),
#        col = "Year",
#        fill = "Year",
#        linetype = "Year")+
# theme_minimal(12)+
#   theme(axis.text.x = element_text(angle = 0),
#         axis.title.y = element_text(angle = 0, vjust = 0.5),
#         legend.position = 'none',
#         strip.text.x = element_text(size=0))+
#   facet_wrap(~lab, scales = 'free')

ggsave(paste0("output/figs/prop_exceedance_",marg_mod,".pdf"), height = 3, width =6)













# ------ PLOT MODEL 0


bts_dat = read_csv("output/simulations/simulation_summary/prop_exceedance_model_bootstrap_mod_2_.csv",
                   col_names = c('temp', 'prop_ex_1942', 'prop_ex_2020', 'bts'))
bts_dat$bts %>% unique

bts_dat$prop_ex_1942[is.na(bts_dat$prop_ex_1942)] = 0
bts_dat$prop_ex_2020[is.na(bts_dat$prop_ex_2020)] = 0


bts_dat = bts_dat %>%
  group_by(temp) %>%
  summarise(prop_ex_2020_u = quantile(prop_ex_2020, 0.9),
            prop_ex_2020_l = quantile(prop_ex_2020, 0.1),
            prop_ex_1942_u = quantile(prop_ex_1942, 0.9),
            prop_ex_1942_l = quantile(prop_ex_1942, 0.1))

# group by temps and plot line?

true_data = read_csv("output/simulations/simulation_summary/prop_exceedance_model_true_mod_2.csv",
                     col_names = c('temp', 'prop_ex_1942', 'prop_ex_2020'))
true_data$prop_ex_1942[is.na(true_data$prop_ex_1942)] = 0
true_data$prop_ex_2020[is.na(true_data$prop_ex_2020)] = 0




# ----- Uncondition
prob_1942 = read_csv("output/prob_observing_T_anywhere_clim_mod_2.csv") %>%
  filter(temps == 27) %>%
  pull(prob_1942)

prob_2020 = read_csv("output/prob_observing_T_anywhere_clim_mod_2.csv") %>%
  filter(temps == 27) %>%
  pull(prob_2020)


bts_unconditioned = bts_dat
true_unconditioned = true_data
bts_unconditioned$prop_ex_2020_u = prob_2020*bts_unconditioned$prop_ex_2020_u
bts_unconditioned$prop_ex_2020_l = prob_2020*bts_unconditioned$prop_ex_2020_l
bts_unconditioned$prop_ex_1942_u = prob_1942*bts_unconditioned$prop_ex_1942_u
bts_unconditioned$prop_ex_1942_l = prob_1942*bts_unconditioned$prop_ex_1942_l
true_unconditioned$prop_ex_2020 = prob_2020*true_unconditioned$prop_ex_2020
true_unconditioned$prop_ex_1942 = prob_1942*true_unconditioned$prop_ex_1942



bts_dat = rbind(bts_dat %>%
  mutate(lab = "conditioned"),
bts_unconditioned %>%
  mutate(lab = "unconditioned")) %>%
  filter(temp < 33)

true_data = rbind(true_data %>%
        mutate(lab = "conditioned"),
      true_unconditioned %>%
        mutate(lab = "unconditioned"))%>%
  filter(temp < 33)





bts_dat %>%
  ggplot()+
  geom_ribbon(aes(temp, ymin = prop_ex_2020_l, ymax = prop_ex_2020_u, fill = '2020'), alpha = 0.25)+
  geom_ribbon(aes(temp, ymin = prop_ex_1942_l, ymax = prop_ex_1942_u, fill = '1942'), alpha = 0.25)+
  geom_line(data =true_data,
            aes(temp, prop_ex_1942, col = '1942', linetype = '1942'), size  =1)+
  geom_line(data = true_data,
            aes(temp, prop_ex_2020, col = '2020', linetype = '2020'), size  =1)+
  scale_x_continuous(limits = c(27, 33),
                     breaks = c(28,  30,  32, 34),
                     label = paste0(c(28,  30,  32,34),"°C"))+
  labs(x = "Temperature",
       y = expression(E[o]^D),
       col = "Year",
       fill = "Year",
       linetype = "Year")+
theme_minimal(12)+
  theme(axis.text.x = element_text(angle = 0),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = 'none',
        strip.text.x = element_text(size=0))+
  facet_wrap(~lab, scales = 'free')

ggsave("output/figs/prop_exceedance_mod_2.pdf", height = 3, width =6)

