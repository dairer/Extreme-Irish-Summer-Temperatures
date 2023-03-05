gc()
rm(list = ls())
library(tidyverse)

setwd("~/Extreme-Irish-Summer-Temperatures/")
source('src/models/marginal_models/gpd_models.R')


obs_sites = read_csv("data/processed/obs_data.csv") %>%
  dplyr::select(Station, scale_9, threshold_9, Long.projected, Lat.projected) %>%
  unique() %>%
  left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>%
              dplyr::select(Station, dist_sea))

obs_grid = read_csv(paste0("data/processed/thresh_exceedance_lambda_num_quantiles_25.csv"))  %>%
  left_join(obs_sites) %>%
  left_join(read_csv("data/processed/obs_data.csv") %>%
              dplyr::select(year, loess_temp_anom) %>% unique)%>%
  filter(year %in% c(1942, 2020))

gpd_pars = read_csv("output/gpd_model_fits/model_2_true.csv")
obs_grid$scale = my_predict_2(unlist(gpd_pars), obs_grid$scale_9, obs_grid$loess_temp_anom)$scale
obs_grid$shape = my_predict_2(unlist(gpd_pars), obs_grid$scale_9, obs_grid$loess_temp_anom)$shape

# ---- Read in simulations
my_simulations_extremes = c()
for(i in seq(1, 100)){
my_simulations_extremes = c(my_simulations_extremes, readRDS(paste0("output/simulations/simulations_on_obs_grid/true/mod_2_run_",i)))
}

my_simulations_extremes = my_simulations_extremes[seq(25000)]

# ---- standardise simulations to have cost = 1
my_simulations_standardised = list()
for (i in 1:length(my_simulations_extremes)) {
  this_cost = mean(my_simulations_extremes[[i]])
  my_simulations_standardised[[i]] = my_simulations_extremes[[i]]/this_cost
}

# --- maximum frechet margin at each of my simulation sites
max_at_each_site = c()
for(s in seq(length(my_simulations_standardised[[1]]))){
  max_at_each_site = c(max_at_each_site, lapply(my_simulations_standardised, "[[", s) %>% unlist %>% max)
}

m = length(my_simulations_standardised)
L = 300

res_1942 = c()
res_2020 = c()

for(temp_i_want in seq(25,36)){
  print(temp_i_want)

  obs_grid$frechet = -1/log(1-obs_grid$thresh_exceedance_9*(1-evd::pgpd(temp_i_want - obs_grid$threshold_9,  scale = obs_grid$scale, shape =  obs_grid$shape[1] )) )
  T_2020 = obs_grid %>% filter(year == 2020) %>% pull(frechet) 
  T_1942 = obs_grid %>% filter(year == 1942) %>% pull(frechet) 
  
  T_1942[T_1942 == -Inf] = Inf
  T_2020[T_2020 == -Inf] = Inf
  
  b_2020 = min(T_2020/max_at_each_site)
  b_1942 = min(T_1942/max_at_each_site)
  m = length(my_simulations_standardised)
  
  above_1942 = 0
  above_2020 = 0
  for(i in seq(L)){
    print(i)
    this_set_scaled = my_simulations_standardised %>%
      map(~{ 
        .x*evd::rgpd(n=1, loc = 1, scale = 1, shape = 1)
      }) %>%
      map(~{ 
        if(mean(.x)>=r_thresh){
          c(sum(.x*b_2020 > T_2020) > 0,
            sum(.x*b_1942 > T_1942) > 0)
        }
      })
    above_2020 = above_2020 + (lapply(this_set_scaled, "[[", 1) %>% unlist %>% sum)
    above_1942 = above_1942 + (lapply(this_set_scaled, "[[", 2) %>% unlist %>% sum)
  }
  
  tibble(temp = temp_i_want, 
         p_1942 = above_1942 / (length(my_simulations_standardised)*L*b_1942),
         p_2020 = above_2020 / (length(my_simulations_standardised)*L*b_2020)) %>%
  write_csv("output/prob_extreme_temp_imp_samp_mod_2.csv", append = T)
}


read_csv("output/prob_extreme_temp_imp_samp_mod_2.csv",
         col_names = c('bts','temp', 'p_1942', 'p_2020')) %>%
  filter(temp<=35, temp>25) %>%
  group_by(temp) %>%
  summarise(p_1942_lower = quantile(p_1942, 0.1, na.rm=T),
            p_1942_upper = quantile(p_1942, 0.9, na.rm=T),
            p_2020_lower = quantile(p_2020, 0.1, na.rm=T),
            p_2020_upper = quantile(p_2020, 0.9, na.rm=T)) %>%
  ggplot()+
  geom_ribbon(aes(temp,  ymin = 1/(92*p_1942_lower), ymax = 1/(92*p_1942_upper), fill = '1942'), alpha = 0.5)+
  geom_ribbon(aes(temp,  ymin = 1/(92*p_2020_lower), ymax = 1/(92*p_2020_upper), fill = '2020'), alpha = 0.5)+
  geom_line(data = read_csv("output/prob_extreme_temp_imp_samp_mod_2.csv",
                            col_names = c('temp', 'p_1942', 'p_2020')) %>%
              filter(temp<=35, temp>25), aes(temp, 1/(92*p_1942), col = '1942'))+
  geom_line(data = read_csv("output/prob_extreme_temp_imp_samp_mod_2.csv",
                            col_names = c('temp', 'p_1942', 'p_2020')) %>%
              filter(temp<=35, temp>25), aes(temp, 1/(92*p_2020), col = '2020'), linetype = 'dashed')+
  theme_minimal()+
  labs(x = "Temperature",
       y = "Return period (Years)",
       fill = 'Year', col = "Year")+
  theme_minimal(12)+
  scale_x_continuous(limits = c(26, 34),
                     breaks = c(26, 28,  30,  32, 34),
                     label = paste0(c(26, 28,  30,  32,34),"°C"))+
  scale_y_log10(breaks = c(0.1, 1, 10,  100,  1000),
                labels = c(0.1, 1, 10,  100,  1000))+
  coord_flip(ylim = c(0.075, 1000), xlim = c(26, 34))+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = 'none')
ggsave(filename = "output/figs/spatial_rl_mod_2.pdf",  width = 5, height = 3)
