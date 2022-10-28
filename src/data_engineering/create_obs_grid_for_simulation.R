rm(list=ls())
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")



# remove sites that are too close to simualte
obs_sites_1 = read_csv("data/processed/obs_data_dist_to_sea.csv") %>% 
  dplyr::select(Station, Long.projected, Lat.projected)

obs_sites_2 = read_csv("data/processed/obs_data_dist_to_sea.csv") %>% 
  dplyr::select(Station, Long.projected, Lat.projected)

dist_h <- function(long1, lat1, long2, lat2) {
  sqrt((long1 - long2)^2 + (lat1 - lat2)^2)
}

# itterate through all stations and find climate grid point
closest_irel = c()
dist_smallest = c()
for(i in seq(nrow(obs_sites_1))){
  smallest_dist = 9999999999
  id_of_smallest_dist_irel = 9999999999
  
  this_dist = 10000
  for(j in seq(nrow(obs_sites_2))){
    this_dist = dist_h(obs_sites_1[i,]$Long.projected,
                       obs_sites_1[i,]$Lat.projected,
                       obs_sites_2[j,]$Long.projected,
                       obs_sites_2[j,]$Lat.projected)
    
    if((this_dist < smallest_dist) & this_dist > 0){
      smallest_dist = this_dist
      id_of_smallest_dist_irel = obs_sites_2$Station[j]
    }
  }
  dist_smallest = c(dist_smallest,smallest_dist)
  
  closest_irel = c(closest_irel, id_of_smallest_dist_irel)
}



sites_to_simulate = tibble(s1 = obs_sites_1$Station, s2 = closest_irel, dist_smallest) %>%
  arrange(desc(dist_smallest)) %>%
  head(86) 

sites_to_keep = c(sites_to_simulate$s1, sites_to_simulate$s2) %>% unique



sites_to_sim = read_csv("data/processed/obs_data_dist_to_sea.csv") %>% 
  filter(Station %in% sites_to_keep)

locs_to_pred = sites_to_sim %>%
  dplyr::select(Long.projected, Lat.projected)%>% 
  as.data.frame()



list(sites_to_sim = sites_to_sim, locs_to_pred = locs_to_pred) %>% saveRDS("data/processed/obs_grid_to_simulate")
