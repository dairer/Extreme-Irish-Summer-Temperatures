rm(list = ls())
library(tidyverse)
setwd("~/Extreme-Irish-Summer-Temperatures/data/processed/")

met_eireann = read.csv("met_eireann_data.csv") %>% 
  as_tibble %>%
  dplyr::select(loc_id,date,maxtp,Station,Lat,Long)

met_eireann[met_eireann$Station == "DUBLIN (MERRION SQUARE)",]$loc_id = 3923 # --- we have 2 "DUBLIN (MERRION SQUARE)" sites in me_sites

ni = read.csv("NI_data.csv") %>% 
  as_tibble %>%
  dplyr::select(loc_id = id,date,maxtp,Station,Lat,Long)

# ----- combine data
obs_data = rbind(met_eireann,cs3,ni)
obs_data = obs_data %>% unique # --- remove mysterious duplicate observations

# --- create tibble describing location and data span of each site
sites = obs_data %>%
  group_by(Station) %>%
  summarise(Long, Lat,loc_id, num_obs = n(),
            min_date = min(date),
            max_date = max(date)) %>%
  unique() %>%
  ungroup()


# --- removing sites with very little data
sites = sites %>%
  filter(num_obs > (365*5))

obs_data = obs_data %>%
  filter(Station %in% sites$Station)

# ------ Project to UTM
my_coords = sites %>% dplyr::select(Long, Lat)
sp::coordinates(my_coords)<- c("Long", "Lat")
sp::proj4string(my_coords) <- sp::CRS("+proj=longlat +datum=WGS84") # add Coordinate Reference System, which is WGS84
proj_cords <- sp::spTransform(my_coords, sp::CRS(paste0("+proj=utm +zone=29 ellps=WGS84"))) # project to utm zone 29
proj_cords = proj_cords %>% sp::coordinates()/100000
sites$Long.projected = proj_cords[,1]
sites$Lat.projected = proj_cords[,2]

dist_h <- function(long1, lat1, long2, lat2) sqrt((long1 - long2)^2 + (lat1 - lat2)^2)

cluster_sites = function(my_sites, distance_epsilon){
  clustered_sites = c()
  for(i in seq(nrow(my_sites) - 1)){
    stations_close = c()
    for(j in seq((i+1),nrow(my_sites))){
      this_dist = dist_h(my_sites[i,]$Long.projected,
                         my_sites[i,]$Lat.projected,
                         my_sites[j,]$Long.projected,
                         my_sites[j,]$Lat.projected)
      
      if(this_dist < distance_epsilon){
        
        stations_close = c(stations_close,my_sites$Station[j])
      }
    }
    clustered_sites = rbind(clustered_sites,
                            tibble(Station = my_sites[i,]$Station, 
                                   stations_close =list(stations_close)))
  }
  clustered_sites %>% drop_na()
}
sites_to_cluster = cluster_sites(sites, 0.025)


# ----- combine sites that are close
for(i in seq(nrow(sites_to_cluster))){
  
  print(paste0("Combining ",sites_to_cluster[i,]$Station," and ",sites_to_cluster[i,]$stations_close[[1]]))

  s1 = sites_to_cluster[i,]$Station
  s2 = sites_to_cluster[i,]$stations_close[[1]]
  
  this_pair = sites %>%
    filter(Station %in% c(s1, s2))
  
  # ---- one site is a sub set of another
  id_former = this_pair$min_date %>% as.Date %>% which.min()
  id_latter = this_pair$max_date %>% as.Date %>% which.max()
  
  if(id_former == id_latter){
    # remove other site (the site that is a subset)
    obs_data = obs_data %>%
      filter(Station != this_pair[c(1,2)[c(1,2) != id_former],]$Station)
    
  }else if((this_pair$max_date[1] < this_pair$min_date[2]) || (this_pair$max_date[2] < this_pair$min_date[1])){
    # ---- sites don't overlap
    # --- check if one max is greater than the other min
    if(this_pair$max_date[1] < this_pair$min_date[2]){
      print(paste0("Rename ", this_pair$Station[1]," to ", this_pair$Station[2]))
      obs_data$Long[obs_data$Station == this_pair$Station[1]] = obs_data %>% filter(Station == this_pair$Station[2]) %>% pull(Long) %>% .[1]
      obs_data$Lat[obs_data$Station == this_pair$Station[1]] = obs_data %>% filter(Station == this_pair$Station[2]) %>% pull(Lat) %>% .[1]
      obs_data$loc_id[obs_data$Station == this_pair$Station[1]] = obs_data %>% filter(Station == this_pair$Station[2]) %>% pull(loc_id) %>% .[1]
      obs_data$Station[obs_data$Station == this_pair$Station[1]] = this_pair$Station[2]
      sites = sites %>% filter(Station != this_pair$Station[1])
    }else{
      print(paste0("Rename ", this_pair$Station[2]," to ", this_pair$Station[1]))
      obs_data$Long[obs_data$Station == this_pair$Station[2]] = obs_data %>% filter(Station == this_pair$Station[1]) %>% pull(Long) %>% .[1]
      obs_data$Lat[obs_data$Station == this_pair$Station[2]] = obs_data %>% filter(Station == this_pair$Station[1]) %>% pull(Lat) %>% .[1]
      obs_data$loc_id[obs_data$Station == this_pair$Station[2]] = obs_data %>% filter(Station == this_pair$Station[1]) %>% pull(loc_id) %>% .[1]
      obs_data$Station[obs_data$Station == this_pair$Station[2]] = this_pair$Station[1]
      sites = sites %>% filter(Station != this_pair$Station[2])
    }
  }else{
    # ---- sites overlap
    # find which site is older
    # which site is more moders
    keep = this_pair[id_latter,]
    remv = this_pair[-id_latter,]
    
    obs_data = obs_data %>% filter(!(Station == remv$Station & (lubridate::as_date(date)>=lubridate::as_date(keep$min_date)))) 
    obs_data$Long[obs_data$Station == remv$Station] = keep$Long
    obs_data$Lat[obs_data$Station == remv$Station] = keep$Lat
    obs_data$loc_id[obs_data$Station == remv$Station] = keep$loc_id
    obs_data$Station[obs_data$Station == remv$Station] = keep$Station
    sites = sites %>% filter(Station != remv$Station)
  }
}

obs_data$Station = str_to_lower(obs_data$Station)

sites = obs_data %>%
  group_by(Station) %>%
  summarise(Long, Lat, loc_id, num_obs = n(),
            min_date = min(date),max_date = max(date)) %>%
  unique() %>%
  ungroup()

my_coords = obs_data %>% dplyr::select(Long, Lat)
sp::coordinates(my_coords)<- c("Long", "Lat")
sp::proj4string(my_coords) <- sp::CRS("+proj=longlat +datum=WGS84") # add Coordinate Reference System, which is WGS84
proj_cords <- sp::spTransform(my_coords, sp::CRS(paste0("+proj=utm +zone=29 ellps=WGS84"))) # project to utm zone 29
proj_cords = proj_cords %>% sp::coordinates()/100000
obs_data$Long.projected = proj_cords[,1]
obs_data$Lat.projected = proj_cords[,2]

obs_data = obs_data %>%
  mutate(date = lubridate::as_date(date))%>%
  filter(
    !(Station == "armagh" &  lubridate::year(date) <=1930), # dont model such old data
    !(Station == "littleton ii b. na m." &  lubridate::year(date) <=2000)) %>% # data too pathcy/suspicious before this date
  filter(!Station %in% c(
    "recess (cloonacartan)", # ---- too few data
    "tuam (airglooney)", # ---- too, suspicious data
    "ballynahinch (for.stn.)" # data too pathcy/suspicious before this date
  )) %>%
  filter(maxtp  < 35)# --- not a believable record

obs_data %>%
  as_tibble() %>%
  unique() %>%
  write_csv("obs_all_data.csv")
