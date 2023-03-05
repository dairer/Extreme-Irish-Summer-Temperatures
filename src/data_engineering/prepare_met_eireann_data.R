# This script reads in raw Met Éireann data and saves max daily temperatures to a csv

rm(list=ls())
library(tidyverse)
setwd("~/Extreme-Irish-Summer-Temperatures/data/met-eireann-data/")

all_data = list.files() %>%
  paste0(.,'/',.,'.csv') %>%
  purrr::map(~{
    
    print(.x)
    x <- readLines(.x)
    nskip = x %>% grep("date",.) %>% .[2]
    
    # read in location data
    my_file_read <- read_csv(.x, skip = nskip - 1)
    
    position_maxt = match("maxt",names(my_file_read))
    
    if(is.na(position_maxt)){
      position_maxt = match("maxtp",names(my_file_read))
    }
    
    if(!is.na(position_maxt)){
      # get code string
      this_id = str_extract(.x, "dly[0-9]{3,}.csv") %>%
        str_remove('.csv')%>%
        str_remove('dly')
      
      tibble(id =this_id,
             date = lubridate::dmy(my_file_read$date), 
             maxtp = as.numeric(unlist(my_file_read[,position_maxt])),
             ind = as.numeric(unlist(my_file_read[,position_maxt-1])))  %>%# indicator string always proceeds variable
        filter(!is.na(maxtp)) %>%
        mutate(ind = as.numeric(ind)) %>%
        filter(ind <= 1) # remove all indicatior variables that are not 0 or 1
    }
  }, .keep = T)%>%
  # save as one big tibble
  plyr::rbind.fill() %>%
  as_tibble()


Stations_meta = read_csv("~/Extreme-Irish-Summer-Temperatures/data/met_eireann_meta.csv") %>%
  rename(loc_id = `Station Number`,
         Station = name,
         Lat = Latitude,
         Long = Longitude) %>%
  dplyr::select(loc_id, Station, Lat, Long)


all_data = all_data %>%
  mutate(loc_id = as.numeric(id)) %>%
  left_join(Stations_meta)


all_data %>% write_csv("~/Extreme-Irish-Summer-Temperatures/data/processed/met_eireann_data.csv")
