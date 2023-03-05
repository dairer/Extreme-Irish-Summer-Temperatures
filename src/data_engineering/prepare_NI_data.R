rm(list=ls())
library(tidyverse)
library(vroom)
setwd("~/Extreme-Irish-Summer-Temperatures/data/NI/")

NI_data = list.files() %>%
  paste0('~/Extreme-Irish-Summer-Temperatures/data/NI/',.) %>%
  purrr::map(~{
    
    site_info = str_extract(.x, "202107.*qcv") %>%
      str_remove("202107_") %>%
      str_remove("_qcv") %>%
      str_split("_")
    
    line_to_start = grep("^data$", readLines(.x)) # csv data starts after "data"
    
    my_file_read = read.csv(.x, skip = line_to_start) %>%
      .[-nrow(.),] # remove last line
    
    
    my_file_read %>%
      dplyr::select(ob_end_time,
                    src_id,
                    max_air_temp, 
                    max_air_temp_j,
                    max_air_temp_q) %>%
      mutate(Station = site_info[[1]][3])
  }) %>%
  plyr::rbind.fill()%>%
  as_tibble() %>%
  left_join(read.csv("../NI_meta.csv"))


# create indicator vars
inds = NI_data$max_air_temp_q %>% as.character()

NI_data$ind_1 = inds %>% substr(start = nchar(inds), stop = nchar(inds)) %>%
  as.numeric()

NI_data$ind_2 = inds %>% substr(start = nchar(inds)-1, stop = nchar(inds)-1) %>%
  as.numeric()

# remove if ind > 0, its either suspect or estimated (https://dap.ceda.ac.uk/badc/ukmo-midas/metadata/doc/QC_J_flags.html)
NI_data$ind_3 = inds %>% substr(start = nchar(inds)-2, stop = nchar(inds)-2) %>%
  as.numeric()

NI_data$ind_4 = inds %>% substr(start = nchar(inds)-3, stop = nchar(inds)-3) %>%
  as.numeric()

NI_data = NI_data %>%
  filter(is.na(ind_3) | (ind_3 == 0)) %>%
  dplyr::select(ob_end_time,
                src_id,
                max_air_temp,
                Station,
                Lat,
                Long)

# where there is more than one obs in a day, take the maximum
NI_data = NI_data %>% 
  group_by(Station) %>%
  group_map(~{
    .x %>%
      mutate(date = lubridate::as_date(ob_end_time)) %>%
      group_by(date) %>%
      slice(which.max(max_air_temp)) %>%
      ungroup()
  },.keep = T)%>%
  # save as one big tibble
  plyr::rbind.fill()%>%
  as_tibble() %>%
  dplyr::select(date, Station, Long, Lat, id = src_id, maxtp = max_air_temp)

# save as csv
NI_data %>%
  write_csv("~/Extreme-Irish-Summer-Temperatures/data/processed/NI_data.csv")
