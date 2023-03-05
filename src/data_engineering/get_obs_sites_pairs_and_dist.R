# --- data used for validating r-paretp process

locs_to_pred = read_csv("data/processed/obs_data_dist_to_sea.csv")
site_pairs = as.data.frame(t(combn(locs_to_pred$Station , 2 )))

dist_h <- function(long1, lat1, long2, lat2) {
  sqrt((long1 - long2)^2 + (lat1 - lat2)^2)
}

dists = c()
for(s in seq(nrow(site_pairs))){
  print(s)
  site_1 = locs_to_pred[locs_to_pred$Station == site_pairs[s,]$V1,]
  site_2 = locs_to_pred[locs_to_pred$Station == site_pairs[s,]$V2,]
  
  dists = c(dists, dist_h(site_1$Long.projected[1],
                          site_1$Lat.projected[1],
                          site_2$Long.projected[1],
                          site_2$Lat.projected[1]))
}

site_pairs$dist = dists

site_pairs %>%
  write_csv("data/processed/obs_pairs_with_dist.csv")
