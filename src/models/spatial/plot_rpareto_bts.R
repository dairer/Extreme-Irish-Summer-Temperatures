
# return level plots
rm(list = ls())
library(tidyverse)
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")


variogram_model_score <- function(h, vario_params){
  alpha <-  (vario_params[1])
  beta <- (vario_params[2])
  nu = 0.2
  res = rep(NA, length(h))
  res[h == 0] = 0
  res[h != 0] = alpha*(1 - ((((h[h != 0] /beta)^nu) * besselK(x = (h[h != 0] / beta), nu = nu))/(gamma(nu)*2^(nu - 1))))
  res
}



emp_chi = read_csv(paste0("~/JRSS_organised_code/corrections/data/processed_data/empirical_chi/chi_obs_bts_25_bins0.95.csv"),
                   col_names = c("i", "bts", "dist", "chi")) %>%
  group_by(dist) %>% 
  summarise(upper = quantile(chi, 0.025, na.rm = T),
            lower = quantile(chi, 0.975, na.rm = T)) 





# ------- MODEL 3
res = read_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_mod_4.csv-Mead.csv"),
               col_names = c("bts", "alpha", "beta")) 

true_vals = read_csv(paste0("output/rpareto_model_fits/true_rpareto_fits_model_mod_4.csv"),
                     col_names = F) %>% 
  .[nrow(.),] %>% # get last row in this file
  unlist %>% as.numeric()





res$alpha %>% hist
res$beta %>% hist
distances = seq(0.005, 4, length.out = 100)

my_mat = c()
for(i in seq(nrow(res))){
  my_mat = rbind(my_mat, 2 * (1 - pnorm(sqrt(variogram_model_score(distances, c(res[i,]$alpha, res[i,]$beta)) / 2))))
}

tibble(distances,
       true =  2 * (1 - pnorm(sqrt(variogram_model_score(distances, true_vals) / 2))),
       lower = my_mat %>% apply(MARGIN = 2, FUN = quantile, 0.025, na.rm = T),
       upper = my_mat %>% apply(MARGIN = 2, FUN = quantile, 0.975, na.rm = T)) %>%
  ggplot()+
  geom_line(aes(distances*100, true))+
  geom_segment(data = emp_chi, aes(x = dist*100, xend = dist*100, y = lower, yend = upper), alpha = 0.5)+
  geom_ribbon(aes(distances*100, ymin = lower, ymax = upper), alpha = 0.25)+
  theme_minimal(12)+
  theme(legend.position = 'none')+
  labs(x = "Distance (KM)",
       y = expression(chi[o]^F))+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

ggsave("output/figs/rparto_model_4_fit.pdf", height = 2.5, width = 4)




# -------- MODEL 1
res = read_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_mod_1_L-BFGS-B.csv"),
               col_names = c("bts", "alpha", "beta")) %>%
  filter(bts>300)

true_vals = read_csv(paste0("output/rpareto_model_fits/true_rpareto_fits_model_mod_1.csv"),
                     col_names = F) %>% 
  .[nrow(.),] %>% # get last row in this file
  unlist %>% as.numeric()


alpha_plot_1 = res %>%
  ggplot()+
  geom_histogram(aes(alpha),fill = 'grey',col = 'black', bins = 20)+
  geom_vline(aes(xintercept = 12.47308), col = 'red', size = 1)+
  theme_minimal(12)+
  labs(x = expression(alpha),
       y = "Count")

my_mat = c()
for(i in seq(nrow(res))){
  my_mat = rbind(my_mat, 2 * (1 - pnorm(sqrt(variogram_model_score(distances, c(res[i,]$alpha, res[i,]$beta)) / 2))))
}


fit_model_1 = tibble(distances,
       true =  2 * (1 - pnorm(sqrt(variogram_model_score(distances, true_vals) / 2))),
       lower = my_mat %>% apply(MARGIN = 2, FUN = quantile, 0.025, na.rm = T),
       upper = my_mat %>% apply(MARGIN = 2, FUN = quantile, 0.975, na.rm = T)) %>%
  ggplot()+
  geom_line(aes(distances*100, true))+
  geom_segment(data = emp_chi, aes(x = dist*100, xend = dist*100, y = lower, yend = upper), alpha = 0.5)+
  geom_ribbon(aes(distances*100, ymin = lower, ymax = upper), alpha = 0.25)+
  theme_minimal(12)+
  theme(legend.position = 'none')+
  labs(x = "Distance (KM)",
       y = expression(chi[o]^F))+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))




# -------- MODEL 2
res = read_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_mod_2_L-BFGS-B.csv"),
               col_names = c("bts", "alpha", "beta")) %>%
  filter(bts>300)

true_vals = read_csv(paste0("output/rpareto_model_fits/true_rpareto_fits_model_mod_2.csv"),
                     col_names = F) %>% 
  .[nrow(.),] %>% # get last row in this file
  unlist %>% as.numeric()

alpha_plot_2 = res %>%
  ggplot()+
  geom_histogram(aes(alpha),fill = 'grey',col = 'black', bins = 20)+
  geom_vline(aes(xintercept = 14.93869), col = 'red', size = 1)+
  theme_minimal(12)+
  labs(x = expression(alpha),
       y = "Count")

my_mat = c()
for(i in seq(nrow(res))){
  my_mat = rbind(my_mat, 2 * (1 - pnorm(sqrt(variogram_model_score(distances, c(res[i,]$alpha, res[i,]$beta)) / 2))))
}

fit_model_2 = tibble(distances,
                     true =  2 * (1 - pnorm(sqrt(variogram_model_score(distances, true_vals) / 2))),
                     lower = my_mat %>% apply(MARGIN = 2, FUN = quantile, 0.025, na.rm = T),
                     upper = my_mat %>% apply(MARGIN = 2, FUN = quantile, 0.975, na.rm = T)) %>%
  ggplot()+
  geom_line(aes(distances*100, true))+
  geom_segment(data = emp_chi, aes(x = dist*100, xend = dist*100, y = lower, yend = upper), alpha = 0.5)+
  geom_ribbon(aes(distances*100, ymin = lower, ymax = upper), alpha = 0.25)+
  theme_minimal(12)+
  theme(legend.position = 'none')+
  labs(x = "Distance (KM)",
       y = expression(chi[o]^F))+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))




mod_1_and_2 = gridExtra::grid.arrange(fit_model_1, fit_model_2, nrow = 1)
ggsave(plot = mod_1_and_2, filename = "output/figs/rparto_mod_1_and_2_fit.pdf", height = 2.5, width = 8)



alpha_plot_1_2 = gridExtra::grid.arrange(alpha_plot_1, alpha_plot_2, nrow = 1)
ggsave(plot = alpha_plot_1_2, filename = "output/figs/rparto_alpha_mod_1_and.pdf", height = 2.5, width = 6)
