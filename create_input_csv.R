#This housekeeping script is specific to an older 4C workflow, which generates input dataset in the format that SimulatePermanence() needs
#from results of an older evaluation code, stored in .Rdata. For the projects used in the paper, this code was already run once and doesn't need to be run again.
rm(list = ls())

#function used to make carbon flux data frame with format needed by SimulatePermanence() input
makeFlux = function(in_df, t0){
  n_sim = unique(in_df$n_sim)
  flux_series = lapply(n_sim, function(x) {
    stock_series = aggregate(class_co2e ~ treatment + year, subset(in_df, n_sim == x & class %in% c(1, 2, 3, 4)), FUN = sum)
    stock_wide = as.data.frame(pivot_wider(stock_series, id_cols = year, names_from = "treatment", values_from = "class_co2e"))
    
    flux_series = data.frame(year = stock_wide$year[-1],
                             project = diff(subset(stock_series, treatment == "treatment")$class_co2e),
                             counterfactual = diff(subset(stock_series, treatment == "control")$class_co2e)) %>%
      mutate(additionality = project - counterfactual) %>%
      pivot_longer(2:4, names_to = "var", values_to = "val") %>%
      mutate(n_sim = x)
  }) %>%
    do.call(rbind, .) %>%
  mutate(started = (year >= t0))
  return(flux_series)
}

project_site = "VCS_934" #Gola_country, CIF_Alto_Mayo, VCS_1396, VCS_934
load(file = paste0("project_input_data/", project_site, ".Rdata")) #t0, eval_classes included

agb_series_df = lapply(seq_along(agb_series_project_sim), function(i) {
  agb_series_project_sim[[i]] %>% mutate(n_sim = i)
}) %>%
  do.call(rbind, .)
flux_series = makeFlux(agb_series_df, t0 = t0)
write.csv(flux_series, paste0("project_input_data/", project_site, ".csv"), row.names = F)

#save a the essential data in a cleaner format in case we need to re-run in the future
saveRDS(list(agb_series_list = agb_series_project_sim, t0 = t0),
        file = paste0("project_input_data/", project_site, ".rds"))

#this part is for future re-runs
project_site = "Gola_country" #Gola_country, CIF_Alto_Mayo, VCS_1396, VCS_934
project_data = readRDS(paste0("project_input_data/", project_site, ".rds"))
agb_series_df = lapply(seq_along(project_data$agb_series_list), function(i) {
  agb_series_project_sim[[i]] %>% mutate(n_sim = i)
}) %>%
  do.call(rbind, .)
flux_series = makeFlux(agb_series_df, t0 = project_data$t0)
write.csv(flux_series, paste0("project_input_data/", project_site, "2.csv"), row.names = F)
