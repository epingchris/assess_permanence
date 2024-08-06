project_site = "VCS_934" #Gola_country, WLT_VNCC_KNT, CIF_Alto_Mayo, VCS_1396, VCS_934

load(file = paste0(file_path, "project_input_data/", project_site, ".Rdata")) #t0, eval_classes included

agb_series_project_sim_df = lapply(seq_along(agb_series_project_sim), function(i) {
  agb_series_project_sim[[i]] %>%
    mutate(n_sim = i)
}) %>%
  do.call(rbind, .) %>%
  mutate(started = (year >= t0))
write.csv(agb_series_project_sim_df, paste0(file_path, "project_input_data/", project_site, ".csv"), row.names = F)
