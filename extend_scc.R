scc_newpipeline = read.csv(file = paste0(file_path, "scc_newpipeline.csv"))
year_range_scc = range(scc_newpipeline$year)
year_min_scc = 2000
year_max_scc = 2500
lm_scc = lm(log(central) ~ year, data = scc_newpipeline)
scc_before = data.frame(year = year_min_scc:(year_range_scc[1] - 1)) %>%
  mutate(central = exp(predict(lm_scc, newdata = ., se.fit = T)$fit)) %>%
  mutate(low = central * 0.5, high = central * 1.5) %>%
  relocate(c(1, 3, 2, 4))
scc_after = data.frame(year = (year_range_scc[2] + 1):year_max_scc) %>%
  mutate(central = exp(predict(lm_scc, newdata = ., se.fit = T)$fit)) %>%
  mutate(low = central * 0.5, high = central * 1.5) %>%
  relocate(c(1, 3, 2, 4))
scc_extended = rbind(scc_before, scc_newpipeline, scc_after)
write.csv(scc_extended, file = paste0(file_path, "scc_extended.csv"))