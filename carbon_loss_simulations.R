absloss_p_init = subset(summ_flux, var == "treatment_proj" & year >= t0 & series != "mean") %>%
  mutate(val = val * (-1), var = NULL, series = NULL)
absloss_c_init = subset(summ_flux, var == "control_proj" & year >= t0 & series != "mean") %>%
  mutate(val = val * (-1), var = NULL, series = NULL)

# additionality = absoloss_c_update$val - absoloss_p_update$val
# samp_additionality = SampByPro(absloss_c_fit$parameters, n = length(absoloss_c_update$val)) -
#   SampByPro(absloss_p_fit$parameters, n = length(absoloss_p_update$val))
# abar = quantile(additionality, 0.05)
# samp_abar = quantile(samp_additionality, 0.05)

##Perform simulations (dynamic a-bar) ----
n_rep = 100 #number of repetitions
H = 50 #evaluation horizon; project duration
H_rel = 10 * H #release horizon
D = 0.03 #discount rate
bp = 5 #buffer period
current_year = lubridate::year(lubridate::now()) #pretend every project starts at 2023
scc_current = scc[which(names(scc) == current_year):length(scc)]

sim_p_loss = matrix(0, H, n_rep)
sim_c_loss = matrix(0, H, n_rep)
sim_additionality = matrix(0, H, n_rep)
sim_credit = matrix(0, H, n_rep)
sim_benefit = matrix(0, H, n_rep)
sim_release = matrix(0, H_rel, n_rep)
sim_damage = matrix(0, H, n_rep)
sim_ep = matrix(0, H, n_rep)
sim_credibility = matrix(1, H, n_rep)
sim_r_sched = vector("list", n_rep)
to_be_released = matrix(0, H, n_rep)

a = Sys.time()
for(j in 1:n_rep){

  r_sched = matrix(0, H, H_rel) #release schedule
  buffer_pool = 0
  for(i in 1:H){
    year_i = t0 + i - 1
    if(year_i <= max(absloss_p_init$year)) {
      sim_p_loss[i, j] = mean(subset(absloss_p_init, year == year_i)$val)
      sim_c_loss[i, j] = mean(subset(absloss_c_init, year == year_i)$val)
      sim_additionality[i, j] = sim_c_loss[i, j] - sim_p_loss[i, j]
      
      if(i <= bp) { #first five years
        if(sim_additionality[i, j] > 0) buffer_pool = buffer_pool + sim_additionality[i, j]
        sim_credit[i, j] = 0
      } else { #sixth year to end of available ex post value
        absloss_p_update = subset(absloss_p_init, year <= year_i)$val
        absloss_c_update = subset(absloss_c_init, year <= year_i)$val
        
        absloss_p_fit = mclust::Mclust(absloss_p_update, 2)
        absloss_c_fit = mclust::Mclust(absloss_c_update, 2)
        
        samp_additionality = SampByPro(absloss_c_fit$parameters, n = length(absloss_c_update)) -
          SampByPro(absloss_p_fit$parameters, n = length(absloss_p_update))
        abar = quantile(samp_additionality, 0.05)
        
        #use buffer pool to fill anticipated releases first
        #only deduct from buffer pool at each year if there is space left for that year
        if(buffer_pool > 0) {
          max_release = ifelse(abar > 0, 0, -abar) #if a-bar is positive, maximum release is zero
          can_be_released = max(0, max_release - sim_release[i, j])
          sim_release[i, j] = sim_release[i, j] + can_be_released
          buffer_pool = buffer_pool - can_be_released
        }

        sim_credit[i, j] = sim_additionality[i, j] + sim_release[i, j]
        if(sim_credit[i, j] > 0){
          to_be_released[i, j] = sim_credit[i, j]
          sim_benefit[i, j] = sim_credit[i, j] * scc_current[i]
          k = i #number of years after i
          while(to_be_released[i, j] > 0 & k < H_rel){
            k = k + 1
            if(k > H) {
              can_be_released = to_be_released[i, j]
              #after project ends, all remaining credits are released the next year
            } else {
              max_release = ifelse(abar > 0, 0, -abar) #if a-bar is positive, maximum release is zero
              can_be_released = max(0, max_release - sim_release[k, j])
            }
            r_sched[i, k] = min(to_be_released[i, j], can_be_released)
            to_be_released[i, j] = to_be_released[i, j] - r_sched[i, k]
            sim_release[k, j] = sim_release[k, j] + r_sched[i, k]
          }
          sim_damage[i, j] = sum(r_sched[i, ] * scc_current[1:length(r_sched[i, ])] / ((1 + D) ^ (1:length(r_sched[i, ]))))
          sim_ep[i, j] = (sim_benefit[i, j] - sim_damage[i, j]) / sim_benefit[i, j]
        } else if(sim_credit[i, j] < 0){
          sim_ep[i, j] = 0
          sim_credibility[i, j] = 0
        }
      }
    } else {
      absloss_p_fit = mclust::Mclust(absloss_p_update, 2)
      absloss_c_fit = mclust::Mclust(absloss_c_update, 2)
      
      sim_p_loss[i, j] = SampByPro(absloss_p_fit$parameters, n = 1)
      sim_c_loss[i, j] = SampByPro(absloss_c_fit$parameters, n = 1)
      absloss_p_update = c(absloss_p_update, sim_p_loss[i, j])
      absloss_c_update = c(absloss_c_update, sim_c_loss[i, j])
      
      sim_additionality[i, j] = sim_c_loss[i, j] - sim_p_loss[i, j]
      sim_credit[i, j] = sim_additionality[i, j] + sim_release[i, j]
      if(sim_credit[i, j] > 0){
        to_be_released[i, j] = sim_credit[i, j]
        sim_benefit[i, j] = sim_credit[i, j] * scc_current[i]
        k = i #number of years after i
        while(to_be_released[i, j] > 0 & k < H_rel){
          k = k + 1
          if(k > H) {
            max_release = mean(sim_release[, j][which(sim_release[, j] != 0)]) #after project ends, release at average release rate during project
          } else {
            max_release = ifelse(abar > 0, 0, -abar) #if a-bar is positive, maximum release is zero
          }
          can_be_released = max(0, max_release - sim_release[k, j])
          r_sched[i, k] = min(to_be_released[i, j], can_be_released)
          to_be_released[i, j] = to_be_released[i, j] - r_sched[i, k]
          sim_release[k, j] = sim_release[k, j] + r_sched[i, k]
        }
        sim_damage[i, j] = sum(r_sched[i, ] * scc_current[1:length(r_sched[i, ])] / ((1 + D) ^ (1:length(r_sched[i, ]))))
        sim_ep[i, j] = (sim_benefit[i, j] - sim_damage[i, j]) / sim_benefit[i, j]
      } else if(sim_credit[i, j] < 0){
        sim_ep[i, j] = 0
        sim_credibility[i, j] = 0
      }
    }
  sim_r_sched[[j]] = r_sched
  }
}
b = Sys.time()
b - a

SummariseSim = function(mat){
  df = mat %>%
    as.data.frame() %>%
    reframe(
      year = row_number(),
      min = apply(., 1, min, na.rm = T),
      max = apply(., 1, max, na.rm = T),
      median = apply(., 1, median, na.rm = T),
      q1 = apply(., 1, function(x) quantile(x, 0.25, na.rm = T)),
      q3 = apply(., 1, function(x) quantile(x, 0.75, na.rm = T))
    )
  return(df)
}

summary_credit = SummariseSim(sim_credit)
summary_release = SummariseSim(sim_release[1:H, ])
summary_ep = SummariseSim(sim_ep %>% replace(., . == 0, NA)) %>%
  replace(., . == Inf| . == -Inf, NA)

summary_cred = sim_ep %>%
  as.data.frame() %>%
  reframe(
    year = row_number(),
    cred = apply(., 1, function(x) length(which(x == 0))))
summary_cred$cred[1:5] = NA
yr_label = c(1, seq(5, 50, by = 5))

ggplot(summary_credit, aes(x = year)) +
  geom_ribbon(data = . %>% dplyr::select(year, min, max), aes(ymin = min, ymax = max), fill = "lightgray", alpha = 0.5) +
  geom_ribbon(data = . %>% dplyr::select(year, q1, q3), aes(ymin = q1, ymax = q3), fill = "darkgray", alpha = 0.5) +
  geom_line(data = . %>% dplyr::select(year, median), aes(y = median), size = 1, color = "red") +
  geom_hline(yintercept = 0, size = 1, color = "black") +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) + 
  scale_y_continuous(name = "Credits (Mg CO2 e)") + 
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, site, "_6f_sim_credit_rep.png"), width = 15, height = 10, unit = "cm")

ggplot(summary_release, aes(x = year)) +
  geom_ribbon(data = . %>% dplyr::select(year, min, max), aes(ymin = min, ymax = max), fill = "lightgray", alpha = 0.5) +
  geom_ribbon(data = . %>% dplyr::select(year, q1, q3), aes(ymin = q1, ymax = q3), fill = "darkgray", alpha = 0.5) +
  geom_line(data = . %>% dplyr::select(year, median), aes(y = median), size = 1, color = "red") +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) + 
  scale_y_continuous(name = "Anticipated releases\n(Mg CO2 e)") + 
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(1, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, site, "_6g_sim_release_rep.png"), width = 15, height = 10, unit = "cm")

ggplot(summary_ep, aes(x = year)) +
  geom_ribbon(data = . %>% dplyr::select(year, min, max), aes(ymin = min, ymax = max), fill = "lightgray", alpha = 0.5) +
  geom_ribbon(data = . %>% dplyr::select(year, q1, q3), aes(ymin = q1, ymax = q3), fill = "darkgray", alpha = 0.5) +
  geom_line(data = . %>% dplyr::select(year, median), aes(y = median), size = 1, color = "red") +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) + 
  scale_y_continuous(name = "eP", limits = c(0, 1)) + 
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, site, "_6h_sim_ep_rep.png"), width = 15, height = 10, unit = "cm")

ggplot(summary_cred, aes(x = year)) +
  geom_line(aes(y = cred), size = 1, color = "red") +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) + 
  scale_y_continuous(name = "Nb. of credibility incidents", limits = c(0, 50)) + 
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, site, "_6i_sim_cred_rep.png"), width = 15, height = 10, unit = "cm")
