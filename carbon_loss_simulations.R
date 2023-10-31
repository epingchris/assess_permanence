rm(list = ls())
library(ggplot2)
library(ggnewscale)
library(tidyverse)
library(magrittr)
library(fitdistrplus)
library(MASS)
library(EnvStats)
library(car)
library(tseries)
library(Matching)
library(mclust)
library(Hmisc)
file_path = "C:/Users/E-Ping Rau/OneDrive - University of Cambridge/carbon_release_pattern/"
source("/Users/E-Ping Rau/OneDrive - University of Cambridge/carbon_release_pattern/divergence.R")

# 1. Load existing data for sites ----
site = "CIF_Alto_Mayo" #Gola_country, WLT_VNCC_KNT, CIF_Alto_Mayo, VCS_1396, VCS_934
site_name = switch(site,
                   Gola = "Gola",
                   WLT_VNCC_KNT = "KNT",
                   CIF_Alto_Mayo = "Alto Mayo",
                   VCS_1396 = "RPA",
                   VCS_934 = "Mai Ndombe")

scc_new = read.csv(file = paste0(file_path, "scc_newpipeline.csv"))
scc_extrap = Hmisc::approxExtrap(scc_new$year, scc_new$central, 2000:2009) %>%
  as.data.frame() %>%
  `colnames<-`(c("year", "central")) %>%
  mutate(low = central * 0.5, high = central * 1.5) %>%
  relocate(c(1, 3, 2, 4)) %>%
  rbind(scc_new)
load(file = paste0(file_path, site, ".Rdata"))

# 2. Calculate carbon stock and flux ----
makeFlux = function(project_series, leakage_series){
  stock_series = aggregate(class_co2e ~ treatment + year, subset(project_series, class %in% eval_classes), FUN = sum)
  stock_wide = as.data.frame(pivot_wider(stock_series, id_cols = year, names_from = "treatment", values_from = "class_co2e"))
  
  flux_series = data.frame(year = stock_wide$year[-1],
                           treatment_proj = diff(subset(stock_series, treatment == "treatment")$class_co2e),
                           control_proj = diff(subset(stock_series, treatment == "control")$class_co2e)) %>%
    mutate(additionality = treatment_proj - control_proj)
  if(!is.null(leakage_series)){
    leakage_aggr = aggregate(class_co2e ~ treatment + year, subset(leakage_series, class %in% eval_classes), FUN = sum)
    flux_series %<>%
      mutate(treatment_leak = diff(subset(leakage_aggr, treatment == "treatment")$class_co2e),
             control_leak = diff(subset(leakage_aggr, treatment == "control")$class_co2e)) %>%
      mutate(leakage = treatment_leak - control_leak)
  }
  return(list(stock = stock_wide, flux = flux_series))
}

stock_series_sim = mapply(function(x, y) makeFlux(project_series = x, leakage_series = y)$stock,
                          x = agb_series_project_sim,
                          y = vector("list", length = length(agb_series_project_sim)),
                          SIMPLIFY = F)

flux_series_sim = mapply(function(x, y) makeFlux(project_series = x, leakage_series = y)$flux,
                         x = agb_series_project_sim,
                         y = vector("list", length = length(agb_series_project_sim)),
                         SIMPLIFY = F)

flux_series_sim_rel = mapply(function(x, y) {
  data.frame(flux_t1_treatment = x$treatment_proj,
             stock_t0_treatment = y[-nrow(y), ]$treatment,
             flux_t1_control = x$control_proj,
             stock_t0_control = y[-nrow(y), ]$control,
             year = x$year) %>%
    mutate(rel_flux_treatment = flux_t1_treatment / stock_t0_treatment,
           rel_flux_control = flux_t1_control / stock_t0_control)
}, x = flux_series_sim, y = stock_series_sim, SIMPLIFY = F)

summariseSeries = function(in_list, sel_col){
  df = as.data.frame(sapply(in_list, function(x) x[, sel_col])) %>%
    `colnames<-`(paste0("V", 1:20)) %>%
    rowwise() %>%
    mutate(mean = mean(V1:V20)) %>%
    ungroup() %>%
    mutate(year = in_list[[1]]$year, var = sel_col) %>%
    pivot_longer(V1:mean, names_to = "series", values_to = "val")
  return(df)
}

summ_stock = rbind(summariseSeries(stock_series_sim, "control"),
                   summariseSeries(stock_series_sim, "treatment"))

summ_flux = rbind(summariseSeries(flux_series_sim, "treatment_proj"),
                  summariseSeries(flux_series_sim, "control_proj"),
                  summariseSeries(flux_series_sim, "additionality"))

absloss_p_init = subset(summ_flux, var == "treatment_proj" & year >= t0 & series != "mean") %>%
  mutate(val = val * (-1), var = NULL, series = NULL)
absloss_c_init = subset(summ_flux, var == "control_proj" & year >= t0 & series != "mean") %>%
  mutate(val = val * (-1), var = NULL, series = NULL)


#test project-counterfactual carbon loss correlation and its impact on a-bar estimation based on different methods
cor.test(absloss_p_init$val, absloss_c_init$val, alternative = "greater", method = "p")
#highly significant for all except Alto Mayo where p-value is 0.06

add_obs = absloss_c_init$val - absloss_p_init$val
abar_obs = quantile(add_obs, 0.05)

absloss_p_fit = mclust::Mclust(absloss_p_init$val, 2)
absloss_c_fit = mclust::Mclust(absloss_c_init$val, 2)

abar_samp = rep(NA, 10000)

for(i in 1:10000){
  add_samp = SampGMM(absloss_c_fit$parameters, n = length(absloss_c_init$val)) -
    SampGMM(absloss_p_fit$parameters, n = length(absloss_p_init$val))
  abar_samp[i] = quantile(add_samp, 0.05)
}

hist(abar_samp, main = site_name, cex.main = 3, cex.axis = 2)
abline(v = abar_obs, lwd = 2, col = "red")
abline(v = median(abar_samp), lwd = 2)

# 2. Perform simulations (dynamic a-bar) ----
#function to generate sampled values from GMM-fitted distributions
SampGMM = function(params, n){
  if(length(params$variance$sigmasq) == 1) params$variance$sigmasq[2] = params$variance$sigmasq[1]
  samp_vec = NULL
  while(length(samp_vec) < n) {
    distr_chosen = ifelse(runif(1) <= params$pro[1], 1, 2)
    val = rnorm(1, mean = params$mean[distr_chosen], sd = sqrt(params$variance$sigmasq[distr_chosen]))
    if(val > 0) samp_vec = c(samp_vec, val)
  }
  return(samp_vec)
}

n_rep = 100 #number of repetitions
H = 50 #evaluation horizon; project duration
D = 0.03 #discount rate
bp = 5 #buffer period
#current_year = lubridate::year(lubridate::now()) #pretend every project starts at 2023
#scc_current = scc[which(names(scc) == current_year):length(scc)]

sim_p_loss = matrix(0, H, n_rep)
sim_c_loss = matrix(0, H, n_rep)
sim_additionality = matrix(0, H, n_rep)
sim_credit = matrix(0, H, n_rep)
sim_benefit = matrix(0, H, n_rep)
sim_abar = matrix(0, H, n_rep)
sim_release = matrix(0, H + 1, n_rep)
sim_damage = matrix(0, H, n_rep)
sim_ep = matrix(0, H, n_rep)
sim_credibility = matrix(1, H, n_rep)
sim_buffer = matrix(0, H, n_rep)
to_be_released = matrix(0, H, n_rep)
sim_r_sched = vector("list", n_rep)

for(j in 1:n_rep){
  r_sched = matrix(0, H, H + 1) #release schedule
  buffer_pool = 0
  for(i in 1:H){
    year_i = t0 + i - 1
    
    #get additionality: use ex post values in years where they are available, sample from fitted distributions otherwise
    if(year_i <= max(absloss_p_init$year)) {
      sim_p_loss[i, j] = sample(subset(absloss_p_init, year == year_i)$val, 1)
      sim_c_loss[i, j] = sample(subset(absloss_c_init, year == year_i)$val, 1)
      sim_additionality[i, j] = sim_c_loss[i, j] - sim_p_loss[i, j]
      
      #calculate a-bar based on carbon loss distributions
      absloss_p_updated = subset(absloss_p_init, year <= year_i)$val
      absloss_c_updated = subset(absloss_c_init, year <= year_i)$val
      
      absloss_p_fit = mclust::Mclust(absloss_p_updated, 2)
      absloss_c_fit = mclust::Mclust(absloss_c_updated, 2)
      
      samp_additionality = SampGMM(absloss_c_fit$parameters, n = length(absloss_c_updated)) -
        SampGMM(absloss_p_fit$parameters, n = length(absloss_p_updated))
      abar = quantile(samp_additionality, 0.05)
      sim_abar[i, j] = abar
    } else if(year_i > max(absloss_p_init$year)){
      sim_p_loss[i, j] = SampGMM(absloss_p_fit$parameters, n = 1)
      sim_c_loss[i, j] = SampGMM(absloss_c_fit$parameters, n = 1)
      sim_additionality[i, j] = sim_c_loss[i, j] - sim_p_loss[i, j]
      sim_abar[i, j] = abar
    }
    
    if(i <= bp) {
      #first five years: no credits/releases; positive additionality added to buffer pool
      if(sim_additionality[i, j] > 0) buffer_pool = buffer_pool + sim_additionality[i, j]
      sim_credit[i, j] = 0
    } else {
      #from sixth year on: get credits and anticipated releases

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
        sim_benefit[i, j] = sim_credit[i, j] * filter(scc_extrap, year == year_i)$central
        k = i #kth year(s), for which we estimate anticipated release
        while(to_be_released[i, j] > 0){
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
        sim_damage[i, j] = sum(r_sched[i, ] * filter(scc_extrap, year %in% t0:(t0 + H))$central / ((1 + D) ^ (1:(H + 1))))
        sim_ep[i, j] = (sim_benefit[i, j] - sim_damage[i, j]) / sim_benefit[i, j]
      } else if(sim_credit[i, j] <= 0){
        sim_ep[i, j] = 0
        sim_credibility[i, j] = 0
      }
    }
    sim_buffer[i, j] = buffer_pool
    sim_r_sched[[j]] = r_sched
  }
}

# 4. Visualise results ----
SummariseSim = function(mat){
  df = mat %>%
    as.data.frame() %>%
    reframe(
      year = row_number(),
      p05 = apply(., 1, function(x) quantile(x, 0.05, na.rm = T)),
      p25 = apply(., 1, function(x) quantile(x, 0.25, na.rm = T)),
      median = apply(., 1, median, na.rm = T),
      p75 = apply(., 1, function(x) quantile(x, 0.75, na.rm = T)),
      p95 = apply(., 1, function(x) quantile(x, 0.95, na.rm = T))
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

ggplot(as.data.frame(sim_credit) %>% mutate(t = row_number()) %>% pivot_longer(V1:V100, names_to = "rep", values_to = "val")) +
  geom_line(aes(x = t, y = val, color = rep), size = 0.5, show.legend = F) +
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
ggsave(paste0(file_path, site, "_6f_sim_credit_real.png"), width = 15, height = 10, unit = "cm")


ggplot(summary_credit, aes(x = year)) +
  geom_ribbon(data = . %>% dplyr::select(year, p05, p95), aes(ymin = p05, ymax = p95), fill = "lightgray", alpha = 0.5) +
  geom_ribbon(data = . %>% dplyr::select(year, p25, p75), aes(ymin = p25, ymax = p75), fill = "darkgray", alpha = 0.5) +
  geom_line(data = . %>% dplyr::select(year, median), aes(y = median), size = 0.5, color = "red") +
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

ggplot(as.data.frame(sim_release[1:H, ]) %>% mutate(t = row_number()) %>% pivot_longer(V1:V100, names_to = "rep", values_to = "val")) +
  geom_line(aes(x = t, y = val, col = rep), size = 0.5, show.legend = F) +
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
ggsave(paste0(file_path, site, "_6g_sim_release_real.png"), width = 15, height = 10, unit = "cm")

ggplot(summary_release, aes(x = year)) +
  geom_ribbon(data = . %>% dplyr::select(year, p05, p95), aes(ymin = p05, ymax = p95), fill = "lightgray", alpha = 0.5) +
  geom_ribbon(data = . %>% dplyr::select(year, p25, p75), aes(ymin = p25, ymax = p75), fill = "darkgray", alpha = 0.5) +
  geom_line(data = . %>% dplyr::select(year, median), aes(y = median), size = 0.5, color = "red") +
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

ggplot(as.data.frame(sim_ep) %>% replace(., . == 0, NA) %>% mutate(t = row_number()) %>% pivot_longer(V1:V100, names_to = "rep", values_to = "val")) +
  geom_line(aes(x = t, y = val, col = rep), size = 0.5, show.legend = F) +
  geom_hline(yintercept = 0, size = 1, color = "black") +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) + 
  scale_y_continuous(name = "eP", limits = c(0, 1)) + 
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, site, "_6h_sim_ep_real.png"), width = 15, height = 10, unit = "cm")

ggplot(summary_ep, aes(x = year)) +
  geom_ribbon(data = . %>% dplyr::select(year, p05, p95), aes(ymin = p05, ymax = p95), fill = "lightgray", alpha = 0.5) +
  geom_ribbon(data = . %>% dplyr::select(year, p25, p75), aes(ymin = p25, ymax = p75), fill = "darkgray", alpha = 0.5) +
  geom_line(data = . %>% dplyr::select(year, median), aes(y = median), size = 0.5, color = "red") +
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
  geom_hline(yintercept = 5, size = 1, color = "black", lty = "dashed") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, site, "_6i_sim_cred_rep.png"), width = 15, height = 10, unit = "cm")
