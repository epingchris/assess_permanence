rm(list = ls())
library(tidyverse)
library(magrittr)
library(mclust)
file_path = "C:/Users/E-Ping Rau/OneDrive - University of Cambridge/carbon_release_pattern/"


# Simulation settings ----

#select site/types (single sites, portfolio or simulated ones)
site = "expo" #Gola_country, WLT_VNCC_KNT, CIF_Alto_Mayo, VCS_1396, VCS_934, expo, portfolio
portfolio_type = "good" #all, good
scale_c = 10  #1.5, 2, 5, 10
use_aomega_theo = T #T, F
save_type = "pdf" #png, pdf

#basic parameters
omega = 0.05

#parameters for the simulations
n_rep = 100 #number of repetitions
H = 50 #evaluation horizon; project duration
D = 0.03 #discount rate
bp = 5 #buffer period


# SCC ----
scc_new = read.csv(file = paste0(file_path, "scc_newpipeline.csv"))
lm_scc = lm(log(central) ~ year, data = scc_new)
scc_extrap1 = data.frame(year = 2000:2009)
scc_extrap2 = data.frame(year = 2101:2500)
scc_extrap1$central = exp(predict(lm_scc, newdata = scc_extrap1, se.fit = T)$fit)
scc_extrap2$central = exp(predict(lm_scc, newdata = scc_extrap2, se.fit = T)$fit)
scc_extrap = rbind(scc_extrap1, dplyr::select(scc_new, c("year", "central")), scc_extrap2) %>%
  mutate(low = central * 0.5, high = central * 1.5) %>%
  relocate(c(1, 3, 2, 4))
year_max = max(scc_extrap$year)


# Functions ----
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

#wrapper function to fit using GMM and generate random samples from GMM-fitted distributions
FitGMM = function(x) {
  if(length(x) > 0) {
    return(mclust::Mclust(x, 2, verbose = F))
  } else {
    return(NULL)
  }
}

SampGMM = function(mclust_obj, n) {
  if(!is.null(mclust_obj)) {
    params = mclust_obj$param
  } else {
    return(NA)
  }
  if(length(params$variance$sigmasq) == 1) params$variance$sigmasq[2] = params$variance$sigmasq[1]
  samp_vec = NULL
  while(length(samp_vec) < n) {
    distr_chosen = ifelse(runif(1) <= params$pro[1], 1, 2)
    val = rnorm(1, mean = params$mean[distr_chosen], sd = sqrt(params$variance$sigmasq[distr_chosen]))
    if(val > 0) samp_vec = c(samp_vec, val)
  }
  return(samp_vec)
}

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


# Core code for the simulation ----
source(paste0(file_path, "carbon_loss_simulations_core.R"))


# Plot results ----
yr_label = c(1, seq(H / 10, H, by = H / 10))

## Figure 3. a-omega distribution ----
expo_limits = case_when(
  site == "expo" & scale_c <= 1.5 ~ c(-10, 10),
  .default = NULL)

ggplot(data = data.frame(var = "add", val = add_samp), aes(x = val)) +
  geom_freqpoly() +
  geom_vline(xintercept = 0, lwd = 0.5, col = "gray") +
  geom_vline(xintercept = ifelse(use_aomega_theo, aomega_theo, aomega_samp), col = "red", lwd = 1) +
  scale_x_continuous(name = "Additionality (Mg CO2e)") + 
  scale_y_continuous(name = "Counts") + 
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, site_name, "_3_aomega_distribution.", save_type), width = 15, height = 15, unit = "cm")


## Figure 4. Credits time series ----
if(site == "expo") {
  range_credit = c(-5, 45)
  range_release = c(0, 5)
} else {
  range_credit = c(min(summ_credit$p05) - 1000, max(summ_credit$p95) + 1000)
  range_release = c(0, max(summ_release$p95) + 1000)
}

# ggplot(sim_credit_long) +
#   geom_line(aes(x = t, y = val, color = rep), size = 0.5, show.legend = F) +
#   geom_hline(yintercept = 0, size = 1, color = "black") +
#   geom_vline(xintercept = year_pres_obs - t0 + 1, lty = "dotted") +
#   scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) + 
#   scale_y_continuous(name = "Credits (Mg CO2 e)", limits = range_credit) + 
#   theme_bw() +
#   theme(axis.line = element_line(linewidth = 0.5),
#         panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 24),
#         axis.text.x = element_text(size = 18, vjust = 0.5),
#         axis.text.y = element_text(size = 18, angle = 45),
#         plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
# ggsave(paste0(file_path, site_name, "_4a_sim_credit_real.", save_type), width = 15, height = 10, unit = "cm")

ggplot(summ_credit, aes(x = year)) +
  geom_ribbon(data = . %>% dplyr::select(year, p05, p95), aes(ymin = p05, ymax = p95), fill = "lightgray", alpha = 0.5) +
  geom_ribbon(data = . %>% dplyr::select(year, p25, p75), aes(ymin = p25, ymax = p75), fill = "darkgray", alpha = 0.5) +
  geom_line(data = . %>% dplyr::select(year, median), aes(y = median), size = 0.5, color = "red") +
  geom_hline(yintercept = 0, size = 1, color = "black") +
  geom_vline(xintercept = year_pres_obs - t0 + 1, lty = "dotted") +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) + 
  scale_y_continuous(name = "Credits (Mg CO2 e)", limits = range_credit) + 
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, site_name, "_4_sim_credit_rep.", save_type), width = 15, height = 10, unit = "cm")


# ## Figure 5. Anticipated release time series ----
# ggplot(subset(sim_release_long, t <= H)) +
#   geom_line(aes(x = t, y = val, col = rep), size = 0.5, show.legend = F) +
#   geom_hline(yintercept = 0, size = 1, color = "black") +
#   geom_vline(xintercept = year_pres_obs - t0 + 1, lty = "dotted") +
#   scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) + 
#   scale_y_continuous(name = "Anticipated releases\n(Mg CO2 e)", limits = range_release) + 
#   theme_bw() +
#   theme(axis.line = element_line(linewidth = 0.5),
#         panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 24),
#         axis.text.x = element_text(size = 18, vjust = 0.5),
#         axis.text.y = element_text(size = 18, angle = 45),
#         plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
# ggsave(paste0(file_path, site_name, "_5a_sim_release_real.", save_type), width = 15, height = 10, unit = "cm")

ggplot(summ_release, aes(x = year)) +
  geom_ribbon(data = . %>% dplyr::select(year, p05, p95), aes(ymin = p05, ymax = p95), fill = "lightgray", alpha = 0.5) +
  geom_ribbon(data = . %>% dplyr::select(year, p25, p75), aes(ymin = p25, ymax = p75), fill = "darkgray", alpha = 0.5) +
  geom_line(data = . %>% dplyr::select(year, median), aes(y = median), size = 0.5, color = "red") +
  geom_hline(yintercept = 0, size = 1, color = "black") +
  geom_vline(xintercept = year_pres_obs - t0 + 1, lty = "dotted") +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) + 
  scale_y_continuous(name = "Anticipated releases\n(Mg CO2 e)", limits = range_release) + 
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(1, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, site_name, "_5_sim_release_rep.", save_type), width = 15, height = 10, unit = "cm")


## Figure 6. EP time series ----
# ggplot(as.data.frame(sim_ep) %>% replace(., . == 0, NA) %>% mutate(t = row_number()) %>% pivot_longer(V1:V100, names_to = "rep", values_to = "val")) +
#   geom_line(aes(x = t, y = val, col = rep), size = 0.5, show.legend = F) +
#   geom_hline(yintercept = 0, size = 1, color = "black") +
#   geom_vline(xintercept = year_pres_obs - t0 + 1, lty = "dotted") +
#   scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) + 
#   scale_y_continuous(name = "EP", limits = c(0, 1)) + 
#   theme_bw() +
#   theme(axis.line = element_line(linewidth = 0.5),
#         panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 24),
#         axis.text.x = element_text(size = 18, vjust = 0.5),
#         axis.text.y = element_text(size = 18, angle = 45),
#         plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
# ggsave(paste0(file_path, site_name, "_6a_sim_ep_real.", save_type), width = 15, height = 10, unit = "cm")

ggplot(summ_ep, aes(x = year)) +
  geom_ribbon(data = . %>% dplyr::select(year, p05, p95), aes(ymin = p05, ymax = p95), fill = "lightgray", alpha = 0.5) +
  geom_ribbon(data = . %>% dplyr::select(year, p25, p75), aes(ymin = p25, ymax = p75), fill = "darkgray", alpha = 0.5) +
  geom_line(data = . %>% dplyr::select(year, median), aes(y = median), size = 0.5, color = "red") +
  geom_vline(xintercept = year_pres_obs - t0 + 1, lty = "dotted") +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) + 
  scale_y_continuous(name = "EP", limits = c(0, 1)) + 
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, site_name, "_6_sim_ep_rep.", save_type), width = 15, height = 10, unit = "cm")


## Figure 7. Credibility time series ----
ggplot(summ_cred, aes(x = year)) +
  geom_line(aes(y = cred), size = 1) +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) + 
  scale_y_continuous(name = "Credibility", limits = c(0.5, 1)) + 
  geom_hline(yintercept = 0.95, color = "red", lty = "dashed", size = 0.5) +
  geom_vline(xintercept = case_when(site != "expo" ~ year_pres_obs - t0 + 1, .default = NULL), lty = "dotted") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, site_name, "_7_sim_cred_rep.", save_type), width = 15, height = 10, unit = "cm")


# Transition of simulation results along 1/lambda_c ----
transition = F

if(transition) {
  scale_c_vec = c(seq(1.1, 5, by = 0.1), seq(6, 15, by = 1))
  
  mean_credit = matrix(NA, n_rep, length(scale_c_vec))
  max_EP = matrix(NA, n_rep, length(scale_c_vec))
  where_max_EP = matrix(NA, n_rep, length(scale_c_vec))
  final_EP = matrix(NA, n_rep, length(scale_c_vec))
  mean_cred = matrix(NA, n_rep, length(scale_c_vec))
  
  for(i_val in 1:length(scale_c_vec)){
    scale_c = scale_c_vec[i_val]
    source(paste0(file_path, "carbon_loss_simulations_core.R"))
    
    #relationship between 1 / lambda_c and mean of [median of credits] from all years
    mean_credit[, i_val] = apply(sim_credit, 2, sum) / (H - bp)
    
    #relationship between 1 / lambda_c and max EP from all years
    max_EP[, i_val] = apply(sim_ep, 2, max)
    
    #relationship between 1 / lambda_c and year where max EP occurs
    where_max_EP[, i_val] = apply(sim_ep, 2, which.max)
    
    #relationship between 1 / lambda_c and final EP
    final_EP[, i_val] = tail(sim_ep, 1)
    
    #relationship between 1 / lambda_c and credibility (proportion of years without credibility incident)
    mean_cred[, i_val] = apply(sim_ep[-c(1:bp), ], 2, function(x) length(which(x > 0)) / (H - bp))
  }
  
  SummariseAcrossTests = function(input, var) {
    n = nrow(input)
    output = input %>% 
      t() %>%
      as.data.frame() %>%
      rowwise() %>% 
      mutate(mean = mean(c_across(1:n_rep)),
             sd = sd(c_across(1:n_rep))) %>%
      ungroup() %>%
      mutate(var = var,
             margin = qt(0.975, df = n_rep - 1) * sd / sqrt(n_rep),
             val_low = mean - margin,
             val_high = mean + margin) %>%
      select(-all_of(c(1:n)))
    return(output)
  }
  
  ## Figure 8.
  data_credit = SummariseAcrossTests(mean_credit, scale_c_vec)
  ggplot(data = data_credit, aes(x = var)) +
    geom_ribbon(aes(ymin = val_low, ymax = val_high), fill = "lightgray") +
    geom_line(aes(y = mean), size = 1) +
    scale_x_continuous(name = expression(1 / lambda[c])) + 
    scale_y_continuous(name = "Mean credits") + 
    theme_bw() +
    theme(axis.line = element_line(linewidth = 0.5),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 18, vjust = 0.5),
          axis.text.y = element_text(size = 18, angle = 45),
          plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
  ggsave(paste0(file_path, "expo_8a_lc_vs_mean_credit.", save_type), width = 15, height = 10, unit = "cm")
  
  data_max_EP = SummariseAcrossTests(max_EP, scale_c_vec)
  ggplot(data = data_max_EP, aes(x = var)) +
    geom_ribbon(aes(ymin = val_low, ymax = val_high), fill = "lightgray") +
    geom_line(aes(y = mean), size = 1) +
    scale_x_continuous(name = expression(1 / lambda[c])) + 
    scale_y_continuous(name = "Maximum EP") + 
    theme_bw() +
    theme(axis.line = element_line(linewidth = 0.5),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 18, vjust = 0.5),
          axis.text.y = element_text(size = 18, angle = 45),
          plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
  ggsave(paste0(file_path, "expo_8b_lc_vs_max_EP.", save_type), width = 15, height = 10, unit = "cm")
  
  data_where_max_EP = SummariseAcrossTests(where_max_EP, scale_c_vec)
  ggplot(data = data_where_max_EP, aes(x = var)) +
    geom_ribbon(aes(ymin = val_low, ymax = val_high), fill = "lightgray") +
    geom_line(aes(y = mean), size = 1) +
    scale_x_continuous(name = expression(1 / lambda[c])) + 
    scale_y_continuous(name = "Year of maximum EP") + 
    theme_bw() +
    theme(axis.line = element_line(linewidth = 0.5),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 18, vjust = 0.5),
          axis.text.y = element_text(size = 18, angle = 45),
          plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
  ggsave(paste0(file_path, "expo_8c_lc_vs_where_max_EP.", save_type), width = 15, height = 10, unit = "cm")
  
  data_final_EP = SummariseAcrossTests(final_EP, scale_c_vec)
  ggplot(data = data_final_EP, aes(x = var)) +
    geom_ribbon(aes(ymin = val_low, ymax = val_high), fill = "lightgray") +
    geom_line(aes(y = mean), size = 1) +
    scale_x_continuous(name = expression(1 / lambda[c])) + 
    scale_y_continuous(name = "EP at project end") + 
    theme_bw() +
    theme(axis.line = element_line(linewidth = 0.5),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 18, vjust = 0.5),
          axis.text.y = element_text(size = 18, angle = 45),
          plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
  ggsave(paste0(file_path, "expo_8d_lc_vs_final_EP.", save_type), width = 15, height = 10, unit = "cm")
  
  data_mean_cred = SummariseAcrossTests(mean_cred, scale_c_vec)
  ggplot(data = data_mean_cred, aes(x = var)) +
    geom_ribbon(aes(ymin = val_low, ymax = val_high), fill = "lightgray") +
    geom_line(aes(y = mean), size = 1) +
    geom_hline(yintercept = 0.95, col = "red", lty = "dashed", size = 0.5) +
    scale_x_continuous(name = expression(1 / lambda[c])) + 
    scale_y_continuous(name = "Mean credibility", limits = c(0.75, 1)) + 
    theme_bw() +
    theme(axis.line = element_line(linewidth = 0.5),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 18, vjust = 0.5),
          axis.text.y = element_text(size = 18, angle = 45),
          plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
  ggsave(paste0(file_path, "expo_8e_lc_vs_no_cred.", save_type), width = 15, height = 10, unit = "cm")

  data_mean_cred = SummariseAcrossTests(mean_cred, scale_c_vec)
  ggplot(data = data_mean_cred, aes(x = var)) +
    geom_ribbon(aes(ymin = val_low, ymax = val_high), fill = "lightgray") +
    geom_line(aes(y = mean), size = 1) +
    geom_hline(yintercept = 0.95, col = "red", lty = "dashed", size = 0.5) +
    scale_x_continuous(name = expression(1 / lambda[c]), limits = c(1, 2)) + 
    scale_y_continuous(name = "Mean credibility", limits = c(0.75, 1)) + 
    theme_bw() +
    theme(axis.line = element_line(linewidth = 0.5),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 18, vjust = 0.5),
          axis.text.y = element_text(size = 18, angle = 45),
          plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
  ggsave(paste0(file_path, "expo_8e_lc_vs_no_cred_short.", save_type), width = 15, height = 10, unit = "cm")
}

  
