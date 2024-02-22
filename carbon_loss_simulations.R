rm(list = ls())
library(tidyverse)
library(magrittr)
library(mclust)
file_path = "C:/Users/epr26/OneDrive - University of Cambridge/carbon_release_pattern/"

# Simulation settings ----

#set simulation type (hypothetical or real-life, single or aggregated)
type = "portfolio" #expo, project, expo_portfolio, portfolio
scale_c = 1.1  #1.1, 2, 5
expo_portfolio_type = "C" #A, B, C
portfolio_type = "good" #five, four, good
project_site = "VCS_934" #Gola_country, WLT_VNCC_KNT, CIF_Alto_Mayo, VCS_1396, VCS_934
use_theo = T #T, F; type == expo_portfolio will override this and set use_theo as F
file_type = "png" #png, pdf
view_snapshot = F
print_real = F #print real trajectories of all 100 repetitions
lc_sensitivity = F
bp_sensitivity = F
ppr_sensitivity = F
H_sensitivity = F

#basic parameters
omega = 0.05

#parameters for the simulations
n_rep = 100 #number of repetitions
H = 50 #evaluation horizon; project duration
D = 0.03 #discount rate
bp = 5 #buffer period
rate_postproj = 2 #rate of post-project release compare to during project


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


# Core code for the simulation ----
source(paste0(file_path, "carbon_loss_simulations_core.R"))


# Plot results ----
yr_label = c(1, seq(H / 10, H, by = H / 10))

## A_omega ----
expo_limits = case_when(
  type == "expo" & scale_c <= 1.5 ~ c(-10, 10),
  .default = NULL)

ggplot(data = data.frame(var = "add", val = add_samp), aes(x = val)) +
  geom_freqpoly() +
  geom_vline(xintercept = 0, lwd = 0.5, col = "gray") +
  geom_vline(xintercept = aomega, col = "red", lwd = 1) +
  scale_x_continuous(name = expression("Additionality (Mg CO"[2]*")")) +
  scale_y_continuous(name = "Counts") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, subfolder, file_pref, "_3_aomega_distribution.", file_type), width = 15, height = 15, unit = "cm")


## Credits time series ----
if(type == "expo") {
  range_credit = c(-5, 45)
  range_release = c(0, 5)
} else if(type == "expo_portfolio") {
  range_credit = c(-5, 110)
  range_release = c(0, 30)
} else {
  range_credit = c(min(summ_credit$p05) - 1000, max(summ_credit$p95) + 1000)
  range_release = c(0, max(summ_release$p95) + 1000)
}

if(print_real) {
  ggplot(sim_credit_long) +
    geom_line(aes(x = t, y = val, color = rep), size = 0.5, show.legend = F) +
    geom_hline(yintercept = 0, size = 1, color = "black") +
    {if(type != "expo" & type != "expo_portfolio") geom_vline(xintercept = year_pres_obs - t0 + 1, size = 0.5, color = "black", lty = "dotted")} +
    scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) +
    scale_y_continuous(name = expression("Credits (Mg CO"[2]*")"), limits = range_credit) +
    theme_bw() +
    theme(axis.line = element_line(linewidth = 0.5),
          panel.grid = element_blank(),
      	axis.title = element_text(size = 24),
      	axis.text.x = element_text(size = 18, vjust = 0.5),
      	axis.text.y = element_text(size = 18),
      	plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
  ggsave(paste0(file_path, subfolder, file_pref, "_4a_sim_credit_real.", file_type), width = 15, height = 10, unit = "cm")
}

ggplot(summ_credit, aes(x = year)) +
  geom_ribbon(data = . %>% dplyr::select(year, p05, p95), aes(ymin = p05, ymax = p95), fill = "lightgray", alpha = 0.5) +
  geom_ribbon(data = . %>% dplyr::select(year, p25, p75), aes(ymin = p25, ymax = p75), fill = "darkgray", alpha = 0.5) +
  geom_line(data = . %>% dplyr::select(year, median), aes(y = median), size = 1, color = "black") +
  geom_hline(yintercept = 0, size = 0.5, color = "black") +
  {if(type != "expo" & type != "expo_portfolio") geom_vline(xintercept = year_pres_obs - t0 + 1, size = 0.5, color = "black", lty = "dotted")} +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) +
  scale_y_continuous(name = expression("Credits (Mg CO"[2]*")"), limits = range_credit) +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, subfolder, file_pref, "_4_sim_credit_rep.", file_type), width = 15, height = 10, unit = "cm")


## EP time series ----

if(print_real) {
  ggplot(as.data.frame(sim_ep) %>% replace(., . == 0, NA) %>% mutate(t = row_number()) %>% pivot_longer(V1:V100, names_to = "rep", values_to = "val")) +
    geom_line(aes(x = t, y = val, col = rep), size = 0.5, show.legend = F) +
    geom_hline(yintercept = 0, size = 0.5, color = "black") +
    geom_vline(xintercept = case_when(type != "expo" & type != "expo_portfolio" ~ year_pres_obs - t0 + 1, .default = NULL), size = 0.5, color = "black", lty = "dotted") +
    scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) +
    scale_y_continuous(name = "EP", limits = c(0, 1)) +
    theme_bw() +
    theme(axis.line = element_line(linewidth = 0.5),
          panel.grid = element_blank(),
      	axis.title = element_text(size = 24),
      	axis.text.x = element_text(size = 18, vjust = 0.5),
      	axis.text.y = element_text(size = 18),
      	plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
  ggsave(paste0(file_path, subfolder, file_pref, "_6a_sim_ep_real.", file_type), width = 15, height = 10, unit = "cm")
}

ggplot(summ_ep, aes(x = year)) +
  geom_ribbon(data = . %>% dplyr::select(year, p05, p95), aes(ymin = p05, ymax = p95), fill = "lightgray", alpha = 0.5) +
  geom_ribbon(data = . %>% dplyr::select(year, p25, p75), aes(ymin = p25, ymax = p75), fill = "darkgray", alpha = 0.5) +
  geom_line(data = . %>% dplyr::select(year, median), aes(y = median), size = 1, color = "black") +
  geom_vline(xintercept = case_when(type != "expo" & type != "expo_portfolio" ~ year_pres_obs - t0 + 1, .default = NULL), size = 0.5, color = "black", lty = "dotted") +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) +
  scale_y_continuous(name = "EP", limits = c(0, 1)) +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, subfolder, file_pref, "_6_sim_ep_rep.", file_type), width = 15, height = 10, unit = "cm")


## Failure risk time series ----
summ_risk = summ_cred %>%
  mutate(risk = 1 - cred)
ggplot(summ_risk, aes(x = year)) +
  geom_line(aes(y = risk), size = 1) +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) + 
  scale_y_continuous(name = "Failure risk", limits = c(0, 0.5)) + 
  geom_hline(yintercept = 0.05, color = "red", lty = "dashed", size = 0.5) +
  {if(type != "expo" & type != "expo_portfolio") geom_vline(xintercept = year_pres_obs - t0 + 1, size = 0.5, color = "black", lty = "dotted")} +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, subfolder, file_pref, "_7_sim_risk_rep.", file_type), width = 15, height = 10, unit = "cm")


# Sensitivity of simulation results to 1/lambda_c ----
lc_sensitivity = T
bp_sensitivity = F
ppr_sensitivity = F
H_sensitivity = F

file_type = "png"

if(lc_sensitivity) {
  type = "expo" #expo, project, expo_portfolio, portfolio
  scale_c_vec = c(seq(1, 10, by = 0.1))
  
  mean_credit = matrix(NA, n_rep, length(scale_c_vec))
  max_EP = matrix(NA, n_rep, length(scale_c_vec))
  mean_pact = matrix(NA, n_rep, length(scale_c_vec))
  risk = matrix(NA, n_rep, length(scale_c_vec))
  
  for(c_i in 1:length(scale_c_vec)){
    scale_c = scale_c_vec[c_i]
    source(paste0(file_path, "carbon_loss_simulations_core.R"))
    
    #relationship between 1 / lambda_c and mean of [median of credits] from all years
    mean_credit[, c_i] = apply(sim_credit, 2, sum) / (H - bp)
    
    #relationship between 1 / lambda_c and max EP from all years
    max_EP[, c_i] = apply(sim_ep, 2, max)
    
    #relationship between 1 / lambda_c and credibility (proportion of years without credibility incident)
    risk[, c_i] = apply(sim_ep[-c(1:bp), ], 2, function(x) length(which(x <= 0)) / (H - bp))
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
             ci_margin = qt(0.975, df = n_rep - 1) * sd / sqrt(n_rep),
             ci_low = mean - ci_margin,
             ci_high = mean + ci_margin) %>%
      select(-all_of(c(1:n)))
    return(output)
  }
  
  ## Figure 2
  data_credit = SummariseAcrossTests(mean_credit, scale_c_vec - 1)
  ggplot(data = data_credit, aes(x = var)) +
#    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "lightgray") +
    geom_line(aes(y = mean)) +
    scale_x_continuous(name = "Drawdown", limits = c(0, 5), breaks = seq(0, 5, by = 1)) + 
    scale_y_continuous(name = "Credits") + 
    theme_bw() +
    theme(axis.line = element_line(linewidth = 0.5),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 18, vjust = 0.5),
          axis.text.y = element_text(size = 18, angle = 45),
          plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
  ggsave(paste0(file_path, "figure2_drawdown_vs_credit.", file_type), width = 10, height = 10, unit = "cm")

  data_max_EP = SummariseAcrossTests(max_EP, scale_c_vec - 1)
  ggplot(data = data_max_EP, aes(x = var)) +
#    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "lightgray") +
    geom_line(aes(y = mean)) +
    scale_x_continuous(name = "Drawdown", limits = c(0, 5), breaks = seq(0, 5, by = 1)) + 
    scale_y_continuous(name = "EP") + 
    theme_bw() +
    theme(axis.line = element_line(linewidth = 0.5),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 18, vjust = 0.5),
          axis.text.y = element_text(size = 18, angle = 45),
          plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
  ggsave(paste0(file_path, "figure2_drawdown_vs_EP.", file_type), width = 10, height = 10, unit = "cm")
  
  data_risk = SummariseAcrossTests(risk, scale_c_vec - 1)
  ggplot(data = data_risk, aes(x = var)) +
#    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "lightgray") +
    geom_line(aes(y = mean)) +
    geom_hline(yintercept = 0.05, col = "red", lty = "dashed", size = 0.5) +
    scale_x_continuous(name = "Drawdown", limits = c(0, 5), breaks = seq(0, 5, by = 1)) + 
    scale_y_continuous(name = "Failure risk", limits = c(0, 0.1)) + 
    theme_bw() +
    theme(axis.line = element_line(linewidth = 0.5),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 18, vjust = 0.5),
          axis.text.y = element_text(size = 18, angle = 45),
          plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
  ggsave(paste0(file_path, "figure2_drawdown_vs_risk.", file_type), width = 10, height = 10, unit = "cm")
}


# Sensitivity to warm-up period length ----
type = "expo" #expo, project, expo_portfolio, portfolio
lc_sensitivity = F
bp_sensitivity = T
ppr_sensitivity = F
H_sensitivity = F
n_rep = 1000

a = Sys.time()

if(bp_sensitivity) {
  scale_c_vec = c(1.1, 1.3, 1.5, 1.7)
  bp_vec = seq(1, 20, by = 1)
  list_cred = vector("list", length(scale_c_vec))
  list_summ_cred = vector("list", length(scale_c_vec))
  
  for(c_i in 1:length(scale_c_vec)){
    scale_c = scale_c_vec[c_i]
    
#    plot_bp = ifelse(scale_c %in% seq(1.1, 1.7, by = 1), T, F)
    
    mean_credit = matrix(NA, n_rep, length(bp_vec))
    max_EP = matrix(NA, n_rep, length(bp_vec))
    final_EP = matrix(NA, n_rep, length(bp_vec))
    cred = matrix(NA, n_rep, length(bp_vec))
    
    for(bp_i in 1:length(bp_vec)){
      bp = bp_vec[bp_i]
      source(paste0(file_path, "carbon_loss_simulations_core.R"))
      
      #relationship between 1 / lambda_c and mean of [median of credits] from all years
      mean_credit[, bp_i] = apply(sim_credit, 2, sum) / (H - bp)
      
      #relationship between 1 / lambda_c and max EP from all years
      max_EP[, bp_i] = apply(sim_ep, 2, max)
      
      #relationship between 1 / lambda_c and credibility (proportion of years without credibility incident)
      cred[, bp_i] = apply(sim_ep[-c(1:bp), ], 2, function(x) length(which(x > 0)) / (H - bp))
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
               ci_margin = qt(0.975, df = n_rep - 1) * sd / sqrt(n_rep),
               ci_low = mean - ci_margin,
               ci_high = mean + ci_margin) %>%
        select(-all_of(c(1:n)))
      return(output)
    }
    
    ## Figures
    data_credit = SummariseAcrossTests(mean_credit, bp_vec)
    data_max_EP = SummariseAcrossTests(max_EP, bp_vec)
    
    list_cred[[c_i]] = cred
    list_summ_cred[[c_i]] = SummariseAcrossTests(cred, bp_vec) %>%
      mutate(scale_c = scale_c)
  }
}
b = Sys.time()
b - a

summ_risk = do.call(rbind, list_summ_cred) %>%
  mutate(ci_low = 1 - ci_low, ci_high = 1 - ci_high, mean = 1 - mean, scale_c = as.factor(scale_c))


#bp length vs failure risk, each curve is a scale_c
ggplot(data = summ_risk, aes(x = var, group = scale_c)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = scale_c), alpha = 0.3) +
  geom_line(aes(y = mean, color = scale_c), size = 0.5) +
  geom_hline(yintercept = 0.05, color = "red", lty = "dashed", size = 1) +
  scale_x_continuous(name = "Warm-up period (years)") + 
  scale_y_continuous(name = "Failure risk", limits = c(0, 0.2)) + 
  scale_color_manual(values = c("red", "pink", "lightblue", "blue")) +
  scale_fill_manual(values = c("red", "pink",  "lightblue", "blue"), guide = NULL) + 
  guides(color = guide_legend(title = "Drawdown rate", title.theme = element_text(size = 18),
                              label.theme = element_text(size = 16), override.aes = list(linewidth = 1))) +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, subfolder, "sensit_wupl_vs_risk_by_scale_c.", file_type), width = 20, height = 15, unit = "cm")


# Sensitivity to post-project release level ----
type = "expo" #expo, project, expo_portfolio, portfolio
lc_sensitivity = F
bp_sensitivity = F
ppr_sensitivity = T
H_sensitivity = F
n_rep = 100

if(ppr_sensitivity) {
  scale_c = 1.3 #5, 10: "good" project, 1.3: "bad" project
  bp = 5
  ppr_vec = seq(1, 5, by = 0.1)

  mean_credit = matrix(NA, n_rep, length(ppr_vec))
  max_EP = matrix(NA, n_rep, length(ppr_vec))
  final_EP = matrix(NA, n_rep, length(ppr_vec))
  mean_cred = matrix(NA, n_rep, length(ppr_vec))
  
  for(ppr_i in 1:length(ppr_vec)){
    a = Sys.time()
    rate_postproj = ppr_vec[ppr_i]
    
    source(paste0(file_path, "carbon_loss_simulations_core.R"))
    
    #relationship between 1 / lambda_c and mean of [median of credits] from all years
    #mean_credit[, ppr_i] = apply(sim_credit, 2, sum) / (H - bp)
    
    #relationship between 1 / lambda_c and max EP from all years
    max_EP[, ppr_i] = apply(sim_ep, 2, max)
    
    #relationship between 1 / lambda_c and credibility (proportion of years without credibility incident)
    #mean_cred[, ppr_i] = apply(sim_ep[-c(1:bp), ], 2, function(x) length(which(x > 0)) / (H - bp))
    b = Sys.time()
    b - a
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
             ci_margin = qt(0.975, df = n_rep - 1) * sd / sqrt(n_rep),
             ci_low = mean - ci_margin,
             ci_high = mean + ci_margin) %>%
      select(-all_of(c(1:n)))
    return(output)
  }
  
  ## Figures
  data_max_EP = SummariseAcrossTests(max_EP, ppr_vec * scale_c)

  ggplot(data = data_max_EP, aes(x = var)) +
#    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "lightgray") +
    geom_line(aes(y = mean), size = 1) +
    scale_y_continuous(limits = c(0, 0.75)) + 
    xlab(label = "Post-project release rate") + 
    ylab(label = "EP") + 
    theme_bw() +
    theme(axis.line = element_line(linewidth = 0.5),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 22),
          axis.text.x = element_text(size = 14, vjust = 0.5),
          axis.text.y = element_text(size = 14),
          plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
  ggsave(paste0(file_path, subfolder, "sensit_ppr_vs_EP.", file_type), width = 15, height = 10, unit = "cm")
}


# Sensitivity to project length ----
type = "expo" #expo, project, expo_portfolio, portfolio
lc_sensitivity = F
bp_sensitivity = F
ppr_sensitivity = F
H_sensitivity = T
n_rep = 100

if(H_sensitivity) {
  scale_c = 1.3 #5, 10: "good" project, 1.3: "bad" project
  bp = 5
  rate_postproj = 2
  H_vec = seq(10, 100, by = 5)
  
  mean_credit = matrix(NA, n_rep, length(H_vec))
  max_EP = matrix(NA, n_rep, length(H_vec))
  final_EP = matrix(NA, n_rep, length(H_vec))
  mean_cred = matrix(NA, n_rep, length(H_vec))
  
  for(H_i in 1:length(H_vec)){
    H = H_vec[H_i]
    
    source(paste0(file_path, "carbon_loss_simulations_core.R"))
    
    #relationship between 1 / lambda_c and mean of [median of credits] from all years
    #mean_credit[, ppr_i] = apply(sim_credit, 2, sum) / (H - bp)
    
    #relationship between 1 / lambda_c and max EP from all years
    max_EP[, H_i] = apply(sim_ep, 2, max)
    
    #relationship between 1 / lambda_c and credibility (proportion of years without credibility incident)
    #mean_cred[, ppr_i] = apply(sim_ep[-c(1:bp), ], 2, function(x) length(which(x > 0)) / (H - bp))
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
             ci_margin = qt(0.975, df = n_rep - 1) * sd / sqrt(n_rep),
             ci_low = mean - ci_margin,
             ci_high = mean + ci_margin) %>%
      select(-all_of(c(1:n)))
    return(output)
  }
  
  ## Figures
  data_max_EP = SummariseAcrossTests(max_EP, H_vec)
  
  ggplot(data = data_max_EP, aes(x = var)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "lightgray") +
    geom_line(aes(y = mean), size = 1) +
    scale_y_continuous(limits = c(0, 0.75)) + 
    xlab(label = "Project duration") + 
    ylab(label = "EP") + 
    theme_bw() +
    theme(axis.line = element_line(linewidth = 0.5),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 22),
          axis.text.x = element_text(size = 14, vjust = 0.5),
          axis.text.y = element_text(size = 14),
          plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
  ggsave(paste0(file_path, subfolder, "sensit_proj_len_vs_EP.", file_type), width = 15, height = 10, unit = "cm")
}