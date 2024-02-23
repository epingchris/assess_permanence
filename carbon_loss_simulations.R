rm(list = ls())
library(tidyverse)
library(magrittr)
library(mclust)
file_path = "C:/Users/epr26/OneDrive - University of Cambridge/carbon_release_pattern/"

# Simulation settings ----

#set simulation type (hypothetical or real-life, single or aggregated)
type = "real_aggr" #hypo, real, hypo_aggr, real_aggr
dd_rate = 2  #1.1, 2, 5
hypo_sensit = "none" #none, dd_rate, warmup, ppr, H
hypo_aggr_type = "C" #A, B, C
real_aggr_type = "three" #five, four, three
project_site = "Gola_country" #Gola_country, WLT_VNCC_KNT, CIF_Alto_Mayo, VCS_1396, VCS_934
use_theo = T #T, F; type != hypo will override this and set use_theo as F
file_type = "png" #png, pdf
view_snapshot = F
print_allrep = F #print real trajectories of all 100 repetitions

#basic parameters
omega = 0.05
n_rep = 100 #number of repetitions
H = 50 #evaluation horizon; project duration
D = 0.03 #discount rate
warmup = 5 #warm-up period
postproject_ratio = 2 #ratio of post-project release compare to during project


# SCC ----
scc = read.csv(file = paste0(file_path, "scc_newpipeline.csv"))
year_range_scc = range(scc$year)
year_min_scc = 2000
year_max_scc = 2500
lm_scc = lm(log(central) ~ year, data = scc)
scc_before = data.frame(year = year_min_scc:(year_range_scc[1] - 1)) %>%
  mutate(central = exp(predict(lm_scc, newdata = ., se.fit = T)$fit)) %>%
  mutate(low = central * 0.5, high = central * 1.5) %>%
  relocate(c(1, 3, 2, 4))
scc_after = data.frame(year = (year_range_scc[2] + 1):year_max_scc) %>%
  mutate(central = exp(predict(lm_scc, newdata = ., se.fit = T)$fit)) %>%
  mutate(low = central * 0.5, high = central * 1.5) %>%
  relocate(c(1, 3, 2, 4))
scc_extended = rbind(scc_before, scc, scc_after)
write.csv(scc_extended, file = paste0(file_path, "scc_extended.csv"))

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

## Time series (summary) ----
if(type == "hypo") {
  range_credit = c(-5, 45)
} else if(type == "hypo_aggr") {
  range_credit = c(-5, 110)
} else {
  range_credit = c(min(summ_credit$p05) - 1000, max(summ_credit$p95) + 1000)
}

#1. Credits
ggplot(summ_credit, aes(x = year)) +
  geom_ribbon(data = . %>% dplyr::select(year, p05, p95), aes(ymin = p05, ymax = p95), fill = "lightgray", alpha = 0.5) +
  geom_ribbon(data = . %>% dplyr::select(year, p25, p75), aes(ymin = p25, ymax = p75), fill = "darkgray", alpha = 0.5) +
  geom_line(data = . %>% dplyr::select(year, median), aes(y = median), size = 1, color = "black") +
  geom_hline(yintercept = 0, size = 0.5, color = "black") +
  {if(type != "hypo" & type != "hypo_aggr") geom_vline(xintercept = year_expost - t0 + 1, size = 0.5, color = "black", lty = "dotted")} +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) +
  scale_y_continuous(name = expression("Credits (Mg CO"[2]*")"), limits = range_credit) +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, subfolder, file_pref, "_summ_credit.", file_type), width = 15, height = 10, unit = "cm")

#2. EP
ggplot(summ_ep, aes(x = year)) +
  geom_ribbon(data = . %>% dplyr::select(year, p05, p95), aes(ymin = p05, ymax = p95), fill = "lightgray", alpha = 0.5) +
  geom_ribbon(data = . %>% dplyr::select(year, p25, p75), aes(ymin = p25, ymax = p75), fill = "darkgray", alpha = 0.5) +
  geom_line(data = . %>% dplyr::select(year, median), aes(y = median), size = 1, color = "black") +
  {if(type != "hypo" & type != "hypo_aggr") geom_vline(xintercept = year_expost - t0 + 1, size = 0.5, color = "black", lty = "dotted")} +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) +
  scale_y_continuous(name = "EP", limits = c(0, 1)) +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, subfolder, file_pref, "_summ_ep.", file_type), width = 15, height = 10, unit = "cm")

# 3. Failure risk
ggplot(summ_risk, aes(x = year)) +
  geom_line(aes(y = risk), size = 1) +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) + 
  scale_y_continuous(name = "Failure risk", limits = c(0, 0.5)) + 
  geom_hline(yintercept = 0.05, color = "red", lty = "dashed", size = 0.5) +
  {if(type != "hypo" & type != "hypo_aggr") geom_vline(xintercept = year_expost - t0 + 1, size = 0.5, color = "black", lty = "dotted")} +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, subfolder, file_pref, "_summ_risk.", file_type), width = 15, height = 10, unit = "cm")


## Time series (all repetitions) ----
sim_credit_long = sim_credit %>%
  as.data.frame() %>%
  mutate(t = row_number()) %>%
  pivot_longer(V1:V100, names_to = "rep", values_to = "val")

sim_ep_long = sim_ep %>%
  as.data.frame() %>%
  replace(., . == 0, NA) %>%
  mutate(t = row_number()) %>%
  pivot_longer(V1:V100, names_to = "rep", values_to = "val")

#1. Credits
ggplot(sim_credit_long) +
  geom_line(aes(x = t, y = val, color = rep), size = 0.5, show.legend = F) +
  geom_hline(yintercept = 0, size = 1, color = "black") +
  {if(type != "hypo" & type != "hypo_aggr") geom_vline(xintercept = year_expost - t0 + 1, size = 0.5, color = "black", lty = "dotted")} +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) +
  scale_y_continuous(name = expression("Credits (Mg CO"[2]*")"), limits = range_credit) +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, subfolder, file_pref, "_allrep_credit.", file_type), width = 15, height = 10, unit = "cm")

#2. EP
ggplot(sim_ep_long) +
  geom_line(aes(x = t, y = val, col = rep), size = 0.5, show.legend = F) +
  geom_hline(yintercept = 0, size = 0.5, color = "black") +
  {if(type != "hypo" & type != "hypo_aggr") geom_vline(xintercept = year_expost - t0 + 1, size = 0.5, color = "black", lty = "dotted")} +
  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) +
  scale_y_continuous(name = "EP", limits = c(0, 1)) +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, subfolder, file_pref, "_allrep_ep.", file_type), width = 15, height = 10, unit = "cm")


## A_omega ----
expo_limits = case_when(
  type == "hypo" & dd_rate <= 1.5 ~ c(-10, 10),
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
ggsave(paste0(file_path, subfolder, file_pref, "_aomega.", file_type), width = 15, height = 15, unit = "cm")


# Sensitivity tests ----
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


## Drawdown rate (1/lambda_c) ----
type = "hypo"
hypo_sensit = "dd_rate" #none, dd_rate, warmup, ppr, H
dd_rate_vec = c(seq(1, 10, by = 0.1))

mean_credit = matrix(NA, n_rep, length(dd_rate_vec))
max_EP = matrix(NA, n_rep, length(dd_rate_vec))
mean_pact = matrix(NA, n_rep, length(dd_rate_vec))
risk = matrix(NA, n_rep, length(dd_rate_vec))

for(dd_i in 1:length(dd_rate_vec)){
  dd_rate = dd_rate_vec[dd_i]
  source(paste0(file_path, "carbon_loss_simulations_core.R"))
  
  #mean of credits from all years (except warm-up period)
  mean_credit[, dd_i] = apply(sim_credit, 2, sum) / (H - warmup)
  
  #max EP from all years
  max_EP[, dd_i] = apply(sim_ep, 2, max)
  
  #per-repetition failure risk (proportion of years without positive credit)
  risk[, dd_i] = apply(sim_failure, 2, sum) / (H - warmup)
}

## Figure 2
drawdown_vec = dd_rate_vec - 1
data_credit = SummariseAcrossTests(mean_credit, drawdown_vec)
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
ggsave(paste0(file_path, subfolder, "drawdown_vs_credit.", file_type), width = 10, height = 10, unit = "cm")

data_max_EP = SummariseAcrossTests(max_EP, drawdown_vec)
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
ggsave(paste0(file_path, subfolder, "drawdown_vs_EP.", file_type), width = 10, height = 10, unit = "cm")

data_risk = SummariseAcrossTests(risk, drawdown_vec)
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
ggsave(paste0(file_path, subfolder, "drawdown_vs_risk.", file_type), width = 10, height = 10, unit = "cm")


## Warm-up period length ----
type = "hypo"
hypo_sensit = "warmup" #none, dd_rate, warmup, ppr, H
n_rep = 1000

a = Sys.time()
dd_rate_vec = c(1.1, 1.3, 1.5, 1.7)
warmup_vec = seq(1, 20, by = 1)
list_summ_risk = vector("list", length(dd_rate_vec))

for(dd_i in 1:length(dd_rate_vec)){
  dd_rate = dd_rate_vec[dd_i]
  
  mean_credit = matrix(NA, n_rep, length(warmup_vec))
  max_EP = matrix(NA, n_rep, length(warmup_vec))
  final_EP = matrix(NA, n_rep, length(warmup_vec))
  risk = matrix(NA, n_rep, length(warmup_vec))
  
  for(w_i in 1:length(warmup_vec)){
    warmup = warmup_vec[w_i]
    source(paste0(file_path, "carbon_loss_simulations_core.R"))
    
    #mean of credits from all years (except warm-up period)
    #mean_credit[, w_i] = apply(sim_credit, 2, sum) / (H - warmup)
    
    #max EP from all years
    #max_EP[, w_i] = apply(sim_ep, 2, max)
    
    #per-repetition failure risk (proportion of years without positive credit)
    risk[, w_i] = apply(sim_failure, 2, sum) / (H - warmup)
  }
  
  ## Figures
  data_credit = SummariseAcrossTests(mean_credit, warmup_vec)
  data_max_EP = SummariseAcrossTests(max_EP, warmup_vec)
  
  list_summ_risk[[dd_i]] = SummariseAcrossTests(risk, warmup_vec) %>%
    mutate(dd_rate = dd_rate)
}
b = Sys.time()
b - a

summ_risk = do.call(rbind, list_summ_risk) %>%
  mutate(dd_rate = as.factor(dd_rate))


#warm-up period length vs failure risk, each curve corresponds to a drawdown rate
ggplot(data = summ_risk, aes(x = var, group = dd_rate)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = dd_rate), alpha = 0.3) +
  geom_line(aes(y = mean, color = dd_rate), size = 0.5) +
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
ggsave(paste0(file_path, subfolder, "warmup_vs_risk_grouped_drawdown.", file_type), width = 20, height = 15, unit = "cm")


## Post-project release rate ----
type = "hypo"
hypo_sensit = "ppr" #none, dd_rate, warmup, ppr, H
n_rep = 100
dd_rate = 1.3
warmup = 5
ppr_vec = seq(1, 5, by = 0.1)

mean_credit = matrix(NA, n_rep, length(ppr_vec))
max_EP = matrix(NA, n_rep, length(ppr_vec))
final_EP = matrix(NA, n_rep, length(ppr_vec))
mean_cred = matrix(NA, n_rep, length(ppr_vec))

for(ppr_i in 1:length(ppr_vec)){
  postproject_ratio = ppr_vec[ppr_i]
  
  source(paste0(file_path, "carbon_loss_simulations_core.R"))
  
  #mean of credits from all years (except the warm-up period)
  #mean_credit[, ppr_i] = apply(sim_credit, 2, sum) / (H - warmup)
  
  #max EP from all years
  max_EP[, ppr_i] = apply(sim_ep, 2, max)
  
  #per-repetition failure risk (proportion of years without positive credit)
  #risk[, ppr_i] = apply(sim_failure, 2, sum) / (H - warmup)
}

## Figures
data_max_EP = SummariseAcrossTests(max_EP, ppr_vec * dd_rate)

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
ggsave(paste0(file_path, subfolder, "ppr_vs_EP.", file_type), width = 15, height = 10, unit = "cm")


## Project duration ----
type = "hypo"
hypo_sensit = "H" #none, dd_rate, warmup, ppr, H
n_rep = 100
dd_rate = 1.3
warmup = 5
postproject_ratio = 2
H_vec = seq(10, 100, by = 5)

mean_credit = matrix(NA, n_rep, length(H_vec))
max_EP = matrix(NA, n_rep, length(H_vec))
final_EP = matrix(NA, n_rep, length(H_vec))
mean_cred = matrix(NA, n_rep, length(H_vec))

for(H_i in 1:length(H_vec)){
  H = H_vec[H_i]
  
  source(paste0(file_path, "carbon_loss_simulations_core.R"))
  
  #mean of credits from all years (except the warm-up period)
  #mean_credit[, ppr_i] = apply(sim_credit, 2, sum) / (H - warmup)
  
  #max EP from all years
  max_EP[, H_i] = apply(sim_ep, 2, max)
  
  #per-repetition failure risk (proportion of years without positive credit)
  #risk[, ppr_i] = apply(sim_failure, 2, sum) / (H - warmup)
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
ggsave(paste0(file_path, subfolder, "H_vs_EP.", file_type), width = 15, height = 10, unit = "cm")