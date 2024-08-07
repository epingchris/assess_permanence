rm(list = ls())
library(tidyverse)
library(magrittr)
library(mclust)

# Settings ----

#read data and functions
source("functions.R")
source("setInput.R") #to be replaced by SimulatePermanence
source("SimulatePermanence.r") #core simulation function
scc = read.csv("scc_extended.csv", row.names = 1) #social cost of carbon
scc = scc %>%
  mutate(value = central) %>%
  dplyr::select(year, value)

#The function SimulatePermanence() takes the following arguments:
# type: character, "theo" for theoretical projects, "real" for real-life projects
#
# mean_drawdown: numerical vector, mean drawdown rate(s) for theoretical project(s)
## its inverse will be used as the lambda parameter of an exponential distribution)
## when aggregate_type is NULL and the length of this argument > 1, it initialises a custom theoretical aggregated project
#
# sites: character vector, name(s) for the real-life project(s)
## the function will look for AGB time series data in "/project_input_data/sites.csv"
## when aggregate_type is NULL and the length of this argument > 1, it initialises a custom real-life aggregated project
#
# aggregate_type: character, specifying default settings for theoretical/real-life aggregated projects
## this argument overrides both mean_drawdown and sites when not NULL:
## "A": theoretical, mean_drawdown = c(1.1, 1.1, 1.1, 5)
## "B": theoretical, mean_drawdown = c(1.1, 1.1, 5, 5)
## "C": theoretical, mean_drawdown = c(1.1, 5, 5, 5)
## "three": real-life, sites = c("Gola_country", "CIF_Alto_Mayo", "VCS_1396")
## "four": real-life, sites = c("Gola_country", "CIF_Alto_Mayo", "VCS_1396", "VCS_934")
#
# verbose: boolean, whether to print basic output at each timestep (default: FALSE)
# runtime: boolean, whether to print (default: TRUE)
# omega: numeric, threshold of acceptable reversal risk (default: 0.05)
# n_rep: numeric, number of repetitions (default: 100)
# H: numeric, project duration (years) (default: 50)
# D: numeric, discount rate (default: 0.03)
# warmup: numeric, warm-up period (years) (default: 5)
# postproject_ratio: numeric, ratio of post-project release compared to during-project additionality accumulation rate (default: 2)
# scc: numeric data frame, containing two columns, year and value (default: scc, a data frame which should be loaded and prepared)

out_path = "C:/Users/epr26/OneDrive - University of Cambridge/assess_permanence_out/"
file_type = "png" #png, pdf
hypo_sensit = "none" #none, dd_rate, warmup, ppr, H


# Run simulation ----
outlist = SimulatePermanence(type = "theo", mean_drawdown = 2) #1.1, 2, 5
outlist = SimulatePermanence(type = "real", sites = "CIF_Alto_Mayo") #Gola_country, CIF_Alto_Mayo, VCS_1396, VCS_934
outlist = SimulatePermanence(type = "theo", aggregate_type = "B") #A, B, C
outlist = SimulatePermanence(type = "real", aggregate_type = "four") #three, four


# Save output objects to RDS objects ----

#unpack output objects
sim_p_loss = outlist$sim_p_loss
sim_c_loss = outlist$sim_c_loss
sim_additionality = outlist$sim_additionality
sim_credit = outlist$sim_credit
sim_benefit = outlist$sim_benefit
sim_aomega = outlist$sim_aomega
sim_release = outlist$sim_release
sim_damage = outlist$sim_damage
sim_ep = outlist$sim_ep
sim_pact = outlist$sim_pact
sim_reversal = outlist$sim_reversal
sim_buffer = outlist$sim_buffer
sim_schedule = outlist$sim_schedule
common_var = outlist$common_var

#summarise time series results
summ_additionality = SummariseSim(sim_additionality)
summ_credit = SummariseSim(sim_credit)
#summ_pact = SummariseSim(sim_pact)
summ_release = SummariseSim(sim_release[1:common_var$H, ])
#summ_buffer = SummariseSim(sim_buffer)
summ_aomega = SummariseSim(sim_aomega)

summ_ep = sim_ep %>%
  replace(., . == 0, NA) %>%
  SummariseSim() %>%
  replace(., . == Inf| . == -Inf, NA)

#per-year reversal risk (proportion of repetitions without positive credits)
summ_risk = data.frame(year = 1:common_var$H, risk = apply(sim_reversal, 1, sum) / ncol(sim_reversal))
summ_risk$risk[1:common_var$warmup] = NA

#gather output objects and save
summ_simulation = list(additionality = summ_additionality,
                       credit = summ_credit,
                       release = summ_release,
                       aomega = summ_aomega,
                       ep = summ_ep,
                       risk = summ_risk)
summ_complete = c(common_var, summ_simulation)

type_label = common_var$type_label
dir.create(paste0(out_path, type_label, "/"))
saveRDS(summ_complete, file = paste0(out_path, type_label, "/", type_label, "_output.rds"))




inpar = setInput(type = "theo", mean_drawdown = 1.1) #1.1, 2, 5
inpar = setInput(type = "real", sites = "Gola_country") #Gola_country, CIF_Alto_Mayo, VCS_1396, VCS_934
inpar = setInput(type = "theo", aggregate_type = "A") #A, B, C
inpar = setInput(type = "real", aggregate_type = "three") #three, four

source("carbon_loss_simulations_core.R")




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
  scale_x_continuous(name = "", breaks = yr_label, labels = yr_label) +
  scale_y_continuous(name = "", limits = range_credit) +
#  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) +
#  scale_y_continuous(name = expression("Credits (Mg CO"[2]*")"), limits = range_credit) +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 22),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
ggsave(paste0(file_path, subfolder, file_pref, "_summ_credit.", file_type), width = 12, height = 9, unit = "cm")

#2. Anticipated releases
summ = readRDS(file = "C:/Users/epr26/OneDrive - University of Cambridge/carbon_release_pattern/real/Gola_test_output.RDS")
ggplot(summ$release, aes(x = year)) +
  geom_ribbon(data = . %>% dplyr::select(year, p05, p95), aes(ymin = p05, ymax = p95), fill = "lightgray", alpha = 0.5) +
  geom_ribbon(data = . %>% dplyr::select(year, p25, p75), aes(ymin = p25, ymax = p75), fill = "darkgray", alpha = 0.5) +
  geom_line(data = . %>% dplyr::select(year, median), aes(y = median), size = 1, color = "black") +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
  {if(type %in% c("hypo", "hypo_aggr") == F) geom_vline(xintercept = year_expost - t0 + 1, linewidth = 0.5, color = "black", lty = "dotted")} +
  scale_x_continuous(name = "", breaks = yr_label, labels = yr_label) +
  scale_y_continuous(name = "", limits = range_credit) +
  #  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) +
  #  scale_y_continuous(name = expression("Credits (Mg CO"[2]*")"), limits = range_credit) +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 22),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
ggsave(paste0(file_path, subfolder, file_pref, "_summ_release.", file_type), width = 12, height = 9, unit = "cm")


#2. EP
ggplot(summ_ep, aes(x = year)) +
  geom_ribbon(data = . %>% dplyr::select(year, p05, p95), aes(ymin = p05, ymax = p95), fill = "lightgray", alpha = 0.5) +
  geom_ribbon(data = . %>% dplyr::select(year, p25, p75), aes(ymin = p25, ymax = p75), fill = "darkgray", alpha = 0.5) +
  geom_line(data = . %>% dplyr::select(year, median), aes(y = median), size = 1, color = "black") +
  {if(type != "hypo" & type != "hypo_aggr") geom_vline(xintercept = year_expost - t0 + 1, size = 0.5, color = "black", lty = "dotted")} +
  scale_x_continuous(name = "", breaks = yr_label, labels = yr_label) +
  scale_y_continuous(name = "", limits = c(0, 1)) +
#  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) +
#  scale_y_continuous(name = "EP", limits = c(0, 1)) +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 22),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
ggsave(paste0(file_path, subfolder, file_pref, "_summ_ep.", file_type), width = 12, height = 9, unit = "cm")

# 3. Non-delivery risk
ggplot(summ_risk, aes(x = year)) +
  geom_line(aes(y = risk), size = 1) +
  geom_hline(yintercept = 0.05, color = "red", lty = "dashed", size = 0.5) +
  {if(type != "hypo" & type != "hypo_aggr") geom_vline(xintercept = year_expost - t0 + 1, size = 0.5, color = "black", lty = "dotted")} +
  scale_x_continuous(name = "", breaks = yr_label, labels = yr_label) +
  scale_y_continuous(name = "", limits = c(0, 0.5)) +
#  scale_x_continuous(name = "Year", breaks = yr_label, labels = yr_label) + 
#  scale_y_continuous(name = "Non-delivery risk", limits = c(0, 0.5)) + 
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 22),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
ggsave(paste0(file_path, subfolder, file_pref, "_summ_risk.", file_type), width = 12, height = 9, unit = "cm")


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
  source("carbon_loss_simulations_core.R")
  
  #mean of credits from all years (except warm-up period)
  mean_credit[, dd_i] = apply(sim_credit, 2, sum) / (H - warmup)
  
  #max EP from all years
  max_EP[, dd_i] = apply(sim_ep, 2, max)
  
  #per-repetition non-delivery risk (proportion of years without positive credit)
  risk[, dd_i] = apply(sim_failure, 2, sum) / (H - warmup)
}

summ_sensit_dd_rate = list(mean_credit = mean_credit,
                           max_EP = max_EP,
                           risk = risk)
saveRDS(summ_sensit_dd_rate, file = paste0(file_path, subfolder, "sensitivity_dd_rate_output.RDS"))
summ_sensit_dd_rate = readRDS(file = paste0(file_path, subfolder, "sensitivity_dd_rate_output.RDS"))
mean_credit = summ_sensit_dd_rate$mean_credit
max_EP = summ_sensit_dd_rate$max_EP
risk = summ_sensit_dd_rate$risk

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
  scale_y_continuous(name = "Non-delivery risk", limits = c(0, 0.1)) + 
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
    source("carbon_loss_simulations_core.R")
    
    #mean of credits from all years (except warm-up period)
    #mean_credit[, w_i] = apply(sim_credit, 2, sum) / (H - warmup)
    
    #max EP from all years
    #max_EP[, w_i] = apply(sim_ep, 2, max)
    
    #per-repetition non-delivery risk (proportion of years without positive credit)
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
  mutate(dd_rate = as.factor(dd_rate - 1))

saveRDS(summ_risk, file = paste0(file_path, subfolder, "sensitivity_warmup_output.RDS"))
summ_risk = readRDS(file = paste0(file_path, subfolder, "sensitivity_warmup_output.RDS"))

#warm-up period length vs non-delivery risk, each curve corresponds to a drawdown rate
ggplot(data = summ_risk, aes(x = var, group = dd_rate)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = dd_rate), alpha = 0.3) +
  geom_line(aes(y = mean, color = dd_rate), size = 0.5) +
  geom_hline(yintercept = 0.05, color = "red", lty = "dashed", size = 1) +
  scale_x_continuous(name = "Warm-up period (years)") + 
  scale_y_continuous(name = "Non-delivery risk", limits = c(0, 0.2)) + 
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
  
  source("carbon_loss_simulations_core.R")
  
  #mean of credits from all years (except the warm-up period)
  #mean_credit[, ppr_i] = apply(sim_credit, 2, sum) / (H - warmup)
  
  #max EP from all years
  max_EP[, ppr_i] = apply(sim_ep, 2, max)
  
  #per-repetition non-delivery risk (proportion of years without positive credit)
  #risk[, ppr_i] = apply(sim_failure, 2, sum) / (H - warmup)
}

data_max_EP = SummariseAcrossTests(max_EP, ppr_vec * dd_rate)
saveRDS(data_max_EP, file = paste0(file_path, subfolder, "sensitivity_ppr_output.RDS"))
data_max_EP = readRDS(file = "C:/Users/epr26/OneDrive - University of Cambridge/carbon_release_pattern/sensit_ppr/sensitivity_ppr_output.RDS")

ggplot(data = data_max_EP, aes(x = var)) +
  #    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "lightgray") +
  geom_line(aes(y = mean), size = 1) +
  scale_y_continuous(limits = c(0, 0.5)) + 
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
  
  source("carbon_loss_simulations_core.R")
  
  #mean of credits from all years (except the warm-up period)
  #mean_credit[, ppr_i] = apply(sim_credit, 2, sum) / (H - warmup)
  
  #max EP from all years
  max_EP[, H_i] = apply(sim_ep, 2, max)
  
  #per-repetition non-delivery risk (proportion of years without positive credit)
  #risk[, ppr_i] = apply(sim_failure, 2, sum) / (H - warmup)
}

data_max_EP = SummariseAcrossTests(max_EP, H_vec)
saveRDS(data_max_EP, file = paste0(file_path, subfolder, "sensitivity_H_output.RDS"))
data_max_EP = readRDS(file = "C:/Users/epr26/OneDrive - University of Cambridge/carbon_release_pattern/sensit_H/sensitivity_H_output.RDS")

ggplot(data = data_max_EP, aes(x = var)) +
#  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "lightgray") +
  geom_line(aes(y = mean), size = 1) +
  scale_y_continuous(limits = c(0, 0.5)) + 
  xlab(label = "Duration (year)") + 
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
