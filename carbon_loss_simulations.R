# Initialise ----

#clean workspace
rm(list = ls())

#load packages
library(tidyverse)
library(magrittr)
library(mclust)

#load functions
source("functions.R") #wrapper functions
source("SimulatePermanence.r") #core simulation function
source("SaveStandard.r") #save standard output to rds objects
source("PlotStandard.r") #plot summary of time series for credit, EP and reversal risk

#load social cost of carbon data
scc = read.csv("scc_extended.csv", row.names = 1)
scc = scc %>%
  mutate(value = central) %>%
  dplyr::select(year, value)

# The core function SimulatePermanence() takes the following arguments:
# 1. type: character, "theo" for theoretical projects, "real" for real-life projects

# 2. mean_drawdown: numerical vector, mean drawdown rate(s) for theoretical project(s)
# # Its inverse will be used as the lambda parameter of an exponential distribution
# # When aggregate_type is NULL and the length of this argument > 1, it initialises a custom theoretical aggregated project

# 3. sites: character vector, name(s) for the real-life project(s)
# # The function will look for a data frame of carbon flux time series at the path "/project_input_data/sites.csv"
# # The data frame must contain the following columns:
# ## year: numeric
# ## var: character, either "project", "counterfactual" or "additionality"
# ## val: numeric, total carbon flux from the last year to this year (Mg CO2e)
# ## n_sim: numeric, index of repetition
# ## started: boolean, whether the year is larger than project start (t0)
# # When aggregate_type is NULL and the length of this argument > 1, it initialises a custom real-life aggregated project

# 4. aggregate_type: character, specifying default settings for theoretical/real-life aggregated projects
# # This argument overrides both mean_drawdown and sites when not NULL:
# # "A": theoretical, mean_drawdown = c(1.1, 1.1, 1.1, 5)
# # "B": theoretical, mean_drawdown = c(1.1, 1.1, 5, 5)
# # "C": theoretical, mean_drawdown = c(1.1, 5, 5, 5)
# # "three": real-life, sites = c("Gola_country", "CIF_Alto_Mayo", "VCS_1396")
# # "four": real-life, sites = c("Gola_country", "CIF_Alto_Mayo", "VCS_1396", "VCS_934")

# 5. verbose: boolean, whether to print basic output at each timestep (default: FALSE)
# 6. runtime: boolean, whether to print (default: TRUE)
# 7. omega: numeric, threshold of acceptable reversal risk (default: 0.05)
# 8. n_rep: numeric, number of repetitions (default: 100)
# 9. H: numeric, project duration (years) (default: 50)
# 10. D: numeric, discount rate (default: 0.03)
# 11. warmup: numeric, warm-up period (years) (default: 5)
# 12. postproject_ratio: numeric, ratio of post-project release compared to during-project additionality accumulation rate (default: 2)
# 13. scc_df: numeric data frame, containing two columns, year and value (default: scc, a data frame which should be loaded and prepared)


# Run simulations ----

#set output path
out_path = "C:/Users/epr26/OneDrive - University of Cambridge/assess_permanence_out/"

#run simulations: one of the four following options
outlist = SimulatePermanence(type = "theo", mean_drawdown = 1.1) #single theoretical projects: 1.1, 2, 5
outlist = SimulatePermanence(type = "real", sites = "VCS_934") #single real-life projects: Gola_country (Gola), CIF_Alto_Mayo (Alto Mayo), VCS_1396 (RPA), VCS_934 (Mai Ndombe)
outlist = SimulatePermanence(type = "theo", aggregate_type = "B") #aggregated theoretical projects: A, B, C
outlist = SimulatePermanence(type = "real", aggregate_type = "four") #aggregated real-life projects: three, four

#save output as RDS objects
summary_out = SaveStandard(outlist, out_path)

#read output from saved RDS object
#summary_out = readRDS(file = paste0(out_path, type_label, "/", type_label, "_output.rds"))

#plotting variable: y-axis range for time series of credits
if(summary_out$type == "theo") {
    range_credit = c(-5, 45)
} else {
    range_credit = c(min(summary_out$credit$p05) * 0.95, max(summary_out$credit$p95) * 1.05)
}

#plot summary of time series for credit, EP and reversal risk
PlotStandard(summary_out, range_credit = range_credit, out_path = out_path, show_axis_title = T)


# Figure 4. Credits/EP/reversal risk vs. mean drawdown rate ----
mean_drawdown_vec = c(seq(1, 6, by = 0.1))
n_rep = 100
file_type = "png"

mean_credit = matrix(NA, n_rep, length(mean_drawdown_vec))
max_EP = matrix(NA, n_rep, length(mean_drawdown_vec))
risk = matrix(NA, n_rep, length(mean_drawdown_vec))

for(i in seq_along(mean_drawdown_vec)){
  mean_drawdown_i = mean_drawdown_vec[i]

  outlist = SimulatePermanence(type = "theo", mean_drawdown = mean_drawdown_i, n_rep = n_rep)
  H = outlist$common_var$H
  warmup = outlist$common_var$warmup

  #mean credits from all years except warm-up period
  mean_credit[, i] = apply(outlist$sim_credit, 2, sum) / (H - warmup)
  
  #max EP from all years
  max_EP[, i] = apply(outlist$sim_ep, 2, max)
  
  #per-repetition reversal risk (proportion of years without positive credit)
  risk[, i] = apply(outlist$sim_reversal, 2, sum) / (H - warmup)
}

summ_sensitivity_drawdown = list(mean_credit = mean_credit,
                                 max_EP = max_EP,
                                 risk = risk)
saveRDS(summ_sensitivity_drawdown, file = paste0(out_path, "output_fig4_sensitivity_drawdown.rds"))
# summ_sensit_dd_rate = readRDS(file = paste0(file_path, subfolder, "sensitivity_dd_rate_output.rds"))
# mean_credit = summ_sensit_dd_rate$mean_credit
# max_EP = summ_sensit_dd_rate$max_EP
# risk = summ_sensit_dd_rate$risk

actual_drawdown_vec = mean_drawdown_vec - 1
#Fig. 4a
data_credit = SummariseAcrossTests(input = mean_credit, var = actual_drawdown_vec, n_rep = n_rep) %>%
  filter(var <= 5)
ggplot(data = data_credit, aes(x = var, y = mean)) +
  geom_line() +
  scale_x_continuous(name = "Drawdown", limits = c(0, 5), breaks = seq(0, 5, by = 1)) + 
  scale_y_continuous(name = "Credits") + 
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
ggsave(paste0(out_path, "output_fig4a_credit_vs_drawdown.", file_type), width = 10, height = 10, unit = "cm")

#Fig. 4b
data_max_EP = SummariseAcrossTests(input = max_EP, var = actual_drawdown_vec, n_rep = n_rep) %>%
  filter(var <= 5)
ggplot(data = data_max_EP, aes(x = var, y = mean)) +
  geom_line() +
  scale_x_continuous(name = "Drawdown", limits = c(0, 5), breaks = seq(0, 5, by = 1)) + 
  scale_y_continuous(name = "EP") + 
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
ggsave(paste0(out_path, "output_fig4b_EP_vs_drawdown.", file_type), width = 10, height = 10, unit = "cm")

#Fig. 4c
data_risk = SummariseAcrossTests(input = risk, var = actual_drawdown_vec, n_rep = n_rep) %>%
  filter(var <= 5)
ggplot(data = data_risk, aes(x = var, y = mean)) +
  geom_line() +
  geom_hline(yintercept = 0.05, col = "red", lty = "dashed") +
  scale_x_continuous(name = "Drawdown", limits = c(0, 5), breaks = seq(0, 5, by = 1)) + 
  scale_y_continuous(name = "Reversal risk") + 
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
ggsave(paste0(out_path, "output_fig4c_risk_vs_drawdown.", file_type), width = 10, height = 10, unit = "cm")


# Figure 5. All time series of mean drawdown rate 1.1, 2, 5
out_path = "C:/Users/epr26/OneDrive - University of Cambridge/assess_permanence_out/"
out_path_complete = paste0(out_path, "theo/")
dir.create(out_path_complete)
range_credit = c(-5, 45)
mean_drawdown_vec = c(1.1, 2, 5)

for(x in mean_drawdown_vec) {
  outlist = SimulatePermanence(type = "theo", mean_drawdown = x)
  type_label = outlist$common_var$type_label
  summary_out = SaveStandard(outlist, out_path = out_path_complete, file_prefix = type_label)
  PlotStandard(summary_out, range_credit = range_credit, show_axis_title = T, out_path = out_path_complete, file_prefix = type_label)
}


# Figure 6. Sensitivity to warm-up period ----
mean_drawdown_vec = c(1.1, 1.3, 1.5, 1.7)
warmup_vec = 1:20
n_rep = 1000
file_type = "png"

summ_risk_list = vector("list", length(mean_drawdown_vec))

for(i in seq_along(mean_drawdown_vec)){
  mean_drawdown_i = mean_drawdown_vec[i]
  
  mean_credit = matrix(NA, n_rep, length(warmup_vec))
  max_EP = matrix(NA, n_rep, length(warmup_vec))
  final_EP = matrix(NA, n_rep, length(warmup_vec))
  risk = matrix(NA, n_rep, length(warmup_vec))
  
  for(j in seq_along(warmup_vec)){
    warmup_j = warmup_vec[j]

    cat("Warm-up period used:", warmup_j, "\n")
    outlist = SimulatePermanence(type = "theo", mean_drawdown = mean_drawdown_i, warmup = warmup_j, n_rep = n_rep)
    H = outlist$common_var$H

    #per-repetition reversal risk (proportion of years without positive credit)
    risk[, j] = apply(outlist$sim_reversal, 2, sum) / (H - warmup_j)
  }
  
  summ_risk_list[[i]] = SummariseAcrossTests(input = risk, var = warmup_vec, n_rep = n_rep) %>%
    mutate(actual_drawdown = mean_drawdown_i  - 1)
}

summ_risk = do.call(rbind, summ_risk_list) %>%
  mutate(actual_drawdown = as.factor(actual_drawdown))
saveRDS(summ_risk, file = paste0(out_path, "output_fig6_sensitivity_warmup.rds"))
#summ_risk = readRDS(file = paste0(out_path, "output_fig6_sensitivity_warmup.rds"))

#reversal risk vs warm-up period length, grouped by drawdown rate
ggplot(data = summ_risk, aes(x = var, group = actual_drawdown)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = actual_drawdown), alpha = 0.3) +
  geom_line(aes(y = mean, color = actual_drawdown), size = 0.5) +
  geom_hline(yintercept = 0.05, color = "red", lty = "dashed", size = 1) +
  scale_x_continuous(name = "Warm-up period (years)") + 
  scale_y_continuous(name = "Reversal risk", limits = c(0, 0.2)) + 
  scale_color_manual(values = c("red", "pink", "lightblue", "blue")) +
  scale_fill_manual(values = c("red", "pink",  "lightblue", "blue"), guide = NULL) + 
  guides(color = guide_legend(title = "Drawdown", title.theme = element_text(size = 18),
                              label.theme = element_text(size = 16), override.aes = list(linewidth = 1))) +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(1.5, 1, 0.7, 0.7, "cm"))
ggsave(paste0(out_path, "output_fig6_risk_vs_warmup_group_drawdown.", file_type), width = 20, height = 15, unit = "cm")


# Figure 7. Sensitivity to post-project release ratio ----
mean_drawdown = 1.3
ppr_ratio_vec = seq(1, 5, by = 0.1)
n_rep = 100
file_type = "png"

mean_credit = matrix(NA, n_rep, length(ppr_ratio_vec))
max_EP = matrix(NA, n_rep, length(ppr_ratio_vec))
final_EP = matrix(NA, n_rep, length(ppr_ratio_vec))
mean_cred = matrix(NA, n_rep, length(ppr_ratio_vec))

for(i in seq_along(ppr_ratio_vec)){
  ppr_ratio_i = ppr_ratio_vec[i]
  
  cat("Post-project release ratio:", ppr_ratio_i, "\n")
  outlist = SimulatePermanence(type = "theo", mean_drawdown = mean_drawdown, ppr_ratio = ppr_ratio_i, n_rep = n_rep)

  #max EP from all years
  max_EP[, i] = apply(outlist$sim_ep, 2, max)
}

data_max_EP = SummariseAcrossTests(input = max_EP, var = ppr_ratio_vec * (mean_drawdown - 1), n_rep = n_rep)
saveRDS(data_max_EP, file = paste0(out_path, "output_fig7_sensitivity_ppr.rds"))
#data_max_EP = readRDS(file = paste0(out_path, "output_fig7_sensitivity_ppr.rds"))

ggplot(data = data_max_EP, aes(x = var, y = mean)) +
  geom_line() +
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
ggsave(paste0(out_path, "output_fig7_risk_vs_ppr.", file_type), width = 15, height = 10, unit = "cm")


# Figure 8. Sensitivity to duration
H_vec = seq(10, 100, by = 5)
mean_drawdown = 1.3
n_rep = 100
file_type = "png"

mean_credit = matrix(NA, n_rep, length(H_vec))
max_EP = matrix(NA, n_rep, length(H_vec))
final_EP = matrix(NA, n_rep, length(H_vec))
mean_cred = matrix(NA, n_rep, length(H_vec))

for(i in seq_along(H_vec)){
  H_i = H_vec[i]

  cat("Project duration:", H_i, "\n")
  outlist = SimulatePermanence(type = "theo", mean_drawdown = mean_drawdown, H = H_i, n_rep = n_rep)

  #max EP from all years
  max_EP[, i] = apply(outlist$sim_ep, 2, max)
}

data_max_EP = SummariseAcrossTests(input = max_EP, var = H_vec, n_rep = n_rep)
saveRDS(data_max_EP, file = paste0(out_path, "output_fig8_sensitivity_H.rds"))
#data_max_EP = readRDS(file = paste0(out_path, "output_fig8_sensitivity_H.rds"))

ggplot(data = data_max_EP, aes(x = var, y = mean)) +
  geom_line() +
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
ggsave(paste0(out_path, "output_fig7_risk_vs_H.", file_type), width = 15, height = 10, unit = "cm")


# Figure 9. All time series of four real-life projects
out_path = "C:/Users/epr26/OneDrive - University of Cambridge/assess_permanence_out/"
out_path_complete = paste0(out_path, "real/")
dir.create(out_path_complete)
site_vec = c("Gola_country", "CIF_Alto_Mayo", "VCS_1396", "VCS_934")

for(x in site_vec) {
  outlist = SimulatePermanence(type = "real", sites = x)
  type_label = outlist$common_var$type_label
  summary_out = SaveStandard(outlist, out_path = out_path_complete, file_prefix = type_label)
  range_credit = c(min(summary_out$credit$p05) * 0.95, max(summary_out$credit$p95) * 1.05)
  PlotStandard(summary_out, range_credit = range_credit, show_axis_title = T, out_path = out_path_complete, file_prefix = type_label)
}

outlist = SimulatePermanence(type = "theo", aggregate_type = "B") #aggregated theoretical projects: A, B, C
outlist = SimulatePermanence(type = "real", aggregate_type = "four") #aggregated real-life projects: three, four


# Figure 10. All time series of aggregated theoretical
out_path = "C:/Users/epr26/OneDrive - University of Cambridge/assess_permanence_out/"
out_path_complete = paste0(out_path, "theo_aggr/")
dir.create(out_path_complete)
range_credit = c(-5, 45)
aggregate_type_vec = c("A", "B", "C")

for(x in aggregate_type_vec) {
  outlist = SimulatePermanence(type = "theo", aggregate_type = x)
  type_label = outlist$common_var$type_label
  summary_out = SaveStandard(outlist, out_path = out_path_complete, file_prefix = type_label)
  PlotStandard(summary_out, range_credit = range_credit, show_axis_title = T, out_path = out_path_complete, file_prefix = type_label)
}


# Figure 11. All time series of aggregated real-life
out_path = "C:/Users/epr26/OneDrive - University of Cambridge/assess_permanence_out/"
out_path_complete = paste0(out_path, "real_aggr/")
dir.create(out_path_complete)
aggregate_type_vec = c("three", "four")

for(x in aggregate_type_vec) {
  outlist = SimulatePermanence(type = "real", aggregate_type = x)
  type_label = outlist$common_var$type_label
  summary_out = SaveStandard(outlist, out_path = out_path_complete, file_prefix = type_label)
  range_credit = c(min(summary_out$credit$p05) * 0.95, max(summary_out$credit$p95) * 1.05)
  PlotStandard(summary_out, range_credit = range_credit, show_axis_title = T, out_path = out_path_complete, file_prefix = type_label)
}