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
file_path = "C:/Users/E-Ping Rau/OneDrive - University of Cambridge/carbon_release_pattern/"
source("/Users/E-Ping Rau/OneDrive - University of Cambridge/carbon_release_pattern/divergence.R")
# 1. Initiatlise ----

## 1a. Generate data for new sites ----
#rmarkdown::render("/Users/E-Ping Rau/OneDrive - University of Cambridge/4C_evaluations/R/Reports/eval_template/evaluation_epingrau.Rmd", clean = FALSE)
#eval_classes = c(1, 2, 3, 4) #only evaluate classes 1-4
## 1b. Load existing data for sites ----
site = "Gola_country" #Gola_country, WLT_VNCC_KNT, CIF_Alto_Mayo, VCS_1396, VCS_934
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


#plot time series
PlotTimeSeries = function(dat, y_label, y_lim, file_suffix){
  ggplot(data = dat, aes(col = var)) +
    geom_line(aes(x = year, y = val, lwd = series, alpha = series), show.legend = F) +
    scale_linewidth_manual("", values = c(1, rep(0.5, 20)), labels = c("Mean", "Simulated", rep(NULL, 19))) +
    scale_alpha_manual("", values = c(1, rep(0.1, 20)), labels = c("Mean", "Simulated", rep(NULL, 19))) +
    ggnewscale::new_scale("lwd") +
    ggnewscale::new_scale("alpha") +
    geom_line(data = subset(dat, series %in% c("mean", "V1")), aes(x = year, y = val, lwd = series, alpha = series)) +
    geom_vline(aes(xintercept = t0)) + 
    scale_color_manual("", values = c("red", "black"), labels = c("Counterfactual", "Project")) +
    scale_linewidth_manual("", values = c(1, 0.5), labels = c("Mean", "Simulated")) +
    scale_alpha_manual("", values = c(1, 0.1), labels = c("Mean", "Simulated")) +
    scale_x_continuous(name = "Year", breaks = seq(1990, 2021, by = 5)) + 
    scale_y_continuous(name = y_label, limits = y_lim) + 
    ggtitle("") +
    theme_bw() +
    theme(axis.line = element_line(linewidth = 0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "none",
          axis.title = element_text(size = 24),
          axis.text = element_text(size = 16))
  ggsave(paste0(file_path, site, file_suffix), width = 15, height = 20, unit = "cm")
}

y_min = switch(site,
               Gola = -4e+05,
               WLT_VNCC_KNT = -4e+05,
               CIF_Alto_Mayo = -8e+05,
               VCS_1396 = -8e+05,
               VCS_934 = -4e+06)

PlotTimeSeries(dat = subset(summ_flux, var %in% c("treatment_proj", "control_proj")),
               #lim_to_use = c(floor(min(dat$val) / (10 ^ 5)) * (10 ^ 5), 0)
               y_label = "C flux (Mg CO2e)",
               y_lim = c(y_min, 25000),
               file_suffix = "_2_time_series.png")


# 3. Describe empirical carbon loss distributions ----

absloss_p = subset(summ_flux, var == "treatment_proj" & year >= t0 & series != "mean") %>%
  mutate(val = val * (-1), var = NULL)
absloss_c = subset(summ_flux, var == "control_proj" & year >= t0 & series != "mean") %>%
  mutate(val = val * (-1), var = NULL)

absloss_p_bfr = subset(summ_flux, var == "treatment_proj" & year < t0 & series != "mean") %>%
  mutate(val = val * (-1), var = NULL)
absloss_c_bfr = subset(summ_flux, var == "control_proj" & year < t0 & series != "mean") %>%
  mutate(val = val * (-1), var = NULL)


## Epsilon (probability of carbon loss) ----
nrow(subset(absloss_p, val > 0)) / nrow(absloss_p)
nrow(subset(absloss_c, val > 0)) / nrow(absloss_c)

nrow(subset(absloss_p_bfr, val > 0)) / nrow(absloss_p_bfr)
nrow(subset(absloss_c_bfr, val > 0)) / nrow(absloss_c_bfr)


## Median ----
median(absloss_p$val)
median(absloss_c$val)

## Divergence ----

#JS divergence (sensitive to number/width of bins!)
#Rice (1944) suggested 2 * (n_samp ^ (1/3))
#https://stats.stackexchange.com/questions/510699/discrete-kl-divergence-with-decreasing-bin-width
#How about PSI (population stability index)? Seems similar to JS divergence
#https://arize.com/blog-course/population-stability-index-psi/

bin_digit = floor(log10(max(absloss_p$val, absloss_c$val) / 12))
break_max = ceiling(max(absloss_p$val, absloss_c$val) / (10 ^ bin_digit)) * (10 ^ bin_digit)

JSTest(subset(absloss_p, val > 0)$val, subset(absloss_p, val > 0)$val, bins = seq(0, break_max, len = 13))


# 4. Fit and plot C loss distributions ----
absloss_list = list(P = subset(absloss_p, val > 0)$val, C = subset(absloss_c, val > 0)$val)
absloss_list_bfr = list(P = subset(absloss_p_bfr, val > 0)$val, C = subset(absloss_c_bfr, val > 0)$val)

## Single parametric distributions: exponential, Gaussian
absloss_fit_exp = lapply(absloss_list, function(x) MASS::fitdistr(x, "exponential"))
absloss_fit_nor = lapply(absloss_list, function(x) MASS::fitdistr(x, "normal"))

#generate sampled values from fitted distributions
samp_absloss_p_exp = rexp(length(absloss_list$P), absloss_fit_exp$P$estimate)
samp_absloss_c_exp = rexp(length(absloss_list$C), absloss_fit_exp$C$estimate)
samp_absloss_df_exp = rbind(data.frame(type = "C", val = samp_absloss_c_exp),
                            data.frame(type = "P", val = samp_absloss_p_exp))

samp_absloss_p_nor = rnorm(length(absloss_list$P), absloss_fit_nor$P$estimate[1], absloss_fit_nor$P$estimate[2])
samp_absloss_c_nor = rnorm(length(absloss_list$C), absloss_fit_nor$C$estimate[1], absloss_fit_nor$C$estimate[2])
samp_absloss_p_nor = samp_absloss_p_nor[which(samp_absloss_p_nor > 0)]
samp_absloss_c_nor = samp_absloss_c_nor[which(samp_absloss_c_nor > 0)]
samp_absloss_df_nor = rbind(data.frame(type = "C", val = samp_absloss_c_nor),
                            data.frame(type = "P", val = samp_absloss_p_nor))

#test goodness of fit
gof_exp_p = EnvStats::gofTest(absloss_list$P, samp_absloss_p_exp)
gof_exp_c = EnvStats::gofTest(absloss_list$C, samp_absloss_c_exp)
paste0(signif(gof_exp_p$statistic, 2),
       " (", ifelse(signif(gof_exp_p$p.value, 2) < 0.001, "< 0.001", signif(gof_exp_p$p.value, 2)), ")")
paste0(signif(gof_exp_c$statistic, 2),
       " (", ifelse(signif(gof_exp_c$p.value, 2) < 0.001, "< 0.001", signif(gof_exp_c$p.value, 2)), ")")

gof_nor_p = EnvStats::gofTest(absloss_list$P, samp_absloss_p_nor)
gof_nor_c = EnvStats::gofTest(absloss_list$C, samp_absloss_c_nor)
paste0(signif(gof_nor_p$statistic, 2),
       " (", ifelse(signif(gof_nor_p$p.value, 2) < 0.001, "< 0.001", signif(gof_nor_p$p.value, 2)), ")")
paste0(signif(gof_nor_c$statistic, 2),
       " (", ifelse(signif(gof_nor_c$p.value, 2) < 0.001, "< 0.001", signif(gof_nor_c$p.value, 2)), ")")

## Gaussian mixture models
absloss_fit_gmm = lapply(absloss_list, function(x) mclust::Mclust(x, 2))
absloss_fit_gmm_bfr = lapply(absloss_list_bfr, function(x) mclust::Mclust(x, 2))

par(mfrow = c(2, 1))
lapply(absloss_fit_gmm, function(x) plot(x, what = "classification", main = "Mclust Classification"))
lapply(absloss_fit_gmm_bfr, function(x) plot(x, what = "classification", main = "Mclust Classification"))
par(mfrow = c(1, 1))

#generate sampled values from fitted distributions
SampByPro = function(params, n){
  if(length(params$variance$sigmasq) == 1) params$variance$sigmasq[2] = params$variance$sigmasq[1]
  samp_vec = NULL
  while(length(samp_vec) < n) {
    distr_chosen = ifelse(runif(1) <= params$pro[1], 1, 2)
    val = rnorm(1, mean = params$mean[distr_chosen], sd = sqrt(params$variance$sigmasq[distr_chosen]))
    if(val > 0) samp_vec = c(samp_vec, val)
  }
  return(samp_vec)
}

samp_absloss_p_gmm = SampByPro(absloss_fit_gmm$P$parameters, n = max(sapply(absloss_list, length)))
samp_absloss_c_gmm = SampByPro(absloss_fit_gmm$C$parameters, n = max(sapply(absloss_list, length)))
samp_absloss_p_bfr_gmm = SampByPro(absloss_fit_gmm_bfr$P$parameters, n = max(sapply(absloss_list_bfr, length)))
samp_absloss_c_bfr_gmm = SampByPro(absloss_fit_gmm_bfr$C$parameters, n = max(sapply(absloss_list_bfr, length)))

#test goodness of fit
gof_gmm_p = EnvStats::gofTest(absloss_list$P, samp_absloss_p_gmm)
gof_gmm_c = EnvStats::gofTest(absloss_list$C, samp_absloss_c_gmm)
paste0(signif(gof_gmm_p$statistic, 2),
       " (", ifelse(signif(gof_gmm_p$p.value, 2) < 0.001, "< 0.001", signif(gof_gmm_p$p.value, 2)), ")")
paste0(signif(gof_gmm_c$statistic, 2),
       " (", ifelse(signif(gof_gmm_c$p.value, 2) < 0.001, "< 0.001", signif(gof_gmm_c$p.value, 2)), ")")

gof_gmm_p_bfr = EnvStats::gofTest(absloss_list_bfr$P, samp_absloss_p_bfr_gmm)
gof_gmm_c_bfr = EnvStats::gofTest(absloss_list_bfr$C, samp_absloss_c_bfr_gmm)
paste0(signif(gof_gmm_p_bfr$statistic, 2),
       " (", ifelse(signif(gof_gmm_p_bfr$p.value, 2) < 0.001, "< 0.001", signif(gof_gmm_p_bfr$p.value, 2)), ")")
paste0(signif(gof_gmm_c_bfr$statistic, 2),
       " (", ifelse(signif(gof_gmm_c_bfr$p.value, 2) < 0.001, "< 0.001", signif(gof_gmm_c_bfr$p.value, 2)), ")")


#visualize mixing proportions
round(absloss_fit_gmm$P$parameters$mean, 1)
round(sqrt(absloss_fit_gmm$P$parameters$variance$sigmasq), 1)
round(absloss_fit_gmm$P$parameters$pro, 2)

round(absloss_fit_gmm$C$parameters$mean, 1)
round(sqrt(absloss_fit_gmm$C$parameters$variance$sigmasq), 1)
round(absloss_fit_gmm$C$parameters$pro, 2)

round(absloss_fit_gmm_bfr$P$parameters$mean, 1)
round(sqrt(absloss_fit_gmm_bfr$P$parameters$variance$sigmasq), 1)
round(absloss_fit_gmm_bfr$P$parameters$pro, 2)

round(absloss_fit_gmm_bfr$C$parameters$mean, 1)
round(sqrt(absloss_fit_gmm_bfr$C$parameters$variance$sigmasq), 1)
round(absloss_fit_gmm_bfr$C$parameters$pro, 2)


#plot
absloss_df = rbind(data.frame(val = absloss_list$P, type = "P"),
                   data.frame(val = absloss_list$C, type = "C"))
absloss_df_bfr = rbind(data.frame(val = absloss_list_bfr$P, type = "P"),
                   data.frame(val = absloss_list_bfr$C, type = "C"))

nbins = 20

max_val = max(sapply(absloss_list, max))
max_plot = ceiling(max_val / 10 ^ floor(log10(max_val))) * 10 ^ floor(log10(max_val))
binwd = max_plot / nbins

ggplot(data = subset(absloss_df, type == "P"), aes(val)) +
  geom_freqpoly(aes(y = stat(density)), lwd = 1, binwidth = binwd, boundary = 0) +
  geom_freqpoly(data = subset(samp_absloss_df_exp, type == "P"), aes(y = stat(density)), lwd = 0.5, col = "skyblue", binwidth = binwd, boundary = 0) +
  geom_freqpoly(data = subset(samp_absloss_df_nor, type == "P"), aes(y = stat(density)), lwd = 0.5, col = "darkblue", binwidth = binwd, boundary = 0) +
  scale_x_continuous(name = "C loss (Mg CO2e)", limits = c(0, max_plot)) +
  scale_y_continuous(name = "Density") +
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18))
ggsave(paste0(file_path, site, "_4_loss_distr_fit_project.png"), width = 15, height = 10, unit = "cm")

ggplot(data = subset(absloss_df, type == "C"), aes(val)) +
  geom_freqpoly(aes(y = stat(density)), lwd = 1, binwidth = binwd, boundary = 0) +
  geom_freqpoly(data = subset(samp_absloss_df_exp, type == "C"), aes(y = stat(density)), lwd = 0.5, col = "skyblue", binwidth = binwd, boundary = 0) +
  geom_freqpoly(data = subset(samp_absloss_df_nor, type == "C"), aes(y = stat(density)), lwd = 0.5, col = "darkblue", binwidth = binwd, boundary = 0) +
  scale_x_continuous(name = "C loss (Mg CO2e)", limits = c(0, max_plot)) +
  scale_y_continuous(name = "Density") +
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18))
ggsave(paste0(file_path, site, "_4_loss_distr_fit_counterfactual.png"), width = 15, height = 10, unit = "cm")


ggplot(data = subset(absloss_df, type == "P"), aes(val)) +
  geom_freqpoly(aes(y = stat(density)), lwd = 1, binwidth = binwd, boundary = 0) +
  geom_freqpoly(data = data.frame(type = "P", val = samp_absloss_p_gmm), aes(y = stat(density)), lwd = 0.5, color = "blue", binwidth = binwd, boundary = 0) +
  scale_x_continuous(name = "C loss (Mg CO2e)", limits = c(0, max_plot)) +
  scale_y_continuous(name = "Density") +
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18))
ggsave(paste0(file_path, site, "_4_loss_distr_gmm_project.png"), width = 15, height = 10, unit = "cm")

ggplot(data = subset(absloss_df, type == "C"), aes(val)) +
  geom_freqpoly(aes(y = stat(density)), lwd = 1, binwidth = binwd, boundary = 0) +
  geom_freqpoly(data = data.frame(type = "C", val = samp_absloss_c_gmm), aes(y = stat(density)), lwd = 0.5, color = "blue", binwidth = binwd, boundary = 0) +
  scale_x_continuous(name = "C loss (Mg CO2e)", limits = c(0, max_plot)) +
  scale_y_continuous(name = "Density") +
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18))
ggsave(paste0(file_path, site, "_4_loss_distr_gmm_counterfactual.png"), width = 15, height = 10, unit = "cm")


# 5. Create step-wise release schedule and calculate credibility and permanence ----

##Calculate a-bar (5% percentile of additionality distribution) ----
samp_mean = rep(NA, 100)
samp_median = rep(NA, 100)
samp_abar = rep(NA, 100)
for(i in 1:100){
  samp_additionality = SampByPro(absloss_fit_gmm$C$parameters, n = 10000) - 
    SampByPro(absloss_fit_gmm$P$parameters, n = 10000)
  samp_mean[i] = mean(samp_additionality)
  samp_median[i] = median(samp_additionality)
  samp_abar[i] = quantile(samp_additionality, 0.05)
}
paste0(round(mean(samp_median)), " [", round(t.test(samp_median)$conf.int[1]), " - ", round(t.test(samp_median)$conf.int[2]), "]")
paste0(round(mean(samp_abar)), " [", round(t.test(samp_abar)$conf.int[1]), " - ", round(t.test(samp_abar)$conf.int[2]), "]")

samp_mean_bfr = rep(NA, 100)
samp_median_bfr = rep(NA, 100)
samp_abar_bfr = rep(NA, 100)
for(i in 1:100){
  samp_additionality_bfr = SampByPro(absloss_fit_gmm_bfr$C$parameters, n = 10000) - 
    SampByPro(absloss_fit_gmm_bfr$P$parameters, n = 10000)
  samp_mean_bfr[i] = mean(samp_additionality_bfr)
  samp_median_bfr[i] = median(samp_additionality_bfr)
  samp_abar_bfr[i] = quantile(samp_additionality_bfr, 0.05)
}

paste0(round(mean(samp_median_bfr)), " [", round(t.test(samp_median_bfr)$conf.int[1]), " - ", round(t.test(samp_median_bfr)$conf.int[2]), "]")
paste0(round(mean(samp_abar_bfr)), " [", round(t.test(samp_abar_bfr)$conf.int[1]), " - ", round(t.test(samp_abar_bfr)$conf.int[2]), "]")


samp_additionality = SampByPro(absloss_fit_gmm$C$parameters, n = 10000) - 
  SampByPro(absloss_fit_gmm$P$parameters, n = 10000)
a_bar = quantile(samp_additionality, 0.05)
samp_additionality_bfr = SampByPro(absloss_fit_gmm_bfr$C$parameters, n = 10000) - 
  SampByPro(absloss_fit_gmm_bfr$P$parameters, n = 10000)
a_bar_bfr = quantile(samp_additionality_bfr, 0.05)
ggplot(data = data.frame(val = samp_additionality), aes(val)) +
  geom_freqpoly(aes(y = stat(density)), lwd = 1, binwidth = binwd, boundary = 0) +
  geom_freqpoly(data = data.frame(val = samp_additionality_bfr), aes(y = stat(density)), lwd = 1, binwidth = binwd, boundary = 0, color = "red") +
  geom_vline(xintercept = 0, lty = "dashed", color = "gray") +
  geom_vline(xintercept = a_bar) +
  geom_vline(xintercept = a_bar_bfr, color = "red") +
  scale_x_continuous(name = "Additionality (Mg CO2e)") +
  scale_y_continuous(name = "Density") +
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18))
ggsave(paste0(file_path, site, "_add_distr_abar.png"), width = 15, height = 20, unit = "cm")


##Perform simulations (static a-bar) ----
mean_ep_vec = rep(NA, 100)
cred_vec = rep(NA, 100)

for(n_rep in 1:100){
  H = 50 #evaluation horizon; project duration
  H_rel = 10 * H #release horizon
  sim_p_loss = rep(0, H)
  sim_c_loss = rep(0, H)
  sim_additionality = rep(0, H)
  sim_credit = rep(0, H)
  sim_benefit = rep(0, H)
  sim_release = rep(0, H_rel)
  sim_damage = rep(0, H)
  sim_ep = rep(0, H)
  sim_credibility = rep(1, H)
  D = 0.03 #discount rate
  current_year = lubridate::year(lubridate::now())
  scc_current = scc[which(names(scc) == current_year):length(scc)]
  for(i in 1:H){
    year = current_year + i - 1
    sim_p_loss[i] = sample(samp_absloss_p_gmm, 1)
    sim_c_loss[i] = sample(samp_absloss_c_gmm, 1)
    sim_additionality[i] = sim_c_loss[i] - sim_p_loss[i]
    sim_credit[i] = sim_release[i] + sim_additionality[i]
    
    if(sim_credit[i] > 0){
      sim_benefit[i] = sim_credit[i] * scc_current[i]
      to_be_released = sim_credit[i]
      j = 0 #number of years after i
      R = rep(0, H_rel - i) #release schedule
      while(to_be_released > 0 & i + j < H_rel){
        j = j + 1
        max_release = ifelse(i + j <= H, max(-a_bar, 0), max(-a_bar_bfr, 0)) #maximum release needed annually
        R[j] = min(to_be_released, max_release - sim_release[i + j])
        to_be_released = to_be_released - R[j]
        sim_release[i + j] = sim_release[i + j] + R[j]
      }
      sim_damage[i] = sum(R * scc_current[(i + 1):H_rel] / ((1 + D) ^ (1:length(R))))
      sim_ep[i] = (sim_benefit[i] - sim_damage[i]) / sim_benefit[i]
    } else if(sim_credit[i] < 0){
      sim_credibility[i] = 0
    }
  }
  
  mean_ep_vec[n_rep] = mean(sim_ep)
  cred_vec[n_rep] = sum(sim_credibility) / H
}

sim_df = rbind(data.frame(var = "Project stock", val = subset(summ_stock, year == t0 & series == "mean" & var == "treatment")$val - cumsum(sim_p_loss)),
               data.frame(var = "Counterfactual stock", val = subset(summ_stock, year == t0 & series == "mean" & var == "control")$val - cumsum(sim_c_loss)),
               data.frame(var = "Project loss", val = sim_p_loss),
               data.frame(var = "Counterfactual loss", val = sim_c_loss),
               data.frame(var = "Additionality", val = sim_additionality),
               data.frame(var = "Credit", val = sim_credit),
               data.frame(var = "Anticipated_release", val = sim_release[1:H]),
               data.frame(var = "Credit_benefit", val = sim_benefit),
               data.frame(var = "Damage", val = sim_damage),
               data.frame(var = "EP", val = sim_ep),
               data.frame(var = "Credibility", val = sim_credibility)) %>%
  mutate(year = rep(current_year:(current_year + H - 1), length(unique(var))))

ggplot(data = subset(sim_df, var %in% c("Project stock", "Counterfactual stock"))) +
  geom_line(aes(x = year, y = val, col = var)) +
  scale_color_manual(values = c("red", "black"), guide = "none") + 
  scale_x_continuous(name = "Year", breaks = seq(current_year, current_year + H, by = 10)) + 
  scale_y_continuous(name = "Carbon stock (Mg CO2e)") + 
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, site, "_6a_sim_stock.png"), width = 15, height = 15, unit = "cm")

ggplot(data = subset(sim_df, var %in% c("Additionality", "Anticipated_release", "Credit"))) +
  geom_line(aes(x = year, y = val, col = var, lty = var)) +
  geom_point(data = subset(sim_df, year %in% subset(sim_df, var == "Credibility" & val != 1)$year & var == "Credit"), aes(x = year, y = val), col = "red") +
  geom_hline(yintercept = 0, lwd = 1.2) +
  scale_color_manual(values = c("black", "red", "blue"), guide = "none") + 
  scale_linetype_manual(values = c(1, 2, 1), guide = "none") + 
  scale_x_continuous(name = "Year", breaks = seq(current_year, current_year + H, by = 10)) + 
  scale_y_continuous(name = "Carbon (Mg CO2e)") + 
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, site, "_6c_sim_credit.png"), width = 15, height = 15, unit = "cm")


ggplot(data = subset(sim_df, var == "EP")) +
  geom_bar(stat = "identity", aes(x = year, y = val)) +
  geom_point(data = subset(sim_df, var == "Credibility"), aes(x = year, y = -val), col = "red") +
  scale_x_continuous(name = "Year", breaks = seq(current_year, current_year + H, by = 10)) + 
  scale_y_continuous(name = "eP", limits = c(0, 1)) + 
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 18, angle = 45),
        plot.margin = margin(0.7, 1, 0.7, 0.7, "cm"))
ggsave(paste0(file_path, site, "_6e_sim_ep.png"), width = 15, height = 15, unit = "cm")

paste0(signif(summary(mean_ep_vec)[3], 2), " [", signif(summary(mean_ep_vec)[2], 2), " - ", signif(summary(mean_ep_vec)[5], 2), "]")
paste0(signif(summary(cred_vec)[3], 2), " [", signif(summary(cred_vec)[2], 2), " - ", signif(summary(cred_vec)[5], 2), "]")
length(which(cred_vec < 0.95)) / length(cred_vec)
