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
rmarkdown::render("/Users/E-Ping Rau/OneDrive - University of Cambridge/4C_evaluations/R/Reports/eval_template/evaluation_epingrau.Rmd", clean = FALSE)
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

summ_flux_rel = rbind(summariseSeries(flux_series_sim_rel, "rel_flux_treatment"),
                      summariseSeries(flux_series_sim_rel, "rel_flux_control")) %>%
  mutate(val = val * -1)

## Plot time series ----
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

#absolute carbon loss
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

#carbon loss rate
# summ_flux_rel_perc = summ_flux_rel %>% mutate(val = val * 100)
# PlotTimeSeries(dat = summ_flux_rel_perc,
#                #lim_to_use = c(floor(min(dat$val) / (10 ^ 5)) * (10 ^ 5), 0)
#                y_label = "C loss rate (%)",
#                y_lim = c(0, 2), #for comparison between the five sites
#                file_suffix = "_2b_rel_time_series.png")


# 3. Describe empirical distributions ----

#absolute carbon loss
absloss_p = -subset(summ_flux, var == "treatment_proj" & year >= t0 & val < 0 & series != "mean")$val
absloss_c = -subset(summ_flux, var == "control_proj" & year >= t0 & val < 0 & series != "mean")$val

## Epsilon (probability of carbon loss) ----
length(absloss_p) / nrow(subset(summ_flux, var == "treatment_proj" & year >= t0 & series != "mean"))
length(absloss_c) / nrow(subset(summ_flux, var == "control_proj" & year >= t0 & series != "mean"))

## Median ----
median(absloss_p)
median(absloss_c)

## Divergence ----

#JS divergence (sensitive to number/width of bins!)
#Rice (1944) suggested 2 * (n_samp ^ (1/3))
#https://stats.stackexchange.com/questions/510699/discrete-kl-divergence-with-decreasing-bin-width
#How about PSI (population stability index)? Seems similar to JS divergence
#https://arize.com/blog-course/population-stability-index-psi/

bin_digit = floor(log10(max(absloss_p, absloss_c) / 12))
break_max = ceiling(max(absloss_p, absloss_c) / (10 ^ bin_digit)) * (10 ^ bin_digit)

JSTest(absloss_p, absloss_c, bins = seq(0, break_max, len = 13))

#KS distance and bootstrap KS test
#https://www.r-bloggers.com/2020/02/monitoring-for-changes-in-distribution-with-resampling-tests/
ks_res = ks.boot(absloss_p, absloss_c, nboots = 10000)
ks_res$ks$statistic
#KS statistic isn't a metric in the formal sense: better not report this
#https://arize.com/blog-course/kolmogorov-smirnov-test/

## Test for stationarity: KS distance of moving-window subsamples ----
absloss_p_long = subset(summ_flux, var == "treatment_proj" & year >= t0 & series != "mean")
absloss_c_long = subset(summ_flux, var == "control_proj" & year >= t0 & series != "mean")

lme_p = lmer(val ~ 1 + (1|year), data = absloss_p_long)
plot(lme_p, resid(., scaled = T) ~ fitted(.) | year, abline = 0) #standardized residuals versus fitted values by year
plot(lme_p, year ~ resid(., scaled = T)) #box-plots of residuals by year
plot(lme_p, sqrt(abs(resid(.))) ~ fitted(.),
     type = c("p", "smooth"),
     par.settings = list(plot.line = list(alpha = 1, col = "red",
                                          lty = 1, lwd = 2))) #scale-location plot

year_vec = t0:2021
ks_vec = rep(NA, length(year_vec))
pval_vec = rep(NA, length(year_vec))
for(i in 1:length(year_vec)){
  ks_res = ks.boot(absloss_p_long$val, subset(absloss_p_long, year == year_vec[i])$val, nboots = 10000)
  ks_vec[i] = ks_res$ks$statistic
  pval_vec[i] = ks_res$ks.boot.pvalue
}
par(mar = c(5, 4, 4, 5) + 0.1)
plot(1:length(year_vec), ks_vec, ylim = c(0, 1), xlab = "Year", type = "l", xaxt = "n")
axis(1, at = 1:length(year_vec), labels = year_vec)
axis(4, at = seq(0, 1, by = 0.2))
mtext("p value", side = 4, line = 3)
lines(1:length(year_vec), pval_vec, ylim = c(0, 1), col = "red")
abline(h = 0.05, col = "red", lty = 2)

## Test for stationarity ----
#trend-stationary: https://en.wikipedia.org/wiki/Trend-stationary_process
#https://www.r-econometrics.com/timeseries/stationarity/
#https://uk.mathworks.com/help/econ/trend-stationary-vs-difference-stationary.html
#https://python.plainenglish.io/time-series-analysis-mastering-the-concepts-of-stationarity-c9fc489893cf

absloss_p_tseries = subset(summ_flux, var == "treatment_proj" & year >= t0 & series != "mean") %>%
  pivot_wider(names_from = "series", values_from = "val") %>%
  mutate(year = NULL, var = NULL)
absloss_c_tseries = subset(summ_flux, var == "control_proj" & year >= t0 & series != "mean") %>%
  pivot_wider(names_from = "series", values_from = "val") %>%
  mutate(year = NULL, var = NULL)

absloss_p_tseries_adf = sapply(absloss_p_tseries, function(x) tseries::adf.test(x)$p.value)
absloss_c_tseries_adf = sapply(absloss_c_tseries, function(x) tseries::adf.test(x)$p.value)
length(which(absloss_p_tseries_adf >= 0.05)) / 20
length(which(absloss_c_tseries_adf >= 0.05)) / 20

absloss_p_tseries_kpss_level = sapply(absloss_p_tseries, function(x) tseries::kpss.test(ts(x), null = "Level")$p.value)
absloss_c_tseries_kpss_level = sapply(absloss_c_tseries, function(x) tseries::kpss.test(x, null = "Level")$p.value)
absloss_p_tseries_kpss_trend = sapply(absloss_p_tseries, function(x) tseries::kpss.test(x, null = "Trend")$p.value)
absloss_c_tseries_kpss_trend = sapply(absloss_c_tseries, function(x) tseries::kpss.test(x, null = "Trend")$p.value)
length(which(absloss_p_tseries_kpss_level < 0.05)) / 20
length(which(absloss_c_tseries_kpss_level < 0.05)) / 20
length(which(absloss_p_tseries_kpss_trend < 0.05)) / 20
length(which(absloss_c_tseries_kpss_trend < 0.05)) / 20

library(Rbeast)
out = beast(absloss_p_tseries$V1, season='none')
plot(out)

library(gam)
absloss_p_df = subset(summ_flux, var == "treatment_proj" & year >= t0 & series != "mean") %>%
  mutate(series = as.factor(series))

gam_out = gam(val ~ s(year) + series, data = absloss_p_df)
par(mfrow = c(1, 2))
plot(gam_out, se = T)
summary(gam_out)

ggplot(data = subset(summ_flux, var == "treatment_proj" & year >= t0 & series != "mean")) +
  geom_line(aes(x = year, y = val, group = series))




# 4. Fit distributions
# relloss_list = list(Project = relloss_p, Counterfactual = relloss_c)
# relloss_fit = list(exp = lapply(relloss_list, function(x) fitdistrplus::fitdist(x, "exp")),
#                    weibull = lapply(relloss_list, function(x) fitdistrplus::fitdist(x, "weibull", lower = c(0, 0))),
#                    gamma = lapply(relloss_list, function(x) fitdistrplus::fitdist(x, "gamma")),
#                    lnorm = lapply(relloss_list, function(x) fitdistrplus::fitdist(x, "lnorm")),
#                    beta = lapply(relloss_list, function(x) fitdistrplus::fitdist(x, "beta")),
#                    norm = lapply(relloss_list, function(x) MASS::fitdistr(x, "normal")))

#additionality_fit_nor = MASS::fitdistr(additionality_series, "normal")

##generate random numbers ----
# relloss_rand_exp = lapply(relloss_fit_exp, function(x) rexp(10000, x$estimate))
# relloss_rand_nor = lapply(relloss_fit_nor, function(x) rnorm(10000, x$estimate[1], x$estimate[2]))
# relloss_rand_nor = lapply(relloss_rand_nor, function(x) x[which(x > 0)])
# 
# absloss_rand_nor = lapply(absloss_fit_nor, function(x) rnorm(10000, x$estimate[1], x$estimate[2]))
# absloss_rand_nor = lapply(absloss_rand_nor, function(x) x[which(x > 0)])
# #exponential, gamma and lognormal results in values larger than 1
# #normal distribution needs to be truncated at zero
# 
# additionality_rand_nor = rnorm(10000, additionality_fit_nor$estimate[1], additionality_fit_nor$estimate[2])
# additionality_rand_nor = additionality_rand_nor[which(additionality_rand_nor > 0)]
# 
#
##examine with p-p plots ----
# relloss_cdf_emp = lapply(relloss_list, ecdf)
# relloss_cdf_fit_exp = lapply(relloss_rand_exp, ecdf)
# relloss_cdf_fit_nor = lapply(relloss_rand_nor, ecdf)
# 
# absloss_cdf_emp = lapply(absloss_list, ecdf)
# absloss_cdf_fit_nor = lapply(absloss_rand_nor, ecdf)
# 
# additionality_cdf_emp = ecdf(additionality_series)
# additionality_cdf_fit_nor = ecdf(additionality_rand_nor)
# 
# 
# plotPP = function(func1, func2, dat, title, wd = 1200, ht = 600, filename){
#   png(paste0(file_path, site, filename), width = wd, height = ht)
#   par(mfrow = c(1, 2), mar = c(6, 6, 6, 2) + 0.1)
#   mapply(function(func1, func2, dat, title) {
#     plot(func1(sort(dat)), func2(sort(dat)),
#          main = title, xlab = "Theoretical cumulative probabilities", ylab = "Sample cumulative probabilities",
#          cex.main = 3, cex.lab = 2.5, cex.axis = 2)
#     lines(c(0, 1), c(0, 1))},
#     func1 = func1, func2 = func2, dat = dat, title = title)
#   dev.off()
# }
#
# plotPP(func1 = relloss_cdf_fit_exp, func2 = relloss_cdf_emp, dat = relloss_list,
#        title = relloss_type, filename = "_relloss_fig3a_fit_pp_exp.png")
# 
# plotPP(func1 = absloss_cdf_fit_nor, func2 = absloss_cdf_emp, dat = absloss_list,
#        title = absloss_type, filename = "_absloss_fig3f_fit_pp_norm.png")
# 
# 
# png(paste0(file_path, site, "_additionality_pp_norm.png"), width = 600, height = 600)
# par(mar = c(6, 6, 6, 2) + 0.1)
# plot(additionality_cdf_fit_nor(sort(additionality_series)), additionality_cdf_emp(sort(additionality_series)),
#      main = "Additionality", xlab = "Theoretical cumulative probabilities", ylab = "Sample cumulative probabilities",
#      cex.main = 3, cex.lab = 2.5, cex.axis = 2)
# lines(c(0, 1), c(0, 1))
# dev.off()
# 
# ##examine with q-q plots ----
# plotQQ = function(x, y, distr, title, wd = 1200, ht = 600, filename){
#   png(paste0(file_path, site, filename), width = wd, height = ht)
#   par(mfrow = c(1, 2), mar = c(6, 6, 6, 2) + 0.1, cex.main = 2.5, cex.lab = 2.25, cex.axis = 2)
#   mapply(function(x, y, distr, title) {
#     params = as.list(y$estimate)
#     do.call(car::qqPlot,
#             c(list(x, distribution = distr, main = title,
#                    xlab = "Theoretical quantiles", ylab = "Sample quantiles"), params))},
#     x = x, y = y, distr = distr, title = title)
#   dev.off()
# }

# plotQQ(x = relloss_list, y = relloss_fit_exp, distr = "exp",
#        title = relloss_type, filename = "_relloss_fig4a_fit_qq_exp.png")
#
# plotQQ(x = absloss_list, y = absloss_fit_nor, distr = "norm",
#        title = absloss_type, filename = "_absloss_fig4f_fit_qq_norm.png")
# 
# 
# png(paste0(file_path, site, "_additionality_fit_qq_norm.png"), width = 600, height = 600)
# par(mfrow = c(1, 2), mar = c(6, 6, 6, 2) + 0.1, cex.main = 2.5, cex.lab = 2.25, cex.axis = 2)
# car::qqPlot(additionality_series, distribution = "norm", mean = additionality_fit_nor$estimate[1], sd = additionality_fit_nor$estimate[2],
#             main = "Additionality", xlab = "Theoretical quantiles", ylab = "Sample quantiles")
# dev.off()
# 
# 
# ##test goodness-of-fit ----
# ChisqGof = function(dat, fit, distr){
#   a = mapply(function(x, y, distr) {
#     EnvStats::gofTest(x, test = "chisq", distribution = distr, param.list = as.list(y$estimate))
#   },
#   x = dat, y = fit, distr = distr, SIMPLIFY = F)
#   return(a)
# }
# 
# gammaScale = function(x){
#   a = list(estimate = c(as.list(x$estimate), scale = 1 / as.list(x$estimate)$rate))
#   a$estimate["rate"] = NULL
#   return(a)
# }
# 
# absloss_chisq_nor = ChisqGof(dat = absloss_list, fit = absloss_fit_nor, distr = "norm")
# 
# sapply(absloss_chisq_nor, function(x) list(stat = x$statistic, pv = x$p.value))
# 
# additionality_chisq_nor = EnvStats::gofTest(additionality_series, test = "chisq", distribution = "norm", param.list = as.list(additionality_fit_nor$estimate))
# (list(stat = additionality_chisq_nor$statistic, pv = additionality_chisq_nor$p.value))
# 
# #weibull and beta are good for project area in Gola, none of them are good for the rest
# relloss_pa_wbl = fitdistrplus::fitdist(relloss_pa * 100, "weibull")
# relloss_pa_bet = fitdistrplus::fitdist(relloss_pa * 100, "beta")


# 4. Fit and plot C loss distributions ----
absloss_list = list(P = absloss_p, C = absloss_c)
absloss_df = rbind(data.frame(type = "P", val = absloss_p),
                   data.frame(type = "C", val = absloss_c))

absloss_fit_exp = lapply(absloss_list, function(x) MASS::fitdistr(x, "exponential"))
samp_absloss_p_exp = rexp(10000, absloss_fit_exp$P$estimate)
samp_absloss_c_exp = rexp(10000, absloss_fit_exp$C$estimate)
samp_absloss_df_exp = rbind(data.frame(type = "C", val = samp_absloss_c_exp),
                            data.frame(type = "P", val = samp_absloss_p_exp))

absloss_fit_nor = lapply(absloss_list, function(x) MASS::fitdistr(x, "normal"))
samp_absloss_p_nor = rnorm(10000, absloss_fit_nor$P$estimate[1], absloss_fit_nor$P$estimate[2])
samp_absloss_c_nor = rnorm(10000, absloss_fit_nor$C$estimate[1], absloss_fit_nor$C$estimate[2])
samp_absloss_df_nor = rbind(data.frame(type = "C", val = samp_absloss_c_nor),
                            data.frame(type = "P", val = samp_absloss_p_nor))


absloss_fit_gmm = lapply(absloss_list, function(x) mclust::Mclust(x, 2))

lapply(absloss_fit_gmm, function(x) x$parameters$pro)

par(mfrow = c(2, 1))
lapply(absloss_fit_gmm, function(x) plot(x, what = "classification", main = "Mclust Classification"))
par(mfrow = c(1, 1))

samp_absloss_p_1 = rnorm(10000 * absloss_fit_gmm$P$parameters$pro[1], mean = absloss_fit_gmm$P$parameters$mean[1], sd = sqrt(absloss_fit_gmm$P$parameters$variance$sigmasq[1]))
samp_absloss_p_2 = rnorm(10000 * absloss_fit_gmm$P$parameters$pro[2], mean = absloss_fit_gmm$P$parameters$mean[2], sd = sqrt(absloss_fit_gmm$P$parameters$variance$sigmasq[2]))
samp_absloss_c_1 = rnorm(10000 * absloss_fit_gmm$C$parameters$pro[1], mean = absloss_fit_gmm$C$parameters$mean[1], sd = sqrt(absloss_fit_gmm$C$parameters$variance$sigmasq[1]))
samp_absloss_c_2 = rnorm(10000 * absloss_fit_gmm$C$parameters$pro[2], mean = absloss_fit_gmm$C$parameters$mean[2], sd = sqrt(absloss_fit_gmm$C$parameters$variance$sigmasq[2]))
samp_absloss_p_1 = samp_absloss_p_1[samp_absloss_p_1 >= 0]
samp_absloss_p_2 = samp_absloss_p_2[samp_absloss_p_2 >= 0]
samp_absloss_p_gmm = c(samp_absloss_p_1, samp_absloss_p_2)
samp_absloss_c_1 = samp_absloss_c_1[samp_absloss_c_1 >= 0]
samp_absloss_c_2 = samp_absloss_c_2[samp_absloss_c_2 >= 0]
samp_absloss_c_gmm = c(samp_absloss_c_1, samp_absloss_c_2)


gof_p = EnvStats::gofTest(absloss_p, samp_absloss_p_gmm)
gof_c = EnvStats::gofTest(absloss_c, samp_absloss_c_gmm)
signif(gof_p$statistic, 2)
signif(gof_p$p.value, 2)
signif(gof_c$statistic, 2)
signif(gof_c$p.value, 2)

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



# PlotDistr(dat = relloss_df,
#           dat_fit = NULL,
#           nbins = 25,
#           x_label = "C loss rate (%)",
#           file_suffix = "_5b_rel_distr.png")

# 5. Fit and plot release distributions ----
additionality_series = subset(summ_flux, var == "additionality" & series != "mean" & year > t0)$val #additionality based on difference of empirical time series

release = additionality_series[which(additionality_series > 0)]
release_fit_exp = MASS::fitdistr(release, densfun = "exponential")
release_fit_wbl = MASS::fitdistr(release, densfun = "weibull")
samp_release_exp = rexp(10000, release_fit_exp$estimate)
samp_release_wbl = rweibull(10000, release_fit_wbl$estimate[1], release_fit_wbl$estimate[2])

max_val = max(abs(range(release)))
max_plot = ceiling(max_val / 10 ^ floor(log10(max_val))) * 10 ^ floor(log10(max_val))
binwd = max_plot / nbins

ggplot(data = data.frame(val = additionality_series), aes(val)) +
  geom_freqpoly(aes(y = stat(density)), lwd = 1, binwidth = binwd, boundary = 0) +
  geom_freqpoly(data = data.frame(val = samp_release_wbl), aes(y = stat(density)), lwd = 0.5, lty = 3, binwidth = binwd, boundary = 0) +
  scale_x_continuous(name = "Additionality (Mg CO2e)", limits = c(-max_plot, max_plot)) +
  scale_y_continuous(name = "Density") +
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18))
ggsave(paste0(file_path, site, "_5_release_distr.png"), width = 15, height = 20, unit = "cm")

# 6. Create step-wise release schedule and calculate credibility and permanence ----
release_factor = 1
release_factor_df = data.frame(rel_factor = NULL, val = NULL)

for(n in 1:10){
  release_factor = 1 * n
  
  median_ep_vec = rep(NA, 100)
  cred_vec = rep(NA, 100)
  
  for(n_rep in 1:100){
    H = 50 #evaluation horizon
    sim_p_loss = rep(0, H)
    sim_c_loss = rep(0, H)
    sim_additionality = rep(0, H)
    sim_credit = rep(0, H)
    sim_credval = rep(0, H)
    sim_release = rep(0, 2 * H)
    sim_damage = rep(0, H)
    sim_ep = rep(0, H)
    sim_credibility = rep(1, H)
    D = 0.03 #discount rate
    current_year = lubridate::year(lubridate::now())
    scc_current = scc[which(names(scc) == current_year):length(scc)]
    for(i in 1:H){
      year = current_year + i - 1
      sim_p_loss[i] = rexp(1, absloss_fit_exp$P$estimate)
      sim_c_loss[i] = rexp(1, absloss_fit_exp$C$estimate)
      sim_additionality[i] = sim_c_loss[i] - sim_p_loss[i]
      sim_credit[i] = sim_release[i] + sim_additionality[i]
      
      if(sim_credit[i] > 0){
        sim_credval[i] = sim_credit[i] * scc_current[i]
        to_be_released = sim_credit[i]
        j = 0 #number of years after i
        R = rep(0, 2 * H - i) #release schedule
        while(to_be_released > 0 & i + j < 2 * H){
          j = j + 1
          R[j] = min(to_be_released, rexp(1, release_fit_exp$estimate * release_factor))
          sim_release[i + j] = sim_release[i + j] + R[j]
          to_be_released = to_be_released - R[j]
        }
        sim_damage[i] = sum(R * scc_current[(i + 1):(2 * H)] / ((1 + D) ^ (1:length(R))))
        sim_ep[i] = (sim_credval[i] - sim_damage[i]) / sim_credval[i]
      } else if(sim_credit[i] < 0){
        sim_credibility[i] = 0
      }
    }
    
    median_ep_vec[n_rep] = mean(sim_ep)
    cred_vec[n_rep] = sum(sim_credibility) / H
  }
  release_factor_df = rbind(release_factor_df, data.frame(rel_factor = release_factor, ep = median_ep_vec, cred = cred_vec))
}

release_factor_df$rel_factor = as.factor(release_factor_df$rel_factor)
ggplot(data = release_factor_df, aes(x = rel_factor)) +
  geom_boxplot(aes(y = ep)) +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 18))

ggplot(data = release_factor_df, aes(x = rel_factor)) +
  geom_boxplot(aes(y = cred)) +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 18))


sim_df = rbind(data.frame(var = "Project stock", val = subset(summ_stock, year == t0 & series == "mean" & var == "treatment")$val - cumsum(sim_p_loss)),
               data.frame(var = "Counterfactual stock", val = subset(summ_stock, year == t0 & series == "mean" & var == "control")$val - cumsum(sim_c_loss)),
               data.frame(var = "Project loss", val = sim_p_loss),
               data.frame(var = "Counterfactual loss", val = sim_c_loss),
               data.frame(var = "Additionality", val = sim_additionality),
               data.frame(var = "Credit", val = sim_credit),
               data.frame(var = "Anticipated_release", val = sim_release[1:H]),
               data.frame(var = "Credit_value", val = sim_credval),
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
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 14, angle = 45),
        plot.margin = margin(0.1, 0.5, 0.1, 0.5, "cm"))
ggsave(paste0(file_path, site, "_6a_sim_stock.png"), width = 15, height = 10, unit = "cm")

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
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 14, angle = 45),
        plot.margin = margin(0.1, 0.5, 0.1, 0.5, "cm"))
ggsave(paste0(file_path, site, "_6c_sim_credit.png"), width = 15, height = 10, unit = "cm")


ggplot(data = subset(sim_df, var == "EP")) +
  geom_bar(stat = "identity", aes(x = year, y = val)) +
  geom_point(data = subset(sim_df, var == "Credibility"), aes(x = year, y = -val), col = "red") +
  scale_x_continuous(name = "Year", breaks = seq(current_year, current_year + H, by = 10)) + 
  scale_y_continuous(name = "eP", limits = c(0, 0.2)) + 
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 14, angle = 45),
        plot.margin = margin(0.1, 0.5, 0.1, 0.5, "cm"))
ggsave(paste0(file_path, site, "_6e_sim_ep.png"), width = 15, height = 10, unit = "cm")


paste0(formatC(median(sim_additionality), format = "e", digits = 1), " [",
       formatC(quantile(sim_additionality, 0.25), format = "e", digits = 1), " - ",
       formatC(quantile(sim_additionality, 0.75), format = "e", digits = 1), "]")
paste0(formatC(median(sim_credit), format = "e", digits = 1), " [",
       formatC(quantile(sim_credit, 0.25), format = "e", digits = 1), " - ",
       formatC(quantile(sim_credit, 0.75), format = "e", digits = 1), "]")
paste0(signif(median(sim_ep), 2), " [",
       signif(quantile(sim_ep, 0.25), 2), " - ",
       signif(quantile(sim_ep, 0.75), 2), "]")
(credibility = sum(sim_credibility) / length(sim_credibility))
