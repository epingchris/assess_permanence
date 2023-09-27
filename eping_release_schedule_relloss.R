library(ggplot2)
library(ggnewscale)
library(tidyverse)
library(magrittr)
library(fitdistrplus)
library(MASS)
library(EnvStats)
library(car)
file_path = "C:/Users/E-Ping Rau/OneDrive - University of Cambridge/carbon_release_pattern/"
source(paste0(file_path, "divergence.R"))

# 1. Initiatlise ----

## 1a. Generate data for new sites ----
rmarkdown::render("/Users/E-Ping Rau/OneDrive - University of Cambridge/4C_evaluations/R/Reports/eval_template/evaluation_epingrau.Rmd", clean = FALSE)
#t0, eval_end, project_end
tmin5 = t0 - 5
years = 1990:eval_end
site_label = site
#eval_classes = c(1, 2, 3, 4) #only evaluate classes 1-4

## 1b. Load existing data for sites ----
site = "WLT_VNCC_KNT" #Gola_country, WLT_VNCC_KNT, CIF_Alto_Mayo, VCS_1396, VCS_934
load(file = paste0(file_path, site, ".Rdata"))

# 2. Calculate carbon stock and flux ----
makeFlux = function(project_series, leakage_series){
  stock_series = aggregate(class_co2e ~ treatment + year, subset(project_series, class %in% eval_classes), FUN = sum)
  stock_wide = as.data.frame(pivot_wider(stock_series, id_cols = year, names_from = "treatment", values_from = "class_co2e"))

  flux_series = data.frame(year = years[-1],
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
PlotTimeSeries = function(dat, lim_to_use, y_label, file_suffix){
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
    scale_y_continuous(name = "C flux (Mg CO2e)", limits = lim_to_use) + 
    ggtitle("") +
    theme_bw() +
    theme(axis.line = element_line(linewidth = 0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "none",
          axis.title = element_text(size = 24),
          axis.text = element_text(size = 16))
  ggsave(paste0(file_path, site_label, file_suffix), width = 15, height = 20, unit = "cm")
}

#absolute carbon loss
PlotTimeSeries(dat = subset(summ_flux, var %in% c("treatment_proj", "control_proj")),
               #lim_to_use = c(floor(min(dat$val) / (10 ^ 5)) * (10 ^ 5), 0)
               lim_to_use = c(-4e+05, 25000), #for comparison between Gola and KNT
               y_label = "C flux (Mg CO2e)",
               file_suffix = "_2a_abs_time_series.png")

#carbon loss rate
summ_flux_rel_perc = summ_flux_rel %>% mutate(val = val * 100)
PlotTimeSeries(dat = summ_flux_rel_perc,
               #lim_to_use = c(floor(min(dat$val) / (10 ^ 5)) * (10 ^ 5), 0)
               lim_to_use = c(-0.2, 1.4), #for comparison between Gola and KNT
               y_label = "C loss rate (%)",
               file_suffix = "_2b_rel_time_series.png")


# 3. Describe empirical distributions ----

#absolute carbon loss
absloss_p = -subset(summ_flux, var == "treatment_proj" & year >= t0 & val < 0)$val
absloss_c = -subset(summ_flux, var == "control_proj" & year >= t0 & val < 0)$val

absloss_df = rbind(data.frame(type = "C", val = absloss_c),
                   data.frame(type = "P", val = absloss_p))

#carbon loss rate
relloss_p = subset(summ_flux_rel, var == "rel_flux_treatment" & year >= t0 & val > 0)$val
relloss_c = subset(summ_flux_rel, var == "rel_flux_control" & year >= t0 & val > 0)$val

relloss_df = rbind(data.frame(type = "C", val = relloss_c),
                   data.frame(type = "P", val = relloss_p)) %>%
  mutate(val = val * 100)

#additionality (based on difference of empirical time series)
additionality_series = subset(summ_flux, var == "additionality" & series != "mean" & year > t0)$val

## Calculate summary statistics ----
#epsilon (probability of carbon loss)
(epsilon_p = length(relloss_p) / nrow(subset(summ_flux_rel, var == "rel_flux_treatment" & year >= t0)))
(epsilon_c = length(relloss_c) / nrow(subset(summ_flux_rel, var == "rel_flux_control" & year >= t0)))

#median, skewness and kurtosis
EmpDistrSumm = function(p, c){
  out_df = cbind(c(median(p), EnvStats::skewness(p), EnvStats::kurtosis(p)),
                          c(median(c), EnvStats::skewness(c), EnvStats::kurtosis(c))) %>%
    `colnames<-`(c("Project", "Counterfactual")) %>%
    `rownames<-`(c("Median", "Skewness", "Kurtosis"))
  return(out_df)
}

(absloss_summary = EmpDistrSumm(p = absloss_p, c = absloss_c))
(relloss_summary = EmpDistrSumm(p = relloss_p, c = relloss_c))

## Plot empirical distributions ----
PlotEmpDistr = function(dat, binwd, x_label, x_lim, file_suffix){
  ggplot(data = dat, aes(val)) +
    geom_freqpoly(aes(col = type), lwd = 1, binwidth = binwd, boundary = 0) +
    scale_color_manual(values = c("red", "black"), labels = c("Counterfactual", "Project")) +
    scale_x_continuous(name = x_label, limits = x_lim) +
    scale_y_continuous(name = "Count") +
    ggtitle("") +
    theme_bw() +
    theme(axis.line = element_line(linewidth = 0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "none",
          axis.title = element_text(size = 24),
          axis.text = element_text(size = 18))
  ggsave(paste0(file_path, site_label, file_suffix), width = 15, height = 20, unit = "cm")
  
}

PlotEmpDistr(dat = absloss_df,
             binwd = 10000,
             x_label = "C loss (Mg CO2e)",
             x_lim = c(0, 4e+05),
             file_suffix = "_3a_abs_emp_distr.png")

PlotEmpDistr(dat = relloss_df,
             binwd = 0.05,
             x_label = "C loss rate (%)",
             x_lim = c(0, 1.5),
             file_suffix = "_3b_rel_emp_distr.png")

ggplot(data = data.frame(type = "add", val = additionality_series), aes(val)) +
  geom_freqpoly(aes(col = type), lwd = 1, binwidth = 15000, boundary = 0) +
  scale_color_manual(values = "black", labels = NULL) +
  scale_x_continuous(name = "Additionality (Mg CO2e)", limits =  c(0, 4e+05)) +
  scale_y_continuous(name = "Count") +
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18))
ggsave(paste0(file_path, site_label, "_3c_add_emp_distr.png"), width = 15, height = 20, unit = "cm")


# 4. Fit distributions ----
relloss_list_full = list(pb = relloss_pb, pa = relloss_pa, cb = relloss_cb, ca = relloss_ca)
relloss_type_full = c("Pre-t0 in project", "Post-t0 in project", "Pre-t0 in counterfactual", "Post-t0 in counterfactual")
relloss_list = list(pa = relloss_pa, ca = relloss_ca)
relloss_type = c("Project", "Counterfactual")

absloss_list = list(pa = absloss_pa, ca = absloss_ca)
absloss_type = c("Project", "Counterfactual")

additionality_type = c("Project", "Counterfactual")


relloss_fit_exp = lapply(relloss_list, function(x) fitdistrplus::fitdist(x, "exp"))
relloss_fit_wbl = lapply(relloss_list, function(x) fitdistrplus::fitdist(x, "weibull", lower = c(0, 0)))
relloss_fit_gam = lapply(relloss_list, function(x) fitdistrplus::fitdist(x, "gamma"))
relloss_fit_lgn = lapply(relloss_list, function(x) fitdistrplus::fitdist(x, "lnorm"))
relloss_fit_bet = lapply(relloss_list, function(x) fitdistrplus::fitdist(x, "beta"))
relloss_fit_nor = lapply(relloss_list, function(x) MASS::fitdistr(x, "normal"))

absloss_fit_nor = lapply(absloss_list, function(x) fitdistrplus::fitdist(x, "norm"))
#gamma and weibull are also able to fit the data properly this time
#since values are bound between 0 and 1, beta distribution is also tested

additionality_fit_nor = MASS::fitdistr(additionality_series, "normal")

##generate random numbers ----
relloss_rand_exp = lapply(relloss_fit_exp, function(x) rexp(10000, x$estimate))
relloss_rand_wbl = lapply(relloss_fit_wbl, function(x) rweibull(10000, x$estimate[1], x$estimate[2]))
relloss_rand_gam = lapply(relloss_fit_gam, function(x) rgamma(10000, shape = x$estimate[1], rate = x$estimate[2]))
relloss_rand_lgn = lapply(relloss_fit_lgn, function(x) rlnorm(10000, x$estimate[1], x$estimate[2]))
relloss_rand_bet = lapply(relloss_fit_bet, function(x) rbeta(10000, x$estimate[1], x$estimate[2]))
relloss_rand_nor = lapply(relloss_fit_nor, function(x) rnorm(10000, x$estimate[1], x$estimate[2]))
relloss_rand_nor = lapply(relloss_rand_nor, function(x) x[which(x > 0)])

absloss_rand_nor = lapply(absloss_fit_nor, function(x) rnorm(10000, x$estimate[1], x$estimate[2]))
absloss_rand_nor = lapply(absloss_rand_nor, function(x) x[which(x > 0)])
#exponential, gamma and lognormal results in values larger than 1
#normal distribution needs to be truncated at zero

additionality_rand_nor = rnorm(10000, additionality_fit_nor$estimate[1], additionality_fit_nor$estimate[2])
additionality_rand_nor = additionality_rand_nor[which(additionality_rand_nor > 0)]


##examine with p-p plots ----
relloss_cdf_emp = lapply(relloss_list, ecdf)
relloss_cdf_fit_exp = lapply(relloss_rand_exp, ecdf)
relloss_cdf_fit_wbl = lapply(relloss_rand_wbl, ecdf)
relloss_cdf_fit_gam = lapply(relloss_rand_gam, ecdf)
relloss_cdf_fit_lgn = lapply(relloss_rand_lgn, ecdf)
relloss_cdf_fit_bet = lapply(relloss_rand_bet, ecdf)
relloss_cdf_fit_nor = lapply(relloss_rand_nor, ecdf)

absloss_cdf_emp = lapply(absloss_list, ecdf)
absloss_cdf_fit_nor = lapply(absloss_rand_nor, ecdf)

additionality_cdf_emp = ecdf(additionality_series)
additionality_cdf_fit_nor = ecdf(additionality_rand_nor)


plotPP = function(func1, func2, dat, title, wd = 1200, ht = 600, filename){
  png(paste0(file_path, site_label, filename), width = wd, height = ht)
  par(mfrow = c(1, 2), mar = c(6, 6, 6, 2) + 0.1)
  mapply(function(func1, func2, dat, title) {
    plot(func1(sort(dat)), func2(sort(dat)),
         main = title, xlab = "Theoretical cumulative probabilities", ylab = "Sample cumulative probabilities",
         cex.main = 3, cex.lab = 2.5, cex.axis = 2)
    lines(c(0, 1), c(0, 1))},
    func1 = func1, func2 = func2, dat = dat, title = title)
  dev.off()
}

plotPP(func1 = relloss_cdf_fit_exp, func2 = relloss_cdf_emp, dat = relloss_list,
       title = relloss_type, filename = "_relloss_fig3a_fit_pp_exp.png")

plotPP(func1 = relloss_cdf_fit_wbl, func2 = relloss_cdf_emp, dat = relloss_list,
       title = relloss_type, filename = "_relloss_fig3b_fit_pp_wbl.png")

plotPP(func1 = relloss_cdf_fit_gam, func2 = relloss_cdf_emp, dat = relloss_list,
       title = relloss_type, filename = "_relloss_fig3c_fit_pp_gam.png")

plotPP(func1 = relloss_cdf_fit_lgn, func2 = relloss_cdf_emp, dat = relloss_list,
       title = relloss_type, filename = "_relloss_fig3d_fit_pp_lnorm.png")

plotPP(func1 = relloss_cdf_fit_bet, func2 = relloss_cdf_emp, dat = relloss_list,
       title = relloss_type, filename = "_relloss_fig3e_fit_pp_beta.png")

plotPP(func1 = relloss_cdf_fit_nor, func2 = relloss_cdf_emp, dat = relloss_list,
       title = relloss_type, filename = "_relloss_fig3f_fit_pp_norm.png")

plotPP(func1 = absloss_cdf_fit_nor, func2 = absloss_cdf_emp, dat = absloss_list,
       title = absloss_type, filename = "_absloss_fig3f_fit_pp_norm.png")


png(paste0(file_path, site_label, "_additionality_pp_norm.png"), width = 600, height = 600)
par(mar = c(6, 6, 6, 2) + 0.1)
plot(additionality_cdf_fit_nor(sort(additionality_series)), additionality_cdf_emp(sort(additionality_series)),
     main = "Additionality", xlab = "Theoretical cumulative probabilities", ylab = "Sample cumulative probabilities",
     cex.main = 3, cex.lab = 2.5, cex.axis = 2)
lines(c(0, 1), c(0, 1))
dev.off()

##examine with q-q plots ----
plotQQ = function(x, y, distr, title, wd = 1200, ht = 600, filename){
  png(paste0(file_path, site_label, filename), width = wd, height = ht)
  par(mfrow = c(1, 2), mar = c(6, 6, 6, 2) + 0.1, cex.main = 2.5, cex.lab = 2.25, cex.axis = 2)
  mapply(function(x, y, distr, title) {
    params = as.list(y$estimate)
    do.call(car::qqPlot,
            c(list(x, distribution = distr, main = title,
                   xlab = "Theoretical quantiles", ylab = "Sample quantiles"), params))},
    x = x, y = y, distr = distr, title = title)
  dev.off()
}

plotQQ(x = relloss_list, y = relloss_fit_exp, distr = "exp",
       title = relloss_type, filename = "_relloss_fig4a_fit_qq_exp.png")

plotQQ(x = relloss_list, y = relloss_fit_wbl, distr = "weibull",
       title = relloss_type, filename = "_relloss_fig4b_fit_qq_wbl.png")

plotQQ(x = relloss_list, y = relloss_fit_gam, distr = "gamma",
       title = relloss_type, filename = "_relloss_fig4c_fit_qq_gam.png")

plotQQ(x = relloss_list, y = relloss_fit_lgn, distr = "lnorm",
       title = relloss_type, filename = "_relloss_fig4d_fit_qq_lnorm.png")

plotQQ(x = relloss_list, y = relloss_fit_bet, distr = "beta",
       title = relloss_type, filename = "_relloss_fig4e_fit_qq_beta.png")

plotQQ(x = relloss_list, y = relloss_fit_nor, distr = "norm",
       title = relloss_type, filename = "_relloss_fig4f_fit_qq_norm.png")

plotQQ(x = absloss_list, y = absloss_fit_nor, distr = "norm",
       title = absloss_type, filename = "_absloss_fig4f_fit_qq_norm.png")


png(paste0(file_path, site_label, "_additionality_fit_qq_norm.png"), width = 600, height = 600)
par(mfrow = c(1, 2), mar = c(6, 6, 6, 2) + 0.1, cex.main = 2.5, cex.lab = 2.25, cex.axis = 2)
car::qqPlot(additionality_series, distribution = "norm", mean = additionality_fit_nor$estimate[1], sd = additionality_fit_nor$estimate[2],
            main = "Additionality", xlab = "Theoretical quantiles", ylab = "Sample quantiles")
dev.off()


##test goodness-of-fit and divergence ----
ChisqGof = function(dat, fit, distr){
  a = mapply(function(x, y, distr) {
    EnvStats::gofTest(x, test = "chisq", distribution = distr, param.list = as.list(y$estimate))
  },
  x = dat, y = fit, distr = distr, SIMPLIFY = F)
  return(a)
}

gammaScale = function(x){
  a = list(estimate = c(as.list(x$estimate), scale = 1 / as.list(x$estimate)$rate))
  a$estimate["rate"] = NULL
  return(a)
}

relloss_chisq_exp = ChisqGof(dat = relloss_list, fit = relloss_fit_exp, distr = "exp")
relloss_chisq_wbl = ChisqGof(dat = relloss_list, fit = relloss_fit_wbl, distr = "weibull")
relloss_chisq_gam = ChisqGof(dat = relloss_list, fit = lapply(relloss_fit_gam, gammaScale), distr = "gamma")
relloss_chisq_lgn = ChisqGof(dat = relloss_list, fit = relloss_fit_lgn, distr = "lnorm")
relloss_chisq_bet = ChisqGof(dat = relloss_list, fit = relloss_fit_bet, distr = "beta")
relloss_chisq_nor = ChisqGof(dat = relloss_list, fit = relloss_fit_nor, distr = "norm")

absloss_chisq_nor = ChisqGof(dat = absloss_list, fit = absloss_fit_nor, distr = "norm")

sapply(relloss_chisq_exp, function(x) list(stat = x$statistic, pv = x$p.value))
sapply(relloss_chisq_wbl, function(x) list(stat = x$statistic, pv = x$p.value))
sapply(relloss_chisq_gam, function(x) list(stat = x$statistic, pv = x$p.value))
sapply(relloss_chisq_lgn, function(x) list(stat = x$statistic, pv = x$p.value))
sapply(relloss_chisq_bet, function(x) list(stat = x$statistic, pv = x$p.value))
sapply(relloss_chisq_nor, function(x) list(stat = x$statistic, pv = x$p.value))

sapply(absloss_chisq_nor, function(x) list(stat = x$statistic, pv = x$p.value))

additionality_chisq_nor = EnvStats::gofTest(additionality_series, test = "chisq", distribution = "norm", param.list = as.list(additionality_fit_nor$estimate))
(list(stat = additionality_chisq_nor$statistic, pv = additionality_chisq_nor$p.value))

#weibull and beta are good for project area in Gola, none of them are good for the rest
relloss_pa_wbl = fitdistrplus::fitdist(relloss_pa * 100, "weibull")
relloss_pa_bet = fitdistrplus::fitdist(relloss_pa * 100, "beta")


##visualise fitted distributions ----
loss_fit_df = rbind(data.frame(type = "C_a", val = rand_exp$ca),
                    data.frame(type = "C_b", val = rand_exp$cb),
                    data.frame(type = "P_a", val = rand_exp$pa),
                    data.frame(type = "P_b", val = rand_exp$pb))

ggplot(data = loss_fit_df, aes(val)) +
  geom_freqpoly(aes(col = type), lwd = 1, position = "dodge", binwidth = 20000, boundary = 0) +
  scale_color_manual(values = c("red", "pink", "blue", "lightblue"),
                     labels = c("Post-t0 counterfactual", "Pre-t0 counterfactual",
                                "Post-t0 project", "Pre-t0 project")) +
  scale_x_continuous(name = "Annual C loss (Mg CO2e)", limits = c(0, 1600000)) +
  scale_y_continuous(name = "Count") +
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18))
ggsave(paste0(file_path, site_label, "_fig2_distr_fit.png"), width = 15, height = 20, unit = "cm")
