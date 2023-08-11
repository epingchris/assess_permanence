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


# 1a. Generate data for new sites ----
rmarkdown::render("/Users/E-Ping Rau/OneDrive - University of Cambridge/4C_evaluations/R/Reports/eval_template/evaluation_epingrau.Rmd", clean = FALSE)
#t0, eval_end, project_end
tmin5 = t0 - 5
years = 1990:eval_end
site_label = site
#eval_classes = c(1, 2, 3, 4) #only evaluate classes 1-4

# 1b. Load existing data for sites ----
site = "WLT_VNCC_KNT" #Gola_country, WLT_VNCC_KNT, CIF_Alto_Mayo, VCS_1396, VCS_934
load(file = paste0(file_path, site, ".Rdata"))

# 2. Calculate carbon flux ----
makeFlux = function(project_series, leakage_series){
  stock_series = aggregate(class_co2e ~ treatment + year, subset(project_series, class %in% eval_classes), FUN = sum)
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
  stock_wide = as.data.frame(pivot_wider(stock_series, id_cols = year, names_from = "treatment", values_from = "class_co2e"))
  return(list(stock = stock_wide, flux = flux_series))
}

flux_series_sim = mapply(function(x, y) makeFlux(project_series = x, leakage_series = y)$flux,
                         x = agb_series_project_sim,
                         y = vector("list", length = length(agb_series_project_sim)),
                         SIMPLIFY = F)

stock_series_sim = mapply(function(x, y) makeFlux(project_series = x, leakage_series = y)$stock,
                         x = agb_series_project_sim,
                         y = vector("list", length = length(agb_series_project_sim)),
                         SIMPLIFY = F)

flux_series_proj_rel = mapply(function(x, y) {
  data.frame(flux_t1_treatment = x$treatment_proj,
             stock_t0_treatment = y[-nrow(y), ]$treatment,
             flux_t1_control = x$control_proj,
             stock_t0_control = y[-nrow(y), ]$control) %>%
    mutate(rel_flux_treatment = flux_t1_treatment / stock_t0_treatment,
           rel_flux_control = flux_t1_control / stock_t0_control)
}, x = flux_series_sim, y = stock_series_sim, SIMPLIFY = F)

summariseSeries = function(in_list, sel_col, flux = T){
  df = as.data.frame(sapply(in_list, function(x) x[, sel_col])) %>%
    `colnames<-`(paste0("V", 1:20)) %>%
    rowwise() %>%
    mutate(mean = mean(V1:V20)) %>%
    ungroup()
  if(flux) {
    df = df %>%
      mutate(years = years[-1], var = sel_col)
  } else {
    df = df %>%
      mutate(years = years, var = sel_col)
  }
  df = df %>%
    pivot_longer(V1:mean, names_to = "series", values_to = "val") %>%
    return(df)
}

summ_data = rbind(summariseSeries(flux_series_sim, "treatment_proj"),
                  summariseSeries(flux_series_sim, "control_proj"),
                  summariseSeries(flux_series_sim, "additionality"))

summ_data_stock = rbind(summariseSeries(stock_series_sim, "control", flux = F),
                  summariseSeries(stock_series_sim, "treatment", flux = F))


summ_data_rel = rbind(summariseSeries(flux_series_proj_rel, "rel_flux_treatment"),
                      summariseSeries(flux_series_proj_rel, "rel_flux_control")) %>%
  mutate(val = val * -1)

#plot relative carbon loss time series
#lim_to_use = c(floor(min(summ_data_rel$val) * 1000) / 10, ceiling(max(summ_data_rel$val) * 1000) / 10)
lim_to_use = c(-0.2, 1.4) #for comparison between Gola and KNT
ggplot(data = summ_data_rel, aes(col = var)) +
  geom_line(aes(x = years, y = val * 100, lwd = series, alpha = series), show.legend = F) +
  scale_linewidth_manual("", values = c(1, rep(0.5, 20)), labels = c("Mean", "Simulated", rep(NULL, 19))) +
  scale_alpha_manual("", values = c(1, rep(0.1, 20)), labels = c("Mean", "Simulated", rep(NULL, 19))) +
  ggnewscale::new_scale("lwd") +
  ggnewscale::new_scale("alpha") +
  geom_line(data = subset(summ_data_rel, series %in% c("mean", "V1")), aes(x = years, y = val * 100, lwd = series, alpha = series)) +
  geom_vline(aes(xintercept = t0)) + 
  scale_color_manual("", values = c("red", "black"), labels = c("Counterfactual", "Project")) +
  scale_linewidth_manual("", values = c(1, 0.5), labels = c("Mean", "Simulated")) +
  scale_alpha_manual("", values = c(1, 0.1), labels = c("Mean", "Simulated")) +
  scale_x_continuous(name = "Year", breaks = seq(1990, 2021, by = 5)) + 
  scale_y_continuous(name = "Relative C loss (%)", limits = lim_to_use) + 
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16))
ggsave(paste0(file_path, site_label, "_relloss_fig1_time_series.png"), width = 15, height = 20, unit = "cm")


# 3. Distribution fitting ----

##visualise absolute losses of each subset ----
absloss_pa = -subset(summ_data, var == "treatment_proj" & years >= t0 & val < 0)$val
absloss_ca = -subset(summ_data, var == "control_proj" & years >= t0 & val < 0)$val

absloss_df = rbind(data.frame(type = "C_a", val = absloss_ca),
                   data.frame(type = "P_a", val = absloss_pa))

ggplot(data = absloss_df, aes(val)) +
  geom_freqpoly(aes(col = type), lwd = 1, binwidth = 10000, boundary = 0) +
  scale_color_manual(values = c("red", "black"), labels = c("Counterfactual", "Project")) +
  scale_x_continuous(name = "C loss (Mg CO2e)", limits = c(0, 4e+05)) +
  scale_y_continuous(name = "Count") +
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18))
ggsave(paste0(file_path, site_label, "_absloss_fig2_distr_emp.png"), width = 15, height = 20, unit = "cm")

#skewness and kurtosis

##visualise relative losses of each subset ----
relloss_pa = subset(summ_data_rel, var == "rel_flux_treatment" & years >= t0 & val > 0)$val
relloss_ca = subset(summ_data_rel, var == "rel_flux_control" & years >= t0 & val > 0)$val
relloss_pb = subset(summ_data_rel, var == "rel_flux_treatment" & years < t0 & val > 0)$val
relloss_cb = subset(summ_data_rel, var == "rel_flux_control" & years < t0 & val > 0)$val

Epsilon = function(var_sel, post = T){
  if(post){
    n_neg = length(subset(summ_data_rel, var == var_sel & years >= t0 & val > 0)$val)
    n_tot = length(subset(summ_data_rel, var == var_sel & years >= t0)$val)
  } else {
    n_neg = length(subset(summ_data_rel, var == var_sel & years < t0 & val > 0)$val)
    n_tot = length(subset(summ_data_rel, var == var_sel & years < t0)$val)
  }
  return(n_neg / n_tot)
}

(epsilon_pa = Epsilon("rel_flux_treatment"))
(epsilon_ca = Epsilon("rel_flux_control"))

relloss_df = rbind(data.frame(type = "C_a", val = relloss_ca),
                   data.frame(type = "P_a", val = relloss_pa)) %>%
  mutate(val = val * 100) #plot in percentage

ggplot(data = relloss_df, aes(val)) +
  geom_freqpoly(aes(col = type), lwd = 1, binwidth = 0.05, boundary = 0) +
  scale_color_manual(values = c("red", "black"), labels = c("Counterfactual", "Project")) +
  scale_x_continuous(name = "Relative C loss (%)", limits = c(0, 1.5)) +
  scale_y_continuous(name = "Count") +
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18))
ggsave(paste0(file_path, site_label, "_relloss_fig2_distr_emp.png"), width = 15, height = 20, unit = "cm")

## visualise empirical additionality distribution ----
flux_pa = subset(summ_data, var == "treatment_proj" & years >= t0)$val
flux_ca = subset(summ_data, var == "control_proj" & years >= t0)$val
additionality_series = flux_pa - flux_ca

png(paste0(file_path, site_label, "_additionality_distribution.png"), width = 600, height = 600)
hist(additionality_series, breaks = seq(-4e+05, 4e+05, by = 10000), main = ifelse(site =="Gola", "Gola", "KNT"), xlab = "Additionality (Mg CO2e)", cex.main = 2, cex.lab = 1.5)
dev.off()


##fit subsets to different distribution families ----
relloss_list_full = list(pb = relloss_pb, pa = relloss_pa, cb = relloss_cb, ca = relloss_ca)
relloss_type_full = c("Pre-t0 in project", "Post-t0 in project", "Pre-t0 in counterfactual", "Post-t0 in counterfactual")
relloss_list = list(pa = relloss_pa, ca = relloss_ca)
relloss_type = c("Project", "Counterfactual")

absloss_list = list(pa = absloss_pa, ca = absloss_ca)
absloss_type = c("Project", "Counterfactual")


relloss_fit_exp = lapply(relloss_list, function(x) fitdistrplus::fitdist(x, "exp"))
relloss_fit_wbl = lapply(relloss_list, function(x) fitdistrplus::fitdist(x, "weibull", lower = c(0, 0)))
relloss_fit_gam = lapply(relloss_list, function(x) fitdistrplus::fitdist(x, "gamma"))
relloss_fit_lgn = lapply(relloss_list, function(x) fitdistrplus::fitdist(x, "lnorm"))
relloss_fit_bet = lapply(relloss_list, function(x) fitdistrplus::fitdist(x, "beta"))
relloss_fit_nor = lapply(relloss_list, function(x) MASS::fitdistr(x, "normal"))

absloss_fit_nor = lapply(absloss_list, function(x) fitdistrplus::fitdist(x, "norm"))
#gamma and weibull are also able to fit the data properly this time
#since values are bound between 0 and 1, beta distribution is also tested

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
