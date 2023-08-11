library(ggplot2)
library(ggnewscale)
library(tidyverse)
library(magrittr)
library(MASS)
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
site = "Gola_country" #Gola_country, WLT_VNCC_KNT, CIF_Alto_Mayo, VCS_1396, VCS_934
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
  stock_wide = pivot_wider(stock_series, id_cols = year, names_from = "treatment", values_from = "class_co2e")
  return(list(stock = stock_wide, flux = flux_series))
}

flux_series_sim = mapply(makeFlux, project_series = agb_series_project_sim, leakage_series = vector("list", length = length(agb_series_project_sim)),
                         SIMPLIFY = F)

summariseSeries = function(in_list, sel_col){
  df = as.data.frame(sapply(in_list, function(x) x[, sel_col])) %>%
    `colnames<-`(paste0("V", 1:20)) %>%
    rowwise() %>%
    mutate(mean = mean(V1:V20)) %>%
    ungroup() %>%
    mutate(years = years[-1], var = sel_col) %>%
    pivot_longer(V1:mean, names_to = "series", values_to = "val") %>%
  return(df)
}

summ_data = rbind(summariseSeries(flux_series_sim, "treatment_proj"),
                  summariseSeries(flux_series_sim, "control_proj"),
                  summariseSeries(flux_series_sim, "additionality"))

flux_series = pivot_wider(subset(summ_data, series == "mean"), names_from = "var", values_from = "val") %>%
  dplyr::select(-one_of(c("series", "treatment_leak", "control_leak", "leakage", "net_additionality")))

write.table(flux_series, file = paste0(file_path, site_label, "_flux_series.csv"), sep = ",", row.names = F)


#plot carbon flux time series
proj_dat = subset(summ_data, var %in% c("treatment_proj", "control_proj"))
#lim_to_use = c(floor(min(proj_dat$val) / (10 ^ 5)) * (10 ^ 5), 0)
lim_to_use = c(-4e+05, 25000) #for comparison between Gola and KNT
ggplot(data = proj_dat, aes(col = var)) +
  geom_line(aes(x = years, y = val, lwd = series, alpha = series), show.legend = F) +
  scale_linewidth_manual("", values = c(1, rep(0.5, 20)), labels = c("Mean", "Simulated", rep(NULL, 19))) +
  scale_alpha_manual("", values = c(1, rep(0.1, 20)), labels = c("Mean", "Simulated", rep(NULL, 19))) +
  ggnewscale::new_scale("lwd") +
  ggnewscale::new_scale("alpha") +
  geom_line(data = subset(proj_dat, series %in% c("mean", "V1")), aes(x = years, y = val, lwd = series, alpha = series)) +
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
ggsave(paste0(file_path, site_label, "_fig1_time_series.png"), width = 15, height = 20, unit = "cm")


# 3. Distribution fitting ----

## section into subsets ----
loss_pb = -subset(summ_data, var == "treatment_proj" & years < t0 & val < 0)$val
loss_pa = -subset(summ_data, var == "treatment_proj" & years >= t0 & val < 0)$val
loss_ca = -subset(summ_data, var == "control_proj" & years >= t0 & val < 0)$val
loss_cb = -subset(summ_data, var == "control_proj" & years < t0 & val < 0)$val
loss_list = list(pb = loss_pb, pa = loss_pa, cb = loss_cb, ca = loss_ca)
loss_type = c("Pre-t0 in project", "Post-t0 in project", "Pre-t0 in counterfactual", "Post-t0 in counterfactual")

##visualise subsets ----
loss_df = rbind(data.frame(type = "C_a", val = loss_ca),
                data.frame(type = "C_b", val = loss_cb),
                data.frame(type = "P_a", val = loss_pa),
                data.frame(type = "P_b", val = loss_pb))

ggplot(data = loss_df, aes(val)) +
  geom_freqpoly(aes(col = type), lwd = 1, position = "dodge", binwidth = 20000, boundary = 0) +
  scale_color_manual(values = c("red", "pink", "blue", "lightblue"),
                     labels = c("Post-t0 counterfactual", "Pre-t0 counterfactual",
                                "Post-t0 project", "Pre-t0 project")) +
  scale_x_continuous(name = "Annual C loss (Mg CO2e)", limits = c(0, 360000)) +
  scale_y_continuous(name = "Count") +
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18))
ggsave(paste0(file_path, site_label, "_fig2_distr_emp.png"), width = 15, height = 20, unit = "cm")

##fit subsets to different distribution families ----

fit_exp = lapply(loss_list, function(x) MASS::fitdistr(x, "exponential"))
#fit_wbl = lapply(loss_list, function(x) MASS::fitdistr(x, "weibull"))
#fit_gam = lapply(loss_list, function(x) MASS::fitdistr(x, "gamma"))
fit_lgn = lapply(loss_list, function(x) MASS::fitdistr(x, "lognormal"))
#gamma: inappropriate scale
#weibull: NaN in parameter estimation for some (flux_ca), suggesting it is not appropriate either

##generate random numbers ----
rand_exp = lapply(fit_exp, function(x) rexp(10000, x$estimate))
rand_lgn = lapply(fit_lgn, function(x) rlnorm(10000, x$estimate[1], x$estimate[2]))

##examine with p-p plots ----
sizeppqq = 1200
cdf_emp = lapply(loss_list, ecdf)
cdf_fit_exp = lapply(rand_exp, ecdf)
cdf_fit_lgn = lapply(rand_lgn, ecdf)

plotPP = function(func1, func2, dat, title){
  plot(func1(sort(dat)), func2(sort(dat)), main = title,
       xlab = "Theoretical quantiles", ylab = "Sample quantiles")
  lines(c(0, 1), c(0, 1))
}

png(paste0(file_path, site_label, "_fig2_fit_pp_exp.png"), width = sizeppqq, height = sizeppqq)
par(mfrow = c(2, 2))
mapply(plotPP, func1 = cdf_fit_exp, func2 = cdf_emp, dat = loss_list, title = loss_type)
dev.off()

png(paste0(file_path, site_label, "_fig2_fit_pp_lnorm.png"), width = sizeppqq, height = sizeppqq)
par(mfrow = c(2, 2))
mapply(plotPP, func1 = cdf_fit_lgn, func2 = cdf_emp, dat = loss_list, title = loss_type)
dev.off()

##examine with q-q plots ----
png(paste0(file_path, site_label, "_fig2_fit_qq_exp.png"), width = sizeppqq, height = sizeppqq)
par(mfrow = c(2, 2))
mapply(function(x, y) car::qqPlot(x, distribution = "exp", rate = y$estimate), x = loss_list, y = fit_exp)
dev.off()

png(paste0(file_path, site_label, "_fig2_fit_qq_lnorm.png"), width = sizeppqq, height = sizeppqq)
par(mfrow = c(2, 2))
mapply(function(x, y) car::qqPlot(x, distribution = "lnorm", meanlog = y$estimate[1], sdlog = y$estimate[2]), x = loss_list, y = fit_lgn)
dev.off()


##test goodness-of-fit and divergence ----
breaks_list = list(pb = c(seq(0, 100000, by = 15000), 10e10),
                   pa = c(seq(0, 120000, by = 15000), 10e10),
                   cb = c(seq(0, 80000, by = 15000), 10e10),
                   ca = c(seq(0, 360000, by = 60000), 10e10))

vec_emp = mapply(function(x, y) hist(x, breaks = y)$count, x = loss_list, y = breaks_list)
vec_exp = mapply(function(x, y, z) hist(rexp(length(x), y$estimate), breaks = z)$count, x = loss_list, y = fit_exp, z = breaks_list)
vec_lgn = mapply(function(x, y, z) hist(rlnorm(length(x), y$estimate[1], y$estimate[2]), breaks = z)$count, x = loss_list, y = fit_lgn, z = breaks_list)

chisq_exp = mapply(function(x, y) chisq.test(x, p = y, rescale.p = T), x = vec_emp, y = vec_exp)
unlist(chisq_exp["p.value", ])
chisq_lgn = mapply(function(x, y) chisq.test(x, p = y, rescale.p = T), x = vec_emp, y = vec_lgn)
unlist(chisq_lgn["p.value", ])
#significantly different from both exponential or lognormal

##test JS divergence ----
JSCalc = function(a, b){
  bin_vec = seq(0, ceiling(max(a, b) / 10 ^ 4) * 10 ^ 4, by = 10000)
  js = compute_js_divergence(a, b, bins = bin_vec)
  return(js)
}

JS_exp = mapply(function(x, y) JSCalc(x, rexp(length(x), y$estimate)), x = loss_list, y = fit_exp)
#[1] 0.016914936 0.018392344 0.008019489 0.114439386
JS_lgn = mapply(function(x, y) JSCalc(x, rlnorm(length(x), y$estimate[1], y$estimate[2])), x = loss_list, y = fit_lgn)
#1] 0.010455531 0.028230478 0.003223241 0.097510087

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


# 3b. Calculate and test JS divergence ----
rexp_p_bfr = rexp(10000, fit_p_bfr$estimate)
rexp_cf_bfr = rexp(10000, fit_cf_bfr$estimate)
rexp_p_aftr = rexp(10000, fit_p_aftr$estimate)
rexp_cf_aftr = rexp(10000, fit_cf_aftr$estimate)

JSTest = function(a, b, alpha = 0.05, n_rep = 10000){
  bin_vec = seq(0, ceiling(max(a, b) / 10 ^ 4) * 10 ^ 4, by = 10000)
  js = compute_js_divergence(a, b, bins = bin_vec)
  pooled = c(a, b)
  js_perm = rep(NA, n_rep)
  for(i in 1:n_rep){
    samp_ind = sample(1:(length(a) + length(b)), length(a))
    samp_a = pooled[samp_ind]
    samp_b = pooled[-samp_ind]
    js_perm[i] = compute_js_divergence(samp_a, samp_b, bins = bin_vec)
  }
  js_crit = quantile(js_perm, 1 - alpha)
  p_val = length(which(js_perm > js)) / n_rep
  if(p_val == 0) {p_val = paste("<", 1 / n_rep)}
  return(list(bin_val = bin_vec, js_obs = js, js_crit = js_crit, p_val = p_val))
}

a0 = Sys.time()
JS_res = list(JS_bfr = JSTest(rexp_p_bfr, rexp_cf_bfr),
              JS_aftr = JSTest(rexp_p_aftr, rexp_cf_aftr),
              JS_cf = JSTest(rexp_cf_bfr, rexp_cf_aftr),
              JS_p = JSTest(rexp_p_bfr, rexp_p_aftr))
a1 = Sys.time()
a1 - a0


#visualise distributions
# ymax = ifelse(site_label %in% c("Gola_country", "WLT_VNCC_KNT"), 4e-05, 2e-05)
# png(paste0(file_path, site_label, "_fig2_distr_comparison.png"), width = 600, height = 800)
# plot(0:1e-05, 0:1e-05, xlim = c(0, 150000), ylim = c(0, ymax), type = "n",
#      xlab = "Annual C loss (Mg CO2e)", ylab = "Prob. density", cex.axis = 2.5, cex.lab = 2.5)
# lines(density(rexp_p_bfr), col = "lightblue", lwd = 3)
# lines(density(rexp_cf_bfr), col = "pink", lwd = 3)
# lines(density(rexp_p_aftr), col = "blue", lwd = 3)
# lines(density(rexp_cf_aftr), col = "red", lwd = 3)
# #legend(x = "topright", legend = c("Project (before t0)", "Counterfactual (before t0)", "Project (after t0)", "Counterfactual (after t0)"),
# #       fill = c("lightblue", "pink", "blue", "red"))
# dev.off()



# 3c. Test time-invariance of distributions using 10-year intervals ----
years_interval = years[which(years < t0)][-1]
n_window = length(years_interval) - 10 + 1
fit_p_param = rep(NA, n_window)
fit_cf_param = rep(NA, n_window)
rexplist_p = vector("list", n_window)
rexplist_cf = vector("list", n_window)
js_div = rep(NA, n_window)

for(i in 1:n_window){
  p_dat = subset(summ_data, var == "treatment_proj" & years >= years_interval[i] & years <= years_interval[i + 9] & val < 0)$val
  fit_p = MASS::fitdistr(-p_dat, "exponential")
  fit_p_param[i] = fit_p$estimate
  
  cf_dat = subset(summ_data, var == "control_proj" & years >= years_interval[i] & years <= years_interval[i + 9] & val < 0)$val
  fit_cf = MASS::fitdistr(-cf_dat, "exponential")
  fit_cf_param[i] = fit_cf$estimate
  
  rexplist_p[[i]] = rexp(1000, fit_p$estimate)
  rexplist_cf[[i]] = rexp(1000, fit_cf$estimate)
  js_div[i] = compute_js_divergence(rexplist_p[[i]], rexplist_cf[[i]], bins = seq(0, ceiling(max(rexplist_p[[i]], rexplist_cf[[i]]) / 10 ^ 4) * 10 ^ 4, by = 10000))
}


#plot density functions through time
sat_val = seq(0.2, 1, len = n_window)
ymax = ifelse(site_label %in% c("Gola_country", "WLT_VNCC_KNT"), 5e-05, 3e-05)

png(paste0(file_path, site_label, "_fig3a_proj_distr_through_time.png"), width = 600, height = 800)
plot(0:1e-05, 0:1e-05, xlim = c(0, 150000), ylim = c(0, ymax), type = "n",
     xlab = "Annual C loss (Mg CO2e)", ylab = "Prob. density", cex.axis = 2, cex.lab = 2)
for(i in 1:(length(years_interval) - 10 + 1)){
  lines(density(rexplist_p[[i]]), col = hsv(2/3, sat_val[i], 1), lwd = 3)
}
dev.off()

png(paste0(file_path, site_label, "_fig3b_cf_distr_through_time.png"), width = 600, height = 800)
plot(0:1e-05, 0:1e-05, xlim = c(0, 150000), ylim = c(0, ymax), type = "n",
     xlab = "Annual C loss (Mg CO2e)", ylab = "Prob. density", cex.axis = 2, cex.lab = 2)
for(i in 1:(length(years_interval) - 10 + 1)){
  lines(density(rexplist_cf[[i]]), col = hsv(0, sat_val[i], 1), lwd = 3)
}
dev.off()

#use linear model to test temporal trend
lm_p_param = lm(param ~ year, data = data.frame(year = 1:n_window,
                                                param = 1/fit_p_param))
lm_cf_param = lm(param ~ year, data = data.frame(year = 1:n_window,
                                                 param = 1/fit_cf_param))
lm_js_div = lm(param ~ year, data = data.frame(year = 1:n_window,
                                               param = js_div))
summary(lm_p_param)
summary(lm_cf_param)
summary(lm_js_div)

#plot rate parameter through time
png(paste0(file_path, site_label, "_fig4a_rate_param_through_time.png"), width = 600, height = 800)
plot(years_interval[1:n_window], years_interval[1:n_window], xlim = c(years_interval[1], years_interval[n_window]),
     ylim = c(floor(min(1/fit_p_param, 1/fit_cf_param) * 10^5), ceiling(max(1/fit_p_param, 1/fit_cf_param) * 10^5))  / (10^5),
     type = "n", xlab = "Time interval (t0)", ylab = "Distribution mean", cex.axis = 2.5, cex.lab = 2.5)
lines(years_interval[1:n_window], 1/fit_p_param, col = "blue")
lines(years_interval[1:n_window], 1/fit_cf_param, col = "red")
#legend(x = "topright", legend = c("Project", "Counterfactual"), fill = c("blue", "red"), cex = 2.5)
# text(x = years_interval[6], y = ceiling(max(fit_p_param, fit_cf_param) * 10^5)  / (10^5) * 0.9,
#      labels = paste0("Project: p = ", round(summary(lm_p_param)$coefficients[2, 4], 7), "\nR2 = ", round(summary(lm_p_param)$r.squared, 2),
#                      "\nCounterfactual: p = ", round(summary(lm_cf_param)$coefficients[2, 4], 7), "\nR2 = ", round(summary(lm_cf_param)$r.squared, 2)),
#      cex = 2)
dev.off()

#plot JS divergence through time
png(paste0(file_path, site_label, "_fig4b_divergence_through_time.png"), width = 600, height = 800)
plot(years_interval[1:n_window], js_div, xlim = c(years_interval[1], years_interval[n_window]),
     ylim = c(0, 0.05), type = "l", xlab = "Time interval (t0)",
     ylab = "JS divergence (project vs. counterfactual)", cex.axis = 2.5, cex.lab = 2.5)
# text(x = years_interval[5], y = 0.04,
#      labels = paste0("p = ", round(summary(lm_js_div)$coefficients[2, 4], 7), "\nR2 = ", round(summary(lm_js_div)$r.squared, 2)),
#      cex = 2)
dev.off()


# 4. Construct release distribution ----

##total carbon drawdown at t_eval ----
eval_end_additionality = sum(subset(flux_series, years >= t0)$additionality)

##make functions to evaluate permanence on an annual basis ----
ep_calc_ann<-function(drawdown, release, SCC_drawdown = scc, SCC_release = scc, r_discount = 0.03){
  # in future we might want to generalise this function to be able to handle a drawdown
  # schedule too
  
  SCC_drawdown = SCC_drawdown[which(names(SCC_drawdown) == as.character(eval_end))]
  project_value = SCC_drawdown * drawdown
  
  # calculate the total damages caused by release
  release_years = names(release)
  release_year_index = as.integer(release_years) - eval_end
  SCC_release = SCC_release[which(names(SCC_release) %in% release_years)]
  SCC_release_discounted <- SCC_release / ((1 + r_discount) ^ release_year_index)
  damage <- SCC_release_discounted * release
  project_damage <- sum(damage)
  
  ep <- (project_value - project_damage) / project_value
  return(list(value = project_value, damage = project_damage, ep = ep))
}

zero_clear<-function(x){
  y<-x[ min( which ( x != 0 )) : max( which( x != 0 )) ]
  y<-y[y>=0]
  y<-c(y,0)
  return(y)
}

##generate random post-t0 additionality samples and distribution ----
diff_aftr = rexp_cf_aftr - rexp_p_aftr

break_seq = seq(-max(abs(range(diff_aftr))) - 1000, max(abs(range(diff_aftr))) + 1000, len = 21)
png(paste0(file_path, site_label, "_fig5_distr_additionality.png"), width = 600, height = 800)
color_list = rep("grey", length(break_seq))
color_list[break_seq < 0] = "black"
hist(diff_aftr, freq = F, breaks = break_seq, col = color_list,
     xlab = "Additionality during project", ylab = "Prob. density", main = "", cex.lab = 2)
dev.off()

#low-risk: carbon release only after project, from the distribution of carbon drawdown during project
#high-risk: carbon release after evaluation, from the distribution of carbon drawdown during project
hrisk = as.numeric(gsub("< ", "", JS_res$JS_aftr$p_val)) > 0.05

if(hrisk) {
  fit_release = MASS::fitdistr(-diff_aftr[which(diff_aftr < 0)], "exponential")
} else {
  fit_release = MASS::fitdistr(diff_aftr[which(diff_aftr > 0)], "exponential")
}


##evaluate permanence (at t_eval) under simulated release schedule ----
release_schedule = vector("list", 100)
len_release = rep(NA, 100)
ep_vec = rep(NA, 100)
for(i in 1:100){
  release_amount = rexp(100, fit_release$estimate) #assumes C drawdown before t_eval = C release after project (low risk) / t_eval (high risk)

  release_cumul = cumsum(release_amount)
  if(release_cumul[1] >= eval_end_additionality){
    release_simul = c(eval_end_additionality, 0)
  } else {
    release_simul = zero_clear(c(eval_end_additionality, eval_end_additionality - release_cumul))
  }
  
  len_release[i] = length(release_simul) - 1
  release_schedule[[i]] = release_simul
  if(hrisk) {
    release_extended = release_simul
  } else {
    release_extended = c(rep(eval_end_additionality, project_end - eval_end), release_simul)
  }
  names(release_extended) = eval_end:(eval_end + length(release_extended) - 1)
  ep_vec[i] = ep_calc_ann(drawdown = eval_end_additionality, release = -diff(release_extended))$ep
}

png(paste0(file_path, site_label, "_fig6a_release_time_", ifelse(hrisk, "hrisk", "lrisk"), ".png"), width = 600, height = 800)
hist(len_release, breaks = 1:45, cex.axis = 2, cex.lab = 2,
     xlab = "Number of years", main = "")
dev.off()

png(paste0(file_path, site_label, "_fig6b_ep_", ifelse(hrisk, "hrisk", "lrisk"), ".png"), width = 600, height = 800)
hist(ep_vec, breaks = seq(0, 0.5, by = 0.025), cex.axis = 2, cex.lab = 2,
     xlab = "eP", main = "")
dev.off()


##evaluate permanence (at t_eval) under fixed release schedule (currently implemented) ----
release_amount_fixed = rep(mean(cf_fit_tmin5$cf_release), 100)
release_cumul_fixed = cumsum(release_amount_fixed)
if(release_cumul_fixed[1] >= eval_end_additionality){
  release_fixed = c(eval_end_additionality, 0)
} else {
  release_fixed = zero_clear(c(eval_end_additionality, eval_end_additionality - release_cumul_fixed))
}

len_release_fixed = length(release_fixed) - 1
release_extended_fixed = c(rep(eval_end_additionality, project_end - eval_end), release_fixed)
names(release_extended_fixed) = eval_end:(eval_end + length(release_extended_fixed) - 1)
ep_fixed = ep_calc_ann(drawdown = eval_end_additionality, release = -diff(release_extended_fixed))$ep

##save data ----
save(scc, t0, tmin5, project_end, eval_end, years, eval_classes, site, site_label,
     agb_series_project_sim, flux_series_sim, summ_data, flux_series,
     fit_p_bfr, fit_cf_bfr, fit_p_aftr, fit_cf_aftr, rexp_p_bfr, rexp_cf_bfr, rexp_p_aftr, rexp_cf_aftr,
     JS_res, rexplist_p, rexplist_cf, years_interval, n_window,
     fit_p_param, fit_cf_param, js_div, lm_p_param, lm_cf_param, lm_js_div,
     additionality, fit_additionality, release, fit_release, hrisk, eval_end_additionality, len_release, release_schedule, ep_vec, file = paste0(file_path, site_label, ".Rdata"))


# 5. Create step-wise release schedule and calculate credibility and permanence ----
ep_calc_stepwise<-function(eval_year, drawdown, drawdown_year, release, SCC_drawdown = scc, SCC_release = scc, r_discount = 0.03){
  # in future we might want to generalise this function to be able to handle a drawdown
  # schedule too
  
  SCC_drawdown = SCC_drawdown[which(names(SCC_drawdown) == as.character(drawdown_year))]
  project_value = SCC_drawdown * drawdown
  
  # calculate the total damages caused by release
  release_years = names(release)
  release_year_index = as.integer(release_years) - eval_year
  SCC_release = SCC_release[which(names(SCC_release) %in% release_years)]
  SCC_release_discounted <- SCC_release / ((1 + r_discount) ^ release_year_index)
  damage <- SCC_release_discounted * release
  project_damage <- sum(damage)
  
  ep = ifelse(project_value != 0, (project_value - project_damage) / project_value, NA)
  return(list(value = project_value, damage = project_damage, ep = ep))
}

cred_lrisk_vec = rep(NA, 100)
cred_hrisk_vec = rep(NA, 100)
ep_sim_vec = rep(NA, 100)
rhat_list = vector("list", 100)
flux_series_sim_list = vector("list", 100)

for(i in 1:100){
  #generate simulated project additionality for another 100 years based on ex post carbon flux distribution
  p_aftr_sim = rexp(100, fit_p_aftr$estimate)
  cf_aftr_sim = rexp(100, fit_cf_aftr$estimate)
  additionality_sim = cf_aftr_sim - p_aftr_sim
  
  flux_series_sim = subset(flux_series, years >= t0) %>%
    dplyr::select(years, additionality) %>%
    rbind(data.frame(years = (eval_end + 1):(eval_end + 100),
                     additionality = additionality_sim)) %>%
    mutate(rhat = 0, credit = 0)
  
  rhat_list[[i]] = vector("list", nrow(flux_series_sim))
  value_vec = rep(NA, length(flux_series_sim$years)) %>%
    `names<-`(flux_series_sim$years)
  damage_vec = rep(NA, length(flux_series_sim$years)) %>%
    `names<-`(flux_series_sim$years)
  
  for(t_eval in flux_series_sim$years){
    #calculate credit as additionality + rhat from previous years (previously anticipated release compensates for present release)
    flux_series_sim %<>% 
      mutate(credit = ifelse(years == t_eval, additionality + rhat, credit))

    #calculate rhat for future years from carbon release distribution (anticipated release to compensate for future release)
    a_t = subset(flux_series_sim, years == t_eval)$additionality
    release_amount = rexp(100, fit_release$estimate)
    release_cumul = cumsum(release_amount)
    if(a_t < 0) {
      rhat = 0
    } else {
      if(release_cumul[1] >= a_t){
        release_series = c(a_t, 0)
      } else {
        release_series = zero_clear(c(a_t, a_t - release_cumul))
      }
      rhat = -diff(release_series)
    }
    
    rhat_df = data.frame(rhat = rhat,
                         i = t_eval,
                         j = (t_eval + 1):(t_eval + length(rhat)))
    rhat_list[[i]][[which(flux_series_sim$years %in% t_eval)]] = rhat_df
    
    years_release = which(flux_series_sim$years %in% rhat_df$j)
    if(length(years_release) > 0){
      flux_series_sim[years_release, "rhat"] = flux_series_sim[years_release, "rhat"] + rhat_df$rhat
    }
    
    #calculate value and damage for eP calculation
    if(t_eval >= 2020){
      names(rhat) = rhat_df$j
      ep_res = ep_calc_stepwise(eval_year = t0,
                                drawdown = max(a_t, 0), drawdown_year = t_eval,
                                release = rhat)
      value_vec[which(names(value_vec) == as.character(t_eval))] = ep_res$value
      damage_vec[which(names(damage_vec) == as.character(t_eval))] = ep_res$damage
    }

  }
  ep_sim_vec[i] = (sum(value_vec, na.rm = T) - sum(damage_vec, na.rm = T)) / sum(value_vec, na.rm = T)

  flux_series_sim_list[[i]] = flux_series_sim
  cred_lrisk_vec[i] = 1 - length(which(flux_series_sim$additionality < 0)) / nrow(flux_series_sim)
  cred_hrisk_vec[i] = 1 - length(which(flux_series_sim$credit < 0)) / nrow(flux_series_sim)
}

ggplot(data = flux_series_sim_list[[1]]) +
  geom_line(aes(x = years, y = additionality)) +
  scale_x_continuous(name = "Year", breaks = seq(t0, eval_end + 100, by = 10)) + 
  scale_y_continuous(name = "Additionality (Mg CO2e)") + 
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        legend.title = element_text(),
        axis.title = element_text(size = 16))
ggsave(paste0(file_path, site_label, "_fig7a_sim_additionality.png"), width = 15, height = 10, unit = "cm")


ggplot(data = flux_series_sim_list[[1]]) +
  geom_line(aes(x = years, y = additionality)) +
  geom_line(aes(x = years, y = rhat), color = "red") +
  geom_line(aes(x = years, y = credit), color = "blue") +
  scale_x_continuous(name = "Year", breaks = seq(t0, eval_end + 100, by = 10)) + 
  scale_y_continuous(name = "Additionality (Mg CO2e)") + 
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        legend.title = element_text(),
        axis.title = element_text(size = 16))
ggsave(paste0(file_path, site_label, "_fig7b_sim_credit.png"), width = 15, height = 10, unit = "cm")

##evaluate permanence (at t_eval)  ----
release_schedule = vector("list", 100)
len_release = rep(NA, 100)
ep_vec = rep(NA, 100)
for(i in 1:100){
  release_amount = rexp(100, fit_release$estimate) #assumes C drawdown before t_eval = C release after project (low risk) / t_eval (high risk)
  
  release_cumul = cumsum(release_amount)
  if(release_cumul[1] >= eval_end_additionality){
    release_simul = c(eval_end_additionality, 0)
  } else {
    release_simul = zero_clear(c(eval_end_additionality, eval_end_additionality - release_cumul))
  }
  
  len_release[i] = length(release_simul) - 1
  release_schedule[[i]] = release_simul
  if(hrisk) {
    release_extended = release_simul
  } else {
    release_extended = c(rep(eval_end_additionality, project_end - eval_end), release_simul)
  }
  names(release_extended) = eval_end:(eval_end + length(release_extended) - 1)
  ep_res = ep_calc_stepwise(eval_year = eval_end,
                               drawdown = eval_end_additionality, drawdown_year = ,
                               release = -diff(release_extended))
  value_vec[i] = ep_res$value
  damage_vec[i] = ep_res$damage
}
ep = (sum(value_vec) - sum(damage_vec)) / sum(value_vec)


##save data ----
save(scc, t0, tmin5, project_end, eval_end, years, eval_classes, site, site_label,
     agb_series_project_sim, flux_series_sim, summ_data, flux_series,
     fit_p_aftr, fit_cf_aftr, rexp_p_aftr, rexp_cf_aftr, diff_aftr,
     p_aftr_sim, cf_aftr_sim, additionality_sim,
     flux_series_sim_list, rhat_list, flux_series_sim, cred_lrisk_vec, cred_hrisk_vec, file = paste0(file_path, site_label, "rhat.Rdata"))

file_path = "C:/Users/E-Ping Rau/OneDrive - University of Cambridge/carbon_release_pattern/"
site_label = "WLT_VNCC_KNT"
load(paste0(file_path, site_label, "rhat.Rdata"))


# 6a. Test how the distributions change based on number of years of observation were used to fit them ----
param_p_part = rep(NA, eval_end - t0 + 1)
param_cf_part = rep(NA, eval_end - t0 + 1)
param_release_part = rep(NA, eval_end - t0 + 1)
JS_p_part = vector("list", eval_end - t0 + 1)
JS_cf_part = vector("list", eval_end - t0 + 1)
JS_release_part = vector("list", eval_end - t0 + 1)

release = rexp_p_aftr - rexp_cf_aftr
fit_release = MASS::fitdistr(-release[which(release < 0)], "exponential")
rand_release = rexp(10000, fit_release$estimate)

for(i in 1:(eval_end - t0 + 1)){
  flux_p_part = subset(summ_data, var == "treatment_proj" & years >= t0 & years < t0 + i & val < 0)$val
  fit_p_part = MASS::fitdistr(-flux_p_part, "exponential")
  rand_p_part = rexp(10000, fit_p_part$estimate)
  param_p_part[i] = fit_p_part$estimate
  JS_p_part[[i]] = JSTest(rand_p_part, rexp_p_aftr, n_rep = 5000)
  
  flux_cf_part = subset(summ_data, var == "control_proj" & years >= t0 & years < t0 + i & val < 0)$val
  fit_cf_part = MASS::fitdistr(-flux_cf_part, "exponential")
  rand_cf_part = rexp(10000, fit_cf_part$estimate)
  param_cf_part[i] = fit_cf_part$estimate
  JS_cf_part[[i]] = JSTest(rand_cf_part, rexp_cf_aftr, n_rep = 5000)
  
  release_part = rand_p_part - rand_cf_part
  fit_release_part = MASS::fitdistr(-release_part[which(release_part < 0)], "exponential")
  rand_release_part = rexp(10000, fit_release_part$estimate)
  param_release_part[i] = fit_release_part$estimate
  JS_release_part[[i]] = JSTest(rand_release_part, rand_release, n_rep = 5000)
}

plot(1:(eval_end - t0 + 1), as.numeric(gsub("< ", "", sapply(JS_p_part, function(x) x$p_val))), type = "l", col = hsv(1/3, 1, 1), lwd = 2, ylim = c(0, 1))
abline(h = 0.05)
plot(1:(eval_end - t0 + 1), as.numeric(gsub("< ", "", sapply(JS_cf_part, function(x) x$p_val))), type = "l", col = hsv(2/3, 1, 1), lwd = 2, ylim = c(0, 1))
abline(h = 0.05)
plot(1:(eval_end - t0 + 1), as.numeric(gsub("< ", "", sapply(JS_release_part, function(x) x$p_val))), type = "l", col = hsv(1, 1, 1), lwd = 2, ylim = c(0, 1))
abline(h = 0.05)


# 6b. Test how the distributions change based on number of replicates of observation were used to fit them ----
summ_data_mod = summ_data %<>%
  filter(series != "mean") %>%
  mutate(series = as.numeric(gsub("V", "", series)))

param_p_part = rep(NA, 20)
param_cf_part = rep(NA, 20)
param_release_part = rep(NA, 20)
JS_p_part = vector("list", 20)
JS_cf_part = vector("list", 20)
JS_release_part = vector("list", 20)

for(i in 1:20){
  flux_p = subset(summ_data_mod, var == "treatment_proj" & years >= t0 & series <= i & val < 0)$val
  fit_p = MASS::fitdistr(-flux_p, "exponential")
  rand_p = rexp(10000, fit_p$estimate)
  param_p_part[i] = fit_p$estimate
  JS_p_part[[i]] = JSTest(rand_p_part, rexp_p_aftr, n_rep = 5000)
  
  flux_cf = subset(summ_data_mod, var == "control_proj" & years >= t0 & series <= i & val < 0)$val
  fit_cf = MASS::fitdistr(-flux_cf, "exponential")
  rand_cf = rexp(10000, fit_cf$estimate)
  param_cf_part[i] = fit_cf$estimate
  JS_cf_part[[i]] = JSTest(rand_cf_part, rexp_cf_aftr, n_rep = 5000)
  
  release = rand_p - rand_cf
  fit_release = MASS::fitdistr(-release[which(release < 0)], "exponential")
  rand_release = rexp(10000, fit_release$estimate)
  param_release_part[i] = fit_release$estimate
  JS_release_part[[i]] = JSTest(rand_release_part, rand_release, n_rep = 5000)
}

plot(1:20, as.numeric(gsub("< ", "", sapply(JS_p_part, function(x) x$p_val))), type = "l", col = hsv(1/3, 1, 1), lwd = 2, ylim = c(0, 1))
abline(h = 0.05)
plot(1:20, as.numeric(gsub("< ", "", sapply(JS_cf_part, function(x) x$p_val))), type = "l", col = hsv(2/3, 1, 1), lwd = 2, ylim = c(0, 1))
abline(h = 0.05)
plot(1:20, as.numeric(gsub("< ", "", sapply(JS_release_part, function(x) x$p_val))), type = "l", col = hsv(1, 1, 1), lwd = 2, ylim = c(0, 1))
abline(h = 0.05)

