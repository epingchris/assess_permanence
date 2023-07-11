library(ggplot2)
library(ggnewscale)
library(tidyverse)
library(magrittr)
library(MASS)
library(stat)
source("divergence.R")


# 1a. Generate data for new sites ----
rmarkdown::render("/Users/E-Ping Rau/OneDrive - University of Cambridge/4C_evaluations/R/Reports/eval_template/evaluation_epingrau.Rmd", clean = FALSE)
#t0, eval_end, project_end
tmin5 = t0 - 5
years = 1990:eval_end
site_label = site
file_path = "C:/Users/E-Ping Rau/OneDrive - University of Cambridge/carbon_release_pattern/"
source(paste0(file_path, "divergence.R"))
#eval_end = 2021
#eval_classes = c(1, 2, 3, 4) #only evaluate classes 1-4

# 1b. Load existing data for sites ----
site = "CIF_Alto_Mayo" #Gola_country, WLT_VNCC_KNT, VCS_1396
load(file = paste0(file_path, site, ".Rdata"))
site_label = site

# 2. Calculate carbon flux ----
makeFlux = function(project_series, leakage_series){
  project_aggr = aggregate(class_co2e ~ treatment + year, subset(project_series, class %in% eval_classes), FUN = sum) 
  flux_series = data.frame(year = years[-1],
                           treatment_proj = diff(subset(project_aggr, treatment == "treatment")$class_co2e),
                           control_proj = diff(subset(project_aggr, treatment == "control")$class_co2e)) %>%
    mutate(additionality = treatment_proj - control_proj)
  if(!is.null(leakage_series)){
    leakage_aggr = aggregate(class_co2e ~ treatment + year, subset(leakage_series, class %in% eval_classes), FUN = sum)
    flux_series %<>%
      mutate(treatment_leak = diff(subset(leakage_aggr, treatment == "treatment")$class_co2e),
             control_leak = diff(subset(leakage_aggr, treatment == "control")$class_co2e)) %>%
      mutate(leakage = treatment_leak - control_leak)
  }
  return(flux_series)
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
lim_to_use = c(floor(min(proj_dat$val) / (10 ^ 5)) * (10 ^ 5), 0)
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
  scale_y_continuous(name = "Flux (Mg CO2e)", limits = lim_to_use) + 
  ggtitle("") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.title = element_text())
ggsave(paste0(file_path, site_label, "_fig1_time_series.png"), width = 15, height = 15, unit = "cm")


# 3. Fit C flux to exponential distributions ----
flux_p_bfr = subset(summ_data, var == "treatment_proj" & years < t0 & val < 0)$val #325
fit_p_bfr = MASS::fitdistr(-flux_p_bfr, "exponential")

flux_cf_bfr = subset(summ_data, var == "control_proj" & years < t0 & val < 0)$val #353
fit_cf_bfr = MASS::fitdistr(-flux_cf_bfr, "exponential")

flux_p_aftr = subset(summ_data, var == "treatment_proj" & years >= t0 & val < 0)$val #209
fit_p_aftr = MASS::fitdistr(-flux_p_aftr, "exponential")

flux_cf_aftr = subset(summ_data, var == "control_proj" & years >= t0 & val < 0)$val #210
fit_cf_aftr = MASS::fitdistr(-flux_cf_aftr, "exponential")


#calculate and test JS divergence
rexp_p_bfr = rexp(1000, fit_p_bfr$estimate)
rexp_cf_bfr = rexp(1000, fit_cf_bfr$estimate)
rexp_p_aftr = rexp(1000, fit_p_aftr$estimate)
rexp_cf_aftr = rexp(1000, fit_cf_aftr$estimate)

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
  if(p_val == 0) {p_val = paste0("< ", 1 / n_rep)}
  return(list(bin_val = bin_vec, js_obs = js, js_crit = js_crit, p_val = p_val))
}

a0 = Sys.time()
JS_res = list(JS_bfr = JSTest(rexp_p_bfr, rexp_cf_bfr, n_rep = 1000),
              JS_aftr = JSTest(rexp_p_aftr, rexp_cf_aftr, n_rep = 1000),
              JS_cf = JSTest(rexp_cf_bfr, rexp_cf_aftr, n_rep = 1000),
              JS_p = JSTest(rexp_p_bfr, rexp_p_aftr, n_rep = 1000))
a1 = Sys.time()
a1 - a0


#visualise distributions
png(paste0(file_path, site_label, "_fig2_distr_comparison.png"), width = 600, height = 600)
plot(0:1e-05, 0:1e-05, xlim = c(0, 200000), ylim = c(0, 4e-05), type = "n",
     xlab = "C flux", ylab = "Density")
lines(density(rexp_p_bfr), col = "lightblue", lwd = 3)
lines(density(rexp_cf_bfr), col = "pink", lwd = 3)
lines(density(rexp_p_aftr), col = "blue", lwd = 3)
lines(density(rexp_cf_aftr), col = "red", lwd = 3)
legend(x = "topright", legend = c("Project (before t0)", "Counterfactual (before t0)", "Project (after t0)", "Counterfactual (after t0)"),
       fill = c("lightblue", "pink", "blue", "red"))
dev.off()


#test time-invariance of distributions using 10-year intervals
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

png(paste0(file_path, site_label, "_fig3_proj_distr_through_time.png"), width = 1000, height = 1000)
plot(0:1e-05, 0:1e-05, xlim = c(0, 200000), ylim = c(0, 4e-05), type = "n",
     xlab = "C flux", ylab = "Density")
for(i in 1:(length(years_interval) - 10 + 1)){
  lines(density(rexplist_p[[i]]), col = hsv(2/3, sat_val[i], 1), lwd = 3)
}
dev.off()

png(paste0(file_path, site_label, "_fig3_cf_distr_through_time.png"), width = 1000, height = 1000)
plot(0:1e-05, 0:1e-05, xlim = c(0, 200000), ylim = c(0, 4e-05), type = "n",
     xlab = "C flux", ylab = "Density")
for(i in 1:(length(years_interval) - 10 + 1)){
  lines(density(rexplist_cf[[i]]), col = hsv(0, sat_val[i], 1), lwd = 3)
}
dev.off()

#plot rate parameter through time
png(paste0(file_path, site_label, "_fig4_rate_param_through_time.png"), width = 600, height = 600)
plot(years_interval[1:n_window], years_interval[1:n_window], xlim = c(years_interval[1], years_interval[n_window]),
     ylim = c(floor(min(fit_p_param, fit_cf_param) * 10^5), ceiling(max(fit_p_param, fit_cf_param) * 10^5))  / (10^5),
     type = "n", xlab = "Time interval (t0)", ylab = "Rate parameter")
lines(years_interval[1:n_window], fit_p_param, col = "blue")
lines(years_interval[1:n_window], fit_cf_param, col = "red")
legend(x = "topright", legend = c("Project", "Counterfactual"), fill = c("blue", "red"))
dev.off()

lm_p_param = lm(param ~ year, data = data.frame(year = 1:n_window,
                                   param = fit_p_param))
lm_cf_param = lm(param ~ year, data = data.frame(year = 1:n_window,
                                   param = fit_cf_param))

#plot rate parameter through time
png(paste0(file_path, site_label, "_fig4_divergence_through_time.png"), width = 600, height = 600)
plot(years_interval[1:n_window], js_div, xlim = c(years_interval[1], years_interval[n_window]), ylim = c(0, 0.05), type = "l",
     xlab = "Time interval (t0)", ylab = "JS divergence project vs. counterfactual")
dev.off()

lm_js_div = lm(param ~ year, data = data.frame(year = 1:n_window,
                                   param = js_div))


# Construct release distributions ----
#Approach 1: the right hand side of distribution before project (additionality release due to random drift)
# release_bfr = rexp_p_bfr - rexp_cf_bfr
# fit_release_bfr = MASS::fitdistr(release_bfr[which(release_bfr > 0)], "exponential")

#Negative side of carbon release distribution after project
#(represents carbon drawdown; conservative estimating that release after project happens at the same rate as drawdown during project)
release_aftr = rexp_p_aftr - rexp_cf_aftr
fit_release_aftr = MASS::fitdistr(-release_aftr[which(release_aftr < 0)], "exponential")

# png(paste0(getwd(), "/", site_label, "_samplediv_fig4a.png"), width = 600, height = 600)
# color_list = rep("grey", length(break_seq))
# color_list[break_seq >= 0] = "black"
# hist(release_bfr, freq = F, breaks = break_seq, col = color_list, xlab = "Release", ylab = "Prob. density", main = "a. Carbon release (before project start)")
# dev.off()

break_seq = seq(-max(abs(range(release_aftr))) - 1000, max(abs(range(release_aftr))) + 1000, len = 15)
png(paste0(file_path, site_label, "_fig5_distr_release_aftr.png"), width = 600, height = 600)
color_list = rep("grey", length(break_seq))
color_list[break_seq < 0] = "black"
hist(release_aftr, freq = F, breaks = break_seq, col = color_list, xlab = "Release", ylab = "Prob. density", main = "Carbon release after project start")
dev.off()


# Estimate carbon drawdown schedule ----
eval_end_additionality = sum(flux_series$additionality)


# Make functions to evaluate permanence on an annual basis ----
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


# Evaluate permanence (at t_eval) under simulated release schedule ----

# #Before-project distribution
# release_schedule_ann_bfr = vector("list", 100)
# len_release_bfr = rep(NA, 100)
# ep_vec_bfr = rep(NA, 100)
# for(i in 1:100){
#   release_amount = rexp(100, fit_release_bfr$estimate)
#   release_cumul = cumsum(release_amount)
#   if(release_cumul[1] >= eval_end_additionality){
#     release_simul = c(eval_end_additionality, 0)
#   } else {
#     release_simul = zero_clear(c(eval_end_additionality, eval_end_additionality - release_cumul))
#   }
#   
#   len_release_bfr[i] = length(release_simul) - 1
#   release_schedule_ann_bfr[[i]] = release_simul
#   release_extended = c(rep(eval_end_additionality, project_end - eval_end), release_simul)
#   names(release_extended) = eval_end:(eval_end + length(release_extended) - 1)
#   ep_vec_bfr[i] = ep_calc_ann(drawdown = eval_end_additionality, release = -diff(release_extended))$ep
# }

#Reversal of after-project distribution
release_schedule_ann_aftr = vector("list", 100)
len_release_aftr = rep(NA, 100)
ep_vec_aftr = rep(NA, 100)
for(i in 1:100){
  release_amount = rexp(100, fit_release_aftr$estimate)
  release_cumul = cumsum(release_amount)
  if(release_cumul[1] >= eval_end_additionality){
    release_simul = c(eval_end_additionality, 0)
  } else {
    release_simul = zero_clear(c(eval_end_additionality, eval_end_additionality - release_cumul))
  }
  
  len_release_aftr[i] = length(release_simul) - 1
  release_schedule_ann_aftr[[i]] = release_simul
  release_extended = c(rep(eval_end_additionality, project_end - eval_end), release_simul)
  names(release_extended) = eval_end:(eval_end + length(release_extended) - 1)
  ep_vec_aftr[i] = ep_calc_ann(drawdown = eval_end_additionality, release = -diff(release_extended))$ep
}


png(paste0(file_path, site_label, "_fig6a_release_time.png"), width = 600, height = 600)
hist(len_release_aftr, breaks = 1:20, cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
     xlab = "Number of years", main = "")
dev.off()

png(paste0(file_path, site_label, "_fig6b_ep.png"), width = 600, height = 600)
hist(ep_vec_aftr, breaks = seq(0, 0.5, by = 0.025), cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
     xlab = "eP", main = "")
dev.off()

save(scc, t0, tmin5, project_end, eval_end, years, eval_classes, site, site_label,
     agb_series_project_sim, flux_series_sim, summ_data, flux_series,
     fit_p_bfr, fit_cf_bfr, fit_p_aftr, fit_cf_aftr, JS_res, rexplist_p, rexplist_cf, lm_p_param, lm_cf_param, lm_js_div,
     fit_release_aftr, len_release_aftr, release_schedule_ann_aftr, ep_vec_aftr, file = paste0(file_path, site_label, ".Rdata"))


# Evaluate permanence (at t_eval) under currently implemented release schedule ----
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

save(site, site_label, t0, eval_end, project_end, eval_classes, scc,
     agb_series_project_sim, agb_series_leakage_sim, summ_data, cf_fit_tmin5, eval_end_additionality,
     release_schedule_ann, len_release, ep_vec, ep_calc_ann, zero_clear, file = paste0(getwd(), "/release_", site_label, ".Rdata"))