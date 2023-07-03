library(ggplot2)
library(ggnewscale)
library(tidyverse)
library(magrittr)
library(MASS)
library(stat)
source("divergence.R")


# Load data ----
load(file = paste0(getwd(), "/release_", "Gola_country", ".Rdata"))
load(file = paste0(getwd(), "/release_", "WLT_VNCC_KNT", ".Rdata"))
load(file = paste0(getwd(), "/release_", "CIF_Alto_Mayo", ".Rdata"))
site_label = site

flux_series = pivot_wider(subset(summ_data, series == "mean" & years >= t0), names_from = "var", values_from = "val") %>%
  dplyr::select(-one_of(c("series", "treatment_leak", "control_leak", "leakage", "net_additionality")))

write.table(flux_series, file = paste0(getwd(), "/", site_label, "_samplediv_flux_series.csv"), sep = ",", row.names = F)


# Graph carbon flux over time ----
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
ggsave(paste0(getwd(), "/", site_label, "_samplediv_fig1.png"), width = 15, height = 15, unit = "cm")


# Fit C flux to exponential distributions ----
flux_p_bfr = subset(summ_data, var == "treatment_proj" & years < t0 & val < 0)$val #325
fit_p_bfr = MASS::fitdistr(-flux_p_bfr, "exponential")

flux_cf_bfr = subset(summ_data, var == "control_proj" & years < t0 & val < 0)$val #353
fit_cf_bfr = MASS::fitdistr(-flux_cf_bfr, "exponential")

flux_p_aftr = subset(summ_data, var == "treatment_proj" & years >= t0 & val < 0)$val #209
fit_p_aftr = MASS::fitdistr(-flux_p_aftr, "exponential")

flux_cf_aftr = subset(summ_data, var == "control_proj" & years >= t0 & val < 0)$val #210
fit_cf_aftr = MASS::fitdistr(-flux_cf_aftr, "exponential")


# Calculate JS divergence and do permutation test ----
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
JS_bfr = JSTest(rexp_p_bfr, rexp_cf_bfr) #N.S.
a1 = Sys.time()
a1 - a0

a0 = Sys.time()
JS_aftr = JSTest(rexp_p_aftr, rexp_cf_aftr) #signif. diff.
a1 = Sys.time()
a1 - a0

a0 = Sys.time()
JS_cf = JSTest(rexp_cf_bfr, rexp_cf_aftr) #signif. diff.
a1 = Sys.time()
a1 - a0

a0 = Sys.time()
JS_p = JSTest(rexp_p_bfr, rexp_p_aftr) #signif. diff.
a1 = Sys.time()
a1 - a0


# Visualise distributions ----

png(paste0(getwd(), "/", site_label, "_samplediv_fig2_new.png"), width = 600, height = 600)
plot(0:1e-05, 0:1e-05, xlim = c(0, 200000), ylim = c(0, 4e-05), type = "n",
     xlab = "C flux", ylab = "Density")
lines(density(rexp_p_bfr), col = "lightblue", lwd = 3)
lines(density(rexp_cf_bfr), col = "pink", lwd = 3)
lines(density(rexp_p_aftr), col = "blue", lwd = 3)
lines(density(rexp_cf_aftr), col = "red", lwd = 3)
legend(x = "topright", legend = c("Project (before t0)", "Counterfactual (before t0)", "Project (after t0)", "Counterfactual (after t0)"),
       fill = c("lightblue", "pink", "blue", "red"))
dev.off()


# Test time-invariance of carbon flux distributions using 10-year intervals ----
years_interval = unique(subset(summ_data, years < t0)$years)
js_div = rep(NA, (length(years_interval) - 10 + 1))
sat_val = seq(0.3, 1, len = length(js_div))

png(paste0(getwd(), "/", site_label, "_samplediv_fig3.png"), width = 600, height = 600)
plot(0:1e-05, 0:1e-05, xlim = c(0, 200000), ylim = c(0, 4e-05), type = "n",
     xlab = "C flux", ylab = "Density")
for(i in 1:(length(years_interval) - 10 + 1)){
  p_dat = subset(summ_data, var == "treatment_proj" & years >= years_interval[i] & years <= years_interval[i + 9] & val < 0)$val
  fit_p = MASS::fitdistr(-p_dat, "exponential")
  rexp_p = rexp(1000, fit_p$estimate)
  
  cf_dat = subset(summ_data, var == "control_proj" & years >= years_interval[i] & years <= years_interval[i + 9] & val < 0)$val
  fit_cf = MASS::fitdistr(-cf_dat, "exponential")
  rexp_cf = rexp(1000, fit_cf$estimate)
  
  js_div[i] = compute_js_divergence(rexp_p, rexp_cf, bins = seq(0, ceiling(max(rexp_p, rexp_cf) / 10 ^ 4) * 10 ^ 4, by = 10000))
  lines(density(rexp_p), col = hsv(2/3, sat_val[i], 1), lwd = 3)
  lines(density(rexp_cf), col = hsv(0, sat_val[i], 1), lwd = 3)
}
legend(x = "topright", legend = c("Project", "Counterfactual"), fill = c("blue", "red"))
dev.off()


# Construct release distributions ----
#Approach 1: the right hand side of distribution before project (additionality release due to random drift)
release_bfr = rexp_p_bfr - rexp_cf_bfr
fit_release_bfr = MASS::fitdistr(release_bfr[which(release_bfr > 0)], "exponential")

#Approach 2: the left hand side of distribution after project (conservative estimate of release at the same rate as accumulation)
release_aftr = rexp_p_aftr - rexp_cf_aftr
fit_release_aftr = MASS::fitdistr(-release_aftr[which(release_aftr < 0)], "exponential")


break_seq = seq(-max(abs(range(release_bfr, release_aftr))) - 1000, max(abs(range(release_bfr, release_aftr))) + 1000, len = 15)

png(paste0(getwd(), "/", site_label, "_samplediv_fig4a.png"), width = 600, height = 600)
color_list = rep("grey", length(break_seq))
color_list[break_seq >= 0] = "black"
hist(release_bfr, freq = F, breaks = break_seq, col = color_list, xlab = "Release", ylab = "Prob. density", main = "a. Carbon release (before project start)")
dev.off()

png(paste0(getwd(), "/", site_label, "_samplediv_fig4b.png"), width = 600, height = 600)
color_list = rep("grey", length(break_seq))
color_list[break_seq < 0] = "black"
hist(release_aftr, freq = F, breaks = break_seq, col = color_list, xlab = "Release", ylab = "Prob. density", main = "b. Carbon release (after project start)")
dev.off()


png(paste0(getwd(), "/", site_label, "_samplediv_fig4c.png"), width = 600, height = 600)
plot(0:1e-05, 0:1e-05, xlim = c(0, 200000), ylim = c(0, 4e-05), type = "n",
     xlab = "C flux", ylab = "Density")
lines(density(rexp(1000, fit_release_bfr$estimate)), col = "blue")
lines(density(rexp(1000, fit_release_aftr$estimate)), col = "red")
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

#Before-project distribution
release_schedule_ann_bfr = vector("list", 100)
len_release_bfr = rep(NA, 100)
ep_vec_bfr = rep(NA, 100)
for(i in 1:100){
  release_amount = rexp(100, fit_release_bfr$estimate)
  release_cumul = cumsum(release_amount)
  if(release_cumul[1] >= eval_end_additionality){
    release_simul = c(eval_end_additionality, 0)
  } else {
    release_simul = zero_clear(c(eval_end_additionality, eval_end_additionality - release_cumul))
  }
  
  len_release_bfr[i] = length(release_simul) - 1
  release_schedule_ann_bfr[[i]] = release_simul
  release_extended = c(rep(eval_end_additionality, project_end - eval_end), release_simul)
  names(release_extended) = eval_end:(eval_end + length(release_extended) - 1)
  ep_vec_bfr[i] = ep_calc_ann(drawdown = eval_end_additionality, release = -diff(release_extended))$ep
}

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


png(paste0(getwd(), "/", site_label, "_samplediv_fig5a.png"), width = 600, height = 600)
hist(len_release_aftr, breaks = 1:(max(len_release_bfr, len_release_aftr) + 1), cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
     xlab = "Number of years", main = "")
dev.off()

png(paste0(getwd(), "/", site_label, "_samplediv_fig5b.png"), width = 600, height = 600)
hist(len_release_bfr, breaks = 1:(max(len_release_bfr, len_release_aftr) + 1), cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
     xlab = "Number of years", main = "")
dev.off()

png(paste0(getwd(), "/", site_label, "_samplediv_fig6a.png"), width = 600, height = 600)
hist(ep_vec_aftr, breaks = seq(0, 1, by = 0.05), cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
     xlab = "eP", main = "")
dev.off()

png(paste0(getwd(), "/", site_label, "_samplediv_fig6b.png"), width = 600, height = 600)
hist(ep_vec_bfr, breaks = seq(0, 1, by = 0.05), cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
     xlab = "eP", main = "")
dev.off()


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
