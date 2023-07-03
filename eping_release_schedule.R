library(MASS)
library(rmarkdown)
library(ggplot2)
library(ggnewscale)
library(EnvStats)
library(extRemes)
library(tidyverse)
library(magrittr)

rmarkdown::render("/Users/E-Ping Rau/OneDrive - University of Cambridge/4C_evaluations/R/Reports/eval_template/evaluation_epingrau.Rmd", clean = FALSE)

#t0, eval_end, project_end
tmin5 = t0 - 5
eval_end = 2021
years = 1990:eval_end
eval_classes = c(1, 2, 3, 4) #only evaluate classes 1-4
if(country_match & site == "Gola"){
  site_label = paste0(site, "_country")
} else {
  site_label = site
}

# Calculate net additionality ----
makeFlux = function(project_series, leakage_series){
  project_aggr = aggregate(class_co2e ~ treatment + year, subset(project_series, class %in% eval_classes), FUN = sum) 
  leakage_aggr = aggregate(class_co2e ~ treatment + year, subset(leakage_series, class %in% eval_classes), FUN = sum)
  treatment_proj = subset(project_aggr, treatment == "treatment")$class_co2e
  control_proj = subset(project_aggr, treatment == "control")$class_co2e
  treatment_leak = subset(leakage_aggr, treatment == "treatment")$class_co2e
  control_leak = subset(leakage_aggr, treatment == "control")$class_co2e
  flux_series = data.frame(year = years[-1],
                           treatment_proj = diff(treatment_proj),
                           control_proj = diff(control_proj),
                           treatment_leak = diff(treatment_leak),
                           control_leak = diff(control_leak))
  flux_series %<>%
    mutate(additionality = treatment_proj - control_proj,
           leakage = treatment_leak - control_leak)
  return(flux_series)
}

flux_series_sim = mapply(makeFlux, agb_series_project_sim, agb_series_leakage_sim,
                         SIMPLIFY = F)

summariseSeries = function(in_list, sel_col){
  df = as.data.frame(
    sapply(in_list, function(x) x[, sel_col]),
    colnames = 1:ncol(df))
  df_mean = data.frame(series = "mean",
                       var = sel_col,
                       years = years[-1],
                       val = apply(df, 1, mean))
  df %<>%
    mutate(var = sel_col, years = years[-1]) %>%
    pivot_longer(V1:V20, names_to = "series", values_to = "val")
  df = rbind(df, df_mean)
  return(df)
}

summ_data = rbind(summariseSeries(flux_series_sim, "treatment_proj"),
                  summariseSeries(flux_series_sim, "control_proj"),
                  summariseSeries(flux_series_sim, "treatment_leak"),
                  summariseSeries(flux_series_sim, "control_leak"),
                  summariseSeries(flux_series_sim, "additionality"),
                  summariseSeries(flux_series_sim, "leakage"))

calcNetAdd = function(add, leak) pmax(ifelse(add > 0, 0, add), add + pmin(0, pmin(0, leak) * add / abs(add)))

netadd_df = data.frame(var = "net_additionality",
                       years = years[-1],
                       series = "mean",
                       val = calcNetAdd(subset(summ_data, var == "additionality" & series == "mean")$val,
                                        subset(summ_data, var == "leakage" & series == "mean")$val))
summ_data = rbind(summ_data, netadd_df)
proj_dat = subset(summ_data, years > t0 & var %in% c("treatment_proj", "control_proj"))
leak_dat = subset(summ_data, years > t0 & var %in% c("treatment_leak", "control_leak"))
add_dat = subset(summ_data, years > t0 & var %in% c("additionality", "leakage", "net_additionality"))

lim0 = sapply(list(proj_dat, leak_dat, add_dat), function(x) floor(range(x$val)[1] / (10^5)) * (10^5))
lim1 = sapply(list(proj_dat, leak_dat, add_dat), function(x) ceiling(range(add_dat$val)[2] / (10^5)) * (10^5))
lim_to_use = c(min(lim0), max(lim1))

lim_to_use = c(-400000, 0)

#fig 1a
ggplot(data = proj_dat, aes(col = var)) +
  geom_line(aes(x = years, y = val, lwd = series, alpha = series), show.legend = F) +
  scale_linewidth_manual("", values = c(1, rep(0.5, 20)), labels = c("Mean", "Simulated", rep(NULL, 19))) +
  scale_alpha_manual("", values = c(1, rep(0.1, 20)), labels = c("Mean", "Simulated", rep(NULL, 19))) +
  ggnewscale::new_scale("lwd") +
  ggnewscale::new_scale("alpha") +
  geom_line(data = subset(proj_dat, series %in% c("mean", "V1")), aes(x = years, y = val, lwd = series, alpha = series)) +
  scale_color_manual("", values = c("red", "black"), labels = c("Project counterfactual", "Project")) +
  scale_linewidth_manual("", values = c(1, 0.5), labels = c("Mean", "Simulated")) +
  scale_alpha_manual("", values = c(1, 0.1), labels = c("Mean", "Simulated")) +
  scale_x_continuous(name = "Year", breaks = (t0 + 1):eval_end) + 
  scale_y_continuous(name = "Flux (Mg CO2e)", limits = lim_to_use) + 
  ggtitle("a. Carbon fluxes (project)") +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        legend.title = element_text())
ggsave(paste0(getwd(), "/", site_label, "_fig1a.pdf"), width = 15, height = 15, unit = "cm")

#fig 1b
ggplot(data = leak_dat, aes(col = var)) +
  geom_line(aes(x = years, y = val, lwd = series, alpha = series), show.legend = F) +
  scale_linewidth_manual("", values = c(1, rep(0.5, 20)), labels = c("Mean", "Simulated", rep(NULL, 19))) +
  scale_alpha_manual("", values = c(1, rep(0.1, 20)), labels = c("Mean", "Simulated", rep(NULL, 19))) +
  ggnewscale::new_scale("lwd") +
  ggnewscale::new_scale("alpha") +
  geom_line(data = subset(leak_dat, series %in% c("mean", "V1")), aes(x = years, y = val, lwd = series, alpha = series)) +
  scale_color_manual("", values = c("red", "black"), labels = c("Leakage counterfactual", "Leakage")) +
  scale_linewidth_manual("", values = c(1, 0.5), labels = c("Mean", "Simulated")) +
  scale_alpha_manual("", values = c(1, 0.1), labels = c("Mean", "Simulated")) +
  scale_x_continuous(name = "Year", breaks = (t0 + 1):eval_end) + 
  scale_y_continuous(name = "Flux (Mg CO2e)", limits = lim_to_use) + 
  ggtitle("b. Carbon fluxes (leakage)") +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        legend.title = element_text())
ggsave(paste0(getwd(), "/", site_label, "_fig1b.pdf"), width = 15, height = 15, unit = "cm")

#fig 2
ggplot(data = add_dat, aes(col = var)) +
  geom_line(aes(x = years, y = val, lwd = series, alpha = series), show.legend = F) +
  scale_linewidth_manual("", values = c(1, rep(0.5, 20)), labels = c("Mean", "Simulated", rep(NULL, 19))) +
  scale_alpha_manual("", values = c(1, rep(0.1, 20)), labels = c("Mean", "Simulated", rep(NULL, 19))) +
  ggnewscale::new_scale("lwd") +
  ggnewscale::new_scale("alpha") +
  geom_line(data = subset(add_dat, series %in% c("mean", "V1")), aes(x = years, y = val, lwd = series, alpha = series)) +
  scale_color_manual("", values = c("blue", "red", "black"), labels = c("Additionality", "Leakage", "Net additionality")) +
  scale_linewidth_manual("", values = c(1, 0.5), labels = c("Mean", "Simulated")) +
  scale_alpha_manual("", values = c(1, 0.1), labels = c("Mean", "Simulated")) +
  scale_x_continuous(name = "Year", breaks = (t0 + 1):eval_end) + 
  scale_y_continuous(name = "Additionality (Mg CO2e)", limits = lim_to_use) + 
  ggtitle("Additionality and leakage") +
  theme(axis.line = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(),
        legend.title = element_text())
ggsave(paste0(getwd(), "/", site_label, "_fig2.pdf"), width = 15, height = 15, unit = "cm")

print(paste("Total net additionality:", round(sum(subset(add_dat, var == "net_additionality")$val), 4), "Mg CO2 e"))


# Fit counterfactual release distribution ----
fitCfRelease = function(from_time, plot_filename, plot_title){
  cf_release = subset(summ_data, years >= from_time & var == "control_proj" & val < 0)$val * -1
  #in years where flux is positive, consider release to be 0: distribution bound on the left (0)
  
  gof_lnorm = EnvStats::gofTest(y = cf_release, distribution = "lnorm")
  gof_gamma = gofgamma::test.AD(cf_release)
  gof_weibull = EWGoF::WLK.test(cf_release)
  
  cf_release_weibull = MASS::fitdistr(cf_release, "weibull") #fit with weibull
  distr_dens_weibull = rweibull(1000, cf_release_weibull$estimate[1], cf_release_weibull$estimate[2])
  
  fit = glm(formula = cf_release ~ 1, family = Gamma)
  cf_release_gamma = MASS::fitdistr(cf_release * coef(fit), "gamma")
  cf_release_gamma$scaling_multiplier = unname(coef(fit))
  distr_dens_gamma = rgamma(1000, cf_release_gamma$estimate[1], cf_release_gamma$estimate[2]) / unname(coef(fit))
  
  hist_val = hist(cf_release, plot = F)
  ylim1 = max(density(distr_dens_weibull)$y, density(distr_dens_gamma)$y, hist_val$density)
  
  png(plot_filename, width = 600, height = 600)
  hist(cf_release, freq = F, ylim = c(0, ylim1), cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
       xlab = "Carbon release (Mg CO2 e)", xaxt = "n",
       main = plot_title)
  axis(1, at = hist_val$breaks, labels = hist_val$breaks)
  lines(density(distr_dens_weibull), col = "blue")
  #lines(density(distr_dens_gamma), col = "red")
  dev.off()
  
  return(list(cf_release = cf_release, pval_lnorm = gof_lnorm$p.value, reject_gamma = gof_gamma$Decision, pval_weibull = gof_weibull$p.value,
              param_weibull = cf_release_weibull$estimate[1:2]))
}

#from t-5
cf_fit_tmin5 = fitCfRelease(from_time = tmin5,
                            plot_filename = paste0(getwd(), "/", site_label, "_fig3.png"),
                            plot_title = "")


# Estimate carbon drawdown schedule ----
eval_end_additionality = sum(subset(summ_data, years > t0 & series == "mean" & var == "net_additionality")$val)


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
release_schedule_ann = vector("list", 100)
len_release = rep(NA, 100)
ep_vec = rep(NA, 100)
for(i in 1:100){
  release_amount = rweibull(100, cf_fit_tmin5$param_weibull[1], cf_fit_tmin5$param_weibull[2])
  release_cumul = cumsum(release_amount)
  if(release_cumul[1] >= eval_end_additionality){
    release_simul = c(eval_end_additionality, 0)
  } else {
    release_simul = zero_clear(c(eval_end_additionality, eval_end_additionality - release_cumul))
  }
  
  len_release[i] = length(release_simul) - 1
  release_schedule_ann[[i]] = release_simul
  release_extended = c(rep(eval_end_additionality, project_end - eval_end), release_simul)
  names(release_extended) = eval_end:(eval_end + length(release_extended) - 1)
  ep_vec[i] = ep_calc_ann(drawdown = eval_end_additionality, release = -diff(release_extended))$ep
}

png(paste0(getwd(), "/", site_label, "_fig4.png"), width = 600, height = 600)
hist(len_release, breaks = 1:15, cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
     xlab = "Number of years", main = "")
dev.off()

png(paste0(getwd(), "/", site_label, "_fig5.png"), width = 600, height = 600)
hist(ep_vec, cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
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

load(file = paste0(getwd(), "/release_", "Gola", ".Rdata"))
