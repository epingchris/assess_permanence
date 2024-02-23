# Set input parameters ----

## A. Hypothetical projects ----
if(type == "hypo") {
  t0 = 2021
  
  year_expost = t0 - 1 #actually no ex post data will be used for hypothetical projects
  
  lambdaP = 1
  lambdaC = 1 / dd_rate
  
  expost_p_loss = rexp(1000, lambdaP)
  expost_c_loss = rexp(1000, lambdaC)
  
  #post-project release rate: a ratio of the counterfactual release rate during project (default to double)
  postproject_release = postproject_ratio / lambdaC
  
  samp_additionality = rexp(1000, lambdaC) - rexp(1000, lambdaP)
  
  aomega = 1 / lambdaP * log(omega * (lambdaP + lambdaC) / lambdaC) #analytical solution

  #input: t0, year_expost, lambdaP/C, expost_lossP/C, postproject_release, aomega
## B. Real-life projects ----
} else if(type == "real") {
  use_theo = F
  load(file = paste0(file_path, "project_input_data/", project_site, ".Rdata")) #t0 included
  
  flux_series_sim = mapply(function(x, y) makeFlux(project_series = x, leakage_series = y)$flux,
                           x = agb_series_project_sim,
                           y = vector("list", length = length(agb_series_project_sim)),
                           SIMPLIFY = F)
  summ_flux = rbind(summariseSeries(flux_series_sim, "treatment_proj"),
                    summariseSeries(flux_series_sim, "control_proj"))
  year_expost = max(summ_flux$year)
  
  absloss_p_init = summ_flux %>%
    subset(var == "treatment_proj" & year >= t0 & series != "mean") %>%
    mutate(val = val * (-1), var = NULL, series = NULL)
  absloss_c_init = summ_flux %>%
    subset(var == "control_proj" & year >= t0 & series != "mean") %>%
    mutate(val = val * (-1), var = NULL, series = NULL)
  
  expost_p_loss = absloss_p_init %>%
    group_by(year) %>%
    summarise(val = mean(val), .groups = "drop") %>%
    pull(val)
  expost_c_loss = absloss_c_init %>%
    group_by(year) %>%
    summarise(val = mean(val), .groups = "drop") %>%
    pull(val)
  
  #post-project release rate: a ratio of the counterfactual release rate during project (default to double)
  absloss_p_fit = FitGMM(absloss_p_init$val)
  absloss_c_fit = FitGMM(absloss_c_init$val)
  postproject_release = mean(SampGMM(absloss_c_fit, n = 1000)) * postproject_ratio
  # 
  # add_samp = SampGMM(absloss_c_fit, n = 1000) - SampGMM(absloss_p_fit, n = 1000)
  # aomega = quantile(add_samp, omega)

  #input: t0, year_expost, absloss_p/c_init, expost_lossP/C, postproject_release, aomega

## C. Aggregated hypothetical projects ----
} else if(type == "hypo_aggr") {
  use_theo = F
  
  t0 = 2021
  
  year_expost = t0 - 1 #actually no ex post data will be used for hypothetical projects
  
  dd_rate_vec = switch(hypo_aggr_type,
                 "A" = c(1.5, 1.5, 1.5, 10),
                 "B" = c(1.5, 1.5, 10, 10),
                 "C" = c(1.5, 10, 10, 10))
  lambdaP_vec = rep(1, length(dd_rate_vec))
  lambdaC_vec = 1 / dd_rate_vec
  
  absloss_p_samp_list = lapply(lambdaP_vec, function(x) rexp(1000, x))
  absloss_c_samp_list = lapply(lambdaC_vec, function(x) rexp(1000, x))
  expost_p_loss = apply(as.data.frame(absloss_p_samp_list), 1, sum)
  expost_c_loss = apply(as.data.frame(absloss_c_samp_list), 1, sum)
  
  #post-project release rate: a ratio of the counterfactual release rate during project (default to double)
  postproject_release = sum(postproject_ratio / lambdaC_vec)
  # 
  # #a-omega: use sampling approach because no analytical a-omega exists yet for portfolio
  # add_samp = mapply(function(x, y) x - y,
  #                   x = absloss_c_samp_list, y = absloss_p_samp_list) %>%
  #   apply(1, sum)
  # aomega = quantile(add_samp, omega)
  # 
  #input: t0, year_expost, lambdaP/C_vec, expost_lossP/C, postproject_release, aomega
  
## D. Aggregated real-life projects ----
} else if(type == "real_aggr") {
  use_theo = F
  
  sites = switch(real_aggr_type,
                 "five" = c("Gola_country", "WLT_VNCC_KNT", "CIF_Alto_Mayo", "VCS_1396", "VCS_934"),
                 "four" = c("Gola_country", "CIF_Alto_Mayo", "VCS_1396", "VCS_934"),
                 "three" = c("Gola_country", "CIF_Alto_Mayo", "VCS_1396"))
  
  summ_flux = vector("list", length(sites))
  absloss_p_init_list = vector("list", length(sites))
  absloss_c_init_list = vector("list", length(sites))
  t0_vec = rep(NA, length(sites))
  
  for(s in sites){
    i = which(sites %in% s)
    load(file = paste0(file_path, "project_input_data/", s, ".Rdata")) #load data
    t0_vec[i] = t0
    
    flux_series_sim = mapply(function(x, y) makeFlux(project_series = x, leakage_series = y)$flux,
                             x = agb_series_project_sim,
                             y = vector("list", length = length(agb_series_project_sim)),
                             SIMPLIFY = F)
    
    summ_flux[[i]] = rbind(summariseSeries(flux_series_sim, "treatment_proj"),
                           summariseSeries(flux_series_sim, "control_proj"))
    
    absloss_p_init_list[[i]] = summ_flux[[i]] %>%
      subset(var == "treatment_proj" & year >= t0 & series != "mean") %>%
      mutate(val = val * (-1), var = NULL, series = NULL)
    absloss_c_init_list[[i]] = summ_flux[[i]] %>%
      subset(var == "control_proj" & year >= t0 & series != "mean") %>%
      mutate(val = val * (-1), var = NULL, series = NULL)
  }
  
  site = "portfolio"
  t0 = min(t0_vec)
  
  year_expost = 2021

  expost_p_loss = mapply(function(x, y) x = x %>% mutate(site = y),
                               x = absloss_p_init_list,
                               y = sites, SIMPLIFY = F) %>%
    do.call(rbind, .) %>%
    group_by(year, site) %>%
    summarise(val = mean(val)) %>%
    ungroup(site) %>%
    summarise(val = sum(val), .groups = "drop") %>%
    pull(val)
  expost_c_loss = mapply(function(x, y) x = x %>% mutate(site = y),
                               x = absloss_c_init_list,
                               y = sites, SIMPLIFY = F) %>%
    do.call(rbind, .) %>%
    group_by(year, site) %>%
    summarise(val = mean(val)) %>%
    ungroup(site) %>%
    summarise(val = sum(val), .groups = "drop") %>%
    pull(val)

  #post-project release rate: a ratio of the counterfactual release rate during project (default to double)
  absloss_p_fit_list = lapply(absloss_p_init_list, function(x) FitGMM(x$val))
  absloss_c_fit_list = lapply(absloss_p_init_list, function(x) FitGMM(x$val))
  postproject_release = absloss_c_fit_list %>%
    sapply(function(x) SampGMM(x, n = 1000)) %>%
    apply(1, sum) %>%
    mean() * postproject_ratio
  
  # absloss_p_samp_list = lapply(absloss_p_fit_list, function(x) SampGMM(x, n = 1000))
  # absloss_c_samp_list = lapply(absloss_c_fit_list, function(x) SampGMM(x, n = 1000))
  # add_samp = mapply(function(x, y) x - y,
  #                   x = absloss_c_samp_list,
  #                   y = absloss_p_samp_list) %>%
  #   as.data.frame() %>%
  #   apply(1, function(x) sum(x, na.rm = T))
  # aomega = quantile(add_samp, omega)
  
  #input: t0, year_expost, absloss_p/c_init_list, expost_lossP/C, postproject_release, aomega
}


# 2. Perform simulations ----
H_max_scc = year_max_scc - t0 + 1

sim_p_loss = matrix(0, H, n_rep)
sim_c_loss = matrix(0, H, n_rep)
sim_additionality = matrix(0, H, n_rep)
sim_credit = matrix(0, H, n_rep)
sim_benefit = matrix(0, H, n_rep)
sim_aomega = matrix(0, H, n_rep)
sim_release = matrix(0, H_max_scc, n_rep)
sim_damage = matrix(0, H, n_rep)
sim_ep = matrix(0, H, n_rep)
sim_pact = matrix(0, H, n_rep)
sim_failure = matrix(F, H, n_rep)
sim_buffer = matrix(0, H, n_rep)
sim_schedule = vector("list", n_rep)

for(j in 1:n_rep){
  schedule = matrix(0, H, H_max_scc) #release schedule
  buffer_pool = 0
  #cat("Buffer at start: ", buffer_pool, "\n")
  for(i in 1:H){
    year_i = t0 + i - 1
    isExPost = year_i <= year_expost

    #get carbon loss values
    #sample from fitted distributions for hypothetical projects or in years where ex post values are not available
    #use ex post values for real-life projects in years where they are available
    if(type == "hypo") {
      sim_p_loss[i, j] = rexp(1, lambdaP)
      sim_c_loss[i, j] = rexp(1, lambdaC)
    } else if(type == "hypo_aggr") {
      sim_p_loss[i, j] = sum(sapply(lambdaP_vec, function(x) rexp(1, x)))
      sim_c_loss[i, j] = sum(sapply(lambdaC_vec, function(x) rexp(1, x)))
    } else if(type == "real") {
      sim_p_loss[i, j] = ifelse(isExPost, expost_p_loss[i], SampGMM(absloss_p_fit, n = 1))
      sim_c_loss[i, j] = ifelse(isExPost, expost_c_loss[i], SampGMM(absloss_c_fit, n = 1))
    } else if(type == "real_aggr") {
      sim_p_loss[i, j] = ifelse(isExPost, expost_p_loss[i], sum(sapply(absloss_p_fit_list, function(x) SampGMM(x, n = 1))))
      sim_c_loss[i, j] = ifelse(isExPost, expost_c_loss[i], sum(sapply(absloss_c_fit_list, function(x) SampGMM(x, n = 1))))
    }

    #calculate a-omega: sample-based unless in single hypothetical project
    if(type == "hypo") {
      samp_additionality = rexp(1000, lambdaC) - rexp(1000, lambdaP)
    } else if(type == "real") {
      absloss_p_fit = FitGMM(subset(absloss_p_init, year <= year_i)$val)
      absloss_c_fit = FitGMM(subset(absloss_c_init, year <= year_i)$val)
      samp_additionality = SampGMM(absloss_c_fit, n = 1000) - SampGMM(absloss_p_fit, n = 1000)
    } else if(type == "hypo_aggr") {
      samp_additionality = mapply(function(x, y) rexp(1000, x) - rexp(1000, y),
                                  x = lambdaC_vec, y = lambdaP_vec) %>%
        apply(1, sum)
    } else if(type == "real_aggr") {
      absloss_p_fit_list = lapply(absloss_p_init_list, function(x) {
        FitGMM(subset(x, year <= year_i & year >= t0)$val)})
      absloss_p_samp_list = lapply(absloss_p_fit_list, function(x) SampGMM(x, n = 1000))
      
      absloss_c_fit_list = lapply(absloss_c_init_list, function(x) {
        FitGMM(subset(x, year <= year_i & year >= t0)$val)})
      absloss_c_samp_list = lapply(absloss_c_fit_list, function(x) SampGMM(x, n = 1000))
      
      samp_additionality = mapply(function(x, y) x - y,
                                  x = absloss_c_samp_list,
                                  y = absloss_p_samp_list) %>%
        as.data.frame() %>%
        apply(1, function(x) sum(x, na.rm = T))
    }
    if(!use_theo) aomega = quantile(samp_additionality, omega)

    sim_additionality[i, j] = sim_c_loss[i, j] - sim_p_loss[i, j]
    sim_aomega[i, j] = aomega
    
    if(i <= warmup) {
      #first five years: releases ignored, drawdown added to project pool
      if(sim_additionality[i, j] > 0) buffer_pool = buffer_pool + sim_additionality[i, j]
      sim_credit[i, j] = 0
      #cat("Buffer at year", i, ": ", buffer_pool, "\n")
    } else {
      #from sixth year on: get credits and anticipated releases
      
      #use buffer pool to fill anticipated releases first
      #only deduct from buffer pool at each year if there is space left for that year
      if(buffer_pool > 0) {
        max_release = ifelse(aomega > 0, 0, -aomega) #if a-omega is positive, maximum release is zero
        can_be_released = min(max(0, max_release - sim_release[i, j]), buffer_pool)
        #if(j < 5) cat("at year", i, ", can be released from buffer =", max_release, "-", sim_release[i, j], "=", can_be_released, "\n")
        sim_release[i, j] = sim_release[i, j] + can_be_released
        buffer_pool = buffer_pool - can_be_released
        #if(j < 5) cat("total release now =", sim_release[i, j], ", left in buffer =", buffer_pool, "\n")
      }
      
      sim_credit[i, j] = sim_additionality[i, j] + sim_release[i, j]
      if(sim_credit[i, j] > 0){
        to_be_released = sim_credit[i, j]
        #cat("credits at year ", i, ": ", to_be_released, "\n")
        sim_benefit[i, j] = sim_credit[i, j] * filter(scc_extended, year == year_i)$central
        k = i #kth year(s), for which we estimate anticipated release
        while(to_be_released > 0 & k < H_max_scc){
          k = k + 1
          if(k > H) {
            max_release = postproject_release #post-project release rate
          } else {
            max_release = ifelse(aomega > 0, 0, -aomega) #if a-omega is positive, maximum release is zero
          }
          #cat("at year", k, ", can be released =", max_release, "-", sim_release[k, j], "=", can_be_released, "\n")
          
          can_be_released = max(0, max_release - sim_release[k, j])
          schedule[i, k] = min(to_be_released, can_be_released)
          to_be_released = to_be_released - schedule[i, k]
          sim_release[k, j] = sim_release[k, j] + schedule[i, k]
          #cat("actually released =", schedule[i, k], ", total release now =", sim_release[k, j], ", left to release =", to_be_released, "\n")
        }
        years_k = (t0 + i):year_max_scc
        sim_damage[i, j] = sum(schedule[i, seq_along(years_k) + i] * filter(scc_extended, year %in% years_k)$central / ((1 + D) ^ seq_along(years_k)))
        sim_ep[i, j] = (sim_benefit[i, j] - sim_damage[i, j]) / sim_benefit[i, j]
        #cat("Credit =", sim_credit[i, j], ", Benefit =", sim_benefit[i, j], ", Damage =", sim_damage[i, j], ", eP =", sim_ep[i, j], "\n")
      } else if(sim_credit[i, j] <= 0){
        sim_ep[i, j] = 0
        sim_failure[i, j] = T
      }
      sim_pact[i, j]  = sim_credit[i, j] * sim_ep[i, j]
    }
    sim_buffer[i, j] = buffer_pool
  }
  sim_schedule[[j]] = schedule
}

#view evolution of a particular iteration
if(view_snapshot) {
  j = 1
  snapshot = data.frame(additionality = sim_additionality[1:50, j],
                        aomega = sim_aomega[1:50, j],
                        release = sim_release[1:50, j],
                        credit = sim_credit[1:50, j],
                        PACT = sim_pact[1:50, j],
                        rsched = apply(sim_schedule[[j]], 1, sum),
                        buffer = sim_buffer[1:50, j])
  View(snapshot)
}

# 3. Summarise results ----
SummariseSim = function(mat){
  df = mat %>%
    as.data.frame() %>%
    reframe(
      year = row_number(),
      p05 = apply(., 1, function(x) quantile(x, 0.05, na.rm = T)),
      p25 = apply(., 1, function(x) quantile(x, 0.25, na.rm = T)),
      median = apply(., 1, median, na.rm = T),
      p75 = apply(., 1, function(x) quantile(x, 0.75, na.rm = T)),
      p95 = apply(., 1, function(x) quantile(x, 0.95, na.rm = T)),
      mean = apply(., 1, mean, na.rm = T),
      sd = apply(., 1, function(x) sd(x, na.rm = T)),
      ci_margin = qt(0.975, df = n_rep - 1) * sd / sqrt(n_rep),
      ci_low = mean - ci_margin,
      ci_high = mean + ci_margin
    )
  return(df)
}

summ_additionality = SummariseSim(sim_additionality)
summ_credit = SummariseSim(sim_credit)
#summ_pact = SummariseSim(sim_pact)
summ_release = SummariseSim(sim_release[1:H, ])
#summ_buffer = SummariseSim(sim_buffer)
summ_aomega = SummariseSim(sim_aomega)

summ_ep = sim_ep %>%
  replace(., . == 0, NA) %>%
  SummariseSim() %>%
  replace(., . == Inf| . == -Inf, NA)

#per-year failure risk (proportion of repetitions without positive credits)
summ_risk = data.frame(year = 1:H,
                       risk = apply(sim_failure, 1, sum) / ncol(sim_failure))
summ_risk$risk[1:warmup] = NA


#4. Set file prefixes and save results ----
if(hypo_sensit != "none") {
  subfolder = paste0("sensitivity_", hypo_sensit, "/")
} else {
  subfolder = paste0(type, "/")
}

if(type == "hypo") {
  dd_rate_text = gsub("\\.", "_", as.character(dd_rate))
  if(hypo_sensit == "dd_rate") {
    file_pref = paste0(hypo_sensit, "_", dd_rate_text)
    summ = list(type = type,
                sensitivity = hypo_sensit,
                dd_rate = dd_rate,
                additionality = summ_additionality,
                credit = summ_credit,
                release = summ_release,
                aomega = summ_aomega,
                ep = summ_ep,
                risk = summ_risk)
  } else if(hypo_sensit == "warmup") {
    file_pref = paste0("dd_rate_", dd_rate_text, "_", hypo_sensit, "_", warmup)
    summ = list(type = type,
                sensitivity = hypo_sensit,
                dd_rate = dd_rate,
                warmup = warmup,
                additionality = summ_additionality,
                credit = summ_credit,
                release = summ_release,
                aomega = summ_aomega,
                ep = summ_ep,
                risk = summ_risk)
  } else if(hypo_sensit == "ppr") {
    file_pref = paste0("dd_rate_", dd_rate_text, "_ppr_", gsub("\\.", "_", as.character(postproject_ratio)))
    summ = list(type = type,
                sensitivity = hypo_sensit,
                dd_rate = dd_rate,
                ppr = postproject_ratio,
                additionality = summ_additionality,
                credit = summ_credit,
                release = summ_release,
                aomega = summ_aomega,
                ep = summ_ep,
                risk = summ_risk)
  } else if(hypo_sensit == "H") {
    file_pref = paste0("dd_rate_", dd_rate_text, "_H_", H)
    summ = list(type = type,
                sensitivity = hypo_sensit,
                dd_rate = dd_rate,
                H = H,
                additionality = summ_additionality,
                credit = summ_credit,
                release = summ_release,
                aomega = summ_aomega,
                ep = summ_ep,
                risk = summ_risk)
  } else {
    subfolder = ifelse(use_theo,
                       "hypo/",
                       "hypo_sampled_aomega/")
    file_pref = paste0("dd_rate_", dd_rate_text, ifelse(use_theo, "", "_sampled_aomega"))
    summ = list(type = type,
                sensitivity = hypo_sensit,
                dd_rate = dd_rate,
                additionality = summ_additionality,
                credit = summ_credit,
                release = summ_release,
                aomega = summ_aomega,
                ep = summ_ep,
                risk = summ_risk)
  }
} else if(type == "real"){
  file_pref = paste0(switch(project_site,
                            "Gola_country" = "Gola",
                            "WLT_VNCC_KNT" = "KNT",
                            "CIF_Alto_Mayo" = "Alto_Mayo",
                            "VCS_1396" = "RPA",
                            "VCS_934" = "Mai_Ndombe"))
  
  summ = list(type = type,
              sensitivity = hypo_sensit,
              project = project_site,
              flux = summ_flux,
              t0 = t0,
              additionality = summ_additionality,
              credit = summ_credit,
              release = summ_release,
              aomega = summ_aomega,
              ep = summ_ep,
              risk = summ_risk)
} else if(type == "hypo_aggr") {
  file_pref = paste0(type, "_", hypo_aggr_type)
  summ = list(type = type,
              sensitivity = hypo_sensit,
              dd_rate = dd_rate_vec,
              additionality = summ_additionality,
              credit = summ_credit,
              release = summ_release,
              aomega = summ_aomega,
              ep = summ_ep,
              risk = summ_risk)
} else if(type == "real_aggr") {
  file_pref = paste0(type, "_", real_aggr_type)
  summ = list(type = type,
              sensitivity = hypo_sensit,
              project = sites,
              flux = summ_flux,
              t0 = t0,
              additionality = summ_additionality,
              credit = summ_credit,
              release = summ_release,
              aomega = summ_aomega,
              ep = summ_ep,
              risk = summ_risk)
}

file_pref = paste0(file_pref, "_test")

saveRDS(summ, file = paste0(file_path, subfolder, file_pref, "_output.RDS"))