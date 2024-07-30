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
    isExPost = (year_i <= t_max)

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
    if(is.null(aomega)) aomega = quantile(samp_additionality, omega)

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
subfolder = paste0(type, "/")
summ_common = list(t0 = t0,
                   t_max = t_max,
                   additionality = summ_additionality,
                   credit = summ_credit,
                   release = summ_release,
                   aomega = summ_aomega,
                   ep = summ_ep,
                   risk = summ_risk)

if(type == "hypo") {
  dd_rate_text = gsub("\\.", "_", as.character(dd_rate))
  ppr_text = gsub("\\.", "_", as.character(postproject_ratio))
  use_theo_text = ifelse(use_theo, "", "_sampled_aomega")

  if(hypo_sensit != "none") {
    subfolder = paste0("sensitivity_", hypo_sensit, "/")
  } else {
    subfolder = paste0("hypo", use_theo_text, "/")
  }

  file_pref = paste0("dd_rate_", dd_rate_text, switch(hypo_sensit,
                                                      "dd_rate" = "",
                                                      "warmup" = paste0("_warmup_", warmup),
                                                      "ppr" = paste0("_ppr_", ppr_text),
                                                      "H" = paste0("_H_", H),
                                                      "none" = use_theo_text))
  
  summ_hypo = list(type = type,
                   sensitivity = hypo_sensit,
                   dd_rate = dd_rate)
  summ = c(summ_hypo, switch(hypo_sensit,
                        "dd_rate" = NULL,
                        "warmup" = list(warmup = warmup),
                        "ppr" = list(ppr = postproject_ratio),
                        "H" = list(H = H),
                        "none" = NULL))
} else if(type == "real"){
  file_pref = switch(project_site,
                     "Gola_country" = "Gola", #1201
                     "WLT_VNCC_KNT" = "KNT",
                     "CIF_Alto_Mayo" = "Alto_Mayo", #944
                     "VCS_1396" = "RPA",
                     "VCS_934" = "Mai_Ndombe")
  
  summ = list(type = type,
              sensitivity = hypo_sensit,
              project = project_site,
              flux = flux_series)
} else if(type == "hypo_aggr") {
  file_pref = paste0(type, "_", hypo_aggr_type)
  summ = list(type = type,
              sensitivity = hypo_sensit,
              dd_rate = dd_rate_vec)
} else if(type == "real_aggr") {
  file_pref = paste0(type, "_", real_aggr_type)
  summ = list(type = type,
              sensitivity = hypo_sensit,
              project = sites,
              flux = flux_series_list)
}
summ = c(summ, summ_common)

saveRDS(summ, file = paste0(file_path, subfolder, file_pref, "_output.RDS"))