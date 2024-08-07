# 1. Unpack input variables ----
type = inpar$type
type_label = inpar$type_label
t0 = inpar$t0
t_max = inpar$t_max
obs_p_loss = inpar$obs_p_loss
obs_c_loss = inpar$obs_c_loss
postproject_release = inpar$postproject_release
aomega = inpar$aomega

if(type == "theo") {
  mean_drawdown = inpar$mean_drawdown
  lambdaP = inpar$lambdaP
  lambdaC = inpar$lambdaC
} else {
  sites = inpar$sites
  c_loss_p_list = inpar$c_loss_p_list
  c_loss_c_list = inpar$c_loss_c_list
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
sim_reversal = matrix(F, H, n_rep)
sim_buffer = matrix(0, H, n_rep)
sim_schedule = vector("list", n_rep)


verbose = T
a = Sys.time()
for(j in 1:n_rep){
  a1 = Sys.time()
  cat("Repetition:", j, "\n")
  schedule = matrix(0, H, H_max_scc) #release schedule
  buffer_pool = 0
  #cat("Buffer at start: ", buffer_pool, "\n")
  for(i in 1:H){
    year_i = t0 + i - 1
    isExPost = (year_i <= t_max)

    #get carbon loss values
    #sample from fitted distributions for hypothetical projects or in years where ex post values are not available
    #use ex post values for real-life projects in years where they are available
    if(type == "theo") {
      sim_p_loss[i, j] = sum(sapply(lambdaP, function(x) rexp(1, x)))
      sim_c_loss[i, j] = sum(sapply(lambdaC, function(x) rexp(1, x)))
    } else if(type == "real") {
      sim_p_loss[i, j] = ifelse(isExPost, obs_p_loss[i], sum(sapply(c_loss_p_fit_list, function(x) SampGMM(x, n = 1))))
      sim_c_loss[i, j] = ifelse(isExPost, obs_c_loss[i], sum(sapply(c_loss_c_fit_list, function(x) SampGMM(x, n = 1))))
    }

    #calculate a-omega: sample-based unless in single hypothetical project
    if(type == "theo") {
      samp_additionality = mapply(function(x, y) rexp(1000, x) - rexp(1000, y),
                                  x = lambdaC, y = lambdaP) %>%
        apply(1, sum)
    } else if(type == "real") {
      c_loss_p_fit_list = lapply(c_loss_p_list, function(x) {
        FitGMM(subset(x, year <= year_i & year >= t0)$val)})
      c_loss_p_samp_list = lapply(c_loss_p_fit_list, function(x) SampGMM(x, n = 1000))
      
      c_loss_c_fit_list = lapply(c_loss_c_list, function(x) {
        FitGMM(subset(x, year <= year_i & year >= t0)$val)})
      c_loss_c_samp_list = lapply(c_loss_c_fit_list, function(x) SampGMM(x, n = 1000))
      
      samp_additionality = mapply(function(x, y) x - y,
                                  x = c_loss_c_samp_list,
                                  y = c_loss_p_samp_list) %>%
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
        sim_reversal[i, j] = T
      }
      sim_pact[i, j]  = sim_credit[i, j] * sim_ep[i, j]
    }
    sim_buffer[i, j] = buffer_pool
    if(verbose) cat("Step", i, ":", "additionality", sim_additionality[i, j], "| anticipated release", sim_release[i, j], "| credit", sim_credit[i, j], "| EP", sim_ep[i, j], "\n")
  }
  sim_schedule[[j]] = schedule
  b1 = Sys.time()
  if(verbose) cat("Repetition", j, "runtime:", b1 - a1, "\n")
}
b = Sys.time()
if(verbose) cat("Total runtime:", b - a, "\n")


# 3. Summarise results (time series) ----
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

#per-year reversal risk (proportion of repetitions without positive credits)
summ_risk = data.frame(year = 1:H, risk = apply(sim_reversal, 1, sum) / ncol(sim_reversal))
summ_risk$risk[1:warmup] = NA

#4. Set file prefixes and save results ----
summ_simulation = list(additionality = summ_additionality,
                       credit = summ_credit,
                       release = summ_release,
                       aomega = summ_aomega,
                       ep = summ_ep,
                       risk = summ_risk)

summ_complete = c(common_var, summ_simulation)
dir.create(paste0(out_path, type_label, "/"))
saveRDS(summ_complete, file = paste0(out_path, type_label, "/", type_label, "_output.rds"))


if(hypo_sensit != "none") {
    ppr_text = gsub("\\.", "_", as.character(postproject_ratio))
    drawdown_text = gsub("\\.", "_", as.character(mean_drawdown))
    subfolder = paste0("sensitivity_", hypo_sensit, "/")
    file_pref = paste0("drawdown_", drawdown_text, switch(hypo_sensit,
                                                          "dd_rate" = "",
                                                          "warmup" = paste0("_warmup_", warmup),
                                                          "ppr" = paste0("_ppr_", ppr_text),
                                                          "H" = paste0("_H_", H),
                                                          "none" = ""))
    summ_hypo = list(type = type,
                     sensitivity = hypo_sensit,
                     mean_drawdown = mean_drawdown,
                     switch(hypo_sensit,
                            "dd_rate" = NULL,
                            "warmup" = list(warmup = warmup),
                            "ppr" = list(ppr = postproject_ratio),
                            "H" = list(H = H),
                            "none" = NULL))
}