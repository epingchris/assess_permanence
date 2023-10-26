absloss_p_init = subset(summ_flux, var == "treatment_proj" & year >= t0 & series != "mean") %>%
  mutate(val = val * (-1), var = NULL, series = NULL)
absloss_c_init = subset(summ_flux, var == "control_proj" & year >= t0 & series != "mean") %>%
  mutate(val = val * (-1), var = NULL, series = NULL)

additionality = absoloss_c_update$val - absoloss_p_update$val
samp_additionality = SampByPro(absloss_c_fit$parameters, n = length(absoloss_c_update$val)) -
  SampByPro(absloss_p_fit$parameters, n = length(absoloss_p_update$val))
abar = quantile(additionality, 0.05)
samp_abar = quantile(samp_additionality, 0.05)

##Perform simulations ----
n_rep = 100 #number of repetitions
H = 50 #evaluation horizon; project duration
H_rel = 10 * H #release horizon
D = 0.03 #discount rate
bp = 5 #buffer period
current_year = lubridate::year(lubridate::now()) #pretend every project starts at 2023
scc_current = scc[which(names(scc) == current_year):length(scc)]

sim_p_loss = matrix(0, H, n_rep)
sim_c_loss = matrix(0, H, n_rep)
sim_additionality = matrix(0, H, n_rep)
sim_credit = matrix(0, H, n_rep)
sim_benefit = matrix(0, H, n_rep)
sim_release = matrix(0, H_rel, n_rep)
sim_damage = matrix(0, H, n_rep)
sim_ep = matrix(0, H, n_rep)
sim_credibility = matrix(1, H, n_rep)
sim_r_sched = vector("list", n_rep)

a = Sys.time()
for(j in 1:n_rep){
  to_be_released = rep(0, H)
  r_sched = matrix(0, H, H_rel) #release schedule
  buffer_pool = 0
  for(i in 1:H){
    year_i = t0 + i - 1
    if(year_i <= max(absloss_p_init$year)) {
      sim_p_loss[i, j] = mean(subset(absloss_p_init, year == year_i)$val)
      sim_c_loss[i, j] = mean(subset(absloss_c_init, year == year_i)$val)
      sim_additionality[i, j] = sim_c_loss[i, j] - sim_p_loss[i, j]
      
      if(i <= bp) { #first five years
        if(sim_additionality[i, j] > 0) buffer_pool = buffer_pool + sim_additionality[i, j]
        sim_credit[i, j] = 0
      } else { #sixth year to end of available ex post value
        absloss_p_update = subset(absloss_p_init, year <= year_i)$val
        absloss_c_update = subset(absloss_c_init, year <= year_i)$val
        
        absloss_p_fit = mclust::Mclust(absloss_p_update, 2)
        absloss_c_fit = mclust::Mclust(absloss_c_update, 2)
        
        samp_additionality = SampByPro(absloss_c_fit$parameters, n = length(absloss_c_update)) -
          SampByPro(absloss_p_fit$parameters, n = length(absloss_p_update))
        abar = quantile(samp_additionality, 0.05)
        
        #use buffer pool to fill anticipated releases first
        if(buffer_pool > 0) k = i
        max_release = ifelse(abar > 0, 0, -abar) #if a-bar is positive, maximum release is zero
        while(buffer_pool > 0){
          sim_release[k, j] = min(buffer_pool, max_release)
          buffer_pool = buffer_pool - sim_release[k, j]
          k = k + 1
        }

        sim_credit[i, j] = sim_additionality[i, j] + sim_release[i, j]
        if(sim_credit[i, j] > 0){
          to_be_released[i] = sim_credit[i, j]
          sim_benefit[i] = sim_credit[i] * scc_current[i]
          k = i #number of years after i
          while(to_be_released[i] > 0 & k < H_rel){
            k = k + 1
            if(k > H) {
              max_release = mean(sim_release[, j][which(sim_release[, j] != 0)]) #after project ends, release at average release rate during project
            } else {
              max_release = ifelse(abar > 0, 0, -abar) #if a-bar is positive, maximum release is zero
            }
            can_be_released = max(0, max_release - sim_release[k, j])
            r_sched[i, k] = min(to_be_released[i], can_be_released)
            to_be_released[i] = to_be_released[i] - r_sched[i, k]
            sim_release[k, j] = sim_release[k, j] + r_sched[i, k]
          }
          sim_damage[i, j] = sum(r_sched[i, ] * scc_current[1:length(r_sched[i, ])] / ((1 + D) ^ (1:length(r_sched[i, ]))))
          sim_ep[i, j] = (sim_benefit[i, j] - sim_damage[i, j]) / sim_benefit[i, j]
        } else if(sim_credit[i, j] < 0){
          sim_ep[i, j] = 0
          sim_credibility[i, j] = 0
        }
      }
    } else {
      absloss_p_fit = mclust::Mclust(absloss_p_update, 2)
      absloss_c_fit = mclust::Mclust(absloss_c_update, 2)
      
      sim_p_loss[i, j] = SampByPro(absloss_p_fit$parameters, n = 1)
      sim_c_loss[i, j] = SampByPro(absloss_c_fit$parameters, n = 1)
      absloss_p_update = c(absloss_p_update, sim_p_loss[i, j])
      absloss_c_update = c(absloss_c_update, sim_c_loss[i, j])
      
      sim_additionality[i, j] = sim_c_loss[i, j] - sim_p_loss[i, j]
      sim_credit[i, j] = sim_additionality[i, j] + sim_release[i, j]
      if(sim_credit[i, j] > 0){
        to_be_released[i] = sim_credit[i, j]
        sim_benefit[i] = sim_credit[i] * scc_current[i]
        k = i #number of years after i
        while(to_be_released[i] > 0 & k < H_rel){
          k = k + 1
          if(k > H) {
            max_release = mean(sim_release[, j][which(sim_release[, j] != 0)]) #after project ends, release at average release rate during project
          } else {
            max_release = ifelse(abar > 0, 0, -abar) #if a-bar is positive, maximum release is zero
          }
          can_be_released = max(0, max_release - sim_release[k, j])
          r_sched[i, k] = min(to_be_released[i], can_be_released)
          to_be_released[i] = to_be_released[i] - r_sched[i, k]
          sim_release[k, j] = sim_release[k, j] + r_sched[i, k]
        }
        sim_damage[i, j] = sum(r_sched[i, ] * scc_current[1:length(r_sched[i, ])] / ((1 + D) ^ (1:length(r_sched[i, ]))))
        sim_ep[i, j] = (sim_benefit[i, j] - sim_damage[i, j]) / sim_benefit[i, j]
      } else if(sim_credit[i, j] < 0){
        sim_ep[i, j] = 0
        sim_credibility[i, j] = 0
      }
    }
  sim_r_sched[[j]] = r_sched
  }
}
b = Sys.time()
b - a

sim_df = rbind(data.frame(var = "Project stock", val = subset(summ_stock, year == t0 & series == "mean" & var == "treatment")$val - cumsum(sim_p_loss)),
               data.frame(var = "Counterfactual stock", val = subset(summ_stock, year == t0 & series == "mean" & var == "control")$val - cumsum(sim_c_loss)),
               data.frame(var = "Project loss", val = sim_p_loss),
               data.frame(var = "Counterfactual loss", val = sim_c_loss),
               data.frame(var = "Additionality", val = sim_additionality),
               data.frame(var = "Credit", val = sim_credit),
               data.frame(var = "Anticipated_release", val = sim_release[1:H]),
               data.frame(var = "Credit_benefit", val = sim_benefit),
               data.frame(var = "Damage", val = sim_damage),
               data.frame(var = "EP", val = sim_ep),
               data.frame(var = "Credibility", val = sim_credibility)) %>%
  mutate(year = rep(current_year:(current_year + H - 1), length(unique(var))))
