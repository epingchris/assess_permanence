setInput = function(type, mean_drawdown = NULL, drawdown_type = NULL, sites = NULL, sites_type = NULL) {
  ## A. Hypothetical projects ----
  if(type == "hypo") {
    if(is.null(mean_drawdown)) {
      stop("Theoretical projects need a mean drawdown value.")
    }
    if(length(mean_drawdown) > 1) {
      stop("Type \"hypo_aggr\" needs to be used for aggregated projects.")
    }

    t0 = 2021
    
    t_max = t0 - 1 #actually no ex post data will be used for hypothetical projects
    
    lambdaP = 1
    lambdaC = 1 / mean_drawdown
    
    obs_p_loss = rexp(1000, lambdaP)
    obs_c_loss = rexp(1000, lambdaC)
    
    #post-project release rate: a ratio of the counterfactual release rate during project (default to double)
    postproject_release = postproject_ratio / lambdaC
    
    aomega = 1 / lambdaP * log(omega * (lambdaP + lambdaC) / lambdaC) #analytical solution
    
    output_list = list(t0 = t0, t_max = t_max, lambdaP = lambdaP, lambdaC = lambdaC,
                       obs_p_loss = obs_p_loss, obs_c_loss = obs_c_loss,
                       postproject_release = postproject_release, aomega = aomega)
    
    ## B. Real-life projects ----
  } else if(type == "real") {
    if(is.null(sites)) {
      stop("Real projects need a site name value.")
    }
    if(length(sites) > 1) {
      stop("Type \"real_aggr\" needs to be used for aggregated projects.")
    }

    agb_series_df = read.csv(paste0(file_path, "project_input_data/", sites, ".csv"), header = T)
    t0 = subset(agb_series_df, started)$year[1]
    t_max = max(agb_series_df$year)
    flux_series = makeFlux(agb_series_df)
    
    absloss_p_init = flux_series %>%
      subset(var == "treatment_proj" & year >= t0) %>%
      mutate(val = val * (-1), var = NULL, series = NULL)
    absloss_c_init = flux_series %>%
      subset(var == "control_proj" & year >= t0) %>%
      mutate(val = val * (-1), var = NULL, series = NULL)
    
    obs_p_loss = absloss_p_init %>%
      group_by(year) %>%
      summarise(val = mean(val), .groups = "drop") %>%
      pull(val)
    obs_c_loss = absloss_c_init %>%
      group_by(year) %>%
      summarise(val = mean(val), .groups = "drop") %>%
      pull(val)
    
    #post-project release rate: a ratio of the counterfactual release rate during project (default to double)
    absloss_p_fit = FitGMM(absloss_p_init$val)
    absloss_c_fit = FitGMM(absloss_c_init$val)
    postproject_release = mean(SampGMM(absloss_c_fit, n = 1000)) * postproject_ratio
    
    output_list = list(t0 = t0, t_max = t_max, absloss_p_init = absloss_p_init, absloss_c_init = absloss_c_init,
                       obs_p_loss = obs_p_loss, obs_c_loss = obs_c_loss,
                       postproject_release = postproject_release, aomega = NULL)
    
  ## A. Theoretical projects (single or aggregated) ----
  } else if(type == "theo") {
    if(is.null(drawdown_type)) {
      if(is.null(mean_drawdown)) {
        stop("Aggregated theoretical projects need a type or a vector of mean drawdown values.")
      } else {
        cat("Drawdown value(s) used:", mean_drawdown)
      }
    } else {
      if(drawdown_type %in% c("A", "B", "C") == F) {
        stop("Drawdown type is not defined.")
      } else {
        mean_drawdown = switch(drawdown_type,
                               "A" = c(1.1, 1.1, 1.1, 5),
                               "B" = c(1.1, 1.1, 5, 5),
                               "C" = c(1.1, 5, 5, 5))
        cat("Drawdown type", drawdown_type, "is used:", mean_drawdown)
      }
    }

    t0 = 2021
    t_max = t0 - 1 #actually no ex post data will be used for hypothetical projects
    lambdaP_vec = rep(1, length(mean_drawdown))
    lambdaC_vec = 1 / mean_drawdown
    
    absloss_p_samp_list = lapply(lambdaP_vec, function(x) rexp(1000, x))
    absloss_c_samp_list = lapply(lambdaC_vec, function(x) rexp(1000, x))
    obs_p_loss = apply(as.data.frame(absloss_p_samp_list), 1, sum)
    obs_c_loss = apply(as.data.frame(absloss_c_samp_list), 1, sum)
    
    #post-project release rate: a ratio of the counterfactual release rate during project (default to double)
    postproject_release = sum(postproject_ratio / lambdaC_vec)
    
    #analytical solution for a_omega (if length of lambdaC is one)
    if(length(lambdaC_vec) == 1) {
      aomega = 1 / lambdaP_vec * log(omega * (lambdaP_vec + lambdaC_vec) / lambdaC_vec)
    } else {
      aomega = NULL
    }

    output_list = list(t0 = t0, t_max = t_max, lambdaP_vec = lambdaP_vec, lambdaC_vec = lambdaC_vec,
                       obs_p_loss = obs_p_loss, obs_c_loss = obs_c_loss,
                       postproject_release = postproject_release, aomega = aomega)
    
    ## D. Aggregated real-life projects ----
  } else if(type == "real_aggr") {
    if(is.null(sites)) {
      if(is.null(sites_type)) {
        stop("For aggregated real-life projects, an aggregated site type or a vector of site names is needed.")
      }
    } else if(sites_type %in% c("four", "three") == F) {
      stop("Aggregated project type not defined.")
    } else {
      sites = switch(sites_type,
                     "four" = c("Gola_country", "CIF_Alto_Mayo", "VCS_1396", "VCS_934"),
                     "three" = c("Gola_country", "CIF_Alto_Mayo", "VCS_1396"))
    }
    
    flux_series_list = vector("list", length(sites))
    absloss_p_init_list = vector("list", length(sites))
    absloss_c_init_list = vector("list", length(sites))
    t0_vec = rep(NA, length(sites))
    
    for(i in seq_along(sites)){
      site_i = sites[i]
      agb_series_df = read.csv(paste0(file_path, "project_input_data/", site_i, ".csv"), header = T)
      t0_vec[i] = subset(agb_series_df, started)$year[1]
      flux_series = makeFlux(agb_series_df)

      flux_series_list[[i]] = lapply(seq_along(agb_series_project_sim), makeFlux) %>%
        do.call(rbind, .)
      absloss_p_init_list[[i]] = flux_series_list[[i]] %>%
        subset(var == "treatment_proj" & year >= t0) %>%
        mutate(val = val * (-1), var = NULL, series = NULL, site = site_i)
      absloss_c_init_list[[i]] = flux_series_list[[i]] %>%
        subset(var == "control_proj" & year >= t0) %>%
        mutate(val = val * (-1), var = NULL, series = NULL, site = site_i)
    }
    
    site = "portfolio"
    t0 = min(t0_vec)
    t_max = 2021
    
    obs_p_loss = absloss_p_init_list %>%
      do.call(rbind, .) %>%
      group_by(year, site) %>%
      summarise(val = mean(val)) %>%
      ungroup(site) %>%
      summarise(val = sum(val), .groups = "drop") %>%
      pull(val)
    obs_c_loss = absloss_c_init_list %>%
      do.call(rbind, .) %>%
      group_by(year, site) %>%
      summarise(val = mean(val)) %>%
      ungroup(site) %>%
      summarise(val = sum(val), .groups = "drop") %>%
      pull(val)
    
    #post-project release rate: a ratio of the counterfactual release rate during project (default to double)
    absloss_p_fit_list = lapply(absloss_p_init_list, function(x) FitGMM(x$val))
    absloss_c_fit_list = lapply(absloss_c_init_list, function(x) FitGMM(x$val))
    postproject_release = absloss_c_fit_list %>%
      sapply(function(x) SampGMM(x, n = 1000)) %>%
      apply(1, sum) %>%
      mean() * postproject_ratio
    
    #input for simulation: t0, t_max, absloss_p/c_init_list, obs_lossP/C, postproject_release
    
    output_list = list(t0 = t0, t_max = t_max, absloss_p_init_list = absloss_p_init_list, absloss_c_init_list = absloss_c_init_list,
                       obs_p_loss = obs_p_loss, obs_c_loss = obs_c_loss,
                       postproject_release = postproject_release, aomega = NULL)
  }
  return(output_list)
}