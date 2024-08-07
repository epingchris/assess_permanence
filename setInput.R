setInput = function(type, mean_drawdown = NULL, sites = NULL, aggregate_type = NULL) {
  # Theoretical projects (single or aggregated) ----
  if(type == "theo") {
    if(is.null(aggregate_type)) {
      if(is.null(mean_drawdown)) {
        stop("An aggregated project type or a mean drawdown value(s) is needed.")
      } else {
        cat("Drawdown value(s) used:", mean_drawdown)
        type_label = paste0(type, "_", paste(gsub("\\.", "_", mean_drawdown), collapse = "_"))
      }
    } else {
      if(aggregate_type %in% c("A", "B", "C") == F) {
        stop("Aggregated project type is not defined.")
      } else {
        mean_drawdown = switch(aggregate_type,
                               "A" = c(1.1, 1.1, 1.1, 5),
                               "B" = c(1.1, 1.1, 5, 5),
                               "C" = c(1.1, 5, 5, 5))
        cat("Aggregate project type", aggregate_type, "is used:", mean_drawdown)
        type_label = paste0(type, "_aggr_", aggregate_type)
      }
    }

    t0 = 2021
    t_max = t0 - 1 #actually no ex post data will be used for hypothetical projects
    lambdaP = rep(1, length(mean_drawdown))
    lambdaC = 1 / mean_drawdown
    
    c_loss_p_samp_list = lapply(lambdaP, function(x) rexp(1000, x))
    c_loss_c_samp_list = lapply(lambdaC, function(x) rexp(1000, x))
    obs_p_loss = apply(as.data.frame(c_loss_p_samp_list), 1, sum)
    obs_c_loss = apply(as.data.frame(c_loss_c_samp_list), 1, sum)
    
    #post-project release rate: a ratio of the counterfactual release rate during project (default to double)
    postproject_release = sum(postproject_ratio / lambdaC)
    
    #analytical solution for a_omega (if length of lambdaC is one)
    if(length(lambdaC) == 1) {
      aomega = 1 / lambdaP * log(omega * (lambdaP + lambdaC) / lambdaC)
    } else {
      aomega = NULL
    }

    output_list = list(type = type, mean_drawdown = mean_drawdown, type_label = type_label,
                       t0 = t0, t_max = t_max, lambdaP = lambdaP, lambdaC = lambdaC,
                       obs_p_loss = obs_p_loss, obs_c_loss = obs_c_loss,
                       postproject_release = postproject_release, aomega = aomega)
    
  # Real-life projects (single or aggregated) ----
  } else if(type == "real") {
    if(is.null(aggregate_type)) {
      if(is.null(sites)) {
        stop("An aggregated project type or a site name(s) is needed.")
      } else {
        cat("Site(s) used:", sites)
        sites_simplified = sapply(sites, function(x) {
          switch(x,
                 "Gola_country" = "Gola", #1201
                 "WLT_VNCC_KNT" = "KNT",
                 "CIF_Alto_Mayo" = "Alto_Mayo", #944
                 "VCS_1396" = "RPA",
                 "VCS_934" = "Mai_Ndombe")
        }) %>% as.vector()
        type_label = paste0(type, "_", paste(sites_simplified, collapse = "_"))
      }
    } else if(aggregate_type %in% c("four", "three") == F) {
      stop("Aggregated project type is not defined.")
    } else {
      sites = switch(aggregate_type,
                     "four" = c("Gola_country", "CIF_Alto_Mayo", "VCS_1396", "VCS_934"),
                     "three" = c("Gola_country", "CIF_Alto_Mayo", "VCS_1396"))
      cat("Aggregated project type", aggregate_type, "is used:", sites)
      type_label = paste0(type, "_aggr_", aggregate_type)
    }
    
    t0_vec = rep(NA, length(sites))
    tmax_vec = rep(NA, length(sites))
    flux_series_list = vector("list", length(sites))
    c_loss_p_list = vector("list", length(sites))
    c_loss_c_list = vector("list", length(sites))
    
    for(i in seq_along(sites)){
      site_i = sites[i]
      agb_series_df = read.csv(paste0("project_input_data/", site_i, ".csv"), header = T)
      t0_vec[i] = subset(agb_series_df, started)$year[1]
      tmax_vec[i] = max(agb_series_df$year)
      flux_series_list[[i]] = makeFlux(agb_series_df)

      c_loss_p_list[[i]] = flux_series_list[[i]] %>%
        subset(var == "treatment_proj" & year >= t0_vec[i]) %>%
        mutate(val = val * (-1), var = NULL, n_sim = NULL, site = site_i)
      c_loss_c_list[[i]] = flux_series_list[[i]] %>%
        subset(var == "control_proj" & year >= t0_vec[i]) %>%
        mutate(val = val * (-1), var = NULL, n_sim = NULL, site = site_i)
    }
    
    t0 = min(t0_vec)
    t_max = 2021
    
    obs_p_loss = c_loss_p_list %>%
      do.call(rbind, .) %>%
      group_by(year, site) %>%
      summarise(val = mean(val)) %>%
      ungroup(site) %>%
      summarise(val = sum(val), .groups = "drop") %>%
      pull(val)
    obs_c_loss = c_loss_c_list %>%
      do.call(rbind, .) %>%
      group_by(year, site) %>%
      summarise(val = mean(val)) %>%
      ungroup(site) %>%
      summarise(val = sum(val), .groups = "drop") %>%
      pull(val)
    
    #post-project release rate: a ratio of the counterfactual release rate during project (default to double)
    c_loss_p_fit_list = lapply(c_loss_p_list, function(x) FitGMM(x$val))
    c_loss_c_fit_list = lapply(c_loss_c_list, function(x) FitGMM(x$val))
    postproject_release = c_loss_c_fit_list %>%
      sapply(function(x) SampGMM(x, n = 1000)) %>%
      apply(1, sum) %>%
      mean() * postproject_ratio
    
    #input for simulation: t0, t_max, c_loss_p/c_init_list, obs_lossP/C, postproject_release
    
    output_list = list(type = type, sites = sites, type_label = type_label,
                       t0 = t0, t_max = t_max, c_loss_p_list = c_loss_p_list, c_loss_c_list = c_loss_c_list,
                       obs_p_loss = obs_p_loss, obs_c_loss = obs_c_loss,
                       postproject_release = postproject_release, aomega = NULL)
  }
  return(output_list)
}