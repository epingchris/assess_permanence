SimulatePermanence = function(type, mean_drawdown = NULL, sites = NULL, aggregate_type = NULL, verbose = F, runtime = T,
                              n_rep = 100, omega = 0.05, H = 50, D = 0.03, warmup = 5, ppr_ratio = 2, scc_df = scc) {
    # A. Get input variables ----
    if(type == "theo") { #theoretical projects (single or aggregated)
        #control input conditions
        if(is.null(aggregate_type)) {
            if(is.null(mean_drawdown)) {
                stop("An aggregated project type or a mean drawdown value(s) is needed.")
            } else {
                cat("Drawdown value(s) used:", mean_drawdown, "\n")
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
                cat("Aggregate project type", aggregate_type, "is used:", mean_drawdown, "\n")
                type_label = paste0(type, "_aggr_", aggregate_type)
            }
        }

        t0 = 2022 #for theoretical projects, simulations are assumed to start at 2022
        t_max = 2021 #for theoretical projects, none of the years will be ex post

        #lambda parameter of the exponential distribution as the inverse of mean drawdown rate
        lambdaP = rep(1, length(mean_drawdown))
        lambdaC = 1 / mean_drawdown

        #observed carbon loss values
        c_loss_p_samp_list = lapply(lambdaP, function(x) rexp(1000, x))
        c_loss_c_samp_list = lapply(lambdaC, function(x) rexp(1000, x))
        obs_p_loss = apply(as.data.frame(c_loss_p_samp_list), 1, sum)
        obs_c_loss = apply(as.data.frame(c_loss_c_samp_list), 1, sum)

        #post-project release rate: a ratio of the counterfactual release rate during project (default to double)
        ppr = sum(ppr_ratio / lambdaC)

        #analytical solution for a_omega (if length of lambdaC is one)
        if(length(mean_drawdown) == 1) {
            aomega = 1 / lambdaP * log(omega * (lambdaP + lambdaC) / lambdaC)
        } else {
            aomega = NULL
        }

    #real-life projects (single or aggregated)
    } else if(type == "real") {
        #control input conditions
        if(is.null(aggregate_type)) {
            if(is.null(sites)) {
                stop("An aggregated project type or a site name(s) is needed.")
            } else {
                cat("Site(s) used:", sites, "\n")
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
            cat("Aggregated project type", aggregate_type, "is used:", sites, "\n")
            type_label = paste0(type, "_aggr_", aggregate_type)
        }

        t0_vec = rep(NA, length(sites))
        t_max_vec = rep(NA, length(sites))
        c_loss_p_list = vector("list", length(sites))
        c_loss_c_list = vector("list", length(sites))

        for(i in seq_along(sites)){
            site_i = sites[i]
            c_flux_df = read.csv(paste0("project_input_data/", site_i, ".csv"), header = T)
            t0_vec[i] = subset(c_flux_df, started)$year[1]
            t_max_vec[i] = max(c_flux_df$year)

            c_loss_p_list[[i]] = c_flux_df %>%
                subset(var == "project" & year >= t0_vec[i]) %>%
                mutate(val = val * (-1), var = NULL, n_sim = NULL, site = site_i)
            c_loss_c_list[[i]] = c_flux_df %>%
                subset(var == "counterfactual" & year >= t0_vec[i]) %>%
                mutate(val = val * (-1), var = NULL, n_sim = NULL, site = site_i)
        }

        t0 = min(t0_vec)
        t_max = min(t_max_vec)
        #min is used to make sure in all years considered to be ex post, there will be available observed values in all projects
        #this may potentially cause some observed values to be discarded, but for now all aggregated projects have same t_max so it isn't an issue

        #GMM-fitted carbon loss distribution, used to generate (similar to lambda parameters for theoretical projects)
        c_loss_p_fit_list = lapply(c_loss_p_list, function(x) FitGMM(x$val))
        c_loss_c_fit_list = lapply(c_loss_c_list, function(x) FitGMM(x$val))

        #observed carbon loss values
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
        ppr = c_loss_c_fit_list %>%
            sapply(function(x) SampGMM(x, n = 1000)) %>%
            apply(1, sum) %>%
            mean() * ppr_ratio

        #no analytical solution for a_omega
        aomega = NULL
    }

    # B. Perform simulations ----
    year_max_scc = max(scc_df$year) #release horizon: number of years into the future for which to calculate release
    H_max_scc = year_max_scc - t0 + 1

    #initialise output matrices
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

    a = Sys.time()
    for(j in 1:n_rep){
        a1 = Sys.time()
        if(verbose) cat("Repetition:", j, "\n")
        schedule = matrix(0, H, H_max_scc) #initialise release schedule
        buffer_pool = 0 #initialise buffer pool

        for(i in 1:H){
            year_i = t0 + i - 1 #get current year
            isExPost = (year_i <= t_max) #check if current year is ex post (observed C loss value is available)

            #get carbon loss value
            if(type == "theo") {
                #for theoretical projects, always sample from exponential distributions
                sim_p_loss[i, j] = sum(sapply(lambdaP, function(x) rexp(1, x)))
                sim_c_loss[i, j] = sum(sapply(lambdaC, function(x) rexp(1, x)))
            } else if(type == "real") {
                #for real-life projects, use ex post (observed) values if available, sample from fitted distributions if not
                sim_p_loss[i, j] = ifelse(isExPost, obs_p_loss[i], sum(sapply(c_loss_p_fit_list, function(x) SampGMM(x, n = 1))))
                sim_c_loss[i, j] = ifelse(isExPost, obs_c_loss[i], sum(sapply(c_loss_c_fit_list, function(x) SampGMM(x, n = 1))))
            }

            #calculate a-omega: sample-based unless in single hypothetical project
            if(type == "theo") {
                #get sampled additionality
                samp_additionality = mapply(function(x, y) rexp(1000, x) - rexp(1000, y),
                                            x = lambdaC, y = lambdaP) %>%
                    as.data.frame() %>%
                    apply(1, function(x) sum(x, na.rm = T))
            } else if(type == "real") {
                #find GMM-fitted C loss distributions using values up to that year
                c_loss_p_fit_list = lapply(c_loss_p_list, function(x) FitGMM(subset(x, year <= year_i & year >= t0)$val))
                c_loss_p_samp_list = lapply(c_loss_p_fit_list, function(x) SampGMM(x, n = 1000))
                
                c_loss_c_fit_list = lapply(c_loss_c_list, function(x) FitGMM(subset(x, year <= year_i & year >= t0)$val))
                c_loss_c_samp_list = lapply(c_loss_c_fit_list, function(x) SampGMM(x, n = 1000))
                
                #get sampled additionality
                samp_additionality = mapply(function(x, y) x - y,
                                            x = c_loss_c_samp_list, y = c_loss_p_samp_list) %>%
                    as.data.frame() %>%
                    apply(1, function(x) sum(x, na.rm = T))
            }
            if(type == "real" | length(mean_drawdown) > 1) aomega = quantile(samp_additionality, omega) #if there is no analytical a_omega, use empirical

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
                    sim_benefit[i, j] = sim_credit[i, j] * filter(scc_df, year == year_i)$value
                    k = i #kth year(s), for which we estimate anticipated release

                    #loop the release process until all creditsa is released
                    while(to_be_released > 0 & k < H_max_scc){
                        k = k + 1
                        if(k > H) {
                            max_release = ppr #post-project release rate
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
                    sim_damage[i, j] = sum(schedule[i, seq_along(years_k) + i] * filter(scc_df, year %in% years_k)$value / ((1 + D) ^ seq_along(years_k)))
                    sim_ep[i, j] = (sim_benefit[i, j] - sim_damage[i, j]) / sim_benefit[i, j]
                    #cat("Credit =", sim_credit[i, j], ", Benefit =", sim_benefit[i, j], ", Damage =", sim_damage[i, j], ", eP =", sim_ep[i, j], "\n")
                } else if(sim_credit[i, j] <= 0){
                    #net credit is negative; reversal happens
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
        if(runtime) cat("Repetition", j, "runtime:", b1 - a1, "\n")
    }
    b = Sys.time()
    if(runtime) cat("Total runtime:", b - a, "\n")

    sim_ep = sim_ep %>%
        replace(., . == Inf| . == -Inf| . == 0, NA)

    #gather common output variables
    if(type == "theo") {
        common_var = list(type = type, type_label = type_label, t0 = t0, t_max = t_max,
                          n_rep = n_rep, omega = omega, H = H, D = D, warmup = warmup, ppr_ratio = ppr_ratio,
                          mean_drawdown = mean_drawdown)
    } else if(type == "real"){
        common_var = list(type = type, type_label = type_label, t0 = t0, t_max = t_max,
                          n_rep = n_rep, omega = omega, H = H, D = D, warmup = warmup, ppr_ratio = ppr_ratio,
                          sites = sites, c_loss_p_list = c_loss_p_list, c_loss_c_list = c_loss_c_list)
    }

    output_list = list(common_var = common_var,
                       sim_p_loss = sim_p_loss, sim_c_loss = sim_c_loss, sim_additionality = sim_additionality, sim_credit = sim_credit,
                       sim_benefit = sim_benefit, sim_aomega = sim_aomega, sim_release = sim_release, sim_damage = sim_damage,
                       sim_ep = sim_ep, sim_pact = sim_pact, sim_reversal = sim_reversal, sim_buffer = sim_buffer, sim_schedule = sim_schedule)
    return(output_list)
}