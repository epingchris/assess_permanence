makeFlux = function(in_df){
  n_sim = unique(in_df$n_sim)
  flux_series = lapply(seq_along(n_sim), function(i) {
    stock_series = aggregate(class_co2e ~ treatment + year, subset(in_df, n_sim == n_sim[i] & class %in% c(1, 2, 3, 4)), FUN = sum)
    stock_wide = as.data.frame(pivot_wider(stock_series, id_cols = year, names_from = "treatment", values_from = "class_co2e"))
    
    flux_series = data.frame(year = stock_wide$year[-1],
                             treatment_proj = diff(subset(stock_series, treatment == "treatment")$class_co2e),
                             control_proj = diff(subset(stock_series, treatment == "control")$class_co2e)) %>%
      mutate(additionality = treatment_proj - control_proj) %>%
      pivot_longer(2:4, names_to = "var", values_to = "val") %>%
      mutate(series = paste0("V", i))
  }) %>%
    do.call(rbind, .)
  return(flux_series)
}

#wrapper function to fit using GMM and generate random samples from GMM-fitted distributions
FitGMM = function(x) {
  if(length(x) == 0) {
    return(NULL)
  }
  
  mclust::Mclust(x, 2, verbose = F)
}

SampGMM = function(mclust_obj, n) {
  if(is.null(mclust_obj)) {
    return(NA)
  }
  
  params = mclust_obj$param
  if(length(params$variance$sigmasq) == 1) params$variance$sigmasq[2] = params$variance$sigmasq[1]
  samp_vec = NULL
  while(length(samp_vec) < n) {
    distr_chosen = ifelse(runif(1) <= params$pro[1], 1, 2)
    val = rnorm(1, mean = params$mean[distr_chosen], sd = sqrt(params$variance$sigmasq[distr_chosen]))
    if(val > 0) samp_vec = c(samp_vec, val)
  }
  return(samp_vec)
}

#create summary of time series where number of row is H
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

