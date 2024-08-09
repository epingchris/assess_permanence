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
SummariseSim = function(mat, n_rep = 100){
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

#function to summarise across simulations with a series of different input variable values, used for sensitivity tests
SummariseAcrossTests = function(input, var, n_rep) {
  n = nrow(input)
  output = input %>%
    t() %>%
    as.data.frame() %>%
    rowwise() %>%
    mutate(mean = mean(c_across(1:n_rep)),
           sd = sd(c_across(1:n_rep))) %>%
    ungroup() %>%
    mutate(var = var,
           ci_margin = qt(0.975, df = n_rep - 1) * sd / sqrt(n_rep),
           ci_low = mean - ci_margin,
           ci_high = mean + ci_margin) %>%
    select(-all_of(c(1:n)))
  return(output)
}
