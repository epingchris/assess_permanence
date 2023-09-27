#Adapted from:
#https://medium.com/datalab-log/measuring-the-statistical-similarity-between-two-samples-using-jensen-shannon-and-kullback-leibler-8d05af514b15

library(stats)

kl_divergence <- function(p, q) {
  probs_df = data.frame(p = p, q = q, intersect = ifelse(p * q == 0, F, T))
  p_intersect = subset(probs_df, intersect)$p
  q_intersect = subset(probs_df, intersect)$q
  
  return(sum(p_intersect * log(p_intersect / q_intersect)))
}

js_divergence <- function(p, q) {
  m <- (1/2) * (p + q)
  return((1/2) * kl_divergence(p, m) + (1/2) * kl_divergence(q, m))
}

compute_kl_divergence <- function(a, b, n_bins = 15) {
  probs_a <- hist(a, breaks = bins, plot = FALSE)$counts / length(a)
  probs_b <- hist(b, breaks = bins, plot = FALSE)$counts / length(b)
  
  return(kl_divergence(probs_a, probs_b))
}

compute_js_divergence <- function(a, b, bins) {
  probs_a <- hist(a, breaks = bins, plot = FALSE)$counts / length(a)
  probs_b <- hist(b, breaks = bins, plot = FALSE)$counts / length(b)
  
  return(js_divergence(probs_a, probs_b))
}

JSTest = function(a, b, bins, alpha = 0.05, nperm = 10000){
  js_obs = compute_js_divergence(a, b, bins = bins)
  pooled = c(a, b)
  js_perm = rep(NA, nperm)
  for(i in 1:nperm){
    samp_ind = sample(1:(length(a) + length(b)), length(a))
    samp_a = pooled[samp_ind]
    samp_b = pooled[-samp_ind]
    js_perm[i] = compute_js_divergence(samp_a, samp_b, bins = bins)
  }
  js_crit = quantile(js_perm, 1 - alpha)
  p_val = length(which(js_perm > js_obs)) / nperm
  if(p_val == 0) {p_val = paste("<", 1 / nperm)}
  return(list(js_obs = js_obs, js_crit = js_crit, p_val = p_val))
}
