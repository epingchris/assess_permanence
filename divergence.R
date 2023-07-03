#Adapted from:
#https://medium.com/datalab-log/measuring-the-statistical-similarity-between-two-samples-using-jensen-shannon-and-kullback-leibler-8d05af514b15

library(stats)

compute_probs <- function(data, n) {
  h <- hist(data, breaks = n, plot = FALSE)
  p <- h$counts / length(data)
  e <- h$breaks
  return(list(e = e, p = p))
}

support_intersection <- function(p, q) {
  sup_int <- Filter(function(x) x[[1]] != 0 & x[[2]] != 0, Map(c, p, q))
  #Map(c, p, q): iteratively take each element in p and in q and apply the c() function to create a list of vectors of pairs of probability values
  #Filter(....., Map(c, p, q)), apply the function to each vector, the functioning being a test to see if both values in the vector are non-zero;
  #only keep the vector if it is true
  return(sup_int)
}

get_probs <- function(list_of_tuples) {
  p <- sapply(list_of_tuples, "[[", 1)
  q <- sapply(list_of_tuples, "[[", 2)
  return(list(p = p, q = q))
}

kl_divergence <- function(p, q) {
  return(sum(p * log(p / q)))
}

js_divergence <- function(p, q) {
  m <- (1/2) * (p + q)
  return((1/2) * kl_divergence(p, m) + (1/2) * kl_divergence(q, m))
}

compute_kl_divergence <- function(train_sample, test_sample, n_bins = 15) {
  probs_train <- compute_probs(train_sample, n = n_bins)
  probs_test <- compute_probs(test_sample, n = probs_train$e)
  
  p <- probs_train$p
  q <- probs_test$p
  
  list_of_tuples <- support_intersection(p, q)
  probs <- get_probs(list_of_tuples)
  
  return(kl_divergence(probs$p, probs$q))
}

compute_js_divergence <- function(train_sample, test_sample, bins) {
  probs_train <- compute_probs(train_sample, n = bins)
  probs_test <- compute_probs(test_sample, n = bins)
  
  p <- probs_train$p
  q <- probs_test$p
  
  list_of_tuples <- support_intersection(p, q) #remove intervals where both p and q are zero
  probs <- get_probs(list_of_tuples) #get the new p and q
  
  return(js_divergence(probs$p, probs$q))
}