#' This function are used for stouffer zscore
#' @param w is the weight of the value
#' @param z is the value itself(the number of value is same as w)
#' @return stoufferzscore value
#' @export
# Stoufferzscore <- function(w,z){
#   w <- w/max(abs(w))
#   weighted_sum <- sum(w*z,na.rm = T)
#   Sqrt_Sumsquares <- sqrt(sum(w^2,na.rm = T))
#   Z <- weighted_sum/Sqrt_Sumsquares
#   return(Z)
# }

Stoufferzscore <- function(w,z){
  w <- w/max(abs(w))
  weighted_sum <- sum(w*z,na.rm = T)
  gene_number <- sum(!is.na(z))
  Z <- weighted_sum/gene_number*100
  return(Z)
}
