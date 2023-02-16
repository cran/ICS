# simulations = function(n = 980, n_outlier = 20,  p,
#                        sigma = diag(c(1,rep(4,p-1))), mean =  rep(0,p),
#                        sigma_outlier = diag(c(1,rep(4,p-1))), mean_outlier =  c(6,rep(0,p-1)) ){
#
#   # We simulate the observations from the majority group
#   x = rmvnorm(n = n, mean = mean, sigma = sigma)
#
#   # We simulate the outliers
#   x_outlier = rmvnorm(n = n_outlier, mean = mean_outlier, sigma = sigma_outlier)
#
#   # We create our dataset with the first n_outlier observations as outliers
#   X_ini = rbind( x, x_outlier)
#
#   return(X_ini)
#
# }
# ## Initialization -----
# alpha_all =  c(-1, -0.5, 0.5, 1)
# set.seed(20212022)
#
#
#
# # Applications of ICS for clustering -----
#
# ## Mixture of two normal distributions ----
# ### Simulation of data ----
# p = 4
# Ntot = 10000
# eps = 0.10
# delta = 6
#
# # Simulation of mixture of two Gaussian distributions
# X_ini <- simulations(n = Ntot*(1-eps), n_outlier = Ntot*eps, p = p,
#                      sigma = diag(1,p), mean =  rep(1,p),
#                      sigma_outlier = diag(1,p),
#                      mean_outlier =  c(delta,rep(1,p-1)))
#
# # Simulations of vectors c_k
# k_all <- 0:30
# mean_all_save <- sapply(k_all, function(k){
#   k_sim(k = k, p=p)
# })
# colnames(mean_all_save) <- paste0("K", k_all)
#
# k_sim <- function(k, p){
#   if(k == 0){
#     val <- rep(1,p)
#   }else if (k=="ex"){
#     val <- c(10^(-12), 10^(-3), 1, 10^6)
#   }else{
#     max_power <- round(k/2,0)
#     min_power <- k - max_power
#     all_val <- seq(max_power, -min_power)
#     ind_val <- sample(2:(length(all_val)-1),p-2)
#     # I want the minimum first
#     val <- 10^all_val[c(length(all_val), ind_val, 1)]
#
#   }
#
#   names(val) <- letters[1:p]
#   val
# }
# library(MASS)
# data(crabs)
# # scale (k=6 OK but not k=7)
# mean_all <- k_sim(k =7, p = ncol(X))
# X = apply(as.matrix(crabs[,-c(1:3)]), 2, log)
# X <-  sweep(X, 2, mean_all, "*")
# print(colMeans(X))
# print(format(kappa(X, exact = TRUE), scientific = TRUE,
#              digits = 2))
# MeanCovW <- function (x, alpha, cf){
#   list(location = colMeans(x),
#        scatter = covW(x, alpha = alpha, cf = cf))
# }
# alpha = 1
# # ICS - ICS_EIGEN
# res_ICS_EIGEN <- ICS::ics2(X, S1 = MeanCov,
#                            S2 = MeanCovW, S2args = list(alpha = alpha, cf = 1))
#
#
# ICS(X, S1 = ICS_cov, S2 = ICS_covW, S2_args = list(alpha = 1, cf = 1/(ncol(X)+2)), QR = FALSE)$gKurt
#
# ICS(X, S1 = ICS_cov(X),
#     S2 = MeanCovW(X, alpha = 1, cf = 1/(ncol(X)+2)),QR = FALSE)$gKurt
#
# ICS(X, S1 = ICS_cov, S2 = ICS_covW, S2_args = list(alpha = 1, cf = 1/(ncol(X)+2)), QR = TRUE)$gKurt
#
