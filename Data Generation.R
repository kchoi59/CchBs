############### code for data generation ###############

r <- c()
for (i in 1:2000) {
  cat("Simulation number ", i, "\n")
  n_pop <- 1000 # full cohort size
  beta <- c(log(2),0)
  lambda.zero <- 1
  rho <- 0.5  # correlation between z1 and z2
  p_z <- 2  # number of covariates
  sigma1 <- matrix(c(1, rho, rho, 1), nrow = p_z, ncol = p_z)  # correlation matrix for (z1, z2)
  ### generating covariates (z1, z2) ###
  z_cts <- mvrnorm(n = n_pop, mu = c(0, 0), Sigma = sigma1)
  #z_cts <- ifelse(z_cts > 0, 1, 0)
  colnames(z_cts) <- c("z1", "z2")
  z1 = z_cts[,1] # covariate possibly missing (z_mis)
  z2 = z_cts[,2] # covariate fully observed (z_obs)
  # cor(z1, z2) # correlation
  ### generating failure time T ###
  # generating unif(0,1)
  u <- runif(n_pop, 0, 1)
  # generating failure time T
  time <- -log(1 - u)/(lambda.zero * exp(z_cts %*% beta)) # failure time following exponential distribution
  ### generating censoring time ###
  # generating censoring time C
  censor <- rexp(n_pop, 10.5)
  # failure indicator delta
  delta <- (time <= censor) * 1
  # observed time
  time_obs <- time * delta + censor * (1 - delta)
  # generating ids
  id <- c(1 : n_pop)
  # data: full cohort data
  mydata_full <- data.frame(id, time, censor, time_obs, delta, z_cts)
  r[i] = sum(mydata_full$delta)/nrow(mydata_full)
}
mean(r)
# censoring = 20% | censor <- rexp(n_pop, 0.225) | cor = 0.8/0.5 | continuous
# censoring = 20% | censor <- rexp(n_pop, 0.35) | cor = 0.8/0.5 | binary
# censoring = 90% | censor <- rexp(n_pop, 10.5) | cor = 0.8/0.5 | continuous
# censoring = 90% | censor <- rexp(n_pop, 13) | cor = 0.8/0.5 | binary
