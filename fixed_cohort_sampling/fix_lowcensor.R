######### Subcohort Sampling | low censoring | fixed cohort | continuous covariates | rho 0.8 #########


rm(list = ls())

library(survival)
library(survey)
library(sampling)
library(addhazard)
library(MASS)
library(gustave)
library(pastecs)
library(psych)
library(writexl)
library(readxl)
library(dplyr)
library(segmented)

seed_num = 2222
set.seed(seed_num)

n_pop <- 1000 # full cohort size 
beta <- c(log(2),0)
lambda.zero <- 1
rho <- 0.8  # correlation between z1 and z2
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
censor <- rexp(n_pop, 0.225)
# failure indicator delta
delta <- (time <= censor) * 1
# observed time
time_obs <- time * delta + censor * (1 - delta)
# generating ids
id <- c(1 : n_pop) 
# data: full cohort data  
mydata_full <- data.frame(id, time, censor, time_obs, delta, z_cts)

beta0 = coef(coxph(Surv(time_obs, delta) ~ z1, 
           data = mydata_full))
beta0

################################ score ######################################

# function for calculating score vector: is it correct?
compute_score_vector <- function(mydata, beta, w) { 
  # mydata: dataset
  # beta: beta 
  # w: sampling weights
  
  score <- 0 # initial score
  # risk set
  for (i in 1 : sum(mydata$delta)) {
    t <- sort(mydata$time_obs[mydata$delta == 1])[i] # sort the observed time
    idx <- which(mydata$time_obs == t) # indices of subjects with observed event at t
    n_tie <- length(idx) # number of subjects with tied event times
    risk_id <- which(mydata$time_obs >= t) # risk set id
    s0 <- sum(w[risk_id] * exp(beta * mydata$z1[risk_id])) # denominator of Z_bar
    s1 <- sum(w[risk_id] * mydata$z1[risk_id] * exp(beta * mydata$z1[risk_id])) # numerator of Z_bar
    
    # calculate the score vector
    for (j in 1 : n_tie) {
      score <- score + w[idx][j] * mydata$delta[idx][j] * (mydata$z1[idx][j] - (s1 / s0))
      # summation over events
    }
  }
  return(score)
}



# read the data
n = 200 # subcohort size
N = nrow(mydata_full) # cohort size
ifmodel = coxph(Surv(time_obs, delta) ~ z2, data = mydata_full) 
inffun = resid(ifmodel, "dfbeta")   # delta-beta: add 1 or not?
mydata_full = cbind(mydata_full, inffun) 
mydata_full$pik = rep(n/N, N); 
mydata_full$wt = 1/mydata_full$pik


sim = 10
beta_srs <- c()
beta_bs <- c()
se_srs <- c()
se_bs <- c()
D1 <- c()
D2 <- c()
est_se_score <- c()
est_se_dfbeta1 <- c()
est_se_dfbeta2 <- c()
vub_bs_dfbeta1 <- c()
vub_bs_dfbeta2 <- c()
vub_bs_score <- c()
vub_bs_dfbeta <- c()
robust_se_bs <- c()
robust_se_srs <- c()
beta_cali <- c()
se_cali <- c()
robust_se_cali <- c()
beta_srs_design <- c()
se_srs_design <- c()
robust_se_srs_design <- c()
vub_bs2 <- c()
est_se_dfbeta <- c()
beta_bs_design <- c()
se_bs_design <- c()
beta_bs_cali <- c()
se_bs_cali <- c()
robust_se_bs_cali <- c()
ci_up_score <- c()
ci_low_score <- c()
ci_ind_score <- c()

for(i in 1 : sim) {
  cat("Simulation number ", i, "\n")
  
  ### srs ###
  srs.subcohort <- with(mydata_full, sample(1 : N, n))
  mydata_full$in.srs.subcohort <- (1:N) %in% srs.subcohort
  in.srs.subcohort <- mydata_full$in.srs.subcohort
  cox_srs = coxph(Surv(time_obs, delta) ~ z1, 
                  data = mydata_full[mydata_full$in.srs.subcohort,],
                  weights = mydata_full$wt[mydata_full$in.srs.subcohort], 
                  robust = T)
  
  beta_srs[i] <- coef(cox_srs)
  se_srs[i] <- summary(cox_srs)$coefficients[, "se(coef)"]
  robust_se_srs[i] <- summary(cox_srs)$coefficients[, "robust se"]
  
  ### bs ###
  s <- samplecube(cbind(mydata_full$inffun,mydata_full$pik), 
                  mydata_full$pik, order = 1, comment = F)
  mydata_full$in.bs.subcohort <- as.logical(s)
  in.bs.subcohort <- mydata_full$in.bs.subcohort
  cox_bs = coxph(Surv(time_obs, delta) ~ z1, 
                 data = mydata_full[mydata_full$in.bs.subcohort,],
                 weights = mydata_full$wt[mydata_full$in.bs.subcohort], 
                 robust = T)
  
  beta_bs[i] <- coef(cox_bs)
  se_bs[i] <- summary(cox_bs)$coefficients[, "se(coef)"] 
  robust_se_bs[i] <- summary(cox_bs)$coefficients[, "robust se"]
  
  ### cali ###
  cch_design <- twophase(id = list(~ 1, ~ 1), 
                         subset = ~ in.srs.subcohort,
                         weights = list(NULL, ~ wt),
                         data= mydata_full, method = "approx")
  cox_cch <- svycoxph(Surv(time_obs, delta) ~ z1, design = cch_design)
  
  beta_srs_design[i] <- coef(cox_cch)
  se_srs_design[i] <- summary(cox_cch)$coefficients[, "se(coef)"] 
  
  cch_cal <- calibrate(cch_design, phase = 2, calfun = "raking", ~ inffun)
  cox_cal_cch <- svycoxph(Surv(time_obs, delta) ~ z1, design = cch_cal)
  
  beta_cali[i] <- coef(cox_cal_cch)
  se_cali[i] <- summary(cox_cal_cch)$coefficients[, "se(coef)"]
  robust_se_cali[i] <- summary(cox_cal_cch)$coefficients[, "robust se"]
  
  ### bs cali ###
  bs_design <- twophase(id = list(~ 1, ~ 1), 
                        subset = ~ in.bs.subcohort,
                        weights = list(NULL, ~ wt),
                        data= mydata_full, method = "approx")
  cox_bs_design <- svycoxph(Surv(time_obs, delta) ~ z1, design = bs_design)
  
  beta_bs_design[i] <- coef(cox_bs_design)
  se_bs_design[i] <- summary(cox_bs_design)$coefficients[, "se(coef)"] 
  
  bs_cal <- calibrate(bs_design, phase = 2, calfun = "raking", ~ inffun + pik)
  cox_cal_bs <- svycoxph(Surv(time_obs, delta) ~ z1, design = bs_cal)
  
  beta_bs_cali[i] <- coef(cox_cal_bs)
  se_bs_cali[i] <- summary(cox_cal_bs)$coefficients[, "se(coef)"]
  robust_se_bs_cali[i] <- summary(cox_cal_bs)$coefficients[, "robust se"]
  
  
  
  ### estimate the bs var
  ## BS -  var est ##
  # BS variance estimator
  var_u_b <- function(Ys, Xs, pik, p, n) {
    c <- (1 - pik) * n / (n - p)
    b <- solve(t(Xs) %*% diag(c / pik^2) %*% Xs, t(Xs) %*% diag(c / pik^2) %*% Ys)
    vub <- c(t(Ys) %*% diag(c / pik^2) %*% Ys - t(Ys) %*% diag(c / pik^2) %*% Xs %*% b - 
               t(b) %*% t(Xs) %*% diag(c / pik^2) %*% Ys + t(b) %*% t(Xs) %*% diag(c / pik^2) %*% Xs %*% b)
    return(vub)
  }
  
  
  mydata_cch <- mydata_full[mydata_full$in.bs.subcohort,]
  mydata_cch$cch_id <- 1 : nrow(mydata_cch)
  cch_sub_id <- mydata_cch[mydata_cch$in.bs.subcohort == T,]$cch_id
  Ys_bs_dfbeta = resid(cox_bs, "dfbeta")[cch_sub_id] * (n/N)
  Xs_bs = mydata_cch$inffun[cch_sub_id] 
  pik_bs = mydata_cch$pik[cch_sub_id] 
  
  
  vub_bs2[i] = var_u_b(Ys_bs_dfbeta, 
                       cbind(mydata_cch$pik[cch_sub_id],
                             mydata_cch$inffun[cch_sub_id]), 
                       mydata_cch$pik[cch_sub_id], 2, sum(s))   # est se using dfbeta
  
  est_se_dfbeta[i] = sqrt(vub_bs2[i]) 
  
  
  ci_up_score[i] <- beta_bs[i] + qnorm(0.975) * est_se_dfbeta[i]
  ci_low_score[i] <- beta_bs[i] - qnorm(0.975) * est_se_dfbeta[i] 
  ci_ind_score[i] <- ifelse((beta0 >= ci_low_score[i]) & (beta0 <= ci_up_score[i]), 1, 0)
  
  
  
}

rr=data.frame(v=c("mean of beta srs","sd of beta srs",
                  "mean of beta srs design","sd of beta srs design",
                  "mean of beta bs","sd of beta bs",
                  "mean of beta bs design","sd of beta bs design",
                  "mean of beta cali", "sd of beta cali",
                  "mean of beta bs cali", "sd of beta bs cali",
                  "mean of se srs","sd of se srs",
                  "mean of se srs design","sd of se srs design",
                  "mean of se bs","sd of se bs",
                  "mean of se cali","sd of se cali",
                  "mean of robust se srs","sd of robust se srs",
                  "mean of robust se bs","sd of robust se bs",
                  "mean of robust se cali","sd of robust se cali",
                  "mean of est se bs", "sd of est se bs",
                  "CI %"),
              n=c(mean(beta_srs),sd(beta_srs), 
                  mean(beta_srs_design),sd(beta_srs_design), 
                  mean(beta_bs),sd(beta_bs),
                  mean(beta_bs_design),sd(beta_bs_design),
                  mean(beta_cali),sd(beta_cali),
                  mean(beta_bs_cali),sd(beta_bs_cali),
                  mean(se_srs),sd(se_srs),
                  mean(se_srs_design),sd(se_srs_design),
                  mean(se_bs),sd(se_bs),
                  mean(se_cali),sd(se_cali),
                  mean(robust_se_srs),sd(robust_se_srs),
                  mean(robust_se_bs),sd(robust_se_bs),
                  mean(robust_se_cali),sd(robust_se_cali),
                  mean(est_se_dfbeta),sd(est_se_dfbeta),
                  sum(ci_ind_score)/length(ci_ind_score)))

(rr=data.frame(n=rr$v,n=round(rr$n,4)))
