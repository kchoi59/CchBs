###########################################################################################################
######### Case Cohort Sampling | high censoring | random cohort | continuous covariates | rho 0.8 #########
###########################################################################################################

library(survival)
library(survey)
library(sampling)
library(addhazard)
library(MASS)
library(gustave)
library(pastecs)
library(writexl)
library(readxl)
library(dplyr)

set.seed(111)
sim = 2000
beta_srs <- c()
beta_bs <- c()
se_srs <- c()
se_bs <- c()
D2 <- c()
vub_bs_dfbeta2 <- c()
est_se_dfbeta2 <- c()
est_se_D2_df <- c()
full_beta <- c()
phase1_var <- c()
phase2_var <- c()
est_se_two <- c()
est_se_two2 <- c()
beta_cali <- c()
beta_bs_cali <- c()
# D1 <- c()
# D3 <- c()
# D4 <- c()
# robust_se_bs <- c()
# robust_se_srs <- c()
# full_se <- c()
# full_robust_se <- c()
# se_cali <- c()
# robust_se_cali <- c()
# se_bs_cali <- c()
# robust_se_bs_cali <- c()
# ci_up_score <-c()
# ci_low_score <-c()
# ci_ind_score <-c()


for(i in 1 : sim) {
  cat("Simulation number ", i, "\n")
  
  ####################################################################
  ####################### generate full cohort #######################
  ####################################################################
  
  n_pop <- 1000 # full cohort size
  beta <- c(log(2), 0)
  lambda.zero <- 1
  rho <- 0.8  # correlation between z1 and z2
  p_z <- 2  # number of covariates
  sigma1 <- matrix(c(1, rho, rho, 1), nrow = p_z, ncol = p_z)  # correlation matrix for (z1, z2)
  ### generating covariates (z1, z2) ###
  z_cts <- mvrnorm(n = n_pop, mu = c(0, 0), Sigma = sigma1)
  # z_cts <- ifelse(z_cts > 0, 1, 0)
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
  
  ####################################################################
  ############################ full cohort ###########################
  ####################################################################
  mydata <- data.frame(id, time, censor, time_obs, delta, z_cts)
  full_cox <- coxph(Surv(time_obs, delta) ~ z1, data = mydata, robust = T)
  full_beta[i] <- coef(full_cox)
  # full_se[i] <- summary(full_cox)$coefficients[, "se(coef)"]
  # full_robust_se[i] <- summary(full_cox)$coefficients[, "robust se"]
  
  
  n <- 100; # subcohort size
  N <- nrow(mydata) # full cohort size
  ifmodel <- coxph(Surv(time_obs, delta) ~ z2, data = mydata)
  inffun <- resid(ifmodel, "dfbeta")
  mydata_if <- cbind(mydata, inffun)
  mydata_if$pik = ifelse(mydata_if$delta == 1, 1, n/(N - sum(mydata_if$delta))) # inclusion probability...!
  mydata_if$wt = 1/mydata_if$pik
  mydata_control <- mydata_if[mydata_if$delta == 0, ]
  
  ####################################################################
  ####################################################################
  ####################### case-cohort sampling #######################
  ####################################################################
  ####################################################################
  
  
  ####################################################################
  ###################### simple random sampling ######################
  ####################################################################
  casectrl <- with(mydata_if,c(which(delta==1),sample(which(delta==0),n)))
  mydata_if$in.ccs <- (1:nrow(mydata_if)) %in% casectrl
  cox_ccs_srs <- coxph(Surv(time_obs, delta) ~ z1, data = mydata_if[mydata_if$in.ccs,],
                       weights = mydata_if$wt[mydata_if$in.ccs], robust = T)
  beta_srs[i] <- coef(cox_ccs_srs)
  # se_srs[i] <- summary(cox_ccs_srs)$coefficients[, "se(coef)"]
  # robust_se_srs[i] <- summary(cox_ccs_srs)$coefficients[, "robust se"]
  
  
  ####################################################################
  ######################## balanced sampling #########################
  ####################################################################
  s = samplecube(cbind(mydata_control$pik,mydata_control$inffun),
                 mydata_control$pik, order = 1, comment = F) # balanced sampling in the control set
  mydata_control$in.bs.subcohort = s # add sample indicator to the control set
  id.bs.subcohort = mydata_control[mydata_control$in.bs.subcohort == 1,]$id # get the id of the samples
  id.cases = mydata_if[mydata_if$delta==1,]$id # get the id of all cases
  id.ccs = c(id.cases,id.bs.subcohort) # combine the id of samples and cases
  bs.ccs <- mydata_if$id %in% id.ccs # case-control sample indicator
  mydata_if$in.bs.ccs <- bs.ccs # case-control sample indicator
  cox_ccs_bs <- coxph(Surv(time_obs, delta) ~ z1,
                      data = mydata_if[mydata_if$in.bs.ccs,],
                      weights = mydata_if$wt[mydata_if$in.bs.ccs], robust = T)
  beta_bs[i] <- coef(cox_ccs_bs)
  # se_bs[i] <- summary(cox_ccs_bs)$coefficients[, "se(coef)"]
  # robust_se_bs[i] <- summary(cox_ccs_bs)$coefficients[, "robust se"]
  # D1[i] <- cox_ccs_bs$var * n/(N-sum(mydata_if$delta))
  D2[i] <- cox_ccs_bs$naive.var  # estimated phase 1 variabilty
  # D3[i] <- coxph(Surv(time_obs, delta) ~ z1,
  #                data = mydata_if[mydata_if$in.bs.ccs,],
  #                weights = mydata_if$wt[mydata_if$in.bs.ccs])$var # should not be robust
  # D4[i] <- cox_ccs_bs$var # robust var
  
  
  ####################################################################
  ############################ calibration ###########################
  ####################################################################
  ccs_design <- twophase(id = list(~ 1, ~ 1),
                         subset = ~ in.ccs,
                         strata = list(NULL, ~ delta),
                         data = mydata_if)
  ccs_cox <- svycoxph(Surv(time_obs, delta) ~ z1, design = ccs_design)
  beta_srs_design[i] <- coef(ccs_cox)
  se_srs_design[i] <- summary(ccs_cox)$coefficients[, "se(coef)"]
  cch_cal <- calibrate(ccs_design, phase = 2, calfun = "raking", ~ inffun)
  cox_cal_cch <- svycoxph(Surv(time_obs, delta) ~ z1, design = cch_cal)
  beta_cali[i] <- coef(cox_cal_cch)
  # se_cali[i] <- summary(cox_cal_cch)$coefficients[, "se(coef)"]
  # robust_se_cali[i] <- summary(cox_cal_cch)$coefficients[, "robust se"]
  
  
  ####################################################################
  ############ re-calibration after balanced sampling ################
  ####################################################################
  bs_design <- twophase(id = list(~ 1, ~ 1),
                        subset = ~ in.bs.ccs,
                        strata = list(NULL, ~ delta),
                        data = mydata_if)
  bs_cox_design <- svycoxph(Surv(time_obs, delta) ~ z1, design = bs_design)
  beta_bs_design[i] <- coef(bs_cox_design)
  se_bs_design[i] <- summary(bs_cox_design)$coefficients[, "se(coef)"]
  bs_cal <- calibrate(bs_design, phase = 2, calfun = "raking", ~ inffun)
  cox_cal_bs <- svycoxph(Surv(time_obs, delta) ~ z1, design = bs_cal)
  beta_bs_cali[i] <- coef(cox_cal_bs)
  # se_bs_cali[i] <- summary(cox_cal_bs)$coefficients[, "se(coef)"]
  # robust_se_bs_cali[i] <- summary(cox_cal_bs)$coefficients[, "robust se"]
  
  ####################################################################
  ############ variance estimation for balanced sampling #############
  ####################################################################
  var_u_b <- function(Ys, Xs, pik, p, n) {
    c <- (1 - pik) * n / (n - p)
    b <- solve(t(Xs) %*% diag(c / pik^2) %*% Xs, t(Xs) %*% diag(c / pik^2) %*% Ys)
    vub <- c(t(Ys) %*% diag(c / pik^2) %*% Ys - t(Ys) %*% diag(c / pik^2) %*% Xs %*% b -
               t(b) %*% t(Xs) %*% diag(c / pik^2) %*% Ys + t(b) %*% t(Xs) %*% diag(c / pik^2) %*% Xs %*% b)
    return(vub)
  }
  
  
  # variance for the balanced sample only
  mydata_if$in.bs.subcohort_v=mydata_if$id %in% id.bs.subcohort
  mydata_ccs_bs <- mydata_if[mydata_if$in.bs.ccs,]
  mydata_ccs_bs$ccs_id <- 1 : nrow(mydata_ccs_bs)
  ccs_sub_id <- mydata_ccs_bs[mydata_ccs_bs$in.bs.subcohort_v == T,]$ccs_id
  Ys_bs_score = resid(cox_ccs_bs, "score")[ccs_sub_id]
  Ys_bs_dfbeta2 = resid(cox_ccs_bs, "dfbeta")[ccs_sub_id] * (n/(N-sum(mydata_if$delta))) 
  Xs_bs = cbind(mydata_ccs_bs$pik[ccs_sub_id],mydata_ccs_bs$inffun[ccs_sub_id])
  pik_bs = mydata_ccs_bs$pik[ccs_sub_id]
  w = (n/n_pop) / pik_bs[1]
 
  vub_bs_dfbeta2[i] = var_u_b(Ys=Ys_bs_dfbeta2, Xs=Xs_bs,
                              pik=pik_bs, p=2, n=sum(s))
  est_se_dfbeta2[i] <- sqrt(vub_bs_dfbeta2[i]) # estimated phase 2 variability
  est_se_D2_df[i] <- sqrt(D2[i] + vub_bs_dfbeta2[i]) # estimated phase 1 + phase 2 variability

  # for the coverage rate of 95% confidence interval
  # ci_up_score[i] <- coef(cox_ccs_bs) + qnorm(0.975) * est_se_D2_df[i]
  # ci_low_score[i] <- coef(cox_ccs_bs) - qnorm(0.975) * est_se_D2_df[i]
  # ci_ind_score[i] <- ifelse((log(2) >= ci_low_score[i]) & (log(2) <= ci_up_score[i]), 1, 0)
}

rr=data.frame(v=c("mean of beta srs","sd of beta srs",
                  "mean of beta bs","sd of beta bs",
                  "mean of beta full","sd of beta full",
                  "mean of beta cali","sd of beta cali",
                  "mean of beta bs cali","sd of beta bs cali",
                  "mean of phase 2 se","sd of phase 2 se",
                  "mean of two phase est se","sd of two phase est se",
                  "mean of phase 1 se","sd of phase 1 se"),
              n=c(mean(beta_srs),sd(beta_srs),
                  mean(beta_bs),sd(beta_bs),
                  mean(full_beta),sd(full_beta),
                  mean(beta_cali),sd(beta_cali),
                  mean(beta_bs_cali),sd(beta_bs_cali),
                  mean(est_se_dfbeta2),sd(est_se_dfbeta2),
                  mean(est_se_D2_df),sd(est_se_D2_df),
                  mean(sqrt(D2)),sd(sqrt(D2))))
(rr=data.frame(n=rr$v,n=round(rr$n,4)))


# mean of sqrt D2: phase 1
# mean of est se (dfbeta2): phase 2
# mean of two phase est se (D2 + vub2): est se
