###########################################################################################################
######### Case Cohort Sampling | high censoring | random cohort | continuous covariates | rho 0.8 #########
###########################################################################################################

rm(list = ls())
setwd("/Users/kc/Library/Mobile Documents/com~apple~CloudDocs/Downloads/thesis_bs")

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
full_beta <- c()
full_se <- c()
full_robust_se <- c()
phase1_var <- c()
phase2_var <- c()
est_se_two <- c()
est_se_two2 <- c()
beta_srs_design <- c()
se_srs_design <- c()
beta_cali <- c()
se_cali <- c()
robust_se_cali <- c()
beta_bs_design <- c()
se_bs_design <- c()
beta_bs_cali <- c()
se_bs_cali <- c()
robust_se_bs_cali <- c()
est_se_D <- c()
est_se_D_pi <- c()
est_se_D_pik <- c()
est_se_D2 <- c()
vub_bs_dfbeta2_pik <- c()
ci_up_se <- c()
ci_low_se <- c()
ci_up_robust <- c()
ci_low_robust <- c()
ci_ind_se <- c()
ci_ind_robust <- c()
ci_ind_se_log2 <- c()
ci_ind_robust_log2 <- c()
vub_bs_dfbeta2_w <- c()
est_se_D_score <- c()
info1 <- c()
D3 <- c()
est_se_D3_pik <- c()
est_se_D3 <- c()
D4 <-c()
ci_up_score <-c()
ci_low_score <-c()
ci_ind_score <-c()
est_se_D_score_w <- c()
est_se_D2_df_w <- c()
est_se_D2_df <- c()
for(i in 1 : sim) {
  cat("Simulation number ", i, "\n")
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
  mydata <- data.frame(id, time, censor, time_obs, delta, z_cts)
  full_cox <- coxph(Surv(time_obs, delta) ~ z1, data = mydata, robust = T)
  full_beta[i] <- coef(full_cox)
  full_se[i] <- summary(full_cox)$coefficients[, "se(coef)"]
  full_robust_se[i] <- summary(full_cox)$coefficients[, "robust se"]
  
  
  n <- 100; # subcohort size
  N <- nrow(mydata) # full cohort size
  ifmodel <- coxph(Surv(time_obs, delta) ~ z2, data = mydata)
  inffun <- resid(ifmodel, "dfbeta")
  mydata_if <- cbind(mydata, inffun)
  mydata_if$pik = ifelse(mydata_if$delta == 1, 1, n/(N - sum(mydata_if$delta))) # inclusion probability...!
  mydata_if$wt = 1/mydata_if$pik
  mydata_control <- mydata_if[mydata_if$delta == 0, ]
  ### srs ###
  casectrl <- with(mydata_if,c(which(delta==1),sample(which(delta==0),n)))
  mydata_if$in.ccs <- (1:nrow(mydata_if)) %in% casectrl
  cox_ccs_srs <- coxph(Surv(time_obs, delta) ~ z1, data = mydata_if[mydata_if$in.ccs,],
                       weights = mydata_if$wt[mydata_if$in.ccs], robust = T)
  beta_srs[i] <- coef(cox_ccs_srs)
  se_srs[i] <- summary(cox_ccs_srs)$coefficients[, "se(coef)"]
  robust_se_srs[i] <- summary(cox_ccs_srs)$coefficients[, "robust se"]
  ### bs ###
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
  se_bs[i] <- summary(cox_ccs_bs)$coefficients[, "se(coef)"]
  robust_se_bs[i] <- summary(cox_ccs_bs)$coefficients[, "robust se"]
  D1[i] <- cox_ccs_bs$var * n/(N-sum(mydata_if$delta))
  # sqrt(cox_ccs_bs$var)*sqrt(sum(mydata_if$in.bs.ccs)/N)
  D2[i] <- cox_ccs_bs$naive.var # ordinary estimate     # cox_ccs_bs$var
  D3[i] <- coxph(Surv(time_obs, delta) ~ z1,
                 data = mydata_if[mydata_if$in.bs.ccs,],
                 weights = mydata_if$wt[mydata_if$in.bs.ccs])$var # should not be robust
  D4[i] <- cox_ccs_bs$var # robust var
  compute_score_vector <- function(mydata, beta, w) {
    # mydata: dataset
    # beta: beta
    # w: sampling weights
    score <- 0 # initial score
    info <- 0
    # risk set
    for (i in 1 : sum(mydata$delta)) {
      t <- sort(mydata$time_obs[mydata$delta == 1])[i] # sort the observed time
      idx <- which(mydata$time_obs == t) # indices of subjects with observed event at t
      n_tie <- length(idx) # number of subjects with tied event times
      risk_id <- which(mydata$time_obs >= t) # risk set id
      s0 <- sum(w[risk_id] * exp(beta * mydata$z1[risk_id])) # denominator of Z_bar
      s1 <- sum(w[risk_id] * mydata$z1[risk_id] * exp(beta * mydata$z1[risk_id])) # numerator of Z_bar
      s11 <- sum(w[risk_id] * (mydata$z1[risk_id])^2 * exp(beta * mydata$z1[risk_id]))
      # calculate the score vector/info
      for (j in 1 : n_tie) {
        score <- score + w[idx][j] * mydata$delta[idx][j] * (mydata$z1[idx][j] - (s1 / s0))
        # summation over events
        info <- info + w[idx][j] * mydata$delta[idx][j] * (s11 * s0 - s1^2)/s0^2
      }
    }
    return(c(score,1/info))
  }
  beta0 = coef(coxph(Surv(time_obs, delta) ~ z1, data = mydata_if))
  w = 1/mydata_if$pik[mydata_if$in.bs.ccs]
  info1[i] <- compute_score_vector(mydata_if[mydata_if$in.bs.ccs,], beta0, w)[2]
  # ci_up_se[i] <- coef(cox_ccs_bs) + qnorm(0.975) * summary(cox_ccs_bs)$coefficients[, "se(coef)"]
  # ci_low_se[i] <- coef(cox_ccs_bs) - qnorm(0.975) * summary(cox_ccs_bs)$coefficients[, "se(coef)"]
  # ci_up_robust[i] <- coef(cox_ccs_bs) + qnorm(0.975) * summary(cox_ccs_bs)$coefficients[, "robust se"]
  # ci_low_robust[i] <- coef(cox_ccs_bs) - qnorm(0.975) * summary(cox_ccs_bs)$coefficients[, "robust se"]
  # ci_up[i] <- confint(cox_ccs_bs)[2]
  # ci_low[i] <- confint(cox_ccs_bs)[1]
  # ci_ind_se[i] <- ifelse((coef(full_cox) >= ci_low_se[i]) & (coef(full_cox) <= ci_up_se[i]), 1, 0)
  # ci_ind_robust[i] <- ifelse((coef(full_cox) >= ci_low_robust[i]) & (coef(full_cox) <= ci_up_robust[i]), 1, 0)
  # ci_ind_se_log2[i] <- ifelse((log(2) >= ci_low_se[i]) & (log(2) <= ci_up_se[i]), 1, 0)
  # ci_ind_robust_log2[i] <- ifelse((log(2) >= ci_low_robust[i]) & (log(2) <= ci_up_robust[i]), 1, 0)
  bs_design <- twophase(id = list(~ 1, ~ 1),
                        subset = ~ bs.ccs,
                        weights = list(NULL, ~ wt),
                        data= mydata_if,
                        method = "approx")
  cox_bs_design <- svycoxph(Surv(time_obs, delta) ~ z1, design = bs_design)
  v1 <- vcov(cox_bs_design)
  phase1_var[i] <- attr(v1, "phase")$phase1
  phase2_var[i] <- attr(v1, "phase")$phase2
  ### cali ###
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
  se_cali[i] <- summary(cox_cal_cch)$coefficients[, "se(coef)"]
  robust_se_cali[i] <- summary(cox_cal_cch)$coefficients[, "robust se"]
  ### cali bs ###
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
  se_bs_cali[i] <- summary(cox_cal_bs)$coefficients[, "se(coef)"]
  robust_se_bs_cali[i] <- summary(cox_cal_bs)$coefficients[, "robust se"]
  ### var est ###
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
  # Ys_bs_dfbeta1 = resid(cox_ccs_bs, "dfbeta")[ccs_sub_id] * (n/N) # * (n/N)
  Ys_bs_dfbeta2 = resid(cox_ccs_bs, "dfbeta")[ccs_sub_id] * (n/(N-sum(mydata_if$delta))) # * (n/N)
  Xs_bs = cbind(mydata_ccs_bs$pik[ccs_sub_id],mydata_ccs_bs$inffun[ccs_sub_id])
  pik_bs = mydata_ccs_bs$pik[ccs_sub_id]
  w = (n/n_pop) / pik_bs[1]
  vub_bs_score[i] = var_u_b(Ys=Ys_bs_score, Xs=Xs_bs,
                            pik=pik_bs, p=2, n=sum(s))
  # vub_bs_dfbeta1[i] = var_u_b(Ys=Ys_bs_dfbeta1, Xs=Xs_bs,
  #                             pik=pik_bs, p=2, n=sum(s))
  vub_bs_dfbeta2[i] = var_u_b(Ys=Ys_bs_dfbeta2, Xs=Xs_bs,
                              pik=pik_bs, p=2, n=sum(s))
  vub_bs_dfbeta2_w[i] <- vub_bs_dfbeta2[i] *  w #(1/N^2) # *((nrow(mydata_if)-sum(mydata_if$delta))/nrow(mydata_if)) # (1/N^2)
  vub_bs_dfbeta2_pik[i] = var_u_b(Ys=Ys_bs_dfbeta2, Xs=Xs_bs,
                                  pik=pik_bs, p=2, n=sum(s)) * pik_bs[1]
  est_se_score[i] <- sqrt(D3[i] * D3[i] * vub_bs_score[i]) # D1
  # est_se_dfbeta1[i] <- sqrt(vub_bs_dfbeta1[i])
  est_se_dfbeta2[i] <- sqrt(vub_bs_dfbeta2[i])
  # * (n/(N-sum(mydata_if$delta)))
  # * n/N
  est_se_D_score[i] <- sqrt(D2[i] + est_se_score[i]^2) # est_se_score[i]^2*w
  est_se_D_score_w[i] <- sqrt(D2[i] + est_se_score[i]^2*w)
  est_se_two[i] <-  sqrt(phase1_var[i] + vub_bs_dfbeta2[i])
  #est_se_D_pik[i] <-  sqrt(D2[i] + vub_bs_dfbeta2_pik[i])
  est_se_D2_df_w[i] <- sqrt(D2[i] + vub_bs_dfbeta2_w[i])
  est_se_D2_df[i] <- sqrt(D2[i] + vub_bs_dfbeta2[i])
  #est_se_D3_pik[i] <-  sqrt(D3[i] + vub_bs_dfbeta2_pik[i])
  est_se_D3[i] <- sqrt(D3[i] + vub_bs_dfbeta2_w[i])
  
  ci_up_score[i] <- coef(cox_ccs_bs) + qnorm(0.975) * est_se_D2_df[i]
  ci_low_score[i] <- coef(cox_ccs_bs) - qnorm(0.975) * est_se_D2_df[i]
  ci_ind_score[i] <- ifelse((log(2) >= ci_low_score[i]) & (log(2) <= ci_up_score[i]), 1, 0)
}

rr=data.frame(v=c("mean of beta srs","sd of beta srs",
                  "mean of beta bs","sd of beta bs",
                  "mean of beta full","sd of beta full",
                  "mean of beta cali","sd of beta cali",
                  "mean of beta bs cali","sd of beta bs cali",
                  "mean of se srs","sd of se srs",
                  "mean of se bs","sd of se bs",
                  "mean of robust se srs","sd of robust se srs",
                  "mean of robust se bs","sd of robust se bs",
                  "mean of est se (score)","sd of est se (score)",
                  "mean of est se (dfbeta2)","sd of est se (dfbeta2)",
                  "mean of two phase est se (full + vub2)","sd of two phase est se (full + vub2)",
                  "mean of two phase est se (svy1 + vub2)","sd of two phase est se (svy1 + vub2)",
                  "mean of two phase est se (D + vub_score)","sd of two phase est se (D + vub_score)",
                  "mean of two phase est se (D + vub_score*w)","sd of two phase est se (D + vub_score*w)",
                  "mean of two phase est se (D + vub2*w)","sd of two phase est se (D + vub2*w)",
                  "mean of two phase est se (D2 + vub2)","sd of two phase est se (D2 + vub2)",
                  "mean of two phase est se (D3 + vub2*w)","sd of two phase est se (D3 + vub2*w)",
                  "mean of sqrt D","sd of sqrt D",
                  "mean of D","sd of D",
                  "mean of sqrt D3","sd of sqrt D3",
                  "mean of D3","sd of D3",
                  "mean of info", "sd of info",
                  "ratio of CI log2"),
              n=c(mean(beta_srs),sd(beta_srs),
                  mean(beta_bs),sd(beta_bs),
                  mean(full_beta),sd(full_beta),
                  mean(beta_cali),sd(beta_cali),
                  mean(beta_bs_cali),sd(beta_bs_cali),
                  mean(se_srs),sd(se_srs),
                  mean(se_bs),sd(se_bs),
                  mean(robust_se_srs),sd(robust_se_srs),
                  mean(robust_se_bs),sd(robust_se_bs),
                  mean(est_se_score),sd(est_se_score),
                  mean(est_se_dfbeta2),sd(est_se_dfbeta2),
                  mean(sqrt(var(full_beta) + vub_bs_dfbeta2*w)),sd(sqrt(var(full_beta) + vub_bs_dfbeta2*w)),
                  mean(est_se_two),sd(est_se_two),
                  mean(est_se_D_score),sd(est_se_D_score),
                  mean(est_se_D_score_w),sd(est_se_D_score_w),
                  mean(est_se_D2_df_w),sd(est_se_D2_df_w),
                  mean(est_se_D2_df),sd(est_se_D2_df),
                  mean(est_se_D3),sd(est_se_D3),
                  mean(sqrt(D2)),sd(sqrt(D2)),
                  mean(D2),sd(D2),
                  mean(sqrt(D3)),sd(sqrt(D3)),
                  mean(D3),sd(D3),
                  mean(info1),sd(info1),
                  sum(ci_ind_score)/length(ci_ind_score)))
(rr=data.frame(n=rr$v,n=round(rr$n,4)))


write_xlsx(rr, paste("c1s1_ccs_conti_N", as.character(n_pop),"_n",as.character(n),"cor",as.character(rho),".xlsx", sep = ""))


write_xlsx(data.frame(beta_srs), paste("random_ccs_conti_beta_srs_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(beta_bs), paste("random_ccs_conti_beta_bs_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(beta_cali), paste("random_ccs_conti_beta_cali_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(beta_bs_cali), paste("random_ccs_conti_beta_cali_bs_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(est_se_D2_df), paste("random_ccs_conti_est_se_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(sqrt(D2)), paste("random_ccs_conti_se1_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(est_se_dfbeta2), paste("random_ccs_conti_se2_N_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(full_beta), paste("random_ccs_conti_beta_full_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))

# mean of sqrt D: phase 1
# mean of est se (dfbeta2): phase 1
# mean of two phase est se (D2 + vub2): est se


######### Case Cohort Sampling | high censoring | random cohort | continuous covariates | rho 0.8 #########


rm(list = ls())
setwd("/Users/kc/Library/Mobile Documents/com~apple~CloudDocs/Downloads/thesis_bs")

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
full_beta <- c()
full_se <- c()
full_robust_se <- c()
phase1_var <- c()
phase2_var <- c()
est_se_two <- c()
est_se_two2 <- c()
beta_srs_design <- c()
se_srs_design <- c()
beta_cali <- c()
se_cali <- c()
robust_se_cali <- c()
beta_bs_design <- c()
se_bs_design <- c()
beta_bs_cali <- c()
se_bs_cali <- c()
robust_se_bs_cali <- c()
est_se_D <- c()
est_se_D_pi <- c()
est_se_D_pik <- c()
est_se_D2 <- c()
vub_bs_dfbeta2_pik <- c()
ci_up_se <- c()
ci_low_se <- c()
ci_up_robust <- c()
ci_low_robust <- c()
ci_ind_se <- c()
ci_ind_robust <- c()
ci_ind_se_log2 <- c()
ci_ind_robust_log2 <- c()
vub_bs_dfbeta2_w <- c()
est_se_D_score <- c()
info1 <- c()
D3 <- c()
est_se_D3_pik <- c()
est_se_D3 <- c()
D4 <-c()
ci_up_score <-c()
ci_low_score <-c()
ci_ind_score <-c()
est_se_D_score_w <- c()
est_se_D2_df_w <- c()
est_se_D2_df <- c()
for(i in 1 : sim) {
  cat("Simulation number ", i, "\n")
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
  mydata <- data.frame(id, time, censor, time_obs, delta, z_cts)
  full_cox <- coxph(Surv(time_obs, delta) ~ z1, data = mydata, robust = T)
  full_beta[i] <- coef(full_cox)
  full_se[i] <- summary(full_cox)$coefficients[, "se(coef)"]
  full_robust_se[i] <- summary(full_cox)$coefficients[, "robust se"]
  
  
  n <- 200; # subcohort size
  N <- nrow(mydata) # full cohort size
  ifmodel <- coxph(Surv(time_obs, delta) ~ z2, data = mydata)
  inffun <- resid(ifmodel, "dfbeta")
  mydata_if <- cbind(mydata, inffun)
  mydata_if$pik = ifelse(mydata_if$delta == 1, 1, n/(N - sum(mydata_if$delta))) # inclusion probability...!
  mydata_if$wt = 1/mydata_if$pik
  mydata_control <- mydata_if[mydata_if$delta == 0, ]
  ### srs ###
  casectrl <- with(mydata_if,c(which(delta==1),sample(which(delta==0),n)))
  mydata_if$in.ccs <- (1:nrow(mydata_if)) %in% casectrl
  cox_ccs_srs <- coxph(Surv(time_obs, delta) ~ z1, data = mydata_if[mydata_if$in.ccs,],
                       weights = mydata_if$wt[mydata_if$in.ccs], robust = T)
  beta_srs[i] <- coef(cox_ccs_srs)
  se_srs[i] <- summary(cox_ccs_srs)$coefficients[, "se(coef)"]
  robust_se_srs[i] <- summary(cox_ccs_srs)$coefficients[, "robust se"]
  ### bs ###
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
  se_bs[i] <- summary(cox_ccs_bs)$coefficients[, "se(coef)"]
  robust_se_bs[i] <- summary(cox_ccs_bs)$coefficients[, "robust se"]
  D1[i] <- cox_ccs_bs$var * n/(N-sum(mydata_if$delta))
  # sqrt(cox_ccs_bs$var)*sqrt(sum(mydata_if$in.bs.ccs)/N)
  D2[i] <- cox_ccs_bs$naive.var # ordinary estimate     # cox_ccs_bs$var
  D3[i] <- coxph(Surv(time_obs, delta) ~ z1,
                 data = mydata_if[mydata_if$in.bs.ccs,],
                 weights = mydata_if$wt[mydata_if$in.bs.ccs])$var # should not be robust
  D4[i] <- cox_ccs_bs$var # robust var
  compute_score_vector <- function(mydata, beta, w) {
    # mydata: dataset
    # beta: beta
    # w: sampling weights
    score <- 0 # initial score
    info <- 0
    # risk set
    for (i in 1 : sum(mydata$delta)) {
      t <- sort(mydata$time_obs[mydata$delta == 1])[i] # sort the observed time
      idx <- which(mydata$time_obs == t) # indices of subjects with observed event at t
      n_tie <- length(idx) # number of subjects with tied event times
      risk_id <- which(mydata$time_obs >= t) # risk set id
      s0 <- sum(w[risk_id] * exp(beta * mydata$z1[risk_id])) # denominator of Z_bar
      s1 <- sum(w[risk_id] * mydata$z1[risk_id] * exp(beta * mydata$z1[risk_id])) # numerator of Z_bar
      s11 <- sum(w[risk_id] * (mydata$z1[risk_id])^2 * exp(beta * mydata$z1[risk_id]))
      # calculate the score vector/info
      for (j in 1 : n_tie) {
        score <- score + w[idx][j] * mydata$delta[idx][j] * (mydata$z1[idx][j] - (s1 / s0))
        # summation over events
        info <- info + w[idx][j] * mydata$delta[idx][j] * (s11 * s0 - s1^2)/s0^2
      }
    }
    return(c(score,1/info))
  }
  beta0 = coef(coxph(Surv(time_obs, delta) ~ z1, data = mydata_if))
  w = 1/mydata_if$pik[mydata_if$in.bs.ccs]
  info1[i] <- compute_score_vector(mydata_if[mydata_if$in.bs.ccs,], beta0, w)[2]
  # ci_up_se[i] <- coef(cox_ccs_bs) + qnorm(0.975) * summary(cox_ccs_bs)$coefficients[, "se(coef)"]
  # ci_low_se[i] <- coef(cox_ccs_bs) - qnorm(0.975) * summary(cox_ccs_bs)$coefficients[, "se(coef)"]
  # ci_up_robust[i] <- coef(cox_ccs_bs) + qnorm(0.975) * summary(cox_ccs_bs)$coefficients[, "robust se"]
  # ci_low_robust[i] <- coef(cox_ccs_bs) - qnorm(0.975) * summary(cox_ccs_bs)$coefficients[, "robust se"]
  # ci_up[i] <- confint(cox_ccs_bs)[2]
  # ci_low[i] <- confint(cox_ccs_bs)[1]
  # ci_ind_se[i] <- ifelse((coef(full_cox) >= ci_low_se[i]) & (coef(full_cox) <= ci_up_se[i]), 1, 0)
  # ci_ind_robust[i] <- ifelse((coef(full_cox) >= ci_low_robust[i]) & (coef(full_cox) <= ci_up_robust[i]), 1, 0)
  # ci_ind_se_log2[i] <- ifelse((log(2) >= ci_low_se[i]) & (log(2) <= ci_up_se[i]), 1, 0)
  # ci_ind_robust_log2[i] <- ifelse((log(2) >= ci_low_robust[i]) & (log(2) <= ci_up_robust[i]), 1, 0)
  bs_design <- twophase(id = list(~ 1, ~ 1),
                        subset = ~ bs.ccs,
                        weights = list(NULL, ~ wt),
                        data= mydata_if,
                        method = "approx")
  cox_bs_design <- svycoxph(Surv(time_obs, delta) ~ z1, design = bs_design)
  v1 <- vcov(cox_bs_design)
  phase1_var[i] <- attr(v1, "phase")$phase1
  phase2_var[i] <- attr(v1, "phase")$phase2
  ### cali ###
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
  se_cali[i] <- summary(cox_cal_cch)$coefficients[, "se(coef)"]
  robust_se_cali[i] <- summary(cox_cal_cch)$coefficients[, "robust se"]
  ### cali bs ###
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
  se_bs_cali[i] <- summary(cox_cal_bs)$coefficients[, "se(coef)"]
  robust_se_bs_cali[i] <- summary(cox_cal_bs)$coefficients[, "robust se"]
  ### var est ###
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
  # Ys_bs_dfbeta1 = resid(cox_ccs_bs, "dfbeta")[ccs_sub_id] * (n/N) # * (n/N)
  Ys_bs_dfbeta2 = resid(cox_ccs_bs, "dfbeta")[ccs_sub_id] * (n/(N-sum(mydata_if$delta))) # * (n/N)
  Xs_bs = cbind(mydata_ccs_bs$pik[ccs_sub_id],mydata_ccs_bs$inffun[ccs_sub_id])
  pik_bs = mydata_ccs_bs$pik[ccs_sub_id]
  w = (n/n_pop) / pik_bs[1]
  vub_bs_score[i] = var_u_b(Ys=Ys_bs_score, Xs=Xs_bs,
                            pik=pik_bs, p=2, n=sum(s))
  # vub_bs_dfbeta1[i] = var_u_b(Ys=Ys_bs_dfbeta1, Xs=Xs_bs,
  #                             pik=pik_bs, p=2, n=sum(s))
  vub_bs_dfbeta2[i] = var_u_b(Ys=Ys_bs_dfbeta2, Xs=Xs_bs,
                              pik=pik_bs, p=2, n=sum(s))
  vub_bs_dfbeta2_w[i] <- vub_bs_dfbeta2[i] *  w #(1/N^2) # *((nrow(mydata_if)-sum(mydata_if$delta))/nrow(mydata_if)) # (1/N^2)
  vub_bs_dfbeta2_pik[i] = var_u_b(Ys=Ys_bs_dfbeta2, Xs=Xs_bs,
                                  pik=pik_bs, p=2, n=sum(s)) * pik_bs[1]
  est_se_score[i] <- sqrt(D3[i] * D3[i] * vub_bs_score[i]) # D1
  # est_se_dfbeta1[i] <- sqrt(vub_bs_dfbeta1[i])
  est_se_dfbeta2[i] <- sqrt(vub_bs_dfbeta2[i])
  # * (n/(N-sum(mydata_if$delta)))
  # * n/N
  est_se_D_score[i] <- sqrt(D2[i] + est_se_score[i]^2) # est_se_score[i]^2*w
  est_se_D_score_w[i] <- sqrt(D2[i] + est_se_score[i]^2*w)
  est_se_two[i] <-  sqrt(phase1_var[i] + vub_bs_dfbeta2[i])
  #est_se_D_pik[i] <-  sqrt(D2[i] + vub_bs_dfbeta2_pik[i])
  est_se_D2_df_w[i] <- sqrt(D2[i] + vub_bs_dfbeta2_w[i])
  est_se_D2_df[i] <- sqrt(D2[i] + vub_bs_dfbeta2[i])
  #est_se_D3_pik[i] <-  sqrt(D3[i] + vub_bs_dfbeta2_pik[i])
  est_se_D3[i] <- sqrt(D3[i] + vub_bs_dfbeta2_w[i])
  
  ci_up_score[i] <- coef(cox_ccs_bs) + qnorm(0.975) * est_se_D2_df[i]
  ci_low_score[i] <- coef(cox_ccs_bs) - qnorm(0.975) * est_se_D2_df[i]
  ci_ind_score[i] <- ifelse((log(2) >= ci_low_score[i]) & (log(2) <= ci_up_score[i]), 1, 0)
}

rr=data.frame(v=c("mean of beta srs","sd of beta srs",
                  "mean of beta bs","sd of beta bs",
                  "mean of beta full","sd of beta full",
                  "mean of beta cali","sd of beta cali",
                  "mean of beta bs cali","sd of beta bs cali",
                  "mean of se srs","sd of se srs",
                  "mean of se bs","sd of se bs",
                  "mean of robust se srs","sd of robust se srs",
                  "mean of robust se bs","sd of robust se bs",
                  "mean of est se (score)","sd of est se (score)",
                  "mean of est se (dfbeta2)","sd of est se (dfbeta2)",
                  "mean of two phase est se (full + vub2)","sd of two phase est se (full + vub2)",
                  "mean of two phase est se (svy1 + vub2)","sd of two phase est se (svy1 + vub2)",
                  "mean of two phase est se (D + vub_score)","sd of two phase est se (D + vub_score)",
                  "mean of two phase est se (D + vub_score*w)","sd of two phase est se (D + vub_score*w)",
                  "mean of two phase est se (D + vub2*w)","sd of two phase est se (D + vub2*w)",
                  "mean of two phase est se (D2 + vub2)","sd of two phase est se (D2 + vub2)",
                  "mean of two phase est se (D3 + vub2*w)","sd of two phase est se (D3 + vub2*w)",
                  "mean of sqrt D","sd of sqrt D",
                  "mean of D","sd of D",
                  "mean of sqrt D3","sd of sqrt D3",
                  "mean of D3","sd of D3",
                  "mean of info", "sd of info",
                  "ratio of CI log2"),
              n=c(mean(beta_srs),sd(beta_srs),
                  mean(beta_bs),sd(beta_bs),
                  mean(full_beta),sd(full_beta),
                  mean(beta_cali),sd(beta_cali),
                  mean(beta_bs_cali),sd(beta_bs_cali),
                  mean(se_srs),sd(se_srs),
                  mean(se_bs),sd(se_bs),
                  mean(robust_se_srs),sd(robust_se_srs),
                  mean(robust_se_bs),sd(robust_se_bs),
                  mean(est_se_score),sd(est_se_score),
                  mean(est_se_dfbeta2),sd(est_se_dfbeta2),
                  mean(sqrt(var(full_beta) + vub_bs_dfbeta2*w)),sd(sqrt(var(full_beta) + vub_bs_dfbeta2*w)),
                  mean(est_se_two),sd(est_se_two),
                  mean(est_se_D_score),sd(est_se_D_score),
                  mean(est_se_D_score_w),sd(est_se_D_score_w),
                  mean(est_se_D2_df_w),sd(est_se_D2_df_w),
                  mean(est_se_D2_df),sd(est_se_D2_df),
                  mean(est_se_D3),sd(est_se_D3),
                  mean(sqrt(D2)),sd(sqrt(D2)),
                  mean(D2),sd(D2),
                  mean(sqrt(D3)),sd(sqrt(D3)),
                  mean(D3),sd(D3),
                  mean(info1),sd(info1),
                  sum(ci_ind_score)/length(ci_ind_score)))
(rr=data.frame(n=rr$v,n=round(rr$n,4)))


write_xlsx(rr, paste("c1s1_ccs_conti_N", as.character(n_pop),"_n",as.character(n),"cor",as.character(rho),".xlsx", sep = ""))


write_xlsx(data.frame(beta_srs), paste("random_ccs_conti_beta_srs_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(beta_bs), paste("random_ccs_conti_beta_bs_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(beta_cali), paste("random_ccs_conti_beta_cali_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(beta_bs_cali), paste("random_ccs_conti_beta_cali_bs_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(est_se_D2_df), paste("random_ccs_conti_est_se_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(sqrt(D2)), paste("random_ccs_conti_se1_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(est_se_dfbeta2), paste("random_ccs_conti_se2_N_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(full_beta), paste("random_ccs_conti_beta_full_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))

# mean of sqrt D: phase 1
# mean of est se (dfbeta2): phase 1
# mean of two phase est se (D2 + vub2): est se

######### Case Cohort Sampling | high censoring | random cohort | continuous covariates | rho 0.8 #########


rm(list = ls())
setwd("/Users/kc/Library/Mobile Documents/com~apple~CloudDocs/Downloads/thesis_bs")

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
full_beta <- c()
full_se <- c()
full_robust_se <- c()
phase1_var <- c()
phase2_var <- c()
est_se_two <- c()
est_se_two2 <- c()
beta_srs_design <- c()
se_srs_design <- c()
beta_cali <- c()
se_cali <- c()
robust_se_cali <- c()
beta_bs_design <- c()
se_bs_design <- c()
beta_bs_cali <- c()
se_bs_cali <- c()
robust_se_bs_cali <- c()
est_se_D <- c()
est_se_D_pi <- c()
est_se_D_pik <- c()
est_se_D2 <- c()
vub_bs_dfbeta2_pik <- c()
ci_up_se <- c()
ci_low_se <- c()
ci_up_robust <- c()
ci_low_robust <- c()
ci_ind_se <- c()
ci_ind_robust <- c()
ci_ind_se_log2 <- c()
ci_ind_robust_log2 <- c()
vub_bs_dfbeta2_w <- c()
est_se_D_score <- c()
info1 <- c()
D3 <- c()
est_se_D3_pik <- c()
est_se_D3 <- c()
D4 <-c()
ci_up_score <-c()
ci_low_score <-c()
ci_ind_score <-c()
est_se_D_score_w <- c()
est_se_D2_df_w <- c()
est_se_D2_df <- c()
for(i in 1 : sim) {
  cat("Simulation number ", i, "\n")
  n_pop <- 3000 # full cohort size
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
  mydata <- data.frame(id, time, censor, time_obs, delta, z_cts)
  full_cox <- coxph(Surv(time_obs, delta) ~ z1, data = mydata, robust = T)
  full_beta[i] <- coef(full_cox)
  full_se[i] <- summary(full_cox)$coefficients[, "se(coef)"]
  full_robust_se[i] <- summary(full_cox)$coefficients[, "robust se"]
  
  
  n <- 300; # subcohort size
  N <- nrow(mydata) # full cohort size
  ifmodel <- coxph(Surv(time_obs, delta) ~ z2, data = mydata)
  inffun <- resid(ifmodel, "dfbeta")
  mydata_if <- cbind(mydata, inffun)
  mydata_if$pik = ifelse(mydata_if$delta == 1, 1, n/(N - sum(mydata_if$delta))) # inclusion probability...!
  mydata_if$wt = 1/mydata_if$pik
  mydata_control <- mydata_if[mydata_if$delta == 0, ]
  ### srs ###
  casectrl <- with(mydata_if,c(which(delta==1),sample(which(delta==0),n)))
  mydata_if$in.ccs <- (1:nrow(mydata_if)) %in% casectrl
  cox_ccs_srs <- coxph(Surv(time_obs, delta) ~ z1, data = mydata_if[mydata_if$in.ccs,],
                       weights = mydata_if$wt[mydata_if$in.ccs], robust = T)
  beta_srs[i] <- coef(cox_ccs_srs)
  se_srs[i] <- summary(cox_ccs_srs)$coefficients[, "se(coef)"]
  robust_se_srs[i] <- summary(cox_ccs_srs)$coefficients[, "robust se"]
  ### bs ###
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
  se_bs[i] <- summary(cox_ccs_bs)$coefficients[, "se(coef)"]
  robust_se_bs[i] <- summary(cox_ccs_bs)$coefficients[, "robust se"]
  D1[i] <- cox_ccs_bs$var * n/(N-sum(mydata_if$delta))
  # sqrt(cox_ccs_bs$var)*sqrt(sum(mydata_if$in.bs.ccs)/N)
  D2[i] <- cox_ccs_bs$naive.var # ordinary estimate     # cox_ccs_bs$var
  D3[i] <- coxph(Surv(time_obs, delta) ~ z1,
                 data = mydata_if[mydata_if$in.bs.ccs,],
                 weights = mydata_if$wt[mydata_if$in.bs.ccs])$var # should not be robust
  D4[i] <- cox_ccs_bs$var # robust var
  compute_score_vector <- function(mydata, beta, w) {
    # mydata: dataset
    # beta: beta
    # w: sampling weights
    score <- 0 # initial score
    info <- 0
    # risk set
    for (i in 1 : sum(mydata$delta)) {
      t <- sort(mydata$time_obs[mydata$delta == 1])[i] # sort the observed time
      idx <- which(mydata$time_obs == t) # indices of subjects with observed event at t
      n_tie <- length(idx) # number of subjects with tied event times
      risk_id <- which(mydata$time_obs >= t) # risk set id
      s0 <- sum(w[risk_id] * exp(beta * mydata$z1[risk_id])) # denominator of Z_bar
      s1 <- sum(w[risk_id] * mydata$z1[risk_id] * exp(beta * mydata$z1[risk_id])) # numerator of Z_bar
      s11 <- sum(w[risk_id] * (mydata$z1[risk_id])^2 * exp(beta * mydata$z1[risk_id]))
      # calculate the score vector/info
      for (j in 1 : n_tie) {
        score <- score + w[idx][j] * mydata$delta[idx][j] * (mydata$z1[idx][j] - (s1 / s0))
        # summation over events
        info <- info + w[idx][j] * mydata$delta[idx][j] * (s11 * s0 - s1^2)/s0^2
      }
    }
    return(c(score,1/info))
  }
  beta0 = coef(coxph(Surv(time_obs, delta) ~ z1, data = mydata_if))
  w = 1/mydata_if$pik[mydata_if$in.bs.ccs]
  info1[i] <- compute_score_vector(mydata_if[mydata_if$in.bs.ccs,], beta0, w)[2]
  # ci_up_se[i] <- coef(cox_ccs_bs) + qnorm(0.975) * summary(cox_ccs_bs)$coefficients[, "se(coef)"]
  # ci_low_se[i] <- coef(cox_ccs_bs) - qnorm(0.975) * summary(cox_ccs_bs)$coefficients[, "se(coef)"]
  # ci_up_robust[i] <- coef(cox_ccs_bs) + qnorm(0.975) * summary(cox_ccs_bs)$coefficients[, "robust se"]
  # ci_low_robust[i] <- coef(cox_ccs_bs) - qnorm(0.975) * summary(cox_ccs_bs)$coefficients[, "robust se"]
  # ci_up[i] <- confint(cox_ccs_bs)[2]
  # ci_low[i] <- confint(cox_ccs_bs)[1]
  # ci_ind_se[i] <- ifelse((coef(full_cox) >= ci_low_se[i]) & (coef(full_cox) <= ci_up_se[i]), 1, 0)
  # ci_ind_robust[i] <- ifelse((coef(full_cox) >= ci_low_robust[i]) & (coef(full_cox) <= ci_up_robust[i]), 1, 0)
  # ci_ind_se_log2[i] <- ifelse((log(2) >= ci_low_se[i]) & (log(2) <= ci_up_se[i]), 1, 0)
  # ci_ind_robust_log2[i] <- ifelse((log(2) >= ci_low_robust[i]) & (log(2) <= ci_up_robust[i]), 1, 0)
  bs_design <- twophase(id = list(~ 1, ~ 1),
                        subset = ~ bs.ccs,
                        weights = list(NULL, ~ wt),
                        data= mydata_if,
                        method = "approx")
  cox_bs_design <- svycoxph(Surv(time_obs, delta) ~ z1, design = bs_design)
  v1 <- vcov(cox_bs_design)
  phase1_var[i] <- attr(v1, "phase")$phase1
  phase2_var[i] <- attr(v1, "phase")$phase2
  ### cali ###
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
  se_cali[i] <- summary(cox_cal_cch)$coefficients[, "se(coef)"]
  robust_se_cali[i] <- summary(cox_cal_cch)$coefficients[, "robust se"]
  ### cali bs ###
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
  se_bs_cali[i] <- summary(cox_cal_bs)$coefficients[, "se(coef)"]
  robust_se_bs_cali[i] <- summary(cox_cal_bs)$coefficients[, "robust se"]
  ### var est ###
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
  # Ys_bs_dfbeta1 = resid(cox_ccs_bs, "dfbeta")[ccs_sub_id] * (n/N) # * (n/N)
  Ys_bs_dfbeta2 = resid(cox_ccs_bs, "dfbeta")[ccs_sub_id] * (n/(N-sum(mydata_if$delta))) # * (n/N)
  Xs_bs = cbind(mydata_ccs_bs$pik[ccs_sub_id],mydata_ccs_bs$inffun[ccs_sub_id])
  pik_bs = mydata_ccs_bs$pik[ccs_sub_id]
  w = (n/n_pop) / pik_bs[1]
  vub_bs_score[i] = var_u_b(Ys=Ys_bs_score, Xs=Xs_bs,
                            pik=pik_bs, p=2, n=sum(s))
  # vub_bs_dfbeta1[i] = var_u_b(Ys=Ys_bs_dfbeta1, Xs=Xs_bs,
  #                             pik=pik_bs, p=2, n=sum(s))
  vub_bs_dfbeta2[i] = var_u_b(Ys=Ys_bs_dfbeta2, Xs=Xs_bs,
                              pik=pik_bs, p=2, n=sum(s))
  vub_bs_dfbeta2_w[i] <- vub_bs_dfbeta2[i] *  w #(1/N^2) # *((nrow(mydata_if)-sum(mydata_if$delta))/nrow(mydata_if)) # (1/N^2)
  vub_bs_dfbeta2_pik[i] = var_u_b(Ys=Ys_bs_dfbeta2, Xs=Xs_bs,
                                  pik=pik_bs, p=2, n=sum(s)) * pik_bs[1]
  est_se_score[i] <- sqrt(D3[i] * D3[i] * vub_bs_score[i]) # D1
  # est_se_dfbeta1[i] <- sqrt(vub_bs_dfbeta1[i])
  est_se_dfbeta2[i] <- sqrt(vub_bs_dfbeta2[i])
  # * (n/(N-sum(mydata_if$delta)))
  # * n/N
  est_se_D_score[i] <- sqrt(D2[i] + est_se_score[i]^2) # est_se_score[i]^2*w
  est_se_D_score_w[i] <- sqrt(D2[i] + est_se_score[i]^2*w)
  est_se_two[i] <-  sqrt(phase1_var[i] + vub_bs_dfbeta2[i])
  #est_se_D_pik[i] <-  sqrt(D2[i] + vub_bs_dfbeta2_pik[i])
  est_se_D2_df_w[i] <- sqrt(D2[i] + vub_bs_dfbeta2_w[i])
  est_se_D2_df[i] <- sqrt(D2[i] + vub_bs_dfbeta2[i])
  #est_se_D3_pik[i] <-  sqrt(D3[i] + vub_bs_dfbeta2_pik[i])
  est_se_D3[i] <- sqrt(D3[i] + vub_bs_dfbeta2_w[i])
  
  ci_up_score[i] <- coef(cox_ccs_bs) + qnorm(0.975) * est_se_D2_df[i]
  ci_low_score[i] <- coef(cox_ccs_bs) - qnorm(0.975) * est_se_D2_df[i]
  ci_ind_score[i] <- ifelse((log(2) >= ci_low_score[i]) & (log(2) <= ci_up_score[i]), 1, 0)
}

rr=data.frame(v=c("mean of beta srs","sd of beta srs",
                  "mean of beta bs","sd of beta bs",
                  "mean of beta full","sd of beta full",
                  "mean of beta cali","sd of beta cali",
                  "mean of beta bs cali","sd of beta bs cali",
                  "mean of se srs","sd of se srs",
                  "mean of se bs","sd of se bs",
                  "mean of robust se srs","sd of robust se srs",
                  "mean of robust se bs","sd of robust se bs",
                  "mean of est se (score)","sd of est se (score)",
                  "mean of est se (dfbeta2)","sd of est se (dfbeta2)",
                  "mean of two phase est se (full + vub2)","sd of two phase est se (full + vub2)",
                  "mean of two phase est se (svy1 + vub2)","sd of two phase est se (svy1 + vub2)",
                  "mean of two phase est se (D + vub_score)","sd of two phase est se (D + vub_score)",
                  "mean of two phase est se (D + vub_score*w)","sd of two phase est se (D + vub_score*w)",
                  "mean of two phase est se (D + vub2*w)","sd of two phase est se (D + vub2*w)",
                  "mean of two phase est se (D2 + vub2)","sd of two phase est se (D2 + vub2)",
                  "mean of two phase est se (D3 + vub2*w)","sd of two phase est se (D3 + vub2*w)",
                  "mean of sqrt D","sd of sqrt D",
                  "mean of D","sd of D",
                  "mean of sqrt D3","sd of sqrt D3",
                  "mean of D3","sd of D3",
                  "mean of info", "sd of info",
                  "ratio of CI log2"),
              n=c(mean(beta_srs),sd(beta_srs),
                  mean(beta_bs),sd(beta_bs),
                  mean(full_beta),sd(full_beta),
                  mean(beta_cali),sd(beta_cali),
                  mean(beta_bs_cali),sd(beta_bs_cali),
                  mean(se_srs),sd(se_srs),
                  mean(se_bs),sd(se_bs),
                  mean(robust_se_srs),sd(robust_se_srs),
                  mean(robust_se_bs),sd(robust_se_bs),
                  mean(est_se_score),sd(est_se_score),
                  mean(est_se_dfbeta2),sd(est_se_dfbeta2),
                  mean(sqrt(var(full_beta) + vub_bs_dfbeta2*w)),sd(sqrt(var(full_beta) + vub_bs_dfbeta2*w)),
                  mean(est_se_two),sd(est_se_two),
                  mean(est_se_D_score),sd(est_se_D_score),
                  mean(est_se_D_score_w),sd(est_se_D_score_w),
                  mean(est_se_D2_df_w),sd(est_se_D2_df_w),
                  mean(est_se_D2_df),sd(est_se_D2_df),
                  mean(est_se_D3),sd(est_se_D3),
                  mean(sqrt(D2)),sd(sqrt(D2)),
                  mean(D2),sd(D2),
                  mean(sqrt(D3)),sd(sqrt(D3)),
                  mean(D3),sd(D3),
                  mean(info1),sd(info1),
                  sum(ci_ind_score)/length(ci_ind_score)))
(rr=data.frame(n=rr$v,n=round(rr$n,4)))


write_xlsx(rr, paste("c1s1_ccs_conti_N", as.character(n_pop),"_n",as.character(n),"cor",as.character(rho),".xlsx", sep = ""))


write_xlsx(data.frame(beta_srs), paste("random_ccs_conti_beta_srs_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(beta_bs), paste("random_ccs_conti_beta_bs_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(beta_cali), paste("random_ccs_conti_beta_cali_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(beta_bs_cali), paste("random_ccs_conti_beta_cali_bs_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(est_se_D2_df), paste("random_ccs_conti_est_se_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(sqrt(D2)), paste("random_ccs_conti_se1_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(est_se_dfbeta2), paste("random_ccs_conti_se2_N_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(full_beta), paste("random_ccs_conti_beta_full_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))

# mean of sqrt D: phase 1
# mean of est se (dfbeta2): phase 1
# mean of two phase est se (D2 + vub2): est se


######### Case Cohort Sampling | high censoring | random cohort | continuous covariates | rho 0.8 #########


rm(list = ls())
setwd("/Users/kc/Library/Mobile Documents/com~apple~CloudDocs/Downloads/thesis_bs")

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
full_beta <- c()
full_se <- c()
full_robust_se <- c()
phase1_var <- c()
phase2_var <- c()
est_se_two <- c()
est_se_two2 <- c()
beta_srs_design <- c()
se_srs_design <- c()
beta_cali <- c()
se_cali <- c()
robust_se_cali <- c()
beta_bs_design <- c()
se_bs_design <- c()
beta_bs_cali <- c()
se_bs_cali <- c()
robust_se_bs_cali <- c()
est_se_D <- c()
est_se_D_pi <- c()
est_se_D_pik <- c()
est_se_D2 <- c()
vub_bs_dfbeta2_pik <- c()
ci_up_se <- c()
ci_low_se <- c()
ci_up_robust <- c()
ci_low_robust <- c()
ci_ind_se <- c()
ci_ind_robust <- c()
ci_ind_se_log2 <- c()
ci_ind_robust_log2 <- c()
vub_bs_dfbeta2_w <- c()
est_se_D_score <- c()
info1 <- c()
D3 <- c()
est_se_D3_pik <- c()
est_se_D3 <- c()
D4 <-c()
ci_up_score <-c()
ci_low_score <-c()
ci_ind_score <-c()
est_se_D_score_w <- c()
est_se_D2_df_w <- c()
est_se_D2_df <- c()
for(i in 1 : sim) {
  cat("Simulation number ", i, "\n")
  n_pop <- 3000 # full cohort size
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
  mydata <- data.frame(id, time, censor, time_obs, delta, z_cts)
  full_cox <- coxph(Surv(time_obs, delta) ~ z1, data = mydata, robust = T)
  full_beta[i] <- coef(full_cox)
  full_se[i] <- summary(full_cox)$coefficients[, "se(coef)"]
  full_robust_se[i] <- summary(full_cox)$coefficients[, "robust se"]
  
  
  n <- 600; # subcohort size
  N <- nrow(mydata) # full cohort size
  ifmodel <- coxph(Surv(time_obs, delta) ~ z2, data = mydata)
  inffun <- resid(ifmodel, "dfbeta")
  mydata_if <- cbind(mydata, inffun)
  mydata_if$pik = ifelse(mydata_if$delta == 1, 1, n/(N - sum(mydata_if$delta))) # inclusion probability...!
  mydata_if$wt = 1/mydata_if$pik
  mydata_control <- mydata_if[mydata_if$delta == 0, ]
  ### srs ###
  casectrl <- with(mydata_if,c(which(delta==1),sample(which(delta==0),n)))
  mydata_if$in.ccs <- (1:nrow(mydata_if)) %in% casectrl
  cox_ccs_srs <- coxph(Surv(time_obs, delta) ~ z1, data = mydata_if[mydata_if$in.ccs,],
                       weights = mydata_if$wt[mydata_if$in.ccs], robust = T)
  beta_srs[i] <- coef(cox_ccs_srs)
  se_srs[i] <- summary(cox_ccs_srs)$coefficients[, "se(coef)"]
  robust_se_srs[i] <- summary(cox_ccs_srs)$coefficients[, "robust se"]
  ### bs ###
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
  se_bs[i] <- summary(cox_ccs_bs)$coefficients[, "se(coef)"]
  robust_se_bs[i] <- summary(cox_ccs_bs)$coefficients[, "robust se"]
  D1[i] <- cox_ccs_bs$var * n/(N-sum(mydata_if$delta))
  # sqrt(cox_ccs_bs$var)*sqrt(sum(mydata_if$in.bs.ccs)/N)
  D2[i] <- cox_ccs_bs$naive.var # ordinary estimate     # cox_ccs_bs$var
  D3[i] <- coxph(Surv(time_obs, delta) ~ z1,
                 data = mydata_if[mydata_if$in.bs.ccs,],
                 weights = mydata_if$wt[mydata_if$in.bs.ccs])$var # should not be robust
  D4[i] <- cox_ccs_bs$var # robust var
  compute_score_vector <- function(mydata, beta, w) {
    # mydata: dataset
    # beta: beta
    # w: sampling weights
    score <- 0 # initial score
    info <- 0
    # risk set
    for (i in 1 : sum(mydata$delta)) {
      t <- sort(mydata$time_obs[mydata$delta == 1])[i] # sort the observed time
      idx <- which(mydata$time_obs == t) # indices of subjects with observed event at t
      n_tie <- length(idx) # number of subjects with tied event times
      risk_id <- which(mydata$time_obs >= t) # risk set id
      s0 <- sum(w[risk_id] * exp(beta * mydata$z1[risk_id])) # denominator of Z_bar
      s1 <- sum(w[risk_id] * mydata$z1[risk_id] * exp(beta * mydata$z1[risk_id])) # numerator of Z_bar
      s11 <- sum(w[risk_id] * (mydata$z1[risk_id])^2 * exp(beta * mydata$z1[risk_id]))
      # calculate the score vector/info
      for (j in 1 : n_tie) {
        score <- score + w[idx][j] * mydata$delta[idx][j] * (mydata$z1[idx][j] - (s1 / s0))
        # summation over events
        info <- info + w[idx][j] * mydata$delta[idx][j] * (s11 * s0 - s1^2)/s0^2
      }
    }
    return(c(score,1/info))
  }
  beta0 = coef(coxph(Surv(time_obs, delta) ~ z1, data = mydata_if))
  w = 1/mydata_if$pik[mydata_if$in.bs.ccs]
  info1[i] <- compute_score_vector(mydata_if[mydata_if$in.bs.ccs,], beta0, w)[2]
  # ci_up_se[i] <- coef(cox_ccs_bs) + qnorm(0.975) * summary(cox_ccs_bs)$coefficients[, "se(coef)"]
  # ci_low_se[i] <- coef(cox_ccs_bs) - qnorm(0.975) * summary(cox_ccs_bs)$coefficients[, "se(coef)"]
  # ci_up_robust[i] <- coef(cox_ccs_bs) + qnorm(0.975) * summary(cox_ccs_bs)$coefficients[, "robust se"]
  # ci_low_robust[i] <- coef(cox_ccs_bs) - qnorm(0.975) * summary(cox_ccs_bs)$coefficients[, "robust se"]
  # ci_up[i] <- confint(cox_ccs_bs)[2]
  # ci_low[i] <- confint(cox_ccs_bs)[1]
  # ci_ind_se[i] <- ifelse((coef(full_cox) >= ci_low_se[i]) & (coef(full_cox) <= ci_up_se[i]), 1, 0)
  # ci_ind_robust[i] <- ifelse((coef(full_cox) >= ci_low_robust[i]) & (coef(full_cox) <= ci_up_robust[i]), 1, 0)
  # ci_ind_se_log2[i] <- ifelse((log(2) >= ci_low_se[i]) & (log(2) <= ci_up_se[i]), 1, 0)
  # ci_ind_robust_log2[i] <- ifelse((log(2) >= ci_low_robust[i]) & (log(2) <= ci_up_robust[i]), 1, 0)
  bs_design <- twophase(id = list(~ 1, ~ 1),
                        subset = ~ bs.ccs,
                        weights = list(NULL, ~ wt),
                        data= mydata_if,
                        method = "approx")
  cox_bs_design <- svycoxph(Surv(time_obs, delta) ~ z1, design = bs_design)
  v1 <- vcov(cox_bs_design)
  phase1_var[i] <- attr(v1, "phase")$phase1
  phase2_var[i] <- attr(v1, "phase")$phase2
  ### cali ###
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
  se_cali[i] <- summary(cox_cal_cch)$coefficients[, "se(coef)"]
  robust_se_cali[i] <- summary(cox_cal_cch)$coefficients[, "robust se"]
  ### cali bs ###
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
  se_bs_cali[i] <- summary(cox_cal_bs)$coefficients[, "se(coef)"]
  robust_se_bs_cali[i] <- summary(cox_cal_bs)$coefficients[, "robust se"]
  ### var est ###
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
  # Ys_bs_dfbeta1 = resid(cox_ccs_bs, "dfbeta")[ccs_sub_id] * (n/N) # * (n/N)
  Ys_bs_dfbeta2 = resid(cox_ccs_bs, "dfbeta")[ccs_sub_id] * (n/(N-sum(mydata_if$delta))) # * (n/N)
  Xs_bs = cbind(mydata_ccs_bs$pik[ccs_sub_id],mydata_ccs_bs$inffun[ccs_sub_id])
  pik_bs = mydata_ccs_bs$pik[ccs_sub_id]
  w = (n/n_pop) / pik_bs[1]
  vub_bs_score[i] = var_u_b(Ys=Ys_bs_score, Xs=Xs_bs,
                            pik=pik_bs, p=2, n=sum(s))
  # vub_bs_dfbeta1[i] = var_u_b(Ys=Ys_bs_dfbeta1, Xs=Xs_bs,
  #                             pik=pik_bs, p=2, n=sum(s))
  vub_bs_dfbeta2[i] = var_u_b(Ys=Ys_bs_dfbeta2, Xs=Xs_bs,
                              pik=pik_bs, p=2, n=sum(s))
  vub_bs_dfbeta2_w[i] <- vub_bs_dfbeta2[i] *  w #(1/N^2) # *((nrow(mydata_if)-sum(mydata_if$delta))/nrow(mydata_if)) # (1/N^2)
  vub_bs_dfbeta2_pik[i] = var_u_b(Ys=Ys_bs_dfbeta2, Xs=Xs_bs,
                                  pik=pik_bs, p=2, n=sum(s)) * pik_bs[1]
  est_se_score[i] <- sqrt(D3[i] * D3[i] * vub_bs_score[i]) # D1
  # est_se_dfbeta1[i] <- sqrt(vub_bs_dfbeta1[i])
  est_se_dfbeta2[i] <- sqrt(vub_bs_dfbeta2[i])
  # * (n/(N-sum(mydata_if$delta)))
  # * n/N
  est_se_D_score[i] <- sqrt(D2[i] + est_se_score[i]^2) # est_se_score[i]^2*w
  est_se_D_score_w[i] <- sqrt(D2[i] + est_se_score[i]^2*w)
  est_se_two[i] <-  sqrt(phase1_var[i] + vub_bs_dfbeta2[i])
  #est_se_D_pik[i] <-  sqrt(D2[i] + vub_bs_dfbeta2_pik[i])
  est_se_D2_df_w[i] <- sqrt(D2[i] + vub_bs_dfbeta2_w[i])
  est_se_D2_df[i] <- sqrt(D2[i] + vub_bs_dfbeta2[i])
  #est_se_D3_pik[i] <-  sqrt(D3[i] + vub_bs_dfbeta2_pik[i])
  est_se_D3[i] <- sqrt(D3[i] + vub_bs_dfbeta2_w[i])
  
  ci_up_score[i] <- coef(cox_ccs_bs) + qnorm(0.975) * est_se_D2_df[i]
  ci_low_score[i] <- coef(cox_ccs_bs) - qnorm(0.975) * est_se_D2_df[i]
  ci_ind_score[i] <- ifelse((log(2) >= ci_low_score[i]) & (log(2) <= ci_up_score[i]), 1, 0)
}

rr=data.frame(v=c("mean of beta srs","sd of beta srs",
                  "mean of beta bs","sd of beta bs",
                  "mean of beta full","sd of beta full",
                  "mean of beta cali","sd of beta cali",
                  "mean of beta bs cali","sd of beta bs cali",
                  "mean of se srs","sd of se srs",
                  "mean of se bs","sd of se bs",
                  "mean of robust se srs","sd of robust se srs",
                  "mean of robust se bs","sd of robust se bs",
                  "mean of est se (score)","sd of est se (score)",
                  "mean of est se (dfbeta2)","sd of est se (dfbeta2)",
                  "mean of two phase est se (full + vub2)","sd of two phase est se (full + vub2)",
                  "mean of two phase est se (svy1 + vub2)","sd of two phase est se (svy1 + vub2)",
                  "mean of two phase est se (D + vub_score)","sd of two phase est se (D + vub_score)",
                  "mean of two phase est se (D + vub_score*w)","sd of two phase est se (D + vub_score*w)",
                  "mean of two phase est se (D + vub2*w)","sd of two phase est se (D + vub2*w)",
                  "mean of two phase est se (D2 + vub2)","sd of two phase est se (D2 + vub2)",
                  "mean of two phase est se (D3 + vub2*w)","sd of two phase est se (D3 + vub2*w)",
                  "mean of sqrt D","sd of sqrt D",
                  "mean of D","sd of D",
                  "mean of sqrt D3","sd of sqrt D3",
                  "mean of D3","sd of D3",
                  "mean of info", "sd of info",
                  "ratio of CI log2"),
              n=c(mean(beta_srs),sd(beta_srs),
                  mean(beta_bs),sd(beta_bs),
                  mean(full_beta),sd(full_beta),
                  mean(beta_cali),sd(beta_cali),
                  mean(beta_bs_cali),sd(beta_bs_cali),
                  mean(se_srs),sd(se_srs),
                  mean(se_bs),sd(se_bs),
                  mean(robust_se_srs),sd(robust_se_srs),
                  mean(robust_se_bs),sd(robust_se_bs),
                  mean(est_se_score),sd(est_se_score),
                  mean(est_se_dfbeta2),sd(est_se_dfbeta2),
                  mean(sqrt(var(full_beta) + vub_bs_dfbeta2*w)),sd(sqrt(var(full_beta) + vub_bs_dfbeta2*w)),
                  mean(est_se_two),sd(est_se_two),
                  mean(est_se_D_score),sd(est_se_D_score),
                  mean(est_se_D_score_w),sd(est_se_D_score_w),
                  mean(est_se_D2_df_w),sd(est_se_D2_df_w),
                  mean(est_se_D2_df),sd(est_se_D2_df),
                  mean(est_se_D3),sd(est_se_D3),
                  mean(sqrt(D2)),sd(sqrt(D2)),
                  mean(D2),sd(D2),
                  mean(sqrt(D3)),sd(sqrt(D3)),
                  mean(D3),sd(D3),
                  mean(info1),sd(info1),
                  sum(ci_ind_score)/length(ci_ind_score)))
(rr=data.frame(n=rr$v,n=round(rr$n,4)))


write_xlsx(rr, paste("c1s1_ccs_conti_N", as.character(n_pop),"_n",as.character(n),"cor",as.character(rho),".xlsx", sep = ""))


write_xlsx(data.frame(beta_srs), paste("random_ccs_conti_beta_srs_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(beta_bs), paste("random_ccs_conti_beta_bs_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(beta_cali), paste("random_ccs_conti_beta_cali_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(beta_bs_cali), paste("random_ccs_conti_beta_cali_bs_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(est_se_D2_df), paste("random_ccs_conti_est_se_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(sqrt(D2)), paste("random_ccs_conti_se1_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(est_se_dfbeta2), paste("random_ccs_conti_se2_N_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))
write_xlsx(data.frame(full_beta), paste("random_ccs_conti_beta_full_N", as.character(n_pop),"_n",as.character(n),"_corr",as.character(rho), ".xlsx", sep = ""))

# mean of sqrt D: phase 1
# mean of est se (dfbeta2): phase 1
# mean of two phase est se (D2 + vub2): est se
