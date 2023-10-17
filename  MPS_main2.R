getwd()
setwd("/Users/yejiko/Downloads/PhD/MPS")

### External requirement
library(dplyr)
library(tidyr)
library(survival)
library(ggplot2)
library(ggpubr)

### Simulation set-up
lambda_0 = 0.03 # parameter for initial time distribution
psi = -0.3 # -0.3, 0, or 0.3
ns = 1000 # sample size
tpoints = 20 # number of time points

msm_est <- c()
cph_est <- c()
cph_est2 <- c()

msm_se <- c()
cph_se <- c()
cph_se2 <- c()

nsim = 500 # number of simulations

set.seed(1234)

while(length(msm_est) < nsim){
  
  time = data.frame(matrix(NA, tpoints, ns))
  time0 = rexp(ns, lambda_0)
  time[1,] = time0
 
  Larray = data.frame(matrix(NA, tpoints, ns)) # time-varying confounder
  Aarray = data.frame(matrix(NA, tpoints, ns)) # exposure
  binAarray = data.frame(matrix(NA, tpoints, ns)) # binary exposure
  cumAarray = data.frame(matrix(NA, tpoints, ns)) # cumulative exposure
  Yarray = data.frame(matrix(NA, tpoints, ns)) # outcome (survival status)
  chkarray = data.frame(matrix(NA, tpoints, ns))
  Ymarray = data.frame(matrix(NA, tpoints, ns)) 
  Am1array = data.frame(matrix(NA, tpoints, ns)) 
  Am1array[1,] = 0
  
  ### j = 1 (baseline time)
  
  ### Generate total duration of employment
  max = ifelse(time0<20, floor(time0),20)
  max[max==0]=1
  worktime0 = runif(ns, 1, max)
  len = worktime0[1]
  
  ### Generate L
  Larray[1,] = 1
  
  ### Generate A
  Aarray[1,] <- rlnorm(ns, mean = log(0.002 - 0.00005*1), 1) # intercept: 0.001 and slope: -0.000015, var= 0.01
                         
  cumAarray[1,] <- Aarray[1,]
  binAarray[1,] <- ifelse(cumAarray[1,]<0.01, 0,1)
  
  ### Generate Y
  chkarray[1,] = exp(psi * binAarray[1,] )
  Yarray[1,] <- ifelse((time0 > chkarray[1,]), 0, 1)
  
  ### Generate T
  col = which(time0 <= chkarray[1,])
  time[1, col] = time0[col] * exp(-(psi * binAarray[1,col] ))
  
  ### j = 2 to tpoints
  
  for(j in (2:tpoints)){

    ### Generate A
    Aarray[j,] <- rlnorm(ns, mean = log(0.002 - 0.00005*j), 1)
    cumAarray <- cumsum(Aarray)
    binAarray[j,] <- ifelse(cumAarray[j,]<0.01, 0,1)
    
    ### Generate Y
    Ymarray[j,] = Yarray[j-1,] 
    Larray[j, which(Ymarray[j,]== 1)] = 0
    binAarray[j, which(Ymarray[j,]== 1)] = 0
    chkarray[j,] = chkarray[j-1,] + exp(psi * binAarray[j,])
    
    for(i in 1:ns){
      ### Generate L
      Larray[j,i] = ifelse(worktime0[i] >= j,j,NA)
      
      Am1array[j,i] = binAarray[j-1,i]
      
      if(time0[i] > chkarray[j,i]){
        Yarray[j,i] = 0
        time[j, i] = time[j-1, i]
      }
      else if (time0[i] <= chkarray[j,i]){
        Yarray[j,i] = 1
        ### Generate T
        # if Ymarray[j,i] == 1
        if(Ymarray[j,i] == 1){
          time[j, i] = time[j-1, i]
        }
        # if Ymarray[j,i] != 1
        else{
          time[j, i] = (j-1) + (time0[i] - chkarray[j-1,i]) * exp(-(psi * binAarray[j,i])) 
        }
      }
    }
  }
  
  ### Verify the simulation data by building a cox MSM
  merge <- data.frame(matrix(NA, ns*tpoints, 2))
  merge[[1]] <- rep(1:tpoints,ns)
  merge[[2]] <- rep(1:ns, each = tpoints)
  names(merge) <- c("time", "id")
  
  Larray <- gather(Larray)
  binAarray <- gather(binAarray)
  Am1array <- gather(Am1array)
  Yarray <- gather(Yarray)
  
  ### Last row of time data frame is the updated time
  survtime <- c(unlist(time[20,]))
  merge$timelen0 <- rep(survtime,each=20)
  merge$timelen0 <- ifelse(merge$timelen0<1, merge$timelen0+1, merge$timelen0)
  
  merge$L <- Larray$value
  merge$A <- binAarray$value
  merge$Am1 <- Am1array$value
  merge$Y <- Yarray$value
  merge <- merge %>% group_by(id) %>% mutate(timelen1 = ifelse(timelen0 > tpoints, tpoints, max(timelen0)))
  
  ### Collapse follow-up time (from retirement to death)
  merge2 <- merge %>% group_by(id) %>% filter(time <= ceiling(timelen1))

  final <- merge2 %>% group_by(id) %>% fill(L, .direction = 'down')
  final <- final %>% group_by(id) %>% mutate(Lm1=lag(L, default=0)) %>% mutate(Lm2=lag(Lm1, default=0))
  
  ### Create time-dependent censoring status (censored if time=tpoints and Y=0) and time entry
  final <- final %>% group_by(id) %>% mutate(is.censored = ifelse((time==tpoints &Y==0),1,0), time0=lag(time, default=0))
  
  
  ### 1. Compute stabilized weights
  mod <- glm(A ~ time, data = final, family = "binomial")

  # calculate propensity score
  final$ps_treatment = predict(mod, type = "response")
  mod_stab <- glm(A ~ L + time, data = final, family = "binomial")
  summary(mod_stab)
  
  # calculate propensity score - stabilized weights
  final$ps_treatment_stabilized = predict(mod, type = "response")
  
  final$prob_ratio_treatment = final$ps_treatment/final$ps_treatment_stabilized
  
  final <- final %>% group_by(id) %>% arrange(id, time) %>% mutate(cumprod_treatment = cumprod(prob_ratio_treatment))

  final$prob_ratio_treatment[which(is.na(final$prob_ratio_treatment))] <- 1
  
  ### 2. Compute censoring weights
  mod <- glm(is.censored ~ time, 
             data = final, family = "binomial")
  censoring_ps = predict(mod, type = "response")
  
  ### Denominator model includes baseline and time-varying covariates 
  mod_stab <- glm(is.censored ~ Lm1 + time, data = final, family = "binomial")
  stabilized_censoring_ps = predict(mod_stab, type="response")
  final$censoring_prob_ratio = censoring_ps/stabilized_censoring_ps
  
  final <- final %>% group_by(id) %>% 
    arrange(id, time) %>% mutate(cumprod_censoring = cumprod(censoring_prob_ratio))
  
  ### 3. Get the final weights
  final$final_weights <- final$cumprod_treatment* final$cumprod_censoring
  
  # truncation at 95% percentile 
  final$final_weights[final$final_weights>quantile(final$final_weights,0.95)] = quantile(final$final_weights,0.95)
    
  ### 4. Fit a cox MSM
  msm <- coxph(Surv(time0, time, Y) ~ A, data = final, weights = final_weights)
  msm_est <- c(msm_est, msm$coefficients)
  msm_se <- c(msm_se, summary(msm)$coef[1,4]) # robust se
  #msm_ci <- rbind(msm_ci, confint(msm))
  
  print(paste(length(msm_est), "th estimate generated!", sep = ""))
  
  ### 5. Fit a standard cox PH
  cph <- coxph(Surv(time0, time, Y) ~ A + L + Am1 + Lm1 + Lm2, data = final)
  cph_est <- c(cph_est, cph$coefficients[1])
  cph_se <- c(cph_se, summary(cph)$coef[1,3])

  ### 6. Compare with a cox PH with no time-dependent variable
  cph2 <- coxph(Surv(time0, time, Y) ~ A + Am1, data = final)
  cph_est2 <- c(cph_est2, cph2$coefficients[1])
  cph_se2 <- c(cph_se2, summary(cph2)$coef[1,3])
}

write.csv(cbind(msm_est[msm_est<1 &msm_est>-1], msm_se[msm_est<1 &msm_est>-1], 
          cph_est[1:length(cph_est)], cph_se[1:length(cph_se)], 
          cph_est2[1:length(cph_est2)], cph_se2[1:length(cph_se2)]), 
          "MPS_sim3_update.csv")

### Bias
round(mean(cph_est2[1:length(cph_est2)]),4)-psi # MC mean for cox PH with no L
round(mean(cph_est[1:length(cph_est)]),4)-psi # MC mean for standard cox
round(mean(msm_est[msm_est<1 &msm_est>-1]),4)-psi # MC mean for cox MSM

### MSE
round(mean((cph_est2[1:length(cph_est2)]-psi)^2),4)
round(mean((cph_est[1:length(cph_est)]-psi)^2),4)
round(mean((msm_est[msm_est<1 &msm_est>-1]-psi)^2),4)

### Coverage rate
# unadjusted
lb =  cph_est2[1:length(cph_est2)] - 1.96 *cph_se2[1:length(cph_est2)]
ub =  cph_est2[1:length(cph_est2)] + 1.96 *cph_se2[1:length(cph_est2)]
cover <- ifelse(lb <psi & ub >psi,1,0)
mean(cover)

# adjusted
lb =  cph_est[1:length(cph_est)] - 1.96 *cph_se[1:length(cph_est)]
ub =  cph_est[1:length(cph_est)] + 1.96 *cph_se[1:length(cph_est)]
cover <- ifelse(lb <psi & ub >psi,1,0)
mean(cover)

# msm
lb =  msm_est[msm_est<1 &msm_est>-1] - 1.96 *msm_se[msm_est<1 &msm_est>-1]
ub =  msm_est[msm_est<1 &msm_est>-1] + 1.96 *msm_se[msm_est<1 &msm_est>-1]
cover <- ifelse(lb <psi & ub >psi,1,0)
mean(cover)

### Create a summary table
sim_res <- data.frame(matrix(NA,9,nsim))

rownames(sim_res) <- c(paste(c("coxPH w/o L: psi", "coxPH with L: psi", 
                               "coxMSM with L: psi"), 0.3), 
                       paste(c("coxPH w/o L: psi", "coxPH with L: psi", 
                               "coxMSM with L: psi"), 0), 
                       paste(c("coxPH w/o L: psi", "coxPH with L: psi", 
                               "coxMSM with L: psi"), -0.3))
sim_res[4,] <- cph_est2
sim_res[5,] <- cph_est
sim_res[6,] <- msm_est

sim_dat <- data.frame(matrix(NA,9,3))
colnames(sim_dat) <- c("mean", "LCI","UCI")
sim_dat[1] <- rowMeans(sim_res)

for(i in 1:9){
  sim_dat[i,2] <- sim_dat[i,1] - 1.96 *sd(sim_res[i,])
  sim_dat[i,3] <- sim_dat[i,1] + 1.96 *sd(sim_res[i,])
}

sim_dat$mod <- rownames(sim_res) 
sim_dat$psi <- c(rep(0.3,3),rep(0,3),rep(-0.3,3))
sim_dat$bias <- (sim_dat$mean - sim_dat$psi) 
sim_dat$bias_lci <- (sim_dat$LCI - sim_dat$psi) 
sim_dat$bias_uci <- (sim_dat$UCI - sim_dat$psi) 

g1 <- ggplot(sim_dat[sim_dat$psi==0.3,], aes(x=mod, y=bias)) +
  ggtitle(expression(paste(psi, " = 0.3"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(NULL, limits = sim_dat[sim_dat$psi==0.3,]$mod,
                   labels = c("Model 3", "Model 2", "Model 1")) +
  scale_y_continuous("Bias") +
  geom_point(shape=15, color="red",size = 3) +
  geom_errorbar(aes(ymin = bias_lci, ymax = bias_uci)) + 
  coord_flip() + 
  geom_text(aes(label=round(bias,4)), hjust=0, vjust=2)+ mytheme()

g2 <- ggplot(sim_dat[sim_dat$psi==0,], aes(x=mod, y=bias)) + 
  ggtitle(expression(paste(psi, " = 0"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(NULL, limits = sim_dat[sim_dat$psi==0,]$mod,
                   labels = c("Model 3", "Model 2", "Model 1")) +
  scale_y_continuous("Bias") +
  geom_point(shape=15, color="red",size = 3) +
  geom_errorbar(aes(ymin = bias_lci, ymax = bias_uci)) + 
  coord_flip() + 
  geom_text(aes(label=round(bias,4)), hjust=0, vjust=2)+ mytheme()

g3 <- ggplot(sim_dat[sim_dat$psi==-0.3,], aes(x=mod, y=bias)) +    
  ggtitle(expression(paste(psi, " = -0.3"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(NULL, limits = sim_dat[sim_dat$psi==-0.3,]$mod, 
                   labels = c("Model 3", "Model 2", "Model 1")) +
  scale_y_continuous("Bias") +
  geom_point(shape=15, color="red",size = 3) +
  geom_errorbar(aes(ymin = bias_lci, ymax = bias_uci)) + 
  coord_flip() + 
  geom_text(aes(label=round(bias,4)), hjust=0, vjust=2) + mytheme()

plot = ggarrange(g1,g2,g3, ncol = 3)
annotate_figure(plot, bottom = text_grob("Figure 1. Comparison of model bias", 
                                      color = "black",size = 14))
