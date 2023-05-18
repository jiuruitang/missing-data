# Created: 04/04/2022
# Modified: 04/04/2022
# Change educ variable to three levels

# read in data
cps18_test_data <- read.csv("cps18_with_edu.csv")
nc_dat <- cps18_test_data[,-1]
summary(nc_dat) 
library(nnet)
library(mvtnorm)
library(survey)
library(gtools)

nc_dat$educ <- nc_dat$educ+1

# Pre-processing
# adjust weights, down weight unit response weights, assign weights to unit nonresponse
nc_adj <- nc_dat
total_W <- sum(nc_dat$weight,na.rm = T)
nc_adj[which(nc_adj$unit == 0),]$weight <- nc_adj[which(nc_adj$unit == 0),]$weight*(1-sum(nc_adj$unit)/nrow(nc_adj))
nc_adj[which(nc_adj$unit == 1),]$weight <- total_W/nrow(nc_adj)
nc_adj$pi <- 1/nc_adj$weight

# add item nonresponse indicator
nc_adj$Rv <-  ifelse(is.na(nc_adj$voted), 1, 0)
nc_adj$Rg <- ifelse(is.na(nc_adj$sex), 1, 0)
nc_adj$Ra <- ifelse(is.na(nc_adj$age), 1, 0)
nc_adj$Re <- ifelse(is.na(nc_adj$race), 1, 0)
nc_adj$Rc <- ifelse(is.na(nc_adj$educ), 1, 0)

nc_adj[which(nc_adj$unit == 1), c("Rv","Rg","Ra","Re","Rc")] <- NA

# Get HT estimators for gender
## For gender, male is 0, female is 1
T_sex <- total_W*0.524 # 3900666
Sd_sex <- sqrt(sum((nc_adj$sex/nc_adj$pi)^2*(1-nc_adj$pi),na.rm = T)) #84927.75

# HT estimators for race
## For race, white non-his 69.9%
T_raceW <- total_W*0.699
T_raceB <- total_W*0.218
T_raceH <- total_W*0.039
T_raceR <- total_W*0.044
raceW <- ifelse(nc_adj$race==1,1,0)
Sd_raceW <- sqrt(sum((raceW/nc_adj$pi)^2*(1-nc_adj$pi),na.rm = T)) #93876.83
## black non-his 21.8%
raceB <- ifelse(nc_adj$race==2,1,0)
Sd_raceB <- sqrt(sum((raceB/nc_adj$pi)^2*(1-nc_adj$pi),na.rm = T)) #60796.61
## hispanic or latino 3.9%
raceH <- ifelse(nc_adj$race==3,1,0)
Sd_raceH <- sqrt(sum((raceH/nc_adj$pi)^2*(1-nc_adj$pi),na.rm = T)) #22426.67
## rest 4.4%
raceR <- ifelse(nc_adj$race==4,1,0)
Sd_raceR <- sqrt(sum((raceR/nc_adj$pi)^2*(1-nc_adj$pi),na.rm = T)) #25918.68

# first assign 2 missing everything except weight as unit non response
nc_adj$unit[which(is.na(nc_adj$sex) & nc_adj$unit == 0)] <- 1

# relevel age
# nc_adj$age <- as.factor(nc_adj$age)
# nc_adj <- within(nc_adj, age <- relevel(age, ref = "(49,59]"))

# impute missing race assuming MAR
m_eth <- multinom(race ~ as.factor(sex) , data = nc_adj[which(nc_adj$Re==0 & nc_adj$unit == 0),])
nc_adj$race[which(nc_adj$Re == 1 & nc_adj$unit == 0)] <- predict(m_eth,newdata = data.frame(sex=nc_adj$sex[which(nc_adj$Re == 1 & nc_adj$unit == 0)]),type = "class")

# impute missing educ assuming MAR
m_educ <- multinom(educ ~ as.factor(sex) + as.factor(race), data = nc_adj[which(nc_adj$Rc==0 & nc_adj$unit == 0),],family=binomial())
# p_educ <- predict(m_educ,newdata = data.frame(sex=nc_adj$sex[which(nc_adj$Rc == 1 & nc_adj$unit == 0)],race=nc_adj$race[which(nc_adj$Rc == 1 & nc_adj$unit == 0)]),type = "response")
nc_adj$educ[which(nc_adj$Rc == 1 & nc_adj$unit == 0)] <- 
  predict(m_educ, newdata = data.frame(sex=nc_adj$sex[which(nc_adj$Rc == 1 & nc_adj$unit == 0)],
                                       race = nc_adj$race[which(nc_adj$Rc == 1 & nc_adj$unit == 0)]),type = "class")

# impute missing age assuming MAR
m_age <- multinom(age ~ as.factor(sex) + as.factor(race) + as.factor(educ), data = nc_adj[which(nc_adj$Ra==0 & nc_adj$unit == 0),])
nc_adj$age[which(nc_adj$Ra == 1 & nc_adj$unit == 0)] <- as.character(predict(m_age,newdata = data.frame(sex=nc_adj$sex[which(nc_adj$Ra == 1 & nc_adj$unit == 0)],
                                                                                                        race=nc_adj$race[which(nc_adj$Ra == 1 & nc_adj$unit == 0)],
                                                                                                        educ=nc_adj$educ[which(nc_adj$Ra == 1 & nc_adj$unit == 0)]),type = "class"))
# impute missing vote assuming MAR
n_mis_V <- sum(is.na(nc_adj$voted[which(nc_adj$unit == 0)]))
m_vote <- glm(voted ~ as.factor(sex) + age + as.factor(race) + as.factor(educ) + as.factor(sex)*as.factor(race)
              + as.factor(sex)*age + as.factor(sex)*as.factor(educ), 
              data = nc_adj[which(nc_adj$Rv==0 & nc_adj$unit == 0),],family=binomial())
p_vote <- predict(m_vote,newdata = data.frame(sex=nc_adj$sex[which(nc_adj$Rv == 1 & nc_adj$unit == 0)],
                                              race=nc_adj$race[which(nc_adj$Rv == 1 & nc_adj$unit == 0)],
                                              age=nc_adj$age[which(nc_adj$Rv == 1 & nc_adj$unit == 0)],
                                              educ=nc_adj$educ[which(nc_adj$Rv == 1 & nc_adj$unit == 0)]),type = "response")
nc_adj$voted[which(nc_adj$Rv == 1 & nc_adj$unit == 0)] <- rbinom(length(p_vote),1,prob = p_vote)

# HT estimators for vote
T_vote <- total_W*0.49 # 3647570
Sd_vote <- sqrt(sum((nc_adj$voted/nc_adj$pi)^2*(1-nc_adj$pi),na.rm = T)) 
# 83631.59 for Rv == 0, U == 0
# 91908.18 for U == 0 after imputing missing vote

## look at the weights of two missing sex obs
#nc_dat[which(nc_dat$unit == 0 & is.na(nc_dat$sex)),]
#hist(nc_dat$weight)
#abline(v = 4402.487, col = "red")
#abline(v = 5805.686, col = "red")

# update
## initialize alpha's
alpha_gs <- rnorm(1)
alpha_es <- rmvnorm(1, mean = rep(0,prod(dim(coef(m_eth)))), sigma = diag(prod(dim(coef(m_eth)))))
alpha_as <- rmvnorm(1, mean = rep(0,prod(dim(coef(m_age)))), sigma = diag(prod(dim(coef(m_age)))))
alpha_vs <- rmvnorm(1, mean = rep(0,length(coef(m_vote))), sigma = diag(length(coef(m_vote))))
alpha_cs <- rmvnorm(1, mean = rep(0,length(coef(m_educ))), sigma = diag(length(coef(m_educ))))

## Missing indicator models
r_eth <- glm(Re~as.factor(sex) + as.factor(educ) + as.factor(age) + as.factor(voted), data = nc_adj[which(nc_adj$unit == 0),],family=binomial())
r_educ <- glm(Rc~as.factor(sex) + as.factor(race) + as.factor(age) +  as.factor(voted),data = nc_adj[which(nc_adj$unit == 0),],family=binomial() )
r_age <- glm(Ra ~ as.factor(sex) + as.factor(race) + as.factor(educ), data = nc_adj[which(nc_adj$unit == 0),],family=binomial())
r_vote <- glm(Rv~as.factor(sex) + as.factor(race) + as.factor(educ) + as.factor(age), data = nc_adj[which(nc_adj$unit == 0),],family=binomial())

## initialize gamma's
gamma_es <- rnorm(length(coef(r_eth)))
gamma_cs <- rnorm(length(coef(r_educ)))
gamma_as <- rnorm(length(coef(r_age)))
gamma_vs <- rnorm(length(coef(r_vote)))

## initialize nu
nus <- rnorm(6)

# MCMC
nc_obs <- nc_adj[which(nc_adj$unit == 0 & nc_adj$Rg == 0),]
ps1 <- sum(nc_obs$W[which(nc_obs$sex == 1 & nc_obs$Y == 1)])/sum(nc_obs$W)

b0_alpha_g <- alpha_gs; b0_alpha_g[] <- 0 
sigma0_alpha_g <- 1
b0_alpha_e <- alpha_es; b0_alpha_e[] <- 0 
sigma0_alpha_e <- diag(c(alpha_es)); diag(sigma0_alpha_e) <- 1
b0_alpha_a <- alpha_as; b0_alpha_a[] <- 0 
sigma0_alpha_a <- diag(c(alpha_as)); diag(sigma0_alpha_a) <- 1
b0_alpha_v <- alpha_vs; b0_alpha_v[] <- 0 
sigma0_alpha_v <- diag(c(alpha_vs)); diag(sigma0_alpha_v) <- 1

b0_gamma_e <- gamma_es; b0_gamma_e[] <- 0 
sigma0_gamma_e <- diag(c(gamma_es)); diag(sigma0_gamma_e) <- 1
b0_gamma_a <- gamma_as; b0_gamma_a[] <- 0 
sigma0_gamma_a <- diag(c(gamma_as)); diag(sigma0_gamma_a) <- 1
b0_gamma_v <- gamma_vs; b0_gamma_v[] <- 0 
sigma0_gamma_v <- diag(c(gamma_vs)); diag(sigma0_gamma_v) <- 1

b0_nu <- nus; b0_nu[] <- 0 
sigma0_nu <- diag(c(nus)); diag(sigma0_nu) <- 1


niter <- 10000
burnin <- 5000
GAMMA_e <- matrix(NA,niter-burnin,length(gamma_es))
GAMMA_c <- matrix(NA,niter-burnin,length(gamma_cs))
GAMMA_a <- matrix(NA,niter-burnin,length(gamma_as))
GAMMA_v <- matrix(NA,niter-burnin,length(gamma_vs))
ALPHA_e <- matrix(NA,niter-burnin,length(alpha_es))
ALPHA_c <- matrix(NA,niter-burnin,length(alpha_cs)) # educ
ALPHA_a <- matrix(NA,niter-burnin,length(alpha_as))
ALPHA_v <- ALPHA_v_origin <- matrix(NA,niter-burnin,length(alpha_vs))
ALPHA_g <- ALPHA_g_origin <- matrix(NA,niter-burnin,length(alpha_gs))

vote_impute <- vote_replicate <- matrix(NA,niter-burnin,length(nc_adj$voted))
age_impute <- age_replicate <- matrix(NA,niter-burnin,length(nc_adj$age))
race_impute <- race_replicate <- matrix(NA,niter-burnin,length(nc_adj$race))
educ_impute <- educ_replicate <- matrix(NA,niter-burnin,length(nc_adj$educ))
sex_impute <- sex_replicate <- matrix(NA,niter-burnin,length(nc_adj$sex))
Ra_replicate <- Rv_replicate <- Re_replicate <- Rc_replicate <- matrix(NA,niter-burnin,sum(nc_adj$unit == 0))

# impute sex and race for U == 1 as a starting point
nc_adj$sex[which(nc_adj$unit == 1)] <- rbinom(sum(nc_adj$unit),1,0.524)
nc_adj$race[which(nc_adj$unit == 1)] <- sample(c(1:4),sum(nc_adj$unit),prob = c(0.699,0.218,0.039,0.044),replace = TRUE)
nc_adj$educ[which(nc_adj$unit == 1)] <- predict(m_educ,newdata = nc_adj[which(nc_adj$unit == 1),],type = "class")
nc_adj$age[which(nc_adj$unit == 1)] <- as.character(predict(m_age,newdata = nc_adj[which(nc_adj$unit == 1),],type = "class"))
# nc_adj$voted[which(nc_adj$unit == 1)] <- rbinom(sum(nc_adj$unit),1,0.49)

vote_unit_prob <- predict(m_vote,newdata = nc_adj[which(nc_adj$unit == 1),],type = "response")
vote_origin <- rbinom(sum(nc_adj$unit),1,vote_unit_prob)
vote_unit_prob_mod <- vote_unit_prob - min(vote_unit_prob)
vote_unit_prob_mod <- ifelse(vote_unit_prob_mod - 0.2 > 0,vote_unit_prob_mod - 0.2,0)
vote_mod <- rbinom(sum(nc_adj$unit),1,vote_unit_prob_mod)
nc_adj$voted[which(nc_adj$unit == 1)] <- vote_mod
sum(nc_adj$weight[which(nc_adj$voted == 1)])

# nc_copy <- nc_adj 
nc_adj <- nc_copy

logit_inv <- function(x){
  return(exp(x)/(1+exp(x)))
}

nc_rep <- nc_adj

for (i in 1:niter){
  print(i)
  {
    n_unit0 <- sum(nc_adj$unit == 0)
    n_unit1 <- sum(nc_adj$unit == 1)
    ### update alpha's
    ## update alpha_g
    m_sex <- glm(sex~1,data = nc_adj[which(nc_adj$unit == 0),],family=binomial())
    sex_mat <- model.matrix(m_sex)
    alpha_g_origin <- alpha_g <- rmvnorm(1,mean = coef(m_sex),sigma = vcov(m_sex))
    
    # replicate data for U==0 cases
    sex_prob <- logit_inv(sex_mat%*%alpha_g)
    sex_rep <- nc_adj$sex
    # sex_rep[which(nc_adj$unit == 0)] <- rbinom(n_unit0,1,prob = sex_prob)
    nc_rep$sex[which(nc_rep$unit == 0)] <- rbinom(n_unit0,1,prob = sex_prob)
    
    # impute sex for U==1 cases
    T_sex_hat <- rnorm(1,T_sex,Sd_sex)
    total_one_sex <- floor((T_sex_hat - sum((nc_adj$sex*nc_adj$weight)[which(nc_adj$unit == 0)]))/mean(nc_adj$weight[which(nc_adj$unit==0)]))
    m_sex1 <- glm(sex~1,data = nc_adj[which(nc_adj$unit == 1),],family=binomial())
    sex_mat1 <- model.matrix(m_sex1)
    # coef in front of sex
    u_sex <- logit(total_one_sex/n_unit1) - mean(sex_mat1%*%t(alpha_g))
    alpha_g[1] <- alpha_g[1] + u_sex
    sex_prob <- logit_inv(sex_mat1%*%alpha_g[1])
    # nc_adj$sex[which(nc_adj$unit == 1)] <- sex_rep[which(nc_adj$unit == 1)] <- rbinom(n_unit1,1,prob = sex_prob)
    nc_adj$sex[which(nc_adj$unit == 1)] <- nc_rep$sex[which(nc_rep$unit == 1)] <- rbinom(n_unit1,1,prob = sex_prob)
    
    ## update alpha_e
    m_eth <- multinom(race ~ as.factor(sex) , data = nc_adj[which(nc_adj$unit == 0),]) # base level is white(1)
    alpha_e <- rmvnorm(1,mean = c(t(coef(m_eth))),sigma = vcov(m_eth))
    eth_mat <- model.matrix(m_eth)
    eth_rep <- nc_adj$race
    alpha_e_mat <- matrix(alpha_e,nrow = 3, byrow = T)
    
    # generate replicate data for eth|sex
    m_eth_rep <- multinom(race ~ as.factor(sex) , data = nc_rep[which(nc_rep$unit == 0),]) # base level is white(1)
    eth_mat_rep <- model.matrix(m_eth_rep)
    eth_prob_rep <-  predict(m_eth_rep,type = "prob")
    eth_prob_alpha <- eth_prob_rep
    eth_prob_alpha[,1] <- 1
    for (l in 1:3){
      eth_prob_alpha[,l+1] <- logit_inv(eth_mat_rep %*%(alpha_e_mat[l,]))
      eth_prob_alpha[,l+1] <- eth_prob_alpha[,l+1]/(1-eth_prob_alpha[,l+1])
    }
    eth_prob_alpha_adj <- eth_prob_alpha/rowSums(eth_prob_alpha)
    nc_rep$race[which(nc_rep$unit == 0)] <- 
      sapply(1:n_unit0,function(x)sample(c(1,2,3,4),size = 1, 
                                         replace = FALSE, prob = eth_prob_alpha_adj[x,]))
    
    # sample unit nonresponse for race
    m_eth1 <- multinom(race ~ as.factor(sex) , data = nc_adj[which(nc_adj$unit == 1),])
    eth_mat1 <- model.matrix(m_eth1)
    black_prob <- logit_inv(eth_mat1 %*%(alpha_e_mat[1,]))
    hist_prob <- logit_inv(eth_mat1 %*%(alpha_e_mat[2,]))
    rest_prob <- logit_inv(eth_mat1 %*%(alpha_e_mat[3,]))
    white_prob <- rep(1,n_unit1)
    ################# Added 01/07 ##################
    black_prob <- black_prob/(1-black_prob)
    hist_prob <- hist_prob/(1-hist_prob)
    rest_prob <- rest_prob/(1-rest_prob)
    ##########################################
    black_prob_adj <- black_prob/(black_prob+hist_prob+rest_prob+white_prob)
    hist_prob_adj <- hist_prob/(black_prob+hist_prob+rest_prob+white_prob)
    rest_prob_adj <- rest_prob/(black_prob+hist_prob+rest_prob+white_prob)
    white_prob_adj <- white_prob/(black_prob+hist_prob+rest_prob+white_prob)
    
    model_white_prob <- mean(white_prob_adj)
    model_black_prob <- mean(black_prob_adj)
    model_hist_prob <- mean(hist_prob_adj)
    model_rest_prob <- mean(rest_prob_adj)
    
    # T_white_hat <- rnorm(1,T_raceW,Sd_raceW)
    T_black_hat <- rnorm(1,T_raceB,Sd_raceB)
    T_hist_hat <- rnorm(1,T_raceH,Sd_raceH)
    T_rest_hat <- rnorm(1,T_raceR,Sd_raceR)
    
    imp_one_rest <- max(floor((T_rest_hat - sum(nc_adj$weight[which(nc_adj$unit == 0 & nc_adj$race == 4)]))/mean(nc_adj$weight[which(nc_adj$unit==0)])),0)
    imp_one_black <- max(floor((T_black_hat - sum(nc_adj$weight[which(nc_adj$unit == 0 & nc_adj$race == 2)]))/mean(nc_adj$weight[which(nc_adj$unit==0)])),0)
    imp_one_hist <- max(floor((T_hist_hat - sum(nc_adj$weight[which(nc_adj$unit == 0 & nc_adj$race == 3)]))/mean(nc_adj$weight[which(nc_adj$unit==0)])),0)
    imp_rest_prob <- imp_one_rest/n_unit1
    imp_black_prob <- imp_one_black/n_unit1
    imp_hist_prob <- imp_one_hist/n_unit1
    imp_white_prob <- (n_unit1 - imp_one_rest - imp_one_black - imp_one_hist)/n_unit1
    
    ################ Edited 01/07 ################
    imp_black_prob_overW <- imp_black_prob/imp_white_prob
    imp_hist_prob_overW <- imp_hist_prob/imp_white_prob
    imp_rest_prob_overW <- imp_rest_prob/imp_white_prob
    
    black_adj <- logit(imp_black_prob_overW/(1+imp_black_prob_overW)) - mean(eth_mat1 %*%(alpha_e_mat[1,]))
    hist_adj <- logit(imp_hist_prob_overW/(1+imp_hist_prob_overW)) - mean(eth_mat1 %*%(alpha_e_mat[2,]))
    rest_adj <- logit(imp_rest_prob_overW/(1+imp_rest_prob_overW)) - mean(eth_mat1 %*%(alpha_e_mat[3,]))
    ################################################
    
    # change intercepts
    alpha_e_mat[1] <- alpha_e_mat[1] + black_adj
    alpha_e_mat[2] <- alpha_e_mat[2] + hist_adj
    alpha_e_mat[3] <- alpha_e_mat[3] + rest_adj
    
    # recalculate
    black_prob <- logit_inv(eth_mat1 %*%(alpha_e_mat[1,]))
    hist_prob <- logit_inv(eth_mat1 %*%(alpha_e_mat[2,]))
    rest_prob <- logit_inv(eth_mat1 %*%(alpha_e_mat[3,]))
    white_prob <- rep(1,n_unit1)
    ################# Added 01/07 ##################
    black_prob <- black_prob/(1-black_prob)
    hist_prob <- hist_prob/(1-hist_prob)
    rest_prob <- rest_prob/(1-rest_prob)
    ##########################################
    black_prob_adj <- black_prob/(black_prob+hist_prob+rest_prob+white_prob)
    hist_prob_adj <- hist_prob/(black_prob+hist_prob+rest_prob+white_prob)
    rest_prob_adj <- rest_prob/(black_prob+hist_prob+rest_prob+white_prob)
    white_prob_adj <- white_prob/(black_prob+hist_prob+rest_prob+white_prob)
    
    # sample race for U == 1
    race_prob <- matrix(NA,n_unit1,4)
    race_prob[,1] <- white_prob_adj
    race_prob[,2] <- black_prob_adj
    race_prob[,3] <- hist_prob_adj
    race_prob[,4] <- rest_prob_adj
    
    nc_adj$race[which(nc_adj$unit == 1)] <- nc_rep$race[which(nc_rep$unit == 1)] <- 
      sapply(1:n_unit1,function(x)sample(c(1,2,3,4),size = 1, replace = FALSE, prob = race_prob[x,]))
    
    ## update alpha_c
    m_educ <- multinom(educ ~ as.factor(sex) + as.factor(race), 
                       data = nc_adj[which(nc_adj$unit == 0),],family=binomial())
    educ_mat <- model.matrix(m_educ)
    alpha_c <- rmvnorm(1,mean = c(t(coef(m_educ))),sigma = vcov(m_educ))
    alpha_c_mat <- matrix(alpha_c,nrow = 2, byrow = T)
    
    # generate replicate educ|sex,race for U=0
    m_educ_rep <- multinom(educ ~ as.factor(sex) + as.factor(race), 
                           data = nc_rep[which(nc_adj$unit == 0),],family=binomial())
    educ_mat_rep <- model.matrix(m_educ_rep)
    educ_prob_rep <- predict(m_educ_rep,type = "prob")
    educ_prob_alpha <- educ_prob_rep
    educ_prob_alpha[,1] <- 1
    for (l in 1:2){
      educ_prob_alpha[,l+1] <- logit_inv(educ_mat_rep %*%(alpha_c_mat[l,]))
      educ_prob_alpha[,l+1] <- educ_prob_alpha[,l+1]/(1-educ_prob_alpha[,l+1])
    }
    educ_prob_alpha_adj <- educ_prob_alpha/rowSums(educ_prob_alpha)
    nc_rep$educ[which(nc_rep$unit == 0)] <- 
      sapply(1:n_unit0,function(x)sample(c(1,2,3),size = 1, replace = FALSE, prob = educ_prob_alpha_adj[x,]))
    
    # generate unit nonresponse for educ
    nc_U1 <- nc_adj[which(nc_adj$unit == 1),]
    educ_mat_U1 <- educ_mat[1:dim(nc_U1)[1],]
    educ_mat_U1[,1] <- 1
    educ_mat_U1[,2] <- nc_U1$sex
    educ_mat_U1[,3] <- ifelse(nc_U1$race == 2,1,0)
    educ_mat_U1[,4] <- ifelse(nc_U1$race == 3,1,0)
    educ_mat_U1[,5] <- ifelse(nc_U1$race == 4,1,0)
    educ_prob <- predict(m_educ,newdata = nc_U1,type = "prob")
    educ_prob_alpha_U1 <- educ_prob
    educ_prob_alpha_U1[,1] <- 1
    for (l in 1:2){
      educ_prob_alpha_U1[,l+1] <- logit_inv(educ_mat_U1%*%(alpha_c_mat[l,]))
      educ_prob_alpha_U1[,l+1] <- educ_prob_alpha_U1[,l+1]/(1-educ_prob_alpha_U1[,l+1])
    }
    educ_prob_alpha_adj_U1 <- educ_prob_alpha_U1/rowSums(educ_prob_alpha_U1)
    nc_adj$educ[which(nc_adj$unit == 1)] <- nc_rep$educ[which(nc_rep$unit == 1)] <- 
      sapply(1:n_unit1,function(x)sample(c(1,2,3),size = 1, replace = FALSE, prob = educ_prob_alpha_adj_U1[x,]))
      
    ## update alpha_a
    m_age <- multinom(age ~ as.factor(sex) + as.factor(race) + as.factor(educ) , data = nc_adj[which(nc_adj$unit == 0),])
    age_mat <- model.matrix(m_age)
    alpha_a <- rmvnorm(1,mean = c(t(coef(m_age))),sigma = vcov(m_age))
    alpha_a_mat <- matrix(alpha_a, nrow = 6, byrow = T)
    
    # generate replicate age|sex, race
    m_age_rep <- multinom(age ~ as.factor(sex) + as.factor(race) + as.factor(educ), 
                          data = nc_rep[which(nc_rep$unit == 0),])
    age_mat_rep <- model.matrix(m_age_rep)
    age_prob_rep <-  predict(m_age_rep,type = "prob")
    age_prob_alpha <- age_prob_rep
    age_prob_alpha[,1] <- 1
    for (l in 1:6){
      age_prob_alpha[,l+1] <- logit_inv(age_mat_rep %*%(alpha_a_mat[l,]))
      # age_prob_alpha[,l+1] <- logit_inv(age_mat_rep %*%(coef(m_age_rep)[l,]))
      age_prob_alpha[,l+1] <- age_prob_alpha[,l+1]/(1-age_prob_alpha[,l+1])
    }
    age_prob_alpha_adj <- age_prob_alpha/rowSums(age_prob_alpha)
    age_char <- sapply(1:n_unit0,function(x)sample(c("(0,29]","(29,39]" , "(39,49]" , "(49,59]" , "(59,69]" , "(69,79]" , "(79,100]"),
                                                   size = 1, replace = FALSE, prob = age_prob_alpha_adj[x,]))
    nc_rep$age[which(nc_rep$unit == 0)] <- age_char
    
    # generate unit nonresponse for age
    # nc_U1 <- nc_adj[which(nc_adj$unit == 1),]
    age_matU1 <- age_mat[1:dim(nc_U1)[1],]
    age_matU1[,2] <- nc_U1$sex
    age_matU1[,3] <- ifelse(nc_U1$race == 2,1,0)
    age_matU1[,4] <- ifelse(nc_U1$race == 3,1,0)
    age_matU1[,5] <- ifelse(nc_U1$race == 4,1,0)
    age_matU1[,6] <- nc_U1$educ
    age_prob <- predict(m_age,newdata = nc_U1,type = "prob")
    age_prob_alpha_U1 <- age_prob
    age_prob_alpha_U1[,1] <- 1
    for (l in 1:6){
      age_prob_alpha_U1[,l+1] <- logit_inv(age_matU1 %*%(alpha_a_mat[l,]))
      # age_prob_alpha_U1[,l+1] <- logit_inv(age_matU1 %*%(coef(m_age)[l,]))
      age_prob_alpha_U1[,l+1] <- age_prob_alpha_U1[,l+1]/(1-age_prob_alpha_U1[,l+1])
    }
    age_prob_alpha_adj_U1 <- age_prob_alpha_U1/rowSums(age_prob_alpha_U1)
    
    age_char <- sapply(1:n_unit1,function(x)sample(c("(0,29]","(29,39]" , "(39,49]" , "(49,59]" , "(59,69]" , "(69,79]" , "(79,100]"),
                                                   size = 1, replace = FALSE, prob = age_prob_alpha_adj_U1[x,]))
    nc_adj$age[which(nc_adj$unit == 1)] <- nc_rep$age[which(nc_rep$unit == 1)] <- age_char
    
    
    ## update alpha_v
    m_vote <- glm(voted ~ as.factor(sex) + age + as.factor(race) + as.factor(educ)
                  + as.factor(sex)*as.factor(race) + as.factor(sex)*age + as.factor(sex)*as.factor(educ),
                  data = nc_adj[which(nc_adj$unit == 0),],family=binomial())
    vote_mat <- model.matrix(m_vote)
    alpha_v_origin <- alpha_v <- rmvnorm(1,mean = c(t(coef(m_vote))),sigma = vcov(m_vote))
    
    # generate vote|sex,race,age
    m_vote_rep <- glm(voted ~ as.factor(sex) + age + as.factor(race) + as.factor(educ) + as.factor(sex)*as.factor(race)
                      + as.factor(sex)*age + as.factor(sex)*as.factor(educ), data = nc_rep[which(nc_rep$unit == 0),],family=binomial())
    vote_mat_rep <- model.matrix(m_vote_rep)
    vote_prob_rep <- logit_inv(vote_mat_rep%*%t(alpha_v))
    nc_rep$voted[which(nc_rep$unit == 0)] <- rbinom(n_unit0,1,prob = vote_prob_rep)
    
    # impute vote for U==1 cases
    T_vote_hat <- rnorm(1,T_vote,Sd_vote)
    total_one_vote <- floor((T_vote_hat - sum((nc_adj$voted*nc_adj$weight)[which(nc_adj$unit == 0)]))/mean(nc_adj$weight[which(nc_adj$unit==0)]))
    vote_mat1 <- vote_mat[1:n_unit1,]
    vote_mat1[,1] <- rep(1,n_unit1)
    vote_mat1[,2] <- nc_adj$sex[which(nc_adj$unit == 1)]
    vote_mat1[,3:13] <- 0
    vote_mat1[,3] <- ifelse(nc_adj$age[which(nc_adj$unit == 1)] == "(29,39]",1,0)
    vote_mat1[,4] <- ifelse(nc_adj$age[which(nc_adj$unit == 1)] == "(39,49]",1,0)
    vote_mat1[,5] <- ifelse(nc_adj$age[which(nc_adj$unit == 1)] == "(49,59]",1,0)
    vote_mat1[,6] <- ifelse(nc_adj$age[which(nc_adj$unit == 1)] == "(59,69]",1,0)
    vote_mat1[,7] <- ifelse(nc_adj$age[which(nc_adj$unit == 1)] == "(69,79]",1,0)
    vote_mat1[,8] <- ifelse(nc_adj$age[which(nc_adj$unit == 1)] == "(79,100]",1,0)
    vote_mat1[,9] <- ifelse(nc_adj$race[which(nc_adj$unit == 1)] == 2,1,0)
    vote_mat1[,10] <- ifelse(nc_adj$race[which(nc_adj$unit == 1)] == 3,1,0)
    vote_mat1[,11] <- ifelse(nc_adj$race[which(nc_adj$unit == 1)] == 4,1,0)
    vote_mat1[,12] <- ifelse(nc_adj$educ[which(nc_adj$unit == 1)] == 2,1,0)
    vote_mat1[,13] <- ifelse(nc_adj$educ[which(nc_adj$unit == 1)] == 3,1,0)
    
    # interactions
    ## sex*race
    vote_mat1[,14] <- ifelse(nc_adj$sex[which(nc_adj$unit == 1)] == 1 & 
                               nc_adj$race[which(nc_adj$unit == 1)] == 2,1,0)
    vote_mat1[,15] <- ifelse(nc_adj$sex[which(nc_adj$unit == 1)] == 1 & 
                               nc_adj$race[which(nc_adj$unit == 1)] == 3,1,0)
    vote_mat1[,16] <- ifelse(nc_adj$sex[which(nc_adj$unit == 1)] == 1 & 
                               nc_adj$race[which(nc_adj$unit == 1)] == 4,1,0)
    ## sex*age
    vote_mat1[,17] <- ifelse(nc_adj$sex[which(nc_adj$unit == 1)] == 1 & 
                               nc_adj$age[which(nc_adj$unit == 1)] == "(29,39]",1,0)
    vote_mat1[,18] <- ifelse(nc_adj$sex[which(nc_adj$unit == 1)] == 1 & 
                               nc_adj$age[which(nc_adj$unit == 1)] == "(39,49]",1,0)
    vote_mat1[,19] <- ifelse(nc_adj$sex[which(nc_adj$unit == 1)] == 1 & 
                               nc_adj$age[which(nc_adj$unit == 1)] == "(49,59]",1,0)
    vote_mat1[,20] <- ifelse(nc_adj$sex[which(nc_adj$unit == 1)] == 1 & 
                               nc_adj$age[which(nc_adj$unit == 1)] == "(59,69]",1,0)
    vote_mat1[,21] <- ifelse(nc_adj$sex[which(nc_adj$unit == 1)] == 1 & 
                               nc_adj$age[which(nc_adj$unit == 1)] == "(69,79]",1,0)
    vote_mat1[,22] <- ifelse(nc_adj$sex[which(nc_adj$unit == 1)] == 1 & 
                               nc_adj$age[which(nc_adj$unit == 1)] == "(79,100]",1,0)
    # sex*educ
    vote_mat1[,23] <- ifelse(nc_adj$sex[which(nc_adj$unit == 1)] == 1 & 
                               nc_adj$educ[which(nc_adj$unit == 1)] == 2,1,0)
    vote_mat1[,24] <- ifelse(nc_adj$sex[which(nc_adj$unit == 1)] == 1 & 
                               nc_adj$educ[which(nc_adj$unit == 1)] == 3,1,0)
    
    # coef in front of vote
    u_adj <- logit(total_one_vote/n_unit1) - mean(vote_mat1%*%t(alpha_v))
    alpha_v[1] <- alpha_v[1] + u_adj
    vote_prob <- logit_inv(vote_mat1%*%t(alpha_v))
    nc_adj$voted[which(nc_adj$unit == 1)] <- nc_rep$voted[which(nc_rep$unit == 1)] <- rbinom(n_unit1,1,prob = vote_prob)
    
    
    ### update gamma's
    # update gamma_e
    r_eth <- glm(Re~as.factor(sex) + as.factor(educ) + as.factor(age) + as.factor(voted), data = nc_adj[which(nc_adj$unit == 0),],family=binomial())
    ethR_mat <- model.matrix(r_eth)
    gamma_e <- rmvnorm(1,mean = c(t(coef(r_eth))),sigma = vcov(r_eth))
    
    # generate replicate of Re|sex,age,voted
    r_eth_rep <- glm(Re~as.factor(sex) + as.factor(educ) + as.factor(age) + as.factor(voted), data = nc_rep[which(nc_rep$unit == 0),],family=binomial())
    ethR_mat_rep <- model.matrix(r_eth_rep)
    # prob_Re_temp <- ifelse(ethR_mat_rep%*%t(gamma_e) > 709,1,ethR_mat_rep%*%t(gamma_e))
    prob_Re_temp <- ethR_mat_rep%*%t(gamma_e)
    prob_Re <- ifelse(prob_Re_temp>709,1,logit_inv(prob_Re_temp))
    nc_rep$Re[which(nc_rep$unit == 0)] <- rbinom(n_unit0,1,prob = prob_Re)
    
    # update gamma_c
    r_educ <- glm(Rc~as.factor(sex) + as.factor(race) + as.factor(age) +  as.factor(voted),
                  data = nc_adj[which(nc_adj$unit == 0),],family=binomial())
    educR_mat <- model.matrix(r_educ)
    gamma_c <- rmvnorm(1,mean = c(t(coef(r_educ))),sigma = vcov(r_educ))
    
    # generate replicate of Rc|sex, race
    r_educ_rep <- glm(Rc~as.factor(sex) + as.factor(race) + as.factor(age) +  as.factor(voted), data = nc_rep[which(nc_rep$unit == 0),],family=binomial())
    educR_mat_rep <- model.matrix(r_educ_rep)
    prob_Rc_temp <- educR_mat_rep%*%t(gamma_c)
    prob_Rc <- ifelse(prob_Rc_temp > 709,1,logit_inv(prob_Rc_temp))
    nc_rep$Rc[which(nc_rep$unit == 0)] <- rbinom(n_unit0,1,prob = prob_Rc)
    
    # update gamma_a
    # r_age <- glm(Ra~as.factor(sex) + as.factor(race) + as.factor(voted), data = nc_adj[which(nc_adj$unit == 0),],family=binomial())
    r_age <- glm(Ra ~ as.factor(sex) + as.factor(race) + as.factor(educ), data = nc_adj[which(nc_adj$unit == 0),],family=binomial())
    ageR_mat <- model.matrix(r_age)
    gamma_a <- rmvnorm(1,mean = c(t(coef(r_age))),sigma = vcov(r_age))
    
    # generate replicate of Ra|sex, race
    r_age_rep <- glm(Ra ~ as.factor(sex) + as.factor(race) + as.factor(educ), data = nc_rep[which(nc_rep$unit == 0),],family=binomial())
    ageR_mat_rep <- model.matrix(r_age_rep)
    prob_Ra <- logit_inv(ageR_mat_rep%*%t(gamma_a))
    nc_rep$Ra[which(nc_rep$unit == 0)] <- rbinom(n_unit0,1,prob = prob_Ra)
    
    # update gamma_v
    r_vote <- glm(Rv~as.factor(sex) + as.factor(race) + as.factor(educ) + as.factor(age),
                  data = nc_adj[which(nc_adj$unit == 0 & nc_adj$Ra == 0),],family=binomial())
    voteR_mat <- model.matrix(r_vote)
    gamma_v <- rmvnorm(1,mean = c(t(coef(r_vote))),sigma = vcov(r_vote))
    
    # generate replicate of Rv|sex, age,race when Ra == 0
    r_vote_rep <- glm(Rv~as.factor(sex) + as.factor(race) + as.factor(educ) + as.factor(age),
                      data = nc_rep[which(nc_rep$unit == 0 & nc_rep$Ra == 0),],family=binomial())
    voteR_mat_rep <- model.matrix(r_vote_rep)
    my_prob_Rv <- logit_inv(voteR_mat_rep%*%t(gamma_v))
    nc_rep$Rv[which(nc_rep$unit == 0 & nc_rep$Ra == 0)] <- rbinom(sum(nc_rep$unit == 0 & nc_rep$Ra == 0),1,prob = my_prob_Rv)
    # introduce dependence
    nc_rep$Rv[which(nc_rep$Ra == 1)] <- 1
    
    ### sample item non-response obs
    ## race
    n_mis_e <- sum(nc_adj[which(nc_adj$unit==0),]$Re == 1)
    pr_eth_miss <- matrix(0,ncol=4,nrow=n_mis_e)
    colnames(pr_eth_miss) <- c("1","2","3","4")
    
    ## subset of missing race obs
    sub_eth <- nc_adj[which(nc_adj$Re == 1 & nc_adj$unit == 0),]
    alpha_e_mat <- matrix(alpha_e, nrow = 3, byrow = T)
    # p(Ei|alpha_e)
    pi_eth <- predict(m_eth,newdata = sub_eth, type = "probs") 
    eth_mis_mat <- eth_mat[1:dim(sub_eth)[1],]
    eth_mis_mat[,2] <- sub_eth$sex
    eth_mis_prob_alpha <- pi_eth
    eth_mis_prob_alpha[,1] <- 1
    for (l in 1:3){
      eth_mis_prob_alpha[,l+1] <- logit_inv(eth_mis_mat %*%(alpha_e_mat[l,]))
      eth_mis_prob_alpha[,l+1] <- eth_mis_prob_alpha[,l+1]/(1-eth_mis_prob_alpha[,l+1])
    }
    pi_eth <- eth_mis_prob_alpha/rowSums(eth_mis_prob_alpha)
    
    
    # p(Ai|Ei)
    age_prob <- pi_eth
    for (j in 1:dim(pi_eth)[1]){
      single_mat <- age_mat[1,]
      single_mat[2] <- sub_eth$sex[j]
      single_mat[6] <- ifelse(sub_eth$educ[j] == 2,1,0)
      single_mat[7] <- ifelse(sub_eth$educ[j] == 3,1,0)
      single_mat[3] <- single_mat[4] <- single_mat[5] <- 0
      for (k in 1:4){
        if (k > 1){
          single_mat[k+1] <- 1
        }
        single_prob <- c(1,logit_inv(single_mat%*%t(alpha_a_mat)))
        single_prob[2:7] <- single_prob[2:7]/(1-single_prob[2:7])
        single_prob_adj <- single_prob/sum(single_prob)
        names(single_prob_adj) <- c("(0,29]","(29,39]","(39,49]","(49,59]","(59,69]","(69,79]", "(79,100]")
        age_prob[j,k] <- single_prob_adj[sub_eth$age[j]]
      }
    }
    
    # p(Ci|Ei)
    educ_prob <- pi_eth
    for (j in 1:dim(pi_eth)[1]){
      single_mat <- educ_mat[1,]
      single_mat[2] <- sub_eth$sex[j]
      single_mat[3] <- single_mat[4] <- single_mat[5] <- 0
      for (k in 1:4){
        if (k > 1){
          single_mat[k+1] <- 1
        }
        single_prob <- c(1,logit_inv(single_mat%*%t(alpha_c_mat)))
        single_prob[2:3] <- single_prob[2:3]/(1-single_prob[2:3])
        single_prob_adj <- single_prob/sum(single_prob)
        educ_prob[j,k] <- single_prob_adj[sub_eth$educ[j]]
      }
    }
    
    # p(Vi|Ei)
    vote_prob <- pi_eth
    vote_mat_mis <- vote_mat[1:dim(sub_eth)[1],]
    vote_mat_mis[,1] <- 1
    for (k in 1:4){
      vote_mat_mis[,2] <- sub_eth$sex
      vote_mat_mis[,3] <- ifelse(sub_eth$age == "(29,39]",1,0)
      vote_mat_mis[,4] <- ifelse(sub_eth$age == "(39,49]",1,0)
      vote_mat_mis[,5] <- ifelse(sub_eth$age == "(49,59]",1,0)
      vote_mat_mis[,6] <- ifelse(sub_eth$age == "(59,69]",1,0)
      vote_mat_mis[,7] <- ifelse(sub_eth$age == "(69,79]",1,0)
      vote_mat_mis[,8] <- ifelse(sub_eth$age == "(79,100]",1,0)
      vote_mat_mis[,9] <- vote_mat_mis[,10] <- vote_mat_mis[,11] <- 0
      vote_mat_mis[,14] <- vote_mat_mis[,15] <- vote_mat_mis[,16] <- 0
      vote_mat_mis[,12] <- ifelse(sub_eth$educ == 2,1,0)
      vote_mat_mis[,13] <- ifelse(sub_eth$educ == 3,1,0)
      if (k==2){
        vote_mat_mis[,9] <- 1
        vote_mat_mis[,14] <- vote_mat_mis[,9]*sub_eth$sex
      }
      if(k==3){
        vote_mat_mis[,10] <- 1
        vote_mat_mis[,15] <- vote_mat_mis[,10]*sub_eth$sex
      }
      if(k==4){
        vote_mat_mis[,11] <- 1
        vote_mat_mis[,16] <- vote_mat_mis[,11]*sub_eth$sex
      }
      
      # interactions
      ## sex*age
      vote_mat_mis[,17] <- ifelse(sub_eth$sex == 1 & sub_eth$age == "(29,39]",1,0)
      vote_mat_mis[,18] <- ifelse(sub_eth$sex == 1 & sub_eth$age == "(39,49]",1,0)
      vote_mat_mis[,19] <- ifelse(sub_eth$sex == 1 & sub_eth$age == "(49,59]",1,0)
      vote_mat_mis[,20] <- ifelse(sub_eth$sex == 1 & sub_eth$age == "(59,69]",1,0)
      vote_mat_mis[,21] <- ifelse(sub_eth$sex == 1 & sub_eth$age == "(69,79]",1,0)
      vote_mat_mis[,22] <- ifelse(sub_eth$sex == 1 & sub_eth$age == "(79,100]",1,0)
      vote_mat_mis[,23] <- ifelse(sub_eth$sex == 1 & sub_eth$educ == 2,1,0)
      vote_mat_mis[,24] <- ifelse(sub_eth$sex == 1 & sub_eth$educ == 3,1,0)
      
      vote_pred <- logit_inv(vote_mat_mis%*%t(alpha_v_origin))
      
      vote_prob[,k] <- ifelse(sub_eth$voted == 1, vote_pred,1 - vote_pred)
    }
    
    # p(R_C|E)
    RC_prob <- pi_eth
    RC_mat_mis <- educR_mat[1:dim(sub_eth)[1],]
    RC_mat_mis[1,] <- 1
    for (k in 1:4){
      RC_mat_mis[,2] <- sub_eth$sex
      RC_mat_mis[,3] <- RC_mat_mis[,4] <- RC_mat_mis[,5] <- 0
      if (k==2){
        RC_mat_mis[,3] <- 1
      }
      if(k==3){
        RC_mat_mis[,4] <- 1
      }
      if(k==4){
        RC_mat_mis[,5] <- 1
      }
      RC_mat_mis[,6] <- ifelse(sub_eth$age == "(29,39]",1,0)
      RC_mat_mis[,7] <- ifelse(sub_eth$age == "(39,49]",1,0)
      RC_mat_mis[,8] <- ifelse(sub_eth$age == "(49,59]",1,0)
      RC_mat_mis[,9] <- ifelse(sub_eth$age == "(59,69]",1,0)
      RC_mat_mis[,10] <- ifelse(sub_eth$age == "(69,79]",1,0)
      RC_mat_mis[,11] <- ifelse(sub_eth$age == "(79,100]",1,0)
      RC_mat_mis[,12] <- sub_eth$educ
      RC_pred <- logit_inv(RC_mat_mis%*%t(gamma_c))
      RC_prob[,k] <- ifelse(sub_eth$Rc == 1, RC_pred,1 - RC_pred)
    }
    
    # p(R_A|E)
    RA_prob <- pi_eth
    RA_mat_mis <- ageR_mat[1:dim(sub_eth)[1],]
    RA_mat_mis[,1] <- 1
    for (k in 1:4){
      RA_mat_mis[,2] <- sub_eth$sex
      RA_mat_mis[,6] <- ifelse(sub_eth$educ==2,1,0)
      RA_mat_mis[,7] <- ifelse(sub_eth$educ==3,1,0)
      RA_mat_mis[,3] <- RA_mat_mis[,4] <- RA_mat_mis[,5] <- 0
      if (k==2){
        RA_mat_mis[,3] <- 1
      }
      if(k==3){
        RA_mat_mis[,4] <- 1
      }
      if(k==4){
        RA_mat_mis[,5] <- 1
      }
      RA_pred <- logit_inv(RA_mat_mis%*%t(gamma_a))
      RA_prob[,k] <- ifelse(sub_eth$Ra == 1, RA_pred,1 - RA_pred)
    }
    
    # p(R_V|E)
    RV_prob <- pi_eth
    RV_mat_mis <- voteR_mat[1:dim(sub_eth)[1],]
    RV_mat_mis[,1] <- 1
    for (k in 1:4){
      RV_mat_mis[,2] <- sub_eth$sex
      RV_mat_mis[,6] <- ifelse(sub_eth$educ==2,1,0)
      RV_mat_mis[,7] <- ifelse(sub_eth$educ==3,1,0)
      RV_mat_mis[,8] <- ifelse(sub_eth$age == "(29,39]",1,0)
      RV_mat_mis[,9] <- ifelse(sub_eth$age == "(39,49]",1,0)
      RV_mat_mis[,10] <- ifelse(sub_eth$age == "(49,59]",1,0)
      RV_mat_mis[,11] <- ifelse(sub_eth$age == "(59,69]",1,0)
      RV_mat_mis[,12] <- ifelse(sub_eth$age == "(69,79]",1,0)
      RV_mat_mis[,13] <- ifelse(sub_eth$age == "(79,100]",1,0)
      RV_mat_mis[,3] <- vote_mat_mis[,4] <- vote_mat_mis[,5] <- 0
      if (k==2){
        RV_mat_mis[,3] <- 1
      }
      if(k==3){
        RV_mat_mis[,4] <- 1
      }
      if(k==4){
        RV_mat_mis[,5] <- 1
      }
      RV_pred <- logit_inv(RV_mat_mis%*%t(gamma_v))
      RV_prob[,k] <- ifelse(sub_eth$Rv == 1, RV_pred,1 - RV_pred)
    }
    
    # put it all together
    pr_eth_miss <- exp(log(pi_eth)+log(age_prob)+log(vote_prob) + log(educ_prob)
                       + log(RC_prob)+log(RA_prob)+log(RV_prob))
    pr_eth_miss_prob <- ifelse(is.na(pr_eth_miss/rowSums(pr_eth_miss)),0,pr_eth_miss/rowSums(pr_eth_miss))
    pr_eth_miss_prob[which(pr_eth_miss_prob == 0)] <- 1e-10
    nc_adj[which(nc_adj$Re == 1 & nc_adj$unit == 0),]$race <- 
      sapply(1:dim(pi_eth)[1],function(x) sample(c(1,2,3,4),1,prob = pr_eth_miss_prob[x,]))
    
    ## Educ
    n_mis_c <- sum(nc_adj[which(nc_adj$unit==0),]$Rc == 1)
    pr_educ_miss <- matrix(0,ncol=3,nrow=n_mis_c)
    ## subset of missing educ obs
    sub_educ <- nc_adj[which(nc_adj$Rc == 1 & nc_adj$unit == 0),]
    
    # p(Ci|alpha_c)
    pi_educ <- predict(m_educ,newdata = sub_educ, type = "probs") 
    educ_mis_mat <- educ_mat[1:dim(sub_educ)[1],]
    educ_mis_mat[,2] <- sub_educ$sex
    educ_mis_mat[,3] <- ifelse(sub_educ$race == 2,1,0)
    educ_mis_mat[,4] <- ifelse(sub_educ$race == 3,1,0)
    educ_mis_mat[,5] <- ifelse(sub_educ$race == 4,1,0)
    educ_mis_prob_alpha <- pi_educ
    educ_mis_prob_alpha[,1] <- 1
    for (l in 1:2){
      educ_mis_prob_alpha[,l+1] <- logit_inv(educ_mis_mat %*%(alpha_c_mat[l,]))
      educ_mis_prob_alpha[,l+1] <- educ_mis_prob_alpha[,l+1]/(1-educ_mis_prob_alpha[,l+1])
    }
    pi_educ <- educ_mis_prob_alpha/rowSums(educ_mis_prob_alpha)
    
    # p(Ai|Ci)
    ageC_prob <- pi_educ
    for (j in 1:dim(pi_educ)[1]){
      single_mat <- age_mat[1,]
      single_mat[2] <- sub_educ$sex[j]
      single_mat[3] <- ifelse(sub_educ$race[j] == 2,1,0)
      single_mat[4] <- ifelse(sub_educ$race[j] == 3,1,0)
      single_mat[5] <- ifelse(sub_educ$race[j] == 4,1,0)
      single_mat[6] <- single_mat[7] <- 0
      for (k in 1:3){
        if (k > 1){
          single_mat[k+4] <- 1
        }
        single_prob <- c(1,logit_inv(single_mat%*%t(alpha_a_mat)))
        single_prob[2:7] <- single_prob[2:7]/(1-single_prob[2:7])
        single_prob_adj <- single_prob/sum(single_prob)
        names(single_prob_adj) <- c("(0,29]","(29,39]","(39,49]","(49,59]","(59,69]","(69,79]", "(79,100]")
        ageC_prob[j,k] <- single_prob_adj[sub_educ$age[j]]
      }
    } 
    
    # p(Vi|Ci)
    voteC_prob <- pi_educ
    voteC_mat_mis <- vote_mat[1:dim(sub_educ)[1],]
    voteC_mat_mis[,1] <- 1
    for (k in 1:3){
      voteC_mat_mis[,2] <- sub_educ$sex
      voteC_mat_mis[,3] <- ifelse(sub_educ$age == "(29,39]",1,0)
      voteC_mat_mis[,4] <- ifelse(sub_educ$age == "(39,49]",1,0)
      voteC_mat_mis[,5] <- ifelse(sub_educ$age == "(49,59]",1,0)
      voteC_mat_mis[,6] <- ifelse(sub_educ$age == "(59,69]",1,0)
      voteC_mat_mis[,7] <- ifelse(sub_educ$age == "(69,79]",1,0)
      voteC_mat_mis[,8] <- ifelse(sub_educ$age == "(79,100]",1,0)
      voteC_mat_mis[,9] <- ifelse(sub_educ$race == 2,1,0)
      voteC_mat_mis[,10] <- ifelse(sub_educ$race == 3,1,0)
      voteC_mat_mis[,11] <- ifelse(sub_educ$race == 4,1,0)
      voteC_mat_mis[,12] <- voteC_mat_mis[,23] <- voteC_mat_mis[,13] <- voteC_mat_mis[,24] <- 0
      if (k == 2){
        voteC_mat_mis[,12] <- 1
        voteC_mat_mis[,23] <- voteC_mat_mis[,12]*sub_educ$sex
      }
      if (k == 3){
        voteC_mat_mis[,13] <- 1
        voteC_mat_mis[,24] <- voteC_mat_mis[,13]*sub_educ$sex
      }
      
      # interactions
      ## sex*age
      voteC_mat_mis[,14] <- ifelse(sub_educ$sex == 1 & sub_educ$race == 2,1,0)
      voteC_mat_mis[,15] <- ifelse(sub_educ$sex == 1 & sub_educ$race == 3,1,0)
      voteC_mat_mis[,16] <- ifelse(sub_educ$sex == 1 & sub_educ$race == 4,1,0)
      voteC_mat_mis[,17] <- ifelse(sub_educ$sex == 1 & sub_educ$age == "(29,39]",1,0)
      voteC_mat_mis[,18] <- ifelse(sub_educ$sex == 1 & sub_educ$age == "(39,49]",1,0)
      voteC_mat_mis[,19] <- ifelse(sub_educ$sex == 1 & sub_educ$age == "(49,59]",1,0)
      voteC_mat_mis[,20] <- ifelse(sub_educ$sex == 1 & sub_educ$age == "(59,69]",1,0)
      voteC_mat_mis[,21] <- ifelse(sub_educ$sex == 1 & sub_educ$age == "(69,79]",1,0)
      voteC_mat_mis[,22] <- ifelse(sub_educ$sex == 1 & sub_educ$age == "(79,100]",1,0)
      
      voteC_pred <- logit_inv(voteC_mat_mis%*%t(alpha_v_origin))
      
      voteC_prob[,k] <- ifelse(sub_educ$voted == 1, voteC_pred,1 - voteC_pred)
    }
    
    # p(Re|Ci)
    REC_prob <- pi_educ
    r_eth_mis_mat <- ethR_mat[1:dim(pi_educ)[1],]
    for (k in 1:3){
      r_eth_mis_mat[,1] <- 1
      r_eth_mis_mat[,2] <- sub_educ$sex
      r_eth_mis_mat[,5] <- ifelse(sub_educ$age == "(29,39]",1,0)
      r_eth_mis_mat[,6] <- ifelse(sub_educ$age == "(39,49]",1,0)
      r_eth_mis_mat[,7] <- ifelse(sub_educ$age == "(49,59]",1,0)
      r_eth_mis_mat[,8] <- ifelse(sub_educ$age == "(59,69]",1,0)
      r_eth_mis_mat[,9] <- ifelse(sub_educ$age == "(69,79]",1,0)
      r_eth_mis_mat[,10] <- ifelse(sub_educ$age == "(79,100]",1,0)
      r_eth_mis_mat[,11] <- sub_educ$voted
      r_eth_mis_mat[,3] <- r_eth_mis_mat[,4] <- 0
      if (k > 1){
        r_eth_mis_mat[,k+1] <- 1
      }
      REC_pred_temp <- r_eth_mis_mat%*%t(gamma_e)
      REC_pred <- ifelse(REC_pred_temp > 709,1,logit_inv(REC_pred_temp))
      REC_prob[,k] <- ifelse(sub_educ$Re==1, REC_pred, 1-REC_pred)
    }
    REC_prob[which(REC_prob == 0)] <- 1e-10
    
    # p(Ra|Ci)
    RAC_prob <- pi_educ
    r_age_mis_mat <- ageR_mat[1:dim(pi_educ)[1],]
    for (k in 1:3){
      r_age_mis_mat[,1] <- 1
      r_age_mis_mat[,2] <- sub_educ$sex
      r_age_mis_mat[,3] <- ifelse(sub_educ$race == 2,1,0)
      r_age_mis_mat[,4] <- ifelse(sub_educ$race == 3,1,0)
      r_age_mis_mat[,5] <- ifelse(sub_educ$race == 4,1,0)
      
      r_age_mis_mat[,6] <- r_age_mis_mat[,7] <- 0
      if (k > 1){
        r_age_mis_mat[,k+4] <- 1
      }
      RAC_pred <- logit_inv(r_age_mis_mat%*%t(gamma_a))
      RAC_prob[,k] <- ifelse(sub_educ$Ra==1, RAC_pred, 1-RAC_pred)
    }
    
    # p(Rv|Ci)
    RVC_prob <- pi_educ
    r_vote_mis_mat <- voteR_mat[1:dim(pi_educ)[1],]
    for (k in 1:3){
      r_vote_mis_mat[,1] <- 1
      r_vote_mis_mat[,2] <- sub_educ$sex
      r_vote_mis_mat[,3] <- ifelse(sub_educ$race == 2,1,0)
      r_vote_mis_mat[,4] <- ifelse(sub_educ$race == 3,1,0)
      r_vote_mis_mat[,5] <- ifelse(sub_educ$race == 4,1,0)
      r_vote_mis_mat[,8] <- ifelse(sub_educ$age == "(29,39]",1,0)
      r_vote_mis_mat[,9] <- ifelse(sub_educ$age == "(39,49]",1,0)
      r_vote_mis_mat[,10] <- ifelse(sub_educ$age == "(49,59]",1,0)
      r_vote_mis_mat[,11] <- ifelse(sub_educ$age == "(59,69]",1,0)
      r_vote_mis_mat[,12] <- ifelse(sub_educ$age == "(69,79]",1,0)
      r_vote_mis_mat[,13] <- ifelse(sub_educ$age == "(79,100]",1,0)
      
      r_vote_mis_mat[,6] <- r_vote_mis_mat[,7] <- 0
      if (k > 1){
        r_vote_mis_mat[,k+4] <- 1
      }
      RVC_pred <- logit_inv(r_vote_mis_mat%*%t(gamma_v))
      RVC_prob[,k] <- ifelse(sub_educ$Rv==1, RVC_pred, 1-RVC_pred)
    }
    
    # put it all together
    pr_educ_miss <- exp(log(pi_educ) + log(ageC_prob) + log(voteC_prob) + log(REC_prob) 
                        + log(RAC_prob) + log(RVC_prob))
    # pr_educ_miss <- pi_educ*ageC_prob*voteC_prob*REC_prob*RAC_prob*RVC_prob
    pr_educ_miss_prob <- ifelse(is.na(pr_educ_miss/rowSums(pr_educ_miss)),0,pr_educ_miss/rowSums(pr_educ_miss))
    nc_adj[which(nc_adj$Rc == 1 & nc_adj$unit == 0),]$educ <- 
      sapply(1:dim(pi_educ)[1],function(x) sample(c(1,2,3),1,prob = pr_educ_miss_prob[x,]))
    
    ## Age
    n_mis_a <- sum(nc_adj[which(nc_adj$unit==0),]$Ra == 1)
    pr_age_miss <- matrix(0,ncol=7,nrow=n_mis_a)
    colnames(pr_age_miss) <- c(row.names(table(nc_adj$age)))
    ## subset of missing age obs
    sub_age <- nc_adj[which(nc_adj$Ra == 1 & nc_adj$unit == 0),]
    
    # p(Ai|alpha_a)
    pi_age <- predict(m_age,newdata = sub_age, type = "probs") 
    age_mis_mat <- age_mat[1:dim(sub_age)[1],]
    age_mis_mat[,2] <- sub_age$sex
    ###### added 03/30, may need to update CPS_latest_debug #####
    age_mis_mat[,3] <- ifelse(sub_age$race == 2,1,0)
    age_mis_mat[,4] <- ifelse(sub_age$race == 3,1,0)
    age_mis_mat[,5] <- ifelse(sub_age$race == 4,1,0)
    age_mis_mat[,6] <- ifelse(sub_age$educ == 2,1,0)
    age_mis_mat[,7] <- ifelse(sub_age$educ == 3,1,0)
    #############################################################
    age_mis_prob_alpha <- pi_age
    age_mis_prob_alpha[,1] <- 1
    for (l in 1:6){
      age_mis_prob_alpha[,l+1] <- logit_inv(age_mis_mat %*%(alpha_a_mat[l,]))
      age_mis_prob_alpha[,l+1] <- age_mis_prob_alpha[,l+1]/(1-age_mis_prob_alpha[,l+1])
    }
    pi_age <- age_mis_prob_alpha/rowSums(age_mis_prob_alpha)
    
    # p(Vi|Ai)
    voteA_prob <- pi_age
    voteA_mat_mis <- vote_mat[1:dim(sub_age)[1],]
    voteA_mat_mis[1,] <- 1
    for (k in 1:7){
      voteA_mat_mis[,2] <- sub_age$sex
      voteA_mat_mis[,9] <- ifelse(sub_age$race == 2,1,0)
      voteA_mat_mis[,10] <- ifelse(sub_age$race == 3,1,0)
      voteA_mat_mis[,11] <- ifelse(sub_age$race == 4,1,0)
      voteA_mat_mis[,12] <- ifelse(sub_age$educ == 2,1,0)
      voteA_mat_mis[,13] <- ifelse(sub_age$educ == 3,1,0)
      voteA_mat_mis[,3:8] <- voteA_mat_mis[,17:22] <- 0
      if (k==2){
        voteA_mat_mis[,3] <- 1
        voteA_mat_mis[,17] <- voteA_mat_mis[,3]*sub_age$sex
      }
      if(k==3){
        voteA_mat_mis[,4] <- 1
        voteA_mat_mis[,18] <- voteA_mat_mis[,4]*sub_age$sex
      }
      if(k==4){
        voteA_mat_mis[,5] <- 1
        voteA_mat_mis[,19] <- voteA_mat_mis[,5]*sub_age$sex
      }
      if(k==5){
        voteA_mat_mis[,6] <- 1
        voteA_mat_mis[,20] <- voteA_mat_mis[,6]*sub_age$sex
      }
      if(k==6){
        voteA_mat_mis[,6] <- 1
        voteA_mat_mis[,21] <- voteA_mat_mis[,6]*sub_age$sex
      }
      if(k==7){
        voteA_mat_mis[,7] <- 1
        voteA_mat_mis[,22] <- voteA_mat_mis[,7]*sub_age$sex
      }
      
      # interactions
      ## sex*race
      voteA_mat_mis[,14] <- ifelse(sub_age$sex == 1 & sub_age$race == 2,1,0)
      voteA_mat_mis[,15] <- ifelse(sub_age$sex == 1 & sub_age$race == 3,1,0)
      voteA_mat_mis[,16] <- ifelse(sub_age$sex == 1 & sub_age$race == 4,1,0)
      
      ## sex*educ
      voteA_mat_mis[,23] <- ifelse(sub_age$sex == 1 & sub_age$educ == 2,1,0)
      voteA_mat_mis[,24] <- ifelse(sub_age$sex == 1 & sub_age$educ == 3,1,0)
      
      voteA_pred <- logit_inv(voteA_mat_mis%*%t(alpha_v_origin))
      
      voteA_prob[,k] <- ifelse(sub_age$voted == 1, voteA_pred,1 - voteA_pred)
    }
    
    # p(R_E|A)
    REA_prob <- pi_age
    REA_mat_mis <- ethR_mat[1:dim(sub_age)[1],]
    REA_mat_mis[,1] <- 1
    for (k in 1:7){
      REA_mat_mis[,2] <- sub_age$sex
      REA_mat_mis[,3] <- ifelse(sub_age$educ == 2,1,0)
      REA_mat_mis[,4] <- ifelse(sub_age$educ == 3,1,0)
      REA_mat_mis[,5:10] <- 0
      if (k==2){
        REA_mat_mis[,5] <-1
      }
      if (k==3){
        REA_mat_mis[,6] <-1
      }
      if (k==4){
        REA_mat_mis[,7] <-1
      }
      if (k==5){
        REA_mat_mis[,8] <-1
      }
      if (k==6){
        REA_mat_mis[,9] <-1
      }
      if (k==7){
        REA_mat_mis[,10] <-1
      }
      REA_mat_mis[,11] <- sub_age$voted
      # prob_Re_temp <- ethR_mat_rep%*%t(gamma_e)
      # prob_Re <- ifelse(prob_Re_temp>709,1,logit_inv(prob_Re_temp))
      REA_pred_temp <- REA_mat_mis%*%t(gamma_e)
      REA_pred <- ifelse(REA_pred_temp > 709,1,logit_inv(REA_pred_temp))
      REA_prob[,k] <- ifelse(sub_age$Re == 1, REA_pred,1 - REA_pred)
    }
    # need a temp
    
    # p(R_V|A)
    RVA_prob <- pi_age
    RVA_mat_mis <- voteR_mat[1:dim(sub_age)[1],]
    RVA_mat_mis[,1] <- 1
    for (k in 1:7){
      RVA_mat_mis[,2] <- sub_age$sex
      RVA_mat_mis[,6] <- ifelse(sub_age$educ == 2,1,0)
      RVA_mat_mis[,7] <- ifelse(sub_age$educ == 3,1,0)
      RVA_mat_mis[,8:13] <- 0
      if (k==2){
        RVA_mat_mis[,8] <-1
      }
      if (k==3){
        RVA_mat_mis[,9] <-1
      }
      if (k==4){
        RVA_mat_mis[,10] <-1
      }
      if (k==5){
        RVA_mat_mis[,11] <-1
      }
      if (k==6){
        RVA_mat_mis[,12] <-1
      }
      if (k==7){
        RVA_mat_mis[,13] <-1
      }
      RVA_mat_mis[,3] <- ifelse(sub_age$race == 2,1,0)
      RVA_mat_mis[,4] <- ifelse(sub_age$race == 3,1,0)
      RVA_mat_mis[,5] <- ifelse(sub_age$race == 4,1,0)
      RVA_pred <- logit_inv(RVA_mat_mis%*%t(gamma_v))
      RVA_prob[,k] <- ifelse(sub_age$Re == 1,RVA_pred,1 - RVA_pred)
    }
    
    # p(R_C|A)
    RCA_prob <- pi_age
    RCA_mat_mis <- educR_mat[1:dim(sub_age)[1],]
    RCA_mat_mis[,1] <- 1
    for (k in 1:7){
      RCA_mat_mis[,2] <- sub_age$sex
      RCA_mat_mis[,12] <- sub_age$voted
      RCA_mat_mis[,6:11] <- 0
      if (k==2){
        RCA_mat_mis[,6] <-1
      }
      if (k==3){
        RCA_mat_mis[,7] <-1
      }
      if (k==4){
        RCA_mat_mis[,8] <-1
      }
      if (k==5){
        RCA_mat_mis[,9] <-1
      }
      if (k==6){
        RCA_mat_mis[,10] <-1
      }
      if (k==7){
        RCA_mat_mis[,11] <-1
      }
      RCA_mat_mis[,3] <- ifelse(sub_age$race == 2,1,0)
      RCA_mat_mis[,4] <- ifelse(sub_age$race == 3,1,0)
      RCA_mat_mis[,5] <- ifelse(sub_age$race == 4,1,0)
      RCA_pred_temp <- RCA_mat_mis%*%t(gamma_c)
      RCA_pred <- ifelse(RCA_pred_temp > 709,1,logit_inv(RCA_pred_temp))
      # RCA_pred <- logit_inv(RCA_mat_mis%*%t(gamma_c))
      RCA_prob[,k] <- ifelse(sub_age$Re == 1,RCA_pred,1 - RCA_pred)
    }
    
    # put it all together
    pr_age_miss <- pi_age*voteA_prob*REA_prob*RVA_prob*RCA_prob
    pr_age_miss_prob <- pr_age_miss/rowSums(pr_age_miss)
    nc_adj[which(nc_adj$Ra == 1 & nc_adj$unit == 0),]$age <- 
      sapply(1:dim(pi_age)[1],function(x)sample(colnames(pr_age_miss),1,prob = pr_age_miss_prob[x,]))
  }
  
  ## vote
  n_mis_v <- sum(nc_adj[which(nc_adj$unit==0),]$Rv == 1)
  pr_vote_miss <- matrix(0,ncol=2,nrow=n_mis_v)
  colnames(pr_vote_miss) <- c("1","2")
  ## subset of missing vote obs
  sub_vote <- nc_adj[which(nc_adj$Rv == 1 & nc_adj$unit == 0),]
  
  # p(Vi|alpha_v)
  m_vote_mis <- glm(voted ~ as.factor(sex) + age + as.factor(race) + as.factor(educ) +
                      as.factor(sex)*as.factor(race)
                    + as.factor(sex)*age + as.factor(sex)*as.factor(educ), 
                    family = binomial(), data = sub_vote) # why NA in sub_vote????
  vote_mat_mis <- model.matrix(m_vote_mis)
  pi_vote <- pr_vote_miss
  pi_vote[,1] <- 1-logit_inv(vote_mat_mis%*%t(alpha_v_origin))
  pi_vote[,2] <- logit_inv(vote_mat_mis%*%t(alpha_v_origin))
  
  # p(R_E|V)
  REV_prob <- pi_vote
  r_eth_mis <- glm(Re~as.factor(sex) + as.factor(educ) + as.factor(age) + as.factor(voted), 
                   family = binomial(),data = sub_vote)
  r_eth_mis_mat <- model.matrix(r_eth_mis)
  r_eth_mis_mat[,11] <- 0
  REV_0_pred <- logit_inv(r_eth_mis_mat%*%t(gamma_e))
  REV_prob[,1] <- ifelse(sub_vote$Re == 1, REV_0_pred, 1-REV_0_pred)
  r_eth_mis_mat[,11] <- 1
  REV_1_pred <- logit_inv(r_eth_mis_mat%*%t(gamma_e))
  REV_prob[,2] <- ifelse(sub_vote$Re == 1, REV_1_pred, 1-REV_1_pred)
  
  # p(R_A|V)
  RAV_prob <- pi_vote
  r_age_mis <- glm(Ra ~ as.factor(sex) + as.factor(race) + as.factor(educ), family = binomial(),data = sub_vote)
  r_age_mis_mat <- model.matrix(r_age_mis)
  RAV_pred <- logit_inv(r_age_mis_mat%*%t(gamma_a))
  RAV_prob[,1] <- ifelse(sub_vote$Ra == 1,RAV_pred, 1-RAV_pred)
  RAV_prob[,2] <- ifelse(sub_vote$Ra == 1,RAV_pred, 1-RAV_pred)
  
  # p(R_C|V)
  RCV_prob <- pi_vote
  r_educ_mis <- glm(Rc~as.factor(sex) + as.factor(race) + as.factor(age) +  as.factor(voted), 
                    family = binomial(),data = sub_vote)
  r_educ_mis_mat <- model.matrix(r_educ_mis)
  r_educ_mis_mat[,12] <- 0
  RCV_0_pred <- logit_inv(r_educ_mis_mat%*%t(gamma_c))
  RCV_prob[,1] <- ifelse(sub_vote$Rc == 1, RCV_0_pred, 1-RCV_0_pred)
  r_educ_mis_mat[,12] <- 1
  RCV_1_pred <- logit_inv(r_educ_mis_mat%*%t(gamma_c))
  RCV_prob[,2] <- ifelse(sub_vote$Rc == 1, RCV_1_pred, 1-RCV_1_pred)
  
  # put it all together
  pr_vote_miss <- pi_vote*REV_prob*RAV_prob*RCV_prob
  pr_vote_miss_prob <- ifelse(is.na(pr_vote_miss/rowSums(pr_vote_miss)),0,pr_vote_miss/rowSums(pr_vote_miss))
  nc_adj[which(nc_adj$Rv == 1 & nc_adj$unit == 0),]$voted <- rbinom(n_mis_v,1,pr_vote_miss_prob[,2])
  
  if (i > burnin){
    GAMMA_e[i-burnin,] <- gamma_e
    GAMMA_a[i-burnin,] <- gamma_a
    GAMMA_v[i-burnin,] <- gamma_v
    GAMMA_c[i-burnin,] <- gamma_c
    ALPHA_e[i-burnin,] <- alpha_e
    ALPHA_c[i-burnin,] <- alpha_c
    ALPHA_a[i-burnin,] <- alpha_a
    ALPHA_v[i-burnin,] <- alpha_v
    ALPHA_v_origin[i-burnin,] <- alpha_v_origin
    ALPHA_g[i-burnin,] <- alpha_g
    ALPHA_g_origin[i-burnin,] <- alpha_g_origin
    vote_impute[i-burnin,] <- nc_adj$voted
    age_impute[i-burnin,] <- nc_adj$age
    educ_impute[i-burnin,] <- nc_adj$educ
    race_impute[i-burnin,] <- nc_adj$race
    sex_impute[i-burnin,] <- nc_adj$sex
    vote_replicate[i-burnin,] <- nc_rep$voted
    age_replicate[i-burnin,] <- nc_rep$age
    educ_replicate[i-burnin,] <- nc_rep$educ
    race_replicate[i-burnin,] <- nc_rep$race
    sex_replicate[i-burnin,] <- nc_rep$sex
    Ra_replicate[i-burnin,] <- nc_rep$Ra[which(nc_rep$unit == 0)]  
    Rv_replicate[i-burnin,] <- nc_rep$Rv[which(nc_rep$unit == 0)]  
    Re_replicate[i-burnin,] <- nc_rep$Re[which(nc_rep$unit == 0)] 
    Rc_replicate[i-burnin,] <- nc_rep$Rc[which(nc_rep$unit == 0)] 
  }
}


###############################################
############ Conditional Margins ##############
###############################################
sv <- matrix(NA,5000,2)
av <- matrix(NA,5000,7)
ev <- matrix(NA,5000,4)
cv <- matrix(NA,5000,3)

vote_tab <- matrix(NA,5000,2)
sex_tab <- matrix(NA,5000,2)
age_tab <- matrix(NA,5000,7)
eth_tab <- matrix(NA,5000,4)
educ_tab <- matrix(NA,5000,3)

sev <- matrix(NA,5000,8)
sav <- matrix(NA,5000,14)
scv <- matrix(NA,5000,6)


for (i in 1:5000){
  Voted <- vote_replicate[i,which(nc_adj$unit == 0)]
  Age <- age_replicate[i,which(nc_adj$unit == 0)]
  Race <- race_replicate[i,which(nc_adj$unit == 0)]
  Sex <- sex_replicate[i,which(nc_adj$unit == 0)]
  Educ <- educ_replicate[i,which(nc_adj$unit == 0)]
  # resample R
  # Ra <- Ra_replicate[i,]
  # Re <- Re_replicate[i,]
  # Rv <- Rv_replicate[i,]
  # Rc <- Rc_replicate[i,]
  # without resample R
  Ra <- nc_adj$Ra[which(nc_adj$unit == 0)]
  Re <- nc_adj$Re[which(nc_adj$unit == 0)]
  Rv <- nc_adj$Rv[which(nc_adj$unit == 0)]
  Rc <- nc_adj$Rc[which(nc_adj$unit == 0)]
  dat_rep <- data.frame(voted = Voted, race = Race, sex = Sex, age = Age, educ = Educ, Ra = Ra, Re = Re, Rv = Rv, Rc = Rc)
  dat_repC <- dat_rep[which(dat_rep$Rv == 0 & dat_rep$Ra == 0 & dat_rep$Re == 0 & dat_rep$Rc == 0),]
  t1 <- prop.table(table(dat_repC$sex,dat_repC$voted),1)[,2]
  sv[i,] <- c(t1)
  t2 <- prop.table(table(dat_repC$age,dat_repC$voted),1)[,2]
  if (length(t2)<7){
    av[i,] <- c(t2,0)
  }
  else{
    av[i,] <- c(t2) 
  }
  t3 <- prop.table(table(dat_repC$race,dat_repC$voted),1)[,2]
  ev[i,] <- c(t3)
  
  t10 <- prop.table(table(dat_repC$educ,dat_repC$voted),1)[,2]
  cv[i,] <- c(t10)
  
  t4 <- prop.table(table(dat_repC$sex))
  sex_tab[i,] <- c(t4)
  t11 <- prop.table(table(dat_repC$educ))
  educ_tab[i,] <- c(t11)
  t5 <- prop.table(table(dat_repC$age))
  if (length(t5) < 7){
    age_tab[i,] <- c(t5,0)
  }
  else{
    age_tab[i,] <- c(t5) 
  }
  t6 <- prop.table(table(dat_repC$race))
  eth_tab[i,] <- c(t6)
  t7 <- prop.table(table(dat_repC$vote))
  vote_tab[i,] <- c(t7)
  
  t8 <- prop.table(table(dat_repC$sex,dat_repC$age,dat_repC$voted),c(1,2))[,,2]
  if (length(t8)<14){
    sav[i,] <- c(t8,0,0)
  }
  else{
    sav[i,] <- c(t8) 
  }
  t9 <- prop.table(table(dat_repC$sex,dat_repC$race,dat_repC$voted),c(1,2))[,,2]
  sev[i,] <- c(t9)
  t12 <- prop.table(table(dat_repC$sex,dat_repC$educ,dat_repC$voted),c(1,2))[,,2]
  scv[i,] <- c(t12)
}

nc_adjC <- nc_copy[which(nc_copy$Rv == 0 & nc_copy$Ra == 0 & nc_copy$Re == 0 & nc_copy$Rc == 0 & nc_copy$unit == 0),]
t1_origin <- prop.table(table(nc_adjC$sex,nc_adjC$voted),1)[,2]
sv_origin <- c(t1_origin)

t2_origin <- prop.table(table(nc_adjC$age,nc_adjC$voted),1)[,2]
av_origin <- c(t2_origin)

t3_origin <- prop.table(table(nc_adjC$race,nc_adjC$voted),1)[,2]
ev_origin <- c(t3_origin)

t9_origin <- prop.table(table(nc_adjC$educ,nc_adjC$voted),1)[,2]
cv_origin <- c(t9_origin)

plot_coef_sv <- matrix(NA,ncol = 2, nrow = 2)
for (i in 1:2){
  plot_coef_sv[i,] = quantile(sv[,i],c(0.025,0.975))
}

G.1 = as.data.frame(cbind(c(sv_origin),plot_coef_sv))
colnames(G.1) = c("truth","inf","sup")

plot_coef_ev <- matrix(NA,ncol = 2, nrow = 4)
for (i in 1:4){
  plot_coef_ev[i,] = quantile(ev[,i],c(0.025,0.975))
}

G.2 = as.data.frame(cbind(c(ev_origin),plot_coef_ev))
colnames(G.2) = c("truth","inf","sup")

plot_coef_av <- matrix(NA,ncol = 2, nrow = 7)
for (i in 1:7){
  plot_coef_av[i,] = quantile(av[,i],c(0.025,0.975))
}

G.3 = as.data.frame(cbind(c(av_origin),plot_coef_av))
colnames(G.3) = c("truth","inf","sup")

plot_coef_cv <- matrix(NA,ncol = 2, nrow = 3)
for (i in 1:3){
  plot_coef_cv[i,] = quantile(cv[,i],c(0.025,0.975))
}

G.9 = as.data.frame(cbind(c(cv_origin),plot_coef_cv))
colnames(G.9) = c("truth","inf","sup")

merge_dat1 <- as.data.frame(rbind(G.1,G.2,G.3,G.9))

# Three way
sav_tab_origin <- prop.table(table(nc_adjC$sex,nc_adjC$age,nc_adjC$voted),c(1,2))[,,2]
sav_origin <- c(sav_tab_origin)

sev_tab_origin <-  prop.table(table(nc_adjC$sex,nc_adjC$race,nc_adjC$voted),c(1,2))[,,2]
sev_origin <- c(sev_tab_origin)

scv_tab_origin <-  prop.table(table(nc_adjC$sex,nc_adjC$educ,nc_adjC$voted),c(1,2))[,,2]
scv_origin <- c(scv_tab_origin)

plot_coef_sev <- matrix(NA,ncol = 2, nrow = 8)
for (i in 1:8){
  plot_coef_sev[i,] = quantile(sev[,i],c(0.025,0.975))
}

plot_coef_sav <- matrix(NA,ncol = 2, nrow = 14)
for (i in 1:14){
  plot_coef_sav[i,] = quantile(sav[,i],c(0.025,0.975),na.rm = T)
}

plot_coef_scv <- matrix(NA,ncol = 2, nrow = 6)
for (i in 1:6){
  plot_coef_scv[i,] = quantile(scv[,i],c(0.025,0.975),na.rm = T)
}

G.5 = as.data.frame(cbind(c(sev_origin,sav_origin,scv_origin),
                          rbind(plot_coef_sev,plot_coef_sav,plot_coef_scv)))
colnames(G.5) = c("truth","inf","sup")

########### Marginal ##############
sex_t1_origin <- prop.table(table(nc_adjC$sex))
sex_tab_origin <- c(sex_t1_origin)

age_t2_origin <- prop.table(table(nc_adjC$age))
age_tab_origin <- c(age_t2_origin)

eth_t3_origin <- prop.table(table(nc_adjC$race))
eth_tab_origin <- c(eth_t3_origin)

vote_t4_origin <- prop.table(table(nc_adjC$voted))
vote_tab_origin <- c(vote_t4_origin)

educ_t9_origin <- prop.table(table(nc_adjC$educ))
educ_tab_origin <- c(educ_t9_origin)


plot_coef_sex <- matrix(NA,ncol = 2, nrow = 2)
for (i in 1:2){
  plot_coef_sex[i,] = quantile(sex_tab[,i],c(0.025,0.975))
}

plot_coef_vote <- matrix(NA,ncol = 2, nrow = 2)
for (i in 1:2){
  plot_coef_vote[i,] = quantile(vote_tab[,i],c(0.025,0.975))
}

plot_coef_educ <- matrix(NA,ncol = 2, nrow = 3)
for (i in 1:3){
  plot_coef_educ[i,] = quantile(educ_tab[,i],c(0.025,0.975))
}

plot_coef_age <- matrix(NA,ncol = 2, nrow = 7)
for (i in 1:7){
  plot_coef_age[i,] = quantile(age_tab[,i],c(0.025,0.975))
}

plot_coef_eth <- matrix(NA,ncol = 2, nrow = 4)
for (i in 1:4){
  plot_coef_eth[i,] = quantile(eth_tab[,i],c(0.025,0.975))
}

G.4 = as.data.frame(cbind(c(vote_tab_origin,sex_tab_origin,eth_tab_origin,age_tab_origin,educ_tab_origin),
                          rbind(plot_coef_vote,plot_coef_sex,plot_coef_eth,plot_coef_age,plot_coef_educ)))
colnames(G.4) = c("truth","inf","sup")
G.4$truth < G.4$sup & G.4$truth > G.4$inf

# Final data
final_dat <- as.data.frame(rbind(G.4, merge_dat1,G.5))
final_dat$genre <- 1:62
final_dat$include <- final_dat$truth < final_dat$sup & final_dat$truth > final_dat$inf

ggplot(final_dat, aes(x = genre, y = truth)) +
  geom_point(shape = 4, size = 0.8) +
  geom_errorbar(aes(ymax = sup, ymin = inf, linetype = include),width = 0.5)+
  geom_vline(xintercept = c(18.5,34.5), col="black", linetype = "dotted") + 
  # scale_colour_manual(name="Error Bars",values=cols) + scale_fill_manual(name="Bar",values=cols) +
  ylab("Probability")+
  xlab("Index") + theme_classic() + theme(legend.position = "none",axis.title.x=element_blank(),
                                          axis.text.x=element_blank(),axis.ticks.x=element_blank())

ggplot(final_dat, aes(x = genre, y = truth, color = include)) +
  geom_point(shape = 4, color = "indianred4", size = 0.8) +
  geom_errorbar(aes(ymax = sup, ymin = inf, color = include, linetype = include),width = 0.5)+
  geom_vline(xintercept = c(18.5,34.5), col="black", linetype = "dotted") + 
  # scale_colour_manual(name="Error Bars",values=cols) + scale_fill_manual(name="Bar",values=cols) +
  ylab("Probability")+
  xlab("Index") + theme_classic() + theme(legend.position = "none",axis.title.x=element_blank(),
                                          axis.text.x=element_blank(),axis.ticks.x=element_blank())

##########################################
########### MI combination rule ##########
##########################################
Vote <- vote_impute[seq(1,(niter-burnin),100),]
Age <- age_impute[seq(1,(niter-burnin),100),]
Race <- race_impute[seq(1,(niter-burnin),100),]
Sex <- sex_impute[seq(1,(niter-burnin),100),]
Educ <- educ_impute[seq(1,(niter-burnin),100),]

n <- dim(Vote)[1]
ans_alpha_G <- ul_a_G <- rep(NA,n)
ans_alpha_E <- ul_a_E <- matrix(NA,n,length(alpha_e))
ans_alpha_A <- ul_a_A <- matrix(NA,n,length(alpha_a))
ans_alpha_V <- ul_a_V <- matrix(NA,n,length(alpha_v))
ans_gamma_E <- ul_g_E <- matrix(NA,n,length(gamma_e))
ans_gamma_A <- ul_g_A <- matrix(NA,n,length(gamma_a))
ans_gamma_V <- ul_g_V <- matrix(NA,n,length(gamma_v))
# ans_nu <- ul_nu <- matrix(NA,n,length(nu))
vote_total <- vote_var <- rep(NA,n)
sex_total <- sex_var <- rep(NA,n)
white_total <- black_total <- hist_total <- rest_total <- rep(NA,n)
white_var <- black_var <- hist_var <- rest_var <- rep(NA,n)

for (i in 1:n){
  test_dat <- as.data.frame(cbind(Vote[i,],Age[i,],Race[i,],Sex[i,],Educ[i,],nc_copy$unit,nc_copy$weight,nc_copy$Ra,nc_copy$Rv,nc_copy$Re,rownames(nc_copy)))
  names(test_dat) <- c("voted","age","race","sex","educ","unit","weight","Ra","Rv","Re","id")
  test_dat$voted <- as.numeric(test_dat$voted)
  test_dat$race <- as.numeric(test_dat$race)
  test_dat$sex <- as.numeric(test_dat$sex)
  test_dat$educ <- as.numeric(test_dat$educ)
  test_dat$unit <- as.numeric(test_dat$unit)
  test_dat$weight <- as.numeric(test_dat$weight)
  test_dat$Ra <- as.numeric(test_dat$Ra)
  test_dat$Re <- as.numeric(test_dat$Re)
  test_dat$id <- as.numeric(test_dat$id)
  mydesign <- svydesign(id = ~id, data = test_dat, weight = ~weight)
  m_vote <- svyglm(voted ~ as.factor(sex) + age + as.factor(race) + as.factor(educ) +
                     as.factor(sex)*as.factor(race)
                   + as.factor(sex)*age + as.factor(sex)*as.factor(educ), 
                   data = test_dat,design = mydesign, family=quasibinomial())
  m_unit <- svyglm(unit ~ as.factor(sex) + as.factor(race) + as.factor(voted) , data = test_dat, design = mydesign, family=quasibinomial()) 
  ans_alpha_V[i,] <- coef(m_vote)
  ul_a_V[i,] <- diag(vcov(m_vote))
  # ans_nu[i,] <- coef(m_unit)
  # ul_nu[i,] <- diag(vcov(m_unit))
  vote_total[i] <- sum(test_dat$voted*test_dat$weight)
  vote_var[i] <- sum((test_dat$voted*test_dat$weight)^2*(1-1/test_dat$weight))
  sex_total[i] <- sum(test_dat$sex*test_dat$weight)
  sex_var[i] <- sum((test_dat$sex*test_dat$weight)^2*(1-1/test_dat$weight))
  white_total[i] <- sum(test_dat$weight[which(test_dat$race == 1)])
  white_var[i] <- sum(test_dat$weight[which(test_dat$race == 1)]^2*(1-1/test_dat$weight)[which(test_dat$race == 1)])
  black_total[i] <- sum(test_dat$weight[which(test_dat$race == 2)])
  black_var[i] <- sum(test_dat$weight[which(test_dat$race == 2)]^2*(1-1/test_dat$weight)[which(test_dat$race == 2)])
  hist_total[i] <- sum(test_dat$weight[which(test_dat$race == 3)])
  hist_var[i] <- sum(test_dat$weight[which(test_dat$race == 3)]^2*(1-1/test_dat$weight)[which(test_dat$race == 3)])
  rest_total[i] <- sum(test_dat$weight[which(test_dat$race == 4)])
  rest_var[i] <- sum(test_dat$weight[which(test_dat$race == 4)]^2*(1-1/test_dat$weight)[which(test_dat$race == 4)])
}

# vote model
u_l_aV <- colMeans(ul_a_V) 
b_l_aV <- apply((scale(ans_alpha_V,scale=FALSE))^2/(n-1),2,sum)
T_l_aV <- (1+1/n)*b_l_aV + u_l_aV # variance of alpha_V
colMeans(ans_alpha_V) # mean of alpha_V

# totals
u_l_Tv <- mean(vote_var)
b_l_Tv <- apply((scale(vote_total,scale=FALSE))^2/(n-1),2,sum)
T_l_Tv <- (1+1/n)*b_l_Tv + u_l_Tv # variance of vote
mean(vote_total)

u_l_Tg <- mean(sex_var)
b_l_Tg <- apply((scale(sex_total,scale=FALSE))^2/(n-1),2,sum)
T_l_Tg <- (1+1/n)*b_l_Tg + u_l_Tg # variance of sex
mean(sex_total)

u_l_TrW <- mean(white_var)
b_l_TrW <- apply((scale(white_total,scale=FALSE))^2/(n-1),2,sum)
T_l_TrW <- (1+1/n)*b_l_TrW + u_l_TrW # variance of white
mean(white_total)

u_l_TrB <- mean(black_var)
b_l_TrB <- apply((scale(black_total,scale=FALSE))^2/(n-1),2,sum)
T_l_TrB <- (1+1/n)*b_l_TrB + u_l_TrB # variance of black
mean(black_total)

u_l_TrH <- mean(hist_var)
b_l_TrH <- apply((scale(hist_total,scale=FALSE))^2/(n-1),2,sum)
T_l_TrH <- (1+1/n)*b_l_TrH + u_l_TrH # variance of black
mean(hist_total)

u_l_TrR <- mean(rest_var)
b_l_TrR <- apply((scale(rest_total,scale=FALSE))^2/(n-1),2,sum)
T_l_TrR <- (1+1/n)*b_l_TrR + u_l_TrR # variance of black
mean(rest_total)

###########################################################
# voter turnout across subgroups
ans_S <- matrix(NA,50,2)
ul_S <- matrix(NA,50,2)
ans_R <- matrix(NA,50,4)
ul_R <- matrix(NA,50,4)
ans_A <- matrix(NA,50,7)
ul_A <- matrix(NA,50,7)
ans_C <- matrix(NA,50,3)
ul_C <- matrix(NA,50,3)
ans <- ul <- rep(NA,50)

for (i in 1:n){
  test_dat <- as.data.frame(cbind(Vote[i,],Age[i,],Race[i,],Sex[i,],Educ[i,],
                                  nc_copy$unit,nc_copy$weight,nc_copy$Ra,nc_copy$Rv,
                                  nc_copy$Re,nc_copy$Rc,rownames(nc_copy)))
  names(test_dat) <- c("voted","age","race","sex","educ","unit","weight","Ra","Rv","Re","Rc","id")
  test_dat$voted <- as.numeric(test_dat$voted)
  test_dat$race <- as.numeric(test_dat$race)
  test_dat$sex <- as.numeric(test_dat$sex)
  test_dat$educ <- as.numeric(test_dat$educ)
  test_dat$unit <- as.numeric(test_dat$unit)
  test_dat$weight <- as.numeric(test_dat$weight)
  test_dat$Ra <- as.numeric(test_dat$Ra)
  test_dat$Re <- as.numeric(test_dat$Re)
  test_dat$Rc <- as.numeric(test_dat$Rc)
  test_dat$id <- as.numeric(test_dat$id)
  mydesign <- svydesign(id = ~id, data = test_dat, weight = ~weight)
  wgt_tbl <-svymean(~interaction(sex,voted), design=mydesign)
  wgt_tblR <- svyby(~as.factor(voted),~as.factor(race),mydesign,svymean)
  wgt_dfR <- data.frame(wgt_tblR)
  ans_R[i,] <- wgt_dfR[,3]
  ul_R[i,] <- wgt_dfR[,5]
  wgt_tblA <- svyby(~as.factor(voted),~as.factor(age),mydesign,svymean)
  wgt_dfA <- data.frame(wgt_tblA)
  ans_A[i,] <- wgt_dfA[,3]
  ul_A[i,] <- wgt_dfA[,5]
  wgt_tblS <- svyby(~as.factor(voted),~as.factor(sex),mydesign,svymean)
  wgt_dfS <- data.frame(wgt_tblS)
  ans_S[i,] <- wgt_dfS[,3]
  ul_S[i,] <- wgt_dfS[,5]
  wgt_tblC <- svyby(~as.factor(voted),~as.factor(educ),mydesign,svymean)
  wgt_dfC <- data.frame(wgt_tblC)
  ans_C[i,] <- wgt_dfC[,3]
  ul_C[i,] <- wgt_dfC[,5]
  wgt_df <- data.frame(svymean(~voted, design=mydesign))
  ans[i] <- wgt_df[,1]
  ul[i] <- wgt_df[,2]
}
colMeans(ans_S)
u_L_S <- colMeans(ul_S^2)
b_L_S <- apply((scale(ans_S,scale=FALSE))^2/(50-1),2,sum)
T_L_S <- (1+1/50)*b_L_S + u_L_S
sqrt(T_L_S)
sqrt(b_L_S/50)
sqrt(u_L_S/50)

colMeans(ans_R)
u_L_R <- colMeans(ul_R^2)
b_L_R <- apply((scale(ans_R,scale=FALSE))^2/(50-1),2,sum)
T_L_R <- (1+1/50)*b_L_R + u_L_R
sqrt(T_L_R)
sqrt(b_L_R/50)
sqrt(u_L_R/50)

colMeans(ans_A)
u_L_A <- colMeans(ul_A^2)
b_L_A <- apply((scale(ans_A,scale=FALSE))^2/(50-1),2,sum)
T_L_A <- (1+1/50)*b_L_A + u_L_A
sqrt(T_L_A)
sqrt(b_L_A/50)
sqrt(u_L_A/50)

colMeans(ans_C)
u_L_C <- colMeans(ul_C^2)
b_L_C <- apply((scale(ans_C,scale=FALSE))^2/(50-1),2,sum)
T_L_C <- (1+1/50)*b_L_C + u_L_C
sqrt(T_L_C)
sqrt(b_L_C/50)
sqrt(u_L_C/50)

mean(ans)
u_L <- mean(ul^2)
b_L <- apply((scale(ans,scale=FALSE))^2/(50-1),2,sum)
T_L <- (1+1/50)*b_L + u_L
sqrt(T_L)
sqrt(b_L/50)
sqrt(u_L/50)

# marginal distributions
ans_S <- matrix(NA,50,2)
ul_S <- matrix(NA,50,2)
ans_R <- matrix(NA,50,4)
ul_R <- matrix(NA,50,4)
ans_C <- matrix(NA,50,3)
ul_C <- matrix(NA,50,3)
ans_A <- matrix(NA,50,7)
ul_A <- matrix(NA,50,7)

for (i in 1:n){
  test_dat <- as.data.frame(cbind(Vote[i,],Age[i,],Race[i,],Sex[i,],Educ[i,],
                                  nc_copy$unit,nc_copy$weight,nc_copy$Ra,nc_copy$Rv,
                                  nc_copy$Re,nc_copy$Rc,rownames(nc_copy)))
  names(test_dat) <- c("voted","age","race","sex","educ","unit","weight","Ra","Rv","Re","Rc","id")
  test_dat$voted <- as.numeric(test_dat$voted)
  test_dat$race <- as.numeric(test_dat$race)
  test_dat$sex <- as.numeric(test_dat$sex)
  test_dat$educ <- as.numeric(test_dat$educ)
  test_dat$unit <- as.numeric(test_dat$unit)
  test_dat$weight <- as.numeric(test_dat$weight)
  test_dat$Ra <- as.numeric(test_dat$Ra)
  test_dat$Re <- as.numeric(test_dat$Re)
  test_dat$Rc <- as.numeric(test_dat$Rc)
  test_dat$id <- as.numeric(test_dat$id)
  mydesign <- svydesign(id = ~id, data = test_dat, weight = ~weight)
  wgt_df_R <- data.frame(svymean(~as.factor(race), design=mydesign))
  ans_R[i,] <- wgt_df_R[,1]
  ul_R[i,] <- wgt_df_R[,2]
  wgt_df_S <- data.frame(svymean(~as.factor(sex), design=mydesign))
  ans_S[i,] <- wgt_df_S[,1]
  ul_S[i,] <- wgt_df_S[,2]
  wgt_df_A <- data.frame(svymean(~as.factor(age), design=mydesign))
  ans_A[i,] <- wgt_df_A[,1]
  ul_A[i,] <- wgt_df_A[,2]
  wgt_df_C <- data.frame(svymean(~as.factor(educ), design=mydesign))
  ans_C[i,] <- wgt_df_C[,1]
  ul_C[i,] <- wgt_df_C[,2]
}
colMeans(ans_S)
u_L_S <- colMeans(ul_S^2)
b_L_S <- apply((scale(ans_S,scale=FALSE))^2/(50-1),2,sum)
T_L_S <- (1+1/50)*b_L_S + u_L_S
sqrt(T_L_S)
sqrt(b_L_S/50)
sqrt(u_L_S/50)

colMeans(ans_R)
u_L_R <- colMeans(ul_R^2)
b_L_R <- apply((scale(ans_R,scale=FALSE))^2/(50-1),2,sum)
T_L_R <- (1+1/50)*b_L_R + u_L_R
sqrt(T_L_R)
sqrt(b_L_R/50)
sqrt(u_L_R/50)

colMeans(ans_A)
u_L_A <- colMeans(ul_A^2)
b_L_A <- apply((scale(ans_A,scale=FALSE))^2/(50-1),2,sum)
T_L_A <- (1+1/50)*b_L_A + u_L_A
sqrt(T_L_A)
sqrt(b_L_A/50)
sqrt(u_L_A/50)

colMeans(ans_C)
u_L_C <- colMeans(ul_C^2)
b_L_C <- apply((scale(ans_C,scale=FALSE))^2/(50-1),2,sum)
T_L_C <- (1+1/50)*b_L_C + u_L_C
sqrt(T_L_C)
sqrt(b_L_C/50)
sqrt(u_L_C/50)

## voter turnout by two subgroups

### 1. voted by sex and race
ans_SR <- matrix(NA,50,8)
ul_SR <- matrix(NA,50,8)
ans_SA <- matrix(NA,50,14)
ul_SA <- matrix(NA,50,14)
ans_SC <- matrix(NA,50,6)
ul_SC <- matrix(NA,50,6)


for (i in 1:n){
  test_dat <- as.data.frame(cbind(Vote[i,],Age[i,],Race[i,],Sex[i,],Educ[i,],
                                  nc_copy$unit,nc_copy$weight,nc_copy$Ra,nc_copy$Rv,
                                  nc_copy$Re,nc_copy$Rc,rownames(nc_copy)))
  names(test_dat) <- c("voted","age","race","sex","educ","unit","weight","Ra","Rv","Re","Rc","id")
  test_dat$voted <- as.numeric(test_dat$voted)
  test_dat$race <- as.numeric(test_dat$race)
  test_dat$sex <- as.numeric(test_dat$sex)
  test_dat$educ <- as.numeric(test_dat$educ)
  test_dat$unit <- as.numeric(test_dat$unit)
  test_dat$weight <- as.numeric(test_dat$weight)
  test_dat$Ra <- as.numeric(test_dat$Ra)
  test_dat$Re <- as.numeric(test_dat$Re)
  test_dat$Rc <- as.numeric(test_dat$Rc)
  test_dat$id <- as.numeric(test_dat$id)
  mydesign <- svydesign(id = ~id, data = test_dat, weight = ~weight)
  wgt_df_SR <- data.frame(svyby(~as.factor(voted),~as.factor(race)+~as.factor(sex),mydesign,svymean))
  ans_SR[i,] <- wgt_df_SR[,4]
  ul_SR[i,] <- wgt_df_SR[,6]
  wgt_df_SA <- data.frame(svyby(~as.factor(voted),~as.factor(age)+~as.factor(sex),mydesign,svymean))
  ans_SA[i,] <- wgt_df_SA[,4]
  ul_SA[i,] <- wgt_df_SA[,6]
  wgt_df_SC <- data.frame(svyby(~as.factor(voted),~as.factor(educ)+~as.factor(sex),mydesign,svymean))
  ans_SC[i,] <- wgt_df_SC[,4]
  ul_SC[i,] <- wgt_df_SC[,6]
}
colMeans(ans_SR)
u_L_SR <- colMeans(ul_SR^2)
b_L_SR <- apply((scale(ans_SR,scale=FALSE))^2/(50-1),2,sum)
T_L_SR <- (1+1/50)*b_L_SR + u_L_SR
sqrt(T_L_SR)
sqrt(b_L_SR/50)
sqrt(u_L_SR/50)

colMeans(ans_SA)
u_L_SA <- colMeans(ul_SA^2)
b_L_SA <- apply((scale(ans_SA,scale=FALSE))^2/(50-1),2,sum)
T_L_SA <- (1+1/50)*b_L_SA + u_L_SA
sqrt(T_L_SA)
sqrt(b_L_SA/50)
sqrt(u_L_SA/50)

colMeans(ans_SC)
u_L_SC <- colMeans(ul_SC^2)
b_L_SC <- apply((scale(ans_SC,scale=FALSE))^2/(50-1),2,sum)
T_L_SC <- (1+1/50)*b_L_SC + u_L_SC
sqrt(T_L_SC)
sqrt(b_L_SC/50)
sqrt(u_L_SC/50)

##############################################
# Estimate coef of U
Vote <- vote_impute[seq(1,(niter-burnin),100),]
Age <- age_impute[seq(1,(niter-burnin),100),]
Race <- race_impute[seq(1,(niter-burnin),100),]
Sex <- sex_impute[seq(1,(niter-burnin),100),]

# calculate beta
ans_vote <- matrix(NA,50,12)
ul_vote <- matrix(NA,50,12)

# calculate gamma using stats package
ans_sex <- matrix(NA,n,2)
ul_sex <- matrix(NA,n,2)

for (i in 1:50){
  data_MI <- as.data.frame(cbind(Vote[i,],Age[i,],Race[i,],Sex[i,],nc_adj$unit,nc_adj$weight))
  names(data_MI) <- c("vote","age","race","sex","U","W")
  data_MI$vote <- as.factor(data_MI$vote)
  data_MI$race <- as.factor(data_MI$race)
  data_MI$sex <- as.factor(data_MI$sex)
  data_MI$U <- as.factor(data_MI$U)
  data_MI$W <- as.numeric(data_MI$W)
  data_MI$id <- row.names(data_MI)
  
  mydesign <- svydesign(id = ~id,data = data_MI,weight = ~W)
  
  m1 <- svyglm(vote~sex+race+age+U,design = mydesign,family = quasibinomial(link = "probit"))
  # m1 <- svyglm(vote~sex+race+age+U,design = mydesign,family = quasibinomial(link = "logit"))
  m2 <- svyglm(sex~U,design = mydesign,family = quasibinomial(link = "probit"))
  # m3 <- svymultinom(race~sex+age+U,design = mydesign,family = quasibinomial(link = "probit"))
  
  ans_vote[i,] <- coef(m1)
  ans_sex[i,] <- coef(m2)
  ul_vote[i,] <- diag(vcov(m1))
  ul_sex[i,] <- diag(vcov(m2))
}

colMeans(ans_vote)
u_L_vote <- colMeans(ul_vote)
b_L_vote <- apply((scale(ans_vote,scale=FALSE))^2/(50-1),2,sum)
T_L_vote <- (1+1/50)*b_L_vote + u_L_vote

colMeans(ans_sex)
u_L_sex <- colMeans(ul_sex)
b_L_sex <- apply((scale(ans_sex,scale=FALSE))^2/(50-1),2,sum)
T_L_sex <- (1+1/50)*b_L_sex + u_L_sex
