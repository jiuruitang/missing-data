load("./Chp4_pop_dat1.RData")


run_index = commandArgs(trailingOnly = TRUE)
V_adj = 1

# libraries
library(mvtnorm)
library(survey)
V_adj <- 1

getSample_Poisson <- function(population = pop_dat_Poisson){
  N <- dim(population)[1]
  pi <- population$pi*100  # increase inclusion probability 100 times to get enough sample
  getSampled <- rbinom(N,1,pi)
  sub_dat <- population[getSampled==1,]
  return(sub_dat)
}

logit_inv <- function(x){
  return(exp(x)/(1+exp(x)))
}

logit <- function(x){
  return(log(x)-log(1-x))
}

ANWC_unit_full<-function(niter = 10000,burnin = 5000, data = sub_dat, alpha = alpha_s, 
                         gamma = gamma_s, pop_mean = pop_mean_HT, pop_sd = pop_sd_HT, 
                         popY_mean = popY_mean_HT, popY_sd = popY_sd_HT){
  
  n_unit0 <- sum(sub_dat$U == 0)
  n_unit1 <- sum(sub_dat$U == 1)
  
  GAMMA <- matrix(NA,niter-burnin,2)
  ALPHA <- matrix(NA,niter-burnin,2)
  # ALPHA_P <- matrix(NA,niter-burnin,length(alpha_p))
  BETA <- matrix(NA,niter-burnin,2)
  X1_impute <- matrix(NA,niter-burnin,length(sub_dat$X1))
  Y_impute <- matrix(NA,niter-burnin,length(sub_dat$Y))
  
  for (i in 1:niter){
    
    # update alpha
    m_Y <- glm(Y~1, data = sub_dat[which(sub_dat$U == 0),],family=binomial())
    Y_mat <- model.matrix(m_Y)
    alpha_Y <- alpha_Y_mod <- rmvnorm(1,mean = coef(m_Y),sigma = vcov(m_Y))
    
    # impute Y for U==1 cases
    T_Y_hat <- rnorm(1,popY_mean,popY_sd)
    total_one_Y <- floor((T_Y_hat - sum((sub_dat$Y*sub_dat$W)[which(sub_dat$U == 0)]))/mean(sub_dat$W[which(sub_dat$U==1)]))
    m_Y1 <- glm(Y~1,data = sub_dat[which(sub_dat$U == 1),],family=binomial())
    Y_mat1 <- model.matrix(m_Y1)
    # coef in front of Y
    if(total_one_Y/n_unit1 < 1 && total_one_Y/n_unit1 > 0){
      u_Y <- logit(total_one_Y/n_unit1) - mean(Y_mat1%*%t(alpha_Y))
      alpha_Y_mod[1] <- alpha_Y_mod[1] + u_Y
      Y_prob <- logit_inv(Y_mat1%*%t(alpha_Y_mod))
      sub_dat$Y[which(sub_dat$U == 1)] <- rbinom(n_unit1,1,prob = Y_prob)
    }
    else{
      if(total_one_Y > n_unit1){
        sub_dat$Y[which(sub_dat$U == 1)] <- rep(1,n_unit1)
      }
      if(total_one_Y <= 0){
        sub_dat$Y[which(sub_dat$U == 1)] <- rep(0,n_unit1)
      }
    }
    
    # update beta
    m_X <- glm(X1~Y, data = sub_dat[which(sub_dat$U == 0),],family=binomial())
    X_mat <- model.matrix(m_X)
    beta_X <- beta_X_mod <- rmvnorm(1,mean = coef(m_X),sigma = vcov(m_X))
    
    # impute X for U==1 cases
    T_X_hat <- rnorm(1,pop_mean,pop_sd)
    total_one_X <- floor((T_X_hat - sum((sub_dat$X1*sub_dat$W)[which(sub_dat$U == 0)]))/mean(sub_dat$W[which(sub_dat$U==1)]))
    m_X1 <- glm(X1~Y, data = sub_dat[which(sub_dat$U == 1),],family=binomial())
    X_mat1 <- model.matrix(m_X1)
    # coef in front of X
    if(total_one_X/n_unit1 < 1 && total_one_X/n_unit1 > 0){
      u_X <- logit(total_one_X/n_unit1) - mean(X_mat1%*%t(beta_X))
      beta_X_mod[1] <- beta_X_mod[1] + u_X
      X_prob <- logit_inv(X_mat1%*%t(beta_X_mod))
      sub_dat$X1[which(sub_dat$U == 1)] <- rbinom(n_unit1,1,prob = X_prob)
    }
    else{
      if(total_one_X > n_unit1){
        sub_dat$X1[which(sub_dat$U == 1)] <- rep(1,n_unit1)
      }
      if(total_one_X <= 0){
        sub_dat$X1[which(sub_dat$U == 1)] <- rep(0,n_unit1)
      }
    }
    
    # update gamma
    m_Rx <- glm(Rx~Y, data = sub_dat[which(sub_dat$U == 0),],family=binomial())
    R_mat <- model.matrix(m_Rx)
    gamma <- rmvnorm(1,mean = c(t(coef(m_Rx))),sigma = vcov(m_Rx))
    
    # sample item non-response obs
    n_mis_X <- sum(sub_dat[which(sub_dat$U==0),]$Rx == 1)
    pr_X_miss <- matrix(0,ncol=2,nrow=n_mis_X)
    colnames(pr_X_miss) <- c("1","2")
    
    # subset of missing X1 in unit respondents
    sub_X <- sub_dat[which(sub_dat$Rx == 1 & sub_dat$U == 0),]
    subX_mat <- model.matrix(glm(X1 ~ Y, family = binomial(), data = sub_X))
    
    pi_X <- pr_X_miss
    #pi_X[,1] <- 1-predict(m_X,newdata = sub_X, type = "response") 
    #pi_X[,2] <- predict(m_X,newdata = sub_X, type = "response") 
    pi_X[,1] <- 1-logit_inv(subX_mat%*%t(beta_X))
    pi_X[,2] <- logit_inv(subX_mat%*%t(beta_X))
    
    pr_X_miss <- pi_X
    pr_X_miss_prob <- pr_X_miss/rowSums(pr_X_miss)
    sub_dat[which(sub_dat$Rx == 1 & sub_dat$U == 0),]$X1 <- rbinom(n_mis_X,1,pr_X_miss_prob[,2])
    
    if (i > burnin){
      GAMMA[i-burnin,] <- gamma
      ALPHA[i-burnin,] <- alpha_Y
      BETA[i-burnin,] <- beta_X
      X1_impute[i-burnin,] <- sub_dat$X1
      Y_impute[i-burnin,] <- sub_dat$Y
    }
  }
  XL <- X1_impute[seq(1,(niter-burnin),100),]
  YL <- Y_impute[seq(1,(niter-burnin),100),]
  return(list(gamma = GAMMA,alpha = ALPHA, beta = BETA,MI_dataX = XL, MI_dataY = YL))
}

ANWC_unitMAR<-function(niter = 10000,burnin = 5000, data = sub_dat, alpha = alpha_s, 
                       gamma = gamma_s, pop_mean = pop_mean_HT, pop_sd = pop_sd_HT, 
                       popY_mean = popY_mean_HT, popY_sd = popY_sd_HT){
  
  n_unit0 <- sum(sub_dat$U == 0)
  n_unit1 <- sum(sub_dat$U == 1)
  
  GAMMA <- matrix(NA,niter-burnin,2)
  ALPHA <- matrix(NA,niter-burnin,2)
  # ALPHA_P <- matrix(NA,niter-burnin,length(alpha_p))
  BETA <- matrix(NA,niter-burnin,2)
  X1_impute <- matrix(NA,niter-burnin,length(sub_dat$X1))
  Y_impute <- matrix(NA,niter-burnin,length(sub_dat$Y))
  
  for (i in 1:niter){
    
    # update alpha
    m_Y <- glm(Y~1, data = sub_dat[which(sub_dat$U == 0),],family=binomial())
    Y_mat <- model.matrix(m_Y)
    alpha_Y <- alpha_Y_mod <- rmvnorm(1,mean = coef(m_Y),sigma = vcov(m_Y))
    
    # impute Y for U==1 cases
    T_Y_hat <- rnorm(1,popY_mean,popY_sd)
    total_one_Y <- floor((T_Y_hat - sum((sub_dat$Y*sub_dat$W)[which(sub_dat$U == 0)]))/mean(sub_dat$W[which(sub_dat$U==1)]))
    m_Y1 <- glm(Y~1,data = sub_dat[which(sub_dat$U == 1),],family=binomial())
    Y_mat1 <- model.matrix(m_Y1)
    
    # coef in front of Y
    Y_prob <- logit_inv(Y_mat1%*%t(alpha_Y))
    sub_dat$Y[which(sub_dat$U == 1)] <- rbinom(n_unit1,1,prob = Y_prob)
    
    # update beta
    m_X <- glm(X1~Y, data = sub_dat[which(sub_dat$U == 0),],family=binomial())
    X_mat <- model.matrix(m_X)
    beta_X <- beta_X_mod <- rmvnorm(1,mean = coef(m_X),sigma = vcov(m_X))
    
    # impute X for U==1 cases
    T_X_hat <- rnorm(1,pop_mean,pop_sd)
    total_one_X <- floor((T_X_hat - sum((sub_dat$X1*sub_dat$W)[which(sub_dat$U == 0)]))/mean(sub_dat$W[which(sub_dat$U==1)]))
    m_X1 <- glm(X1~Y, data = sub_dat[which(sub_dat$U == 1),],family=binomial())
    X_mat1 <- model.matrix(m_X1)
    
    # coef in front of X
    X_prob <- logit_inv(X_mat1%*%t(beta_X))
    sub_dat$X1[which(sub_dat$U == 1)] <- rbinom(n_unit1,1,prob = X_prob)
    
    # update gamma
    m_Rx <- glm(Rx~Y, data = sub_dat[which(sub_dat$U == 0),],family=binomial())
    R_mat <- model.matrix(m_Rx)
    gamma <- rmvnorm(1,mean = c(t(coef(m_Rx))),sigma = vcov(m_Rx))
    
    # sample item non-response obs
    n_mis_X <- sum(sub_dat[which(sub_dat$U==0),]$Rx == 1)
    pr_X_miss <- matrix(0,ncol=2,nrow=n_mis_X)
    colnames(pr_X_miss) <- c("1","2")
    
    # subset of missing X1 in unit respondents
    sub_X <- sub_dat[which(sub_dat$Rx == 1 & sub_dat$U == 0),]
    subX_mat <- model.matrix(glm(X1 ~ Y, family = binomial(), data = sub_X))
    
    pi_X <- pr_X_miss
    #pi_X[,1] <- 1-predict(m_X,newdata = sub_X, type = "response") 
    #pi_X[,2] <- predict(m_X,newdata = sub_X, type = "response") 
    pi_X[,1] <- 1-logit_inv(subX_mat%*%t(beta_X))
    pi_X[,2] <- logit_inv(subX_mat%*%t(beta_X))
    
    pr_X_miss <- pi_X
    pr_X_miss_prob <- pr_X_miss/rowSums(pr_X_miss)
    sub_dat[which(sub_dat$Rx == 1 & sub_dat$U == 0),]$X1 <- rbinom(n_mis_X,1,pr_X_miss_prob[,2])
    
    if (i > burnin){
      GAMMA[i-burnin,] <- gamma
      ALPHA[i-burnin,] <- alpha_Y
      BETA[i-burnin,] <- beta_X
      X1_impute[i-burnin,] <- sub_dat$X1
      Y_impute[i-burnin,] <- sub_dat$Y
    }
  }
  XL <- X1_impute[seq(1,(niter-burnin),100),]
  YL <- Y_impute[seq(1,(niter-burnin),100),]
  return(list(gamma = GAMMA,alpha = ALPHA, beta = BETA,MI_dataX = XL, MI_dataY = YL))
}

getResults_unit <- function(dataMI_X,dataMI_Y,alpha_len, gamma_len,beta_len, nu_len,sub_dat){
  n <- dim(dataMI_X)[1]
  # calculate alpha
  ans <- matrix(NA,n,alpha_len)
  ul <- matrix(NA,n,alpha_len)
  total_varX <- total_varY <- rep(NA,n)
  
  # calculate beta
  ans_b <- matrix(NA,n,beta_len)
  ul_b <- matrix(NA,n,beta_len)
  
  # calculate gamma using stats package
  ans_g <- matrix(NA,n,gamma_len)
  ul_g <- matrix(NA,n,gamma_len)
  
  total_X <- total_Y <- rep(NA,n)
  
  # conditional margin
  ans_X0 <- ans_X1 <- ans_Y0 <- ans_Y1 <- matrix(NA,n,2)
  ul_X0 <- ul_X1 <- ul_Y0 <- ul_Y1 <- matrix(NA,n,2)
  
  for (i in 1:n){
    test_dat <- as.data.frame(cbind(dataMI_X[i,],rownames(sub_dat),sub_dat$W,dataMI_Y[i,],sub_dat$Rx,sub_dat$U),stringsAsFactors = FALSE)
    names(test_dat) <- c("X1","id","W","Y","Rx","U")
    total_varX[i] <- sum((as.numeric(test_dat$X1)*sub_dat$W)^2*(1-1/sub_dat$W))
    total_varY[i] <- sum((as.numeric(test_dat$Y)*sub_dat$W)^2*(1-1/sub_dat$W))
    test_dat$X1 <- as.factor(test_dat$X1)
    test_dat$id <- as.numeric(test_dat$id)
    test_dat$W <- as.numeric(test_dat$W)
    test_dat$pi <- 1/test_dat$W
    test_dat$Y <- as.factor(test_dat$Y)
    test_dat$Rx <- as.factor(test_dat$Rx)
    test_dat$U <- as.factor(test_dat$U)
    mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
    m1 <- svyglm(X1~Y+U,design = mydesign,family = quasibinomial(link = "probit"))
    ans_b[i,] <- coef(m1)
    ul_b[i,] <- diag(vcov(m1))
    m3 <- svyglm(Y~U,design = mydesign,family = quasibinomial(link = "probit"))
    ans[i,] <- coef(m3)
    ul[i,] <- diag(vcov(m3))
    
    # use only unit respondents for Rx model
    m2 <- glm(Rx~Y,data = test_dat[which(test_dat$U==0),], family = binomial(link = "probit"))
    ans_g[i,] <- coef(m2)
    ul_g[i,] <- diag(vcov(m2))
    total_X[i] <- sum(dataMI_X[i,]*sub_dat$W)
    total_Y[i] <- sum(dataMI_Y[i,]*sub_dat$W)
    
    wgt_tblY <- svyby(~as.factor(X1),~as.factor(Y),mydesign,svymean)
    wgt_dfY <- data.frame(wgt_tblY)
    ans_Y0[i,] <- c(t(wgt_dfY[1,2:3]))
    ul_Y0[i,] <- c(t(wgt_dfY[1,4:5]))
    ans_Y1[i,] <- c(t(wgt_dfY[2,2:3]))
    ul_Y1[i,] <- c(t(wgt_dfY[2,4:5]))
    
    wgt_tblX <- svyby(~as.factor(Y),~as.factor(X1),mydesign,svymean)
    wgt_dfX <- data.frame(wgt_tblX)
    ans_X0[i,] <- c(t(wgt_dfX[1,2:3]))
    ul_X0[i,] <- c(t(wgt_dfX[1,4:5]))
    ans_X1[i,] <- c(t(wgt_dfX[2,2:3]))
    ul_X1[i,] <- c(t(wgt_dfX[2,4:5]))
  }
  # colMeans(ans)
  u_L_bar <- colMeans(ul)
  b_L <- apply((scale(ans,scale=FALSE))^2/(n-1),2,sum)
  T_L <- (1+1/n)*b_L + u_L_bar
  
  # colMeans(ans_g)
  u_L_g <- colMeans(ul_g)
  b_L_g <- apply((scale(ans_g,scale=FALSE))^2/(n-1),2,sum)
  T_L_g <- (1+1/n)*b_L_g + u_L_g
  
  # colMeans(ans_b)
  u_L_b <- colMeans(ul_b)
  b_L_b <- apply((scale(ans_b,scale=FALSE))^2/(n-1),2,sum)
  T_L_b <- (1+1/n)*b_L_b + u_L_b
  
  u_L_t <- mean(total_varX)
  b_L_t <- apply((scale(total_X,scale=FALSE))^2/(n-1),2,sum)
  T_L_t <- (1+1/n)*b_L_t + u_L_t
  
  u_L_tY <- mean(total_varY)
  b_L_tY <- apply((scale(total_Y,scale=FALSE))^2/(n-1),2,sum)
  T_L_tY <- (1+1/n)*b_L_tY + u_L_tY
  
  # conditional margin
  u_L_X0 <- colMeans(ul_X0^2)
  b_L_X0 <- apply((scale(ans_X0,scale=FALSE))^2/(50-1),2,sum)
  T_L_X0 <- (1+1/50)*b_L_X0 + u_L_X0
  
  u_L_X1 <- colMeans(ul_X1^2)
  b_L_X1 <- apply((scale(ans_X1,scale=FALSE))^2/(50-1),2,sum)
  T_L_X1 <- (1+1/50)*b_L_X1 + u_L_X1
  
  u_L_Y0 <- colMeans(ul_Y0^2)
  b_L_Y0 <- apply((scale(ans_Y0,scale=FALSE))^2/(50-1),2,sum)
  T_L_Y0 <- (1+1/50)*b_L_Y0 + u_L_Y0
  
  u_L_Y1 <- colMeans(ul_Y1^2)
  b_L_Y1 <- apply((scale(ans_Y1,scale=FALSE))^2/(50-1),2,sum)
  T_L_Y1 <- (1+1/50)*b_L_Y1 + u_L_Y1
  
  return(list(alpha = ans, gamma = ans_g, beta = ans_b, total_mean_X = mean(total_X), total_mean_Y = mean(total_Y) ,
              total_varX = T_L_t, total_varY = T_L_tY,
              alpha_var = T_L, alpha_mean = colMeans(ans),gamma_mean = colMeans(ans_g),beta_mean = colMeans(ans_b),
              gamma_var = T_L_g, beta_var = T_L_b,
              mar_varb = c(b_L,b_L_b,b_L_g,b_L_t,b_L_tY),
              cond_pt = c(colMeans(ans_X0), colMeans(ans_X1), colMeans(ans_Y0), colMeans(ans_Y1)),
              cond_var = c(T_L_X0, T_L_X1, T_L_Y0, T_L_Y1),
              var_b = c(b_L_X0,b_L_X1,b_L_Y0,b_L_Y1)))
}

########### End of Functions ########
sub_dat <- getSample_Poisson(population = pop_dat)
sub_dat$origin_W <- sub_dat$W
sub_dat$W <- sub_dat$W/100
sub_dat$pi <- sub_dat$pi*100

N_pop <- nrow(pop_dat)

# weight adjustment for X=1
w_adj1 <- sum(sub_dat[which(sub_dat$X1 == 1),]$W)/sum(sub_dat[which(sub_dat$U == 0 & sub_dat$X1 == 1),]$W)

# weight adjustment for X=0
w_adj0 <- sum(sub_dat[which(sub_dat$X1 == 0),]$W)/sum(sub_dat[which(sub_dat$U == 0 & sub_dat$X1 == 0),]$W)

# apply the nonresponse adjustment
sub_dat[which(sub_dat$U == 0 & sub_dat$X1 == 1),]$W <- sub_dat[which(sub_dat$U == 0 & sub_dat$X1 == 1),]$W*w_adj1
sub_dat[which(sub_dat$U == 1 & sub_dat$X1 == 1),]$W <- NA
sub_dat[which(sub_dat$U == 0 & sub_dat$X1 == 0),]$W <- sub_dat[which(sub_dat$U == 0 & sub_dat$X1 == 0),]$W*w_adj0
sub_dat[which(sub_dat$U == 1 & sub_dat$X1 == 0),]$W <- NA

# adjust weights and inclusion probabilities
# total_W <- sum(sub_dat$W,na.rm = T)
total_W = nrow(pop_dat)
sub_dat[which(sub_dat$U == 0),]$W <- sub_dat[which(sub_dat$U == 0),]$W*(1-sum(sub_dat$U)/nrow(sub_dat))
sub_dat[which(sub_dat$U == 1),]$W <- total_W/nrow(sub_dat)
sub_dat$pi <- 1/sub_dat$W

# HT estimator of population total
X1 <- sub_dat$X1
W <- sub_dat$W
Rx <- sub_dat$Rx
Y <- sub_dat$Y
U <- sub_dat$U
# pbeta <- 1/sub_dat$W
n_mis <- sum(sub_dat$Rx == 1)
m_X1 <- glm(X1~Y,data = sub_dat[which(Rx ==0 & U == 0),],family=binomial(probit))
p_X1 <- predict(m_X1,newdata = data.frame(Y=Y[which(Rx == 1||U == 1)]),type = "response")
X1[which(Rx == 1 || U == 1)] <- rbinom(length(p_X1),1,prob = p_X1)
sub_dat$X1 <- X1
pop_sd_HT <- sqrt(sum((X1/sub_dat$pi)^2*(1-sub_dat$pi)))/V_adj 
pop_mean_HT <- sum(pop_dat$X1)
popY_sd_HT <- sqrt(sum((Y/sub_dat$pi)^2*(1-sub_dat$pi)))/V_adj 
popY_mean_HT <- sum(pop_dat$Y)

pop_mean = pop_mean_HT
pop_sd = pop_sd_HT
popY_mean = popY_mean_HT
popY_sd = popY_sd_HT

alpha_s <- rnorm(2)
gamma_s <- rnorm(2)
beta_s <- rnorm(2)
nu_s <- rnorm(3)

#######################################
# n_unit0 <- sum(sub_dat$U == 0)
# n_unit1 <- sum(sub_dat$U == 1)
# niter <- 10000
# burnin <- 5000
# 
# GAMMA <- matrix(NA,niter-burnin,2)
# ALPHA <- matrix(NA,niter-burnin,2)
# # ALPHA_P <- matrix(NA,niter-burnin,length(alpha_p))
# BETA <- matrix(NA,niter-burnin,2)
# X1_impute <- matrix(NA,niter-burnin,length(sub_dat$X1))
# Y_impute <- matrix(NA,niter-burnin,length(sub_dat$Y))
# 
# for (i in 1:niter){
#   
#   # update alpha
#   m_Y <- glm(Y~1, data = sub_dat[which(sub_dat$U == 0),],family=binomial())
#   Y_mat <- model.matrix(m_Y)
#   alpha_Y <- alpha_Y_mod <- rmvnorm(1,mean = coef(m_Y),sigma = vcov(m_Y))
#   
#   # impute Y for U==1 cases
#   T_Y_hat <- rnorm(1,popY_mean,popY_sd)
#   total_one_Y <- floor((T_Y_hat - sum((sub_dat$Y*sub_dat$W)[which(sub_dat$U == 0)]))/mean(sub_dat$W[which(sub_dat$U==1)]))
#   m_Y1 <- glm(Y~1,data = sub_dat[which(sub_dat$U == 1),],family=binomial())
#   Y_mat1 <- model.matrix(m_Y1)
#   # coef in front of Y
#   if(total_one_Y/n_unit1 < 1 && total_one_Y/n_unit1 > 0){
#     u_Y <- logit(total_one_Y/n_unit1) - mean(Y_mat1%*%t(alpha_Y))
#     alpha_Y_mod[1] <- alpha_Y_mod[1] + u_Y
#     Y_prob <- logit_inv(Y_mat1%*%t(alpha_Y_mod))
#     sub_dat$Y[which(sub_dat$U == 1)] <- rbinom(n_unit1,1,prob = Y_prob)
#   }
#   else{
#     if(total_one_Y > n_unit1){
#       sub_dat$Y[which(sub_dat$U == 1)] <- rep(1,n_unit1)
#     }
#     if(total_one_Y <= 0){
#       sub_dat$Y[which(sub_dat$U == 1)] <- rep(0,n_unit1)
#     }
#   }
#   
#   # update beta
#   m_X <- glm(X1~Y, data = sub_dat[which(sub_dat$U == 0),],family=binomial())
#   X_mat <- model.matrix(m_X)
#   beta_X <- beta_X_mod <- rmvnorm(1,mean = coef(m_X),sigma = vcov(m_X))
#   
#   # impute X for U==1 cases
#   T_X_hat <- rnorm(1,pop_mean,pop_sd)
#   total_one_X <- floor((T_X_hat - sum((sub_dat$X1*sub_dat$W)[which(sub_dat$U == 0)]))/mean(sub_dat$W[which(sub_dat$U==1)]))
#   m_X1 <- glm(X1~Y, data = sub_dat[which(sub_dat$U == 1),],family=binomial())
#   X_mat1 <- model.matrix(m_X1)
#   # coef in front of X
#   if(total_one_X/n_unit1 < 1 && total_one_X/n_unit1 > 0){
#     u_X <- logit(total_one_X/n_unit1) - mean(X_mat1%*%t(beta_X))
#     beta_X_mod[1] <- beta_X_mod[1] + u_X
#     X_prob <- logit_inv(X_mat1%*%t(beta_X_mod))
#     sub_dat$X1[which(sub_dat$U == 1)] <- rbinom(n_unit1,1,prob = X_prob)
#   }
#   else{
#     if(total_one_X > n_unit1){
#       sub_dat$X1[which(sub_dat$U == 1)] <- rep(1,n_unit1)
#     }
#     if(total_one_X <= 0){
#       sub_dat$X1[which(sub_dat$U == 1)] <- rep(0,n_unit1)
#     }
#   }
#   
#   # update gamma
#   m_Rx <- glm(Rx~Y, data = sub_dat[which(sub_dat$U == 0),],family=binomial())
#   R_mat <- model.matrix(m_Rx)
#   gamma <- rmvnorm(1,mean = c(t(coef(m_Rx))),sigma = vcov(m_Rx))
#   
#   # sample item non-response obs
#   n_mis_X <- sum(sub_dat[which(sub_dat$U==0),]$Rx == 1)
#   pr_X_miss <- matrix(0,ncol=2,nrow=n_mis_X)
#   colnames(pr_X_miss) <- c("1","2")
#   
#   # subset of missing X1 in unit respondents
#   sub_X <- sub_dat[which(sub_dat$Rx == 1 & sub_dat$U == 0),]
#   subX_mat <- model.matrix(glm(X1 ~ Y, family = binomial(), data = sub_X))
#   
#   pi_X <- pr_X_miss
#   #pi_X[,1] <- 1-predict(m_X,newdata = sub_X, type = "response") 
#   #pi_X[,2] <- predict(m_X,newdata = sub_X, type = "response") 
#   pi_X[,1] <- 1-logit_inv(subX_mat%*%t(beta_X))
#   pi_X[,2] <- logit_inv(subX_mat%*%t(beta_X))
#   
#   pr_X_miss <- pi_X
#   pr_X_miss_prob <- pr_X_miss/rowSums(pr_X_miss)
#   sub_dat[which(sub_dat$Rx == 1 & sub_dat$U == 0),]$X1 <- rbinom(n_mis_X,1,pr_X_miss_prob[,2])
#   
#   if (i > burnin){
#     GAMMA[i-burnin,] <- gamma
#     ALPHA[i-burnin,] <- alpha_Y
#     BETA[i-burnin,] <- beta_X
#     X1_impute[i-burnin,] <- sub_dat$X1
#     Y_impute[i-burnin,] <- sub_dat$Y
#   }
# }
# XL <- X1_impute[seq(1,(niter-burnin),100),]
# YL <- Y_impute[seq(1,(niter-burnin),100),]
# 
# dataMI_X = XL
# dataMI_Y = YL
##############################################
testListANWC <- ANWC_unit_full()
resultListANWC <- getResults_unit(dataMI_X = testListANWC$MI_dataX, dataMI_Y = testListANWC$MI_dataY,
                                  alpha_len=2, gamma_len = 2,beta_len = 3, sub_dat = sub_dat)
testListMAR <- ANWC_unitMAR()
resultListMAR <- getResults_unit(dataMI_X = testListMAR$MI_dataX, dataMI_Y = testListMAR$MI_dataY,
                                 alpha_len=2, gamma_len = 2,beta_len = 3, sub_dat = sub_dat)

MI_dataANWC <- data.frame(X = testListANWC$MI_dataX, Y = testListANWC$MI_dataY)
MI_dataMAR <- data.frame(X = testListMAR$MI_dataX, Y = testListMAR$MI_dataY)

save(resultListANWC,resultListMAR,sub_dat,MI_dataANWC,MI_dataMAR,pop_sd_HT,file = paste("./Chp4/pop1_N/Mis_",run_index,".RData",sep=""))



