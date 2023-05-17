# Generates populations for simulation studies regarding 
# handling item and unit nonresponse with survey weights
# and auxiliary information

# libraries
library(mvtnorm)
library(survey)
V_adj <- 1

load("~/Missing/cps2012.rda")
dat12 <- da36383.0001
weights <- dat12$HWHHWGT/10000
# there are some zero weights entries (vacant)
nonzero_W <- weights[which(weights != 0)]


# this function generates population as follow:
## pattern mixture model for U
## selection model for item nonresponse
gen_pop_unitMix <- function(N=length(nonzero_W),HTW = nonzero_W,alpha = c(0.3,-0.5),
                            gamma = c(-1.05,0.2),beta=c(0.2,0.4,-0.5),nu = -1.2){
  # pi <- rbeta(N,2,5) 
  W <- HTW
  pi <- 1/W
  
  # sample U
  pi_U <- pnorm(nu)
  U <- rbinom(N,1,p=pi_U)
  
  # sample Y|U
  Y <- U
  pi_Y <- pnorm(model.matrix(Y~U)%*%alpha)
  Y <- rbinom(N,1,p=pi_Y)
  
  # sample X|Y,U
  X <- Y
  pi_X <- pnorm(model.matrix(X~Y+U)%*%beta)
  X1 <- rbinom(N,1,p=pi_X)
  
  # sample Rx|X,Y,U
  Rx <- X1
  pi_Rx <- pnorm(model.matrix(Rx~Y)%*%gamma)
  Rx <- rbinom(N,1,p = pi_Rx)
  
  # get a dataframe for population
  return(as.data.frame(cbind(Y,X1,Rx,W,pi,U)))
}

# this function generates population as follow:
## pattern mixture model for U
## selection model for item nonresponse
plogit <- function(X){
  return(1/(1+exp(-1*X)))
}

gen_pop_unitMix_logit <- function(N=length(nonzero_W),HTW = nonzero_W,alpha = c(0.3,-0.5),
                            gamma = c(-1.05,0.2),beta=c(0.2,0.4,-0.5),nu = -1.2){
  # pi <- rbeta(N,2,5) 
  W <- HTW
  pi <- 1/W
  
  # sample U
  pi_U <- plogit(nu)
  U <- rbinom(N,1,p=pi_U)
  
  # sample Y|U
  Y <- U
  pi_Y <- plogit(model.matrix(Y~U)%*%alpha)
  Y <- rbinom(N,1,p=pi_Y)
  
  # sample X|Y,U
  X <- Y
  pi_X <- plogit(model.matrix(X~Y+U)%*%beta)
  X1 <- rbinom(N,1,p=pi_X)
  
  # sample Rx|X,Y,U
  Rx <- X1
  pi_Rx <- plogit(model.matrix(Rx~Y)%*%gamma)
  Rx <- rbinom(N,1,p = pi_Rx)
  
  # get a dataframe for population
  return(as.data.frame(cbind(Y,X1,Rx,W,pi,U)))
}

pop_dat <- gen_pop_unitMix(alpha = c(0.3,-0.5),gamma = c(-1.05,0.2),beta=c(0.2,0.4,-0.5),nu = -1.2)
save(pop_dat,file = "Chp4_pop_dat1.RData")
pop_dat <- gen_pop_unitMix(alpha = c(0.3,-2.0),gamma = c(-1.05,0.2),beta=c(0.2,0.4,-0.5),nu = -1.2)
save(pop_dat,file = "Chp4_pop_dat2.RData")
pop_dat <- gen_pop_unitMix(alpha = c(0.3,-0.5),gamma = c(-1.05,0.2),beta=c(0.2,0.4,-2.0),nu = -1.2)
save(pop_dat,file = "Chp4_pop_dat3.RData")
pop_dat <- gen_pop_unitMix(alpha = c(0.3,-2.0),gamma = c(-1.05,0.2),beta=c(0.2,0.4,-2.0),nu = -1.2)
save(pop_dat,file = "Chp4_pop_dat4.RData")

pop_dat <- gen_pop_unitMix(alpha = c(0.3,-0.5),gamma = c(-1.05,0.2),beta=c(0.2,2,-0.5),nu = -1.2)
save(pop_dat,file = "Chp4_pop_dat5.RData")
pop_dat <- gen_pop_unitMix(alpha = c(0.3,-2.0),gamma = c(-1.05,0.2),beta=c(0.2,2,-0.5),nu = -1.2)
save(pop_dat,file = "Chp4_pop_dat6.RData")
pop_dat <- gen_pop_unitMix(alpha = c(0.3,-0.5),gamma = c(-1.05,0.2),beta=c(0.2,2,-2.0),nu = -1.2)
save(pop_dat,file = "Chp4_pop_dat7.RData")
pop_dat <- gen_pop_unitMix(alpha = c(0.3,-2.0),gamma = c(-1.05,0.2),beta=c(0.2,2,-2.0),nu = -1.2)
save(pop_dat,file = "Chp4_pop_dat8.RData")


# logit data generation
pop_dat <- gen_pop_unitMix_logit(alpha = c(0.3,-2.0),gamma = c(-1.05,0.2),beta=c(0.2,0.4,-0.5),nu = -1.2)
save(pop_dat,file = "Chp4_pop_dat2_logit.RData")
pop_dat <- gen_pop_unitMix_logit(alpha = c(0.3,-2.0),gamma = c(-1.05,0.2),beta=c(0.2,2,-2.0),nu = -1.2)
save(pop_dat,file = "Chp4_pop_dat8_logit.RData")

pop_dat <- gen_pop_unitMix_logit(alpha = c(0.3,-0.5),gamma = c(-1.05,0.2),beta=c(0.2,0.4,-0.5),nu = -1.2)
save(pop_dat,file = "Chp4_pop_dat1_logit.RData")
pop_dat <- gen_pop_unitMix_logit(alpha = c(0.3,-0.5),gamma = c(-1.05,0.2),beta=c(0.2,0.4,-2.0),nu = -1.2)
save(pop_dat,file = "Chp4_pop_dat3_logit.RData")
pop_dat <- gen_pop_unitMix_logit(alpha = c(0.3,-2.0),gamma = c(-1.05,0.2),beta=c(0.2,0.4,-2.0),nu = -1.2)
save(pop_dat,file = "Chp4_pop_dat4_logit.RData")

pop_dat <- gen_pop_unitMix_logit(alpha = c(0.3,-0.5),gamma = c(-1.05,0.2),beta=c(0.2,2,-0.5),nu = -1.2)
save(pop_dat,file = "Chp4_pop_dat5_logit.RData")
pop_dat <- gen_pop_unitMix_logit(alpha = c(0.3,-2.0),gamma = c(-1.05,0.2),beta=c(0.2,2,-0.5),nu = -1.2)
save(pop_dat,file = "Chp4_pop_dat6_logit.RData")
pop_dat <- gen_pop_unitMix_logit(alpha = c(0.3,-0.5),gamma = c(-1.05,0.2),beta=c(0.2,2,-2.0),nu = -1.2)
save(pop_dat,file = "Chp4_pop_dat7_logit.RData")

# library(LaplacesDemon)
# ESS(GAMMA[,1])