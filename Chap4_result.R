
# 09/12
# Process Chp4/pop1
n_sim <- sim_n <- 1000
ANWC_mean <- ANWC_var <- MAR_mean <- MAR_var <- matrix(NA,sim_n,length(resultListANWC$alpha_mean) + length(resultListANWC$beta_mean)+
                                  length(resultListANWC$gamma_mean))
ANWC_Tx <- ANWC_Tx_var <- ANWC_Ty <- ANWC_Ty_var <- rep(NA,sim_n)
MAR_Tx <- MAR_Tx_var <- MAR_Ty <- MAR_Ty_var <- rep(NA,sim_n)
alpha_premiss <- matrix(NA,sim_n,length(resultListANWC$alpha_mean))
beta_premiss <- matrix(NA,sim_n,length(resultListANWC$beta_mean))
gamma_premiss <- matrix(NA,sim_n,length(resultListANWC$gamma_mean))

ANWC_b_mat <- matrix(NA,sim_n,length(resultListANWC$mar_varb))
MAR_b_mat <- matrix(NA,sim_n,length(resultListMAR$mar_varb))

ANWC_cond_pt <- MAR_cond_pt <- matrix(NA,sim_n,8)
ANWC_cond_var <- MAR_cond_var <- matrix(NA,sim_n,8)
ANWC_varb <- MAR_varb <-  matrix(NA,sim_n,8)

condY0_prop <- condY1_prop <- rep(NA,sim_n)
condX1_prop <- condX0_prop <- rep(NA,sim_n)

ans_X0 <- ans_X1 <- ans_Y0 <- ans_Y1 <- matrix(NA, sim_n, 2)
ul_X0 <- ul_X1 <- ul_Y0 <- ul_Y1 <- matrix(NA, sim_n, 2)

Tx_premiss <- Ty_premiss <- rep(NA,sim_n)

for (i in 1:n_sim){
  print(i)
  load(paste("./Chp4/pop2_designW_logit/Mis_",i,".RData",sep=""))
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean,resultListANWC$beta_mean, resultListANWC$gamma_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$beta_var, resultListANWC$gamma_var)
  ANWC_Tx[i] <- resultListANWC$total_mean_X
  ANWC_Tx_var[i] <- resultListANWC$total_varX
  ANWC_Ty[i] <- resultListANWC$total_mean_Y
  ANWC_Ty_var[i] <- resultListANWC$total_varY
  
  MAR_mean[i,] <- c(resultListMAR$alpha_mean,resultListMAR$beta_mean, resultListMAR$gamma_mean)
  MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$beta_var, resultListMAR$gamma_var)
  MAR_Tx[i] <- resultListMAR$total_mean_X
  MAR_Tx_var[i] <- resultListMAR$total_varX
  MAR_Ty[i] <- resultListMAR$total_mean_Y
  MAR_Ty_var[i] <- resultListMAR$total_varY
  # 
  ANWC_cond_pt[i,] <- resultListANWC$cond_pt
  MAR_cond_pt[i,] <- resultListMAR$cond_pt
  ANWC_cond_var[i,] <- resultListANWC$cond_var
  MAR_cond_var[i,] <- resultListMAR$cond_var
  ANWC_varb[i,] <- resultListANWC$var_b
  MAR_varb[i,] <- resultListMAR$var_b
  ANWC_b_mat[i,] <- resultListANWC$mar_varb
  MAR_b_mat[i,] <- resultListMAR$mar_varb
  
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y,sub_dat$U,sub_dat$Rx),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y","U","Rx")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$pi <- 1/test_dat$W
  test_dat$Y <- as.factor(test_dat$Y)
  test_dat$U <- as.factor(test_dat$U)
  test_dat$Rx <- as.factor(test_dat$Rx)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m2 <- svyglm(X1~Y+U,design = mydesign,family = quasibinomial(link = "probit"))
  m4 <- svyglm(Y~U,design = mydesign,family = quasibinomial(link = "probit"))
  m6 <- glm(Rx~Y,data = test_dat[which(test_dat$U==0),], family = binomial(link = "probit"))
  alpha_premiss[i,] <- coef(m4)
  beta_premiss[i,] <- coef(m2)
  gamma_premiss[i,] <- coef(m6)
  
  sub_dat_Y1 <- sub_dat[which(sub_dat$Y == 1),]
  condY1_prop[i] <- sum(sub_dat_Y1$X1*sub_dat_Y1$W)/sum(sub_dat_Y1$W)
  
  sub_dat_Y0 <- sub_dat[which(sub_dat$Y == 0),]
  condY0_prop[i] <- sum(sub_dat_Y0$X1*sub_dat_Y0$W)/sum(sub_dat_Y0$W)
  
  sub_dat_X0 <- sub_dat[which(sub_dat$X1 == 0),]
  condX0_prop[i] <- sum(sub_dat_X0$Y*sub_dat_X0$W)/sum(sub_dat_X0$W)
  
  sub_dat_X1 <- sub_dat[which(sub_dat$X1 == 1),]
  condX1_prop[i] <- sum(sub_dat_X1$Y*sub_dat_X1$W)/sum(sub_dat_X1$W)
  
  ##########
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
  
  Tx_premiss[i] <- sum(sub_dat$X1*sub_dat$W)
  Ty_premiss[i] <- sum(sub_dat$Y*sub_dat$W)
}

# ANWC result
# coverage of CI
unit_ANWC1 <- data.frame(ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_Tx, ANWC_TV = ANWC_Tx_var,ANWC_Y = ANWC_Ty,ANWC_Y_var = ANWC_Ty_var)
1-sum(unit_ANWC1$ANWC_mean.1 - 2.009*sqrt(unit_ANWC1$ANWC_var.1) > colMeans(alpha_premiss)[1]|unit_ANWC1$ANWC_mean.1 + 2.009*sqrt(unit_ANWC1$ANWC_var.1) < colMeans(alpha_premiss)[1])/sim_n
1-sum(unit_ANWC1$ANWC_mean.2 - 2.009*sqrt(unit_ANWC1$ANWC_var.2) > colMeans(alpha_premiss)[2]|unit_ANWC1$ANWC_mean.2 + 2.009*sqrt(unit_ANWC1$ANWC_var.2) < colMeans(alpha_premiss)[2])/sim_n

1-sum(unit_ANWC1$ANWC_mean.3 - 2.009*sqrt(unit_ANWC1$ANWC_var.3) > colMeans(beta_premiss)[1]|unit_ANWC1$ANWC_mean.3 + 2.009*sqrt(unit_ANWC1$ANWC_var.3) < colMeans(beta_premiss)[1])/sim_n
1-sum(unit_ANWC1$ANWC_mean.4 - 2.009*sqrt(unit_ANWC1$ANWC_var.4) > colMeans(beta_premiss)[2]|unit_ANWC1$ANWC_mean.4 + 2.009*sqrt(unit_ANWC1$ANWC_var.4) < colMeans(beta_premiss)[2])/sim_n
1-sum(unit_ANWC1$ANWC_mean.5 - 2.009*sqrt(unit_ANWC1$ANWC_var.5) > colMeans(beta_premiss)[3]|unit_ANWC1$ANWC_mean.5 + 2.009*sqrt(unit_ANWC1$ANWC_var.5) < colMeans(beta_premiss)[3])/sim_n

1-sum(unit_ANWC1$ANWC_mean.6 - 2.009*sqrt(unit_ANWC1$ANWC_var.6) > colMeans(gamma_premiss)[1]|unit_ANWC1$ANWC_mean.6 + 2.009*sqrt(unit_ANWC1$ANWC_var.6) < colMeans(gamma_premiss)[1])/sim_n
1-sum(unit_ANWC1$ANWC_mean.7 - 2.009*sqrt(unit_ANWC1$ANWC_var.7) > colMeans(gamma_premiss)[2]|unit_ANWC1$ANWC_mean.7 + 2.009*sqrt(unit_ANWC1$ANWC_var.7) < colMeans(gamma_premiss)[2])/sim_n

1-sum(unit_ANWC1$ANWC_T - 2.009*sqrt(unit_ANWC1$ANWC_TV) >  sum(pop_dat$X1) | unit_ANWC1$ANWC_T + 2.009*sqrt(unit_ANWC1$ANWC_TV) <  sum(pop_dat$X1))/sim_n
1-sum(unit_ANWC1$ANWC_Y - 2.009*sqrt(unit_ANWC1$ANWC_Y_var) >  sum(pop_dat$Y) | unit_ANWC1$ANWC_Y + 2.009*sqrt(unit_ANWC1$ANWC_Y_var) <  sum(pop_dat$Y))/sim_n

colMeans(alpha_premiss)
colMeans(beta_premiss)
colMeans(gamma_premiss)
sum(pop_dat$Y)
sum(pop_dat$X1)

# colmeans
colMeans(ANWC_mean)
mean(ANWC_Ty)
mean(ANWC_Tx)

# pre-miss sd
apply(alpha_premiss,2,sd)
apply(beta_premiss,2,sd)
apply(gamma_premiss,2,sd)
sd(Ty_premiss)
sd(Tx_premiss)

# sd of mean
apply(ANWC_mean,2,sd)
sd(unit_ANWC1$ANWC_Y)
sd(unit_ANWC1$ANWC_T)

# mean of sd
sqrt(colMeans(ANWC_var))
sqrt(mean(unit_ANWC1$ANWC_Y_var))
sqrt(mean(unit_ANWC1$ANWC_TV))


# sqrt(b/50)
sqrt(colMeans(ANWC_b_mat)/50)

# MAR result
# coverage of CI
unit_MAR <- data.frame(MAR_mean = MAR_mean, MAR_var = MAR_var, MAR_T = MAR_Tx, MAR_TV = MAR_Tx_var,MAR_Y = MAR_Ty,MAR_Y_var = MAR_Ty_var)
1-sum(unit_MAR$MAR_mean.1 - 2.009*sqrt(unit_MAR$MAR_var.1) > colMeans(alpha_premiss)[1]|unit_MAR$MAR_mean.1 + 2.009*sqrt(unit_MAR$MAR_var.1) < colMeans(alpha_premiss)[1])/sim_n
1-sum(unit_MAR$MAR_mean.2 - 2.009*sqrt(unit_MAR$MAR_var.2) > colMeans(alpha_premiss)[2]|unit_MAR$MAR_mean.2 + 2.009*sqrt(unit_MAR$MAR_var.2) < colMeans(alpha_premiss)[2])/sim_n

1-sum(unit_MAR$MAR_mean.3 - 2.009*sqrt(unit_MAR$MAR_var.3) > colMeans(beta_premiss)[1]|unit_MAR$MAR_mean.3 + 2.009*sqrt(unit_MAR$MAR_var.3) < colMeans(beta_premiss)[1])/sim_n
1-sum(unit_MAR$MAR_mean.4 - 2.009*sqrt(unit_MAR$MAR_var.4) > colMeans(beta_premiss)[2]|unit_MAR$MAR_mean.4 + 2.009*sqrt(unit_MAR$MAR_var.4) < colMeans(beta_premiss)[2])/sim_n
1-sum(unit_MAR$MAR_mean.5 - 2.009*sqrt(unit_MAR$MAR_var.5) > colMeans(beta_premiss)[3]|unit_MAR$MAR_mean.5 + 2.009*sqrt(unit_MAR$MAR_var.5) < colMeans(beta_premiss)[3])/sim_n

1-sum(unit_MAR$MAR_mean.6 - 2.009*sqrt(unit_MAR$MAR_var.6) > colMeans(gamma_premiss)[1]|unit_MAR$MAR_mean.6 + 2.009*sqrt(unit_MAR$MAR_var.6) < colMeans(gamma_premiss)[1])/sim_n
1-sum(unit_MAR$MAR_mean.7 - 2.009*sqrt(unit_MAR$MAR_var.7) > colMeans(gamma_premiss)[2]|unit_MAR$MAR_mean.7 + 2.009*sqrt(unit_MAR$MAR_var.7) < colMeans(gamma_premiss)[2])/sim_n

1-sum(unit_MAR$MAR_T - 2.009*sqrt(unit_MAR$MAR_TV) >  sum(pop_dat$X1) | unit_MAR$MAR_T + 2.009*sqrt(unit_MAR$MAR_TV) <  sum(pop_dat$X1))/sim_n
1-sum(unit_MAR$MAR_Y - 2.009*sqrt(unit_MAR$MAR_Y_var) >  sum(pop_dat$Y) | unit_MAR$MAR_Y + 2.009*sqrt(unit_MAR$MAR_Y_var) <  sum(pop_dat$Y))/sim_n

# colmeans
colMeans(MAR_mean)
mean(MAR_Ty)
mean(MAR_Tx)

# sd of mean
apply(MAR_mean,2,sd)
sd(unit_MAR$MAR_Y)
sd(unit_MAR$MAR_T)

# mean of sd
sqrt(colMeans(MAR_var))
sqrt(mean(unit_MAR$MAR_Y_var))
sqrt(mean(unit_MAR$MAR_TV))


# sqrt(b/50)
sqrt(colMeans(MAR_b_mat)/50)

#######
# conditional margins
colMeans(ANWC_cond_pt) #X0,X1,Y0,Y1
colMeans(MAR_cond_pt)

# truth
c(colMeans(ans_X0),colMeans(ans_X1),colMeans(ans_Y0),colMeans(ans_Y1))

# pre-miss sd
sqrt(colMeans(ul_X0))
sqrt(colMeans(ul_X1))
sqrt(colMeans(ul_Y0))
sqrt(colMeans(ul_Y1))

# estimate from ANWC
colMeans(ANWC_cond_pt)[c(1,3,5,7)]

# CI coverage 
1-sum(ANWC_cond_pt[,1] - 2.009*sqrt(ANWC_cond_var[,1]) > 
        colMeans(ans_X0)[1]|ANWC_cond_pt[,1] + 2.009*sqrt(ANWC_cond_var[,1]) 
< colMeans(ans_X0)[1])/sim_n
1-sum(ANWC_cond_pt[,3] - 2.009*sqrt(ANWC_cond_var[,3]) > 
        colMeans(ans_X1)[1]|ANWC_cond_pt[,3] + 2.009*sqrt(ANWC_cond_var[,3]) 
      < colMeans(ans_X1)[1])/sim_n
1-sum(ANWC_cond_pt[,5] - 2.009*sqrt(ANWC_cond_var[,5]) > 
        colMeans(ans_Y0)[1]|ANWC_cond_pt[,5] + 2.009*sqrt(ANWC_cond_var[,5]) 
      < colMeans(ans_Y0)[1])/sim_n
1-sum(ANWC_cond_pt[,7] - 2.009*sqrt(ANWC_cond_var[,7]) > 
        colMeans(ans_Y1)[1]|ANWC_cond_pt[,7] + 2.009*sqrt(ANWC_cond_var[,7]) 
      < colMeans(ans_Y1)[1])/sim_n

# sd
apply(ANWC_cond_pt,2,sd)

# avg. est. sd
sqrt(colMeans(ANWC_cond_var))

# var_b
sqrt(colMeans(ANWC_varb)/50)

# estimate from MAR
colMeans(MAR_cond_pt)[c(1,3,5,7)]

# CI coverage 
1-sum(MAR_cond_pt[,1] - 2.009*sqrt(MAR_cond_var[,1]) > 
        colMeans(ans_X0)[1]|MAR_cond_pt[,1] + 2.009*sqrt(MAR_cond_var[,1]) 
      < colMeans(ans_X0)[1])/sim_n
1-sum(MAR_cond_pt[,3] - 2.009*sqrt(MAR_cond_var[,3]) > 
        colMeans(ans_X1)[1]|MAR_cond_pt[,3] + 2.009*sqrt(MAR_cond_var[,3]) 
      < colMeans(ans_X1)[1])/sim_n
1-sum(MAR_cond_pt[,5] - 2.009*sqrt(MAR_cond_var[,5]) > 
        colMeans(ans_Y0)[1]|MAR_cond_pt[,5] + 2.009*sqrt(MAR_cond_var[,5]) 
      < colMeans(ans_Y0)[1])/sim_n
1-sum(MAR_cond_pt[,7] - 2.009*sqrt(MAR_cond_var[,7]) > 
        colMeans(ans_Y1)[1]|MAR_cond_pt[,7] + 2.009*sqrt(MAR_cond_var[,7]) 
      < colMeans(ans_Y1)[1])/sim_n

# sd
apply(MAR_cond_pt,2,sd)

# avg. est. sd
sqrt(colMeans(MAR_cond_var))

# var_b
sqrt(colMeans(MAR_varb)/50)

# 09/13/2021
# get u and b for each simulation
u_mat <- matrix(NA,n_sim,9)
b_mat <- matrix(NA,n_sim,9)
Tx_premiss <- Ty_premiss <- rep(NA,sim_n)

for (i in 1:n_sim){
  print(i)
  load(paste("./Chp4/pop1/Mis_",i,".RData",sep=""))
  X_input <- MI_dataANWC[,which(grepl("X", names(MI_dataANWC), fixed = TRUE))]
  Y_input <- MI_dataANWC[,which(grepl("Y", names(MI_dataANWC), fixed = TRUE))]
  
  ans <- matrix(NA,50,2)
  ul <- matrix(NA,50,2)
  total_varX <- total_varY <- rep(NA,50)
  
  # calculate beta
  ans_b <- matrix(NA,50,3)
  ul_b <- matrix(NA,50,3)
  
  # calculate gamma using stats package
  ans_g <- matrix(NA,50,2)
  ul_g <- matrix(NA,50,2)
  
  total_X <- total_Y <- rep(NA,50)
  
  for (j in 1:50){
    test_dat <- as.data.frame(cbind(t(X_input[j,]),rownames(sub_dat),sub_dat$W,t(Y_input[j,]),sub_dat$U,sub_dat$Rx),stringsAsFactors = FALSE)
    names(test_dat) <- c("X1","id","W","Y","U","Rx")
    total_varX[j] <- sum((as.numeric(test_dat$X1)*sub_dat$W)^2*(1-1/sub_dat$W))
    total_varY[j] <- sum((as.numeric(test_dat$Y)*sub_dat$W)^2*(1-1/sub_dat$W))
    test_dat$X1 <- as.factor(test_dat$X1)
    test_dat$id <- as.numeric(test_dat$id)
    test_dat$W <- as.numeric(test_dat$W)
    test_dat$pi <- 1/test_dat$W
    test_dat$Y <- as.factor(test_dat$Y)
    test_dat$U <- as.factor(test_dat$U)
    test_dat$Rx <- as.factor(test_dat$Rx)
    mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
    
    m1 <- svyglm(X1~Y+U,design = mydesign,family = quasibinomial(link = "probit"))
    ans_b[j,] <- coef(m1)
    ul_b[j,] <- diag(vcov(m1))
    m3 <- svyglm(Y~U,design = mydesign,family = quasibinomial(link = "probit"))
    ans[j,] <- coef(m3)
    ul[j,] <- diag(vcov(m3))
    m2 <- glm(Rx~Y,data = test_dat, family = binomial(link = "probit"))
    ans_g[j,] <- coef(m2)
    ul_g[j,] <- diag(vcov(m2))
    
    total_X[j] <- sum(X_input[j,]*sub_dat$W)
    total_Y[j] <- sum(Y_input[j,]*sub_dat$W)
  }
  
  Tx_premiss[i] <- sum(sub_dat$X1*sub_dat$W)
  Ty_premiss[i] <- sum(sub_dat$Y*sub_dat$W)
  # colMeans(ans)
  u_L_a <- colMeans(ul)
  b_L_a <- apply((scale(ans,scale=FALSE))^2/(50-1),2,sum)
  T_L_a <- (1+1/50)*b_L_a + u_L_a
  
  # colMeans(ans_g)
  u_L_g <- colMeans(ul_g)
  b_L_g <- apply((scale(ans_g,scale=FALSE))^2/(50-1),2,sum)
  T_L_g <- (1+1/50)*b_L_g + u_L_g
  
  # colMeans(ans_b)
  u_L_b <- colMeans(ul_b)
  b_L_b <- apply((scale(ans_b,scale=FALSE))^2/(50-1),2,sum)
  T_L_b <- (1+1/50)*b_L_b + u_L_b
  
  u_L_t <- mean(total_varX)
  b_L_t <- apply((scale(total_X,scale=FALSE))^2/(50-1),2,sum)
  T_L_t <- (1+1/50)*b_L_t + u_L_t
  
  u_L_tY <- mean(total_varY)
  b_L_tY <- apply((scale(total_Y,scale=FALSE))^2/(50-1),2,sum)
  T_L_tY <- (1+1/50)*b_L_tY + u_L_tY
  
  u_mat[i,] <- c(u_L_a,u_L_b,u_L_g,u_L_t,u_L_tY)
  b_mat[i,] <- c(b_L_a,b_L_b,b_L_g,b_L_t,b_L_tY)
}

sd(Tx_premiss)
sd(Ty_premiss)



####################
X_input <- MI_dataANWC[,which(grepl("X", names(MI_dataANWC), fixed = TRUE))]
Y_input <- MI_dataANWC[,which(grepl("Y", names(MI_dataANWC), fixed = TRUE))]
test_dat <- as.data.frame(cbind(t(X_input[j,]),rownames(sub_dat),sub_dat$W,t(Y_input[j,]),sub_dat$U,sub_dat$Rx),stringsAsFactors = FALSE)
names(test_dat) <- c("X1","id","W","Y","U","Rx")
#total_varX[j] <- sum((as.numeric(test_dat$X1)*sub_dat$W)^2*(1-1/sub_dat$W))
#total_varY[j] <- sum((as.numeric(test_dat$Y)*sub_dat$W)^2*(1-1/sub_dat$W))
test_dat$X1 <- as.factor(test_dat$X1)
test_dat$id <- as.numeric(test_dat$id)
test_dat$W <- as.numeric(test_dat$W)
test_dat$pi <- 1/test_dat$W
test_dat$Y <- as.factor(test_dat$Y)
test_dat$U <- as.factor(test_dat$U)
test_dat$Rx <- as.factor(test_dat$Rx)

# Tx_premiss_origin <- Ty_premiss_origin <- rep(NA,n_sim)
ans_X0_origin <- ans_X1_origin <- ans_Y0_origin <- ans_Y1_origin <- matrix(NA, sim_n, 2)
ul_X0_origin <- ul_X1_origin <- ul_Y0_origin <- ul_Y1_origin <- matrix(NA, sim_n, 2)

alpha_premiss_origin <- matrix(NA,sim_n,length(resultListANWC$alpha_mean))
beta_premiss_origin <- matrix(NA,sim_n,length(resultListANWC$beta_mean))
gamma_premiss_origin <- matrix(NA,sim_n,length(resultListANWC$gamma_mean))

for (i in 1:n_sim){
  print(i)
  load(paste("./Chp4/pop3/Mis_",i,".RData",sep=""))
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$origin_W,sub_dat$Y,sub_dat$U,sub_dat$Rx),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y","U","Rx")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)/100
  # test_dat$pi <- 1/test_dat$W
  test_dat$Y <- as.factor(test_dat$Y)
  test_dat$U <- as.factor(test_dat$U)
  test_dat$Rx <- as.factor(test_dat$Rx)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  # Tx_premiss_origin[i] <- sum(sub_dat$X1*sub_dat$origin_W/100)
  # Ty_premiss_origin[i] <- sum(sub_dat$Y*sub_dat$origin_W/100)
  
  m2 <- svyglm(X1~Y+U,design = mydesign,family = quasibinomial(link = "probit"))
  m4 <- svyglm(Y~U,design = mydesign,family = quasibinomial(link = "probit"))
  m6 <- glm(Rx~Y,data = test_dat[which(test_dat$U==0),], family = binomial(link = "probit"))
  alpha_premiss_origin[i,] <- coef(m4)
  beta_premiss_origin[i,] <- coef(m2)
  gamma_premiss_origin[i,] <- coef(m6)
  
  # wgt_tblY <- svyby(~as.factor(X1),~as.factor(Y),mydesign,svymean)
  # wgt_dfY <- data.frame(wgt_tblY)
  # ans_Y0_origin[i,] <- c(t(wgt_dfY[1,2:3]))
  # ul_Y0_origin[i,] <- c(t(wgt_dfY[1,4:5]))
  # ans_Y1_origin[i,] <- c(t(wgt_dfY[2,2:3]))
  # ul_Y1_origin[i,] <- c(t(wgt_dfY[2,4:5]))
  # 
  # wgt_tblX <- svyby(~as.factor(Y),~as.factor(X1),mydesign,svymean)
  # wgt_dfX <- data.frame(wgt_tblX)
  # ans_X0_origin[i,] <- c(t(wgt_dfX[1,2:3]))
  # ul_X0_origin[i,] <- c(t(wgt_dfX[1,4:5]))
  # ans_X1_origin[i,] <- c(t(wgt_dfX[2,2:3]))
  # ul_X1_origin[i,] <- c(t(wgt_dfX[2,4:5]))
}

pop_dat_Y1 <- pop_dat[which(pop_dat$Y == 1),]
ans_Y1_true <- 1-mean(pop_dat_Y1$X1)

pop_dat_Y0 <- pop_dat[which(pop_dat$Y == 0),]
ans_Y0_true <- 1-mean(pop_dat_Y0$X1)

pop_dat_X0 <- pop_dat[which(pop_dat$X1 == 0),]
ans_X0_true <- 1-mean(pop_dat_X0$Y)

pop_dat_X1 <- pop_dat[which(pop_dat$X1 == 1),]
ans_X1_true <- 1-mean(pop_dat_X1$Y)

# CI coverage for ANWC of conditional proportions, based on truth
1-sum(ANWC_cond_pt[,1] - 2.009*sqrt(ANWC_cond_var[,1]) > 
        ans_X0_true|ANWC_cond_pt[,1] + 2.009*sqrt(ANWC_cond_var[,1]) 
      < ans_X0_true)/sim_n
1-sum(ANWC_cond_pt[,3] - 2.009*sqrt(ANWC_cond_var[,3]) > 
        ans_X1_true|ANWC_cond_pt[,3] + 2.009*sqrt(ANWC_cond_var[,3]) 
      < ans_X1_true)/sim_n
1-sum(ANWC_cond_pt[,5] - 2.009*sqrt(ANWC_cond_var[,5]) > 
        ans_Y0_true|ANWC_cond_pt[,5] + 2.009*sqrt(ANWC_cond_var[,5]) 
      < ans_Y0_true)/sim_n
1-sum(ANWC_cond_pt[,7] - 2.009*sqrt(ANWC_cond_var[,7]) > 
        ans_Y1_true|ANWC_cond_pt[,7] + 2.009*sqrt(ANWC_cond_var[,7]) 
      < ans_Y1_true)/sim_n

# CI coverage for MAR of conditional proportions, based on truth
1-sum(MAR_cond_pt[,1] - 2.009*sqrt(MAR_cond_var[,1]) > 
        ans_X0_true|MAR_cond_pt[,1] + 2.009*sqrt(MAR_cond_var[,1]) 
      < ans_X0_true)/sim_n
1-sum(MAR_cond_pt[,3] - 2.009*sqrt(MAR_cond_var[,3]) > 
        ans_X1_true|MAR_cond_pt[,3] + 2.009*sqrt(MAR_cond_var[,3]) 
      < ans_X1_true)/sim_n
1-sum(MAR_cond_pt[,5] - 2.009*sqrt(MAR_cond_var[,5]) > 
        ans_Y0_true|MAR_cond_pt[,5] + 2.009*sqrt(MAR_cond_var[,5]) 
      < ans_Y0_true)/sim_n
1-sum(MAR_cond_pt[,7] - 2.009*sqrt(MAR_cond_var[,7]) > 
        ans_Y1_true|MAR_cond_pt[,7] + 2.009*sqrt(MAR_cond_var[,7]) 
      < ans_Y1_true)/sim_n

# CI coverage of ANWC prop (original weights)
1-sum(ANWC_cond_pt[,1] - 2.009*sqrt(ANWC_cond_var[,1]) > 
        colMeans(ans_X0_origin)[1]|ANWC_cond_pt[,1] + 2.009*sqrt(ANWC_cond_var[,1]) 
      < colMeans(ans_X0_origin)[1])/sim_n
1-sum(ANWC_cond_pt[,3] - 2.009*sqrt(ANWC_cond_var[,3]) > 
        colMeans(ans_X1_origin)[1]|ANWC_cond_pt[,3] + 2.009*sqrt(ANWC_cond_var[,3]) 
      < colMeans(ans_X1_origin)[1])/sim_n
1-sum(ANWC_cond_pt[,5] - 2.009*sqrt(ANWC_cond_var[,5]) > 
        colMeans(ans_Y0_origin)[1]|ANWC_cond_pt[,5] + 2.009*sqrt(ANWC_cond_var[,5]) 
      < colMeans(ans_Y0_origin)[1])/sim_n
1-sum(ANWC_cond_pt[,7] - 2.009*sqrt(ANWC_cond_var[,7]) > 
        colMeans(ans_Y1_origin)[1]|ANWC_cond_pt[,7] + 2.009*sqrt(ANWC_cond_var[,7]) 
      < colMeans(ans_Y1_origin)[1])/sim_n

# CI coverage of MAR prop (original weights)
1-sum(MAR_cond_pt[,1] - 2.009*sqrt(MAR_cond_var[,1]) > 
        colMeans(ans_X0_origin)[1]|MAR_cond_pt[,1] + 2.009*sqrt(MAR_cond_var[,1]) 
      < colMeans(ans_X0_origin)[1])/sim_n
1-sum(MAR_cond_pt[,3] - 2.009*sqrt(MAR_cond_var[,3]) > 
        colMeans(ans_X1_origin)[1]|MAR_cond_pt[,3] + 2.009*sqrt(MAR_cond_var[,3]) 
      < colMeans(ans_X1_origin)[1])/sim_n
1-sum(MAR_cond_pt[,5] - 2.009*sqrt(MAR_cond_var[,5]) > 
        colMeans(ans_Y0_origin)[1]|MAR_cond_pt[,5] + 2.009*sqrt(MAR_cond_var[,5]) 
      < colMeans(ans_Y0_origin)[1])/sim_n
1-sum(MAR_cond_pt[,7] - 2.009*sqrt(MAR_cond_var[,7]) > 
        colMeans(ans_Y1_origin)[1]|MAR_cond_pt[,7] + 2.009*sqrt(MAR_cond_var[,7]) 
      < colMeans(ans_Y1_origin)[1])/sim_n

# proportion truth
ans_X0_true
ans_X1_true
ans_Y0_true
ans_Y1_true

# premiss sd based on samples with original weights
apply(ans_X0_origin,2,sd)
apply(ans_X1_origin,2,sd)
apply(ans_Y0_origin,2,sd)
apply(ans_Y1_origin,2,sd)

# CI coverage for ANWC of conditional proportions, based on truth
1-sum(ANWC_cond_pt[,1] - 2.009*sqrt(ANWC_cond_var[,1]) > 
        ans_X0_true|ANWC_cond_pt[,1] + 2.009*sqrt(ANWC_cond_var[,1]) 
      < ans_X0_true)/sim_n
1-sum(ANWC_cond_pt[,3] - 2.009*sqrt(ANWC_cond_var[,3]) > 
        ans_X1_true|ANWC_cond_pt[,3] + 2.009*sqrt(ANWC_cond_var[,3]) 
      < ans_X1_true)/sim_n
1-sum(ANWC_cond_pt[,5] - 2.009*sqrt(ANWC_cond_var[,5]) > 
        ans_Y0_true|ANWC_cond_pt[,5] + 2.009*sqrt(ANWC_cond_var[,5]) 
      < ans_Y0_true)/sim_n
1-sum(ANWC_cond_pt[,7] - 2.009*sqrt(ANWC_cond_var[,7]) > 
        ans_Y1_true|ANWC_cond_pt[,7] + 2.009*sqrt(ANWC_cond_var[,7]) 
      < ans_Y1_true)/sim_n


# CI coverage for MAR of conditional proportions, based on truth
1-sum(MAR_cond_pt[,1] - 2.009*sqrt(MAR_cond_var[,1]) > 
        ans_X0_true|MAR_cond_pt[,1] + 2.009*sqrt(MAR_cond_var[,1]) 
      < ans_X0_true)/sim_n
1-sum(MAR_cond_pt[,3] - 2.009*sqrt(MAR_cond_var[,3]) > 
        ans_X1_true|MAR_cond_pt[,3] + 2.009*sqrt(MAR_cond_var[,3]) 
      < ans_X1_true)/sim_n
1-sum(MAR_cond_pt[,5] - 2.009*sqrt(MAR_cond_var[,5]) > 
        ans_Y0_true|MAR_cond_pt[,5] + 2.009*sqrt(MAR_cond_var[,5]) 
      < ans_Y0_true)/sim_n
1-sum(MAR_cond_pt[,7] - 2.009*sqrt(MAR_cond_var[,7]) > 
        ans_Y1_true|MAR_cond_pt[,7] + 2.009*sqrt(MAR_cond_var[,7]) 
      < ans_Y1_true)/sim_n

# sd
apply(ANWC_cond_pt,2,sd)

# avg. est. sd
sqrt(colMeans(ANWC_cond_var))

# var_b
sqrt(colMeans(ANWC_varb)/50)

# sd
apply(MAR_cond_pt,2,sd)

# avg. est. sd
sqrt(colMeans(MAR_cond_var))

# var_b
sqrt(colMeans(MAR_varb)/50)


#############################################
########## 09/20 ############################
#############################################
# Summary of simulation results #
n_sim <- sim_n <- 1000
ANWC_mean <- ANWC_var <- MAR_mean <- MAR_var <- matrix(NA,sim_n,length(resultListANWC$alpha_mean) + length(resultListANWC$beta_mean)+
                                                         length(resultListANWC$gamma_mean))
ANWC_Tx <- ANWC_Tx_var <- ANWC_Ty <- ANWC_Ty_var <- rep(NA,sim_n)
MAR_Tx <- MAR_Tx_var <- MAR_Ty <- MAR_Ty_var <- rep(NA,sim_n)
alpha_premiss <- matrix(NA,sim_n,length(resultListANWC$alpha_mean))
beta_premiss <- matrix(NA,sim_n,length(resultListANWC$beta_mean))
gamma_premiss <- matrix(NA,sim_n,length(resultListANWC$gamma_mean))

ANWC_b_mat <- matrix(NA,sim_n,length(resultListANWC$mar_varb))
MAR_b_mat <- matrix(NA,sim_n,length(resultListMAR$mar_varb))

ANWC_cond_pt <- MAR_cond_pt <- matrix(NA,sim_n,8)
ANWC_cond_var <- MAR_cond_var <- matrix(NA,sim_n,8)
ANWC_varb <- MAR_varb <-  matrix(NA,sim_n,8)

Tx_premiss_origin <- Ty_premiss_origin <- rep(NA,sim_n)

ans_X0_origin <- ans_X1_origin <- ans_Y0_origin <- ans_Y1_origin <- matrix(NA, sim_n, 2)
ul_X0_origin <- ul_X1_origin <- ul_Y0_origin <- ul_Y1_origin <- matrix(NA, sim_n, 2)

for (i in 1:n_sim){
  print(i)
  load(paste("./Chp4/pop8_designW_logit/Mis_",i,".RData",sep=""))
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean,resultListANWC$beta_mean, resultListANWC$gamma_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$beta_var, resultListANWC$gamma_var)
  ANWC_Tx[i] <- resultListANWC$total_mean_X
  ANWC_Tx_var[i] <- resultListANWC$total_varX
  ANWC_Ty[i] <- resultListANWC$total_mean_Y
  ANWC_Ty_var[i] <- resultListANWC$total_varY

  MAR_mean[i,] <- c(resultListMAR$alpha_mean,resultListMAR$beta_mean, resultListMAR$gamma_mean)
  MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$beta_var, resultListMAR$gamma_var)
  MAR_Tx[i] <- resultListMAR$total_mean_X
  MAR_Tx_var[i] <- resultListMAR$total_varX
  MAR_Ty[i] <- resultListMAR$total_mean_Y
  MAR_Ty_var[i] <- resultListMAR$total_varY

  ANWC_cond_pt[i,] <- resultListANWC$cond_pt
  MAR_cond_pt[i,] <- resultListMAR$cond_pt
  ANWC_cond_var[i,] <- resultListANWC$cond_var
  MAR_cond_var[i,] <- resultListMAR$cond_var
  ANWC_varb[i,] <- resultListANWC$var_b
  MAR_varb[i,] <- resultListMAR$var_b
  ANWC_b_mat[i,] <- resultListANWC$mar_varb
  MAR_b_mat[i,] <- resultListMAR$mar_varb

  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y,sub_dat$U,sub_dat$Rx,sub_dat$origin_W),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y","U","Rx","Origin_W")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$pi <- 1/test_dat$W
  test_dat$Y <- as.factor(test_dat$Y)
  test_dat$U <- as.factor(test_dat$U)
  test_dat$Rx <- as.factor(test_dat$Rx)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m2 <- svyglm(X1~Y+U,design = mydesign,family = quasibinomial(link = "probit"))
  m4 <- svyglm(Y~U,design = mydesign,family = quasibinomial(link = "probit"))
  m6 <- glm(Rx~Y,data = test_dat[which(test_dat$U==0),], family = binomial(link = "probit"))
  alpha_premiss[i,] <- coef(m4)
  beta_premiss[i,] <- coef(m2)
  gamma_premiss[i,] <- coef(m6)
  
  test_dat$W <- as.numeric(test_dat$Origin_W)/100
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  
  wgt_tblY <- svyby(~as.factor(X1),~as.factor(Y),mydesign,svymean)
  wgt_dfY <- data.frame(wgt_tblY)
  ans_Y0_origin[i,] <- c(t(wgt_dfY[1,2:3]))
  ul_Y0_origin[i,] <- c(t(wgt_dfY[1,4:5]))
  ans_Y1_origin[i,] <- c(t(wgt_dfY[2,2:3]))
  ul_Y1_origin[i,] <- c(t(wgt_dfY[2,4:5]))

  wgt_tblX <- svyby(~as.factor(Y),~as.factor(X1),mydesign,svymean)
  wgt_dfX <- data.frame(wgt_tblX)
  ans_X0_origin[i,] <- c(t(wgt_dfX[1,2:3]))
  ul_X0_origin[i,] <- c(t(wgt_dfX[1,4:5]))
  ans_X1_origin[i,] <- c(t(wgt_dfX[2,2:3]))
  ul_X1_origin[i,] <- c(t(wgt_dfX[2,4:5]))
  
  Tx_premiss_origin[i] <- sum(sub_dat$X1*sub_dat$origin_W/100)
  Ty_premiss_origin[i] <- sum(sub_dat$Y*sub_dat$origin_W/100)
}

# truth
colMeans(alpha_premiss)
colMeans(beta_premiss)
colMeans(gamma_premiss)
sum(pop_dat$Y)
sum(pop_dat$X1)

# pre-miss sd
apply(alpha_premiss,2,sd)
apply(beta_premiss,2,sd)
apply(gamma_premiss,2,sd)
sd(Ty_premiss_origin)
sd(Tx_premiss_origin)

# colmeans
colMeans(ANWC_mean)
mean(ANWC_Ty)
mean(ANWC_Tx)

# ANWC result
# coverage of CI
unit_ANWC2 <- data.frame(ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_Tx, ANWC_TV = ANWC_Tx_var,ANWC_Y = ANWC_Ty,ANWC_Y_var = ANWC_Ty_var)
1-sum(unit_ANWC2$ANWC_mean.1 - 2.009*sqrt(unit_ANWC2$ANWC_var.1) > colMeans(alpha_premiss)[1]|unit_ANWC2$ANWC_mean.1 + 2.009*sqrt(unit_ANWC2$ANWC_var.1) < colMeans(alpha_premiss)[1])/sim_n
1-sum(unit_ANWC2$ANWC_mean.2 - 2.009*sqrt(unit_ANWC2$ANWC_var.2) > colMeans(alpha_premiss)[2]|unit_ANWC2$ANWC_mean.2 + 2.009*sqrt(unit_ANWC2$ANWC_var.2) < colMeans(alpha_premiss)[2])/sim_n

1-sum(unit_ANWC2$ANWC_mean.3 - 2.009*sqrt(unit_ANWC2$ANWC_var.3) > colMeans(beta_premiss)[1]|unit_ANWC2$ANWC_mean.3 + 2.009*sqrt(unit_ANWC2$ANWC_var.3) < colMeans(beta_premiss)[1])/sim_n
1-sum(unit_ANWC2$ANWC_mean.4 - 2.009*sqrt(unit_ANWC2$ANWC_var.4) > colMeans(beta_premiss)[2]|unit_ANWC2$ANWC_mean.4 + 2.009*sqrt(unit_ANWC2$ANWC_var.4) < colMeans(beta_premiss)[2])/sim_n
1-sum(unit_ANWC2$ANWC_mean.5 - 2.009*sqrt(unit_ANWC2$ANWC_var.5) > colMeans(beta_premiss)[3]|unit_ANWC2$ANWC_mean.5 + 2.009*sqrt(unit_ANWC2$ANWC_var.5) < colMeans(beta_premiss)[3])/sim_n

1-sum(unit_ANWC2$ANWC_mean.6 - 2.009*sqrt(unit_ANWC2$ANWC_var.6) > colMeans(gamma_premiss)[1]|unit_ANWC2$ANWC_mean.6 + 2.009*sqrt(unit_ANWC2$ANWC_var.6) < colMeans(gamma_premiss)[1])/sim_n
1-sum(unit_ANWC2$ANWC_mean.7 - 2.009*sqrt(unit_ANWC2$ANWC_var.7) > colMeans(gamma_premiss)[2]|unit_ANWC2$ANWC_mean.7 + 2.009*sqrt(unit_ANWC2$ANWC_var.7) < colMeans(gamma_premiss)[2])/sim_n

1-sum(unit_ANWC2$ANWC_Y - 2.009*sqrt(unit_ANWC2$ANWC_Y_var) >  sum(pop_dat$Y) | unit_ANWC2$ANWC_Y + 2.009*sqrt(unit_ANWC2$ANWC_Y_var) <  sum(pop_dat$Y))/sim_n
1-sum(unit_ANWC2$ANWC_T - 2.009*sqrt(unit_ANWC2$ANWC_TV) >  sum(pop_dat$X1) | unit_ANWC2$ANWC_T + 2.009*sqrt(unit_ANWC2$ANWC_TV) <  sum(pop_dat$X1))/sim_n

# sd of mean
apply(ANWC_mean,2,sd)
sd(unit_ANWC2$ANWC_Y)
sd(unit_ANWC2$ANWC_T)

# mean of sd
sqrt(colMeans(ANWC_var))
sqrt(mean(unit_ANWC2$ANWC_Y_var))
sqrt(mean(unit_ANWC2$ANWC_TV))


# sqrt(b/50)
sqrt(colMeans(ANWC_b_mat)/50)

# MAR result
# colmeans
colMeans(MAR_mean)
mean(MAR_Ty)
mean(MAR_Tx)

# coverage of CI
unit_MAR2 <- data.frame(MAR_mean = MAR_mean, MAR_var = MAR_var, MAR_T = MAR_Tx, MAR_TV = MAR_Tx_var,MAR_Y = MAR_Ty,MAR_Y_var = MAR_Ty_var)
1-sum(unit_MAR2$MAR_mean.1 - 2.009*sqrt(unit_MAR2$MAR_var.1) > colMeans(alpha_premiss)[1]|unit_MAR2$MAR_mean.1 + 2.009*sqrt(unit_MAR2$MAR_var.1) < colMeans(alpha_premiss)[1])/sim_n
1-sum(unit_MAR2$MAR_mean.2 - 2.009*sqrt(unit_MAR2$MAR_var.2) > colMeans(alpha_premiss)[2]|unit_MAR2$MAR_mean.2 + 2.009*sqrt(unit_MAR2$MAR_var.2) < colMeans(alpha_premiss)[2])/sim_n

1-sum(unit_MAR2$MAR_mean.3 - 2.009*sqrt(unit_MAR2$MAR_var.3) > colMeans(beta_premiss)[1]|unit_MAR2$MAR_mean.3 + 2.009*sqrt(unit_MAR2$MAR_var.3) < colMeans(beta_premiss)[1])/sim_n
1-sum(unit_MAR2$MAR_mean.4 - 2.009*sqrt(unit_MAR2$MAR_var.4) > colMeans(beta_premiss)[2]|unit_MAR2$MAR_mean.4 + 2.009*sqrt(unit_MAR2$MAR_var.4) < colMeans(beta_premiss)[2])/sim_n
1-sum(unit_MAR2$MAR_mean.5 - 2.009*sqrt(unit_MAR2$MAR_var.5) > colMeans(beta_premiss)[3]|unit_MAR2$MAR_mean.5 + 2.009*sqrt(unit_MAR2$MAR_var.5) < colMeans(beta_premiss)[3])/sim_n

1-sum(unit_MAR2$MAR_mean.6 - 2.009*sqrt(unit_MAR2$MAR_var.6) > colMeans(gamma_premiss)[1]|unit_MAR2$MAR_mean.6 + 2.009*sqrt(unit_MAR2$MAR_var.6) < colMeans(gamma_premiss)[1])/sim_n
1-sum(unit_MAR2$MAR_mean.7 - 2.009*sqrt(unit_MAR2$MAR_var.7) > colMeans(gamma_premiss)[2]|unit_MAR2$MAR_mean.7 + 2.009*sqrt(unit_MAR2$MAR_var.7) < colMeans(gamma_premiss)[2])/sim_n

1-sum(unit_MAR2$MAR_Y - 2.009*sqrt(unit_MAR2$MAR_Y_var) >  sum(pop_dat$Y) | unit_MAR2$MAR_Y + 2.009*sqrt(unit_MAR2$MAR_Y_var) <  sum(pop_dat$Y))/sim_n
1-sum(unit_MAR2$MAR_T - 2.009*sqrt(unit_MAR2$MAR_TV) >  sum(pop_dat$X1) | unit_MAR2$MAR_T + 2.009*sqrt(unit_MAR2$MAR_TV) <  sum(pop_dat$X1))/sim_n


# sd of mean
apply(MAR_mean,2,sd)
sd(unit_MAR2$MAR_Y)
sd(unit_MAR2$MAR_T)

# mean of sd
sqrt(colMeans(MAR_var))
sqrt(mean(unit_MAR2$MAR_Y_var))
sqrt(mean(unit_MAR2$MAR_TV))


# sqrt(b/50)
sqrt(colMeans(MAR_b_mat)/50)

# proportion truth
pop_dat_Y1 <- pop_dat[which(pop_dat$Y == 1),]
ans_Y1_true <- 1-mean(pop_dat_Y1$X1)

pop_dat_Y0 <- pop_dat[which(pop_dat$Y == 0),]
ans_Y0_true <- 1-mean(pop_dat_Y0$X1)

pop_dat_X0 <- pop_dat[which(pop_dat$X1 == 0),]
ans_X0_true <- 1-mean(pop_dat_X0$Y)

pop_dat_X1 <- pop_dat[which(pop_dat$X1 == 1),]
ans_X1_true <- 1-mean(pop_dat_X1$Y)

ans_X0_true
ans_X1_true
ans_Y0_true
ans_Y1_true

# premiss sd based on samples with original weights
apply(ans_X0_origin,2,sd)
apply(ans_X1_origin,2,sd)
apply(ans_Y0_origin,2,sd)
apply(ans_Y1_origin,2,sd)

# estimate
colMeans(ANWC_cond_pt)
colMeans(MAR_cond_pt)

# CI coverage for ANWC of conditional proportions, based on truth
1-sum(ANWC_cond_pt[,1] - 2.009*sqrt(ANWC_cond_var[,1]) > 
        ans_X0_true|ANWC_cond_pt[,1] + 2.009*sqrt(ANWC_cond_var[,1]) 
      < ans_X0_true)/sim_n
1-sum(ANWC_cond_pt[,3] - 2.009*sqrt(ANWC_cond_var[,3]) > 
        ans_X1_true|ANWC_cond_pt[,3] + 2.009*sqrt(ANWC_cond_var[,3]) 
      < ans_X1_true)/sim_n
1-sum(ANWC_cond_pt[,5] - 2.009*sqrt(ANWC_cond_var[,5]) > 
        ans_Y0_true|ANWC_cond_pt[,5] + 2.009*sqrt(ANWC_cond_var[,5]) 
      < ans_Y0_true)/sim_n
1-sum(ANWC_cond_pt[,7] - 2.009*sqrt(ANWC_cond_var[,7]) > 
        ans_Y1_true|ANWC_cond_pt[,7] + 2.009*sqrt(ANWC_cond_var[,7]) 
      < ans_Y1_true)/sim_n


# CI coverage for MAR of conditional proportions, based on truth
1-sum(MAR_cond_pt[,1] - 2.009*sqrt(MAR_cond_var[,1]) > 
        ans_X0_true|MAR_cond_pt[,1] + 2.009*sqrt(MAR_cond_var[,1]) 
      < ans_X0_true)/sim_n
1-sum(MAR_cond_pt[,3] - 2.009*sqrt(MAR_cond_var[,3]) > 
        ans_X1_true|MAR_cond_pt[,3] + 2.009*sqrt(MAR_cond_var[,3]) 
      < ans_X1_true)/sim_n
1-sum(MAR_cond_pt[,5] - 2.009*sqrt(MAR_cond_var[,5]) > 
        ans_Y0_true|MAR_cond_pt[,5] + 2.009*sqrt(MAR_cond_var[,5]) 
      < ans_Y0_true)/sim_n
1-sum(MAR_cond_pt[,7] - 2.009*sqrt(MAR_cond_var[,7]) > 
        ans_Y1_true|MAR_cond_pt[,7] + 2.009*sqrt(MAR_cond_var[,7]) 
      < ans_Y1_true)/sim_n

# sd
apply(ANWC_cond_pt,2,sd)

# avg. est. sd
sqrt(colMeans(ANWC_cond_var))

# var_b
sqrt(colMeans(ANWC_varb)/50)

# sd
apply(MAR_cond_pt,2,sd)

# avg. est. sd
sqrt(colMeans(MAR_cond_var))

# var_b
sqrt(colMeans(MAR_varb)/50)

sample_total_W <- rep(NA,500)
for (i in 1:500){
  load(paste("./Chp4/pop4_designW/Mis_",i,".RData",sep=""))
  sample_total_W[i] <- sum(sub_dat$origin_W/100)
}

#######################################
# 2023/01/05 Update
# Summary of simulation results #
#######################################
{
n_sim <- sim_n <- 1000
ANWC_mean <- ANWC_var <- MAR_mean <- MAR_var <- matrix(NA,sim_n,length(resultListANWC$alpha_mean) + length(resultListANWC$beta_mean)+
                                                         length(resultListANWC$gamma_mean))
ANWC_Tx <- ANWC_Tx_var <- ANWC_Ty <- ANWC_Ty_var <- rep(NA,sim_n)
MAR_Tx <- MAR_Tx_var <- MAR_Ty <- MAR_Ty_var <- rep(NA,sim_n)
alpha_premiss <- matrix(NA,sim_n,length(resultListANWC$alpha_mean))
beta_premiss <- matrix(NA,sim_n,length(resultListANWC$beta_mean))
gamma_premiss <- matrix(NA,sim_n,length(resultListANWC$gamma_mean))

ANWC_b_mat <- matrix(NA,sim_n,length(resultListANWC$mar_varb))
MAR_b_mat <- matrix(NA,sim_n,length(resultListMAR$mar_varb))

ANWC_cond_pt <- MAR_cond_pt <- matrix(NA,sim_n,8)
ANWC_cond_var <- MAR_cond_var <- matrix(NA,sim_n,8)
ANWC_varb <- MAR_varb <-  matrix(NA,sim_n,8)

Tx_premiss_origin <- Ty_premiss_origin <- rep(NA,sim_n)

ans_X0_origin <- ans_X1_origin <- ans_Y0_origin <- ans_Y1_origin <- matrix(NA, sim_n, 2)
ul_X0_origin <- ul_X1_origin <- ul_Y0_origin <- ul_Y1_origin <- matrix(NA, sim_n, 2)

for (i in 1:n_sim){
  print(i)
  load(paste("./Chp4/pop1_N/Mis_",i,".RData",sep=""))
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean,resultListANWC$beta_mean, resultListANWC$gamma_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$beta_var, resultListANWC$gamma_var)
  ANWC_Tx[i] <- resultListANWC$total_mean_X
  ANWC_Tx_var[i] <- resultListANWC$total_varX
  ANWC_Ty[i] <- resultListANWC$total_mean_Y
  ANWC_Ty_var[i] <- resultListANWC$total_varY
  
  MAR_mean[i,] <- c(resultListMAR$alpha_mean,resultListMAR$beta_mean, resultListMAR$gamma_mean)
  MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$beta_var, resultListMAR$gamma_var)
  MAR_Tx[i] <- resultListMAR$total_mean_X
  MAR_Tx_var[i] <- resultListMAR$total_varX
  MAR_Ty[i] <- resultListMAR$total_mean_Y
  MAR_Ty_var[i] <- resultListMAR$total_varY
  
  ANWC_cond_pt[i,] <- resultListANWC$cond_pt
  MAR_cond_pt[i,] <- resultListMAR$cond_pt
  ANWC_cond_var[i,] <- resultListANWC$cond_var
  MAR_cond_var[i,] <- resultListMAR$cond_var
  ANWC_varb[i,] <- resultListANWC$var_b
  MAR_varb[i,] <- resultListMAR$var_b
  ANWC_b_mat[i,] <- resultListANWC$mar_varb
  MAR_b_mat[i,] <- resultListMAR$mar_varb
  
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y,sub_dat$U,sub_dat$Rx,sub_dat$origin_W),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y","U","Rx","Origin_W")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$pi <- 1/test_dat$W
  test_dat$Y <- as.factor(test_dat$Y)
  test_dat$U <- as.factor(test_dat$U)
  test_dat$Rx <- as.factor(test_dat$Rx)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m2 <- svyglm(X1~Y+U,design = mydesign,family = quasibinomial(link = "probit"))
  m4 <- svyglm(Y~U,design = mydesign,family = quasibinomial(link = "probit"))
  m6 <- glm(Rx~Y,data = test_dat[which(test_dat$U==0),], family = binomial(link = "probit"))
  alpha_premiss[i,] <- coef(m4)
  beta_premiss[i,] <- coef(m2)
  gamma_premiss[i,] <- coef(m6)
  
  test_dat$W <- as.numeric(test_dat$Origin_W)/100
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  
  wgt_tblY <- svyby(~as.factor(X1),~as.factor(Y),mydesign,svymean)
  wgt_dfY <- data.frame(wgt_tblY)
  ans_Y0_origin[i,] <- c(t(wgt_dfY[1,2:3]))
  ul_Y0_origin[i,] <- c(t(wgt_dfY[1,4:5]))
  ans_Y1_origin[i,] <- c(t(wgt_dfY[2,2:3]))
  ul_Y1_origin[i,] <- c(t(wgt_dfY[2,4:5]))
  
  wgt_tblX <- svyby(~as.factor(Y),~as.factor(X1),mydesign,svymean)
  wgt_dfX <- data.frame(wgt_tblX)
  ans_X0_origin[i,] <- c(t(wgt_dfX[1,2:3]))
  ul_X0_origin[i,] <- c(t(wgt_dfX[1,4:5]))
  ans_X1_origin[i,] <- c(t(wgt_dfX[2,2:3]))
  ul_X1_origin[i,] <- c(t(wgt_dfX[2,4:5]))
  
  Tx_premiss_origin[i] <- sum(sub_dat$X1*sub_dat$origin_W/100)
  Ty_premiss_origin[i] <- sum(sub_dat$Y*sub_dat$origin_W/100)
}

# truth
colMeans(alpha_premiss)
colMeans(beta_premiss)
colMeans(gamma_premiss)
sum(pop_dat$Y)
sum(pop_dat$X1)

t2_truth <- c(colMeans(alpha_premiss),colMeans(beta_premiss),colMeans(gamma_premiss))

# pre-miss sd
apply(alpha_premiss,2,sd)
apply(beta_premiss,2,sd)
apply(gamma_premiss,2,sd)
sd(Ty_premiss_origin)
sd(Tx_premiss_origin)

t2_presd <- c(apply(alpha_premiss,2,sd),apply(beta_premiss,2,sd),apply(gamma_premiss,2,sd))

# colmeans
colMeans(ANWC_mean)
mean(ANWC_Ty)
mean(ANWC_Tx)

t2_MDAM_est <- colMeans(ANWC_mean)

# ANWC result
# coverage of CI
unit_ANWC2 <- data.frame(ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_Tx, ANWC_TV = ANWC_Tx_var,ANWC_Y = ANWC_Ty,ANWC_Y_var = ANWC_Ty_var)
1-sum(unit_ANWC2$ANWC_mean.1 - 2.009*sqrt(unit_ANWC2$ANWC_var.1) > colMeans(alpha_premiss)[1]|unit_ANWC2$ANWC_mean.1 + 2.009*sqrt(unit_ANWC2$ANWC_var.1) < colMeans(alpha_premiss)[1])/sim_n
1-sum(unit_ANWC2$ANWC_mean.2 - 2.009*sqrt(unit_ANWC2$ANWC_var.2) > colMeans(alpha_premiss)[2]|unit_ANWC2$ANWC_mean.2 + 2.009*sqrt(unit_ANWC2$ANWC_var.2) < colMeans(alpha_premiss)[2])/sim_n

1-sum(unit_ANWC2$ANWC_mean.3 - 2.009*sqrt(unit_ANWC2$ANWC_var.3) > colMeans(beta_premiss)[1]|unit_ANWC2$ANWC_mean.3 + 2.009*sqrt(unit_ANWC2$ANWC_var.3) < colMeans(beta_premiss)[1])/sim_n
1-sum(unit_ANWC2$ANWC_mean.4 - 2.009*sqrt(unit_ANWC2$ANWC_var.4) > colMeans(beta_premiss)[2]|unit_ANWC2$ANWC_mean.4 + 2.009*sqrt(unit_ANWC2$ANWC_var.4) < colMeans(beta_premiss)[2])/sim_n
1-sum(unit_ANWC2$ANWC_mean.5 - 2.009*sqrt(unit_ANWC2$ANWC_var.5) > colMeans(beta_premiss)[3]|unit_ANWC2$ANWC_mean.5 + 2.009*sqrt(unit_ANWC2$ANWC_var.5) < colMeans(beta_premiss)[3])/sim_n

1-sum(unit_ANWC2$ANWC_mean.6 - 2.009*sqrt(unit_ANWC2$ANWC_var.6) > colMeans(gamma_premiss)[1]|unit_ANWC2$ANWC_mean.6 + 2.009*sqrt(unit_ANWC2$ANWC_var.6) < colMeans(gamma_premiss)[1])/sim_n
1-sum(unit_ANWC2$ANWC_mean.7 - 2.009*sqrt(unit_ANWC2$ANWC_var.7) > colMeans(gamma_premiss)[2]|unit_ANWC2$ANWC_mean.7 + 2.009*sqrt(unit_ANWC2$ANWC_var.7) < colMeans(gamma_premiss)[2])/sim_n

1-sum(unit_ANWC2$ANWC_Y - 2.009*sqrt(unit_ANWC2$ANWC_Y_var) >  sum(pop_dat$Y) | unit_ANWC2$ANWC_Y + 2.009*sqrt(unit_ANWC2$ANWC_Y_var) <  sum(pop_dat$Y))/sim_n
1-sum(unit_ANWC2$ANWC_T - 2.009*sqrt(unit_ANWC2$ANWC_TV) >  sum(pop_dat$X1) | unit_ANWC2$ANWC_T + 2.009*sqrt(unit_ANWC2$ANWC_TV) <  sum(pop_dat$X1))/sim_n

t2_MDAM_CI <- c(
  1-sum(unit_ANWC2$ANWC_mean.1 - 2.009*sqrt(unit_ANWC2$ANWC_var.1) > colMeans(alpha_premiss)[1]|unit_ANWC2$ANWC_mean.1 + 2.009*sqrt(unit_ANWC2$ANWC_var.1) < colMeans(alpha_premiss)[1])/sim_n,
  1-sum(unit_ANWC2$ANWC_mean.2 - 2.009*sqrt(unit_ANWC2$ANWC_var.2) > colMeans(alpha_premiss)[2]|unit_ANWC2$ANWC_mean.2 + 2.009*sqrt(unit_ANWC2$ANWC_var.2) < colMeans(alpha_premiss)[2])/sim_n,
  1-sum(unit_ANWC2$ANWC_mean.3 - 2.009*sqrt(unit_ANWC2$ANWC_var.3) > colMeans(beta_premiss)[1]|unit_ANWC2$ANWC_mean.3 + 2.009*sqrt(unit_ANWC2$ANWC_var.3) < colMeans(beta_premiss)[1])/sim_n,
  1-sum(unit_ANWC2$ANWC_mean.4 - 2.009*sqrt(unit_ANWC2$ANWC_var.4) > colMeans(beta_premiss)[2]|unit_ANWC2$ANWC_mean.4 + 2.009*sqrt(unit_ANWC2$ANWC_var.4) < colMeans(beta_premiss)[2])/sim_n,
  1-sum(unit_ANWC2$ANWC_mean.5 - 2.009*sqrt(unit_ANWC2$ANWC_var.5) > colMeans(beta_premiss)[3]|unit_ANWC2$ANWC_mean.5 + 2.009*sqrt(unit_ANWC2$ANWC_var.5) < colMeans(beta_premiss)[3])/sim_n,
  1-sum(unit_ANWC2$ANWC_mean.6 - 2.009*sqrt(unit_ANWC2$ANWC_var.6) > colMeans(gamma_premiss)[1]|unit_ANWC2$ANWC_mean.6 + 2.009*sqrt(unit_ANWC2$ANWC_var.6) < colMeans(gamma_premiss)[1])/sim_n,
  1-sum(unit_ANWC2$ANWC_mean.7 - 2.009*sqrt(unit_ANWC2$ANWC_var.7) > colMeans(gamma_premiss)[2]|unit_ANWC2$ANWC_mean.7 + 2.009*sqrt(unit_ANWC2$ANWC_var.7) < colMeans(gamma_premiss)[2])/sim_n
)

# sd of mean
apply(ANWC_mean,2,sd)
sd(unit_ANWC2$ANWC_Y)
sd(unit_ANWC2$ANWC_T)

t2_MDAM_sd <- apply(ANWC_mean,2,sd)

# mean of sd
sqrt(colMeans(ANWC_var))
sqrt(mean(unit_ANWC2$ANWC_Y_var))
sqrt(mean(unit_ANWC2$ANWC_TV))

t2_MDAM_AEsd <- sqrt(colMeans(ANWC_var))

# sqrt(b/50)
sqrt(colMeans(ANWC_b_mat)/50)

# MAR result
# colmeans
colMeans(MAR_mean)
mean(MAR_Ty)
mean(MAR_Tx)

t2_ICIN_est <- colMeans(MAR_mean)

# coverage of CI
unit_MAR2 <- data.frame(MAR_mean = MAR_mean, MAR_var = MAR_var, MAR_T = MAR_Tx, MAR_TV = MAR_Tx_var,MAR_Y = MAR_Ty,MAR_Y_var = MAR_Ty_var)
1-sum(unit_MAR2$MAR_mean.1 - 2.009*sqrt(unit_MAR2$MAR_var.1) > colMeans(alpha_premiss)[1]|unit_MAR2$MAR_mean.1 + 2.009*sqrt(unit_MAR2$MAR_var.1) < colMeans(alpha_premiss)[1])/sim_n
1-sum(unit_MAR2$MAR_mean.2 - 2.009*sqrt(unit_MAR2$MAR_var.2) > colMeans(alpha_premiss)[2]|unit_MAR2$MAR_mean.2 + 2.009*sqrt(unit_MAR2$MAR_var.2) < colMeans(alpha_premiss)[2])/sim_n

1-sum(unit_MAR2$MAR_mean.3 - 2.009*sqrt(unit_MAR2$MAR_var.3) > colMeans(beta_premiss)[1]|unit_MAR2$MAR_mean.3 + 2.009*sqrt(unit_MAR2$MAR_var.3) < colMeans(beta_premiss)[1])/sim_n
1-sum(unit_MAR2$MAR_mean.4 - 2.009*sqrt(unit_MAR2$MAR_var.4) > colMeans(beta_premiss)[2]|unit_MAR2$MAR_mean.4 + 2.009*sqrt(unit_MAR2$MAR_var.4) < colMeans(beta_premiss)[2])/sim_n
1-sum(unit_MAR2$MAR_mean.5 - 2.009*sqrt(unit_MAR2$MAR_var.5) > colMeans(beta_premiss)[3]|unit_MAR2$MAR_mean.5 + 2.009*sqrt(unit_MAR2$MAR_var.5) < colMeans(beta_premiss)[3])/sim_n

1-sum(unit_MAR2$MAR_mean.6 - 2.009*sqrt(unit_MAR2$MAR_var.6) > colMeans(gamma_premiss)[1]|unit_MAR2$MAR_mean.6 + 2.009*sqrt(unit_MAR2$MAR_var.6) < colMeans(gamma_premiss)[1])/sim_n
1-sum(unit_MAR2$MAR_mean.7 - 2.009*sqrt(unit_MAR2$MAR_var.7) > colMeans(gamma_premiss)[2]|unit_MAR2$MAR_mean.7 + 2.009*sqrt(unit_MAR2$MAR_var.7) < colMeans(gamma_premiss)[2])/sim_n

1-sum(unit_MAR2$MAR_Y - 2.009*sqrt(unit_MAR2$MAR_Y_var) >  sum(pop_dat$Y) | unit_MAR2$MAR_Y + 2.009*sqrt(unit_MAR2$MAR_Y_var) <  sum(pop_dat$Y))/sim_n
1-sum(unit_MAR2$MAR_T - 2.009*sqrt(unit_MAR2$MAR_TV) >  sum(pop_dat$X1) | unit_MAR2$MAR_T + 2.009*sqrt(unit_MAR2$MAR_TV) <  sum(pop_dat$X1))/sim_n

t2_ICIN_CI <- c(
  1-sum(unit_MAR2$MAR_mean.1 - 2.009*sqrt(unit_MAR2$MAR_var.1) > colMeans(alpha_premiss)[1]|unit_MAR2$MAR_mean.1 + 2.009*sqrt(unit_MAR2$MAR_var.1) < colMeans(alpha_premiss)[1])/sim_n,
  1-sum(unit_MAR2$MAR_mean.2 - 2.009*sqrt(unit_MAR2$MAR_var.2) > colMeans(alpha_premiss)[2]|unit_MAR2$MAR_mean.2 + 2.009*sqrt(unit_MAR2$MAR_var.2) < colMeans(alpha_premiss)[2])/sim_n,
  1-sum(unit_MAR2$MAR_mean.3 - 2.009*sqrt(unit_MAR2$MAR_var.3) > colMeans(beta_premiss)[1]|unit_MAR2$MAR_mean.3 + 2.009*sqrt(unit_MAR2$MAR_var.3) < colMeans(beta_premiss)[1])/sim_n,
  1-sum(unit_MAR2$MAR_mean.4 - 2.009*sqrt(unit_MAR2$MAR_var.4) > colMeans(beta_premiss)[2]|unit_MAR2$MAR_mean.4 + 2.009*sqrt(unit_MAR2$MAR_var.4) < colMeans(beta_premiss)[2])/sim_n,
  1-sum(unit_MAR2$MAR_mean.5 - 2.009*sqrt(unit_MAR2$MAR_var.5) > colMeans(beta_premiss)[3]|unit_MAR2$MAR_mean.5 + 2.009*sqrt(unit_MAR2$MAR_var.5) < colMeans(beta_premiss)[3])/sim_n,
  1-sum(unit_MAR2$MAR_mean.6 - 2.009*sqrt(unit_MAR2$MAR_var.6) > colMeans(gamma_premiss)[1]|unit_MAR2$MAR_mean.6 + 2.009*sqrt(unit_MAR2$MAR_var.6) < colMeans(gamma_premiss)[1])/sim_n,
  1-sum(unit_MAR2$MAR_mean.7 - 2.009*sqrt(unit_MAR2$MAR_var.7) > colMeans(gamma_premiss)[2]|unit_MAR2$MAR_mean.7 + 2.009*sqrt(unit_MAR2$MAR_var.7) < colMeans(gamma_premiss)[2])/sim_n
)

# sd of mean
apply(MAR_mean,2,sd)
sd(unit_MAR2$MAR_Y)
sd(unit_MAR2$MAR_T)

t2_ICIN_sd <- apply(MAR_mean,2,sd)

# mean of sd
sqrt(colMeans(MAR_var))
sqrt(mean(unit_MAR2$MAR_Y_var))
sqrt(mean(unit_MAR2$MAR_TV))

t2_ICIN_AEsd <- sqrt(colMeans(MAR_var))


# sqrt(b/50)
sqrt(colMeans(MAR_b_mat)/50)

# proportion truth
pop_dat_Y1 <- pop_dat[which(pop_dat$Y == 1),]
ans_Y1_true <- 1-mean(pop_dat_Y1$X1)

pop_dat_Y0 <- pop_dat[which(pop_dat$Y == 0),]
ans_Y0_true <- 1-mean(pop_dat_Y0$X1)

pop_dat_X0 <- pop_dat[which(pop_dat$X1 == 0),]
ans_X0_true <- 1-mean(pop_dat_X0$Y)

pop_dat_X1 <- pop_dat[which(pop_dat$X1 == 1),]
ans_X1_true <- 1-mean(pop_dat_X1$Y)

ans_X0_true
ans_X1_true
ans_Y0_true
ans_Y1_true

t1_truth <- c(sum(pop_dat$Y),
              sum(pop_dat$X1),
              ans_X0_true,
              ans_X1_true,
              ans_Y0_true,
              ans_Y1_true)

# premiss sd based on samples with original weights
apply(ans_X0_origin,2,sd)
apply(ans_X1_origin,2,sd)
apply(ans_Y0_origin,2,sd)
apply(ans_Y1_origin,2,sd)

t1_presd <- c(
  sd(Ty_premiss_origin),
  sd(Tx_premiss_origin),
  apply(ans_X0_origin,2,sd)[1],
  apply(ans_X1_origin,2,sd)[1],
  apply(ans_Y0_origin,2,sd)[1],
  apply(ans_Y1_origin,2,sd)[1]
)

# estimate
colMeans(ANWC_cond_pt)
colMeans(MAR_cond_pt)

t1_MDAM_est <- c(mean(ANWC_Ty),mean(ANWC_Tx),colMeans(ANWC_cond_pt)[c(1,3,5,7)])
t1_ICIN_est <- c(mean(MAR_Ty),mean(MAR_Tx),colMeans(MAR_cond_pt)[c(1,3,5,7)])

# CI coverage for ANWC of conditional proportions, based on truth
1-sum(ANWC_cond_pt[,1] - 2.009*sqrt(ANWC_cond_var[,1]) > 
        ans_X0_true|ANWC_cond_pt[,1] + 2.009*sqrt(ANWC_cond_var[,1]) 
      < ans_X0_true)/sim_n
1-sum(ANWC_cond_pt[,3] - 2.009*sqrt(ANWC_cond_var[,3]) > 
        ans_X1_true|ANWC_cond_pt[,3] + 2.009*sqrt(ANWC_cond_var[,3]) 
      < ans_X1_true)/sim_n
1-sum(ANWC_cond_pt[,5] - 2.009*sqrt(ANWC_cond_var[,5]) > 
        ans_Y0_true|ANWC_cond_pt[,5] + 2.009*sqrt(ANWC_cond_var[,5]) 
      < ans_Y0_true)/sim_n
1-sum(ANWC_cond_pt[,7] - 2.009*sqrt(ANWC_cond_var[,7]) > 
        ans_Y1_true|ANWC_cond_pt[,7] + 2.009*sqrt(ANWC_cond_var[,7]) 
      < ans_Y1_true)/sim_n

t1_MDAM_CI <- c(
  1-sum(unit_ANWC2$ANWC_Y - 2.009*sqrt(unit_ANWC2$ANWC_Y_var) >  sum(pop_dat$Y) | unit_ANWC2$ANWC_Y + 2.009*sqrt(unit_ANWC2$ANWC_Y_var) <  sum(pop_dat$Y))/sim_n,
  1-sum(unit_ANWC2$ANWC_T - 2.009*sqrt(unit_ANWC2$ANWC_TV) >  sum(pop_dat$X1) | unit_ANWC2$ANWC_T + 2.009*sqrt(unit_ANWC2$ANWC_TV) <  sum(pop_dat$X1))/sim_n,
  1-sum(ANWC_cond_pt[,1] - 2.009*sqrt(ANWC_cond_var[,1]) > 
          ans_X0_true|ANWC_cond_pt[,1] + 2.009*sqrt(ANWC_cond_var[,1]) 
        < ans_X0_true)/sim_n,
  1-sum(ANWC_cond_pt[,3] - 2.009*sqrt(ANWC_cond_var[,3]) > 
          ans_X1_true|ANWC_cond_pt[,3] + 2.009*sqrt(ANWC_cond_var[,3]) 
        < ans_X1_true)/sim_n,
  1-sum(ANWC_cond_pt[,5] - 2.009*sqrt(ANWC_cond_var[,5]) > 
          ans_Y0_true|ANWC_cond_pt[,5] + 2.009*sqrt(ANWC_cond_var[,5]) 
        < ans_Y0_true)/sim_n,
  1-sum(ANWC_cond_pt[,7] - 2.009*sqrt(ANWC_cond_var[,7]) > 
          ans_Y1_true|ANWC_cond_pt[,7] + 2.009*sqrt(ANWC_cond_var[,7]) 
        < ans_Y1_true)/sim_n
)

# CI coverage for MAR of conditional proportions, based on truth
1-sum(MAR_cond_pt[,1] - 2.009*sqrt(MAR_cond_var[,1]) > 
        ans_X0_true|MAR_cond_pt[,1] + 2.009*sqrt(MAR_cond_var[,1]) 
      < ans_X0_true)/sim_n
1-sum(MAR_cond_pt[,3] - 2.009*sqrt(MAR_cond_var[,3]) > 
        ans_X1_true|MAR_cond_pt[,3] + 2.009*sqrt(MAR_cond_var[,3]) 
      < ans_X1_true)/sim_n
1-sum(MAR_cond_pt[,5] - 2.009*sqrt(MAR_cond_var[,5]) > 
        ans_Y0_true|MAR_cond_pt[,5] + 2.009*sqrt(MAR_cond_var[,5]) 
      < ans_Y0_true)/sim_n
1-sum(MAR_cond_pt[,7] - 2.009*sqrt(MAR_cond_var[,7]) > 
        ans_Y1_true|MAR_cond_pt[,7] + 2.009*sqrt(MAR_cond_var[,7]) 
      < ans_Y1_true)/sim_n

t1_ICIN_CI <- c(
  1-sum(unit_MAR2$MAR_Y - 2.009*sqrt(unit_MAR2$MAR_Y_var) >  sum(pop_dat$Y) | unit_MAR2$MAR_Y + 2.009*sqrt(unit_MAR2$MAR_Y_var) <  sum(pop_dat$Y))/sim_n,
  1-sum(unit_MAR2$MAR_T - 2.009*sqrt(unit_MAR2$MAR_TV) >  sum(pop_dat$X1) | unit_MAR2$MAR_T + 2.009*sqrt(unit_MAR2$MAR_TV) <  sum(pop_dat$X1))/sim_n,
  1-sum(MAR_cond_pt[,1] - 2.009*sqrt(MAR_cond_var[,1]) > 
          ans_X0_true|MAR_cond_pt[,1] + 2.009*sqrt(MAR_cond_var[,1]) 
        < ans_X0_true)/sim_n,
  1-sum(MAR_cond_pt[,3] - 2.009*sqrt(MAR_cond_var[,3]) > 
          ans_X1_true|MAR_cond_pt[,3] + 2.009*sqrt(MAR_cond_var[,3]) 
        < ans_X1_true)/sim_n,
  1-sum(MAR_cond_pt[,5] - 2.009*sqrt(MAR_cond_var[,5]) > 
          ans_Y0_true|MAR_cond_pt[,5] + 2.009*sqrt(MAR_cond_var[,5]) 
        < ans_Y0_true)/sim_n,
  1-sum(MAR_cond_pt[,7] - 2.009*sqrt(MAR_cond_var[,7]) > 
          ans_Y1_true|MAR_cond_pt[,7] + 2.009*sqrt(MAR_cond_var[,7]) 
        < ans_Y1_true)/sim_n
)

# sd
apply(ANWC_cond_pt,2,sd)

t1_MDAM_sd <- c(sd(unit_ANWC2$ANWC_Y),sd(unit_ANWC2$ANWC_T),apply(ANWC_cond_pt,2,sd)[c(1,3,5,7)])

# avg. est. sd
sqrt(colMeans(ANWC_cond_var))

t1_MDAM_AEsd <- c(sqrt(mean(unit_ANWC2$ANWC_Y_var)),sqrt(mean(unit_ANWC2$ANWC_TV)),sqrt(colMeans(ANWC_cond_var))[c(1,3,5,7)])

# var_b
sqrt(colMeans(ANWC_varb)/50)
t1_MDAM_bm <- c(sqrt(colMeans(ANWC_b_mat)/50)[c(9,8)],sqrt(colMeans(ANWC_varb)/50)[c(1,3,5,7)])

# sd
apply(MAR_cond_pt,2,sd)
t1_ICIN_sd <- c(sd(unit_MAR2$MAR_Y),sd(unit_MAR2$MAR_T),apply(MAR_cond_pt,2,sd)[c(1,3,5,7)])

# avg. est. sd
sqrt(colMeans(MAR_cond_var))

t1_ICIN_AEsd <- c(sqrt(mean(unit_MAR2$MAR_Y_var)),sqrt(mean(unit_MAR2$MAR_TV)),sqrt(colMeans(MAR_cond_var))[c(1,3,5,7)])


# var_b
sqrt(colMeans(MAR_varb)/50)

########################################
# convert to dataframe
t1 <- data.frame(cbind(t1_truth,t1_MDAM_est,t1_ICIN_est,t1_MDAM_CI,t1_ICIN_CI,t1_presd,t1_MDAM_sd,t1_ICIN_sd,t1_MDAM_AEsd,t1_ICIN_AEsd,t1_MDAM_bm))
t2 <- data.frame(cbind(t2_truth,t2_MDAM_est,t2_ICIN_est,t2_MDAM_CI,t2_ICIN_CI,t2_presd,t2_MDAM_sd,t2_ICIN_sd,t2_MDAM_AEsd,t2_ICIN_AEsd))
t1$t1_MDAM_CI <- t1$t1_MDAM_CI*100
t1$t1_ICIN_CI <- t1$t1_ICIN_CI*100

t2$t2_MDAM_CI <- t2$t2_MDAM_CI*100
t2$t2_ICIN_CI <- t2$t2_ICIN_CI*100

}
xtable(t1,digits = c(3,3,3,3,1,1,3,3,3,3,3,3))
xtable(t2,digits = c(1,3,3,3,1,1,3,3,3,3,3))

round(t1$t1_MDAM_sd,digits = 3)
round(t1$t1_ICIN_sd,digits = 3)
