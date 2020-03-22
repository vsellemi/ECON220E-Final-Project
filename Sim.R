## simulating returns of test assets and factors

## for applications to double selection lasso approach



library("mvtnorm")

library("pracma")

library("MASS")


# set dimensions

n <- 100 # Number of assets

p <- 25

T <- 240

d <- 3

# Load rwal world asset data (to be used in calibration)
Ri <- port_3x2b # test asset
Ri = Ri[,2:ncol(Ri)]


set.seed(69) # for reproducability

# Load tune center

tune_center <- matrix(c(0.000000018,  1.8e-8, 1.8e-8, 1.80E-08, 1.80E-08, 1.80E-08, 1.80E-08,  1.80E-08, 1.80E-08,  1.80E-08,  1.80E-08, 1.80E-08, 1.80E-08, 1.80E-08,  1.80E-08, 0.00000000000176,4.30E-13, 8.11E-12, 3.25E-11,6.28E-12, 2.94E-11 
                        , 2.07E-12, 2.44E-11, 8.23E-12,1.01E-12,1.42E-12
                        , 2.45E-11, 2.52E-12, 2.44E-10, 7.04E-11  ), nrow=15, ncol =2)



# (1) simulate Ce (nxd) and Ch1 (nx4) independently from multivariate normals

# must calibrate parameters using Fama-French 5 factors:

# calibrate: chi, eta, lambda, Sigmaz, Ce, Ch1, and Sigmah

temp3 = cov(cbind(allfactors[c("MktRf","HML", "SMB", "UMD", "cash", "HML_Devil", "gma")], Ri)) # Using the cov matrix for the FF factors
Ch1 = temp3[8: 107, 1:4]
Ce = temp3[8:107, 5:7]

#mean_Ce  <- as.matrix(rep(0,d))

#cov_Ce   <- diag(rep(1,d))

#mean_Ch1 <- as.matrix(rep(0,4))

#cov_Ch1   <- diag(rep(1,4))

#Ce  <- mvrnorm(n, mean_Ce, cov_Ce)

#Ch1 <- mvrnorm(n, mean_Ch1, cov_Ch1)



# (2) calculate Ch2, initialize theta0 (p-4)x1, theta1 (p-4)x4, and Ceps nx(p-4) ~ N(m,S) 

theta0 <- matrix(0, nrow = p-4, ncol = 1) # Using a matrix of ones
theta1 <- matrix(1, nrow = p-4, ncol = 4)

mean_Ceps  <- as.matrix(rep(0,p-4))

cov_Ceps   <- diag(rep(1,p-4))

Ceps   <-  mvrnorm(n, mean_Ceps, cov_Ceps)

Ch2    <- matrix(1, nrow = n, ncol = 1) %*% t(theta0) + Ch1 %*% t(theta1) + Ceps



# (3) Cg

#xi  <- matrix(1, nrow = 1, ncol = d)

#chi <- matrix(0, nrow = d, ncol = p)

Ch  <- cbind(Ch1, Ch2)  # (nxp)

xi = t(matrix(cbind(as.matrix(1), t(as.matrix(rep(0,d-1)))), nrow = d, ncol = 1)) # no loadings on redundant factors
#chi = t(rbind(matrix(1, nrow = 4, ncol = d), matrix(0, nrow = (p-4), ncol = d)))
chi = t(cbind(rbind(as.matrix(rep(1,4)), as.matrix(rep(0,p-4))), matrix(0, nrow = p, ncol = 1),matrix(1, nrow = p, ncol = 1))) # no loading of the useless factor on Ch

Cg  <- matrix(1, nrow = n, ncol = 1) %*% xi + Ch %*% t(chi) + Ce  # Cg ~ (nxd)



# (4) Cz

#eta <- matrix(0, nrow = d, ncol = p)
#eta = rbind(matrix(1, nrow= 1, ncol = p), matrix(0, nrow = 2, ncol = p)) 

eta = cbind(matrix(1, nrow = 3, ncol = 4), matrix(0, nrow = 3, ncol = p-4)) # No loading of g on h2
Cz  <- Cg - Ch %*% t(eta) # (nxd) - (nxp)(pxd) 


# (5) Er

gamma0  <- matrix(1, nrow = 1, ncol =1)

lambdag <- matrix(c(1,0,0), nrow = 3, ncol = 1)  # one useful g, one useless, one redundant

lambdah <- rbind(matrix(c(1,1,1,1),nrow = 4, ncol = 1), matrix(0, nrow = p-4, ncol = 1)) # 4 useful h, p-4 useless

Ert  <- matrix(1, nrow = n, ncol = 1) %*% gamma0 + Cg %*% lambdag + Ch %*% lambdah



# (6) calculate betas

#Sigmaz <- matrix(1,nrow = d, ncol = d)

Sigmaz = diag(rep(1,d))
Sigmah = diag(rep(1,p))

Betag  <- Cz %*% inv(Sigmaz) # (nxd)(dxd)  #use pracma for inverse!!

#Sigmah <- matrix(1,nrow = p, ncol = p) 

Betah  <- Ch %*% inv(Sigmah) # (nxp)(pxp)




# Monte Carlo Simulations (Repeat 2000 times)

nsim <- 100

estlambda <- data.frame(matrix(0, ncol = 3, nrow = nsim))
tstatlambda <- data.frame(matrix(0, ncol = 3, nrow = nsim))

for (m in 1: nsim){
  
  disp(m)

  Sigmau <- matrix(1,nrow = n, ncol = n)   #variance of sigmau disturbances
  
  
  Rt = matrix(nrow = 100)
  Ht = matrix(ncol = p)
  Gt = matrix(ncol = d)
  
  # Draw T data points to create dataset 
  
  
  for (i in 1:T){
  
  
    ut     <- t(rmvt(1, sigma = diag(100), df = 5)) #??    # draw (nx1) ut from student t distribution with 5 deg of freedom and Sigmau var
    
    
    # (7) generate ht, zt -- >
    
    mean_ht  <- as.matrix(rep(0,p))
    
    ht  <- as.matrix(mvrnorm(1, mean_ht, Sigmah)) 
    #ht = t(temp2)                        
    
    mean_zt <- as.matrix(rep(0,d))
    
    zt <- as.matrix(mvrnorm(1, mean_zt, Sigmaz))
    
    gt <- as.matrix(eta)%*%as.matrix(ht) + zt
    
    rt     <- Ert + Betag %*% gt + Betah %*% ht  + as.matrix(ut)
    
    Rt = cbind(Rt, rt)
    Ht = rbind(Ht, t(ht))
    Gt = rbind(Gt, t(gt))
  
  }
            
  Rt = Rt[, 2:ncol(Rt)]
  Ht = Ht[2:nrow(Ht),]
  Gt = Gt[2:nrow(Gt),]
  
  
  # Enter data in DS model:
  
  # test factor individually
  
  for (j in 1:3) {
    
    gt <- t(Gt[,j]) # test factor
    
    ht <- t(Ht)  # control factor
    
    
    # use the average tuning parameter from 200 random seeds
    
    model_ds  <- DS(Rt, gt, ht, -log(tune_center[j,1]), -log(tune_center[j,2]),1,seed_num)
    
    tstat_ds  <- model_ds$lambdag_ds/model_ds$se_ds
    
    lambda_ds <- model_ds$gamma_ds[1]
  
    # combine the results in a table (data frame)
    
    #temp <- data.frame(tstat_ds, lambda_ds)
    
    #result[,3:ncol(result)] <- rbind(result[,3:ncol(result)],temp)
    
    #result$tstat_ds[j]   <- tstat_ds
    estlambda[m,j]  =  lambda_ds
    tstatlambda[m,j]  =  tstat_ds
  
  }
}

lambda_ds1 = colMeans(estlambda)
tstat_ds1 = colMeans(estlambda)

hist(estlambda[,1])

#simresult$lambda_ds = simresult$lambda_ds/nsim 
