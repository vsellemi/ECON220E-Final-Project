
## simulating returns of test assets and factors
## for applications to double selection lasso approach

library("mvtnorm")
library("pracma")

# must calibrate parameters using Fama-French 5 factors:
# calibrate: chi, eta, lambda, Sigmaz, Ce, Ch1, and Sigmah
mean_Ce  <- as.matrix(rep(0,p))
cov_Ce   <- diag(rep(1,p))
mean_Ch1 <- as.matrix(rep(0,4))


# set dimensions
n <- 100
p <- 25
T <- 240

set.seed(69) # for reproducability
# (1) simulate Ce (nxd) and Ch1 (nx4) independently from multivariate normals
Ce  <- dmvnorm()
Ch1 <- dmvnorm()

# (2) calculate Ch2, initialize theta0 (p-4)x1, theta1 (p-4)x4, and Ceps nx(p-4) ~ N(m,S) 
theta0 <- matrix(0, nrow = p-4, ncol = 1) # initialized but should fill
theta1 <- matrix(0, nrow = p-4, ncol = 4)
Ceps   <- dmvnorm()
Ch2    <- matrix(1, nrow = n, ncol = 1) %*% t(theta0) + Ch1 %*% t(theta_1) + Ceps

# (3) Cg
xi  <- matrix(1, nrow = 1, ncol = 1)
chi <- matrix(0, nrow = d, ncol = p)
Ch  <- cbind(Ch1, Ch2)  # (nxp)

Cg  <- matrix(1, nrow = n, ncol = 1) %*% xi + C_h %*% t(chi) + Ce  # Cg ~ (nxd)

# (4) Cz
eta <- matrix(0, nrow = d, ncol = p)
Cz  <- Cg - Ch %*% t(eta) # (nxd) - (nxp)(pxd) 

# (5) Er
gamma0  <- matrix(1, nrow = 1, ncol =1)
lambdag <- matrix(c(1,0,0), nrow = 3, ncol = 1)  # one useful g, one useless, one redundant
lambdah <- rbind(matrix(c(1,1,1,1),nrow = 4, ncol = 1), matrix(0, nrow = p-4, ncol = 1)) # 4 useful h, p-4 useless
Ert  <- matrix(1, nrow = n, ncol = 1) %*% gamma0 + Cg %*% lambdag + Ch %*% lambdah

# (6) calculate betas
Sigmaz <- matrix(1,nrow = d, ncol = d)
Betag  <- C_z %*% inv(Sigmaz) # (nxd)(dxd)  #use pracma for inverse!!
Sigmah <- matrix(1,nrow = p, ncol = p) 
Betah  <- C_h %*% inv(Sigmah) # (nxp)(pxp)
Sigmau <- matrix(1,nrow = n, ncol = n)   #variance of sigmau disturbances
ut     <- dt(sample(1:1000, n, replace = TRUE)/1000,df = 5) #??    # draw (nx1) ut from student t distribution with 5 deg of freedom and Sigmau var
rt     <- Ert + Betag %*% as.matrix(gt) + Betah %*% as.matrix(ht) + ut

# (7) generate ht, zt -- >




