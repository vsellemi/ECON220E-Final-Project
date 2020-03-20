library("glmnet")
library("rapportools")
library("pracma")
library("matlib")

PriceRisk_OLS <- function(Ri,gt,ht) {
  # for the factors gt by controling ht
  # using the Feng, Giglio and Xiu (2017) approach.
  
  # Input:
  # Ri, nxT stock excess returns
  # gt, dxT testing factors
  # ht, pxT control factors
  
  # Output:
  # output for the No Selection OLS Estimation
  # lambdag_ols, se_ols, and lambda_ols
  
  # robust to missing data
  if (sum(is.nan(Ri)) > 0) {print("missing data in returns Ri - will append to zero")}
  if (sum(is.nan(gt)) > 0) {print("missing data in factors gt - will append to zero")}
  if (sum(is.nan(ht)) > 0) {print("missing data in factors ht - will append to zero")}
  
  # impute zero for missing observations
  Ri[is.nan(Ri)] <- 0
  gt[is.nan(gt)] <- 0
  ht[is.nan(ht)] <- 0
  
  # data information
  n <- dim(Ri)[1]  
  T <- dim(Ri)[2]
  d <- dim(gt)[1]
  p <- dim(ht)[1]
  
  cov_h <- matrix(NaN, nrow = n, ncol = p)
  for (nn in 1:n){
    temp <- cov( cbind((as.matrix(Ri[nn,])), t(ht)) )
    cov_h[nn,] <- temp[1,2:ncol(temp)]
  }
  cov_g <- matrix(NaN, nrow = n, ncol = d)
  for (nn in 1:n) {
    temp <- cov( cbind(as.matrix(Ri[nn,]), t(as.matrix(gt))) )
    cov_g[nn,] <- temp[1,2:ncol(temp)]
  }
  
  ER <- rowMeans(Ri)
  
  # for later use
  ones <- as.matrix(rep(1,n))
  M0   <- eye(n) - ones%*% inv(t(ones)%*%ones)%*%t(ones)
  
  # for no selection OLS
  # no selection estimate
  X <- cbind(cov_g,cov_h)
  X_zero <- cbind(rep(1,n),cov_g,cov_h)
  lambda_full_zero <- inv( t(X_zero) %*% X_zero )%*%(t(X_zero)%*%ER)
  lambda_full      <- inv(t(X)%*%M0%*%X)%*%(t(X)%*%M0%*%ER)
  lambdag_ols      <- lambda_full[1:d]
  rm(X)
  nomissing <- which(colSums(is.nan(rbind(ht,gt))) == 0)
  Lnm <- length(nomissing)
  
  # calculate avar_ols
  zthat3 <- matrix(NaN,nrow = d, ncol = Lnm)
  for (i in 1:d) {
    M_mdl <- eye(Lnm) - t(ht[,nomissing])%*%inv(ht[,nomissing]%*%t(ht[,nomissing]))%*%ht[,nomissing]
    zthat3[i,] <- M_mdl%*%as.matrix(gt[i,nomissing])
    rm(M_mdl)
  }
  Sigmazhat3 <- zthat3%*%t(zthat3)/Lnm
  
  vt <- rbind(gt[,nomissing],ht[,nomissing])
  temp3 <- matrix(0,nrow = d, ncol = d)
  i <- 0
  for (l in nomissing) {
    i <- i+1
    mt <- 1 - t(lambda_full)%*%rbind(as.matrix(gt[,l]),as.matrix(ht[,l]))
    temp3 <- temp3 + mt^2* (inv(Sigmazhat3)%*%zthat3[,i]%*%t(zthat3[,i])%*%inv(Sigmazhat3) )
  }
  avar_lambdag3 <- diag(temp3) / Lnm
  se3 <- sqrt(avar_lambdag3/Lnm)
  rm(temp3)
  
  # scaled lambda
  V_bar <- vt - as.matrix(rowMeans(vt))%*%t(as.matrix(rep(1,Lnm)))
  var_v <- V_bar%*%t(V_bar) / Lnm
  
  lambda_ols <- diag(var_v)%*%lambda_full
  #rm(list = c("X", "vt", "V_bar", "var_v", "lambda_full"))
  
  # output
  
  # output for OLS, No Selection
  result <- list("lambdag_ols" = lambdag_ols,
                 "se_ols" = se3,
                 "lambda_ols" = lambda_ols,
                 "lambda_ols_zero" = lambda_full_zero)
  
  ###################
  # new avar estimate
  ###################
  avar   <- matrix(0,nrow = p+d, ncol = p+d)
  meanBt <- matrix(0,nrow = p+d, 1)
  
  sigmavhat <- var_v
  vtbar <- V_bar
  Gammahat <- sigmavhat%*%lambda_full
  
  for (t in 1:Lnm){
    Bt <- inv(sigmavhat)%*%vtbar[,t]%*%t(vtbar[,t])%*%inv(sigmavhat)%*%Gammahat
    At <- inv(sigmavhat)%*%vtbar[,t]%*%t(vtbar[,t])%*%inv(sigmavhat)*drop(t(Gammahat)%*%inv(sigmavhat)%*%vtbar[,t])
    avar   <- avar - 2*At/Lnm + Bt%*%t(Bt) / Lnm
    meanBt <- meanBt + Bt/Lnm
  }
  avarhatLFM <- diag(inv(sigmavhat) + avar - (meanBt)%*%t(meanBt)/Lnm ) # diag(inv(sigmavhat))/T
  result1 <- list("se_new" = sqrt(avarhatLFM))
  return(c(result,result1))
}
