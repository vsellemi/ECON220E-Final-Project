rm(list=ls())

library("glmnet")
library("rapportools")
library("pracma")
library("matlib")


DS <- function(Ri, gt, ht, tune1, tune2, alpha, seednum){
  ### The main function for Double-Selection implementation
  ### no cross-validation for 1st and 2nd selections
  
  # INPUT DIMENSIONS: Ri (nxT), ht (pxT), gt (dxT)
  
  # dim of g is 1!
  if (is.empty(alpha)) {alpha <- 1}
  
  # robust to missing data
  if (sum(is.nan(Ri)) > 0) {print("missing data in returns Ri - will append to zero")}
  if (sum(is.nan(gt)) > 0) {print("missing data in factors gt - will append to zero")}
  if (sum(is.nan(ht)) > 0) {print("missing data in factors ht - will append to zero")}
  
  # impute zero for missing observations
  Ri[is.nan(Ri)] <- 0
  gt[is.nan(gt)] <- 0
  ht[is.nan(ht)] <- 0
  
  # data information
  if (is.null(dim(rf))) {n <- length(Ri)} else {n <- dim(Ri)[1]}
  if (is.null(dim(ht))) {p <- length(ht)} else {p <- dim(ht)[1]}
  if (is.null(dim(gt))) {d <- length(gt)} else {d <- dim(gt)[1]}
  
  tmp1  <- cov(cbind(t(gt), t(Ri)))
  cov_g <- tmp1[(d+1):nrow(tmp1), 1:d]
  tmp2  <- cov(cbind(t(ht),t(Ri)))
  cov_h <- tmp2[(p+1):nrow(tmp2), 1:p]
  
  ER    <- rowMeans(Ri)
  
  beta  <- matrix(nrow = n, ncol = p)
  for (i in 1:p) {
    beta[,i] <- cov_h[,i] / var(ht[i,])
  }
  penalty <- colMeans(beta^2)
  penalty <- penalty/mean(penalty)   # check dim
  
  lambda0 <- exp(-seq(0,35,35/100))
  
  # 1st selection in cross-sectional regression
  model1 <- glmnet(x = cov_h*diag(penalty),
                   y = ER,
                   family = 'gaussian',
                   standardize = FALSE,
                   lambda = exp(-tune1),
                   alpha = alpha)
  model1_est <- coef.glmnet(model1)
  sel1 <- t( which(model1_est[2:(p+1)] != 0) )
  
  err1 <- mean( (ER - cbind(t(rep(1,n),cov_h*diag(penalty)))*model1_est)*2 )
  
  # 2nd selection
  sel2 <- vector()
  err2 <- matrix(nrow = d, ncol = 1)
  for (i in 1:d) {
    model2 <- glmnet(x = cov_h*diag(penalty),
                     y = cov_g[,i],
                     family = 'gaussian',
                     standardize = FALSE,
                     lambda = exp(-tune2),
                     alpha = alpha)
    model2_est <- coef.glmnet(model2)
    sel2 <- rbind(sel2, which(model2_est[2:length(model2_est)]!=0))
    err2[i] <- mean( (cov_g[,i]- cbind(t(rep(1,n)), cov_h*diag(penalty))*model2_est)^2)
  }
  sel2 <- t(unique(sel2))
  
  # 3rd selection for avar zt
  sel3 <- vector()
  for (i in 1:d) {
    TCSVout <- TCSV(Ri, gt[i,], ht, lambda0, 10, 1, alpha, seednum)
    sel3 <- rbind(sel3, TCSVout$sel3_1se)
  }
  sel3 <- unique(sel3)
  
  # post-selection estimation and inference
  dsout <- infer(Ri, gt, ht, sel1, sel2, sel3)
  ssout <- infer(Ri, gt, ht, sel1, sel2 = NULL, sel3) # why sel2 null?
  
  result <- list("lambdag_ds" = dsout$lambdag, # output for Double Selection
                 "se_ds" = dsout$se,
                 "gamma_ds" = dsout$gamma,
                 "lambdag_ss" = ssout$lambdag, # output for Single Selection
                 "se_ss" = ssout$se,
                 "gamma_ss" = ssout$gamma,
                 "sel1" = sel1,                # selection results
                 "sel2" = sel2,
                 "sel3" = sel3,
                 "select" = union(sel1,sel2),
                 "err1" = err1,
                 "err2" = err2)
}

# --------------------------------------------------------------- #

TCSV <- function(Ri, gt, ht, lambda, Kfld, Jrep, alpha, seednum) {
  ### the function for cross-validation over time
  ### only for 3rd selection in the DS function
  if (is.empty(seednum)) {seednum <- 101}
  
  # data information
  p <- dim(ht)[1]
  T <- dim(ht)[2]
  
  L <- length(lambda)
  
  cvm3  <- array(NaN, dim = c(L,Kfld,Jrep)) 
  cvm33 <- vector()
  
  nomissing <- t(colSums(is.nan(rbind(ht,gt))) == 0)
  
  for (j in 1:Jrep) {
    set.seed(seednum + j)
    indices <- sample(1:Kfld,T)
    for (k in 1:Kfld) {
      test  <- (indices == k)
      train <- (indices != k)
      
      ht_train <- ht[,intersect(train,nomissing)]
      gt_train <- gt[,intersect(train,nomissing)]
      
      ht_test <- ht[,test]
      gt_test <- gt[,test]
      model3  <- glmnet(x = t(ht_train),
                        y = t(gt_train),
                        family = 'gaussian',
                        intercept = FALSE, #why? 
                        standardize = TRUE, 
                        lambda = lambda,
                        alpha = alpha)
      gt_pred <- t(ht_test)*model3$beta
      
      LL3  <- length(model3$lambda)
      temp <- (repmat(t(gt_test),n=1,m=LL3) - gt_pred)^2
      temp[is.na(temp)] <- 0
      cvm[1:LL3,k,j] <- colMeans(temp)
    }
    cvm33 <- cbind[cvm33,cvm3[,,j]]
  }
  cv_sd3  <- std(t(cvm33)) / sqrt(Kfld*Jrep)
  cvm333  <- rowMeans(cvm33)
  l_sel3  <- which.min(cvm333)
  cvm33ub <- cvm333[l_sel3] + cv_sd3[l_sel3]
  l3_1se  <- which(cvm333[1:l_sel3] >= cvm33ub)
  l3_1se  <- l3_1se[length(l3_1se)]
  if (is.empty((l3_1se))) {l3_1se <- l_sel3}
  
  # to reestimate the model with all data, refit the model
  model3 <- glmnet(x = t(ht[,nomissing]),
                   y = t(gt[nomissing]),
                   family = 'gaussian',
                   intercept = FALSE,
                   standardize = TRUE,
                   lambda = lambda[c(l3_1se,l_sel3)],
                   alpha = alpha)
  sel3 <- which(model3$beta[,2] != 0)
  sel3_1se <- which(model3.beta[,1] != 0)
  output <- list("sel3" = sel3, "lambda3" = lambda[l_sel3],
                 "sel3_1se" = sel3_1se, "lambda3_1se" = lambda[l3_1se])
  return(output)
}


# -------------------------------------------------------------- #

infer <- function(Ri, gt, ht, sel1, sel2, sel3){
   ### the function for estimation and inference
  if (is.null(dim(Ri))) {n <- length(Ri)} else {n <- dim(Ri)[1]}
  if (is.null(dim(ht))) {p <- length(ht)} else {p <- dim(ht)[1]}
  if (is.null(dim(gt))) {d <- length(gt)} else {d <- dim(gt)[1]}
  
  tmp1  <- cov(cbind(t(gt),t(Ri)))
  cov_g <- tmp1[(d+1):nrow(tmp1), 1:d]
  tmp2  <- cov(cbind(t(ht),t(Ri)))
  cov_h <- tmp2[(p+1):nrow(tmp2),1:p]
  
  ER    <- rowMeans(Ri)
  M0    <- eye(n) - t(rep(1,n))* inv(rep(1,n)*t(rep(1,n)))*rep(1,n)
  
  nomissing <- which(colSums(is.nan(rbind(ht,gt))) == 0)
  Lnm       <- length(nomissing)
  select    <- union(sel1,sel2)
  
  X <- cbind(cov_g, cov_h[,select]) 
  lambda_full <- inv(t(X)*M0*X)*(t(X)*M0*ER)
  lambdag     <- lambda_full[1:d]
  rm(X)
  
  # for double selection inference: AVAR
  zthat = matrix(NaN, nrow = d, ncol = Lnm)
  for (i in 1:d){
    M_mdl <- eye(Lnm) - t(ht[sel3,nomissing])*inv(ht[sel3,nomissing]*t(ht[sel3,nomissing]))*ht[sel3,nomissing]
    zthat[i,] = M_mdl*t(gt[i,nomissing])
    rm(M_mdl)
  }
  Sigmazhat <- zthat*t(zthat) / Lnm
  temp2     <- matrix(0,nrow = d, ncol = d)
  ii <- 0 
  for (l in nomissing) {
    ii <- ii + 1
    mt    <- 1 - t(lambda_full)*rbind(gt[1:d,l],ht[select,l])
    temp2 <- (mt^2) * (inv(Sigmazhat)*zthat[,ii]*t(zthat[,ii])*inv(Sigmazhat))
  }
  avar_lambdag <- diag(temp2)/Lnm
  se <- sqrt(avar_lambdag/Lnm)
  rm(temp2)
  
  # scaled lambda for DS
  vt <- rbind(gt[,nomissing], ht[select,nomissing])
  V_bar <- vt - rowMeans(vt,2)*rep(1,Lnm)
  var_v <- V_bar*t(V_bar) / Lnm
  gamma <- diag(var_v)*lambda_full
  rm(list = c("X", "vt", "V_bar", "var_v", "lambda_full"))
  
  output <- list("lambdag" = lambdag,
                 "se" = se,
                 "gamma" = gamma)
}

