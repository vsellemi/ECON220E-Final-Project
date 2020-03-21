#%%

#----------------------------------------------------------------------------- 

#DOUBLE SELECTION LASSO                                

#----------------------------------------------------------------------------- 

#Author: Victor Sellemi

#*Based on 2020 Matlab Code of Guanhao Feng, Stefano Giglio and Dacheng Xiu*

#Date  : March 2020

#----------------------------------------------------------------------------- 

#%%



library("glmnet")

library("rapportools")

library("pracma")

#library("matlib")



# main results



main_dir <- "C:\\Users\\anind\\Google Drive\\Coursework\\Econ_220E_Metrics\\Proj\\Our code\\Final"

data_dir <- "C:\\Users\\anind\\Google Drive\\Coursework\\Econ_220E_Metrics\\Proj\\Our code\\Final\\data"



setwd(data_dir)



# fix the random seed for any additional CV

seed_num <- 100



# factor data

allfactors <- read.csv("factors.csv")                  

date       <- allfactors$Date

rf         <- allfactors$RF

factors    <- allfactors[,3:ncol(allfactors)]



L <- length(date)

P <- ncol(factors)



# test porfolios

port_5x5 <- read.csv("port_5x5.csv", header = FALSE) 

port_5x5 <- port_5x5[,2:ncol(port_5x5)]                                   # remove date col

port_5x5 <- port_5x5 - rf                                                 # excess returns



port_3x2 <- read.csv("port_3x2.csv", header = FALSE)

port_3x2 <- port_3x2[,2:ncol(port_3x2)]

port_3x2 <- port_3x2 - rf



port_202 <- read.csv("port202.csv", header = FALSE)

port_202 <- port_202[,2:ncol(port_202)]/100

port_202 <- port_202 - rf



# other information

summary         <- read.csv('summary.csv')

factorname      <- summary$Row

factorname_full <- summary$Description

year_pub        <- summary$Year

year_end        <- summary$Year_end



port_5x5_id <- read.csv("port_5x5_id.csv") 

port_3x2_id <- read.csv("port_3x2_id.csv")



mkt_ind <- which(is.element(factorname,"MktRf"))

smb_ind <- which(is.element(factorname,"SMB"))

hml_ind <- which(is.element(factorname,"HML"))



#%% form a smaller set of portfolios for bivariate sorted porfolios

kk = 10   # minimum number of stocks in a portfolio



include_3x2 <- which(port_3x2_id$min_stk6 >= kk)

port_3x2b   <- data.frame(matrix(nrow = nrow(port_3x2)))



for (i in 1:P) {
  
  if (is.element(i,include_3x2)) {
    
    port_3x2b <- cbind(port_3x2b, port_3x2[,(i*6-5):(i*6)])
    
  }
  
}



Ri <- port_3x2b # test asset
Ri = Ri[,2:ncol(Ri)]


# load tune_center for 200 randome seeds selected by cross-validations

# For each random seed in 1:200, we run a cross-validation and find the

# best tuning parameter.

# Then we take the average for the 200 selected tuning parameters at the

# log scale as tune center, because we plot heat maps on log(lambda).



# load cross validation results from xiu, giglio, feng (2020)

tune_center <- matrix(c(0.000000018,  1.8e-8, 1.8e-8, 1.80E-08, 1.80E-08, 1.80E-08, 1.80E-08,  1.80E-08, 1.80E-08,  1.80E-08,  1.80E-08, 1.80E-08, 1.80E-08, 1.80E-08,  1.80E-08, 0.00000000000176,4.30E-13, 8.11E-12, 3.25E-11,6.28E-12, 2.94E-11 
                        , 2.07E-12, 2.44E-11, 8.23E-12,1.01E-12,1.42E-12
                        , 2.45E-11, 2.52E-12, 2.44E-10, 7.04E-11  ), nrow=15, ncol =2)

  


# choose control factors before 2012

ContrlList <- which(year_pub < 2012)

ControlFactor <- factors[,ContrlList]

FF3 <- t(factors[,c(mkt_ind,smb_ind,hml_ind)])



# test factors since 2012

TestList <- which(year_pub >= 2012)

TestFactor <- factors[,TestList]



result <- data.frame(matrix(0,nrow = length(TestList),ncol = 10))

names(result) <- c("tstat_ds", "lambda_ds", "tstat_ss", "lambda_ss", "avg", "tstat_avg",
                   
                   "lambda_ols", "tstat_ols", "lambda_FF3", "tstat_FF3")



# test factor individually

#for (j in 1:length(TestList)) {

for (j in 1:2) {
  
  disp(j)
  
  gt <- t(TestFactor[,j]) # test factor
  
  ht <- t(ControlFactor)  # control factor
  
  
  
  # robust to missing data
  
  if (sum(is.na(Ri)) > 0) {print("missing data in returns Ri - will append to zero")}
  
  if (sum(is.nan(gt)) > 0) {print("missing data in factors gt - will append to zero")}
  
  if (sum(is.nan(ht)) > 0) {print("missing data in factors ht - will append to zero")}
  
  
  
  # impute zero for missing observations
  
  Ri[is.na(Ri)] <- 0
  
  gt[is.nan(gt)] <- 0
  
 # ht[is.nan(ht)] <- 0
  
  
  
  # use the average tuning parameter from 200 random seeds
  
  model_ds  <- DS(t(Ri), gt, ht, -log(tune_center[j,1]), -log(tune_center[j,2]),1,seed_num)
  
  tstat_ds  <- model_ds$lambdag_ds/model_ds$se_ds
  
  lambda_ds <- model_ds$gamma_ds[1]
  
  
  
  # Single-Selection results, replace with a huge tune2
  
  model_ss  <- DS(t(Ri), gt, ht, -log(tune_center[j,1]), -log(1),1,seed_num)
  
  tstat_ss  <- model_ss$lambdag_ds/model_ss$se_ds
  
  lambda_ss <- model_ss$gamma_ds[1]
  
  
  # controlling everything, no selection, OLS
  
  model_ols  <- PriceRisk_OLS(t(Ri), gt, ht)
  
  tstat_ols  <- model_ols$lambdag_ols/model_ols$se_ols
  
  lambda_ols <- model_ols$lambda_ols[1]
  
  
  
  # only control FF3 by OLS
  
  model_FF3 <- PriceRisk_OLS(t(Ri), gt, FF3)
  
  tstat_FF3 <- model_FF3$lambdag_ols/model_FF3$se_ols
  
  lambda_FF3 <- model_FF3$lambda_ols[1]
  
  
  
  # time series average
  
  avg <- mean(gt)
  
  tstat_avg <- avg/std(gt)*sqrt(sum(!is.nan(gt)))
  
  
  
  # combine the results in a table (data frame)
  
  #temp <- data.frame(tstat_ds, lambda_ds, tstat_ss, lambda_ss, avg, tstat_avg,
  
  #                   lambda_ols, tstat_ols, lambda_FF3, tstat_FF3)
  
  #result[,3:ncol(result)] <- rbind(result[,3:ncol(result)],temp)
  
  result$tstat_ds[j]   <- tstat_ds
  
  result$lambda_ds[j]  <- lambda_ds
  
  result$tstat_ss[j]   <- tstat_ss
  
  result$lambda_ss[j]  <- lambda_ss
  
  result$avg[j]        <- avg
  
  result$tstat_avg[j]  <- tstat_avg
  
  result$lambda_ols[j] <- lambda_ols
  
  result$tstat_ols[j]  <- tstat_ols
  
  result$lambda_FF3[j] <- lambda_FF3  
  
  result$tstat_FF3[j]  <- tstat_FF3  
  
  
  
  
  
}



show(factorname_full[model_ds$sel1])



# extract outputs

lambda_ds  <- result$lambda_ds * 10000 # bp

tstat_ds   <- result$tstat_ds



lambda_ss  <- result$lambda_ss*10000   # bp

tstat_ss   <- result$tstat_ss



lambda_FF3 <- result$lambda_FF3*10000  # bp

tstat_FF3  <- result$tstat_FF3



avg        <- result$avg*10000         # bp

tstat_avg  <- result$tstat_avg



lambda_ols <- result$lambda_ols*10000  # bp

tstat_ols  <- result$tstat_ols



# factor names for those factors since 2012

factornames <- factorname_full[TestList]



result <- data.frame(TestList,factornames,lambda_ds,tstat_ds,lambda_ss,tstat_ss,lambda_FF3,
                     
                     tstat_FF3,lambda_ols,tstat_ols,avg,tstat_avg)



# display the table

show(result)



# output Table as a CSV file

out_file <- "C:\\Users\\anind\\Google Drive\\Coursework\\Econ_220E_Metrics\\Proj\\Our code\\Final\\main.csv"


write.csv(result, file = out_file)
