rm(list=ls())

main_dir = "/Users/victorsellemi/Downloads/DoubleLasso"
data_dir = "/Users/victorsellemi/Downloads/DoubleLasso/data"

# --------------------------------------------------------------------------- #

#%% data import and cleaning

#getwd() 
setwd(data_dir)

# factor data
allfactors <- read.csv("factors.csv")                  
date       <- allfactors$Date
rf         <- allfactors$RF
factors    <- allfactors[,3:ncol(allfactors)]

L <- length(date)
P <- ncol(factors)

# test porfolios
port_5x5 <- read.csv("port_5x5.csv", header = FALSE) 
port_3x2 <- port_3x2[,2:ncol(port_5x5)]                                   # remove date col
port_5x5 <- port_5x5 - rf                                                 # excess returns

port_3x2 <- read.csv("port_3x2.csv", header = FALSE)
port_3x2 <- port_3x2[,2:ncol(port_3x2)]
port_3x2 <- port_3x2 - rf

#port_5x5_seq <- read.csv("port_5x5_seq.csv", header = FALSE)
#port_5x5_seq <- port_5x5_seq[,2:ncol(port_5x5_seq)]
#port_5x5_seq <- port_5x5_seq - rf

#port_3x2_seq <- read.csv("port_3x2_seq.csv", header = FALSE)
#port_3x2_seq <- port_3x2_seq[,2:ncol(port_3x2_seq)]
#port_3x2_seq <- port_3x2_seq - rf

port_202 <- read.csv("port202.csv", header = FALSE)
port_202 <- port_202[0,2:ncol(port_202)]/100
port_202 <- port_202 - rf

# other information
summary         <- read.csv('summary.csv')
factorname      <- summary$Row
factorname_full <- summary$Description
year_pub        <- summary$Year
year_end        <- summary$Year_end

port_5x5_id <- read.csv("port_5x5_id.csv") 
port_3x2_id <- read.csv("port_3x2_id.csv")

#%% form a smaller set of portfolios for bivariate sorted porfolios
kk = 10   # minimum number of stocks in a portfolio

include_3x2 <- which(port_3x2_id$min_stk6 >= kk)
port_3x2b   <- data.frame(matrix(nrow = nrow(port_3x2)))

for (i in 1:P) {
  if (is.element(i,include_3x2)) {
    port_3x2b <- cbind(port_3x2b, port_3x2[,(i*6-5):(i*6)])
  }
}

include_5x5 <- which(port_5x5_id$min_stk >= kk)
port_5x5b   <- data.frame(matrix(nrow = nrow(port_5x5)))

for (i in 1:P) {
  if (is.element(i,include_5x5))
    port_5x5b <- cbind(port_5x5b, port_5x5[,(i*25-24):(i*25)])
}

