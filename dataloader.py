#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DATALOADER
@author: victorsellemi
"""

#%% data import and management

# specify directories
main_dir = "/Users/victorsellemi/Downloads/DoubleLasso"
data_dir = main_dir + "/data"

#call directory
os.chdir(data_dir)                                       

#----------------------
# factor data
allfactors = pd.read_csv("factors.csv")
date    = allfactors.Date                                # date
rf      = allfactors.RF                                  # avg mkt return
factors = allfactors.iloc[:,2:-1]

L = len(date)
P = factors.shape[1]

#----------------------
# test portfolios
port_5x5 = pd.read_csv("port_5x5.csv", header = None).iloc[:,1:-1]
port_5x5 = port_5x5.sub(rf,axis=0)                      # excess returns

port_3x2 = pd.read_csv("port_3x2.csv", header = None).iloc[:,1:-1]
port_3x2 = port_3x2.sub(rf,axis=0)  

#port_5x5_seq = pd.read_csv("port_5x5_seq.csv", header = None).iloc[:,1:-1]
#port_5x5_seq = port_5x5_seq.sub(rf,axis=1)

#port_3x2_seq = pd.read_csv("port_3x2_seq.csv", header = None).iloc[:,1:-1]
#port_3x2_seq = port_3x2_seq.sub(rf,axis=1)  

port_202 = (pd.read_csv("port202.csv", header = None).iloc[:,1:-1]).div(100)
port_202 = port_202.sub(rf,axis=0)  

#---------------------
# other information
summary         = pd.read_csv("summary.csv")
factorname      = summary.Row
factorname_full = summary.Description
year_pub        = summary.Year
year_end        = summary.Year_end

port_5x5_id = pd.read_csv("port_5x5_id.csv")
port_3x2_id = pd.read_csv("port_3x2_id.csv")


#---------------------
# form a smaller set of portfolios for bivariate sorted porfolios
kk = 10                 # minimun number of stocks in a portfolio

include_3x2 = port_3x2_id.loc[port_3x2_id['min_stk6'] >= kk]
port_3x2b   = np.array([])

for i in range(P):
    if i in include_3x2.index:
        np.append(port_3x2, port_3x2.iloc[:,(i*6-5):(i*6)])

include_5x5 = port_5x5_id.loc[port_5x5_id['min_stk'] >= kk]
port_5x5b   = []

for i in range(P):
    if i in include_5x5.index:
        np.append(port_5x5, port_5x5.iloc[:,(i*6-5):(i*6)])