import pandas as pd
import numpy as np
from scipy import stats
import glob
import math
import sys
import csv
import ROOT as rt
import fittools as ft
pd.set_option('display.max_colwidth',None)
pd.set_option('display.max_columns', None)
pd.set_option('expand_frame_repr', False)
pd.set_option('display.max_rows', None)


##### Input FD csv file and important params#####
filename = "B_20-7.csv" #"B_20-7_ZFD.csv"#"B_20-3.csv"
countThreshold = 5 #Chisq min stat criteria
#n_nuisances=186 #hardcoded number of nuisances
#outfile = "B_20-7.root"#"B_20-7_ZFD.root"
filename = sys.argv[1]
outfile = sys.argv[3]
n_nuisances = float(sys.argv[2])
DO_POISSON= True
#############################



df = pd.read_csv(filename, index_col=None, header=0, delimiter=" ")
#df manipulation for g/s/b separation
df['RegionSplit'] = df['RegionName'].str.split('_')
df['CriteriaLen'] = df['RegionSplit'].apply(lambda x:len(x))
df['Rank'] = df['RegionSplit'].apply( lambda x:x[2] )#0L doesnt have rank, Njet will fill this criteria
df['NLep'] = df['RegionSplit'].apply( lambda x:x[0] )
#print(df)

histFile = rt.TFile(outfile,"RECREATE")
rt.gStyle.SetOptFit(1)

print("------------------Evaluating all bins------------------")
ft.analyzedf(df,countThreshold,n_nuisances, histFile, "all", DO_POISSON)

DO_POISSON= False

print('\n')
print("------------------all 0L bins------------------")
L0 = df.loc[ df['NLep'] == 'Ch0L' ]
ft.analyzedf(L0,countThreshold,n_nuisances, histFile, "0L", DO_POISSON)

print('\n')
print("------------------all 1L bins------------------")
L1 = df.loc[ df['NLep'] == 'Ch1L' ]
ft.analyzedf(L1,countThreshold,n_nuisances, histFile, "1L", DO_POISSON)

print('\n')
print("------------------all 2L bins------------------")
L2 = df.loc[ df['NLep'] == 'Ch2L' ]
ft.analyzedf(L2,countThreshold,n_nuisances, histFile, "2L", DO_POISSON)


print('\n')
print("------------------all 3L bins------------------", DO_POISSON)
L3 = df.loc[ df['NLep'] == 'Ch3L' ]
ft.analyzedf(L3,countThreshold,n_nuisances, histFile, "3L", DO_POISSON)


print('\n')
print("------------------Evaluating all gold bins------------------")
dfgold = df.loc[ df['Rank'] == 'gold' ]
ft.analyzedf(dfgold,countThreshold,n_nuisances, histFile,"gold", DO_POISSON)

print('\n')
print("Evaluating all silver bins")
dfslvr = df.loc[ df['Rank'] == 'slvr' ]
ft.analyzedf(dfslvr,countThreshold,n_nuisances, histFile, "slvr", DO_POISSON)

print('\n')
print("Evaluating all bronze bins")
dfbron = df.loc[ df['Rank'] == 'bron' ]
pulldf = ft.analyzedf(dfbron,countThreshold,n_nuisances, histFile, "bron", DO_POISSON)

#hack to evaluate bronze
#pull = pulldf.loc[ pulldf['Pre_O-E/sqrt(E)'] > 6 ]

#print(pull.shape)
#print(pull)
#histFile.Write()
histFile.Close()

