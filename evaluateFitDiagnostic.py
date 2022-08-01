import pandas as pd
import numpy as np
from scipy import stats
import glob
import math
import sys
import csv
import fittools as ft
pd.set_option('display.max_colwidth',None)
pd.set_option('display.max_columns', None)
pd.set_option('expand_frame_repr', False)
pd.set_option('display.max_rows', None)


##### Input FD csv file and important params#####
filename = "B_20-3.csv"
countThreshold = 5 #Chisq min stat criteria
n_nuisances=186 #hardcoded number of nuisances
#############################

df = pd.read_csv(filename, index_col=None, header=0, delimiter=" ")
#df manipulation for g/s/b separation
df['RegionSplit'] = df['RegionName'].str.split('_')
df['CriteriaLen'] = df['RegionSplit'].apply(lambda x:len(x))
df['Rank'] = df['RegionSplit'].apply( lambda x:x[2] )#0L doesnt have rank, Njet will fill this criteria
df['NLep'] = df['RegionSplit'].apply( lambda x:x[0] )
#print(df)

print("------------------Evaluating all bins------------------")
ft.analyzedf(df,countThreshold,n_nuisances)

