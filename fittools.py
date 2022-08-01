import pandas as pd
import numpy as np
from scipy import stats
import glob
import math
import sys

import csv


def getStatus(df,countThreshold,n_nuisances):
	totalMCpre=0
	totalMCpost=0
	totalD=0
	totalMCpre = df['bprefit'].sum()
	totalMCpost = df['bpostfit'].sum()
	totalD = df['data'].sum()
	print("Prefit: Total Expected:",round(totalMCpre),"Total Observed:",totalD,"O-E=",round(totalD-totalMCpre),"nuisances= 0")
	print("Postfit: Total Expected:",round(totalMCpost),"Total Observed:",totalD,"O-E=",round(totalD-totalMCpost),"nuisances= ",n_nuisances)
    
def getChisq(df, countThreshold,n_nuisances):
	chidf = pd.DataFrame()
	chidf['observed'] = df['data']
	chidf['expected_pre'] = df['bprefit']
	chidf['expected_post'] = df['bpostfit']
	chidf['variance_pre'] = df['bprefit_err']
	chidf['variance_post'] = df['bpostfit_err']
	chidf['variance_pre'] = chidf['variance_pre'].apply(lambda x:x*x)
	chidf['variance_post'] = chidf['variance_post'].apply(lambda x:x*x)
    
    
    #chipre = chidf.loc[ chidf['expected_pre'] > (10e-5) ]
	TotalBins = chidf.shape[0]
	chipre = chidf.loc[ chidf['expected_pre'] > countThreshold]
	chipost = chidf.loc[ chidf['expected_post'] > countThreshold]
	SelectedBins_pre = chipre.shape[0]
	SelectedBins_post = chipost.shape[0]
    
    #calculate columns for pre chisq
	chipre['diff'] = chipre['observed'] - chipre['expected_pre']
	#chipre['sqdiff'] = chipre.apply(lambda x:(x*x))
	chipre['sqdiff'] = chipre['diff']*chipre['diff']
	#chipre['ratio'] = chipre['sqdiff']/chipre['variance_pre']
	chipre['ratio'] = chipre['sqdiff']/chipre['expected_pre']
	chi2pre = chipre['ratio'].sum()
	chi2preNDF = chi2pre/(float(SelectedBins_pre))
	pvalue_pre = 1- stats.chi2.cdf(chi2pre,SelectedBins_pre)
    
    #calculate columns for post chisq
	chipost['diff'] = chipost['observed'] - chipost['expected_post']
	#chipost['sqdiff'] = chipost.apply(lambda x:(x*x))
	chipost['sqdiff'] = chipost['diff']*chipost['diff']
	#chipost['ratio'] = chipost['sqdiff']/(chidf['expected_post'] - chidf['variance_post'] )
	chipost['ratio'] = chipost['sqdiff']/chipost['expected_post']
	#chipost['ratio'] = chipost['sqdiff']/chidf['variance_post']
	chi2post = chipost['ratio'].sum()
	chi2postNDF = chi2post/(float(SelectedBins_post))
	pvalue_post = 1- stats.chi2.cdf(chi2post,SelectedBins_post)
    
	print("Prefit Chi2  (O-E)^2 / E")
	print('{0: <12}'.format("Chisq"), '{0: <12}'.format("Chisq/NDF"), '{0: <12}'.format("P-val"),end='')
	print('{0: <12}'.format("Selected Bins"), '{0: <12}'.format("Total Bins"))
	print('{0: <12}'.format(round(chi2pre,1)), '{0: <12}'.format(round(chi2preNDF,1)), '{0: <12}'.format(round(pvalue_pre,2)), '{0: <12}'.format(SelectedBins_pre), '{0: <12}'.format(TotalBins))
	print("Postfit Chi2 (O-E)^2 / E")
	print('{0: <12}'.format("Chisq"), '{0: <12}'.format("Chisq/NDF"), '{0: <12}'.format("P-val"),end='')
	print('{0: <12}'.format("Selected Bins"), '{0: <12}'.format("Total Bins"))
	print('{0: <12}'.format(round(chi2post,1)), '{0: <12}'.format(round(chi2postNDF,1)), '{0: <12}'.format(round(pvalue_post,2)), '{0: <12}'.format(SelectedBins_post), '{0: <12}'.format(TotalBins))

def getPLike(df):
	poissondf = pd.DataFrame()
	poissondf['observed'] = df['data']
	poissondf['expected'] = df['bprefit']
	poissondf['loglam'] = np.log(poissondf['expected'])
	poissondf['kloglam'] = (-1)* poissondf['observed']*poissondf['loglam']    
	poissondf['logL'] =  poissondf['kloglam'] + poissondf['expected']
	sumLogLpre = poissondf['logL'].sum()
	
	poissondf = pd.DataFrame()
	poissondf['observed'] = df['data']
	poissondf['expected'] = df['bpostfit']
	poissondf['loglam'] = np.log(poissondf['expected'])
	poissondf['kloglam'] = (-1)* poissondf['observed']*poissondf['loglam']    
	poissondf['logL'] =  poissondf['kloglam'] + poissondf['expected']
	sumLogLpost = poissondf['logL'].sum()
	
	print("Prefit -LogLikelihood")
	print('{0: <12}'.format(round(-sumLogLpre,1)))
	print("Postfit -LogLikelihood")
	print('{0: <12}'.format(round(-sumLogLpost,1)))
	
def getPull(df):
	res = df
	res = res.loc[ res['bprefit'] > (10e-5) ]#remove extraneous bins created by diagnostic 
	res['Pre_O-E/sqrt(E)'] = (res['data'] - res['bprefit'])/ res['bprefit']**(1./2.)
	res['Pull_O-E/sqrt(E-VF)'] = (res['data'] - res['bpostfit'])/ ( res['bpostfit'] - res['bpostfit_err']**2 )**(1./2.)
	res = res[['RegionName','BinNumber','Pre_O-E/sqrt(E)','Pull_O-E/sqrt(E-VF)','data','bprefit','bpostfit']]
    
	print('Top 10 Leading normalized residuals sorted by pull') 
	res = res.sort_values(by='Pull_O-E/sqrt(E-VF)', key=pd.Series.abs, ascending=False)
	print(res[:10])
	print('Top 10 Leading normalized residuals sorted by pre')    
	res = res.sort_values(by='Pre_O-E/sqrt(E)',key=pd.Series.abs, ascending=False)
	print(res[:10])
    
def analyzedf(df, countThreshold, n_nuisances):
	getStatus(df,countThreshold,n_nuisances)
	getChisq(df,countThreshold,n_nuisances)
	getPLike(df)
	getPull(df)
    
