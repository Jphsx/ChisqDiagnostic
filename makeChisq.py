import pandas as pd
import numpy as np
from scipy import stats
import glob
import math
import sys

import csv
#import ROOT

#filename = "prepost.csv"
#filename = "pre18.csv"
#filename = "M4_27_22_allyearFD.csv"
#filename = "M5_31_22_2016FD.csv"
#filename = "M_6-2-22_nobtag2016.csv"
#filename = "M_6-2-22_bshape2016.csv"
#filename = "M_6-2-22_bnorm2016.csv"
#filename = "B19-5.csv"
#filename = "B19-6.csv"
#filename = "B19-8.csv"
#filename = "B19-9.csv"
#filename = "B20-1.csv"
#filename = "B20-2.csv"
#filename = "B20-1A.csv"
#filename = "M_6-27.csv"
filename = "B_20-3.csv"
#filename = "B_20-7.csv"
#filename = "B_20-7_ZFD.csv"
#filename = str( sys.argv[1] )
pd.set_option('display.max_colwidth',None)
pd.set_option('display.max_columns', None)
pd.set_option('expand_frame_repr', False)
pd.set_option('display.max_rows', None)

df = pd.read_csv(filename, index_col=None, header=0, delimiter=" ")
#print(df)



#create a column with a list of criteria strings
df['RegionSplit'] = df['RegionName'].str.split('_')
#create a column that stores the length of the criteria string list
df['CriteriaLen'] = df['RegionSplit'].apply(lambda x:len(x))

#print(df)
countThreshold = 5 # minimum events allowed in Chisq calculation
#exit()
#pass in the columns you need to calculate chi2

def reportChisq_mc(data,mc,mc_err,countThreshold,nuisances): #run w.r.t MC error
    #data['observed'] = data[(data.columns)[0]]
    chidf = pd.DataFrame()
    chidf['observed'] = data
    chidf['expected'] = mc
    chidf['variance'] = mc_err.apply(lambda x:(x*x)) 
    #define error= sqrt( VE + VEPoisson)
    chidf['variance'] = chidf['variance'] + chidf['expected']
    chidf['diff'] = chidf['observed'] - chidf['expected']
    chidf['sqdiff'] = chidf['diff'].apply(lambda x:(x*x))
    chidf['ratio'] = chidf['sqdiff']/chidf['variance']
    #select data within threshold
    #remove extraneous bins
    chidf = chidf.loc[ chidf['expected'] > (10e-5)  ]
    TotalBins = chidf.shape[0]
    chidf = chidf.loc[ chidf['expected'] > countThreshold]
    SelectedBins = chidf.shape[0]
    chi2 = chidf['ratio'].sum() 
    chi2ndf = chi2/(float(SelectedBins))
    pvalue = 1- stats.chi2.cdf(chi2,SelectedBins)
    #round everything 
    print('{0: <12}'.format(round(chi2,1)), '{0: <12}'.format(round(chi2ndf,1)), '{0: <12}'.format(round(pvalue,2)), '{0: <12}'.format(SelectedBins), '{0: <12}'.format(TotalBins))

def reportChisq_data(data,data_err,mc,countThreshold): #run w.r.t data error
    reportChisq_mc(data,mc,data_err,countThreshold)

def report_PLike(data, mc):
    poissondf = pd.DataFrame()
    poissondf['observed'] = data
    poissondf['expected'] = mc
    poissondf['loglam'] = np.log(poissondf['expected'])
    poissondf['kloglam'] = (-1)* poissondf['observed']*poissondf['loglam']    
    poissondf['logL'] =  poissondf['kloglam'] + poissondf['expected']
    sumLogL = poissondf['logL'].sum()
    print('{0: <12}'.format(round(-sumLogL,1)))

def getResiduals(RegionName,BinNumber,data,data_error,expected,expected_error,countThreshold,Mode=1,pltname=0, hist=0) :
    res = pd.DataFrame()
    res['RegionName'] = RegionName
    res['BinNumber'] = BinNumber
    res['id'] = res['RegionName']+res['BinNumber'].astype(str)
    res['observed'] = data
    res['expected'] = expected
    res['data_err'] = data_error
    res['exp_err'] = expected_error
    #res = res.loc[ res['observed'] <= countThreshold]
    res = res.loc[ res['expected'] <= countThreshold]
    res = res.loc[ res['expected'] > (10e-5) ]#remove extraneous bins created by diagnostic 
    res['diff'] = res['observed'] - res['expected']
    res['O-E/E'] = res['diff']/res['exp_err']
    res['O-E/O'] = res['diff']/res['data_err']
    res['exp_var'] = res['exp_err'].apply(lambda x:(x*x))
    res['exp_var'] = res['exp_var'] + res['expected']
    res['obs_var'] = res['data_err'].apply(lambda x:(x*x))
    #res['O+E'] = res['exp_var'] + res['obs_var']
    res['O+E'] = res['exp_var'].apply(lambda x:math.sqrt(x))
    res['O-E/sqrt(VE)'] = res['diff']/res['O+E']
    res['absresE'] = res['O-E/E'].abs()
    res['absresO'] = res['O-E/O'].abs()
    res['absresVE'] = res['O-E/sqrt(VE)'].abs()
    #if(sortMode == 0):
    res = res.sort_values(by='absresVE',ascending=False)
    #if(sortMode == 1):
    #    res = res.sort_values(by='absresO',ascending=False)
    res = res[['RegionName','BinNumber','O-E/sqrt(VE)','observed','expected']]
    if(Mode ==1):
        print(res[:20])

    if(Mode == 2):
        #res['Global Bin Number'] = res.index
        #plt = res.plot.scatter(x='Global Bin Number',y='O-E/sqrt(VE)')
        #plt.set_ylim(-5.7,5.7)
        #plt = plt.get_figure()
        #plt.savefig(pltname)
        #print(res)
        #hist = rt.TH1D(pltname,pltname,40,-6,6)
        for pull in res['O-E/sqrt(VE)'].to_numpy():
        	hist.Fill(pull)
        	#print(pull)
        #hist.Draw()
        #print(res)
		

    	

def getStatus(RegionName,BinNumber,data,expected,countThreshold,n_nuisances):
    totalMC=0
    totalD=0
    status = pd.DataFrame()
    status['RegionName'] = RegionName
    status['BinNumber'] = BinNumber
    status['observed'] = data
    status['expected'] = expected
    totalMC = status['expected'].sum()
    totalD = status['observed'].sum()
    print("countThreshold=",countThreshold,"Total Expected:",round(totalMC),"Total Observed:",totalD,"O-E=",round(totalD-totalMC),"nuisances=", n_nuisances)

def makePlot(RegionName,BinNumber,data,data_error,expected,expected_error,countThreshold,Mode=1,pltname=0, hist=0):
    getResiduals(RegionName,BinNumber,data,data_error,expected,expected_error,countThreshold,Mode,pltname,hist)


	

def evaluate_Chisqs(df,mode,countThreshold,n_nuisances):
    if(mode==0):
        print('{0:<25}'.format(" "),'{0: <12}'.format("Chisq"), '{0: <12}'.format("Chisq/NDF"), '{0: <12}'.format("P-val"),end='')
        print('{0: <12}'.format("Selected Bins"), '{0: <12}'.format("Total Bins"))
        print('{0: <25}'.format("Background only prefit"),end='')
        reportChisq_mc(df['data'],df['bprefit'],df['bprefit_err'],countThreshold,0)
        print('{0: <25}'.format("Background only postfit"),end='')
        reportChisq_mc(df['data'],df['bpostfit'],df['bpostfit_err'],countThreshold,n_nuisances)
        #do poisson stuff
        print('{0:<25}'.format(" "),'{0: <12}'.format("-LogLikelihood"))
        print('{0: <25}'.format("Background only prefit"),end='')
        report_PLike(df['data'],df['bprefit'])
        print('{0: <25}'.format("Background only postfit"),end='')
        report_PLike(df['data'],df['bpostfit'])
     #   print('{0: <25}'.format("S+B prefit"),end='')
     #   reportChisq_mc(df['data'],df['sbprefit'],df['sbprefit_err'],countThreshold)
     #   print('{0: <25}'.format("S+B postfit"),end='')
     #   reportChisq_mc(df['data'],df['sbpostfit'],df['sbpostfit_err'],countThreshold)
    #if(mode==1):#removing this for now
    #    reportChisq_mc(df['data'],df['data_err'],df['bprefit'],countThreshold)
    #    reportChisq_mc(df['data'],df['data_err'],df['bpostfit'],countThreshold)
    #    reportChisq_mc(df['data'],df['data_err'],df['sbprefit'],countThreshold)
    #    reportChisq_mc(df['data'],df['data_err'],df['sbpostfit'],countThreshold)
    print("Full status Prefit")
    getStatus(df['RegionName'],df['BinNumber'],df['data'],df['bprefit'],countThreshold,0)
    print("Full status Postfit")
    getStatus(df['RegionName'],df['BinNumber'],df['data'],df['bpostfit'],countThreshold,n_nuisances)
    print('Top 20 Leading Residuals w.r.t bpostfit')
    getResiduals(df['RegionName'],df['BinNumber'],df['data'],df['data_err'],df['bpostfit'],df['bpostfit_err'],999999999,1)
    print('Top 20 Leading Low Stat Residuals w.r.t bpostfit')
    getResiduals(df['RegionName'],df['BinNumber'],df['data'],df['data_err'],df['bpostfit'],df['bpostfit_err'],countThreshold)
    print('Top 20 Leading Residuals w.r.t bprefit')
    getResiduals(df['RegionName'],df['BinNumber'],df['data'],df['data_err'],df['bprefit'],df['bprefit_err'],999999999,1)
    print('Top 20 Leading Low Stat Residuals w.r.t bprefit')
    getResiduals(df['RegionName'],df['BinNumber'],df['data'],df['data_err'],df['bprefit'],df['bprefit_err'],countThreshold)

    #generating plots if enabled
    #makePlot(pltname,prefit,postfit,

#####################################################################start
#DEFINE N NUISANCE
n_nuisances=186
#####

print("------------------Evaluating all bins------------------")
evaluate_Chisqs(df,0,countThreshold,n_nuisances)
df['RegionName'],df['BinNumber'],df['data'],df['data_err'],df['bprefit'],df['bprefit_err'],countThreshold


#make a 1d plot of pulls prefit and postfit
import ROOT as rt
#outfile = rt.TFile("canvasTemp.root","RECREATE")
prehist = rt.TH1D("prefitPull","prefitPull",40,-6,6)
makePlot(df['RegionName'],df['BinNumber'],df['data'],df['data_err'],df['bprefit'],df['bprefit_err'],99999999,2,"prefitPull",prehist)
#prehist.Draw()
posthist = rt.TH1D("postfitPull","postfitPull",40,-6,6)

#posthist.SetLineColor(rt.kRed)
posthist.SetLineWidth(2)
posthist.SetTitle("Post-fit Pull per Analysis Bin ;(Expected-Observed)/#sigma_{Expected}; N bins")
makePlot(df['RegionName'],df['BinNumber'],df['data'],df['data_err'],df['bpostfit'],df['bpostfit_err'],99999999,2,"all2Lpost.pdf",posthist)


#c1 = rt.TCanvas("c","c")
posthist.Draw("E")
posthist.Fit("gaus")
#fit = posthist.GetFunction("gaus")
#posthist.Draw("SAME")
#fit.SetLineColor(rt.kRed)
#fit.Draw("SAME")
#outfile.WriteTObject(c1)

#outfile.Close()



#fit.Draw()
#slice up our df and look at smaller pieces
#create rank column
df['Rank'] = df['RegionSplit'].apply( lambda x:x[2] )
df['NLep'] = df['RegionSplit'].apply( lambda x:x[0] )

print('\n')
print("------------------all 0L bins------------------")
L0 = df.loc[ df['NLep'] == 'Ch0L' ]
evaluate_Chisqs(L0,0,countThreshold,n_nuisances)

print('\n')
print("------------------all 1L bins------------------")
L1 = df.loc[ df['NLep'] == 'Ch1L' ]
evaluate_Chisqs(L1,0,countThreshold,n_nuisances)


print('\n')
print("------------------all 2L bins------------------")
L2 = df.loc[ df['NLep'] == 'Ch2L' ]
evaluate_Chisqs(L2,0,countThreshold,n_nuisances)


print('\n')
print("------------------all 3L bins------------------")
L3 = df.loc[ df['NLep'] == 'Ch3L' ]
evaluate_Chisqs(L3,0,countThreshold,n_nuisances)


print('\n')
print("------------------Evaluating all gold bins------------------")
dfgold = df.loc[ df['Rank'] == 'gold' ]
evaluate_Chisqs(dfgold,0,countThreshold,n_nuisances)

print('\n')
print("Evaluating all silver bins")
dfslvr = df.loc[ df['Rank'] == 'slvr' ]
evaluate_Chisqs(dfslvr,0,countThreshold,n_nuisances)

print('\n')
print("Evaluating all bronze bins")
dfbron = df.loc[ df['Rank'] == 'bron' ]
evaluate_Chisqs(dfbron,0,countThreshold,n_nuisances)




