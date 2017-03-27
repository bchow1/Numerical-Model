#!/usr/bin/python
import os
import sys
import socket
import time
import numpy as np
#import numpy.ma as ma
import sqlite3
import matplotlib.pyplot as plt
import fileinput
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import spline

compName = socket.gethostname()
env      = os.environ.copy()

print compName

# Local modules
if compName == 'sm-bnc' or compName == 'sage-d600':
  sys.path.append('C:\\Users\\sid\\python')
  sys.path.append('D:\\TestSCICHEM\\Scripts')
  sys.path.append('D:\\TestSCICHEM\\Scripts\\Chemistry')
  sys.path.append('D:\\TestSCICHEM\\Scripts\\AERMOD')

if  compName == 'pj-linux4':
  sys.path.append('/home/user/bnc/python')
  sys.path.append('/home/user/bnc/TestSCICHEM/Scripts')
  sys.path.append('/home/user/bnc/TestSCICHEM/Scripts/Chemistry')
  sys.path.append('/home/user/bnc/TestSCICHEM/Scripts/AERMOD')

if compName == 'Durga':
  sys.path.append('C:\\Users\\Bishusunita\\BNC\\python')
  sys.path.append('C:\\Users\\Bishusunita\\BNC\\TestSCICHEM\\Scripts')
  sys.path.append('C:\\Users\\Bishusunita\\BNC\\TestSCICHEM\\Scripts\\Chemistry')
  sys.path.append('C:\\Users\\Bishusunita\\BNC\\TestSCICHEM\\Scripts\\AERMOD')

if compName == 'sm-bnc':  
  sys.path.append('C:\\Users\\sid\\python')
  sys.path.append('D:\\TestSCICHEM\\Scripts')
  sys.path.append('D:\\TestSCICHEM\\Scripts\\Chemistry')
  sys.path.append('D:\\TestSCICHEM\\Scripts\\AERMOD')
  
# Local libs
import utilDb
import setSCIparams as SCI
import run_cmd
import measure

if compName == 'sm-bnc':
  os.chdir('D:\\SCIPUFF\\runs\\EPRI\\DoletHill')
if compName == 'Durga':
  os.chdir('C:\\Users\\Bishusunita\\BNC\\Publications\\SCICHEM\\AtmEnv\\DoletHills')
if compName == 'pj-linux4':
  #os.chdir('/home/user/bnc/scipuff/runs/EPRI/DoletHill')
  os.chdir('/home/user/bnc/EPRI/SCICHEM3.0B3/PPTCases/DoletHills')

isColor = True
isTotal = True
print os.getcwd()
samFile = 't1345678.sam' 
if sys.argv.__len__() > 1:
  prjName = sys.argv[1]
  ans = raw_input('Use sampler file for project %s? (y/Y): '%prjName)
  if ans[0].upper() != 'Y':
    sys.exit()
  utilDb.createSmpDb([prjName],samFiles=[samFile])
  preDbName = '%s.smp.db'%prjName
else:
  preDbName = 'dolethills_ppm.smp.db'

mySciFiles = SCI.Files(preDbName)
(preConn,preCur,preCreateDb) = utilDb.Smp2Db(preDbName,mySciFiles)
  
obsDbName = 'aztec_20050908_some.csv.db'
obsConn = sqlite3.connect(obsDbName)
obsConn.row_factory = sqlite3.Row
obsCur = obsConn.cursor()

# Load all samplers from all traverses

allSmps = np.loadtxt(samFile,skiprows=1,usecols=(0,1,2))

# Number of sampler in each traverse
smpNos   = {1:30, 3:54, 4:45, 5:17, 6:94, 7:82, 8:210}

# Distance to midpoint of the traverses
# Only using values for traverse 3,6 and 8
srcLoc = [94.839,-38.8915]
trDist = {1:[111.,222.],3:[80.2538,-50.2851],4:[111.,222.],5:[111.,222.],\
          6:[50.0423,-72.4951],7:[111.,222.],8:[6.6158,-80.5116]}
figNo  = {1:11, 3:12, 4:12, 5:12, 6:13, 7:13, 8:14}

def createArraysAndPlots(preDbName,samFile):

  preConn = sqlite3.connect(preDbName)
  preConn.row_factory = sqlite3.Row
  preCur = preConn.cursor()
    
  obsDbName = 'aztec_20050908_some.csv.db'
  obsConn = sqlite3.connect(obsDbName)
  obsConn.row_factory = sqlite3.Row
  obsCur = obsConn.cursor()
  
  # Load all samplers from all traverses
  
  allSmps = np.loadtxt(samFile,skiprows=1,usecols=(0,1,2))
  
  # Number of sampler in each traverse
  smpNos   = {1:30, 3:54, 4:45, 5:17, 6:94, 7:82, 8:210}
  varNames = ['NOx','NOy','SO2','O3'] 
  
  if 'ppm' in preDbName:
    fug2ppm = [1.,1., 1., 1.] # SO2,O3,NOx,NOy
    print 'Using ppm conc. units'
  else:
    fug2ppm = [0.075/196.,0.075/147., 0.053/100., 0.053/100.] # SO2,O3,NOx,NOy
    print 'Using ug/m3 conc. units'
  time.sleep(1)
  
  trList = [1,3,6,8] #[1,3,4,5,6,7,8] #Note: Need 1 to create TrTbl
  
  locDict = {}

  print 'Create TrTbl'
  ans = raw_input('Press y|Y to continue')

  # Set locList and create TrTbl for all traverses
  for trNo in trList:      
      # Increment starting row
      startNo = 0 
      for trn,trl in smpNos.items():    
        endNo = startNo + trl - 1   
        if trn == trNo:
          break 
        startNo = endNo + 1
        
      print 'For traverse ',trNo,', # samplers = ',smpNos[trNo]
      print 'Start = ',startNo,allSmps[startNo,0],allSmps[startNo,1]
      print 'End   = ',endNo,allSmps[endNo,0],allSmps[endNo,0]
         
      locList = []    
      for allNo in range(startNo,endNo):
    
        #print allNo, allSmps[allNo,:]
        xSmp = allSmps[allNo,0]
        ySmp = allSmps[allNo,1]
        #print xSmp,ySmp
        locList.append([xSmp,ySmp])
        
        if False:
          # plot location
          plt.clf()
          plt.hold(True)
          tSmps = np.loadtxt('Alltraverse%d.sam'%trNo,skiprows=1,usecols=(0,1))
          locA = np.array(locList)
          plt.plot(locA[:,0],locA[:,1],color='red',marker='o')
          plt.plot(tSmps[:,0],tSmps[:,1]+.1,color='green',marker='+')
          plt.xlim([-95.8,-93.3])
          plt.ylim([31.6,33.16])
          plt.hold(False)
          plt.savefig('Traverse%d_loc.png'%trNo)
          
      locDict.update({trNo:locList})
      
      # Create TrTbl
      if ans == 'y' or ans == 'Y':
        createTrTbl(obsCur,preCur,trNo,locList)
        obsConn.commit()
      else:
        print 'Skip creating TrTbl'
            
  #print locDict
    
  # Loop over the traverses
  for trNo in trList:
    
    fnameStr = str(trNo) + "test_stat.csv"
    statFile = open(fnameStr,"w")
    #Normalized Mean Bias nmb, Mean Bias Error mbe, Mean Normalized Bias Error mnbe, Mean Normalized Gross Error mnge
    statFile.write("TrNo, Distance, Species, Mean Observed (ppb), Mean Predicted(ppb), NMSE, MFB, r, "\
                   " NMB, MBE(ppb)\n")    
    
    obsQry = "Select * From TrTbl where trNo = %d"%trNo
    #  0     1   2   3   4    5   6  7   8      9      10     11     12      13      14     15    16
    # trNo xobs yobs tHr NOx NOy O3 SO2 preNO2 preNO preNO3 preN2O5 preHNO3 preHONO prePAN preO3 preSO2
    obsList = utilDb.db2List(obsCur,obsQry) 
    #for obs in obsList:
    #  print obs
      
    obsArray = np.array(obsList)
    xObs    = obsArray[:,1]
    yObs    = obsArray[:,2]
    tObs    = obsArray[:,3]
    NOx_obs = obsArray[:,4]
    NOy_obs = obsArray[:,5]
    SO2_obs = obsArray[:,6] 
    O3_obs  = obsArray[:,7]
    NOx_pre = obsArray[:,8] + obsArray[:,9]
    NOy_pre = NOx_pre +  obsArray[:,10] + obsArray[:,11] + obsArray[:,12] + obsArray[:,13]  + obsArray[:,14]
    O3_pre =  obsArray[:,15] 
    SO2_pre = obsArray[:,16]   
    #var_pre  = ((tpre2-tHr)*var_pre1 + (tHr-tpre1)*var_pre2)/(tpre2-tpre1)
   
    dist = np.sqrt( (trDist[trNo][0] - srcLoc[0])**2 + \
                    (trDist[trNo][1] - srcLoc[1])**2 )
    print '\nGetting stats for Traverse %d at %d km'%(trNo,dist)
    statPrx = '%d , %d '%(trNo,dist)
    
    # NOx 
    mask = NOx_obs < -9900
    NOx_obs = np.ma.array(NOx_obs,mask=mask).compressed()
    NOx_pre = np.ma.array(NOx_pre,mask=mask).compressed()
    NOx_obs[NOx_obs < 0.] = 0.
    NOx_pre = NOx_pre*fug2ppm[1]*1e3  # ppb
    if isTotal:
      NOx_bgd = min(NOx_obs[:3].max(),NOx_obs[-3:].max())
      NOx_pbgd = min(NOx_pre[:3].max(),NOx_pre[-3:].max()) 
      print '\nNOx',NOx_bgd,NOx_pbgd
      NOx_diff = NOx_bgd - NOx_pbgd
      calcStats(NOx_obs,NOx_pre+NOx_diff,'NOx',statFile=statFile,statPrx=statPrx )
      NOx_obs = NOx_obs - NOx_bgd
      NOx_obs[NOx_obs < 0.] = 0.
      NOx_pre = NOx_pre - NOx_pbgd
      NOx_pre[NOx_pre < 0.] = 0.
    else:
      NOx_bgd = min(NOx_obs[:3].min(),NOx_obs[-3:].min())
      print '\nNOx',NOx_bgd
      calcStats(NOx_obs,NOx_pre+NOx_bgd,'NOx',statFile=statFile,statPrx=statPrx )
      NOx_obs = NOx_obs - NOx_bgd
    xNOx    = np.ma.array(xObs,mask=mask).compressed()
    

    # NOy
    mask = NOy_obs < -9900
    NOy_obs = np.ma.array(NOy_obs,mask=mask).compressed()
    NOy_pre = np.ma.array(NOy_pre,mask=mask).compressed()
    xNOy    = np.ma.array(xObs,mask=mask).compressed()
    NOy_obs[NOy_obs < 0.] = 0.
    NOy_pre = NOy_pre*fug2ppm[1]*1e3 # ppb
    if isTotal:
      NOy_bgd = min(NOy_obs[:3].max(),NOy_obs[-3:].max())
      NOy_pbgd = min(NOy_pre[:3].max(),NOy_pre[-3:].max())
      print '\nNOy',NOy_bgd,NOy_pbgd
      NOy_diff = NOy_bgd - NOy_pbgd
      calcStats(NOy_obs,NOy_pre+NOy_diff,'NOy',statFile=statFile,statPrx=statPrx )
      NOy_obs = NOy_obs - NOy_bgd
      NOy_obs[NOy_obs < 0.] = 0.
      NOy_pre = NOy_pre - NOy_pbgd
      NOy_pre[NOy_pre < 0.] = 0.
    else:
      NOy_bgd = min(NOy_obs[:3].min(),NOy_obs[-3:].min())
      print '\nNOy',NOy_bgd
      calcStats(NOy_obs,NOy_pre + NOy_bgd,'NOy',statFile=statFile,statPrx=statPrx )
      NOy_obs = NOy_obs - NOy_bgd
    
    # O3 
    mask = O3_obs < -9900
    O3_obs = np.ma.array(O3_obs,mask=mask).compressed()
    O3_pre = np.ma.array(O3_pre,mask=mask).compressed()
    xO3    = np.ma.array(xObs,mask=mask).compressed()
    O3_obs[O3_obs < 0.] = 0.
    O3_pre = O3_pre*fug2ppm[1]*1e3  # ppb
    if isTotal:
      O3_bgd = min(O3_obs[:3].min(),O3_obs[-3:].min())
      O3_pbgd = min(O3_pre[:3].min(),O3_pre[-3:].min()) 
      print '\nO3',O3_bgd,O3_pbgd
      O3_diff = O3_bgd - O3_pbgd
      calcStats(O3_obs,O3_pre+O3_diff,'O3',statFile=statFile,statPrx=statPrx )
      O3_obs = O3_obs - O3_bgd
      O3_obs[O3_obs < 0.] = 0.
      O3_pre = O3_pre - O3_pbgd
      O3_pre[O3_pre < 0.] = 0.
    else:
      O3_bgd = min(O3_obs[:3].min(),O3_obs[-3:].min())
      print '\nO3',O3_bgd
      calcStats(O3_obs,O3_pre +  O3_bgd,'O3',statFile=statFile,statPrx=statPrx )
      O3_obs = O3_obs - O3_bgd
    
    
    # SO2
    mask = SO2_obs < -9900
    SO2_obs = np.ma.array(SO2_obs,mask=mask).compressed()
    SO2_pre = np.ma.array(SO2_pre,mask=mask).compressed()
    xSO2    = np.ma.array(xObs,mask=mask).compressed()
    SO2_obs[SO2_obs < 0.] = 0.
    SO2_pre = SO2_pre*fug2ppm[1]*1e3 # ppb
    if isTotal:
      SO2_bgd = 0.
      SO2_bgd = min(SO2_obs[:3].min(),SO2_obs[-3:].min())
      SO2_pbgd = min(SO2_pre[:3].min(),SO2_pre[-3:].min()) 
      print '\nSO2',SO2_bgd,SO2_pbgd
      SO2_diff = SO2_bgd - SO2_pbgd
      calcStats(SO2_obs,SO2_pre+SO2_diff,'SO2',statFile=statFile,statPrx=statPrx )
      SO2_obs = SO2_obs - SO2_bgd
      SO2_obs[SO2_obs < 0.] = 0.
      SO2_pre = SO2_pre - SO2_pbgd
      SO2_pre[SO2_pre < 0.] = 0.
    else:
      SO2_bgd = min(SO2_obs[:3].min(),SO2_obs[-3:].min())
      print '\nSO2',SO2_bgd
      calcStats(SO2_obs,SO2_pre +  SO2_bgd,'SO2',statFile=statFile,statPrx=statPrx )
      SO2_obs = SO2_obs - SO2_bgd
        
    #print NOx_obs
    #print NOx_pre
    params1 = {'axes.labelsize': 10, 'text.fontsize': 10, 'xtick.labelsize': 10,
                'ytick.labelsize': 10, 'legend.pad': 0.1,  
                'legend.fontsize': 8, 'lines.markersize': 6, 'lines.width': 2.0,
                'font.size': 10, 'text.usetex': False}
    plt.rcParams.update(params1)
    
    fig = plt.figure(1)
    plt.clf()
    plt.hold(True)
    
    fig.subplots_adjust(bottom=0.2)
    fig.subplots_adjust(top=0.9)
    figCaption1 = ''#'Figure %d:  Comparison of the mean observed and predicted concentrations for Dolet Hills plume '%figNo[trNo]
    figCaption2 = ''#'traverse study at %d Km for Traverse %s.'%(dist,trNo)                 
    plt.text(0.45,-0.09,'Longitude',fontsize=10)
    plt.text(-0.04,-0.16,figCaption1,fontsize=10)
    plt.text(-0.04,-0.20,figCaption2,fontsize=10)
    plt.text(-.09,0.7,'Concentration Perturbation (ppb)',rotation='vertical',fontsize=10)
    
    plt.setp(plt.gca(), frame_on=False, xticks=(), yticks=())
    
    xlim = None
    if trNo == 3:
      xlim = [93.55,93.85]
    if trNo == 6:
      xlim = [93.85,94.2]
    if trNo == 8:
      xlim = [93.9,94.9]
    
    print 'xlim = ',xlim
    msize = 5
    lsty  = '-'
    msty  = '^'
    if isColor:
      ObsClr = 'red'
      PreClr = 'green'
    else:
      ObsClr = '0.25'
      PreClr = '0.5'
    
    #NOx            
    ax = fig.add_subplot(2,2,1)
    plt.text(0.07,0.85,'NOx',transform=ax.transAxes)

    LhO, = plt.plot(xNOx,NOx_obs,linestyle='None',marker='o',markersize=msize, markevery=2, markerfacecolor=ObsClr)
    LhP, = plt.plot(xNOx,NOx_pre,linestyle=lsty,marker=msty,markersize=msize, markevery=4,color=PreClr,markerfacecolor=PreClr)

    ax.set_xticklabels([])
    if xlim is not None:
      ax.set_xlim(xlim)
    ax.xaxis.set_major_locator(MaxNLocator(6))
    
    #NOy
    ax = fig.add_subplot(2,2,2)
    plt.text(0.07,0.85,'NOy',transform=ax.transAxes)

    LhO, = plt.plot(xNOy,NOy_obs,linestyle='None',marker='o',markersize=msize, markevery=2,markerfacecolor=ObsClr)
    LhP, = plt.plot(xNOy,NOy_pre,linestyle=lsty,marker=msty,markersize=msize, markevery=4,color=PreClr,markerfacecolor=PreClr)

    ax.set_xticklabels([])
    if xlim is not None:
      ax.set_xlim(xlim)
    ax.xaxis.set_major_locator(MaxNLocator(6))
    
    #O3
    ax = fig.add_subplot(2,2,3)
    plt.text(0.07,0.85,'O3',transform=ax.transAxes)

    LhO, = plt.plot(xO3,O3_obs,linestyle='None',marker='o',markersize=msize, markevery=2,markerfacecolor=ObsClr)
    LhP, = plt.plot(xO3,O3_pre,linestyle=lsty,marker=msty,markersize=msize, markevery=4,color=PreClr,markerfacecolor=PreClr)

    if xlim is not None:
      ax.set_xlim(xlim)
    ax.xaxis.set_major_locator(MaxNLocator(6))
    
    #SO2
    ax = fig.add_subplot(2,2,4)
    plt.text(0.07,0.85,'SO2',transform=ax.transAxes)

    LhO, = plt.plot(xSO2,SO2_obs,linestyle='None',marker='o',markersize=msize, markevery=2,markerfacecolor=ObsClr)
    LhP, = plt.plot(xSO2,SO2_pre,linestyle=lsty,marker=msty,markersize=msize, markevery=4,color=PreClr,markerfacecolor=PreClr)

    if xlim is not None:
      ax.set_xlim(xlim)
    ax.xaxis.set_major_locator(MaxNLocator(6))
    fig.legend([LhO,LhP],['OBSERVED','SCICHEM-3.0'],ncol=2, bbox_to_anchor=(0.9,0.96))
    
    plt.hold(False)
    plt.savefig('Traverse%d.png'%(trNo))
    print 'Created Traverse%d.png in %s'%(trNo,os.getcwd())
           
  return

def calcStats(obsArray, preArray,varName,statFile=None,statPrx=''):
  if varName == 'O3':
    pass
  obsMean =  obsArray.mean()
  predMean = preArray.mean()
  print 'Min:',obsArray.min(),preArray.min()
  print 'Max:',obsArray.max(),preArray.max()
  if min(obsArray)< 0.:
    pass
  NMSE = measure.nmse_1(preArray,obsArray)
  #NME = measure.nme(preArray[:,1],obsArray[:,1], cutoff=0.0)
  #pearCoeff = pearsonr(obsArray,preArray[:,1])
  biasFac = measure.bf(preArray,obsArray)
  MFB = measure.mfbe(preArray,obsArray,cutoff=0.)
  mnbe = measure.mnbe(preArray,obsArray,cutoff=0.)
  correlation = measure.correlation(preArray,obsArray)
  mage = measure.mage(preArray,obsArray)
  nmb = measure.nmb(preArray,obsArray)
  mbe = measure.mbe(preArray,obsArray)
  rmse = measure.rmse(preArray,obsArray)
  mnge = measure.mnge(preArray,obsArray,cutoff=0.)
  nmaef = measure.nmaef(preArray,obsArray)
  nmbf  = measure.nmbf(preArray,obsArray)
  print 'stats: ',obsMean, predMean,';',NMSE, biasFac, correlation
  if statFile is not None:
    partStr = statPrx + "," + varName + ","
    statFile.write(partStr)
    statFile.write("%8.2f, %8.2f, %07.3f, %07.3f,%07.2f, %07.3f,%07.2f\n"%\
                  (obsMean, predMean, NMSE, MFB, correlation, nmb, mbe))
  
  return

# Create trTable if required 
def createTrTbl(obsCur,preCur,trNo,locList):

  if trNo == 1 :  
    obsCur.execute('DROP table if exists TrTbl')  
    createStr = 'CREATE table TrTbl (trNo integer, xobs real, yobs real, tHr real, '\
                'NOx real, NOy real, O3 real, SO2 real,  preNO2 real  ,preNO real '\
                '  ,preNO3 real   ,preN2O5 real   ,preHNO3 real   ,preHONO real  '\
                ' ,prePAN real   ,preO3 real   ,preSO2 real    )'
    print createStr
    obsCur.execute(createStr)    
 
  for xSmp,ySmp in locList:
        
    # Find xSmp and ySmp from observed db (aztec_20050908_some.csv.db)
    tHr = 0.
    
    initStr = "select  ABS(gps_lon_deg_5s) as xobs "
    initStr += " ,ABS(gps_lat_deg_5s) as yobs , date_DDMMYY, time_cdt , no_ppbv_5s + no2_ppbv_5s  "
    initStr += " ,noy_ppbv_5s, so2_ppbv_5s, o3_ppbv_5s "
    initStr += " from datatable where abs(xobs + %7.4f) < 0.0005 and abs(yobs - %7.4f) < 0.0005 "%( xSmp,ySmp)
    obsList1 = utilDb.db2List(obsCur,initStr) 
    tObs   = obsList1[0][3]
    tHr = 24. + run_cmd.hms2hr(tObs)
    xobs,yobs,date,time,nox,noy,so2,o3 = obsList1[0]
    
    #Now update the predicted values
    varQry =  'select time,xSmp,ySmp, varName, Avg(Value) from samTable a,smpTable p where a.colNo=p.colNo'
    varQry += " and varName in ('NO2','NO','NO3','N2O5','HNO3','HONO','PAN', 'O3','SO2')"
    preQry =  varQry + " and ABS(xSmp- %7.4f) < 0.0005 and ABS(ySmp- %7.4f) < 0.0005 "%(xSmp,ySmp)\
                     + " and ABS(time-%8.4f) < 0.165 group by varName "%(tHr)  
    print preQry
    preList = utilDb.db2List(preCur,preQry) 

    insertStr = "INSERT into TrTbl VALUES (%d,%f,%f,%8.2f,%8.3f,%8.3f,%8.3f,%8.3f,"% \
                (trNo,xobs,yobs,tHr,nox,noy,so2,o3)
    
    cPre = {}            
    for i in range(len(preList)): 
      cPre.update({preList[i][3]:preList[i][4]})
    print cPre
      
    for sp in ('NO2','NO','NO3','N2O5','HNO3','HONO','PAN', 'O3','SO2'): 
      insertStr += '%7.4f '%cPre[sp]
      if sp == 'SO2':
        insertStr += ')'
      else:                               
        insertStr += ','
    obsCur.execute(insertStr)
    print insertStr
  
  return
  
# Main program
if __name__ == '__main__':
  print 'Starting'
  createArraysAndPlots(preDbName,samFile)
  print 'Done:-)'
