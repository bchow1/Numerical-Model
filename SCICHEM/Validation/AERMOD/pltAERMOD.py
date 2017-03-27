#!/bin/env python
import os
import sys
import shutil
import socket
import numpy as np
import sqlite3
import matplotlib.pyplot as plt
import pandas as pd

compName = socket.gethostname()
env      = os.environ.copy()

# Local modules
if compName == 'sm-bnc' or compName == 'sage-d600':
  sys.path.append('C:\\Users\\sid\\python')
  sys.path.append('D:\\TestSCICHEM\\Scripts')
if  compName == 'pj-linux4':
  sys.path.append('/home/user/bnc/python:/home/user/bnc/TestSCICHEM/Scripts')

# Local lib files
import setSCIparams as SCI
import run_cmd
import utilDb
import measure
import pltBarAermod

def getMaxConc(prjName,iHr=1):
  sciConn,sciCur = getSmpDb(prjName)
  smpFac         = 1.e9  # kg/m3 to ug/m3
  sciArray       = smpDbMax(sciCur,iHr,smpFac)
  return sciArray

def setAerModDir(valDir,ymdhm=None):
  
  if ymdhm is None:
    ymdhm  = time.strftime("%y%m%d.%H%M")
  inpDir = os.path.join(valDir,'Inputs','AERMOD')
  outDir = os.path.join(valDir,'Outputs',ymdhm,'AERMOD')
  pltDir = os.path.join(valDir,'Outputs',ymdhm,'Plots')
  
  return(inpDir,outDir,pltDir)

def getObsArray(prjName,obsDir=None,obsName=None,iHr=1):
  
  if obsName is None:
    obsNameList = {'kinso2':'kinso2','kinsf6':'KINSF6.SUM','bowline':'bow','clifty':'CCR75'}
    obsNameList.update({'baldwin':'bal','martin':'MCR','tracy':'TRACY','pgrass':'PGRSF6.SUM'})
    try:
      obsName = obsNameList[prjName.lower()]
    except KeyError:
      if prjName.upper().startswith('KSF6'):
        obsName = 'KINSF6.SUM'      
  if obsName is None:
    obsArray = []
    return obsArray
  if obsDir is None:
    obsDir = os.path.join('..','Observations')
  # Read observations from obs directory from .obs files
  if prjName == 'pgrass':
    obsFile = os.path.join(obsDir,obsName)
    colList = (i for i in range(5,7))
    QCArray    = np.loadtxt(obsFile,skiprows=4,usecols=colList)
    obsArray   = QCArray[:,0]*QCArray[:,1]*1.e6 # ug/m3
    print np.shape(obsArray)
  elif prjName.upper().startswith('KSF6'):
    obsFile = os.path.join(obsDir,obsName+'.db')
    YY = prjName[9:11]
    MM = prjName[5:6]
    DD = prjName[6:8]
    qryStr  = 'select max(CHI) from dataTable where YY = %s and MM = %s and DD = %s group by HH'%(YY,MM,DD)
    obsArray = utilDb.db2Array(obsFile,qryStr,dim=1)
    obsArray = obsArray/167.  # Convert from ppt to ug/m3 for SF6
  elif prjName.lower().startswith('martin'):
    obsFile    = os.path.join(obsDir,obsName + '%02d.OBS'%iHr)
    useCols = [i for i in range(11)]  # Extra sampler compared to AERMOD RE list !!
    obsHrArray = np.loadtxt(obsFile,usecols=useCols)
    nHr        = len(obsHrArray)
    obsArray   = np.zeros(nHr) - 999.
    for hr in range(nHr):
      obsArray[hr] = max(obsHrArray[hr,3:])
  else:
    obsFile    = os.path.join(obsDir,(obsName + '%02d.obs'%iHr).upper())
    obsHrArray = np.loadtxt(obsFile)
    nHr        = len(obsHrArray)
    obsArray   = np.zeros(nHr) - 999.
    for hr in range(nHr):
      obsArray[hr] = max(obsHrArray[hr,4:])
  obsArray  = np.sort(obsArray)[::-1]
  print 'OBS max =', obsArray[0],obsArray[min(25,len(obsArray))-1],min(25,len(obsArray))
  #for iSmp in range(nSmp):
  #  obsArray[:,iSmp] = np.sort(obsArray[:,iSmp])[::-1]
  #  print 'OBS max =', obsArray[0,iSmp],obsArray[25,iSmp]
  return obsArray

def getAerArray(prjName,aerDir=None,aerName=None,iHr=1):
  
  if aerName is None:
    aerNameList = {'kinso2':'KS2AER','kinsf6':None,'bowline':'BOWAER','clifty':'CCRAER'}
    aerNameList.update({'baldwin':'BALAER','martin':'MCRAT2','tracy':'TRAER','pgrass':'PGRASS'})
    try:
      aerName = aerNameList[prjName.lower()]
    except KeyError:
      if prjName.upper().startswith('KSF6'):
        aerName = None
    
  if aerDir is None:
    aerDir = os.path.join('..','Aermod')
 
  if prjName.upper().startswith('KSF6'):
    aerName = os.path.join(aerDir,prjName.upper()).replace('I','R').replace('_','.')
    useCols = [i for i in range(8)]
    aerConc = np.loadtxt(aerName,skiprows=8,usecols=useCols)
    aerArray = aerConc[:,1]
    aerArray = np.sort(aerArray)[::-1]
  elif prjName.upper().startswith('BALDWIN'):
    aerArray = np.loadtxt("/home/user/bnc/AERMOD/v15181/runs/baldwin/AERMOD/balaer%02d.rnk"%iHr,skiprows=8,usecols=[1])
  else:
    aerFile = os.path.join(aerDir,aerName + '%02d.PST.db'%iHr)
    aerConn = sqlite3.connect(aerFile)
    aerConn.row_factory = sqlite3.Row
    aerCur = aerConn.cursor()    
    #aerQry = 'select Cavg from datatable order by Cavg desc'
    # Remove duplicate data period similar to RANK-FILE
    aerQry = 'select max(Cavg) as cMax from datatable group by date order by cMax desc'
    aerArray = utilDb.db2Array(aerCur,aerQry,dim=1)
    print 'AERMOD max =', aerArray[0],aerArray[min(25,len(aerArray))-1],min(25,len(aerArray))
    aerConn.close()
  
  return aerArray

def getSmpDb(prjName):
  mySciFiles = SCI.Files(prjName)
  smpDb = '%s.smp.db'%(prjName)
  (smpDbConn,smpDbCur,smpCreateDb) = utilDb.Smp2Db(smpDb,mySciFiles)
  return (smpDbConn,smpDbCur)

def countSmpDb(smpDbCur):
  smpIds = map(int,utilDb.db2Array(smpDbCur,'select distinct(smpid) from samTable'))
  nSmp   = len(smpIds)
  nTimes = utilDb.db2Array(smpDbCur,"select count(value) from samTable a, smptable p where a.colno=p.colno \
                          and varname='C' and smpId =1",dim=0)
  print 'nSmp,nTimes = ',nSmp,nTimes
  return (nTimes,nSmp,smpIds)

def getPstMax(iHr,prjName,sciFac):
  
  #pstFile  = prjName + '.' + '%dhr'%iHr + '.pst'
  pstFile  = 'so2' + '.' + '%dhr'%iHr + '.pst'
  
  skipRows = 8
  colNames = ['x', 'y', 'avgC', 'zlev', 'zhill', 'zflag', 'ave' ,'grp', 'date']
  useCols  = [0, 1, 2]

  print 'pstFile = ',pstFile

  df = pd.read_table(pstFile,skiprows=skipRows,sep=r'\s*',names=colNames)
  #print df['avgC'].max(),df['avgC'][0],df['avgC'][-1]
  print df.describe()
  
  sciArray = df.sort('avgC',axis=0,ascending=False)['avgC'].values
  
  sciArray = sciArray*sciFac
  print 'SCICHEM max = 1:',sciArray[0],', ',sciArray[min(25,len(sciArray))-1],':',min(25,len(sciArray))
  
  return sciArray
  
  

def smpDbMax(smpDbCur,iHr,smpFac,smpIds=None,nTimes=None):
  
  if nTimes is None:
    (nTimes,nSmp,smpIds) = countSmpDb(smpDbCur)
  else:
    nSmp = len(smpIds)
    
  if iHr > 1:
    preQry  = "select value from samTable a, smptable p where a.colno=p.colno "
    preQry += "and varname='C' and smpId ="
    sciMax  = np.zeros((nSmp,nTimes/iHr),float)
    for j,smpId in enumerate(smpIds):
      sciQry = preQry + str(smpId) + ' order by time'
      hrMax  =  utilDb.db2Array(smpDbCur,sciQry)
      # avgArray[j,:] = np.convolve(hrMax, np.ones(iHr)/iHr)
      for i in range(0,len(hrMax)-iHr+1,iHr):
        sciMax[j,i/iHr] = np.mean(hrMax[i:i+iHr])    
    sciMax   = np.reshape(sciMax,np.size(sciMax))  
    sciArray = np.sort(sciMax)[::-1]
  else:
    sciQry  = "select max(value) as maxVal from samTable a, smptable p where a.colno=p.colno "
    sciQry += "and varname='C' and (time -round(time)) = 0. group by time order by maxVal desc"
    sciArray = utilDb.db2Array(smpDbCur,sciQry,dim=1)
  sciArray = sciArray*smpFac
  print 'SCICHEM max = 1:',sciArray[0],', ',sciArray[min(25,len(sciArray))-1],':',min(25,len(sciArray))
  
  return sciArray

def plotData(obsArray, aerArray, sciArray, figName, figTitle, cutoff=0.0, units=None):
  minLen = min(len(obsArray),len(aerArray),len(sciArray))
  aerArray = aerArray[:minLen][obsArray[:minLen] > cutoff]
  sciArray = sciArray[:minLen][obsArray[:minLen] > cutoff]
  obsArray = obsArray[:minLen][obsArray[:minLen] > cutoff] 
  maxVal = max(max(obsArray),max(aerArray),max(sciArray))
  plt.clf()
  plt.hold(True)
  ax = plt.subplot(111)
  LSCI = plt.scatter(obsArray,sciArray,marker='o',color='r')
  LAER = plt.scatter(obsArray,aerArray,marker='d',color='b')
  plt.xlim(0,maxVal)
  plt.ylim(0,maxVal)
  plt.plot([0,maxVal],[0,maxVal],'k-')
  plt.plot([0,maxVal],[0,maxVal*0.5],'r-')
  plt.plot([0,maxVal],[0,maxVal*2],'r-')
  if units is None:
    units = r'($\mu g/m^3$)'
  plt.xlabel(r'Observed('  + units + ')')
  plt.ylabel(r'Predicted(' + units + ')')
  #ax.set_xscale('log')
  #ax.set_yscale('log')
  plt.hold(False)
  plt.title(figTitle)
  plt.legend([LSCI,LAER],['SCICHEM', 'AERMOD'],loc='upper left')
  plt.savefig(figName)   
  return

# Code for SCICHEM plots
def cmpAerMod(inpDir,prjName,dstDir=None,color=True):
 
  sciFac = 1.e+9 # kg/m3 to microg/m3
  MaxVals = ['1 hr', '3hr', '24hr', 'All']
  MaxObs  = [0,0,0,0]  # OBS
  MaxPre1 = [0,0,0,0]  # SCICHEM-2012
  MaxPre2 = [0,0,0,0]  # AERMOD
  
  prjName = prjName.lower()
  if prjName == 'kinso2':
    yMax = 2500
  if prjName == 'kinsf6':
    yMax = 2500
  if prjName == 'bowline':
    yMax = 850
  if prjName == 'baldwin':
    yMax = 3100
  if prjName == 'clifty':
    yMax = 2500
  if prjName == 'martin':
    yMax = 2500
  if prjName == 'tracy':
    yMax = 2500
  if prjName == 'pgrass':
    yMax = 2500
  
  statFile = open(prjName+"_stat.csv","w")
  statFile.write("Case, Dur, maxObs, maxAER, maxSCI, RHC_AER, RHC_SCI\n")

  # SCICHEM predictions
  sciConn,sciCur       = getSmpDb(prjName)
  (nTimes,nSmp,smpIds) = countSmpDb(sciCur)
  
  sciBcData = []
  aerBcData = []
  obsBcData = []
  
  for iHr in [1,3,24]:
    
    if iHr > 1 and (prjName == 'tracy' or prjName == 'pgrass'):
      continue
    
    # Read observations from obs directory from .obs files
    obsDir = os.path.join(inpDir,'Observations')
    obsArray = getObsArray(prjName,obsDir=obsDir,iHr=iHr)
    
    # Get AERMOD max concentration array
    aerDir = os.path.join(inpDir,'AERMOD')
    aerArray = getAerArray(prjName,aerDir=aerDir,iHr=iHr)

    # Get SCIPUFF max concentration array
    #sciArray = getPstMax(iHr,prjName,1.0)

    sciArray = smpDbMax(sciCur,iHr,sciFac,smpIds=smpIds,nTimes=nTimes)
    # Get maximum 1 hr avg SCIPUFF concentration
    #sciArray = getMaxConc(prjName)*cFac
   
    #
    # Print top 10 values
    #
    for val in zip(obsArray[:10], aerArray[:10],sciArray[:10]):
      print '%5.1f %5.1f %5.1f'%(val[0],val[1],val[2])
    
    sciBcData.append(sciArray[0])
    aerBcData.append(aerArray[0])
    obsBcData.append(obsArray[0])
    if iHr == 1:
      sciAll = np.mean(sciArray)
      aerAll = np.mean(aerArray)
      obsAll = np.mean(obsArray)
    #continue
  
  
    #    
    #
    # Q-Q plot
    # 
    figName = prjName + str(iHr) +'.png'
    figTitle = '%s: Max Concentration for %02d Hr Average Concentrations'%(prjName.upper(),iHr)
    plotData(obsArray, aerArray, sciArray, figName, figTitle)
    
    obsRHC = measure.rhc(obsArray)
    sciRHC = measure.rhc(sciArray)
    aerRHC = measure.rhc(aerArray)
    print aerRHC/obsRHC,sciRHC/obsRHC
    #calcStats(obsMax1Hr,pre1Max1Hr,statFile)
  
    if statFile is not None:
      # Case, Dur, maxObs, maxAER, maxSCI, RHC_AER, RHC_SCI
      statFile.write("%s, %02d hr, %10.3f, %10.3f, %10.3f, %6.3f, %6.3f\n"%\
                    (prjName,iHr,obsArray[0],aerArray[0],sciArray[0],aerRHC/obsRHC,sciRHC/obsRHC))
    if dstDir is not None:
        newName = os.path.join(dstDir,figName)
        try:
          shutil.move(figName,newName)
          print 'Moved %s to plots dir'%newName
        except EnvironmentError:
          print 'Failed to move %s to plots dir'%newName
  #
  # Bar plot
  #
  obsBcData.append(obsAll)
  sciBcData.append(sciAll)
  aerBcData.append(aerAll)
  #pltName = prjName +'_bar.eps'
  pltName = prjName +'_bar.png'

  if prjName == 'kinso2':
    title = 'Figure 3: Comparison of AERMOD and SCICHEM predicted maximum concentrations'
    title2 = 'averaged over different time periods for Kincaid SO2 study'
  elif prjName == 'baldwin':
    title = 'Figure 4: Comparison of AERMOD and SCICHEM predicted maximum concentrations'
    title2 = 'averaged over different time periods for Baldwin study'
  elif prjName == 'bowline':
    title = 'Figure 5: Comparison of AERMOD and SCICHEM predicted maximum concentrations'
    title2 = 'averaged over different time periods for Bowline study'
  elif prjName == 'clifty':
    title = 'Figure 6: Comparison of AERMOD and SCICHEM predicted maximum concentrations'
    title2 = 'averaged over different time periods for Clifty study'    
  else:
    title = 'Figure 7: Comparison of AERMOD and SCICHEM predicted maximum concentrations'
    title2 = 'averaged over different time periods for % study'%prjName 
  pltBarAermod.plotBarGraph(obsBcData,sciBcData, aerBcData,title, pltName,yMax=yMax, xTicLab=['1Hr','3Hr','24hr','All'],title2=title2)
  
# Main program
if __name__ == '__main__':

  compName = socket.gethostname()
  env = os.environ.copy()
  
  # Local modules
  if  compName == 'pj-linux4':
    sys.path.append('/home/user/bnc/python')
    runDir = ''
  if compName == 'sm-bnc':
    outDir          = 'D:\\TestSCICHEM\\Outputs\\2014_03_20'
    SCIPUFF_BASEDIR = "D:\\SCIPUFF\\Repository\\workspace\\EPRI\\bin\\intel\\Win32\\Release"
    INIfile         = "D:\\SCIPUFF\\Repository\\Workspace\\EPRI\\scipuff.ini"
    
  prjName = 'Clifty'
  runDir  = os.path.join(outDir,'AERMOD',prjName,'SCICHEM')
  inpDir  = os.path.join(os.path.dirname(outDir).replace('Outputs','Inputs'),'AERMOD',prjName)
  os.chdir(runDir)
  
  if 'kinsf6' in runDir.upper():
    units   = 'PPT'
    matName = 'SF6'
  else:
    units = 'ug/m3'
  
  if units == 'PPT' and matName == 'SF6':
    cFac = 167. # Convert from ug/m3 to ppt for SF6
  else:
    cFac = 1.
    
  cmpAerMod(inpDir,prjName)       
