
#!/bin/env python
import os
import sys
import time
import socket
import shutil
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import sqlite3

compName = socket.gethostname()

# Local modules
if compName == 'sm-bnc' or compName == 'sage-d600':
  sys.path.append('C:\\Users\\sid\\python')
if  compName == 'pj-linux4':
  sys.path.append('/home/user/bnc/python')
  
#import measure
import utilDb
import setSCIparams as SCI

# Code for SCICHEM 2012 plots

def cmpTva(inpDir,prjName,prePfx2=None):
  
  obsPfx,preCur1,preCur2 = getDBConn(inpDir,prjName)
  
  varNames = ["SO2", "O3", "NOx",  "NOy"]   #["SO2", "O3", "NOx",  "NOy"]
  
  preCreateStr = utilDb.db2List(preCur1,'select sql FROM sqlite_master where name = "samTable"')
  
  ug2ppm = False
  if 'mcUnits' in preCreateStr[0][0]:
    mcUnits = utilDb.db2List(preCur1,'select mcUnits from samTable limit 1')[0][0]
    if mcUnits == 'ug/m3':
      ug2ppm  = True
      fug2ppm = [0.075/196.,0.075/147.,0.053/100.,0.053/100.]

  if "tva_980825" in prjName:
    distance = [20, 55, 110]
    times    = [11.5, 12.75, 14.75]
    zSmp     = [520, 600, 620]
    zSmp2    = [520, 600, 620]
    
  if "tva_980826" in prjName:
    distance = [18, 27, 86, 59, 93]#, 141
    times    = [10.25, 10.75, 12.25, 12.5, 12.75]#, 13.25
    zSmp     = [465, 500, 659, 910, 819, 662]   
    zSmp2    = [465, 500, 659, 910, 819, 662]   

  if "tva_990706" in prjName:
    distance = [11, 31, 65]  #, 89]
    times    = [12.5, 13.25, 16.0] #[12.0, 13.0, 16.0, 16.75]
    #times    = [10.5, 10.0, 10.0] Test for o3 wings
    zSmp     = [501, 505, 500] #[505, 491, 448, 533]
    zSmp2    = [501, 505, 500] #, 533]
  
  if "tva_990715" in prjName:
    distance = [16, 62, 106]
    times    = [11.75, 12.5, 17]#[11.5, 12.25, 16.5] #
    zSmp     = [411, 727, 438] #[445, 589, 587]
    zSmp2    = [415, 584, 582]
  
  statFile = open(prjName+"_stat.csv","w")
  statFile.write("Case, Distance, PlumeNo, Species, Mean Obs, Mean Pre, NMSE, MFB, r\n")
  
  for idt,dist in enumerate(distance): 

    if preCur2 is None:
      if prePfx2 is None:
        print "Skipping plots for SCICHEM-99 prediction"
      else:
        # Predictions SCICHEM-01
        pre2DbName = prePfx2+'.smp.db' #prePfx2 + '_' + str(dist) + 'km' + '.csv.db'       
        print 'Predicted DB SCICHEM99', pre2DbName
        preConn2 = sqlite3.connect(pre2DbName)
        preConn2.row_factory = sqlite3.Row
        preCur2 = preConn2.cursor()
        print preConn2
    
    lArc = 2.*dist*np.pi/180.
   
    for varName in varNames:

      print '\n',dist,' km, ',varName

      # Prediction query
      if varName == "NOx":
        preQry1 = 'select xSmp,Sum(Value) from samTable a,smpTable p where a.colNo=p.colNo'
        preQry1 += " and varName in ('NO2', 'NO' ) "
        preQry1 += " and zSmp = %f and time = %3f group by smpId"%(zSmp[idt],times[idt])
      elif varName == 'NOy':
        preQry1 = "select xSmp,Sum(Value)from samTable a,smpTable p where a.colNo=p.colNo"
        preQry1 += " and varName in ('NO2','NO','NO3','N2O5','HNO3','HONO','PAN'  ) "
        preQry1 += " and zSmp = %f and time = %3f group by smpId"%(zSmp[idt],times[idt])
      else:  
        preQry1  = "select xSmp,Value from samTable a,smpTable p where a.colNo=p.colNo and "
        preQry1 += "varName = '%s' and zSmp = %f and time = %3f order by smpId"%(varName,zSmp[idt],times[idt])      
      print preQry1
      
      # Predictions #1 (2012)
      preArray1 = utilDb.db2Array(preCur1,preQry1)     
      if ug2ppm:
        preArray1[:,1] = preArray1[:,1]*fug2ppm[varNames.index(varName)]
      
      if varName == 'SO2':
        # Set x index where SO2 is max. Same for all obs plumes
        iMax1 = np.where(preArray1[:,1] == preArray1[:,1].max())[0][0]
          
      # Predictions #2 (v2100)
      if varName == "NOx":      
        preQry2 = 'select xSmp,Sum(Value) from samTable a,smpTable p where a.colNo=p.colNo'
        preQry2 += " and varName in ('NO2', 'NO' ) "
        preQry2 += " and zSmp = %f and time = %3f group by smpId"%(zSmp2[idt],times[idt])
      elif varName == 'NOy':
        preQry2 = "select xSmp,Sum(Value)  from samTable a,smpTable p where a.colNo=p.colNo"
        preQry2 += " and varName in ('NO2','NO','NO3','N2O5','HNO3','HONO','PAN'  ) "
        preQry2 += " and zSmp = %f and time = %3f group by smpId"%(zSmp2[idt],times[idt])
      else:  
        preQry2  = "select xSmp,Value from samTable a,smpTable p where a.colNo=p.colNo and "
        preQry2 += "varName = '%s' and zSmp = %f and time = %3f order by smpId"%(varName,zSmp2[idt],times[idt])      
      
      preArray2 = utilDb.db2Array(preCur2,preQry2)
      #print preQry2  
      if varName == 'SO2':
        # Set x index where SO2 is max. Same for all obs plumes
        iMax2 = np.where(preArray2[:,1] == preArray2[:,1].max())[0][0]
          
      # Observed data
      if "tva_980825" in prjName:
        if dist == 20:
          pls = [3,4,5]
        if dist == 55:
          pls = [6,7,8]
        if dist == 110:
          pls = [9,10,11]
          
      if "tva_980826" in prjName:
        if dist == 18:
          pls = [1,2]
        if dist == 27:
          pls = [4,5,6]
        if dist == 86:
          pls = [8]
        if dist == 59:
          pls = [9] 
        if dist ==93:
          pls = [10]
        if dist == 141:
          pls = [11]    

      if "tva_990706" in prjName:
        if dist == 11:
          pls = [1,2,3,4,5]
        if dist == 31:
          pls = [6,7,8]
        if dist == 65:
          pls = [9,10,11]
        if dist == 89:
          pls = [12]
          
      if "tva_990715" in prjName:
        if dist == 16:
          pls = [1,2,3,4]
        if dist == 62:
          pls = [6,7,8]
        if dist == 106:
          pls = [9,10,11]

      # Set Omax to zero for max plume no.
      if varName == 'SO2':
        oMax  = [0 for i in range(12)]
      
      for ipl in pls:
        
        obsDbName = obsPfx + str(dist) + 'km_obs' + str(ipl) + '.csv.db'
        print '\n obsDbName = ',obsDbName
      
        obsConn = sqlite3.connect(obsDbName)
        obsConn.row_factory = sqlite3.Row
        obsCur = obsConn.cursor()

        # Observations
        if "tva_990706" in prjName:
          obsQry = 'select CAST(plumeKM as real) as x, avg(' + varName + '/1000.) from dataTable where plumeKM > -9999.0 group by x'
        else:
          obsQry = 'select CAST(plumeKM as real)as x, avg(' + varName + ') from dataTable  where plumeKM > -9999.0 group by x' 

        print obsQry
        
        obsArray = utilDb.db2Array(obsCur,obsQry)
        if varName == 'SO2':
          oMax[ipl] = np.where(obsArray[:,1] == obsArray[:,1].max())[0][0]
          print 'Max Obs SO2 at x index = ', oMax[ipl],' with x, value = ',obsArray[oMax[ipl],0],obsArray[oMax[ipl],1]

        # Set x value to 0. where SO2 conc is maximum for plume ipl
        obsArray[:,0] = obsArray[:,0] - obsArray[oMax[ipl],0]
        print ipl,obsArray[:,0]
        statFile.write("%s, %02d, %d, %s,"%(prjName,dist,ipl,varName))

        print '\nFor varName %s, plume no. %d,x index = %d, x = %f ,conc = %f\n'%\
              (varName,ipl,oMax[ipl],obsArray[oMax[ipl],0],obsArray[oMax[ipl],1])
        
        # Align the prediction1 and observed centerlines
        for i in range(len(preArray1)):
          # Set x values at iMax1(index where SO2 is max) to 0 for plume ipl. 
          preArray1[i,0] = float(i-iMax1)*lArc 
          if i == iMax1:
            print 'Adjusted preArray1 x,c = ',preArray1[i,:] 
          
        # Align the prediction2 and observed centerlines
        for i in range(len(preArray2)):
          # Set x values where SO2 max to obs x value for plume ip1
          preArray2[i,0] = float(i-iMax2)*lArc 
          if i == iMax2:
            print 'Adjusted preArray2 x,c = ',preArray2[i,:]
           
        # Save npz for later plotting
        npzFile = str(dist) +'km_'+ str(ipl) + '_' + varName + '_'+ str(times[idt]) + 'hr' + '.npz'
        figTitle = '%s (@ %d km)'%(varName,dist)
        print npzFile
        np.savez(npzFile,obsArray,preArray1,preArray2)
           
        # Get statistics
        pre1 = np.copy(obsArray)
        pre1[:,1] = np.interp(obsArray[:,0],preArray1[:,0],preArray1[:,1])
        print obsArray[:,1].max(),pre1[:,1].max()
        calcStats(obsArray,pre1,statFile=statFile)
        if varName == varNames[-1] and ipl == pls[-1]:
          statFile.write("\n")

        obsConn.close()

    if prePfx2 is not None:
      print 'preConn2 is still open'
      #preConn2.close()
          
  return

def getDBConn(inpDir,prjName):
  
  # Predicted SCICHEM-3.0 data
  prjName1 = os.path.join(prjName)
  print '**********\n',prjName1
  preConn1,preCur1 = getSmpDb(prjName1)

  # Observed data 
  obsPfx = os.path.join(inpDir,prjName,'OBSERVATIONS',prjName + '_')
  print obsPfx
  
  # Use prePfx2 + '_' + str(dist) + 'km' + '.csv.db'
  prjName2 = os.path.join(inpDir,prjName,'SCICHEM-99',prjName)
  print prjName2
  preConn2,preCur2 = getSmpDb(prjName2)
  
  return obsPfx,preCur1,preCur2

def createSubPlots(prjName,dstDir=None,skip99=True,color=True):
  print 'Creating sub plots'
  params1 = {'axes.labelsize': 10, 'text.fontsize': 10, 'xtick.labelsize': 10,
                'ytick.labelsize': 10, 'legend.pad': 0.1,  
                'legend.fontsize': 8, 'lines.markersize': 6, 'lines.width': 2.0,
                'font.size': 10, 'text.usetex': False}
  plt.rcParams.update(params1)
  varName = ['NOx','NOy','O3','SO2']
  isPPB = True  

  fList = os.listdir('.')
  for fName in fList:
    if fName.endswith('.npz'):
      if 'SO2' in fName:
        fig = plt.figure()
        plt.clf()
        fig.hold(True)        
        dist = fName.split('_')[0].replace('km','')
        plt.title('Downwind distance %s km for %s'%(dist,prjName))
        ax = fig.add_subplot(111)
        plt.text(.3,-0.09,'Cross plume distance(km)',transform=ax.transAxes,fontsize=14)
        if isPPB:
          cUnits = 'ppb'
          plt.text(-.1,0.8,'Perturbation Concentration (%s)'%cUnits,transform=ax.transAxes,rotation='vertical',fontsize=14)
        else:
          cUnits = 'ppm' 
          plt.text(-.15,0.8,'Perturbation Concentration (%s)'%cUnits,transform=ax.transAxes,rotation='vertical',fontsize=14)
        plt.setp(plt.gca(), frame_on=False, xticks=(), yticks=())
        figName = prjName + '_'  + fName.replace('_SO2_','_').replace('.npz','.png')
        for i in range(4):
          fig,LhO,LhP1,LhP2 = subplotConc(fig,fName,dist,i,varName[i],skip99=skip99,color=color,isPPB=isPPB)
        if skip99:
          fig.legend([LhO,LhP1,LhP2],['OBSERVED','SCICHEM'])
        else:
          fig.legend([LhO,LhP1,LhP2],['OBSERVED','SCICHEM','SCICHEM-99'])
        fig.hold(False)
        plt.savefig(figName)
        if dstDir is not None:
          newName = os.path.join(dstDir,figName)
          try:
            shutil.move(figName,newName)
            print 'Moved %s to plots dir'%newName
          except EnvironmentError:
            print 'Failed to move %s to plots dir'%newName
            
def subplotConc(fig,fName,dist,i,varName,skip99=False,color=False,isPPB=False):
  
  rmAmb = False
  dist  = float(fName.split('_')[0].replace('km',''))  
  fName = fName.replace('SO2',varName)  
  print fName,', Remove ambient = ',rmAmb
  
  
  npzfile = np.load(fName)
  obsData = npzfile['arr_0'] 
  preData1 = npzfile['arr_1'] 
  preData2 = npzfile['arr_2']
  
  ax = fig.add_subplot(2,2,i+1)
  #print obsData[:,1]
  if varName == 'O3':
    if dist == 110 or dist == 55:
      cO = min(obsData[-5:-1,1].min(),obsData[0:5,1].min())
    else:
      cO = max(obsData[-1,1],obsData[0,1])
    if rmAmb:
      c1 = max(preData1[-1,1],preData1[0,1])
    else:
      c1 = 0.
    c2 = max(preData2[-1,1],preData2[0,1])
    plt.text(0.85,0.1,'%s'%varName,transform=ax.transAxes)
  else:
    cO = min(obsData[-1,1],obsData[0,1])
    if rmAmb:
      c1 = min(preData1[-1,1],preData1[0,1])
    else:
      c1 = 0.
    c2 = min(preData2[-1,1],preData2[0,1])
    #plt.text(0.1,0.8,'%s'%varName,transform=ax.transAxes)
    plt.text(0.85,0.9,'%s'%varName,transform=ax.transAxes)
        
  print 'Plume Min/Max =',varName,cO,c1,c2
  
  if color:
    obsColor = 'red'
    sci3Color = 'green'
    sci1Color = 'blue'
  else:
    obsColor = '0.25'
    sci3Color = '0.5'
    sci1Color = '0.75'    
  # Obs
  C = obsData[:,1] - cO
  if isPPB:
    C = C*1.e3
  LhO, = plt.plot(obsData[:,0],C,linestyle='None',marker='o',markersize=6,markerfacecolor=obsColor) 

  # SCICHEM 3.0
  #C = ma.masked_where(preData1[:,1]<0.,preData1[:,1]) - c1
  C = preData1[:,1] - c1
  if isPPB:
    C = C*1.e3
  LhP1, = plt.plot(preData1[:,0],C,linestyle='-',color=sci3Color,marker='s',markersize=6,markerfacecolor=sci3Color)

  # SCICHEM 1.0 
  if skip99:
    LhP2 = None
  else:
    # Skip plots for SCICHEM-99
    C = ma.masked_where(preData2[:,1]<0.,preData2[:,1]) - c2
    LhP2, = plt.plot(preData2[:,0],C,linestyle='-',color=sci1Color,marker='^',markersize=6,markerfacecolor=sci1Color)

  if int(dist) < 30:
    plt.xlim([-10,10])
  else:
    plt.xlim([-30,30])
  
  #if varName == 'O3':
  #  plt.ylim([0,0.1])
  if i == 0 or i == 1:
    ax.set_xticklabels([])
    
  return fig,LhO,LhP1,LhP2
  
          
def pltCmpConc(dist, varName, obsData, preData1, preData2, figTitle, figName):
  #import pdb; pdb.set_trace()
    # Set up plot parameters
  params1 = {'axes.labelsize': 10, 'text.fontsize': 10, 'xtick.labelsize': 10,
                'ytick.labelsize': 10, 'legend.pad': 0.1,  
                'legend.fontsize': 8, 'lines.markersize': 6, 'lines.width': 2.0,
                'font.size': 10, 'text.usetex': False}
  plt.rcParams.update(params1)
  
  fig = plt.figure()
  plt.clf()
  fig.hold(True)
  ax = fig.add_subplot(111)
  #print obsData[:,1]
  if varName == 'O3':
    fac = 1
    cO = max(obsData[-1,1],obsData[0,1])
    c1 = max(preData1[-1,1],preData1[0,1])
    c2 = max(preData2[-1,1],preData2[0,1])
  else:
    fac = 1
    cO = min(obsData[-1,1],obsData[0,1])
    c1 = min(preData1[-1,1],preData1[0,1])
    c2 = min(preData2[-1,1],preData2[0,1])
  C = obsData[:,1] - fac*cO
  LhO  = plt.plot(obsData[:,0],C,linestyle='None',marker='o',markersize=6,markerfacecolor='0.25') 
  LkO  = 'OBS'
  C = ma.masked_where(preData1[:,1]<0.,preData1[:,1]) - fac*c1
  LhP1 = plt.plot(preData1[:,0],C,linestyle='-',color='0.5',marker='s',markersize=6,markerfacecolor='0.5')
  LkP1 = 'SCICHEM-3.0'
  C = ma.masked_where(preData2[:,1]<0.,preData2[:,1]) - fac*c2
  LhP2 = plt.plot(preData2[:,0],C,linestyle='-',color='0.75',marker='^',markersize=6,markerfacecolor='0.75')
  LkP2 = 'SCICHEM-1.0'
  plt.ylabel('Perturbation Concentration (ppm)')
  plt.xlabel('Cross plume distance (km)')
  plt.title(figTitle)
  if dist < 21:
    plt.xlim([-20,20])
  else:
    plt.xlim([-50,50])
  #if varName == 'O3':
  #  plt.ylim([0,0.1])
  plt.legend([LkO,LkP1,LkP2],bbox_to_anchor=(0.02,0.98),loc=2,borderaxespad=0.)
  lgnd  = plt.gca().get_legend()
  ltext = lgnd.get_texts()
  plt.setp(ltext,fontsize=9)
  fig.hold(False)
  plt.savefig(figName)
  #sys.exit()
  #plt.show()
  return

def calcStats(obsArray, preArray,statFile=None):
  obsMean =  np.mean(obsArray[:,1])
  predMean = np.mean(preArray[:,1])
  print 'Min:',min(obsArray[:,1]),min(preArray[:,1])
  if min(obsArray[:,1])< 0.:
    pass
  #NMSE = measure.nmse_1(preArray[:,1],obsArray[:,1], cutoff=0.0)
  #biasFac = measure.bf(preArray[:,1],obsArray[:,1], cutoff=0.0)
  #MFB = measure.mfbe(preArray[:,1],obsArray[:,1], cutoff=0.0)
  #mnbe = measure.mnbe(preArray[:,1],obsArray[:,1], cutoff=0.0)
  #correlation = measure.correlation(preArray[:,1],obsArray[:,1], cutoff=0.0 )
  #mage = measure.mage(preArray[:,1],obsArray[:,1])
  #print 'stats: ',obsMean, predMean, NMSE, biasFac, correlation
  #if statFile is not None:
  #  statFile.write("%8.3f, %8.3f, %8.3f, %8.3f,%8.3f\n"%(obsMean,predMean,NMSE,MFB,correlation))

def getSmpDb(prjName):
    mySciFiles = SCI.Files(prjName)
    smpDb = '%s.smp.db'%(prjName)
    print '%%%%%%%%%%%%', smpDb
    # Create database for calculated data 
    ## print 'Create smpDb ',smpDb,' in ',os.getcwd()
    (smpDbConn,smpDbCur,smpCreateDb) = utilDb.Smp2Db(smpDb,mySciFiles)
    return (smpDbConn,smpDbCur)

def setChemDir(valDir,ymdhm=None):
  
  if ymdhm is None:
    ymdhm  = time.strftime("%y%m%d.%H%M")
  inpDir = os.path.join(valDir,'Inputs','Chemistry')
  outDir = os.path.join(valDir,'Outputs',ymdhm,'Chemistry')
  pltDir = os.path.join(valDir,'Outputs',ymdhm,'Plots')
  
  return(inpDir,outDir,pltDir)
            
# Main program
if __name__ == '__main__':

  if compName == 'sm-bnc' or compName == 'sage-d600':
    #valDir = 'D:\\TestSCICHEM'
    valDir = 'D:\\TestSCICHEM'
  if compName == 'pj-linux4':
    valDir = '/home/user/bnc/TestSCICHEM'
  
  #typeList = 'v3b2  v3b2_TN3 v3b2_TN10 v3b2_HiRes v3b2_FixAmb v3b2_FixAmb_HiISOP v3b2_FixAmb_HiOHHO2 v3b2_FixAmb_LoISOP v3b2_FixAmb_LoOHHO2'
  typeList = 'v3b2_TN3_noLSV_utc'
  
  for ymdhm in typeList.split():
     
    inpDir,outDir,pltDir = setChemDir(valDir,ymdhm=ymdhm)
    #130511_EPRI_git for 98 cases
    #loLSV_HiVOC for 99 cases
    
    for prjName in  ['tva_980825','tva_990706']: # ['tva_980825','tva_980826', 'tva_990706', 'tva_990715']
      runDir = os.path.join(outDir,prjName,'SCICHEM')
      os.chdir(runDir)
      print 'runDir = ',runDir
      cmpTva(inpDir,prjName)
      #createSubPlots(prjName,dstDir=pltDir)
      #createSubPlots(prjName)
  
  print 'Done'
