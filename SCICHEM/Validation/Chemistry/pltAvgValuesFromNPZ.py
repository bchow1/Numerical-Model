#!/bin/env python
import os
import sys
import socket
import numpy as np
import numpy.ma as ma
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import sqlite3

compName = socket.gethostname()

# Local modules
if  compName == 'linux':
  sys.path.append('/home/user/bnc/python')

import measure
import utilDb
import setSCIparams as SCI

# Code for SCICHEM 2012 plots
def pltAvgData(prjName,preDirs=None,preLbls=None):
  rmAmb = True
  print 'Plot averaged data with rmAmb set to ',rmAmb
  params1 = {'axes.labelsize': 8, 'text.fontsize': 8,'text.family':'Helvetica' ,'xtick.labelsize': 8,
            'ytick.labelsize': 8, 'legend.pad': 0.1,  
            'legend.fontsize': 8, 'lines.markersize': 5, 'lines.width': 2.0,
            'font.size': 8, 'text.usetex': False}
  plt.rcParams.update(params1)
  varName = ['NOx','NOy','O3','SO2']  
  lblForManus = ['(a)','(b)','(c)','(d)','(e)','(f)' ] # Using lblForManus[lblCt]
  lblCt = 0

  if "tva_980825" in prjName:
    distance = [20, 55, 110]
    times    = [11.5, 12.75, 14.75]
    zSmp     = [520, 600, 620]
    zSmp2    = [520, 600, 620]
    xLim     = [[-15,15],[-25,25],[-30,30]]
    #yLim     = {'NOx':[[0,200],[0,80],[0,25]],'NOy':[[0,210],[0,90],[0,30]],\
    #            'O3':[[-50,10],[-30,20],[-10,40]],'SO2':[[0,20],[0,10],[0,5]]}

    yLim     = {'NOx':[[0,200],[0,80],[0,25]],'NOy':[[0,210],[0,90],[0,30]],\
                'O3':[[-50,10],[-30,20],[-10,40]],'SO2':[[0,20],[0,10],[0,5]]}
    '''
    # For total concentration
    yLim     = {'NOx':[[0,200],[0,80],[0,25]],'NOy':[[0,210],[0,90],[0,30]],\
                'O3':[[-50,70],[-30,70],[-10,90]],'SO2':[[0,20],[0,10],[0,5]]}
    '''
    # Background from flt_summ.xls 
    c0       = {'NOx':[0.7,0.7,1.0],'NOy':[4.,6.,8.],'O3':[62.,70.,80.],'SO2':[1.,1.,1.]}
   
  if "tva_980826" in prjName:
    distance = [18, 27, 59] # 86, 93, 141]
    times    = [10.25, 10.75, 12.25, 12.5, 12.75]
    zSmp     = [465, 500, 659, 910, 819, 662]       
    zSmp2    = [465, 500, 659, 910, 819, 662]   
    # Background from flt_summ.xls ( Using tva_980826 values currently)
    c0        = {'NOx':[0.7,0.7,1.0],'NOy':[4.,6.,8.],'O3':[62.,68.,80.],'SO2':[1.,1.,1.]}
    xLim     = [[-15.,15],[-20,20],[-30,30]]
    yLim     = {'NOx':[[0,160],[0,80],[0,25]],'NOy':[[0,160],[0,90],[0,30]],\
                'O3':[[-50,10],[-30,15],[-10,40]],'SO2':[[0,15],[0,10],[0,5]]}    

  if "tva_990706" in prjName:
    distance = [11, 31, 65]#, 89]
    times    = [12.25, 13.0, 16.0, 16.75]#[12.5, 13.25, 16.0]#
    zSmp     = [501, 505, 500]#[505, 491, 448, 533]
    zSmp2     = [501, 505, 500]#, 533]
    xLim      = [[-10.,10],[-20,20],[-30,30]]
    yLim      = {'NOx':[[0,160],[0,80],[0,25]],'NOy':[[0,160],[0,90],[0,30]],\
                  'O3':[[-50,10],[-30,15],[-10,40]],'SO2':[[0,15],[0,6],[0,2.5]]}
    # Background from flt_summ.xls 
    #c0        = {'NOx':[0.7,0.7,1.0],'NOy':[4.,6.,8.],'O3':[62.,68.,80.],'SO2':[1.,1.,1.]}
    c0       = {'NOx':[0.7,0.7,1.0],'NOy':[4.,6.,8.],'O3':[55.,55.,60.],'SO2':[1.,1.,1.]}
    lblCt = 3

  if "tva_990715" in prjName:
    distance = [16, 62, 106]
    times    = [11.5, 12.25, 16.5]#[11.75, 12.5, 17]#
    zSmp     = [411, 727, 438]#[445, 589, 587]#
    zSmp2    = [415, 584, 582]
    # Background from flt_summ.xls 
    c0       = {'NOx':[0.7,0.7,1.0],'NOy':[4.,6.,8.],'O3':[62.,68.,80.],'SO2':[1.,1.,1.]}
    xLim     = [[-10.,10],[-20,20],[-30,30]]
    yLim     = {'NOx':[[0,160],[0,80],[0,25]],'NOy':[[0,160],[0,90],[0,30]],\
                'O3':[[-50,10],[-30,15],[-10,40]],'SO2':[[0,15],[0,10],[0,5]]}    
  
  statFile = open(prjName+"_stat.csv","w")
  #Normalized Mean Biasnmb,Mean Bias Error mbe, Mean Normalized Bias Error mnbe,Mean Normalized Gross Error (mnge) 
  statFile.write("Case, Distance, Species, Mean Observed (ppb), Mean Predicted(ppb), NMSE, MFB, r,  \
       NMB, MBE(ppb), MNBE, MNGE, NMAEF, NMBF\n")

  fList = os.listdir('.')
  for id,dist in enumerate(distance):
    pls = getTraverseList(prjName, dist)
    print 'Creating plots for %s at distance %s in %s'%(prjName,dist,os.getcwd())
    maxRatio =[]
    for iv,varName in enumerate(['NOx','NOy','O3','SO2']):
      print 'varName = %s'%varName       
      obsArray = []
      xMin =  1.e+20
      xMax = -1.e+20
      for ipl in pls:
        npzFileName = str(dist) +'km_' +  str(ipl) + '_' + varName
        print npzFileName
        for fName in fList:
          if fName.endswith('.npz'):
            if npzFileName in fName:
              npzfile = np.load(fName)#load(fName)
              obsData = npzfile['arr_0']
              obsData = obsData[obsData[:,0] > -9998.]
              xMin = min(xMin,min(obsData[:,0]))  
              xMax = max(xMax,max(obsData[:,0]))
              #print obsData[:,0]
              print 'Plume No. %d: xMin,Xmax = %g %g'%(ipl,min(obsData[:,0]),max(obsData[:,0]))
              obsArray.append(obsData)
              preData1 = npzfile['arr_1']
              #
              if preDirs is None:
                preData2  = npzfile['arr_2']
              else:
                preData2L = []
                for preDir in preDirs:
                  #print 'Reading %s'%os.path.join(preDir,prjName,'SCICHEM',fName)
                  npzfile2 = np.load(os.path.join(preDir,prjName,'SCICHEM',fName))
                  preData2L.append(npzfile2['arr_1'])                
               
              # Convert concentration from ppm to ppb
              obsData[:,1] = obsData[:,1]*1e3 
              preData1[:,1] = preData1[:,1]*1e3
              if preDirs is None:
                preData2[:,1] = preData2[:,1]*1e3
              else:
                for preData2 in preData2L:          
                  preData2[:,1] = preData2[:,1]*1e3
      
      if len(obsArray) == len(pls):
        xAvg = np.linspace(xMin, xMax, 40)
        #print 'xAvg = ',xAvg[0],xAvg[-1]
        avgData = np.zeros(len(xAvg))
        colors = ['b','g','c','k','y']
 
        avgList = []
        for ipl in range(len(pls)):
          f = interp1d(obsArray[ipl][:,0],obsArray[ipl][:,1],bounds_error=False,fill_value=-9999.)
          y = f(xAvg)
          avgList.append(y)
          if varName == 'NOx' and dist==18:
            y = ma.masked_where(y<-9998.0,y)
            #plt.plot(xAvg,y)
            #plt.show()
            pass
        avgList = np.array(avgList)
        avgList = ma.masked_where(avgList<-9998.0,avgList)
        npl,nx  = np.shape(avgList)
        avgData = np.zeros(nx) - 9999.0
        avgData = ma.masked_where(avgData<-9998.0,avgData)
        errData = avgData.copy()
        for ix in range(nx):
          nc = avgList[:,ix].count()
          if nc <= 1:
            continue
          else: 
            avgData[ix] = np.ma.mean(avgList[:,ix])
            errData[ix] = np.ma.std(avgList[:,ix])#/np.sqrt(nc)
        if rmAmb:           
          if varName == 'O3':
            c1 = max(preData1[-1,1],preData1[0,1])
            if preDirs is None:
              c2 = max(preData2[-1,1],preData2[0,1])
            else:
              c2L = []
              for preData2 in preData2L:
                c2L.append(max(preData2[-1,1],preData2[0,1]))
          else:
            c1 = min(preData1[-1,1],preData1[0,1])
            if preDirs is None:
              c2 = min(preData2[-1,1],preData2[0,1])
            else:
              c2L = []
              for preData2 in preData2L:
                c2L.append(min(preData2[-1,1],preData2[0,1]))
         
          #print varName, dist, c0[varName][id],avgData[:].max(),preData1[:,1].max()
        else:
          c1 = 0.
          if preDirs is None:
            c2 = 0.
          else:
            c2L = []
            for preData2 in preData2L:
              c2L.append(0.)
          #c0[varName][id] = avgData[0]+avgData[-1])
          
        oData         = ma.masked_all((nx,2))
        oData[:,0]    = xAvg
        oBgd          = c0[varName][id]
        oData[:,1]    = avgData[:] -  oBgd
        preData1[:,1] = preData1[:,1] - c1
        
 
        
        if preDirs is None:
          preData2[:,1] = preData2[:,1] - c2
        else:
          for i in range(len(preData2L)):
            preData2L[i][:,1] = preData2L[i][:,1] - c2L[i]
                      
        pData         = ma.masked_all((nx,2))
        pData[:,0]    = xAvg
        f = interp1d(preData1[:,0],preData1[:,1],bounds_error=False,fill_value=-9999.)
        pData[:,1]    = f(xAvg)
        pData[:,1]    = ma.masked_where(pData[:,1]<-9998.,pData[:,1])         
        if False:  # Single plots
          plt.clf()
          plt.hold(True)
          lgd = []          
          lgd.append('Avg. Plume')
          plt.plot(xAvg,avgData, '-',color='r',)
          plt.text(0.1,0.8,'%s'%varName,transform=ax.transAxes)
          plt.errorbar(xAvg, Cobs, yerr=errData, fmt='o',color='r')
          lgd.append('SCICHEM 3.0')
          plt.plot(preData1[:,0],preData1[:,1], 's')
          lgd.append('SCICHEM 1.0')
          plt.plot(preData2[:,0],preData2[:,1], '+')
          plt.xlim(xLim[id])
          plt.ylim(yLim[varName][id])
          plt.legend(lgd)
          plt.hold(False)
          #plt.show()
          print os.getcwd()
          plt.savefig(varName + '_' + str(dist) + 'km.tif', dpi=300)
        else:
          # Create subplots
          if iv == 0:
            print 'Start fig'
            fig = plt.figure()
            plt.clf()
            fig.hold(True)
            fig.subplots_adjust(bottom=0.2)
            ax = fig.add_subplot(111)
            #plt.title(' Downwind distance %s km '% dist,fontsize=11)
            plt.text(.35,-0.08,'Cross plume distance(km)',transform=ax.transAxes,fontsize=10)
            #plt.text(0.48,-0.14,lblForManus[lblCt],fontsize=13,transform=ax.transAxes)# Labels for Manuscript plots
            lblCt = lblCt+1
            runDate = prjName.replace('tva_','')
            figNo   = 99
            if "tva_980825" in prjName:
              runDate = '08/25/98'
              if dist == 20:
                 figNo = 6
              if dist == 55:
                 figNo = 7
            if "tva_990706" in prjName:
              runDate = '07/06/99'
              if dist == 11:
                 figNo = 8
              if dist == 31:
                 figNo = 9
              if dist == 65:
                 figNo = 10              
            figCaption1 = 'Figure %d:  Comparison of the mean observed and predicted concentrations for TVA cumberland plume '%figNo
            figCaption2 = 'traverse study at %d Km on %s. Error bars represent one standard deviation from the mean.'%(dist,runDate)                 
            #plt.text(-0.05,-0.17,figCaption1,transform=ax.transAxes,fontsize=10)
            #plt.text(-0.05,-0.23,figCaption2,transform=ax.transAxes,fontsize=10)
            print figCaption1
            print figCaption2
            plt.text(-.07,0.7,'Perturbation Concentration (ppb)',transform=ax.transAxes,rotation='vertical',fontsize=10)
            plt.setp(plt.gca(), frame_on=False, xticks=(), yticks=())
            
        
          #print varName, dist, oData[:,1].max(),preData1[:,1].max()  

          if figNo == 10 and varName=='O3':
            print oData
            print preData1
          if preDirs is None:
            fig,LhO,LhP1 = subplotConc(fig,oData,preData1,iv,varName,xLim=xLim[id],\
                                        yLim=yLim[varName][id],errData=errData)
            fig.legend([LhO,LhP1],['OBSERVED','SCICHEM-3.0'],ncol=2, bbox_to_anchor=(0.9,0.96))
          else:
            fig,LhList = subplotConc(fig,oData,preData1,iv,varName,xLim=xLim[id],\
                                        yLim=yLim[varName][id],errData=errData,preData2L=preData2L)
            LbList = ['Observed']
            for lbl in preLbls:
              LbList.append(lbl)
            fig.legend(LhList,LbList,ncol=4) #, bbox_to_anchor=(0.9,0.96))
        
          if iv == 3:
            
            
            fig.hold(False)
            figName = prjName + '_'  + str(dist) + '_km' + '_ResSen.png'
            #figName = 'Figure%d.eps'%figNo
            plt.savefig(figName,dpi=300, figsize=(8, 6))
            print 'Created %s in %s'%(figName,os.getcwd())
            #plt.show()
          statFile.write("%s, %02d, %s,"%(prjName,dist,varName))
          calcStats(oData+oBgd,pData+oBgd,varName,statFile=statFile)
          if varName == 'NOy':
             maxRatio = []
             maxRatio.append(oData[:,1].max()+oBgd)
             maxRatio.append(pData[:,1].max()+oBgd)
          if varName == 'SO2':      
            print  '***Ratio at' , dist ,': Obs, pred ',maxRatio[0]/(oData[:,1].max()+oBgd), maxRatio[1]/(pData[:,1].max()+oBgd)
            
            


            
def calcStats(obsArray, preArray,varName,statFile=None):
  if varName == 'O3':
    pass
  obsMean =  obsArray[:,1].mean()
  predMean = preArray[:,1].mean()
  print 'Min:',obsArray[:,1].min(),preArray[:,1].min()
  print 'Max:',varName, obsArray[:,1].max(),preArray[:,1].max()
  if min(obsArray[:,1])< 0.:
    pass
  NMSE = measure.nmse_1(preArray[:,1],obsArray[:,1])
  #NME = measure.nme(preArray[:,1],obsArray[:,1], cutoff=0.0)
  #pearCoeff = pearsonr(obsArray[:,1],preArray[:,1])
  biasFac = measure.bf(preArray[:,1],obsArray[:,1])
  MFB = measure.mfbe(preArray[:,1],obsArray[:,1],cutoff= -9998)
  mnbe = measure.mnbe(preArray[:,1],obsArray[:,1],cutoff=-9998)
  correlation = measure.correlation(preArray[:,1],obsArray[:,1])
  mage = measure.mage(preArray[:,1],obsArray[:,1])
  nmb = measure.nmb(preArray[:,1],obsArray[:,1])
  mbe = measure.mbe(preArray[:,1],obsArray[:,1])
  rmse = measure.rmse(preArray[:,1],obsArray[:,1])
  mnge = measure.mnge(preArray[:,1],obsArray[:,1],cutoff=-9998)
  nmaef = measure.nmaef(preArray[:,1],obsArray[:,1])
  nmbf  = measure.nmbf(preArray[:,1],obsArray[:,1])
  print 'stats: ',obsMean, predMean,';',NMSE, biasFac, correlation
  if statFile is not None:
    statFile.write("%8.2f, %8.2f, %8.2f, %8.2f,%8.2f, %8.2f,%8.2f,%8.2f,%8.2f,%8.2f,%8.2f,%8.2f,%8.2f\n"%\
                   (obsMean, predMean, NMSE, MFB, correlation, nmb, mbe, mnbe, mnge, nmaef, nmbf, obsArray[:,1].max(),preArray[:,1].max()))
  
def subplotConc(fig,obsData,preData1,iVar,varName,xLim=None,yLim=None,errData=None,preData2L=None):
  
  from scipy.integrate import simps
  
  isColor = True
  
  if isColor:
    clrs = ['r','g','b','m','c','k','y']
  else:
    clrs = ['0.25',  '0.5' ,  '0.75' ,'0.75', '0.5' ,'0.25'] 

  ax = fig.add_subplot(2,2,iVar+1)
  if varName == 'O3':
    plt.text(0.85,0.1,'%s'%varName,transform=ax.transAxes)
  else:
    plt.text(0.85,0.9,'%s'%varName,transform=ax.transAxes)
    
  if errData is not None:
    plt.errorbar(obsData[:,0],obsData[:,1], yerr=errData, fmt='o', color='k')

  
  LhO, = plt.plot(obsData[:,0],obsData[:,1],linestyle='None',marker='o',markersize=5,markerfacecolor=clrs[0])
  
  obsMask = obsData[:,1].mask
  obsX    = np.ma.array(obsData[:,0],mask=obsMask).compressed()
  obsY    = np.ma.array(obsData[:,1],mask=obsMask).compressed()
  IntObs = simps(obsY, obsX)
  
  if preData2L is not None:
    IntList = [IntObs]
    LhList = [LhO]
    lSty = ['-','--','-.',':','__','---','-']
    mSty = ['^','>','<','d','p','*','h']

    for i,preData2 in enumerate(preData2L): 
      LhP2, = plt.plot(preData2[:,0],preData2[:,1],linestyle=lSty[i],color=clrs[i+1],marker=mSty[i],markersize=5,markerfacecolor=clrs[i+1])
      LhList.append(LhP2)
      IntList.append(simps(preData2[:,1], preData2[:,0] + max(-preData2[0,0],0.) + 0.01))
  else:
     LhP1, = plt.plot(preData1[:,0],preData1[:,1],linestyle='None',color=clrs[2],marker='s',markersize=5,markerfacecolor=clrs[2])

  if xLim is not None:
    plt.xlim(xLim)
  if yLim is not None:
    plt.ylim(yLim)
  if iVar == 0 or iVar == 1:
    ax.set_xticklabels([])
    
  if preData2L is None:  
    return fig,LhO,LhP1
  else:
    print 'CrossWind Integrated concentration = ',IntList
    return fig,LhList
  
def getTraverseList(prjName, dist):
  pls =[]
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
      
  return pls

        
# Main program
if __name__ == '__main__':

  preDirs  = [None]
  preLbls  = [None]
  runTypes = [None]
  
  if compName == 'bnc' or compName == 'sage-d600':
    runDir = "D:\\TestSCICHEM\\Outputs\\2017_03_27\\Chemistry"  
    runDir = "D:\\TestSCICHEM\\Outputs"
    typeList = 'v3b2_TN3_noLSV_utc'
    #typeList = 'v3b2 v3b2_FixAmb v3b2_FixAmb_LoISOP v3b2_FixAmb_HiISOP v3b2_FixAmb_HiOHHO2  v3b2_FixAmb_LoOHHO2' 
    
    # Figure 6-11
    #typeList = 'v3b2_TN3'
    # Figure 12
    #typeList = 'v3b2_TN3 v3b2_TN10'
    # Figure 13
    #typeList = 'v3b2_FixAmb_HiOHHO2 v3b2_FixAmb v3b2_FixAmb_LoOHHO2'
    # Figure 14
    typeList = 'v3b2_TN3 v3b2_HiRes'
    
  if compName == 'pj-linux4':
    runDir = '/home/user/bnc/TestSCICHEM/Outputs'
    typeList = 'v3b3'
  if compName == 'Durga':
    runDir = "C:\\Users\\Bishusunita\\BNC\\TestSCICHEM\\Outputs"
    # Figure 6
    typeList = 'v3b2_TN3'
    #typeList = 'v3b2_TN3_eps'#'v3b2_FixAmb_HiOHHO2 v3b2_FixAmb v3b2_FixAmb_LoOHHO2 ' # v3b2  v3b2_TN10  v3b2_TN3 
    # HiLoOH v3b2_FixAmb_HiOHHO2  v3b2_FixAmb v3b2_FixAmb_LoOHHO2 
    # Layer v3b2_TN3 v3b2_TN10
    # Hi    v3b2_HiRes
  runTypes = typeList.split() 
  if runTypes[0] is not None:
    preDirs = []
    preLbls = []
    
    for runType in runTypes:
    
      preDirs.append(os.path.join(runDir,'%s'%runType,'Chemistry'))
  
      if runType == 'v3b2':
        preLbl = 'SCICHEM-3.0(VISTAS_WEST)'
      elif runType == 'v3b2_TN3':
        preLbl = 'SCICHEM-3.0' #(Tennessee Region with 3 Layer)'
      elif runType == 'v3b2_TN10':
        preLbl = 'SCICHEM-3.0(Tennessee Region with 10 Layer)'
      elif 'HiRes' in runType: # == 'v3b2_HiRes':
        preLbl = 'SCICHEM-3.0(High Resolution)'
      elif runType == 'v3b2_FixAmb':
        preLbl = ' SCICHEM-3.0 (Local Amb)'#Medium Range  OH & HO2 '#
      elif runType == 'v3b2_FixAmb_HiOHHO2':
        preLbl = 'High OH & HO2'
      elif runType == 'v3b2_FixAmb_LoOHHO2':
        preLbl = 'Low OH & HO2'
      else:    
        preLbl = 'SCICHEM-3.0(Local Amb ' + runType.replace('v3b2_FixAmb_','') + ')' 
      preLbls.append(preLbl)
    
  print preDirs
  print preLbls
  
  import time
  time.sleep(3)
  for prjName in ['tva_980825','tva_990706']: #,'tva_980825','tva_990706']:
    if runTypes[0] is not None:
      prjDir = os.path.join(preDirs[0],prjName,'SCICHEM')
    else:
      prjDir = os.path.join(runDir,prjName,'SCICHEM')
    os.chdir(prjDir)
    print 'prjDir = ',prjDir
    pltAvgData(prjName,preDirs=preDirs,preLbls=preLbls)
    
  print 'Done'
