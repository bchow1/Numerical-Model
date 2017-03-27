#!/bin/python

#C:\Windows\System32\cmd.exe /K C:\Python27\python.exe $(FULL_CURRENT_PATH) 
import os
import sys
import socket
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sqlite3


sys.path.append('/home/user/bnc/python')


compName = socket.gethostname()
print compName

if compName == 'Durga': 
  #os.chdir('C:\\Users\\Bishusunita\\BNC\\Publications\\SCICHEM\\AtmEnv\\Reviewer')
  sys.path.append('C:\\Users\\sid\\python')
  # C:\Users\Bishusunita\BNC\TestSCICHEM\Outputs\150309_HiRes\Chemistry\tva_990706\SCICHEM
  valDir = 'C:\\Users\\Bishusunita\\BNC\\TestSCICHEM'#C:\Users\Bishusunita\BNC\TestSCICHEM\Inputs\Chemistry\
  ymdhm='v3b2_TN3'
  #prjName in  ['tva_980825']: # ['tva_980825','tva_980826', 'tva_990706', 'tva_990715']
  prjName = 'tva_980825' #'tva_990706'
  samFile = '%s.sam'%prjName
  smpFile =   '%s.smp'%prjName
  print os.getcwd()
#  ymdhm=os.getcwd().split('\')[7]
#  prjName=os.getcwd().split('\')[9]
  
if compName == 'sm-bnc' or compName == 'sage-d600':
  sys.path.append('C:\\Users\\sid\\python')
  valDir = 'V:\\TestSCICHEM'
  ymdhm='TVA_O3'
  #prjName in  ['tva_980825']: # ['tva_980825','tva_980826', 'tva_990706', 'tva_990715']
  
  
  #prjName = 'tva_990706_v3b2_TN3'
  #prjDir = 'tva_990706'
  

  prjName = 'tva_980825_v3b2_TN3'
  prjDir = 'tva_980825'
      
  samFile = '%s.sam'%prjName
  smpFile =   '%s.smp'%prjName
  print os.getcwd()  
  
if  compName == 'pj-linux4':
  sys.path.append('/home/user/bnc/python')
  valDir = '/home/user/bnc/TestSCICHEM'
  #ymdhm=os.getcwd().split('/')[7]
  #prjName=os.getcwd().split('/')[9]
  ymdhm='v3b3'
  prjName = 'tva_980825'
  prjDir = 'tva_980825'
  samFile = '%s.sam'%prjName
  smpFile = prjName + '.smp'
  
  #smpFile = sys.argv[1] + '.smp'
  if not os.path.exists(smpFile):
    print 'Error: Cannot find smpFile %s'%smpFile
    sys.exit()
ambFile = smpFile.replace('.smp','.asmp')
if not os.path.exists(ambFile):
  ambFile = None

# Local modules    
import convertUTM
import utilDb
sys.path.append(os.path.join(valDir,'Scripts','Chemistry'))
import pltSO2Slice

def getPlConc(slatlon,prjName,dist,ipl,varName):

  myConvertUTM = convertUTM.convertUTM()
  myConvertUTM.latd,myConvertUTM.lngd = slatlon
  print myConvertUTM.latd,myConvertUTM.lngd
  myConvertUTM.GeogToUTM()
  srcX,srcY = myConvertUTM.x*1e-3,myConvertUTM.y*1e-3
  print srcX,srcY

  dbDir   = os.path.join(cheDir,prjDir,'OBSERVATIONS')
  obsPfx  = os.path.join(dbDir,'%s_')%prjDir
  print prjDir

  obsList = []
  for pNo in ipl:
    obsDbName = obsPfx + str(dist) + 'km_obs' + str(pNo) + '.csv.db'
    print '\n obsDbName = ',obsDbName
        
    obsConn = sqlite3.connect(obsDbName)
    obsConn.row_factory = sqlite3.Row
    obsCur = obsConn.cursor()

    # Observations
    if "tva_990706" in prjName:
      #obsQry = 'select CAST(plumeKM as real) as x, avg(' + varName + '/1000.) from dataTable where plumeKM > -9999.0 group by x'
      obsQry = 'select Lat,Long,'+varName+'/1000. from dataTable'
    else:
      obsQry = 'select Lat,Long,'+varName+' from dataTable where Lat != -9999.0'
    print obsQry
    obsArray = utilDb.db2Array(obsCur,obsQry)
    print np.shape(obsArray),obsArray[0,0],obsArray[0,1],obsArray[0,2]
    obsList.append(obsArray)

  for pNo,obsArray in enumerate(obsList):
    for i in range(len(obsArray)):
      myConvertUTM.latd =  obsArray[i,0]
      myConvertUTM.lngd = -obsArray[i,1]
      #print  myConvertUTM.latd, myConvertUTM.lngd
      myConvertUTM.GeogToUTM()
      x = myConvertUTM.x*1e-3 - srcX
      y = myConvertUTM.y*1e-3 - srcY
      obsArray[i,0] = x
      obsArray[i,1] = y
    obsList[pNo] = obsArray
    print pNo,obsList[pNo][0,:]
    print pNo,obsList[pNo][-1,:]
  obsCur.close()
  obsConn.close()
  
  return obsList

cheDir = os.path.join(valDir,'Inputs','Chemistry')
outDir = os.path.join(valDir,'Outputs',ymdhm,'Chemistry')
pltDir = os.path.join(valDir,'Outputs',ymdhm,'Plots')

#
inpDir = os.path.join(cheDir,prjName,'SCICHEM')
runDir = os.path.join(outDir,prjDir,'SCICHEM')
# samFile = os.path.join(inpDir,samFile)
samFile = '%s.sam'%prjDir
#
os.chdir(runDir)
print os.getcwd()
myFile = open(prjName+"_Plume.csv","w")
#Normalized Mean Biasnmb,Mean Bias Error mbe, Mean Normalized Bias Error mnbe,Mean Normalized Gross Error (mnge) 

# Step 1: Find start and end line number for Arc@31km (78-?) from sam file
# Step 2: Find column numbers of C004_077 to C004_??? from smp file
# Step 3: Get the values for a given time i.e row (T = Column 1 in smp) for these columns
# Step 4: get x for each sampler . find values for a particular T.
#Step 5: Plot

srcLoc  = [0.,0.]
lineStyles = ['-','--',':','__','---']
markStyles = ['o','s','^','>','<']


if "tva_990706" in prjName:
  slatlon = (36.39,-87.6523)
  distList = [11,31,65]
  xlims   = [-20.,100.]
  ylims   = [-60.,60.]
  xFigCap = 82
  yFigCap = 52
  
  # Background from flt_summ.xls 
  c0   = {'NOx':[0.7,0.7,1.0],'NOy':[4.,6.,8.],'O3':[55.,55.,60.],'SO2':[0.631,1.,1.]}
  cLim = {'NOx':[[0,160],[0,80],[0,25]],'NOy':[[0,160],[0,90],[0,30]],\
                  'O3':[[-50,50,6],[-50,50,6],[-50,50,6]],'SO2':[[0,15],[0,6],[0,2.5]]}
  #cLim = {'O3':[[-50,50,6]]}  

if "tva_980825" in prjName:
  slatlon = (36.39,-87.6523)
  distList = [20,55,110]
  xlims   = [-100.,130.]
  ylims   = [-100.,100.]
  xFigCap = 95
  yFigCap = 85
  # Background from flt_summ.xls 
  c0   = {'NOx':[0.7,0.7,1.0],'NOy':[4.,6.,8.],'O3':[62.,68.,80.],'SO2':[1.,1.,1.]}
  cLim = {'NOx':[[0,160],[0,80],[0,25]],'NOy':[[0,160],[0,90],[0,30]],\
                'O3':[[-50,50,6],[-50,50,6],[-50,50,6]],'SO2':[[0,15],[0,10],[0,5]]}  

  verList   = [['65. -150.','65. 150.','400. 700.','200'],
               ['65. -150.','65. 150.','400. 700.','200'],
               ['65. -150.','65. 150.','400. 700.','200']]

for dNo,dist in enumerate(distList):
  # tva_990706 
  if dist == 11:
    hr      = 12.5
    dstr    = 'Arc@11.0km'
    ipl     = [1,2,3,4,5]
    xLims =   [-10.,10.]
    o3Lims =  [ -50.,10.]
    so2Lims = [ 0.,14]
    tString  = '12_00'  # For hSlice 
    zH       = 501. 
    figCapS  = 'a'
    figCapO = 'd'
  elif dist == 31:
    hr      = 13.25
    dstr    = 'Arc@31.0km'
    ipl     = [6,7,8]
    xLims =   [-20.,20.]
    o3Lims =  [ -30.,15.]
    so2Lims = [ 0.,6.]
    tString  = '13_00'  # For hSlice 
    zH       = 505. 
    figCapS  = 'b'
    figCapO  = 'e'
  elif dist == 65:
    hr      = 16.0
    dist    = 65
    dstr    = 'Arc@65.0km'
    ipl     = [9,10,11]
    xLims =   [-30.,30.]
    o3Lims =  [ -10.,40.]
    so2Lims = [ 0.,2.5]
    tString  = '16_00'  # For hSlice 
    zH       = 500. 
    figCapS  = 'c'
    figCapO  = 'f'
  #"tva_980825" 
  elif dist == 20:
    dstr    = 'Arc@20.0km'
    xLims =   [-15.,15.]
    o3Lims =  [ -50.,10.]
    so2Lims = [ 0.,20.]
    hr      = 11.5
    zH      = 520. 
    ipl     = [3,4,5]
    tString  = '11_30'  # For hSlice    
    figCapS  = 'a'
    figCapO = 'd'
  elif dist == 55:
    dstr    = 'Arc@55.0km'
    xLims =   [-25.,25.]
    o3Lims =  [ -30.,20.]
    so2Lims = [ 0.,10.]
    hr      = 12.75
    zH      = 600.
    ipl     = [6,7,8]
    # Start point, end point, vertical range, no. of points
    tString  = '12_45'  # For hSlice 
    figCapS  = 'b'
    figCapO = 'e'
  elif dist == 110:
    dstr    = 'Arc@110.0km'
    xLims =   [-30.,30.]
    o3Lims =  [ -10.,40.]
    so2Lims = [ 0.,5.]
    hr      = 14.75
    zH      = 620.
    ipl     = [9,10,11]
    tString  = '14_45'  # For hSlice
    figCapS  = 'c'
    figCapO = 'f'
  else:
    print 'Error in dist value'
    sys.exit()

  varName = 'SO2'
  obsSO2List = getPlConc(slatlon,prjName,dist,ipl,varName)

  varName = 'O3'
  obsO3List = getPlConc(slatlon,prjName,dist,ipl,varName)

  if '990706' in prjName:
    samVals = pd.read_table(samFile,skiprows=1,sep=r'\s*',names=['x','y','z','mc','spList','Dist'],dtype={'x':np.float,'y':np.float,'z':np.float})
    sams    = samVals[samVals['Dist']== dstr]
  else:
    samVals = pd.read_table(samFile,skiprows=1,sep=r'\s*',names=['x','y','z','mc','spList'],dtype={'x':np.float,'y':np.float,'z':np.float})
    sams    = samVals[abs(samVals['z']-zH) < 1.]

  #
  print 'samVals: %d,(%g,%g),(%g,%g)'%(len(samVals),samVals['x'].iloc[0],samVals['y'].iloc[0],samVals['x'].iloc[-1],samVals['y'].iloc[-1])
  samRows = sams.index
  print samRows[0],samRows[-1]
  print samVals.iloc[samRows[0]],samVals.iloc[samRows[-1]]

  spList = samVals['spList'][0].split('(')[1].split(')')[0].split(',')
  print 'Species list: ',spList

  # Add 1 for starting from 0 
  beg_smp = samRows[0]  + 1  
  end_smp = samRows[-1] + 1

  print 'Start and end sam: %d,%d,%d,(%g,%g),(%g,%g)'%(beg_smp,end_smp,len(sams),sams['x'].iloc[0],sams['y'].iloc[0],sams['x'].iloc[-1],sams['y'].iloc[-1])

  sR = len(samVals) + 2
  smpVals = pd.read_table(smpFile,skiprows=sR,sep=r'\s*',header=0)
  print smpVals.columns.values
  if ambFile is not None:
    ambVals = pd.read_table(ambFile,skiprows=sR,sep=r'\s*',header=0)
    print ambVals.columns.values
    

  spNos = {}
  spNo  = 0
  for cNo in smpVals.columns:
    print cNo
    if '_001' in cNo:
      print spNo, spList[spNo]
      spNos.update({spList[spNo]:cNo.replace('_001','')})
      spNo += 1
    if '_002' in cNo:
      break
  print spNos

  print 'Time  = [',smpVals['T'].values[0],', ',smpVals['T'].values[-1],']'
  tRow = smpVals[smpVals['T'] == 3600.*hr].index[0]
  print 'Row No. = %d in smpFile for hr %g'%(tRow,smpVals['T'][tRow]/3600.)

  SO2Cols = [ '%s_%03d'%(spNos['SO2'],smpNo) for smpNo in range(beg_smp,end_smp + 1) ]
  SO2 = smpVals[SO2Cols].values[tRow]*1e3  # ppb?
  if ambFile is not None:
    aSO2 = ambVals[SO2Cols].values[tRow]*1e3  # ppb?

  O3Cols = [ '%s_%03d'%(spNos['O3'],smpNo) for smpNo in range(beg_smp,end_smp + 1) ]
  O3 = smpVals[O3Cols].values[tRow]*1e3 

  # Find index for maximum SO2 conc
  print
  print 'Predicted max:'
  mxId = np.where(SO2 == SO2.max())[0][0]
  xPre = sams['x'].iloc[mxId]
  yPre = sams['y'].iloc[mxId]
  sPre = np.sqrt((sams['x']-xPre)**2 + (sams['y']-yPre)**2)*np.sign(sams['y']-yPre)
  print 'nSmp, len(SO2), SO2.max, maxId, xPre, yPre = ',len(sams['x']),len(SO2),SO2.max(),mxId,xPre,yPre
  print
  print 'X   Y  Dist SO2' 
  for i in range(len(SO2)):
    if abs(SO2[i]) > 0. :
      print sams['x'].iloc[i],sams['y'].iloc[i],sPre.iloc[i],SO2[i]

  # Find index for maximum obs conc
  xObs = []
  yObs = []
  for ip,pNo in enumerate(ipl):
    print
    print 'Observed max for plume No:',pNo
    obsSO2Array = obsSO2List[ip]
    obsO3Array = obsO3List[ip]
    mxObsId = np.where(obsSO2Array[:,2] == obsSO2Array[:,2].max())[0][0]
    xo = obsSO2Array[mxObsId,0]
    yo = obsSO2Array[mxObsId,1]
    xObs.append(xo)
    yObs.append(yo)
    print 'len(Obs), CObs.max,Obs(maxId) = ',len(obsSO2Array),obsSO2Array[:,2].max(),mxObsId,obsSO2Array[mxObsId,:]
    print 'Diff between observed and predicted max:',np.sqrt((xPre-xo)**2 + (yPre-yo)**2)

  # Plot SO2 and O3 plume centerlines
  if True:

    # Plot for SO2
    isLog = True
    plt.figure()
    plt.clf()
    plt.hold(True)
    
    # Contour plots
    ntvFile = prjName + '_hSlice_' + 'SO2_' + tString +'hr.ntv'
    print os.getcwd()
    print ntvFile,obsSO2Array[:,2].max()*1.e3,SO2.max() # Set obs unit to ppb
    lnorm,levels,vmin,vmax,cbar = pltSO2Slice.pltNtv(ntvFile,hSlice=int(zH),figHold=True,isLog=isLog)
    
    # Obs points
    cBck = obsSO2Array[:,2].min()
    cSO2 = (obsSO2Array[:,2] - cBck)*1e3 # ppb
    sax = plt.scatter(obsSO2Array[:,0],obsSO2Array[:,1],c=cSO2,edgecolors='none',norm=lnorm,\
                      cmap=plt.cm.jet,marker='s',s=20)
                      
    # Source
    plt.scatter((srcLoc[0],),(srcLoc[1],),color='k',marker='o',s=20)
  
    if not isLog:
      sax.set_clim(vmin,vmax)
    else:
      sax.set_clim(10.**vmin,10.**vmax)
  
    lhSets = []
    lgSets = []
    
    # Plume center line for observed plume
    for ip,pNo in enumerate(ipl[-1:]):
      obsSet, = plt.plot((srcLoc[0],xObs[ip]),(srcLoc[1],yObs[ip]),color='r',linestyle='-')
  
    # Plume center line for predicted plume
    preSet, = plt.plot((srcLoc[0],xPre),(srcLoc[1],yPre),color='g')
    
    lhSets.append(obsSet)
    lhSets.append(preSet)
    lgSets.append('Observed SO2')
    lgSets.append('Predicted  SO2')
    plt.xlim(xlims)
    plt.ylim(ylims)
    plt.legend(lhSets,lgSets,ncol=1,bbox_to_anchor=(0.02,0.98),loc=2,borderaxespad=0.)
    plt.xlabel('X(km)')
    plt.ylabel('Y(km)')
    plt.text(xFigCap, yFigCap, '( %s )'%figCapS, fontsize =14)
    plt.title('Plume centerlines at %d km'%dist)
    cbar.ax.text(0.1,1.03,'ppb',fontsize=11)
    plt.hold(False)
    plt.savefig('%s_Cline_SO2_%dkm.png'%(prjName,dist))
    plt.show()

    # Plot for O3
    isLog = False
    plt.figure()
    plt.clf()
    plt.hold(True)
    
    # Contour plots
    ntvFile = prjName + '_hSlice_' + 'O3_' + tString +'hr.ntv'
    print os.getcwd()
    print ntvFile,obsO3Array[:,2].max()*1.e3,O3.max() # Set obs unit to ppb
    lnorm,levels,vmin,vmax,cbar = pltSO2Slice.pltNtv(ntvFile,hSlice=int(zH),figHold=True,isLog=isLog,cLim=cLim['O3'][dNo])
    
    # Obs points
    cBck = c0['O3'][dNo]
    cO3 = obsO3Array[:,2]*1e3 - cBck # ppb
    sax = plt.scatter(obsO3Array[:,0],obsO3Array[:,1],c=cO3,edgecolors='none',norm=lnorm,\
                      cmap=plt.cm.jet,marker='s',s=20)
                      
    # Source
    plt.scatter((srcLoc[0],),(srcLoc[1],),color='k',marker='o',s=20)
  
    if not isLog:
      sax.set_clim(vmin,vmax)
    else:
      sax.set_clim(10.**vmin,10.**vmax)
  
    lhSets = []
    lgSets = []
    
    # Plume center line for observed plume
    for ip,pNo in enumerate(ipl[-1:]):
      obsSet, = plt.plot((srcLoc[0],xObs[ip]),(srcLoc[1],yObs[ip]),color='r',linestyle='-')
  
    # Plume center line for predicted plume
    preSet, = plt.plot((srcLoc[0],xPre),(srcLoc[1],yPre),color='g')
    
    lhSets.append(obsSet)
    lhSets.append(preSet)
    lgSets.append('Observed O3')
    lgSets.append('Predicted O3')
    plt.xlim(xlims)
    plt.ylim(ylims)
    plt.legend(lhSets,lgSets,ncol=1,bbox_to_anchor=(0.02,0.98),loc=2,borderaxespad=0.)
    plt.xlabel('X(km)')
    plt.ylabel('Y(km)')
    plt.text(xFigCap, yFigCap, '( %s )'%figCapO, fontsize =14)
    plt.title('Plume centerlines at %d km'%dist)
    cbar.ax.text(0.1,1.03,'ppb',fontsize=11)
    plt.hold(False)
    plt.savefig('%s_Cline_O3_%dkm.png'%(prjName,dist))
    plt.show()


    '''
    plt.figure()
    plt.clf()
    plt.hold(True)
    
    # Vertical slice contour plots
    ntvFile = prjName + '_v3b2_TN3_vSlice_' + tString +'hr.ntv'
    print ntvFile,obsArray[:,2].max()*1.e3,SO2.max() # Set obs unit to ppb
    lnorm,levels,vmin,vmax = pltSO2Slice.pltNtv(ntvFile,figHold=True,isLog=not isLinear)
    #plt.xlim(xlims)
    #plt.ylim(ylims)
    plt.title('Vertical concentration slice for %s at %d km'%(prjName,dist))
    #ntvFile = prjName + '_vSlice_' + tString +'hr.ntv'
    plt.xlabel('X(km)')
    plt.ylabel('Z(m)')
    plt.hold(False)
    plt.savefig('%s_vSlice_%dkm.png'%(prjName,dist))
    plt.show()
    '''

  # Plot SO2 and O3 plume cross section
  if False:
    O3Cols = [ '%s_%03d'%(spNos['O3'],smpNo) for smpNo in range(beg_smp,end_smp + 1) ]
    O3 = smpVals[O3Cols].values[tRow]*1e3
    print 'Predicted O3 = ',O3
    if ambFile is not None:
      aO3 = ambVals[O3Cols].values[tRow]*1e3  # ppb?
      print 'Predicted amb O3 = ',aO3
      O3  = O3 - aO3
      SO2 = SO2 - aSO2

    figName = smpFile.replace('.smp','')
    fig = plt.figure()
    plt.clf()
    plt.hold(True)
    plt.setp(plt.gca(),frame_on=False,xticks=(),yticks=())
    ax = fig.add_subplot(1,2,2)
    plt.scatter(sPre,SO2,c='g',marker='o')
    print 'Predicted SO2 = ',SO2
    plt.xlim(xLims)
    plt.ylim(so2Lims)
    plt.text(0.85,0.1,'SO2',transform=ax.transAxes) 
    ax = fig.add_subplot(1,2,1)
    plt.scatter(sPre,O3,c='r',marker='o')
    plt.xlim(xLims)
    plt.ylim(o3Lims)
    plt.text(0.85,0.1,'O3',transform=ax.transAxes) 
    #plt.legend([so2set,o3set],['SO2','O3'],bbox_to_anchor=(0.02,0.98),loc=2,borderaxespad=0.)
    plt.title('Cross section for %s at %d km'%(figName,dist))
    plt.hold(False)
    plt.savefig('%s_crs_%dkm.png'%(figName,dist))
    plt.show()

print ('DONE')
