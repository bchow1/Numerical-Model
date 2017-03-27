#!/bin/python
import os
import sys
import socket
import fileinput
import subprocess
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.cm as cm
from matplotlib import colors

compName = socket.gethostname()

# Local modules
if  compName == 'pj-linux4':
  pyDir = '/home/user/bnc/python'
  testDir = '/home/user/bnc/TestSCICHEM'
else:
  pyDir = 'C:\\Users\\sid\\python'
  testDir = 'V:\\TestSCICHEM\\Script\\Chemistry'
scriptDir = os.path.join(testDir,'script')

sys.path.append(pyDir)
sys.path.append(scriptDir)
    
import utilDb

def crtNtv(env,scipp,iniFile,prjName,tString,rmOut=True,ntvFile='temp.ntv',hSlice=-1,vSlice=None,spName='SO2'):

  if os.path.exists(ntvFile):

    print '***************************************'
    print 'File %s already exists. Skipping crtNtv'%ntvFile
    print '***************************************'
    return ntvFile

  else:

    inSCIpp = open('scipp.input','w') 
    inSCIpp.write('%s\n'%iniFile)
    inSCIpp.write('KE\n%s\n'%prjName)
    if hSlice > 0.:
      inSCIpp.write('Concentration\nComponents\n%s\nHorizontal Slice\nMean Value\n'%spName)
      # slice lower left; upper right; Slice height
      inSCIpp.write('%s\n\n\n%s\n'%(tString,hSlice))
      inSCIpp.write('NG\n')
      inSCIpp.write('%s\n'%ntvFile)
    else:
      inSCIpp.write('Concentration\nComponents\n%s\nVertical Slice\nMean Value\n'%spName)
      inSCIpp.write('%s\n'%tString)
      # slice starting point
      print 'slice starting point %s'%vSlice[0]
      inSCIpp.write('%s\n'%vSlice[0])
      # slice ending point
      print 'slice ending point %s'%vSlice[1]
      inSCIpp.write('%s\n'%vSlice[1])
      # slice  vertical range
      print 'slice vertical range %s'%vSlice[2]
      inSCIpp.write('%s\n'%vSlice[2])
      # slice no. of vertical points
      print 'slice no. of vertical points %s'%vSlice[3]
      inSCIpp.write('%s\n'%vSlice[3])
      inSCIpp.write('NG\n')
      inSCIpp.write('%s\n'%ntvFile)
    inSCIpp.close()
    scipp_inp = open('scipp.input','r')
    scipp_out = open('scipp.output','w')
    scipp_err = open('scipp.error','w')
    h = subprocess.Popen(scipp, env=env, bufsize=0, shell=False,stdin=scipp_inp, stdout=scipp_out,stderr=scipp_err)
    h.communicate()
    scipp_inp.close()
    scipp_out.close()
    scipp_err.close()
    #os.remove('scipp.input')
    #os.remove('scipp.error')
    if rmOut:
      os.remove('scipp.output')
  
  return ntvFile

def getCtlist(env,scipp,iniFile,prjName,timeList=None):

  if timeList is None:
   
    tString = ':00'
    ntvFile = crtNtv(env,scipp,iniFile,prjName,tString,hSlice=1000.,rmOut=False)
    usrTfile = '%s_UserTimes.txt'%prjName
    try:
      shutil.move('scipp.output',usrTfile)
    except EnvironmentError:
      raise
  
    if os.path.exists(usrTfile):
      timeList   = []
      timeLStart = 9999
      for line in fileinput.input(usrTfile):
        if 'Available Field Times' in line:
          timeLStart = fileinput.lineno() + 1        
        if fileinput.lineno() > timeLStart:
          if 'User time' in line or len(line.strip()) == 0:
            break
          timeList.append([line.split('(')[1].split(')')[0].strip(),float(line.split('(')[0].strip())])
      fileinput.close()
    
    if os.path.exists(ntvFile):
      os.remove(ntvFile)

  return timeList

def pltNtv(ntvFile,nSkip=13,isLog=True,hSlice=-1,isLatLon=False,figHold=False,cLim=None):

  inFile = ntvFile.replace('.ntv','')

  nodeData = np.loadtxt(ntvFile, skiprows=nSkip, usecols=(0, 1, 2, 5), dtype={'names':('nId','x','y','cmean'),\
                        'formats':('int','float','float','float')})
  triFile = ntvFile.replace('.ntv','.tri')
  triData  = np.loadtxt(inFile + '.tri', skiprows=1, usecols=(0, 1, 2, 3), dtype={'names':('tId','nIda','nIdb','nIdc'),\
                        'formats':('int','int','int','int')})
  
  nodeData = np.sort(nodeData,order='nId')
  triData  = np.sort(triData,order='tId')
  
  triangles = np.array([triData['nIda']-1,triData['nIdb']-1,triData['nIdc']-1])
  
  # Rotate for LatLon
  if isLatLon:
    y = nodeData['x']
    x = nodeData['y']
  else:
    x = nodeData['x']
    y = nodeData['y']
  c = nodeData['cmean'] * 1e3 # ppb
  
  #if not isLog:
  c0 = c[0]
  c  = c - c0
  
  npts = len(x)
  maxc = max(c)
  minc = min(c)
  
  print 'Min, Max Value = ',minc,maxc  
  print 'No. of points = ',npts

  if False:
    c = c/maxc
    scaled = True
  else:
    scaled = False

  if isLog:
    c = np.ma.masked_where(c<1e-30,c)
    c = np.ma.filled(c,1e-30)
  
  if isLog:
    logBase = 10.
    clrmax = int(np.log10(maxc)) - 1 # 0
    clrmin = clrmax - 4
    clrlev = (clrmax - clrmin + 1)
    levels = np.logspace(clrmin,clrmax,num=clrlev,base=logBase)
    clrmap = cm.get_cmap('jet',clrlev-1)
    lnorm  = colors.LogNorm(levels,clip=False)
  else:
    if cLim is None:
      clrmin = minc
      clrmax = maxc
      clrlev = 6
    else:
      clrmin = cLim[0]
      clrmax = cLim[1]
      clrlev = cLim[2]
    
    levels = np.linspace(clrmin,clrmax,num=clrlev)
    clrmap = cm.get_cmap('jet',clrlev-1)
    lnorm  = colors.Normalize(levels,clip=False)
  
  print levels

  # tricontour.
  if not figHold:
    fig = plt.figure()
    plt.clf()
    plt.hold(True)
  if isLog:
    if hSlice < 0:
      cax = plt.tricontourf(y,x,c, triangles=triangles, norm= lnorm, levels = levels, cmap=plt.cm.jet, vmin = 10**clrmin)
    else:
      cax = plt.tricontourf(x,y,c, triangles=triangles, norm= lnorm, levels = levels, cmap=plt.cm.jet, vmin = 10**clrmin)
  else:
    cax = plt.tricontourf(x,y,c, triangles=triangles, norm= lnorm, levels = levels, cmap=plt.cm.jet, vmin=clrmin)
  cax.cmap.set_under('white')
  
  if isLog:
    cax.set_clim(10**clrmin,10**clrmax)
    cbar = plt.colorbar(ticks=levels,format="%5.3e")
    cbar.ax.set_yticklabels(levels)
  else:
    cax.set_clim(clrmin,clrmax)
    cbar = plt.colorbar(ticks=levels,format="%3.1f")
    cbar.ax.set_yticklabels(levels)
  plt.tricontour(x,y,c, triangles=triangles, norm= lnorm, levels = levels, colors='k')
  
  if isLatLon:
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
  else:
    if hSlice < 0:
      plt.xlabel('Y')
      plt.ylabel('Z')
    else:
      plt.xlabel('X')
      plt.ylabel('Y')
      
  if not figHold:
    plt.hold(False)
    if scaled:
      plt.title('Concentrations for '+ inFile + ' scaled with %8.3e'%maxc)
    else:
      plt.title('Concentrations for '+ inFile)
    fig.savefig(inFile+'.png')
    return
  else:
    return lnorm,levels,clrmin,clrmax,cbar
    

if __name__ == '__main__':
  
  env = os.environ.copy()
  if os.sys.platform != 'win32':
    binDir      = "/home/user/bnc/EPRI/SCICHEM3.0B3/bin/linux"
    iniFile     = binDir #+ "/scipuff.ini"
    env["PATH"] = "%s:%s" % (binDir,env["PATH"])
    scipp       = os.path.join(binDir,'scipp')
    testDir   = '/home/user/bnc/TestSCICHEM'
  else:
    #binDir      = "D:\\EPRI\\SCICHEM3.0b2\\bin\\win32"
    binDir      = "D:\\SCIPUFF\\EPRI\\workspace\\EPRI\\vs2012\\bin\\intel\\x64\\Debug"
    vendir      = ""
    iniFile     = "D:\\EPRI\\SCICHEM3.0b2\\bin\\win32\\scipuff.ini"
    env["PATH"] = "%s;%s" % (binDir,vendir)
    scipp       = os.path.join(binDir,'scipp.exe')
    runsci      = os.path.join(binDir,'runSCI.exe')
    testDir   = 'V:\\TestSCICHEM'
 
  print 'Using scipp from ',binDir
  print scipp
  print env["PATH"]

  caseName = 'tva_980825'
  #caseName = 'tva_990706'
        
  matList   = {}
  ntvFile   = ''
  if sys.argv.__len__() > 1:
    runDir = os.getcwd()
  else:
    runDir = os.path.join(testDir,'Outputs','TVA_O3','Chemistry',caseName,'SCICHEM')
  os.chdir(runDir)
  #iniFile  = os.getcwd() # + "/scipuff.ini"
  print runDir

  # Start point, end point, vertical range, no. of points 
  verList   = [['65. -150.','65. 150.','400. 700.','200'],
               ['65. -150.','65. 150.','400. 700.','200'],
               ['65. -150.','65. 150.','400. 700.','200']]

  if sys.argv.__len__() > 2:
    prjName = sys.argv[2]
  else:
    prjName = None

  if caseName == 'tva_980825':

    if prjName is None:
      # Set project names for creating ntv files
      prjName  = 'tva_980825_v3b2_TN3'

    cScale   = 1. # Conc scaling
    isLatLon = False
    hrList   = ['11:30','12:45','14:45']
    htList   = ['520.','600','620']
    # Start point, end point, vertical range, no. of points 
    xlabel   = "X(Km)"
    ylabel   = "Y(Km)"

  if caseName == 'tva_990706':

    if prjName is None:
      prjName  = 'tva_990706_v3b2_TN3'
    hrList   = ['12:00','13:00','16:00']
    htList   = ['501.','505','500']

  # Get time concentration output time list for forward project
  timeList = getCtlist(env,scipp,iniFile,prjName,timeList=None)
  print 'Times in %s file:'%prjName
  for tl in timeList:
    print tl
  print

  # Create ntv files 
  for iHr in range(len(hrList)):
 
      tString = hrList[iHr]
      
      for spName in ['SO2','O3']:
      
        if spName == 'O3':
          isLog = False
        else:
          isLog = True
        # Horizontal Slice
        hSlice  = htList[iHr]
        ntvFile = prjName + '_hSlice_' + spName + '_' + tString.replace(':','_') +'hr.ntv'
        crtNtv(env,scipp,iniFile,prjName,tString,rmOut=False,ntvFile=ntvFile,hSlice=hSlice,spName=spName)
        pltNtv(ntvFile,hSlice=hSlice,isLog=isLog)
        
        '''
        # Vertical slice
        vSlice  = verList[iHr]
        print vSlice
        print 'Starting plot of %s ...'%ntvFile
        ntvFile = prjName + '_vSlice_' + spName + '_' + tString.replace(':','_') +'hr.ntv'
        crtNtv(env,scipp,iniFile,prjName,tString,rmOut=False,ntvFile=ntvFile,vSlice=vSlice,spName=spName)      
        pltNtv(ntvFile,hSlice=hSlice)
        '''
        
        print
        print ntvFile
        print

      
  print 'Done:-)'
