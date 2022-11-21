#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Run script for non-isothermal simulations of flow, DCR and enthalpy bilantion in macro-scale CF  
# -- Versions -- 
    # 1.0 -- run simulation with the given inlet temperature list (steady state) 

import os
import time
import sys
import numpy as np
import matplotlib.pyplot as plt

svVer = 13                                              # solver version
postProcesList  = ['CO','O2','CO2','T','p']					        # list of considered gas components
# outFile     = 'avgFielOut.dat'      					# outlet file for results
isOnWallIn0org =1 							            # 1 - isOnWallCoating and isInWallCoating fields are prepared in the 0.org folder, 0 - not prepared
nIterFlow    = 1001							            # number of iterations for the first run of flow simulation
# ~nIterFlow    = 5
nIterFlowThen = 800							            # number of iterations for other runs of flow simulations

# patchName = 'average'                                 	# different patch names for postProcesing - for kraken
patchName = 'areaAverage'                   		    # different patch names for postProcesing - my PC                     

# -- tempLst with solved temperatures
tempLst = np.linspace(350,850,5)
startTemp = tempLst[0]

# -- number of used processors					
nProc = 8					

# -- my pc or kraken?
computer = 'kraken'

# -- auxiliary func for parse of folder names - very ugly but working, 
#see code below for usage
def isInt(value):
  try:
    int(value)
    return True
  except ValueError:
    return False    

# -- time the code
startTime = time.time()

grphsOut = np.zeros((len(tempLst),5))
# -- foreach temp
for ind2 in range(len(tempLst)):
    temp = tempLst[ind2]						# current temperature

    # -- first run -- blockMesh etc. 
    if startTemp == temp:
        # -- clean of old cases
        #	os.system('./Allclean')
        # -- compile solver
        #os.system('./DCR_thermV%d/Allwclean'%svVer)
        #os.system('./DCR_thermV%d/Allwmake'%svVer)
        # --set number of proc
        print ('Updating decomposeParDict numberOfSubdomains to %d'%nProc)
        fileHn  = open('system/decomposeParDict','r')
        lineLst = fileHn.readlines()
        for ind in range(len(lineLst)):
            if lineLst[ind].find('numberOfSubdomains') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
                lineLst[ind] = 'numberOfSubdomains\t%d;\n'%nProc
        with open('system/decomposeParDict','w') as fileHn:
            fileHn.writelines(lineLst)
        os.system('./Allrun')
        # os.system('rm -r tempRes')
        # os.system('mkdir tempRes')
		
        # -- set number of iterations for first flow simulation in controlDict
        print('Updating controlDict endtime to %d'%(int(nIterFlow)))
        fileHn  = open('system/controlDict','r')
        lineLst = fileHn.readlines()
        for ind in range(len(lineLst)):
            if lineLst[ind].find('endTime') >= 0 and lineLst[ind].find('//') == -1 and lineLst[ind].find('stopAt') == -1:#!!do not put comment on executed temperature
                lineLst[ind] = 'endTime\t%d;\n'%int(nIterFlow)
            if lineLst[ind].find('writeInterval') >= 0 and lineLst[ind].find('//') == -1 and lineLst[ind].find('stopAt') == -1:#!!do not put comment on executed temperature
                spLine = lineLst[ind].split('\t')
                intSpLine = int(spLine[1].replace(';\n',''))
                lineLst[ind] = 'writeInterval \t %d;\n' %int(nIterFlow)
        with open('system/controlDict','w') as fileHn:
            fileHn.writelines(lineLst)
        dirsI = [0]

    # -- if not first run get current highest dir and set new endtime
    else:
        dirsI = []
        dirsS = os.listdir('./processor0/')
        for Dir in dirsS:
            if isInt(Dir):
                dirsI.append(int(Dir))
        print('Updating controlDict endtime to %d'%(int(nIterFlowThen+max(dirsI))))
        fileHn  = open('system/controlDict','r')
        lineLst = fileHn.readlines()
        for ind in range(len(lineLst)):
            if lineLst[ind].find('endTime') >= 0 and lineLst[ind].find('//') == -1 and lineLst[ind].find('stopAt') == -1:#!!do not put comment on executed temperature
                lineLst[ind] = 'endTime\t%d;\n'%int(nIterFlowThen+max(dirsI))
            if lineLst[ind].find('writeInterval') >= 0 and lineLst[ind].find('//') == -1 and lineLst[ind].find('stopAt') == -1:#!!do not put comment on executed temperature
                spLine = lineLst[ind].split('\t')
                intSpLine = int(spLine[1].replace(';\n',''))
                lineLst[ind] = 'writeInterval \t %d;\n' %int(nIterFlowThen+max(dirsI))
        with open('system/controlDict','w') as fileHn:
            fileHn.writelines(lineLst)

    # -- set temperature
    os.system('cp -r 0 %d'%max(dirsI))
    fileHn  = open('%d/T'%max(dirsI),'r')
    lineLst = fileHn.readlines()
    for ind in range(len(lineLst)):
        if lineLst[ind].find('internalField   uniform') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
            lineLst[ind] = 'internalField   uniform %g;\n'%temp
        if lineLst[ind].find('value           uniform') >= 0 and lineLst[ind].find('//') == -1:#!!do not put comment on executed temperature
            lineLst[ind] = 'value           uniform %g;\n'%temp
    with open('%d/T'%max(dirsI),'w') as fileHn:
        fileHn.writelines(lineLst)   
    os.system('decomposePar -time %d -fields >> log.decomposePar%g'%(max(dirsI),tempLst[ind2]))
    os.system('rm -r %d'%max(dirsI))


    # -- run simulation 
    os.system('./Allrun-parallel2')       
    os.system('mv log.flowDCRentV  log.flowDCRentV%g'%tempLst[ind2])

    # -- postProcess
    os.system("reconstructPar -fields '(CO T O2 U CO2 p)' -latestTime >> log.reconstructPar%g"%tempLst[ind2])         # reconstruct Case
    dirsI = []
    dirsS = os.listdir('./processor0')
    for Dir in dirsS:
        if isInt(Dir):
            dirsI.append(int(Dir))
    os.system('cp -r %d %g'%(max(dirsI),tempLst[ind2]))	
    os.system('rm -r %d'%max(dirsI))
    
    # --postprocessing reading patches
    os.system("postProcess -func 'patchAverage(name=outlet,CO)' -latestTime > log.patchAverageOutletCO")
    os.system("postProcess -func 'patchAverage(name=outlet,T)' -latestTime > log.patchAverageOutletT")
    os.system("postProcess -func 'patchAverage(name=outlet,O2)' -latestTime > log.patchAverageOutletO2")
    os.system("postProcess -func 'patchAverage(name=outlet,CO2)' -latestTime > log.patchAverageOutletCO2")
    os.system("postProcess -func 'patchAverage(name=inlet,p)' -latestTime > log.patchAverageOutletp")
    # os.system("cp %d/T tempRes/%g_T_field"%(max(dirsI),temp))

    outConcLst = []
    # ~ fileNm = 'log.patchAverageInletP'
    # ~ fileHn = open(fileNm,'r')
    # ~ lineLst= fileHn.readlines()
    # ~ for ind in range(len(lineLst)):
        # ~ if lineLst[ind].find('%s'%patchName) >= 0:
            # ~ outConcLst.append(float(lineLst[ind].split(' ')[-1]))
    # ~ fileHn.close()        
    for ind3 in range(len(postProcesList)):
        fileNm = 'log.patchAverageOutlet%s'%(postProcesList[ind3])
        fileHn = open(fileNm,'r')
        lineLst= fileHn.readlines()
        for ind in range(len(lineLst)):
            if lineLst[ind].find('%s'%patchName) >= 0:
                grphsOut[ind2,ind3] = (float(lineLst[ind].split(' ')[-1]))
        fileHn.close()

print(grphsOut)
np.save('grphsOut.npy',grphsOut)
plt.plot(tempLst-273,1e6*grphsOut[:,0],label='CO')
plt.plot(tempLst-273,1e6*grphsOut[:,1],label='O2')
plt.plot(tempLst-273,1e6*grphsOut[:,2],label='CO2')
plt.legend()
plt.title('Outlet yi in dependence on inlet temperature')
plt.show()
plt.close()
plt.plot(tempLst-273,grphsOut[:,3]-tempLst)
plt.title('Temperature difference in dependece on inlet temperature')
plt.show()
plt.close()
plt.plot(tempLst-273,grphsOut[:,4])
plt.title('Inlet pressure in dependece on inlet temperature')
plt.plot()
plt.show()

#os.system('reconstructPar -newTimes')

endTime = time.time()
print('Total program execution time: %7.2f s\n'%(endTime-startTime))

