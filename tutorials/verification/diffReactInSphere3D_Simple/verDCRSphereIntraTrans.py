# -- Script to do verification simulation of isothermal reaction inside sphere

# -- imports
from re import X
from time import time
import numpy as np
import os 
import shutil as sh
import matplotlib.pyplot as plt

# -- auxiliary function to determine if val is float
def isFloat(val):
    try:
        float(val)
        return True
    except:
        return False


# -- setup case
baseCaseDir = 'baseCase'
outFolder = 'ZZ_cases'

# -- case data
yInf = 0.01     # molar fraction farfar from sphere
p = 101325      # presure
Runiv = 8.314   # universal gas constant
R = 1           # sphere radius
Rinf = 1.1         # infinite radius 
cellSizeLst = [0.2,0.1,0.05]  # cell Size
TLst = [500]  # temperature 
k0 = 5e7      # reaction pre-exponential factor
EA = 90e3     # reaction activation energy
tort = 3      # tortuosity
solver = 'reactingHetCatSimpleFoam' # used solver

# -- create case
for T in TLst:
    for cellSize in cellSizeLst:
        # -- create caseFolder based on baseCase
        caseName = 'intraTrasn_yInf_%g_R_%g_Rinf_%g_T_%g_tort_%g_cS_%g'%(yInf,R,Rinf,T,tort,cellSize)
        print('Preparing case %s'%caseName)
        caseDir = '%s/%s/'%(outFolder,caseName)

        if os.path.isdir(caseDir):                                          #ensure, that the caseDir is clear
            sh.rmtree(caseDir)

        sh.copytree(baseCaseDir,caseDir)


        # -- change the case files
        # NOTE TH: I only change temperature and domain size, much more can be studied 
        print('Changing temperature to %g'%T)
        with open(caseDir + '/0.org/T', 'r') as file:
            data = file.readlines()

        for ind in range(len(data)):
            if data[ind].find('isoT') >= 0:
                data[ind] = data[ind].replace('isoT',str(T))

        with open(caseDir + '/0.org/T', 'w') as file:
            file.writelines(data)
        
        # NOTE TH: I only change temperature, much more can be studied 
        print('Changing blockMeshDict cube size to (%gx%gx%g)'%(Rinf*2,Rinf*2,Rinf*2))
        with open(caseDir + '/system/blockMeshDict', 'r') as file:
            data = file.readlines()

        for ind in range(len(data)):
            if data[ind].find('dSN') >= 0:
                data[ind] = data[ind].replace('dSN',str(Rinf))
                
            if data[ind].find('nDisc') >= 0:
                data[ind] = data[ind].replace('nDisc',str(int(Rinf/cellSize)*2))

        with open(caseDir + '/system/blockMeshDict', 'w') as file:
            file.writelines(data)
        
        # -- run the simulation
        os.chdir(caseDir)
        os.system('./AllrunIntraSphere')

        # -- save results
        # -- load Dfree, Deff, and k from log file
        print('Reading parameters from log file.')
        # NOTETH: this can be slow (in that case maybe use something like tail?)
        with open('log.%s'%solver, 'r') as fl:
            lines = fl.readlines()
        # -- go from back to get latest result
        #DFree
        for lineInd in range(len(lines)-1,0,-1):
            if lines[lineInd].find('DFree') >= 0:
                DFree = float(lines[lineInd].split(' ')[-1].replace('\n',''))
                print('DFree = %g'%DFree)
                break
            if lineInd == 0:
                print('DFree not found.')
        #DEff
        for lineInd in range(len(lines)-1,0,-1):
            if lines[lineInd].find('DEff') >= 0:
                DEff = float(lines[lineInd].split(',')[0].split(' ')[-1].replace('\n',''))
                print('DEff = %g'%DEff)
                break
            if lineInd == 0:
                print('DEff not found.')
        #k
        for lineInd in range(len(lines)-1,0,-1):
            if lines[lineInd].find('max(k)') >= 0:
                k = float(lines[lineInd].split(' ')[-1].replace(').\n',''))
                print('k = %g'%k)
                break
            if lineInd == 0:
                print('k not found.')

        # -- compute analytical results
        # NOTETH: this is probably wrong
        kL = DFree/(Rinf-R)         # mass transfer coefficient (absolutely not sure about this)

        thiele = R*(k/DEff)**(0.5)  # thiele modulus
        # BiM = kL*R/DEff             # Biot number
        # etaAnal = 3./(thiele**2) * (thiele*1./np.tanh(thiele)-1)/(1+(thiele*1./np.tanh(thiele)-1)/BiM)                      # analytical effectivness factor
        etaAnal = 3./(thiele**2) * (thiele*1./np.tanh(thiele) - 1)                      # analytical effectivness factor
        print('\nThiele modulus = %g\nAnalytical effectivness factor is %g'%(thiele,etaAnal))

        # -- compute simulation results
        # -- read real source
        with open('log.intSrcSphere', 'r') as fl:
            lines = fl.readlines()
        #reaction source
        for lineInd in range(len(lines)-1,0,-1):
            if lines[lineInd].find('reaction source') >= 0:
                rS = float(lines[lineInd].split(' ')[-1].replace('\n',''))*k
                break
            if lineInd == 0:
                print('Reaction source not found.')
        rSqIdeal = 4./3*np.pi*R**3*k*yInf*p/Runiv/T
        print('reaction source = %g'%rS)
        etaSim = rS/rSqIdeal
        print('Simulation effectivness factor is %g, relative error = %g'%(etaSim,(etaAnal-etaSim)/etaAnal))
        
        # -- load concetration profile to compare with analytical solution
        timeLst = os.listdir('./')
        timeLst = sorted([time for time in timeLst if isFloat(time)])
        for timeInd in range(len(timeLst)):
            with open('postProcessing/graphCellFace(start=(000),end=(100),fields=(CO))/%s/line.xy'%timeLst[timeInd],'r') as fl:
                lines = fl.readlines()
            simData = np.zeros((len(lines)-1,2))
            for lineInd in range(len(lines)-1):
                simData[lineInd,:] = np.array([num for num in np.array(lines[lineInd+1].replace('\n','').split(' ')) if isFloat(num)]).astype(float)
            plt.plot(simData[:,0],simData[:,1],label='sim, time = %s'%timeLst[timeInd])

        # -- analytical solution
        cinf = yInf
        rAnal = np.linspace(1e-6,R,100)
        cAnal = cinf * R/rAnal * np.sinh(rAnal*(k/DEff)**0.5)/np.sinh(thiele)
        with open('anal.csv','w')as fl:
            fl.writelines('r,yinf\n')
            for i in range(len(rAnal)):
                fl.writelines('%g,%g\n'%(rAnal[i],cAnal[i]))
        plt.plot(rAnal,cAnal,label='analytical')
        plt.legend()
        plt.show()
        # print('Norm (analytical-simulation)')
        
        os.chdir('../../')