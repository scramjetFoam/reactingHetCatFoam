# -- Script to do verification simulation of isothermal reaction inside and outside sphere

# -- imports
import numpy as np
import os 
import shutil as sh

# -- setup case
baseCaseDir = 'baseCase'
outFolder = 'ZZ_cases'

# -- case data
yInf = 0.01     # molar fraction farfar from sphere
p = 101325      # presure
Runiv = 8.314   # universal gas constant
R = 1           # sphere radius
Rinf = 5         # infinite radius 
TLst = [500,550] # temperature 
k0 = 5e7      # reaction pre-exponential factor
EA = 90e3     # reaction activation energy
tort = 3      # tortuosity
solver = 'reactingHetCatSimpleFoam' # used solver

# -- create case
for T in TLst:
    # -- create caseFolder based on baseCase
    caseName = 'yInf_%g_R_%g_Rinf_%g_T_%g_tort_%g'%(yInf,R,Rinf,T,tort)
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
    print('Changing blockMeshDict cube size to (%gx%gx%g)'%(Rinf,Rinf,Rinf))
    with open(caseDir + '/system/blockMeshDict', 'r') as file:
        data = file.readlines()

    for ind in range(len(data)):
        if data[ind].find('dSN') >= 0:
            data[ind] = data[ind].replace('dSN',str(Rinf))
            
        if data[ind].find('nDisc') >= 0:
            data[ind] = data[ind].replace('nDisc',str(Rinf*10))

    with open(caseDir + '/system/blockMeshDict', 'w') as file:
        file.writelines(data)
    
    # -- run the simulation
    os.chdir(caseDir)
    os.system('./Allrun')

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
    BiM = kL*R/DEff             # Biot number
    etaAnal = 3./(thiele**2) * (thiele*1./np.tanh(thiele)-1)/(1+(thiele*1./np.tanh(thiele)-1)/BiM)                      # analytical effectivness factor
    print('\nBiot number = %g\nThiele modulus = %g\nAnalytical effectivness factor is %g'%(BiM,thiele,etaAnal))

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
    

    os.chdir('../../')