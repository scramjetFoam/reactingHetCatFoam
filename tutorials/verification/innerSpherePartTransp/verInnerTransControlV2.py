# verInnerTransControlV2.py
# -- Script for verification cases with heat & mass transfer in a spherical particle
# -- NOTE: SCRIPT USES ARGUMENTS:
#       -- "makeMesh"
#       -- "runSim"
#       -- "showPlots"
# -- TODO: 
#       (1) edit [if isothermal + else] in the loop
#       (2) for non-isoT: edit plot making
#       (3) add csv export
#       (4) for non-isoT: incorporate theoretical solution


# -- imports
import numpy as np
import os
import re
import shutil as sh
import matplotlib.pyplot as plt
import sys
from auxiliarFuncs import *
# from shootChandra import *


# -- set solver to be used
solverLst = ['reactingHetCatSimpleFoam']
solver = solverLst[0]

# -- set number of the enthalpy corrections
numOfTCorr = 1
isothermal = (True if numOfTCorr == 0 else False)  # isothermal logic

# -- script arguments logic
args = sys.argv
runSim = (True if 'runSim' in args else False)
showPlots = (True if 'showPlots' in args else False)
makeMesh = (True if 'makeMesh' in args else False)
getCsv = (True if 'getCsv' in args else False)

# -- directory naming
baseCaseDir = 'baseCase'
outFolder = 'ZZ_cases'
if isothermal: outFolder += 'isoT'
ZZZ_path = 'ZZZ_res'
if isothermal: ZZZ_path += 'isoT'
# ZZZ_file = 'mult_steadySt'
# ZZZ_filepath = ZZZ_path+'/'+ZZZ_file

# -- set case parameters
yInf = 0.02      # molar fraction farfar from sphere
p = 101325      # presure	
Runiv = 8.314   # universal gas constant
R = 0.1           # sphere radius

# -- set geometry
domainSize = 1.1*R
tort = 2                # tortuosity
kappaEff = 2            # mass transfer coefficient
DFreeZ = 1e-5           # set diffusivity in fluid

# -- list parameters
thieleLst = [0.5,0.75,1.,2,4]   # Thiele modulus
TLst = [300]                    # temperature
gammaLst = [20]                 # gamma - Arrhenius number
betaLst = [0.6]                 # beta - Prater number
cellSizeLst = [0.5*R]           # NOTE: The mesh will be much more refined inside the sphere (5 5)

# == ARCHIVED SETTINGS: 
# -- 12/11/2022 khyrm@multipede: [5 cases for multSteadySt with (5 5) refinement]
thieleLst = [0.5,0.75,1.,2,4]
TLst = [300]
gammaLst = [20]
betaLst = [0.6]
cellSizeLst = [0.35*R]

# -- 14/11/2022 khyrm@multipede: [2 cases for multSteadySt with (8 8) refinement]
baseCaseDir += '_88'
outFolder += '_88'
# thieleLst = [0.2,4]
thieleLst = [0.2]
TLst = [300]
gammaLst = [20]
betaLst = [0.6]
cellSizeLst = [0.35*R]


# -- prepare prototype mesh for each cellSize
if makeMesh:
    for cellSize in cellSizeLst:
        if not os.path.isdir('%s'%outFolder): os.mkdir('%s'%outFolder)
        meshDir = '%s/protoMesh/%g'%(outFolder,cellSize)
        print('Preparing mesh %s',meshDir)
        # -- check that meshDir is clean
        if os.path.isdir(meshDir): sh.rmtree(meshDir)
        # -- copy files
        sh.copytree(baseCaseDir,meshDir)
        changeInCaseFolders(meshDir,'system/blockMeshDict',['dSN', 'nDisc'],[str(domainSize),str(int(domainSize/cellSize*2))])
        changeInCaseFolders(meshDir,'system/snappyHexMeshDict',['spR'],[str(R)])
        changeInCaseFolders(meshDir,'system/snappyHexMeshDictIntraTrans',['spR'],[str(R)])
        os.chdir(meshDir)
        os.system('chmod u=rwx All*') # Just to make sure.
        os.system('./AllmeshIntraSphere')
        os.chdir('../../../')

# -- numpy array for results
if isothermal:
    resNp = np.zeros((len(cellSizeLst)))
    resNp2 = np.zeros((len(cellSizeLst)))
else:
    resNp = np.zeros((2,len(thieleLst)+1))

# -- create cases for
cases = [(T,thiele,cellSize,beta,gamma) for T in TLst for thiele in thieleLst for cellSize in cellSizeLst for beta in betaLst for gamma in gammaLst ]
for case in cases:
    # -- parameters
    T,thiele,cellSize,beta,gamma = case
    EA = gamma*Runiv*T
    DFree = DFreeZ
    DEff = DFree/tort*0.5
    sHr = -beta/yInf/p*Runiv*T/DEff*kappaEff*T
    k0 = (thiele/R)**2 * DEff/(np.exp(-gamma))

    caseName = 'intraTrans_phi_%g_beta_%g_cellSize_%g'%(thiele,beta,cellSize)
    caseDir = '%s/%s/'%(outFolder,caseName)
    meshDir = '%s/protoMesh/%g'%(outFolder,cellSize)

    if runSim:
        print('Preparing case %s'%caseName)
        # -- check that caseDir is clean
        if os.path.isdir(caseDir): sh.rmtree(caseDir)
        # -- copy files
        sh.copytree(meshDir,caseDir)
        # -- change the case files
        changeInCaseFolders(caseDir,'0.org/T',['isoT'],[str(T)])
        changeInCaseFolders(caseDir,'0.org/CO',['yCOSet'],[str(yInf)])
        changeInCaseFolders(caseDir,'system/controlDict',['customSolver'],[solver])
        changeInCaseFolders(caseDir,'system/fvSolution',['nTCorr'],[str(numOfTCorr)])
        changeInCaseFolders(caseDir,'constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(k0),str(EA),str(sHr)])
        changeInCaseFolders(caseDir,'constant/transportProperties',['kappaEffSet','tortSet','DSet'],[str(kappaEff),str(tort),str(DFreeZ)])
        # -- run simulation
        os.chdir(caseDir)
        os.system('chmod u=rwx AllrunIntraSphere')
        os.system('./AllrunIntraSphere')
    else:
        os.chdir(caseDir)
    
    k = k0*np.exp(-gamma)       # ???
    k0Art = k0*np.exp(-gamma)   # ???
    rSqIdeal = 4/3*np.pi*R**3*k0Art*yInf*p/Runiv/T     # ideal reaction source
    rS = read_real_source()                            # simulation reaction source

    print('Case with thiele = %g'%thiele)
    
    if isothermal: 
        # -- analytical effectivness factor
        etaAnal = 3./(thiele**2)*(thiele*1./np.tanh(thiele)-1)
        print('Reaction source = %g, simulation effectivness factor is %g, relative error = %g'%(rS,etaSim,(etaAnal-etaSim)/etaAnal))
        print('\nThiele modulus = %g\nAnalytical effectivness factor is %g'%(thiele,etaAnal))   
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
        rAnal = simData[:,0]
        yAnal = yInf * R/rAnal * np.sinh(rAnal*(k/DEff)**0.5)/np.sinh(thiele)
        resNp[cellSizeInd] = abs(etaAnal-etaSim)
        resNp2[cellSizeInd] = np.linalg.norm(yAnal-simData[:,1])
        with open('anal.csv','w') as fl:
            fl.writelines('r,yinf\n')
            fl.writelines(['%g,%g\n'%(rAnal[i],yAnal[i]) for i in range(length(rAnal))])
        plt.plot(rAnal,yAnal,label='analytical')
        plt.legend()
        plt.title('Concentration profile comparison for k0=%d'%k0)
        plt.savefig('resHere.png')
        if showPlots:
            plt.show()
    else:
        etaSim = rS/rSqIdeal  # simulation effectivness factor
        print('beta',beta,'thiele',thiele,'etaSim',etaSim,'R',R,'k0Art',k0Art,'Deff',DEff)
        resNp[:,thieleLst.index(thiele)] = np.array([thiele,etaSim])
    os.chdir('../../')


# -- create plots
if showPlots:
    if isothermal:
        # -- mesh dependence plot
        plt.plot(cellSizeLst,resNp,label='eta diff (my)')
        # plt.plot(cellSizeLst,resNp[1],label='eta diff (STF)')
        plt.plot(cellSizeLst,resNp2,label='whole sol diff (my)')
        # plt.plot(cellSizeLst,resNp2[1],label='whole sol diff (STF)')
        plt.plot(cellSizeLst,cellSizeLst,label='slope = 1')
        plt.plot(cellSizeLst,np.array(cellSizeLst)**2,label='slope = 2')
        title = 'Dependence of the error on the mesh.'
        fileName = 'error_mesh.png'
        plt.title(title)
        plt.yscale('log')
        plt.xscale('log')
        plt.legend()
        plt.show()
    else:
        # -- eta-thiele diagram
        with open('baseCase/analRes06.csv' ,'r') as fl:
            lines = fl.readlines()
        analRes = np.zeros((len(lines)-1,2))
        for i in range(1,len(lines)):
            analRes[i-1,:] = np.array(lines[i].split(',')) 
        plt.plot(analRes[:,0], analRes[:,1],'r',label='anal. res')
        plt.plot(resNp[0,:],resNp[1,:],'x',label='sim res.')
        title = 'Multiple steady states.'
        fileName = 'mult_steadySt.png'
        plt.title(title)
        plt.yscale('log')
        plt.xscale('log')
        plt.legend()
        # plt.savefig('%s/%s'%(ZZZ_path,fileName))
        plt.show()