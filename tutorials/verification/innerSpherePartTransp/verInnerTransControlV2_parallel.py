# -- verInnerTransControlV2_parallel.py
# -- Script for verification cases with heat & mass transfer in a spherical particle.
# -- V2_parallel:
# --    uses parallel snappHexMesh instead of protoMeshes
# --    only for meshing + running the simulation

# -- imports
import numpy as np
import os
import shutil as sh
import matplotlib.pyplot as plt
import sys
from auxiliarFuncs import *

# -- set solver to be used
solver = 'reactingHetCatSimpleFoam'

# -- set number of the enthalpy corrections
numOfTCorr = 1

# -- directory naming
baseCaseDir = 'baseCase'
outFolder = 'ZZ_cases'
ZZZ_path = 'ZZZ_res'
if not os.path.isdir(outFolder): os.mkdir(outFolder)
if not os.path.isdir(ZZZ_path): os.mkdir(ZZZ_path)

# -- set case parameters
yInf = 0.02         # molar fraction farfar from sphere
p = 101325          # presure	
Runiv = 8.314       # universal gas constant
R = 0.1             # sphere radius
# T0Lst = [300, 800]
T0Lst = [300]

# -- set geometry
domainSize = 1.1*R
tort = 2            # tortuosity
kappaEff = 2        # mass transfer coefficient
DFreeZ = 1e-5       # set diffusivity in fluid

# -- list parameters
TLst = [300]
gammaLst = [20]
betaLst = [0.6]
# thieleLst = [4.0, 0.5]
thieleLst = [4.0]
# cellSizeLst = [0.2*R, 0.1*R]
cellSizeLst = [1.6*R]

# -- numpy array for results
resNp = np.zeros((2,len(thieleLst)+1))

# -- create cases for
cases = [(T,thiele,cellSize,beta,gamma) for T in TLst for thiele in thieleLst for cellSize in cellSizeLst for beta in betaLst for gamma in gammaLst]
for case in cases:
    # -- parameters
    T,thiele,cellSize,beta,gamma = case
    EA = gamma*Runiv*T
    DFree = DFreeZ
    DEff = DFree/tort*0.5
    sHr = -beta/yInf/p*Runiv*T/DEff*kappaEff*T
    k0 = (thiele/R)**2 * DEff/(np.exp(-gamma))
    T0 = T0Lst[thieleLst.index(thiele)]

    caseName = 'intraTrans_phi_%g_beta_%g_cellSize_%g_T_%g'%(thiele,beta,cellSize,T)
    caseDir = '%s/%s'%(outFolder,caseName)

    print('Preparing case %s'%caseName)
    # -- check that caseDir is clean
    if os.path.isdir(caseDir): sh.rmtree(caseDir)
    # -- copy files
    sh.copytree(baseCaseDir,caseDir)
    # -- change the case files
    changeInCaseFolders(caseDir,'system/blockMeshDict',['dSN', 'nDisc'],[str(domainSize),str(int(domainSize/cellSize*2))])
    changeInCaseFolders(caseDir,'system/snappyHexMeshDict',['spR'],[str(R)])
    changeInCaseFolders(caseDir,'system/snappyHexMeshDictIntraTrans',['spR'],[str(R)])
    changeInCaseFolders(caseDir,'0.org/T',['initT','boundT'],[str(T0),str(T)])
    changeInCaseFolders(caseDir,'0.org/CO',['yCOSet'],[str(yInf)])
    changeInCaseFolders(caseDir,'system/controlDict',['customSolver'],[solver])
    changeInCaseFolders(caseDir,'system/fvSolution',['nTCorr'],[str(numOfTCorr)])
    changeInCaseFolders(caseDir,'constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(k0),str(EA),str(sHr)])
    changeInCaseFolders(caseDir,'constant/transportProperties',['kappaEffSet','tortSet','DSet'],[str(kappaEff),str(tort),str(DFreeZ)])

    # -- make mesh && run simulation        
    os.chdir(caseDir)
    os.system('chmod u=rwx All*')        
    os.system('python3 Allrun.py')
    os.chdir('../../')