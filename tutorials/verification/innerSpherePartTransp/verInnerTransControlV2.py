# verInnerTransControlV2.py
# -- Script for verification cases with heat & mass transfer in a spherical particle
# -- NOTE: SCRIPT USES ARGUMENTS:
#       -- "makeMesh"
#       -- "runSim"
#       -- "showPlots"

# -- imports
import numpy as np
import os
import re
import shutil as sh
import matplotlib.pyplot as plt
import sys
from auxiliarFuncs import *


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
# ZZZ_file = 'res.csv'
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
# 12/11/2022
cellSizeLst = [0.2*R]

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

# # -- create cases for
# cases = [(T,thiele,cellSize,gamma,beta) for T in TLst for thiele in thieleLst for cellSize in cellSizeLst for gamma in gammaLst for beta in betaLst]
# for case in cases:
#     # -- parameters
#     T,thiele,cellSize = case
#     EA = gamma*Runiv*T
#     DEff = DFreeZ/tort*0.5
#     k0 = (thiele/R)**2 * DEff/(np.exp(-gamma))
#     sHr = -beta/yInf/p*Runiv*T/DEff*kappaEff*T

#     caseName = 'intraTrans_phi_%g_beta_%g_cellSize_%g'%(thiele,beta,cellSize)
#     caseDir = '%s/%s/'%(outFolder,caseName)

#     if runSim:
#         print('Preparing case')