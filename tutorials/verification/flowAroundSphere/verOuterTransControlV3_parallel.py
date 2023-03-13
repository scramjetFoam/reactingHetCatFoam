# verOuterTransControlV3.py
# -- Script for verification simulation of conjugated mass transport
# -- V3_parallel: uses parallel snappHexMesh instead of protoMeshes

# -- imports
import numpy as np
import os
import shutil as sh
import matplotlib.pyplot as plt
import sys
from auxiliarFuncsV3 import *

# -- set solver to be used
solver = 'reactingHetCatSimpleFoam'

# -- script arguments logic
args = sys.argv

# -- directory naming
baseCaseDir = 'baseCase'
outFolder = 'ZZ_cases'
ZZZ_path = 'ZZZ_res'
ZZZ_file = 'flow.csv'
ZZZ_filepath = ZZZ_path+'/'+ZZZ_file
if not os.path.isdir('%s'%outFolder): os.mkdir('%s'%outFolder)

# -- set number of enthalpy correctors
numOfTCorr = 0

# -- case parameters
yInf = 0.01         # molar fraction farfar from sphere
p = 101325          # presure
Runiv = 8.314       # universal gas constant
sHr = -283e3        # standard reaction enthalpy (physical 283e3)	
T = 500             # temperature
DFreeZ = 1e-5       # Diffusivity in fluid
kappaEff = 1        # mass transfer coefficient
eps = 0.5           # material porosity
nu = 5.58622522e-05 # kinematic viscosity
EA = 90e3           # activation energy

# -- geometry
R = 0.01            # sphere radius
length1 = 15*R      # inlet <-> sphere centre
length2 = 45*R      # sphere centre <-> outlet
width = 15*R        # top|bottom wall <-> sphere centre

# -- list parameters [ORIGINAL]
ReLst = [10, 40, 80, 160]            # Reynolds number
invLst = [round(Re*nu/2/R,4) for Re in ReLst] # inlet velocity
thieleLst = [2,6]       # Thiele modulus
cellSizeLst = [0.4*R]  # FV cell Size
# cellSizeLst = [0.4*R, 0.2*R, 0.1*R]
tortLst = [1]           # tortuosity

# -- create cases for:
cases = [(inv,cellSize,tort,thiele) for inv in invLst for cellSize in cellSizeLst for tort in tortLst for thiele in thieleLst]
for case in cases:
    # -- parameters
    inv,cellSize,tort,thiele = case
    Re = ReLst[invLst.index(inv)]
    DFree = DFreeZ
    DEff = DFree*eps/tort
    k0Art = DEff*(thiele/R)**2 
    k0 = k0Art*np.exp(EA/(Runiv*T))
    # thiele tort re cS
    caseName = 'flow_phi_%g_tort_%g_Re_%g_cS_%g'%(thiele,tort,Re,cellSize)
    caseDir = '%s/%s/'%(outFolder,caseName)
    

    print('Preparing case %s'%caseName)
    # -- check that caseDir is clean
    if os.path.isdir(caseDir): sh.rmtree(caseDir)   # ensure the caseDir is clear
    # -- copy files
    sh.copytree(baseCaseDir,caseDir)
    # -- write parameters
    changeInCaseFolders(caseDir,'system/blockMeshDict',['length1', 'length2', 'width','nDiscX','nDiscYZ'],[str(length1),str(length2),str(width),str(int((length1+length2)/cellSize)),str(int(2*width/cellSize))])
    changeInCaseFolders(caseDir,'system/snappyHexMeshDict',['spR'],[str(R)])
    changeInCaseFolders(caseDir,'0.org/T',['isoT'],[str(T)])
    changeInCaseFolders(caseDir,'0.org/CO',['yCOSet'],[str(yInf)])
    changeInCaseFolders(caseDir,'0.org/U', ['inv'],[str(inv)])
    changeInCaseFolders(caseDir,'system/controlDict',['customSolver'],[solver])
    changeInCaseFolders(caseDir,'system/fvSolution',['nTCorr'],[str(numOfTCorr)])
    changeInCaseFolders(caseDir,'constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(k0),str(EA),str(sHr)])
    changeInCaseFolders(caseDir,'constant/transportProperties',['kappaEffSet','tortSet','DSet'],[str(kappaEff),str(tort),str(DFreeZ)])
    # -- run simulation
    os.chdir(caseDir)
    os.system('chmod u=rwx All*')
    os.system('python3 Allrun.py')

    os.chdir('../../')