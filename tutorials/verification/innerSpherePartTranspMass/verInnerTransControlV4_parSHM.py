# -- verInnerTransControlV4_parSHM.py
# -- Script for verification cases with heat & mass transfer in a spherical particle.
# -- NOTE: V4_parSHM:
#   -- uses parallel snappHexMesh instead of protoMeshes
#   -- only for meshing + running the simulation
#   -- for evaulation, use the regular V4 script
# -- NOTE SCRIPT ARGUMENTS:
#   -- "onlyMesh"

# -- imports
import numpy as np
import os
import shutil as sh
import matplotlib.pyplot as plt
import sys
sys.path.append('../../../ZZ_pythonCtrlScripts')
from OF_caseClass import OpenFOAMCase
from auxiliarFuncs import *

# -- obtained from shootChandraVM.py
nProcs = 60  # number of processors
# nProcs = 2

# -- set solver to be used
solverLst = ['reactingHetCatSimpleFoam']
solver = solverLst[0]

# -- set number of the enthalpy corrections
numOfTCorr = 1
isothermal = (True if numOfTCorr == 0 else False)  # isothermal logic

# -- script arguments logic
args = sys.argv
onlyMesh = (False if 'onlyMesh' not in args else True)

# -- directory naming
baseDir = 'baseCase'
outFolder = 'ZZ_cases'
if isothermal: outFolder += '_isoT'
ZZZ_path = 'ZZZ_res'
if isothermal: ZZZ_path += '_isoT'
if not os.path.exists(ZZZ_path): os.mkdir(ZZZ_path)
if not os.path.exists(outFolder): os.mkdir(outFolder)

# -- set case parameters
yInf = 0.02      # molar fraction farfar from sphere
p = 101325      # presure	
Runiv = 8.314   # universal gas constant
R = 0.1           # sphere radius
T0 = False     # used for multiple steady state (800), comment out or set to False if not used
T0 = 800

# -- set geometry
domainSize = 1.1*R
tort = 2                # tortuosity
kappaEff = 2            # mass transfer coefficient
DFreeZ = 1e-5           # set diffusivity in fluid

# -- list parameters [ORIGINAL]
thieleLst = [0.2, 0.4, 0.5, 0.75, 1, 2, 4]
TLst = [300]
gammaLst = [20]
betaLst = [0.6]
cellSizeLst = [0.4*R]  # NOTE: The mesh will be much more refined inside the sphere: (5 5)

# -- chemical species
specieNames = np.array(["CO", "prod", "N2"])
species = ' '.join(specieNames) 
molMass = np.array([ 28e-3, 28e-3, 28e-3])
sigmaVs = np.array([ 18.0, 18.0, 18.0])
nuVec = np.array([-1,1,0])
alphaVec = np.array([1, 0, 0])
yIn = np.array([yInf, 0, 1-yInf])
MgIn = np.sum(yIn * molMass)
wIn = yIn * molMass / MgIn


# -- mesh tests
# thieleLst = [0.5, 4.0]
# thieleLst = [0.5]
thieleLst = [4.0]
# cellSizeLst = [1.6*R, 0.8*R, 0.4*R, 0.2*R, 0.1*R]
# cellSizeLst = [1.6*R, 0.8*R]
cellSizeLst = [0.4*R, 0.2*R, 0.1*R]

# T0 = 800 # for mult. st. states, set to False to ignore
T0 = False

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
    if not T0: T0 = T

    caseName = 'intraTrans_phi_%g_beta_%g_cellSize_%g_T_%g'%(thiele,beta,cellSize,T)
    caseDir = '%s/%s/'%(outFolder,caseName)
    meshDir = '%s/protoMesh/%g'%(outFolder,cellSize)

    caseHere = OpenFOAMCase()
    caseHere.loadOFCaseFromBaseCase(baseDir)
    caseHere.changeOFCaseDir(caseDir)
    caseHere.copyBaseCase()
    # -- mesh parameters
    caseHere.replace([
        ["system/decomposeParDict",['np'],[str(nProcs)]],
        ["system/blockMeshDict",['dSN','nDisc'],[str(domainSize),str(int(domainSize/cellSize*2))]],
        ["system/snappyHexMeshDict",['spR'],[str(R)]]
    ])
    # -- case parameters
    caseHere.replace([
        ['0.org/T',['initT','boundT'],[str(T0),str(T)]],
        ['0.org/CO',['yCOSet'],[str(yInf)]],
        ['system/controlDict',['customSolver'],[solver]],
        ['system/fvSolution',['nTCorr'],[str(numOfTCorr)]],
        ['constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(k0),str(EA),str(sHr)]],
        ['constant/transportProperties',['kappaEffSet','tortSet','DSet'],[str(kappaEff),str(tort),str(DFreeZ)]]
    ])

    # -- meshing:
    caseHere.runCommands([
        'rm -rf 0',
        'mkdir 0', 
        'cp -rf 0.org/* 0', 
        'mkdir dynamicCode',
        'blockMesh > log.blockMesh',
        # 'paraFoam -touch',
        'decomposePar > log1.decomposePar',
        'foamJob -parallel -screen snappyHexMesh -overwrite > log.snappyHexMesh',
        'reconstructParMesh -time 0 > log1.reconstructPar1',
        'reconstructParMesh -constant > log1.reconstructPar2',
        'topoSet > log.topoSet',
        'rm -rf */cellLevel',
        'rm -rf */pointLevel',
    ])
    
    # -- running the simulation:
    if not onlyMesh:
        caseHere.runCommands([
            'rm -rf 0',
            'mkdir 0',
            'cp -rf 0.org/* 0',
            'decomposePar -force > log2.decomposePar',
            'foamJob -parallel renumberMesh -overwrite > log.renumberMesh', 
            'foamJob -parallel -screen %s > log1.%s' % (solver, solver),
        ])

        # -- rerun if the first run fails
        rerun = False
        with open('%slog1.%s'%(caseDir,solver), 'r') as f:
            if (f.readlines()[-3][:-1] != 'End'): rerun = True
        if rerun: caseHere.runCommands(['foamJob -parallel -screen %s > log2.%s' % (solver, solver)])

        caseHere.runCommands([
            'reconstructPar -latestTime > log2.reconstructPar',
            'intSrcSphere > log.intSrcSphere',
            # 'postProcess -func integrace -latestTime > log.integrace',
            "postProcess -func 'graphCellFace(start = (0 0 0), end = (1 0 0), fields=(CO))' > log.postProcess"
        ])
