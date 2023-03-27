# -- verInnerTransControlMassV1_parSHM.py
# -- Script for verification cases with heat & mass transfer in a spherical particle

# -- imports
import numpy as np
import os
import shutil as sh
import matplotlib.pyplot as plt
import sys
sys.path.append('../../../ZZ_pythonCtrlScripts')
from OF_caseClass import OpenFOAMCase
from auxiliarFuncsMassV1 import *

# -- obtained from shootChandraVM.py
# nonisoT_etaAnal = 42.094    # phi = 0.5
nonisoT_etaAnal =  7.516    # phi = 4.0

# -- set solver to be used
solverLst = ['reactingHetCatSimpleFoamM']
solver = solverLst[0]

# -- simulation parameters
endTime = 500
nProc = 8
numOfCCorr = 1
# -- set number of the enthalpy corrections
numOfTCorr = 1
isothermal = (True if numOfTCorr == 0 else False)  # isothermal logic

# -- directory naming
baseDir = '../baseCaseMass_local_intraTransp'
outFolder = 'ZZ_cases'
if isothermal: outFolder += '_isoT'
ZZZ_path = 'ZZZ_res'
if isothermal: ZZZ_path += '_isoT'
# ZZZ_file = 'mult_steadySt'
# ZZZ_filepath = ZZZ_path+'/'+ZZZ_file
if not os.path.exists(ZZZ_path): os.mkdir(ZZZ_path)
if not os.path.exists(outFolder): os.mkdir(outFolder)

# -- set case parameters
yInf = 0.02      # molar fraction farfar from sphere
p = 101325      # presure	
Runiv = 8.314   # universal gas constant
R = 0.1           # sphere radius
T0 = False     # used for multiple steady state (800), comment out or set to False if not used
# T0 = 800

# -- set geometry
domainSize = 1.1*R
tort = 2                # tortuosity
kappaEff = 2            # mass transfer coefficient
DFreeZ = 1e-5           # set diffusivity in fluid

# -- list parameters [ORIGINAL]
thieleLst = [0.2, 0.4, 0.5, 0.75, 1.0, 2.0, 4.0]
# thieleLst = [0.4, 0.5, 0.75, 1.0, 2.0, 4.0]
thieleLst = [0.2]
TLst = [300]
gammaLst = [20]
betaLst = [0.6]
# cellSizeLst = [0.4*R]  # NOTE: The mesh will be much more refined inside the sphere: (5 5)

# -- thiele lst
# thieleLst = [0.5]
cellSizeLst = [0.2*R]

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

# -- tests
# thieleLst = [0.5]
# thieleLst = [4.0]
# cellSizeLst = [1.6*R, 0.8*R, 0.4*R, 0.2*R, 0.1*R]
# cellSizeLst = [0.4*R, 0.2*R, 0.1*R]
# cellSizeLst = [1.6*R, 0.8*R]


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
    if not T0: T0 = T

    caseName = 'intraTrans_phi_%g_beta_%g_cellSize_%g_T_%g'%(thiele,beta,cellSize,T)
    caseDir = '%s/%s/'%(outFolder,caseName)

    caseHere = OpenFOAMCase()
    caseHere.loadOFCaseFromBaseCase(baseDir)
    caseHere.changeOFCaseDir(caseDir)
    caseHere.copyBaseCase()
    # ---
    # -- create and update boundary conditions
    for chSpI in range(len(specieNames)):
        name = specieNames[chSpI]
        caseHere.runCommands(['cp 0.org/bsChemSp 0.org/%sMass' % name])
        caseHere.addToDictionary([['0.org/%sMass' % name, '\n\t"(sphere|infinity)"\n\t{\n\t\ttype\t\t\tfixedValue;\n\t\tvalue\t\t\tuniform wChSpSet;\n\t}\n', 'boundaryField']])
        caseHere.replace([["0.org/%sMass"% (name), ['wChSpSetInit','wChSpSet','nameSet'], ['%.5g'%(wIn[chSpI]), '%.5g'%(wIn[chSpI]), str(name)]]])
    caseHere.addToDictionary([['0.org/U', '\n\t"(sphere|infinity)"\n\t{\n\t\ttype\t\t\tfixedValue;\n\t\tvalue\t\t\tuniform (0 0 0.0);\n\t}\n', 'boundaryField']])
    caseHere.addToDictionary([['0.org/T', '\n\t"(sphere|infinity)"\n\t{\n\t\ttype\t\t\tfixedValue;\n\t\tvalue\t\t\tuniform boundT;\n\t}\n', 'boundaryField']])
    caseHere.addToDictionary([['0.org/p', '\n\tinfinity\n\t{\n\t\ttype\t\t\ttotalPressure;\n\t\tp0\t\t\t\tuniform 101325;\n\t\tvalue\t\t\tuniform 101325;\n\n\t\tU\t\t\t\tU;\n\t\tphi\t\t\t\tphi;\n\t}\n', 'boundaryField']])            
    
    caseHere.runCommands(['rm 0.org/bsChemSp'])
    
    # -- update transport properties
    caseHere.addToDictionary([['constant/transportProperties', 'species (%s);\n' % species, '']])

    for nameInd in range(len(specieNames)):
        name = specieNames[nameInd]
        caseHere.addToDictionary([
            ['constant/transportProperties', '%s\n{\n}\n' % name, ''],
            ['constant/transportProperties', 'D  D\t[0 2 -1 0 0 0 0] DSet;\n', name ],
            ['constant/transportProperties', 'sigmaV\t%g;\n' % sigmaVs[nameInd], name ],
            ['constant/transportProperties', 'molM\t%g;\n' % molMass[nameInd], name ],
            ['constant/transportProperties', 'nuVec\t(0 %g);\n' % nuVec[nameInd], name ],
            ['constant/transportProperties', 'alphaVec\t(0 %g);\n' % alphaVec[nameInd], name ],
        ])

    # -- update blockMeshDict
    caseHere.replace([
        ['system/blockMeshDict', ['length1', 'length2', 'width', 'width'], [str(domainSize), str(domainSize), str(domainSize), str(domainSize)]],
        ['system/blockMeshDict', ['nDiscX', 'nDiscYZ'], [str(int(domainSize/cellSize*2)), str(int(domainSize/cellSize*2))]]
    ])
    # -- update others
    caseHere.replace([
        ["system/decomposeParDict",['np'],[str(nProc)]],
        ["system/snappyHexMeshDict",['spR'],[str(R)]]
    ])
    caseHere.replace([
        ['0.org/U',['inv'],[str(0)]]
    ])
    # -- setup verification case from the base case
    # -- creation and updates of the boundary conditions
    caseHere.setParameters([
        ['system/controlDict', 'endTime', str(endTime), ''],
        ['system/decomposeParDict', 'numberOfSubdomains', str(nProc), ''],
        ['system/fvSolution', 'nConcCorrectors', str(numOfCCorr), ''],
        ['system/fvSolution', 'nTempCorrectors', str(numOfTCorr), ''],
        ['system/fvSolution', 'p', str(0.18), 'fields'],
        ['constant/transportProperties', 'D', str(DFreeZ), 'CO']
    ])
    # ---
    caseHere.replace([
        ['0.org/T',['initT','boundT'],[str(T0),str(T)]],
        # ['0.org/CO',['yCOSet'],[str(yInf)]],
        ['system/controlDict',['customSolver'],[solver]],
        ['system/fvSolution',['customSolver'],['customSolver|%s' %species.replace(' ','Mass|')]],
        ['constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(k0),str(EA),str(sHr)]],
        ['constant/transportProperties',['kappaEffSet','tortSet','DSet'],[str(kappaEff),str(tort),str(DFreeZ)]]
    ])
    caseHere.runCommands([
        'rm -rf 0',
        'mkdir 0', 
        'cp -rf 0.org/* 0', 
        'mkdir dynamicCode',
        'blockMesh > log.blockMesh',
        'paraFoam -touch',
        'decomposePar > log1.decomposePar',
        'foamJob -parallel -screen snappyHexMesh -overwrite > log.snappyHexMesh',
        'reconstructParMesh -time 0 > log1.reconstructPar1',
        'reconstructParMesh -constant > log1.reconstructPar2',
        'rm -rf processor*',
        'topoSet > log.topoSet',
        'rm -rf */cellLevel',
        'rm -rf */pointLevel',
    ])
    caseHere.runCommands([
        'decomposePar -force > log2.decomposePar',
        'foamJob -parallel -screen renumberMesh -overwrite > log.renumberMesh', 
        'foamJob -parallel -screen %s > log1.%s' % (solver, solver),
    ])

    # -- rerun if the first run fails
    rerun = False
    with open('%slog1.%s'%(caseDir,solver), 'r') as f:
        if (f.readlines()[-3][:-1] != 'End'): rerun = True
    if rerun: caseHere.runCommands(['foamJob -parallel -screen %s > log2.%s' % (solver, solver)])

    caseHere.runCommands([
        'reconstructPar -latestTime > log2.reconstructPar',
        'intSrcSphereM > log.intSrcSphereM',
        # 'postProcess -func integrace -latestTime > log.integrace',
        "postProcess -func 'graphCellFace(start = (0 0 0), end = (1 0 0), fields=(COMass))' > log.postProcess",
        'rm -rf processor*'
    ])

    os.chdir(caseDir)
    k = k0*np.exp(-gamma)       # ???
    k0Art = k0*np.exp(-gamma)   # ???
    rSqIdeal = 4/3*np.pi*R**3*k0Art*yInf*p/Runiv/T     # ideal reaction source
    rS = read_real_source()                            # simulation reaction source
    etaSim = rS/rSqIdeal                               # simulation effectivness factor

    print('Case with thiele = %g'%thiele)
    
    # nonisothermal
    print('beta',beta,'thiele',thiele,'etaSim',etaSim,'R',R,'k0Art',k0Art,'Deff',DEff)

    os.chdir('../../')