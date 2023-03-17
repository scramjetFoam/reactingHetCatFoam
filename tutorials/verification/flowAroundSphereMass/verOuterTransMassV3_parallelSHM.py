# -- Script for verification simulation of conjugated mass transport
# -- parallel snappyHexMesh version

import sys
import numpy as np
import os
sys.path.append('../../../ZZ_pythonCtrlScripts')
from OF_caseClass import OpenFOAMCase
from auxiliarFuncsV3 import *

# -- set solver to be used
solverLst = ['reactingHetCatSimpleFoamM','scalarTransportFoamCO']
solver = solverLst[0]

# -- set number of the enthalpy corrections
numOfTCorr = 0
isothermal = (True if numOfTCorr==0 else False)  # isothermal logic

# -- script arguments logic
args = sys.argv
runSim = (True if 'runSim' in args else False)
showPlots = (True if 'showPlots' in args else False)
makeMesh = (True if 'makeMesh' in args else False)
getCsv = (True if 'getCsv' in args else False)
errMesh = (True if 'errMesh' in args else False)

# -- baseCase directory 
verTestDir = "../tutorials/verification/flowAroundSphereMass"
# -- directory naming
baseCaseDir = '../baseCaseMass_local'
outFolder = 'ZZ_cases/'
ZZZ_path = 'ZZZ_res'
ZZZ_file = 'flow.csv'
ZZZ_filepath = ZZZ_path+'/'+ZZZ_file

# -- case parameters
yInf = 0.01         # molar fraction of CO farfar from sphere
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
ReLst = [10, 40, 80, 160]                       # Reynolds number
invLst = [round(Re*nu/2/R,4) for Re in ReLst]   # inlet velocity
thieleLst = [2, 6]                              # Thiele modulus
cellSizeLst = [0.4*R]                           # FV cell Size
# tortLst = [0.5, 1.0, 2.5, 5.0]                  # tortuosity         
tortLst = [0.5, 2.5, 5.0]  
endTime = 500
numOfCCorr = 2
nProc = 12

# -- TESTING
# ReLst = [10]
# tortLst = [0.5]  
# thieleLst = [2]
# cellSizeLst = [0.4*R]
# ----------

# -- chemical species
specieNames = np.array(["CO", "prod", "N2"])
species = ' '.join(specieNames) 
molMass = np.array([28e-3, 28e-3, 28e-3])
sigmaVs = np.array([18.0, 18.0, 18.0])
nuVec = np.array([-1,1,0])
alphaVec = np.array([1, 0, 0])
yIn = np.array([yInf, 0, 1-yInf])
MgIn = np.sum(yIn * molMass)
wIn = yIn * molMass / MgIn

# -- different cellSize meshes baseCase list
meshesCaseLst = []

# -- numpy array for results:
if errMesh:
    emdNp = np.zeros((2,len(cellSizeLst)))

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
    
    # -- mesh parameters
    caseHere = OpenFOAMCase()
    caseHere.loadOFCaseFromBaseCase(baseCaseDir)
    caseHere.changeOFCaseDir(caseDir)
    caseHere.copyBaseCase()

    # -- creation and updates of the boundary conditions
    for chSpI in range(len(specieNames)):
        name = specieNames[chSpI]
        caseHere.runCommands(['cp 0.org/bsChemSp 0.org/%sMass' % name])
        caseHere.replace([["0.org/%sMass"% (name), ['wChSpSetInit', 'wChSpSet', 'nameSet'], ['%.5g'%(wIn[chSpI]), '%.5g'%(wIn[chSpI]), str(name)]]])
        caseHere.addToDictionary([['0.org/%sMass' % name, '\n\tsides\n\t{\n\t\ttype zeroGradient;\n\t}\n\n', 'boundaryField']])
    caseHere.runCommands(['rm 0.org/bsChemSp'])

    # -- transport properties updates
    caseHere.addToDictionary([['constant/transportProperties', 'species (%s);\n' % species, '']])
    for nameInd in range(len(specieNames)):
        name = specieNames[nameInd]
        caseHere.addToDictionary( 
            [
                ['constant/transportProperties', '%s\n{\n}\n' % name, ''],
                ['constant/transportProperties', 'D  D\t[0 2 -1 0 0 0 0] DSet;\n', name],
                ['constant/transportProperties', 'sigmaV\t%g;\n' % sigmaVs[nameInd], name],
                ['constant/transportProperties', 'molM\t%g;\n' % molMass[nameInd], name],
                ['constant/transportProperties', 'nuVec\t(0 %g);\n' % nuVec[nameInd], name],
                ['constant/transportProperties', 'alphaVec\t(0 %g);\n' % alphaVec[nameInd], name],
            ]
        )

    # -- replace in blockMeshDict    
    caseHere.replace(
        [
            [
                "system/blockMeshDict",
                ['length1', 'length2', 'width','nDiscX','nDiscYZ'],
                [str(length1),str(length2),str(width),str(int((length1+length2)/cellSize)),str(int(2*width/cellSize))]
            ], 
            [
                "system/snappyHexMeshDict", 
                ['spR'], 
                [str(R)]
            ]
        ]
    )
    caseHere.setParameters([
        ['system/controlDict', 'endTime', str(endTime), ''],
        # ['system/decomposeParDict', 'numberOfSubdomains', str(nProc), ''],
        ['system/fvSolution', 'nConcCorrectors', str(numOfCCorr), ''],
        ['system/fvSolution', 'nTempCorrectors', str(numOfTCorr), ''],
        ['system/fvSolution', 'p', str(0.18), 'fields'],
        # [, 'D', str(DFreeZ), 'CO'],
    ])
    caseHere.replace([
        ['0.org/T',['isoT'],[str(T)]],
        ['0.org/U', ['inv'],[str(inv)]],
        ['system/controlDict',['customSolver'],[solver]],
        ['system/fvSolution',['customSolver'],['customSolver|%s' %species.replace(' ','Mass|')]],
        ['constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(k0),str(EA),str(sHr)]],
        ['constant/transportProperties',['kappaEffSet','tortSet','DSet'],[str(kappaEff),str(tort),str(DFreeZ)]]
    ])
    caseHere.setParameters([
        ['system/decomposeParDict', 'numberOfSubdomains', str(nProc), '']
    ])

    if makeMesh:
        caseHere.runCommands([
            'rm -rf 0',
            'mkdir 0', 
            'cp -rf 0.org/* 0', 
            'mkdir dynamicCode',
            'blockMesh > log.blockMesh',
            # 'paraFoam -touch',
            'decomposePar > log1.decomposePar',
            'foamJob -parallel -screen snappyHexMesh -overwrite > log.snappyHexMesh',
            'reconstructParMesh -constant -time 0 > log1.regconstructParMesh',
            # 'rm -rf */cellLevel',
            # 'rm -rf */pointLevel',
        ])
    if runSim:
        caseHere.runCommands([
            # 'rm -rf 0',
            # 'mkdir 0',
            # 'cp -rf 0.org/* 0',
            'decomposePar -force > log2.decomposePar',
            'foamJob -parallel -screen renumberMesh -overwrite > log.renumberMesh', 
            'foamJob -parallel -screen %s > log.%s' % (solver, solver),
        ])
        if not os.path.isdir('processor0/100'): caseHere.runCommands(['foamJob -parallel -screen %s > log2.%s' % (solver, solver)])
        caseHere.runCommands([
            'reconstructPar -latestTime > log2.reconstructPar',
            'intSrcSphereM > log.intSrcSphereM',
            'postProcess -func integrace -latestTime > log.integrace',
            'rm -rf processor*'
        ])
    
