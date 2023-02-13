# -- Script for verification simulation of conjugated mass transport
# -- NOTE SCRIPT ARGUMENTS:
#       -- "makeMesh": prepare mesh for each cellSize 
#           (this must be done manually if protoMesh doesn't exist yet)
#       -- "runSim": run the simulation
#       -- "showPlots": shows Plots
#       -- "getCsv": generate csv files for TeX plots 
#       -- "errMesh": evaluate error mesh dependency 
#           (for otherwise identical parameters) 
# -- NOTE V3 changelog: 
#       -- thieleLst + tortLst --> k0
#       -- ReLst + R --> invLst
#       -- flow.csv stores Thiele modulus
#       -- removed eta_anal
#       -- uses auxiliarFuncsV3
#       -- added error mesh dependency tests

# -- using OF_caseClass for the description see OF_caseClass.py
# -- imports 

from OF_caseClass import OpenFOAMCase
import sys
import numpy as np
import os
# from auxiliarFuncsV3 import *

# -- set solver to be used
solverLst = ['reactingHetCatSimpleFoamM','scalarTransportFoamCO']
solver = solverLst[0]

# -- script arguments logic
# args = sys.argv
# runSim = (True if 'runSim' in args else False)
# showPlots = (True if 'showPlots' in args else False)
# makeMesh = (True if 'makeMesh' in args else False)
# getCsv = (True if 'getCsv' in args else False)
# errMesh = (True if 'errMesh' in args else False)

# -- baseCase directory 
# verTestDir = "../tutorials/verification/flowAroundSphereMass"

# -- directory naming
baseCaseDir = '../tutorials/tested/baseCaseMassRash'
outFolder = '../ZZ_cases/testRashV1'

# -- inlet parameters
specieNames = np.array(["ethylene", "O2", "HCl", "N2", "C2H4Cl2", "H2O"])
nIn = np.array([ 6.521E-03, 2.588E-03, 1.281E-02, 1.947E-02, 5.175E-03, 5.175E-03 ])
yIn = nIn / np.sum(nIn)
for yInInd in range(len(specieNames)):
    print( "Inlet molar fraction of the specie %s is %g" % (specieNames[yInInd], yIn[yInInd]))

MolMass = np.array([ 28.05e-3, 32.0e-3, 36.46e-3, 28.01e-3, 98.96e-3, 18.01e-3])
MgIn = np.sum( yIn * MolMass )
print( "Inlet molar mass is %g." %MgIn ) 

wIn = yIn * MolMass / MgIn
for yInInd in range(len(specieNames)):
    print( "Inlet mass fraction of the specie %s is %g" % (specieNames[yInInd], wIn[yInInd]))

TIn = 200 + 273     # inlet temperature in Kelvins
UIn = 0.65          # inlet velocity
pOut = 5.1e5          # outlet pressure

# -- reaction parameters
sHr = -245274       # standard reaction enthalpy (physical 283e3)	
A1 = 2.78e5         # pre-exponential factor 1
A2 = 2e8            # pre-exponential factor 2
EA1 = 4.53e4        # activation energy 1
EA2 = 5.19e4        # activation energy 2
KEq = 0.2
cCuCl = 0.012

Runiv = 8.314       # universal gas constant

DFreeZ = 1e-5       # Diffusivity in fluid
kappaEff = 1        # mass transfer coefficient
eps = 0.5           # material porosity
EA = 90e3           # activation energy

# -- geometry
R = 0.01            # sphere radius
length1 = 15*R      # inlet <-> sphere centre
length2 = 45*R      # sphere centre <-> outlet
width = 15*R        # top|bottom wall <-> sphere centre

# -- create new folder from basecase folder
case = OpenFOAMCase()
case.loadOFCaseFromBaseCase(baseCaseDir)
case.changeOFCaseDir(outFolder)
case.copyBaseCase()

# -- set parameters:
# 1) inlet parameters
for nameInd in range(len(specieNames)):
    name = specieNames[nameInd]
    case.runCommands( [ 'cp 0.org/bsChemSp 0.org/%sMass' % name] )
    case.replace( [ [ "0.org/%sMass" % ( name ), [ 'wChSpSet', 'nameSet' ], [ '%.5g' % (wIn[nameInd]), str(name) ] ] ] )
case.setParameters( [ 
                        [ '0.org/T', 'internalField', 'uniform %.5g' % TIn, '' ],
                        [ '0.org/U', 'internalField', 'uniform (%.5g 0 0)' % UIn, ''],
                        [ '0.org/p', 'internalField', 'uniform %.5g' % pOut, ''],
                    ] )

# 2) transportProperties parameters
for nameInd in range(len(specieNames)):
    name = specieNames[nameInd]
    case.setParameters( [ [ 'constant/transportProperties', 'molM', '%.5g' % MolMass[nameInd], name ] ])

# 3) reactionProperties parameters
case.setParameters( [
                        [ 'constant/reactiveProperties', 'k0', '%.5g' % A1, 'reaction02' ],
                        [ 'constant/reactiveProperties', 'k2', '%.5g' % A2, 'reaction02' ],
                        [ 'constant/reactiveProperties', 'EA', '%.5g' % EA1, 'reaction02' ],
                        [ 'constant/reactiveProperties', 'E2', '%.5g' % EA2, 'reaction02' ],
                        [ 'constant/reactiveProperties', 'KEq', '%.5g' % KEq, 'reaction02' ], 
                        [ 'constant/reactiveProperties', 'cCuCl', '%.5g' % cCuCl, 'reaction02' ],
                        [ 'constant/reactiveProperties', 'sHr', '%.5g' % sHr, 'reaction02' ],
                    ] )

# 3) thermophysicalProperties parameters
case.setParameters( [
                        [ 'constant/thermophysicalProperties', 'molWeight', '%.5g' % MgIn, '' ],
                    ] )

# 4) fvSchemes parameters
case.setParameters( [
                        [ 'system/fvSchemes', 'default', 'bounded Gauss SFCD' , 'divSchemes' ], 
                    ] )
                
# 5) 

# # -- different cellSize meshes baseCase list
# meshesCaseLst = []

# # -- prepare prototype mesh for each cellSize
# if makeMesh:
#     for cellSize in cellSizeLst:
#         meshCase = OpenFOAMCase()
#         meshCase.loadOFCaseFromBaseCase(baseCaseDir)
#         meshDir = '%s/protoMesh/%g'%(outFolder,cellSize)
#         meshCase.changeOFCaseDir(meshDir)
#         meshCase.copyBaseCase()
#         replaceInBlockMesh = ["system/blockMeshDict",['length1', 'length2', 'width','nDiscX','nDiscYZ'],[str(length1),str(length2),str(width),str(int((length1+length2)/cellSize)),str(int(2*width/cellSize))]]
#         replaceInSnappyHexMesh = ["system/snappyHexMeshDict", ['spR'], [str(R)]]             
#         meshCase.replace([replaceInBlockMesh, replaceInSnappyHexMesh])
#         meshCase.setParameter(['system/controlDict', 'endTime', str(endTime), ''])
#         meshCase.setParameter(['system/decomposeParDict', 'numberOfSubdomains', str(nProc), ''])
#         meshCase.runCommands(['chmod u=rwx All*', './Allmesh'])
#         meshesCaseLst.append(meshCase)
#     if not runSim: sys.exit()

        
     