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
specieNames = np.array(["ethylene", "O2", "HCl", "C2H4Cl2", "H2O", "N2"])
namesStr = ''
namesSolver = ''
for nameInd in range(len(specieNames)):
    namesStr += '%sMass ' % specieNames[nameInd]
nIn = np.array([ 6.521E-03, 2.588E-03, 1.281E-02, 5.175E-03, 5.175E-03, 1.947E-02 ])
yIn = nIn / np.sum(nIn)
for yInInd in range(len(specieNames)):
    print( "Inlet molar fraction of the specie %s is %g" % (specieNames[yInInd], yIn[yInInd]))

MolMass = np.array([ 28.05e-3, 32.0e-3, 36.46e-3, 98.96e-3, 18.01e-3, 28.01e-3])
MgIn = np.sum( yIn * MolMass )
print( "Inlet molar mass is %g." %MgIn )

sigmaVs = np.array([ 36.6, 16.3, 23.3, 83.2, 13.1, 18.5 ])

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

# -- geometry generation parameters
nCellsBetweenLevels = 4 # 4
rashLvl = '(1 2)' # (1 2)
cylLvl = '(1 1)' # (1 1)
nX = 110
nY = 32

# -- rashigs zone parameters
porEps = 0.42
tort = 5
dP = 5e-9
kappa = 5

# -- numerics and computing
nConc = 2
nTemp = 2
nProc = 12
endTime = 10

# ---------------------------------------------------------------------------------------------------------------------
# -- set the parameters

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
case.runCommands( [ 'rm 0.org/bsChemSp' ] )
case.setParameters( [ 
                        [ '0.org/T', 'internalField', 'uniform %.5g' % TIn, '' ],
                        [ '0.org/U', 'internalField', 'uniform (%.5g 0 0)' % UIn, ''],
                        [ '0.org/p', 'internalField', 'uniform %.5g' % pOut, ''],
                    ] )

# 2) geometry genereation parameters
case.setParameters( [
                        [ 'system/snappyHexMeshDict', 'nCellsBetweenLevels', '%d' % nCellsBetweenLevels, '' ],
                        [ 'system/snappyHexMeshDict', 'level', rashLvl, 'rashigs' ],
                        [ 'system/snappyHexMeshDict', 'level', cylLvl, 'cylinder' ],
                    ] )
case.replace (  [ 
                    [ 'system/blockMeshDict', [ 'nX', 'nY'], [ '%d' % nX, '%d' % nY ] ]
                ] )

# 3) transportProperties parameters
# -- chemical species parameters

for nameInd in range(len(specieNames)):
    name = specieNames[nameInd]
    case.setParameters( [ 
                            [ 'constant/transportProperties', 'molM', '%.5g' % MolMass[nameInd], name ],
                            [ 'constant/transportProperties', 'sigmaV', '%.5g' % sigmaVs[nameInd], name ],
                        ])
case.setParameters( [
                        [ 'constant/transportProperties', 'species', '(%s)' % namesStr[:-1].replace('Mass',''), '' ]
                    ] )
# -- porousZone parameters
case.setParameters( [
                        [ 'constant/transportProperties', 'porEps', '%.5g' % porEps, '' ],
                        [ 'constant/transportProperties', 'tort', '%.5g' % tort, '' ],
                        [ 'constant/transportProperties', 'dP', '%.5g' % dP, '' ],
                        [ 'constant/transportProperties', 'kappa', '%.5g' % kappa, 'coatZone' ],
                    ] )

# 4) reactionProperties parameters
case.setParameters( [
                        [ 'constant/reactiveProperties', 'k0', '%.5g' % A1, 'reaction02' ],
                        [ 'constant/reactiveProperties', 'k2', '%.5g' % A2, 'reaction02' ],
                        [ 'constant/reactiveProperties', 'EA', '%.5g' % EA1, 'reaction02' ],
                        [ 'constant/reactiveProperties', 'E2', '%.5g' % EA2, 'reaction02' ],
                        [ 'constant/reactiveProperties', 'KEq', '%.5g' % KEq, 'reaction02' ], 
                        [ 'constant/reactiveProperties', 'cCuCl', '%.5g' % cCuCl, 'reaction02' ],
                        [ 'constant/reactiveProperties', 'sHr', '%.5g' % sHr, 'reaction02' ],
                    ] )

# 5) thermophysicalProperties parameters
case.setParameters( [
                        [ 'constant/thermophysicalProperties', 'molWeight', '%.5g' % MgIn, '' ],
                    ] )

# 6) fvSchemes parameters
case.setParameters( [
                        [ 'system/fvSchemes', 'default', 'bounded Gauss SFCD' , 'divSchemes' ], 
                    ] )
                
# 7) fvSolution parameters
case.replace(   [
                    [ 'system/fvSolution', ['customSolver'], [namesStr[:-1].replace(' ', '|')] ], 
                ] )
case.setParameters( [
                        [ 'system/fvSolution', 'nConcCorrectors', '%d' %nConc , 'SIMPLE' ], 
                        [ 'system/fvSolution', 'nTempCorrectors', '%d' %nTemp , 'SIMPLE' ], 
                    ] )

# 8) decomposeParDict parameters
case.setParameters( [
                        [ 'system/decomposeParDict', 'numberOfSubdomains', '%d' %nProc , '' ], 
                    ] )
# 9) residuals
case.replace(   [
                    [ 'system/residuals', [ 'flsLst' ], [ namesStr[:-1] ] ], 
                ] ) 

# 10) controlDict
case.setParameters( [
                        [ 'system/controlDict', 'application', solver, '' ], 
                        [ 'system/controlDict', 'endTime', '%d' % endTime, '' ], 
                    ] ) 



# -------------------------------------------------------------------------------------------
# -- prepare mesh and run the simulation
case.runCommands(   [
                        'chmod 755 -R ./*',
                        './Allrun-parallel'
                    ] )

        
     