# python script to create and run test simulation for trasnient simulation using reactHetCat

from OF_caseClass import OpenFOAMCase
import sys
import numpy as np
import os

solvers = ['rhoPimpleFoam']

TIn = 400
UIn = 1
pOut = 101325

porEps = 0.4
tort = 2
dP = 1e-5
kappa = 5

specieNames = np.array(["N2"])
yIn = np.array([1.0])
MolMass = np.array([ 28.01e-3])
MgIn = np.sum( yIn * MolMass )

sigmaVs = np.array([ 18.5 ])

for solver in solvers:
    caseOutFolder = '../ZZ_cases/test%s' % solver

    tutorialFolder = '../baseCases/baseCaseMass/'
    
    testCase = OpenFOAMCase()
    testCase.loadOFCaseFromBaseCase(tutorialFolder)
    testCase.changeOFCaseDir(caseOutFolder)
    testCase.copyBaseCase()

    # -- set parameters:
    # 1) inlet parameters
    testCase.runCommands( [ 'rm 0.org/bsChemSp' ] )
    testCase.setParameters( [ 
        [ '0.org/T', 'internalField', 'uniform %.5g' % TIn, '' ],
        [ '0.org/U', 'internalField', 'uniform (%.5g 0 0)' % UIn, ''],
        [ '0.org/p', 'internalField', 'uniform %.5g' % pOut, ''],
    ] )

    # 2) geometry genereation parameters
    # testCase.replace (  [ 
    #     [ 'system/blockMeshDict', [ 'nX', 'nY'], [ '%d' % nX, '%d' % nY ] ]
    # ] )

    # 3) transportProperties parameters
    # -- chemical species parameters
    # for nameInd in range(len(specieNames)):
    #     name = specieNames[nameInd]
    #     case.setParameters( [ 
    #     [ 'constant/transportProperties', 'molM', '%.5g' % MolMass[nameInd], name ],
    #     [ 'constant/transportProperties', 'sigmaV', '%.5g' % sigmaVs[nameInd], name ],
    # ])
    # case.setParameters( [
    #     [ 'constant/transportProperties', 'species', '(%s)' % namesStr[:-1].replace('Mass',''), '' ]
    # ] )
    # -- porousZone parameters
    testCase.setParameters( [
        [ 'constant/transportProperties', 'porEps', '%.5g' % porEps, '' ],
        [ 'constant/transportProperties', 'tort', '%.5g' % tort, '' ],
        [ 'constant/transportProperties', 'dP', '%.5g' % dP, '' ],
        [ 'constant/transportProperties', 'kappa', '%.5g' % kappa, 'coatZone' ],
    ] )

    # 4) reactionProperties parameters
    # testCase.setParameters( [
    #     [ 'constant/reactiveProperties', 'k0', '%.5g' % A1, 'reaction02' ],
    #     [ 'constant/reactiveProperties', 'k2', '%.5g' % A2, 'reaction02' ],
    #     [ 'constant/reactiveProperties', 'EA', '%.5g' % EA1, 'reaction02' ],
    #     [ 'constant/reactiveProperties', 'E2', '%.5g' % EA2, 'reaction02' ],
    #     [ 'constant/reactiveProperties', 'KEq', '%.5g' % KEq, 'reaction02' ], 
    #     [ 'constant/reactiveProperties', 'cCuCl', '%.5g' % cCuCl, 'reaction02' ],
    #     [ 'constant/reactiveProperties', 'sHr', '%.5g' % sHr, 'reaction02' ],
    # ] )

    # 5) thermophysicalProperties parameters
    testCase.setParameters( [
        [ 'constant/thermophysicalProperties', 'molWeight', '%.5g' % (MgIn*1000), '' ],
        [ 'constant/thermophysicalProperties', 'Cp', '%.5g' % (CpInMass), '' ],
    ] )

    # 6) fvSchemes parameters
    case.setParameters( [
        [ 'system/fvSchemes', 'default', divScheme , 'divSchemes' ], 
        # [ 'system/fvSchemes', 'default', 'bounded Gauss upwind phi' , 'divSchemes' ], 
    ] )
                    
    # 7) fvSolution parameters
    if not meshDone:
        case.replace([
            [ 'system/fvSolution', ['customSolver'], [namesStr[:-1].replace(' ', '|')] ], 
        ] )
    case.setParameters([
        [ 'system/fvSolution', 'nConcCorrectors', '%d' %nConc , 'SIMPLE' ], 
        [ 'system/fvSolution', 'nTempCorrectors', '%d' %nTemp , 'SIMPLE' ], 
    ] )

    # 8) decomposeParDict parameters
    case.setParameters([
        [ 'system/decomposeParDict', 'numberOfSubdomains', '%d' %nProc , '' ], 
    ] )

    # 9) residuals
    if not meshDone:
        case.replace([
            [ 'system/residuals', [ 'flsLst' ], [ 'U p T %s' % namesStr[:-1] ] ], 
        ]) 

    # 10) controlDict
    case.setParameters([
        [ 'system/controlDict', 'application', solver, '' ], 
        [ 'system/controlDict', 'endTime', '%d' % endTime, '' ], 
        [ 'system/controlDict', 'writeInterval', '%d' % wrInt, '' ], 
    ]) 

    # 11) sample
    case.setParameters([
        [ 'system/sample', 'fields', "(U p T %s)" % namesStr[:-1] , '' ], 
    ] )