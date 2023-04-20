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

meshDone = True
# meshDone = False

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
# baseCaseDir = '../ZZ_cases/testRashV1'
# baseCaseDir = '../ZZ_cases/testRashV4'
baseCaseDir = '../ZZ_cases/V2meshV1'
baseCaseDir = '../ZZ_cases/V2meshV1_reform'
baseCaseDir = '../ZZ_cases/V2meshV1_vagon'

# baseCaseDir = '../ZZ_cases/V2meshV16_betNum'

outFolder = '../ZZ_cases/V2meshV3_moreIt'
outFolder = '../ZZ_cases/V2meshV4_testNewRepl'
outFolder = '../ZZ_cases/V2meshV5_EqRelaxUpdate'
outFolder = '../ZZ_cases/V2meshV6_highHTcoeff'
outFolder = '../ZZ_cases/V2meshV7_highHTcoeffHK'
outFolder = '../ZZ_cases/V2meshV8_6but'
outFolder = '../ZZ_cases/V2meshV9_8DiffKin'
outFolder = '../ZZ_cases/V2meshV11_moreIter'
outFolder = '../ZZ_cases/V2meshV12_betterNumerics'
outFolder = '../ZZ_cases/V2meshV13_noInialGuess'
outFolder = '../ZZ_cases/V2meshV17_diffSolver'

outFolder = '../ZZ_cases/V2reformV3'
outFolder = '../ZZ_cases/V2vagonV3'
# outFolder = '../ZZ_cases/testRashV5'
# outFolder = '../ZZ_cases/testRashV6'
# outFolder = '../ZZ_cases/testRashV7'

# -- inlet parameters
specieNames = np.array(["ethylene", "O2", "HCl", "C2H4Cl2", "H2O", "N2"])
scalarfields = ['T']
vectorfields = ['U']
namesStr = ''
namesSolver = ''
for nameInd in range(len(specieNames)):
    namesStr += '%sMass ' % specieNames[nameInd]
    scalarfields.append('%sMass' % specieNames[nameInd])
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

CpIn = np.array([61.07, 29.15, 29.64, 87.43, 35.13, 29.49]) # J mol**-1 K**-1
CpInM = np.sum(CpIn * yIn)  # J mol**-1 K**-1
CpInMass = CpInM / MgIn

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
cCuCl = 0.34
cCuCl = 0.4

# -- geometry generation parameters
nCellsBetweenLevels = 5 # 4
rashLvl = '(1 2)' # (1 2)
cylLvl = '(1 2)' # (1 1)
nX = 143
nY = 42

# -- rashigs zone parameters
porEps = 0.42
tort = 5
dP = 5e-9
kappa = 5

# -- numerics and computing
nConc = 0
nTemp = 0
nBoth = 1
nProc = 16
eqTRelx = "0.9999"
# endTime = 3000
# endTime = 250 
endTime = 200 
# endTime = 10
# wrInt = 3000
# wrInt = 100
wrInt = 200
divScheme = 'bounded Gauss SFCD'
# divScheme = 'bounded Gauss upwind'
# divScheme = 'bounded Gauss linearUpwind grad(U);'

# -- parameters for the whole reactor
N = 33
Nzacatek = 0
fields = "'(T U p %s)'" % namesStr[:-1]

# -- update of the kinetics during the whole reactor simulation
whenChange = [0, 7, 22, 26]
whenChange = [0, 7, 26]
howMuch = [1, 1.36, 1.67, 1.81]
howMuch = [1, 1.2, 1.9]

# ---------------------------------------------------------------------------------------------------------------------
# -- set the parameters

# -- create new folder from basecase folder
case = OpenFOAMCase()
case.loadOFCaseFromBaseCase(baseCaseDir)
case.changeOFCaseDir(outFolder)
case.copyBaseCase()

# -- set parameters:
# 1) inlet parameters
if not meshDone:
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
    # case.setParameters( [
        # [ 'system/snappyHexMeshDict', 'nCellsBetweenLevels', '%d' % nCellsBetweenLevels, '' ],
        # [ 'system/snappyHexMeshDict', 'level', rashLvl, 'rashigs' ],
        # [ 'system/snappyHexMeshDict', 'level', cylLvl, 'cylinder' ],
    # ] )
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
    # [ 'system/fvSolution', 'nTogCorrectors', '%d' %nBoth , 'SIMPLE' ], 
    # [ 'system/fvSolution', 'T', eqTRelx , 'equations' ], 
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


# sys.exit()
# -------------------------------------------------------------------------------------------
# -- prepare mesh
if not meshDone:
    case.runCommands(   [
                            'chmod 755 -R ./*',
                            './Allrun-parallel'
                        ] )
# else:
#     case.runCommands(   [
#                             'decomposePar > log.decomposePar',
#                             'foamJob -parallel -screen %s' % solver,
#                         ] )
if Nzacatek == 0:
    case.runCommands([
        'rm -rf constant/boundaryData',
    ])
# -- the whole reactor
case.runCommands([
    # 'rm -rf processor*',
    'rm -rf processor*',
    'decomposePar > log.decomposePar',
    # 'foamJob -parallel -screen %s' % solver,
])

if Nzacatek > 0:
    case.updateTimes()
    lT = case.latestTime
    case.setParameters([
        [ 'system/controlDict', 'endTime', '%d' % (lT + endTime), '' ], 
        [ 'system/controlDict', 'writeInterval', '%d' % (lT + endTime), '' ], 
    ]) 
else:
    case.setParameters([
        [ 'system/controlDict', 'endTime', '%d' % endTime, '' ], 
        [ 'system/controlDict', 'writeInterval', '%d' % endTime, '' ], 
    ]) 

for sim in range(Nzacatek,N):
    print('Running sim. %d %g -- %g m'%(sim,(sim)*0.1,(sim+1)*0.1))
    for j in range(len(whenChange)):
        if whenChange[j] == sim:
            case.setParameters([
                [ 'constant/reactiveProperties', 'cCuCl', '%.5g' % (cCuCl * howMuch[j]), 'reaction02' ], 
            ])
    case.runCommands([
        'foamJob -parallel -screen %s > log.rash%d' % (solver, sim),
        "reconstructPar -latestTime -fields %s > log.reconstruct%d" % ( fields, sim ),
        "postProcess -func 'sample' > log.sample%d" % sim, 
    ])
    case.updateTimes()
    lT = case.latestTime
    case.runCommands(
        [
            'cp -r %d 0.%d' %(lT, lT)
        ]
    )
    print('Result in 0.%d folder' % lT)
    # case.runCommands([
    #     'cp -r 0 %d' % (lT+1)
    # ])
    # case.updateTimes()
    lT = case.latestTime
    if sim == 0:
        case.runCommands([
            'mkdir constant/boundaryData',
            'mkdir constant/boundaryData/inlet',
        ])
    case.runCommands([
        'mkdir constant/boundaryData/inlet/%d' % lT,
        'cp postProcessing/sample/%d/outlet_field/points constant/boundaryData/inlet/' % lT,
        # 'rm -f constant/boundaryData/inlet/%d/*' % lT,
    ])
    for field in scalarfields:
        case.runCommands([
            'cp postProcessing/sample/%d/outlet_field/scalarField/%s constant/boundaryData/inlet/%d/' % ( lT, field, lT ) 
        ])
    case.runCommands([
        "foamJob -parallel -screen postProcess -func 'patchAverage(p,name=inlet,patch=inlet)' > log.postProcessP%d" % sim 
    ])
    with open ('%s/postProcessing/patchAverage(p,name=inlet,patch=inlet)/0/surfaceFieldValue.dat' % case.dir, 'r') as fl:
        lines = fl.readlines()
        print(lines[-1].split("\t")[-1].replace('\n',''))
        dp = float(lines[-1].split("\t")[-1].replace('\n','')) - pOut
    
    # case.setParameters( [
    #     [ '%d/p' % lT, 'p0', 'uniform %.5g' % (pOut - dp), 'outlet'],
    #     [ '%d/p' % lT, 'value', 'uniform %.5g' % (pOut - dp), 'outlet'],
    # ] )

    # pOut = pOut - dp
    print('POut = %g' % pOut)

    # with open('%s/%d/p' % (case.dir, lT), 'r') as fl:
    #     lineBC = fl.readlines()
    #     indInlet = -1
    #     indZavorka = -1
    #     uz = 0
    #     for i in range(len(lineBC)):
    #         if 'outlet' in lineBC[i]:
    #             indInlet = i
    #             uz = 1
    #         if uz == 1 and '}' in lineBC[i]:
    #             if sim == 0:
    #                 indZavorka = i
    #                 uz = 0
    #             else:
    #                 uz = 2
    #         elif uz == 2 and '}' in lineBC[i]:
    #             if sim != 0:
    #                 indZavorka = i
    #             uz = 0
    #     newBC = []
    #     print('inlet at', indInlet,indZavorka)
    #     for i in range(indInlet+2):
    #         newBC.append(lineBC[i])
    #     newBC.append("\ntype\ttotalPressure;\np0\tuniform %g;\nvalue\tuniform %g;\nU\tU;\nphi\tphi;\n" % (pOut, pOut))
    #     for i in range(indZavorka,len(lineBC)):
    #         newBC.append(lineBC[i])
    #     with open('%s/%d/p' % (case.dir, lT), 'w') as fl:
    #         for i in range(len(newBC)):
    #             fl.writelines(newBC[i])
    
    for field in vectorfields:
        case.runCommands([
            'cp postProcessing/sample/%d/outlet_field/vectorField/%s constant/boundaryData/inlet/%d/'%(lT,field,lT) 
        ])

    for field in scalarfields:
        with open('%s/%d/%s' % (case.dir, lT, field), 'r') as fl:
            lineBC = fl.readlines()
        indInlet = -1
        indZavorka = -1
        uz = 0
        for i in range(len(lineBC)):
            if 'inlet' in lineBC[i]:
                indInlet = i
                uz = 1
            if uz == 1 and '}' in lineBC[i]:
                if sim == 0:
                    indZavorka = i
                    uz = 0
                else:
                    uz = 2
            elif uz == 2 and '}' in lineBC[i]:
                if sim != 0:
                    indZavorka = i
                uz = 0
        newBC = []
        print('inlet at', indInlet,indZavorka)
        for i in range(indInlet+2):
            newBC.append(lineBC[i])
        newBC.append("\t\t\t\ntype\ttimeVaryingMappedFixedValue;\t\t\t\noffset 0;\t\t\t\nsetAverage off;\t\t\t\nmapMethod nearest;\n")
        for i in range(indZavorka,len(lineBC)):
            newBC.append(lineBC[i])
        with open('%s/%d/%s' % (case.dir, lT,field), 'w') as fl:
            for i in range(len(newBC)):
                fl.writelines(newBC[i])
        
    for field in vectorfields:
        with open('%s/%d/%s' % (case.dir, lT, field), 'r') as fl:
            lineBC = fl.readlines()
        indInlet = -1
        indZavorka = -1
        uz = 0
        for i in range(len(lineBC)):
            if 'inlet' in lineBC[i] and uz == 0:
                indInlet = i
                uz = 1
            if uz == 1 and '}' in lineBC[i]:
                if sim == 0:
                    indZavorka = i
                    uz = 3
                else:
                    uz = 2
            elif uz == 2 and '}' in lineBC[i]:
                if sim != 0:
                    indZavorka = i
                uz = 3
        print('inlet at', indInlet,indZavorka)
        newBC = []
        for i in range(indInlet+2):
            newBC.append(lineBC[i])
        newBC.append("\t\t\t\ntype\ttimeVaryingMappedFixedValue;\t\t\t\noffset (0 0 0);\t\t\t\nsetAverage off;\t\t\t\nmapMethod nearest;\n")
        for i in range(indZavorka,len(lineBC)):
            newBC.append(lineBC[i])
        with open('%s/%d/%s' % (case.dir, lT, field), 'w') as fl:
            for i in range(len(newBC)):
                fl.writelines(newBC[i])

    # for field in scalarfields:
        # # repaire here
        # case.runCommands([ "foamJob -parallel -screen postProcess -func 'patchAverage(%s,name=outlet,patch=outlet)' > log.postProcess%s%d" % (field,field,sim) ])
        # with open ('%s/postProcessing/patchAverage(%s,name=outlet,patch=outlet)/0/surfaceFieldValue.dat' % (case.dir, field), 'r') as fl:
        #     lines = fl.readlines()
        #     fCorr = float(lines[-1].split("\t")[-1].replace('\n',''))
        #     # print(fCorr)

        # case.runCommands([ "postProcess -func 'patchAverage(%s,name=inlet,patch=inlet)' > log.postProcess_o%s%d" % (field,field,sim) ])
        # with open ('%s/postProcessing/patchAverage(%s,name=inlet,patch=inlet)/0/surfaceFieldValue.dat' % (case.dir, field), 'r') as fl:
        #     lines = fl.readlines()
        #     # print(lines[-1].split("\t")[-1].replace('\n',''))
        #     fsouc = float(lines[-1].split("\t")[-1].replace('\n',''))
        # print(fCorr,fsouc,fCorr/fsouc) 

        # with open ('%s/constant/boundaryData/inlet/%d/%s' % (case.dir, lT, field) , 'r') as fl:
        #     lines = fl.readlines()
        # for i in range(3,len(lines)-1):
        #     # f = float(lines[i].replace('\n',''))
        #     f = float(lines[i].replace('\n',''))*float(fCorr)/float(fsouc)
        #     lines[i] = '%.15g\n'%(f)
        # with open ('%s/constant/boundaryData/inlet/%d/%s' % (case.dir, lT, field) , 'w') as fl:
        #     for i in range(len(lines)):
        #         fl.writelines(lines[i])
        # case.runCommands([ "postProcess -func 'patchAverage(%s,name=inlet,patch=inlet)' > log.postProcess_n%s%d"%(field,field,sim) ])

    # case.runCommands([ "postProcess -func 'patchAverage(U,name=outlet,patch=outlet)' > log.postProcessU%d" % sim ])
    # case.runCommands([ "foamJob -parallel -screen postProcess -func 'patchAverage(U,name=outlet,patch=outlet)' > log.postProcessU%d" % sim ])
    # with open ('%s/postProcessing/patchAverage(U,name=outlet,patch=outlet)/0/surfaceFieldValue.dat' % case.dir, 'r') as fl:
    #     lines = fl.readlines()
    # Ucorr = (lines[-1].split("\t")[-1].replace('\n','').replace('(','').replace(')','').split(' '))

    # case.runCommands([ "postProcess -func 'patchAverage(U,name=inlet,patch=inlet)' > log.postProcessU%d" % sim ])
    # with open ('%s/postProcessing/patchAverage(U,name=inlet,patch=inlet)/0/surfaceFieldValue.dat' % case.dir, 'r') as fl:
    #     lines = fl.readlines()
    # # print(lines[-1].split("\t")[-1].replace('\n',''))
    # Usouc = (lines[-1].split("\t")[-1].replace('\n','').replace('(','').replace(')','').split(' '))
    # print(Ucorr,Usouc,float(Ucorr[0])/float(Usouc[0]),float(Ucorr[1])/float(Usouc[1]),float(Ucorr[2])/float(Usouc[2]))


    # for field in scalarfields:
    #     with open('%s/%d/%s' % (case.dir, lT, field), 'r') as fl:
    #         lineBC = fl.readlines()
    #     indInlet = -1
    #     indZavorka = -1
    #     uz = 0
    #     for i in range(len(lineBC)):
    #         if 'inlet' in lineBC[i]:
    #             indInlet = i
    #             uz = 1
    #         if uz == 1 and '}' in lineBC[i]:
    #             if sim == 0:
    #                 indZavorka = i
    #                 uz = 0
    #             else:
    #                 uz = 2
    #         elif uz == 2 and '}' in lineBC[i]:
    #             if sim != 0:
    #                 indZavorka = i
    #             uz = 0
    #     newBC = []
    #     print('inlet at', indInlet,indZavorka)
    #     for i in range(indInlet+2):
    #         newBC.append(lineBC[i])
    #     newBC.append("\t\t\t\ntype\ttimeVaryingMappedFixedValue;\t\t\t\noffset 0;\t\t\t\nsetAverage off;\t\t\t\nmapMethod nearest;\n")
    #     for i in range(indZavorka,len(lineBC)):
    #         newBC.append(lineBC[i])
    #     with open('%s/%d/%s' % (case.dir, lT,field), 'w') as fl:
    #         for i in range(len(newBC)):
    #             fl.writelines(newBC[i])
        
    # for field in vectorfields:
    #     with open('%s/%d/%s' % (case.dir, lT, field), 'r') as fl:
    #         lineBC = fl.readlines()
    #     indInlet = -1
    #     indZavorka = -1
    #     uz = 0
    #     for i in range(len(lineBC)):
    #         if 'inlet' in lineBC[i] and uz == 0:
    #             indInlet = i
    #             uz = 1
    #         if uz == 1 and '}' in lineBC[i]:
    #             if sim == 0:
    #                 indZavorka = i
    #                 uz = 3
    #             else:
    #                 uz = 2
    #         elif uz == 2 and '}' in lineBC[i]:
    #             if sim != 0:
    #                 indZavorka = i
    #             uz = 3
    #     print('inlet at', indInlet,indZavorka)
    #     newBC = []
    #     for i in range(indInlet+2):
    #         newBC.append(lineBC[i])
    #     newBC.append("\t\t\t\ntype\ttimeVaryingMappedFixedValue;\t\t\t\noffset (0 0 0);\t\t\t\nsetAverage off;\t\t\t\nmapMethod nearest;\n")
    #     for i in range(indZavorka,len(lineBC)):
    #         newBC.append(lineBC[i])
    #     with open('%s/%d/%s' % (case.dir, lT, field), 'w') as fl:
    #         for i in range(len(newBC)):
    #             fl.writelines(newBC[i])



    # for field in scalarfields:
    #     # with open('%s/%d/%s' % (case.dir, lT, field), 'r') as fl:
    #     with open('%s/0/%s' % (case.dir, field), 'r') as fl:
    #         lineBC = fl.readlines()
    #     indInlet = -1
    #     indZavorka = -1
    #     uz = 0
    #     for i in range(len(lineBC)):
    #         if 'inlet' in lineBC[i]:
    #             indInlet = i
    #             uz = 1
    #         if uz == 1 and '}' in lineBC[i]:
    #             # if sim == 0:
    #             indZavorka = i
    #             uz = 0
    #             # else:
    #             #     uz = 2
    #         elif uz == 2 and '}' in lineBC[i]:
    #             # if sim != 0:
    #             #     indZavorka = i
    #             uz = 0
    #     newBC = []
    #     print('inlet at', indInlet,indZavorka)
    #     for i in range(indInlet+2):
    #         newBC.append(lineBC[i])
    #     newBC.append("\t\t\t\ntype\ttimeVaryingMappedFixedValue;\t\t\t\noffset 0;\t\t\t\nsetAverage off;\t\t\t\nmapMethod nearest;\n")
    #     for i in range(indZavorka,len(lineBC)):
    #         newBC.append(lineBC[i])
    #     with open('%s/%d/%s' % (case.dir, lT,field), 'w') as fl:
    #         for i in range(len(newBC)):
    #             fl.writelines(newBC[i])
        
    #     # case.replace(
    #     #     [
    #     #         ['%s/%d/%s' % (case.dir, lT,field), ['HODNOTA'], ['%.6g' % ]]
    #     #     ]
    #     # )
        
    # for field in vectorfields:
    #     with open('%s/0/%s' % (case.dir, field), 'r') as fl:
    #         lineBC = fl.readlines()
    #     indInlet = -1
    #     indZavorka = -1
    #     uz = 0
    #     for i in range(len(lineBC)):
    #         if 'inlet' in lineBC[i] and uz == 0:
    #             indInlet = i
    #             uz = 1
    #         if uz == 1 and '}' in lineBC[i]:
    #             # if sim == 0:
    #             indZavorka = i
    #             uz = 3
    #             # else:
    #             #     uz = 2
    #         elif uz == 2 and '}' in lineBC[i]:
    #             # if sim != 0:
    #             #     indZavorka = i
    #             uz = 3
    #     print('inlet at', indInlet,indZavorka)
    #     newBC = []
    #     for i in range(indInlet+2):
    #         newBC.append(lineBC[i])
    #     newBC.append("\t\t\t\ntype\ttimeVaryingMappedFixedValue;\t\t\t\noffset (0 0 0);\t\t\t\nsetAverage off;\t\t\t\nmapMethod nearest;\n")
    #     for i in range(indZavorka,len(lineBC)):
    #         newBC.append(lineBC[i])
    #     with open('%s/%d/%s' % (case.dir, lT, field), 'w') as fl:
    #         for i in range(len(newBC)):
    #             fl.writelines(newBC[i])

    # with open ('%s/constant/boundaryData/inlet/%d/U' % (case.dir, lT) , 'r') as fl:
    #     lines = fl.readlines()
    # for i in range(3,len(lines)-1):
    #     Ux = float(lines[i].split(' ')[0].replace('(',''))*float(Ucorr[0])/float(Usouc[0])#*pOut/(pOut+dp)
    #     Uy = float(lines[i].split(' ')[1])*float(Ucorr[1])/float(Usouc[1])*pOut/(pOut+dp)*0
    #     Uz = float(lines[i].split(' ')[2].replace(')\n',''))*float(Ucorr[2])/float(Usouc[2])*pOut/(pOut+dp)*0
    #     # Ux = float(lines[i].split(' ')[0].replace('(',''))#*float(Ucorr[0])/float(Usouc[0])#*pOut/(pOut+dp)
    #     # Uy = float(lines[i].split(' ')[1])
    #     # Uz = float(lines[i].split(' ')[2].replace(')\n',''))
    #     lines[i] = '(%.15g %.15g %.15g)\n'%(Ux,Uy, Uz)
    # with open ('%s/constant/boundaryData/inlet/%d/U' % (case.dir, lT) , 'w') as fl:
    #     for i in range(len(lines)):
    #         fl.writelines(lines[i])

    case.runCommands([ 'decomposePar -latestTime -fields > log.decompose%d' % (sim) ])
    case.setParameters([
        [ 'system/controlDict', 'endTime', '%d' % (lT + endTime), '' ], 
        [ 'system/controlDict', 'writeInterval', '%d' % (lT + endTime), '' ], 
    ]) 
             
        
     