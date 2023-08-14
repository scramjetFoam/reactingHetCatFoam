# -- python script to create test case for reactor pipe simulations

from OF_caseClass import OpenFOAMCase
import numpy as np
import math
from blockMeshDictClassV8 import *

# PARAMETERS ################################################
# -- openfoam base case directory
# baseCaseDir = '../ZZ_cases/testSOFC'
baseCaseDir = '../ZZ_cases/testSOFC'
outFolder = '../ZZ_cases/testSOFCV3'
# outFolder = '../ZZ_cases/testSOFCV4'

# -- geometry infromation
createMesh = False
createMesh = True
L = 1.2e-4                                                  # -- length of the geometry stl
LBuff = 30e-6                                               # -- length of the buffer infront of the geometry
r = 1.1e-4                                                  # -- side of the geometry 
dA = 5e-6
dX, dY, dZ = dA, dA, dA                                     # -- computational cell dimensions
x0 = y0 = z0 = 0.0                                          # -- coordinate system origin
grX = grY = grZ = "1.0"                                     # -- grading

# -- information about chemical species
specieNames = np.array(["H2", "O2", "H2O", "N2"])           # -- what reacts
# specieNames = np.array(["H2", "N2", "H2O", "O2"])           # -- what reacts
species = ''
for specieName in specieNames:
    species += specieName + ' '
# molMass = np.array([ 2e-3, 28e-3, 18e-3, 32e-3])           # -- molar masses
molMass = np.array([ 2e-3, 32e-3, 18e-3, 28e-3])           # -- molar masses
sigmaVs = np.array([ 6.0, 16.3, 13.0, 18.5])               # -- diffusion volumes (see Fuller in skripta Chemické inženýrství II)
sigmaVs = np.array([ 6.0, 18.5, 13.0, 16.3])               # -- diffusion volumes (see Fuller in skripta Chemické inženýrství II)
# nuVec = np.array([0, -0.5, 0, 0])                          # -- stoichiometric coefficients
nuVec = np.array([-1, -0.5, 1, 0])                          # -- stoichiometric coefficients
# nuVec = np.array([-1, 0, 1, 0])                          # -- stoichiometric coefficients
alphaVec = np.array([0.5, 0.5, 0, 0])                           # -- discus this with me
# alphaVec = np.array([0, 0, 0, 0])                           # -- discus this with me
yIn = np.array([0.5, 0.5, 1e-16, 1-0.5-0.5])              # -- inlet molar fractions
# yIn = np.array([0.3, 0.3, 0, 1-0.3-0.3])              # -- inlet molar fractions
yInInit = np.array([1e-16, 1e-16, 1e-16, 1-3*1e-16])              # -- inlet molar fractions

# -- other BC
Tin = 470
Tg = 470
alpha = 0
pOut = 101326
UIn = 0.0

# -- information about reaction zone
reactZones = ['solid']                             # -- name of the reaction zone
reacZonesProp = np.array(
    [
        ['porEps', '[0 0 0 0 0 0 0]', 0.0],                 # -- porosity         
        ['tort', '[0 0 0 0 0 0 0]', 4],                     # -- tortuosity
        ['dP', '[0 1 0 0 0 0 0]', 4e-13],                   # -- pore diameter
        ['kappa', '[1 1 -3 -1 0 0 0]', 2],                  # -- heat conductivity
        ['coatFrac', '[0 0 0 0 0 0 0]', 1],                 # -- coating fraction
        ['sigmaScInfEl', '[-1 -3 3 0 0 2 0]', 2e5],        # -- coeficients in electric conduictivity
        ['sigmaScEAEl', '[0 0 0 1 0 0 0]', -1150],          
        ['sigmaScInfIon', '[-1 -3 3 0 0 2 0]', 6.92e4],     
        ['sigmaScEAIon', '[0 0 0 1 0 0 0]', -9681],      
    ]
)

# -- reactions
reactions = ['reaction02']
kinSwitch = [4]
reactionsInfo = np.array(
    [
        ['j0A', '[0 -3 0 0 0 1 0]', 1e5],                 
        ['EA', '[1  2 -2 0 -1 0 0]', 10],                 
        ['m', '[0 0 0 0 0 0 0]', 0.5],                 
        ['F', '[0 0 1 0 -1 1 0]', 96485.0],                 
        ['alphaA', '[0 0 0 0 0 0 0]', 0.3],         
        # ['j0K', '[0 -3 0 0 0 1 0]', 1000],                 
        # ['EA', '[1  2 -2 0 -1 0 0]', 10],                 
        # ['m', '[0 0 0 0 0 0 0]', 0.5],                 
        # ['F', '[0 0 1 0 -1 1 0]', 96485.0],         
    ]
)

# -- gas dynamic viscosity 
As = 2.08e-06                                               # -- argon -- sutherland equation in openfoam implementation
Ts = 110.0  

# -- permeability of the packed bed
kappa = 5e-16

# -- heat conductivity
Cp = 520

# -- reaction
k0 = 1e12
# k0 = 0
dH = 0

# --numerics
relEqMass = 0.95
relEqT = 0.95

# CODE ITSELF ###############################################
# -- basic calculations
MgIn = np.sum(yIn * molMass)                                # -- inlet molar mass
wIn = yIn * molMass / MgIn                                  # -- inlet mass fractions

# -- create OpenFOAMCase object to change values in dictionaries
baseCase = OpenFOAMCase()
baseCase.loadOFCaseFromBaseCase(baseCaseDir)
baseCase.changeOFCaseDir(outFolder)
baseCase.copyBaseCase()

# -- preparation of the blockMeshDictFile for geometry generation
fvMesh = mesh()                                             # -- create mesh

### BLOCKS ###
# -- first block
xC, yC, zC = x0, y0-LBuff, z0
xE, yE, zE = xC+r, y0+L, zC+r
# xE, yE, zE = xC+r, y0+L, zC+dA

# vertices
vertices = [
        [xC, yC, zC],
        [xE, yC, zC],
        [xE, yE, zC],
        [xC, yE, zC],
        [xC, yC, zE],
        [xE, yC, zE],
        [xE, yE, zE],
        [xC, yE, zE],
    ]

# neighbouring blocks
neighbours = []

# number of cells
nCX = int(round(abs(xE-xC)/dX))
nCY = int(round(abs(yE-yC)/dY))
nCZ = int(round(abs(zE-zC)/dZ))
nCells = [nCX, nCY, nCZ]

# grading
grading = [grX, grY, grZ]

# create the block
firstMid = fvMesh.addBlock(vertices, neighbours, nCells, grading)

## PATCHES ###

# -- inlet
inlet = list()
inlet.append(firstMid.retFXZ0())

fvMesh.addPatch("inletY", "patch", inlet)

# -- outlet
outlet = list()
outlet.append(firstMid.retFXZE())

fvMesh.addPatch("outletY", "patch", outlet)

# -- walls
walls = list()
walls.append(firstMid.retFXY0())
walls.append(firstMid.retFXYE())

fvMesh.addPatch("wall", "wall", walls)

inletX = list()
inletX.append(firstMid.retFYZ0())

fvMesh.addPatch("inletX", "patch", inletX)

outletX = list()
outletX.append(firstMid.retFYZE())

fvMesh.addPatch("outletX", "patch", outletX)

### WRITE ###
fvMesh.writeBMD("%s/system/" % baseCase.dir)

baseCase.replace(
    [
        ['system/snappyHexMeshDict', ['coatZone'], [reactZones[0]]]
    ]
)

# 1) BOUNDARY CONDITIONS
# -- creation of the boundary conditions for chemical species
for chSpI in range(len(specieNames)):
    name = specieNames[chSpI]
    baseCase.runCommands(
        [
            'cp 0.org/bsChemSp 0.org/%sMass' % name         # -- create file with BC
        ]
    )
    baseCase.replace(                                       # -- change values in BC
        [ 
            # ["0.org/%sMass" % ( name ), ['nameSet', 'wChSpSetInit', 'wChSpSet'], [ str(name),  '%.5g' % (wIn[chSpI]), '%.5g' % (wIn[chSpI])] ]
            ["0.org/%sMass" % ( name ), ['nameSet', 'wChSpSetInit', 'wChSpSet'], [ str(name),  '%.5g' % (yInInit[chSpI]), '%.5g' % (wIn[chSpI])] ]
        ] 
    )
baseCase.runCommands(                                                           
    [
        'rm 0.org/bsChemSp'                                # -- remove baseBC
    ]
)

# -- update BC for velocity, pressure and temperature
baseCase.replace(
    [
        ['0.org/T', ['isoT'] ,['%g' % Tin]],
        ['0.org/p', ['101325'] ,['%g' % pOut]],
    ]
)

# 2) constant/porosityProperties
# -- set permeability
baseCase.replace(
    [
         ['constant/porosityProperties', ["reactingCellZone"], [reactZones[0]]]
    ]
)
baseCase.setParameters(
    [
        ["constant/porosityProperties", "d", "(%g %g %g)" % (1/kappa,1/kappa,1/kappa), "DarcyForchheimerCoeffs"]
    ]
)

# 3) constant/transportProperties
# -- add chemical species to transportProperties 
baseCase.addToDictionary( 
    [
        [ 'constant/transportProperties', 'species (%s);\n' % species, ''],
    ]
)
for nameInd in range(len(specieNames)):
    name = specieNames[nameInd]
    baseCase.addToDictionary( 
        [
            [ 'constant/transportProperties', '%s\n{\n}\n' % name, ''],
            [ 'constant/transportProperties', 'sigmaV\t%g;\n' % sigmaVs[nameInd], name ],
            [ 'constant/transportProperties', 'molM\t%g;\n' % molMass[nameInd], name ],
            # [ 'constant/transportProperties', 'nuVec\t(0 0 %g);\n' % nuVec[nameInd], name ],
            # [ 'constant/transportProperties', 'alphaVec\t(0 0 %g);\n' % alphaVec[nameInd], name ],
            [ 'constant/transportProperties', 'nuVec\t(0 %g);\n' % nuVec[nameInd], name ],
            [ 'constant/transportProperties', 'alphaVec\t(0 %g);\n' % alphaVec[nameInd], name ],
        ]
    )

# -- define zone parameters in transport properties
reacZones = '('
for i in range(len(reactZones)):
    reacZones += '%s ' % reactZones[i]
    baseCase.addToDictionary(
        [
            ['constant/transportProperties', '\n%s\n{\n}\n' % reactZones[i], ''],
        ]
    )
reacZones = reacZones[:-1] + ')'

baseCase.setParameters(
    [
        ['constant/transportProperties', 'zones', reacZones, ''],
    ]
)

for rowI in range(reacZonesProp.shape[0]):
    for colI in range(reacZonesProp.shape[1]-2):
        baseCase.addToDictionary(
            [
                ['constant/transportProperties', '%s %s %s %g;\n' % (reacZonesProp[rowI, 0], reacZonesProp[rowI, 0], reacZonesProp[rowI,1], float(reacZonesProp[rowI,2+colI])), reactZones[colI]]
            ]
        )

# 4) constant/thermophysicalProperties
baseCase.setParameters( 
    [
        [ 'constant/thermophysicalProperties', 'molWeight', '%.5g' % (MgIn*1000), '' ],
        [ 'constant/thermophysicalProperties', 'Cp', '%g' % Cp, '' ],
        [ 'constant/thermophysicalProperties', 'As', '%g' % As, '' ],
        [ 'constant/thermophysicalProperties', 'Ts', '%g' % Ts, '' ],
    ] 
)

# 5) constant/reactiveProperties
baseCase.setParameters(
    [
        # ['constant/reactiveProperties', 'activeReacts', '(0 0 1)', ''],
        ['constant/reactiveProperties', 'activeReacts', '(0 1)', ''],
        ['constant/reactiveProperties', 'k0', '%g' % k0, 'reaction01'],
        ['constant/reactiveProperties', 'sHr', '%g' % dH, 'reaction01'],
    ]
)
# -- create reactions 
for i in range(len(reactions)):
    baseCase.addToDictionary(
        [
            ['constant/reactiveProperties', '\n\n\t%s\n\t{\n\n\t}\n' % reactions[i], 'reactProps'],
            ['constant/reactiveProperties', '\tkinSetupSwitch %d;\n' % kinSwitch[i], '%s'%reactions[i]],
        ]
    )

for rowI in range(reactionsInfo.shape[0]):
    for colI in range(reactionsInfo.shape[1]-2):
        baseCase.addToDictionary(
            [
                ['constant/reactiveProperties', '\t%s %s %s %g;\n' % (reactionsInfo[rowI, 0], reactionsInfo[rowI, 0], reactionsInfo[rowI,1], float(reactionsInfo[rowI,2+colI])), reactions[colI]]
            ]
        )

# 6) system/fvSolution
baseCase.replace(
    [
        ['system/fvSolution',['customSolver'],['customSolver|%s' %species.replace(' ','Mass|')]],
    ]
)
baseCase.setParameters(
    [
        ['system/fvSolution', 'nConcCorrectors', '1', ''],
        ['system/fvSolution', 'nTempCorrectors', '0', ''],
        ['system/fvSolution', 'nTogCorrectors', '0', ''],
        ['system/fvSolution', 'nNonOrthogonalCorrectors', '0', ''],
    ]
)
baseCase.replace(
    [
        # ['system/fvSolution', ['"(customSolver|COMass|O2Mass|CO2Mass|ArMass|)" 	 0.9;'], ['"(customSolver|COMass|O2Mass|CO2Mass|ArMass|)"  %g;' % relEqMass]],
        ['system/fvSolution', ['"(T)" 0.9995;'], ['"(T)" %g;' % relEqT]]
    ]
)

# 7) system/fvSchemes
# baseCase.replace(
#     [
#         ['system/fvSchemes', ['corrected 0.5'], ['uncorrected']],
#     ]
# )

# # others
baseCase.replace(
    [
        [ 'system/residuals', [ 'customSolver' ], [ 'U p T elPhi elIon %s' %species.replace(' ','Mass\n') ] ], 
    ]
) 

if createMesh:
    baseCase.runCommands(
        [
            'chmod 755 ./* -R',
            'blockMesh > log.blockMesh',
            # 'rm -rf 0',
            # 'cp -r 0.org 0',
            # 'decomposePar > log.decomposePar',
            'rm -rf 0',
            'cp -r 0.org 0',
            'paraFoam -touch',
            'topoSet > log.topoSet',
            'createPatch -overwrite > log.createPatch',
            'snappyHexMesh -overwrite > log.snappyHexMesh',
            # 'setFields > log.setFields',
        ]
    )

# -- run the simulation
baseCase.runCommands(
    [
        'chmod 755 ./* -R',

        # 'topoSet > log.topoSet',
        # 'reactingHetCatSimpleFoamM'
    ]
)
