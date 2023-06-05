# -- python script to create test case for reactor pipe simulations

from OF_caseClass import OpenFOAMCase
import numpy as np
import math
from blockMeshDictClassV8 import *

# PARAMETERS ################################################
# -- openfoam base case directory
baseCaseDir = '../baseCases/baseCaseMass/'
outFolder = '../ZZ_cases/testForLucka1'

# -- geometry infromation
L = 1                                                       # -- length of the tube
r = 0.05                                                    # -- radius of the tube 
wAng = 1.0                                                  # -- wedge angle
dX, dY, dZ = 3e-3, 3e-3, 3e-3                               # -- computational cell dimensions
x0 = y0 = z0 = 0.0                                          # -- coordinate system origin
grX = grY = grZ = "1.0"                                     # -- grading

# -- information about chemical species
specieNames = np.array(["CO", "O2", "CO2", "Ar"])           # -- what reacts
species = ''
for specieName in specieNames:
    species += specieName + ' '
molMass = np.array([ 28e-3, 32e-3, 44e-3, 40e-3])           # -- molar masses
sigmaVs = np.array([ 18.0, 16.3, 26.9, 16.2])               # -- diffusion volumes (see Fuller in skripta Chemické inženýrství II)
nuVec = np.array([-1, -0.5, 1, 0])                          # -- stoichiometric coefficients
alphaVec = np.array([0, 0, 0, 0])                           # -- discus this with me
yIn = np.array([0.01, 0.005, 0, 1-0.01-0.005])              # -- inlet molar fractions

# -- other BC
Tin = 470
Tg = 470
alpha = 10
pOut = 101326
UIn = 0.5

# -- information about reaction zone
reactZones = ['reactZoneLucka']                             # -- name of the reaction zone
reacZonesProp = np.array(
    [
        ['porEps', '[0 0 0 0 0 0 0]', 0.4],                 # -- porosity         
        ['tort', '[0 0 0 0 0 0 0]', 4],                     # -- tortuosity
        ['dP', '[0 1 0 0 0 0 0]', 4e-6],                    # -- pore diameter
        ['kappa', '[1 1 -3 -1 0 0 0]', 1],                  # -- heat conductivity
        ['coatFrac', '[0 0 0 0 0 0 0]', 0.5],               # -- coating fraction
    ]
)

# -- gas dynamic viscosity 
As = 2.08e-06                                               # -- argon -- sutherland equation in openfoam implementation
Ts = 110.0  

# -- permeability of the packed bed
kappa = 5e-8

# -- heat conductivity
Cp = 520

# -- reaction
k0 = 1e15
KInh = 100

# --numerics
relEqMass = 0.9
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
yMax = dZ/math.atan(wAng/180*math.pi*0.5)                   # -- yMax from wedge angle
fvMesh = mesh()                                             # -- create mesh
nCZ = 1                                                     # -- wedge

### BLOCKS ###
# -- first block
xC, yC = x0, y0
xE, yE = xC+L, yC+r

# vertices
vertices = [
        [xC, yC, z0-yC/yMax*dZ],
        [xE, yC, z0-yC/yMax*dZ],
        [xE, yE, z0-yE/yMax*dZ],
        [xC, yE, z0-yE/yMax*dZ],
        [xC, yC, z0+yC/yMax*dZ],
        [xE, yC, z0+yC/yMax*dZ],
        [xE, yE, z0+yE/yMax*dZ],
        [xC, yE, z0+yE/yMax*dZ],
    ]

# neighbouring blocks
neighbours = []

# number of cells
nCX = int(round(abs(xE-xC)/dX))
nCY = int(round(abs(yE-yC)/dY))
nCells = [nCX, nCY, nCZ]

# grading
grading = [grX, grY, grZ]

# create the block
firstMid = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = reactZones[0])

### PATCHES ###
# -- wedge patches
wedgeZ0 = list()
for block in fvMesh.blocks:
    wedgeZ0.append(block.retFXY0())

fvMesh.addPatch("wedgeZ0", "wedge", wedgeZ0)

wedgeZE = list()
for block in fvMesh.blocks:
    wedgeZE.append(block.retFXYE())

fvMesh.addPatch("wedgeZE", "wedge", wedgeZE)

# -- inlet
inlet = list()
inlet.append(firstMid.retFYZ0())

fvMesh.addPatch("inlet", "patch", inlet)

# -- outlet
outlet = list()
outlet.append(firstMid.retFYZE())

fvMesh.addPatch("outlet", "patch", outlet)

# -- walls
walls = list()
walls.append(firstMid.retFXZE())

fvMesh.addPatch("wall", "wall", walls)

### WRITE ###
fvMesh.writeBMD("%s/system/" % baseCase.dir)

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
            ["0.org/%sMass" % ( name ), ['nameSet', 'wChSpSetInit', 'wChSpSet'], [ str(name),  '%.5g' % (wIn[chSpI]), '%.5g' % (wIn[chSpI])] ]
        ] 
    )
    baseCase.addToDictionary(                              
        [
            [ '0.org/%sMass' % name, '\n\t"(wedgeZ0|wedgeZE)"\n\t{\n\t\ttype wedge;\n\t}\n', 'boundaryField'],
            [ '0.org/%sMass' % name, '\n\t"(wall)"\n\t{\n\t\ttype zeroGradient;\n\t}\n', 'boundaryField'],
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
        ['0.org/U', ['inv'] ,['%g' % UIn]],
        ['0.org/T', ['isoT'] ,['%g' % Tin]],
        ['0.org/p', ['101325'] ,['%g' % pOut]],
    ]
)
baseCase.addToDictionary(                              
    [
        [ '0.org/U', '\n\t"(wedgeZ0|wedgeZE)"\n\t{\n\t\ttype wedge;\n\t}\n', 'boundaryField'],
        [ '0.org/p', '\n\t"(wedgeZ0|wedgeZE)"\n\t{\n\t\ttype wedge;\n\t}\n', 'boundaryField'],
        [ '0.org/T', '\n\t"(wedgeZ0|wedgeZE)"\n\t{\n\t\ttype wedge;\n\t}\n', 'boundaryField'],
        [ '0.org/U', '\n\t"(wall)"\n\t{\n\t\ttype noSlip;\n\t}\n', 'boundaryField'],
        [ '0.org/p', '\n\t"(wall)"\n\t{\n\t\ttype zeroGradient;\n\t}\n', 'boundaryField'],
        [ '0.org/T', '\n\t"(wall)"\n\t{\n\t\ttype zeroGradient;\n\t}\n', 'boundaryField'],
        [ '0.org/T', '\n\t"(wall)"\n\t{\n\t\ttype codedMixed;\n\trefValue uniform %g;\n\trefGradient uniform 0;\n\tvalueFraction uniform 0;\n\tvalue uniform %g;\n\tname robinFront;\ncode #{\n\t\tconst fvPatch& boundaryPatch = patch();\n\t\tconst fvBoundaryMesh& boundaryMesh = boundaryPatch.boundaryMesh();\n\t\tconst fvMesh& mesh = boundaryMesh.mesh();\n\t\tconst scalarField& delta = patch().deltaCoeffs();\n\t\tconst volScalarField& kappa((mesh.lookupObject<volScalarField>("kappaEff")));\n\t\tconst volScalarField& epsEff((mesh.lookupObject<volScalarField>("epsEff")));\n\t\tconst volScalarField& kappaG((mesh.lookupObject<volScalarField>("thermo:kappa")));\n\t\tthis->refValue() = %g;\n\t\tthis->refGrad() = 0.;\n\t\tthis->valueFraction() = 1.0/(1.0+(kappa+kappaG*epsEff)/(%g/delta));\n\t#};\n\t}\n'% (Tg,Tg,Tg,alpha), 'boundaryField'],
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
            [ 'constant/transportProperties', 'nuVec\t(%g 0);\n' % nuVec[nameInd], name ],
            [ 'constant/transportProperties', 'alphaVec\t(%g 0);\n' % alphaVec[nameInd], name ],
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
        ['constant/reactiveProperties', 'activeReacts', '(1 0)', ''],
        ['constant/reactiveProperties', 'k0', '%g' % k0, 'reaction00'],
        ['constant/reactiveProperties', 'kin', '%g' % KInh , 'reaction00'],
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
        ['system/fvSolution', 'nConcCorrectors', '0', ''],
        ['system/fvSolution', 'nTempCorrectors', '0', ''],
        ['system/fvSolution', 'nTogCorrectors', '1', ''],
        ['system/fvSolution', 'nNonOrthogonalCorrectors', '0', ''],
    ]
)
baseCase.replace(
    [
        ['system/fvSolution', ['"(customSolver|COMass|O2Mass|CO2Mass|ArMass|)" 	 0.9;'], ['"(customSolver|COMass|O2Mass|CO2Mass|ArMass|)"  %g;' % relEqMass]],
        ['system/fvSolution', ['"(T)" 0.9995;'], ['"(T)" %g;' % relEqT]]
    ]
)

# 7) system/fvSchemes
baseCase.replace(
    [
        ['system/fvSchemes', ['corrected 0.5'], ['uncorrected']],
    ]
)

# others
baseCase.replace(
    [
        [ 'system/residuals', [ 'customSolver' ], [ 'U p T %s' %species.replace(' ','Mass\n') ] ], 
    ]
) 

# -- run the simulation
baseCase.runCommands(
    [
        'chmod 755 ./* -R',
        'blockMesh > log.blockMeshDict',
        'rm -rf 0',
        'cp -r 0.org 0',
        'paraFoam -touch',
        'reactingHetCatSimpleFoamM'
    ]
)
