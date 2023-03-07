from OF_caseClass import OpenFOAMCase
import numpy as np

baseCase = OpenFOAMCase()
baseCase.loadOFCaseFromBaseCase('../baseCases/baseCaseMass')
baseCase.changeOFCaseDir('../ZZ_cases/fromKrak/toMult/10_baseCaseVSpaciV7')
baseCase.copyBaseCase()

# -- chemical species
specieNames = np.array(["CO", "O2", "CO2", "Ar"])
species = ''
species2 = ''
for specieName in specieNames:
    species += specieName + ' '
    species2 += specieName + 'Mass '
molMass = np.array([ 28e-3, 32e-3, 44e-3, 40e-3])
sigmaVs = np.array([ 18.0, 16.3, 26.9, 16.2])
nuVec = np.array([-1, -0.5, 1, 0])
alphaVec = np.array([0, 0, 0, 0])
yIn = np.array([0.002, 0.001, 0, 1-0.002-0.001])
MgIn = np.sum(yIn * molMass)
wIn = yIn * molMass / MgIn

reacZonesProp = np.array(
    [
        ['porEps', '[0 0 0 0 0 0 0]', 0.4, 0.3],
        ['tort', '[0 0 0 0 0 0 0]', 20, 4],
        ['dP', '[0 1 0 0 0 0 0]', 4e-6, 1e-6],
        ['kappa', '[1 1 -3 -1 0 0 0]', 5, 6],
        ['coatFrac', '[0 0 0 0 0 0 0]', 0.0, 1.0],
    ]
)

reactZones = ['wall', 'reactionZone']

for chSpI in range(len(specieNames)):
    name = specieNames[chSpI]
    baseCase.runCommands(['cp 0.org/bsChemSp 0.org/%sMass' % name])
    baseCase.replace( [ [ "0.org/%sMass" % ( name ), [ 'wChSpSet', 'nameSet' ], [ '%.5g' % (wIn[chSpI]), str(name) ] ] ] )
    baseCase.addToDictionary( 
        [
            [ '0.org/%sMass' % name, '\n\t"(outlet|walls|inletCyl|cylinder)"\n\t{\n\t\ttype zeroGradient;\n\t}\n\n', 'boundaryField'],
            [ '0.org/%sMass' % name, '\n\t"(symmetry.*)"\n\t{\n\t\ttype symmetryPlane;\n\t}\n\n', 'boundaryField'],
        ]
    )
baseCase.runCommands(['rm 0.org/bsChemSp'])

baseCase.replace(
    [
        ['0.org/U', ['(inv 0 0.0)'] ,['(0 0 UIN_intF)']],
    ]
)
baseCase.addToDictionary(
    [
        ['0.org/U', '\n\tinletCyl\n\t{\n\t\ttype\tfixedValue;\n\t\tvalue\tuniform (0 0 UIN_cyl);\n\t}\n', 'boundaryField'],
        ['0.org/U', '\n\t"(symmetry).*"\n\t{\n\t\ttype\tsymmetryPlane;\n\t}\n', 'boundaryField'],
        ['0.org/p', '\n\t"(symmetry).*"\n\t{\n\t\ttype\tsymmetryPlane;\n\t}\n', 'boundaryField'],
        ['0.org/T', '\n\t"(symmetry).*"\n\t{\n\t\ttype\tsymmetryPlane;\n\t}\n', 'boundaryField'],
        ['0.org/T', '\n\t"(walls|inletCyl|cylinder)"\n\t{\n\t\ttype\tzeroGradient;\n\t}\n', 'boundaryField'],
        ['0.org/p', '\n\t"(walls|inletCyl|cylinder)"\n\t{\n\t\ttype\tzeroGradient;\n\t}\n', 'boundaryField'],
        ['0.org/U', '\n\t"cylinder"\n\t{\n\t\ttype\tnoSlip;\n\t}\n', 'boundaryField'],
    ]
)

baseCase.replace(
    [
        ['0.org/T', ['isoT'], [str(433)]]
    ]
)

baseCase.addToDictionary( 
    [
        [ 'constant/transportProperties', 'species (%s);\n' % species, ''],
    ]
)

baseCase.setParameters(
    [
        [ 'system/decomposeParDict', 'numberOfSubdomains', '16' , '' ], 
    ] 
)

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
        ['system/fvSolution', 'nConcCorrectors', str(0), ''],
        ['system/fvSolution', 'nTempCorrectors', str(0), ''],
        ['system/fvSolution', 'nTogCorrectors', str(2), ''],
        # ['system/fvSolution', 'nTempCorrectors', str(7), ''],
    ]
)

baseCase.replace([
    [ 'system/residuals', [ 'customSolver' ], [ 'U p T %s' % species2 ] ], 
]) 

for rowI in range(reacZonesProp.shape[0]):
    for colI in range(reacZonesProp.shape[1]-2):
        baseCase.addToDictionary(
            [
                ['constant/transportProperties', '%s %s %s %g;\n' % (reacZonesProp[rowI, 0], reacZonesProp[rowI, 0], reacZonesProp[rowI,1], float(reacZonesProp[rowI,2+colI])), reactZones[colI]]
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

baseCase.setParameters( 
    [
        [ 'constant/thermophysicalProperties', 'molWeight', '%.5g' % (MgIn*1000), '' ],
        [ 'constant/thermophysicalProperties', 'Cp', '520.6768', '' ],
    ] 
)

baseCase.replace(
    [
        ['system/fvSolution',['customSolver'],['customSolver|%s' %species.replace(' ','Mass|')]],
    ]
)

baseCase.setParameters(
    [
        ['constant/reactiveProperties', 'k0', str(5e18), 'reaction00'],
        ['constant/reactiveProperties', 'kin', str(280), 'reaction00'],
        ['constant/reactiveProperties', 'activeReacts', '(1 0)', ''],
        ['system/fvSchemes', 'default', 'bounded Gauss vanLeer01', 'divSchemes'],
        ['system/fvSchemes', 'default', 'cellLimited Gauss linear 0.5', 'gradSchemes'],
        ['system/controlDict', 'endTime', '200', ''],
        ['system/controlDict', 'writeInterval', '100', ''],
    ]
)
baseCase.addToDictionary(
    [
        ['system/fvSchemes', '\tdiv((phi*interpolate(Cp)),T)  bounded Gauss linearUpwind grad(T);\n', 'divSchemes']
    ]
)

baseCase.runCommands(
    [
        'cp ../10_baseCaseVSpaciV6/system/snappyHexMeshDict ./system',
        'cp ../10_baseCaseVSpaciV6/Allrun-parallel ./',
        'cp ../10_baseCaseVSpaciV6/stitchMeshSc.sh ./',
        'cp ../10_baseCaseVSpaciV6/system/sample ./system',
        'cp ../10_baseCaseVSpaciV6/system/topoSetDict ./system',
        'cp ../10_baseCaseVSpaciV6/system/topoSetDict.2 ./system',
        'cp ../10_baseCaseVSpaciV6/system/topoSetDict.3 ./system',
        'cp ../10_baseCaseVSpaciV6/system/createPatchDict ./system',
        'cp ../10_baseCaseVSpaciV6/system/createPatchDict.2 ./system',
        'chmod 775 ./* -R'
    ]
)

# not used
baseCase.replace(
    [
        ['constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(1),str(1),str(0)]],
    ]
)