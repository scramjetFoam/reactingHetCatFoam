# -- using OF_caseClass for the description see OF_caseClass.py
# -- imports 

from OpenFoamData import OpenFoamData

from OF_caseClass import OpenFOAMCase
import sys
import numpy as np
import os

caseOutFolder = '../ZZ_cases/testAda'
numpyOutDir = '%s/ZZ_testNumpy' % caseOutFolder
fieldsOfInt = ['U', 'p']

testCase = OpenFOAMCase()
testCase.loadOFCaseFromBaseCase('/opt/openfoam10/tutorials/incompressible/icoFoam/cavity/cavity')
testCase.changeOFCaseDir(caseOutFolder)
testCase.copyBaseCase()
testCase.setParameters(
    [
        ['system/controlDict', 'endTime', str(3), ''],
        ['system/fvSchemes', 'default', 'Gauss SFCD', 'divSchemes'],
    ]
)
testCase.replace(
    [
        ['constant/physicalProperties', ['0.01'], ['0.02']]
    ]
)
testCase.addToDictionary(
    [
        ['system/fvSchemes','div(phi,U) Gauss upwind phi;\n', 'divSchemes']
    ]
)
testCase.runCommands(
    [
        'blockMesh > log.blockMesh',
        'icoFoam > log.icoFoam',
    ]
)

oFData = OpenFoamData(caseOutFolder, 1, 3, 'case', fieldsOfInt, (numpyOutDir))
oFData.loadTimeLst()
oFData.loadYsFromOFData()
oFData.saveYsFromNPField()

# oFData.loadYsFromNPField()

oFData.calcAvgY()

for i in range(len(fieldsOfInt)):
    oFData.writeField(
        oFData.avgs[i],
        fieldsOfInt[i],
        'avg_%s'%(fieldsOfInt[i]),
        caseOutFolder,
        outDir='%s/1000/' % (caseOutFolder)
    )