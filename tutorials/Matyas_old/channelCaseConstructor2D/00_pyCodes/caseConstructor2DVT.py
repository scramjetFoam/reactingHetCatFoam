#!/usr/bin/python

#########DO NOT EDIT####################################################

# -- runtime monitoring ------------------------------------------------
from time import time
start_time = time()
# ----------------------------------------------------------------------

import os
import math
import shutil as sh
import math

#IMPORT BLOCK-CUSTOM====================================================
from blockMeshDictClassV5nicer import *
from topoSetDictClassV1 import *
from caseConstructorAuxFuncsVM import *
from fileHeaderFooter import *


#=======================================================================
#							EDITABLE
#=======================================================================

# -- dictionary with different media info===============================
# -- name = list with used porous zones
# -- perm = list with permeabilities of the media
# -- inPt = list of interpolation points files 
# -- notUse = list of lists which blocks not use (plugs)
# -- zone = list of lists, where the zones are active ([[inlet],[outlet]])
# -- set used mediaType as mediaName

mediaInfoDict = {
    # 'CF1':          {'name':['wall'],
    #                  'perm':[12.8612e-12]},
    # 'CF2':          {'name':['wall','coat'],
    #                  'perm':[4.76617e-12,4.76617e-12]},
    # 'CF3':          {'name':['wall','coat'],
    #                  'perm':[0.623756e-12,4.76617e-12]},
    'uncurved1layer': {'name':['wall','coat'],
                     'perm':[5.0e-13,2.7600003091200346e-13],
                     'inPt':['SpaciPt1'],
                     'notUse': [[[0,0],[4,4]],[[0,4],[4,0]]],
                     'zone':[[[2,0],[2,1],[0,2],[1,2],[2,2],[3,2],[4,2],[2,3],[3,2],[4,2],[2,4]],
                             [[1,0],[3,0],[0,1],[1,1],[3,1],[4,1],[0,3],[1,3],[3,3],[4,3],[1,4],[3,4]]],
    },
    # 'uncurved2layer': {'name':['wall','coat','lOncoat'],
    #                  'perm':[5.0e-13,2.7600003091200346e-15,2.7600003091200346e-15],
    #                  'inPt':['SpaciPt1','SpaciPt2'],
    #                  'notUse': [[[0,0],[6,6]],[[0,6],[6,0]]],
    #                  'zone':[[[3,0],[3,1],[3,2],[0,3],[1,3],[2,3],[3,3],[4,3],[5,3],[6,3],[3,4],[3,5],[3,6]],
    #                          [[2,0],[4,0],[2,1],[4,1],[0,2],[1,2],[2,2],[4,2],[5,2],[6,2],[0,4],[1,4],[2,4],[4,4],[5,4],[6,4],[2,4],[4,4],[2,5],[4,5],[2,6],[4,6]],
    #                          [[1,0],[5,0],[0,1],[1,1],[5,1],[6,1],[0,5],[1,5],[5,5],[6,5],[1,6],[5,6]]]
    #                 },
}

#-----------------------------------------------------------------------
# I/O DATA
#-----------------------------------------------------------------------

# -- media type 
mediaName = 'uncurved1layer'

# -- base openfoam folder for simulations
baseCase = '../10_baseCasePimple/'

# -- out openfoam folder
outFolder= '../testTomas/'

# -- name of the simulation folder
name = 'pimple_test_%s_2D_V5'%mediaName

#-----------------------------------------------------------------------
# CASE/RUN SETTINGS
#-----------------------------------------------------------------------
nCores      = 6                                                         # number of cores to run the case on
startTime   = 0                                                         # simulation startTime
endTime     = 2000                                                      # simulation endTime
wrInt       = 0.01                                                        # simulation write interval

#-----------------------------------------------------------------------
# GEOMETRY DATA
#-----------------------------------------------------------------------
# -- width dimensions
WWl = 137e-6                                                            # wall thickness
WCh = 1.13e-3*0.5                                                       # channel width (square)                                                            
nCh = 4                                                                 # number of channels

# --length dimensions
LBf1 = 5.0e-3   # 7.0e-3                                                # length of buffer in front of channels
LBf2 = 40.0e-3   #                                                       # length of buffer behind channels
LPl = 1.0e-3   # 1.0e-3                                                # plug length
LCh = 20.0e-3  # 76.2e-3, 30.2e-3                                      # channel length (total, including the plugs, but WITHOUT the buffers)

LCh = LCh - 2*LPl
print('Geometry info: WWl = %g, WCh = %g, LBf1 = %g, LBf2 = %g, LPl = %g, LCh = %g'%(WWl, WCh, LBf1, LBf2, LPl, LCh))

# -- list with length dimensions
LLst = [LBf1, LPl, LCh, LPl, LBf2]

#-----------------------------------------------------------------------
# MESH DATA
#-----------------------------------------------------------------------
dX, dY, dZ   = 30e-6, 30e-6, 30e-6

# -- multiplication factors for number of cells in coating and wall
nTimesCoat = 2  
nTimesWall = 2

#-----------------------------------------------------------------------
# CASE PARAMETERS
#-----------------------------------------------------------------------
spVel    = 1000000                                                         # space velocity

print('Mesh discretization: %s'%(str([dX, dY, dZ])))

#-----------------------------------------------------------------------
# POROSITY DATA
#-----------------------------------------------------------------------
porNames = mediaInfoDict[mediaName]['name']
perms = mediaInfoDict[mediaName]['perm']
print('Creating geometry with %d zones: %s, with permeabilities: %s' % (len(porNames), str(porNames), str(perms)))

dVals  = [1.0/val for val in perms]                                     # go from porosity do d coefficient
fVals  = [0.0 for val in dVals]                                         # f values (not used at the time)

#########PREFERABLY DO NOT EDIT#########################################
# -- create outFolder
os.system('mkdir %s' % outFolder)

# -- folder with source data
pyFolder = '../00_pyCodes/'                                             #source directory for python codes
pyList   = [                                                            #python files
    'caseConstructor2D',
    'blockMeshDictClassV5nicer',
    'topoSetDictClassV1',
    'caseConstructorAuxFuncsVM',
]


#########DO NOT EDIT####################################################

# -- calculate inlet velocity acording to space velocity
uIn = LCh/((1.0/spVel)*3600)

#-------------------------------------------------------------------
# DATA GATHERING
#-------------------------------------------------------------------
# -- copy the data from the baseCase
caseName = '%s/'%(name)
caseDir = outFolder + caseName

if os.path.isdir(caseDir):                                          #ensure, that the caseDir is clear
    sh.rmtree(caseDir)

sh.copytree(baseCase, caseDir)

# -- copy the data from the pyFolder
for pyCode in pyList:
    sh.copyfile(pyFolder + pyCode + '.py', caseDir + 'ZZ_genScripts/' + pyCode + '.py')

# -- copy stitchMeshSc.sh from the pyFolder
sh.copyfile(pyFolder + 'stitchMeshSc.sh', caseDir + 'stitchMeshSc.sh')
    
#-------------------------------------------------------------------
# BLOCKMESHDICT PREPARATION
#-------------------------------------------------------------------
# -- Load interpolation points
edgS = []
for edgsPts in mediaInfoDict[mediaName]['inPt']:
    edgS.append(importInterpPoints(edgsPts))
print('Loaded %d points files: %s'%(len(edgS), edgS))

# -- create width lists
# -- if there is coating, substract its latest verPoint, then add layers
if len(edgS) >= 1:       
    WChY = WCh - edgS[-1].verPoint2[1] + WWl/2
else:
    WChY = WCh

# -- old implementation 
WLst = WCh, WWl

# -- widths list 
# -- free channel
WLstX, WLstY = [dX], [WChY]
infoBLock = ['free']
# -- layers between layers
for wInd in range(len(edgS)-1):
    WLstY.append(edgS[-1-wInd].verPoint2[1]-edgS[-2-wInd].verPoint2[1])
    infoBLock.append(mediaInfoDict[mediaName]['name'][len(edgS)-wInd])
# -- last layer on the wall
try:
    WLstY.append(edgS[0].verPoint2[1]-WWl/2)
    infoBLock.append(mediaInfoDict[mediaName]['name'][1])
except:
    pass
# -- Wall
WLstY.append(WWl)
infoBLock.append(mediaInfoDict[mediaName]['name'][0])
# -- last layer on the wall
try:
    WLstY.append(edgS[0].verPoint2[1]-WWl/2)
    infoBLock.append(mediaInfoDict[mediaName]['name'][1])
except:
    pass
# -- layers between layers
for wInd in range(len(edgS)-1):
    WLstY.append(edgS[1+wInd].verPoint2[1]-edgS[wInd].verPoint2[1])
    infoBLock.append(mediaInfoDict[mediaName]['name'][wInd+2])
# -- free channel
WLstY.append(WChY)
infoBLock.append('free')

print('Created lists WLstX %s and WLstY %s'%(str(WLstX), str(WLstY)))
print('Infoblock:', infoBLock)

# -- number of blocks at X and Y direction
nBlocksX = 1
nBlocksY = len(WLstY)
nBlocksL = len(LLst)

# -- number of cells in X, Y & Z directions
nCXLst = [1]
        
nCYLst = []
for i in range(nBlocksY):   
    r = 1
    try:
        info = infoBLock[i%nBlocksY]
        if info == 'wall':   r = nTimesWall
        elif 'coat' in info: r = nTimesCoat
    except:
        pass
    nCYLst.append(WLstY[i]//dY*r)
    if nCYLst[-1] == 0:
        nCYLst[-1] = 1

nCZLst = [length//dZ for length in LLst]

print('Created disretization arrays: X, Y, Z', (nCXLst, nCYLst, nCZLst))

# mesh grading (basic)
grX, grY, grZ = 1, 1, 1

# mesh scale
mScale = 1

#-------------------------------------------------------------------
# -- recompute dX,dY,dZ to be of the actually used size
dXLst = dX
dYLst = [WLstY[ind]/float(nCYLst[ind]) for ind in range(nBlocksY)]
dZLst = [LLst[ind]/float(nCZLst[ind]) for ind in range(nBlocksL)]

# width of big channel
WBCh = WWl + WCh*2

#-------------------------------------------------------------------
# -- define which blocks are to be used
totBlocksX = range(nBlocksX)
totBlocksY = range(nBlocksY*nCh)  # range of total number of blocks in Ydir
totBlocksL = range(nBlocksL)       # range of total number of blocks in Zdir

blockUse = [[[True] for j in totBlocksY] for l in range(5)]  # True-filled template
for l in (1, 3):  # layers with plugs -> change plugs to False
    for j in totBlocksY:
        for k in totBlocksX:
            nBChx = k//nBlocksX                                         # number of channel in X direction
            nBChy = j//nBlocksY                                         # number of channel in Y direction
            if [k%nBlocksX, j%nBlocksY] in mediaInfoDict[mediaName]['notUse'][(nBChx + nBChy + (l-1)//2)%2]:
                blockUse[l][j][0] = False

zones = [[[None] for j in totBlocksY] for l in range(5)]  # None-filled template
for l in (1, 2, 3):  # layers with zones -> change zones to zoneName
    for zInd in range(len(mediaInfoDict[mediaName]['name'])):
        zoneName = mediaInfoDict[mediaName]['name'][zInd]
        for j in totBlocksY:
            for k in totBlocksX:
                if [k%nBlocksX,j%nBlocksY] in mediaInfoDict[mediaName]['zone'][zInd]:
                    zones[l][j][0] = zoneName

print('Prepared zones:', zones)

#-------------------------------------------------------------------
# -- generate the block objects
blockLst   = [[[[] for k in totBlocksX] for j in totBlocksY] for i in totBlocksL]  # blank list of lists of lists for block objects
x0, y0, z0 = 0.0, 0.0, 0.0
blockCount = 0

zC = z0
for i in totBlocksL:          # layers along Zdir
    yC = y0
    auxLstY = []
    for j in totBlocksY:      # Ydir
        xC = x0
        auxLstX = []
        for k in totBlocksX:  # Xdir
            if blockUse[i][j][k]:
                cZone = zones[i][j][k]
                lcLst = xC, yC, zC
                grLst = grX, grY, grZ
                dmLst = WLstX[k%nBlocksX], WLstY[j%nBlocksY], LLst[i]
                nCLst = nCXLst, nCYLst, nCZLst
                blockPars = funBlockPars(  # from caseConstructorAuxFuncs
                    blockCount,
                    lcLst,  # location list
                    dmLst,  # dimensions list
                    nCLst,  # nCells list
                    grLst,  # grading list
                    cZone,  # zone
                    i, j, k, 
                    WLst, WBCh, 
                    nBlocksX, nBlocksY)  # blockPars = [nPointsBefBlock, coordinates_list, nCells_list, grading_list, cZone, edges_list]
                blockLst[i][j][k] = hexBlockClass(blockPars) 
                blockCount += 1
            xC += WLstX[k%nBlocksX]
        yC += WLstY[j%nBlocksY]
    zC += LLst[i]

#-------------------------------------------------------------------
# -- prepare boundaries
boundStrLst = []

# - empty
emptyLeftBoundary = []
for i in totBlocksL:
    for j in totBlocksY:
        if isinstance(blockLst[i][j][-1], hexBlockClass):  # check if not plug
            emptyLeftBoundary.append(blockLst[i][j][-1].retFYZE())
emptyLeftBoundaryStr = retBoundString('emptyLeft', 'empty', emptyLeftBoundary)
boundStrLst.append(emptyLeftBoundaryStr)

emptyRightBoundary = []
for i in totBlocksL:
    for j in totBlocksY:
        if isinstance(blockLst[i][j][0], hexBlockClass):
            emptyRightBoundary.append(blockLst[i][j][0].retFYZ0())
emptyRightBoundaryStr = retBoundString('emptyRight', 'empty', emptyRightBoundary)
boundStrLst.append(emptyRightBoundaryStr)

# - symmetries
symTopBoundary = []
for i in totBlocksL:
    for k in totBlocksX:
        if isinstance(blockLst[i][-1][k], hexBlockClass): 
            symTopBoundary.append(blockLst[i][-1][k].retFXZE())    
symTopBoundaryStr = retBoundString('symmetryTop', 'symmetryPlane', symTopBoundary)
boundStrLst.append(symTopBoundaryStr)

symBottomBoundary = []
for i in totBlocksL:
    for k in totBlocksX:
        if isinstance(blockLst[i][0][k], hexBlockClass):
            symBottomBoundary.append(blockLst[i][0][k].retFXZ0())
symBottomBoundaryStr = retBoundString('symmetryBottom', 'symmetryPlane', symBottomBoundary)
boundStrLst.append(symBottomBoundaryStr)

# - wall boundary
# -- surrounding walls: for each used blocks checks for unused neighbours: defines wall
wallBoundary = []
for i in totBlocksL:
    for j in totBlocksY:
        for k in totBlocksX:
            if blockUse[i][j][k]:
                # if k > 0 and not blockUse[i][j][k-1]:
                #     wallBoundary.append(blockLst[i][j][k].retFYZ0())
                # if k < nBlocksX*nCh-1 and not blockUse[i][j][k+1]:
                #     wallBoundary.append(blockLst[i][j][k].retFYZE())
                if j > 0 and not blockUse[i][j-1][k]:
                    wallBoundary.append(blockLst[i][j][k].retFXZ0())
                if j < nBlocksY*nCh-1 and not blockUse[i][j+1][k]:
                    wallBoundary.append(blockLst[i][j][k].retFXZE())
                if i > 0 and not blockUse[i-1][j][k]:
                    wallBoundary.append(blockLst[i][j][k].retFXY0())
                if i < nBlocksL-1 and not blockUse[i+1][j][k]:
                    wallBoundary.append(blockLst[i][j][k].retFXYE())
wallBoundaryStr = retBoundString('walls', 'wall', wallBoundary)
boundStrLst.append(wallBoundaryStr)

# - inlet boundary
inletBoundary = []
i = 0
for j in totBlocksY:
    for k in totBlocksX:
        if blockUse[i][j][k]:
            inletBoundary.append(blockLst[i][j][k].retFXY0())
inletBoundaryStr = retBoundString('inlet', 'patch', inletBoundary)
boundStrLst.append(inletBoundaryStr)

# - outlet boundary
outletBoundary = []
i = -1
for j in totBlocksY:
    for k in totBlocksX:
        if blockUse[i][j][k]:
            outletBoundary.append(blockLst[i][j][k].retFXYE())
outletBoundaryStr = retBoundString('outlet', 'patch', outletBoundary)
boundStrLst.append(outletBoundaryStr)

# - internal faces to be removed
masterPerpBoundary = []
for i in range(nBlocksL-1):
    for j in totBlocksY:
        for k in totBlocksX:
            if blockUse[i][j][k] and blockUse[i+1][j][k]:
                masterPerpBoundary.append(blockLst[i][j][k].retFXYE())
masterPerpBoundaryStr = retBoundString('masterPerp', 'patch', masterPerpBoundary)
boundStrLst.append(masterPerpBoundaryStr)
        
slavePerpBoundary = []
for i in range(1,nBlocksL):
    for j in totBlocksY:
        for k in totBlocksX:
            if blockUse[i][j][k] and blockUse[i-1][j][k]:
                slavePerpBoundary.append(blockLst[i][j][k].retFXY0())
slavePerpBoundaryStr = retBoundString('slavePerp', 'patch', slavePerpBoundary)
boundStrLst.append(slavePerpBoundaryStr)

masterHorBoundary = []
slaveHorBoundary = []
masterVerBoundary = []
slaveVerBoundary = []
for i in totBlocksL:
    for j in totBlocksY:
        for k in totBlocksX:
            if blockUse[i][j][k]:
                # if k < nBlocksX*nChX-1 and blockUse[i][j][k+1]:
                #     masterVerBoundary.append(blockLst[i][j][k].retFYZE())
                #     slaveVerBoundary.append(blockLst[i][j][k+1].retFYZ0())
                if j < nBlocksY*nCh-1 and blockUse[i][j+1][k]:
                    masterHorBoundary.append(blockLst[i][j][k].retFXZE())
                    slaveHorBoundary.append(blockLst[i][j+1][k].retFXZ0())

masterHorBoundaryStr = retBoundString('masterHor', 'patch', masterHorBoundary)
boundStrLst.append(masterHorBoundaryStr)
slaveHorBoundaryStr = retBoundString('slaveHor', 'patch', slaveHorBoundary)
boundStrLst.append(slaveHorBoundaryStr)
masterVerBoundaryStr = retBoundString('masterVer', 'patch', masterVerBoundary)
boundStrLst.append(masterVerBoundaryStr)
slaveVerBoundaryStr = retBoundString('slaveVer', 'patch', slaveVerBoundary)
boundStrLst.append(slaveVerBoundaryStr)

mergePairs = []

stitchPairs = [
    ['masterPerp', 'slavePerp'],
    ['masterHor', 'slaveHor'],
    #['masterVer', 'slaveVer'],
]

writeStitching('stitchMeshSc.sh', stitchPairs)

#-------------------------------------------------------------------
# -- file generation
bMD = open(caseDir + 'system/blockMeshDict', 'w')		                # open file for writing

#-------------------------------------------------------------------
# write the file header
headStr = retFileHeaderStr()
for line in headStr:
    bMD.write(line)
    
#-------------------------------------------------------------------
# convert to metres
bMD.write('convertToMeters \t' + repr(mScale) + '; \n\n')

#-------------------------------------------------------------------
# write vertices
bMD.write('vertices \n( \n')
l = 0
for i in totBlocksL:
    for j in totBlocksY:
        for k in totBlocksX:
            if blockUse[i][j][k]:
                for vert in blockLst[i][j][k].vCoords:
                    bMD.write('\t (' + ' '.join(str(e) for e in vert) + ')\t//' + ' %03d'%l + '\n')
                    l += 1
bMD.write('); \n\n')

#-------------------------------------------------------------------
# write edges
bMD.write('edges \n( \n')
for i in totBlocksL:                                                                                              
    for j in totBlocksY:
        for k in totBlocksX:
            if blockUse[i][j][k]:     
                for line in blockLst[i][j][k].retPolyLine():		    # NOTE: .retPolyLine() or .retSpline()
                    bMD.write(str(line))
bMD.write('\n); \n\n')

#-------------------------------------------------------------------
# write blocks
bMD.write('blocks \n( \n')
for i in totBlocksL:
    for j in totBlocksY:
        for k in totBlocksX:
            if blockUse[i][j][k]:
                for line in blockLst[i][j][k].retBlockString():
                    bMD.write(line)
bMD.write('); \n\n')

#-------------------------------------------------------------------
# write boundaries
bMD.write('boundary \n( \n')
for boundStr in boundStrLst:
    for line in boundStr:
        bMD.write(line)
bMD.write('); \n\n')

#-------------------------------------------------------------------
# write patch pairs to be merged
bMD.write('mergePatchPairs \n( \n')
for line in mergePairs:
    bMD.write('\t' + line + '\n')
bMD.write('); \n\n')

#-------------------------------------------------------------------
# write the file footer
footStr = retFileFooterStr()
for line in footStr:
    bMD.write(line)

#-------------------------------------------------------------------
# close file
bMD.close()
   
#-------------------------------------------------------------------
# POROSITY PROPERTIES UPDATE
#-------------------------------------------------------------------
# -- open the file and write the header
pPD = open(caseDir +'constant/porosityProperties', 'w')

porClassMember = porosityPropsWriterClass()                             # prepare class member for easy writing

headStr = porClassMember.retFileHeaderStr()
for line in headStr:
    pPD.write(line)
    
# -- write the rest of the porosity classes
for ind in range(len(porNames)):
    currD = [dVals[ind] for dInd in range(3)]
    currF = [fVals[ind] for dInd in range(3)]
        
    currPorName = porNames[ind]
    currActive  = 'yes'
    currZoneName= porNames[ind]
    
    porProps    = [currPorName, currActive, currZoneName, currD, currF]
    
    porStr = porClassMember.retPorEntryStr(porProps)
    for line in porStr:
        pPD.write(line)
                
# -- write the file footer
footStr = porClassMember.retFileFooterStr()
for line in footStr:
    pPD.write(line)

#-------------------------------------------------------------------
# close file
pPD.close()
    
#-------------------------------------------------------------------
# CASE RUN PARAMETERS UPDATE
#-------------------------------------------------------------------

# decomposeParDict
idStr = ['numberOfSubdomains ']

pVals = [repr(nCores)]

# write everything to the file
with open(caseDir + './system/decomposeParDict', 'r') as file:
    # read a list of lines into data
    data = file.readlines()
    
for j in range(len(idStr)):
    for i in range(len(data)):
        fInd = data[i].find(idStr[j])
        if fInd>-1:
            data[i] = data[i][:fInd] + idStr[j] + '\t' + pVals[j] + ';\n'

with open(caseDir + './system/decomposeParDict', 'w') as file:
    file.writelines(data)

#-------------------------------------------------------------------
# controlDict
idStr = ['startTime ', 'endTime ', 'writeInterval ']
pVals = [repr(startTime), repr(endTime), repr(wrInt)]

# write everything to the file
with open(caseDir + './system/controlDict', 'r') as file:
    # read a list of lines into data
    data = file.readlines()

for j in range(len(idStr)):
    for i in range(len(data)):
        fInd = data[i].find(idStr[j])
        if fInd>-1:
            data[i] = data[i][:fInd] + idStr[j] + '\t' + pVals[j] + ';\n'

with open(caseDir + './system/controlDict', 'w') as file:
    file.writelines(data)

# -- boundary velocity update
with open(caseDir + '/0.org/U', 'r') as file:
    data = file.readlines()

# uIn_cyl = Vcap/60*1e-6 / (math.pi*radius[0]**2)
for ind in range(len(data)):
    if data[ind].find('UIN_intF') >= 0:#!!do not put comment on executed temperature
        data[ind] = 'internalField\tuniform\t(0 0 %g);\n'%(uIn)
    if data[ind].find('UIN_w') >= 0:#!!do not put comment on executed temperature
        data[ind] = '\t\t\tvalue\tuniform (0 0 %g);\n'%(uIn)
    # if data[ind].find('UIN_cyl') >= 0:#!!do not put comment on executed temperature
    #     data[ind] = '\t\t\tvalue\tuniform (0 0 %g);\n'%(uIn_cyl)

with open(caseDir + '/0.org/U', 'w') as file:
    file.writelines(data)

# -- runtime monitoring ------------------------------------------------
print('run %ss'%(time() - start_time))
# ----------------------------------------------------------------------
