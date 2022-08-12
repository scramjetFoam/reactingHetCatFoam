#!/usr/bin/python

#FILE DESCRIPTION=======================================================
# Python script to automatically generate cases for meso-scale model of
# flow in monolitic catalytic filter
#
# IO PARAMETERS
# 1. baseCase   - where are the input data saved (this is copied and
#                 modified 
# 2. outFolder  - where do I want the case to be created
# 3. caseName   - what is the name of the ouput folder
#
# CASE PARAMETERS
# 1. inlet velocity
# 2. channel dimensions
# 3. number of porosity "tiles"
# 4. length of porosity "tiles" - expressed in number of cells
#
# MESH PARAMETERS
# 1. dimensions of a single cell, dX,dY,dZ (these are slightly altered
#    during the script execution) 
#
# NOTE:
# - the blockMeshDict generation as well as topoSetDict generation are
#   directly included into the code. this is not very nice as any changes
#   in one of these scripts would enforce a new version of this
#   caseConstructor. however, these two may be separated in the future
#   to clean up the code. this file is mostly a proof of concept
#
# CHANGELOG:
# 2019-03-19 (TH)
# current state:
# - added coating layer with different permeability
# - fixed problem with more channel in hor. or vert. direction
# 2018-11-14 (MI)
# current state:
# - uses Tomas's curved channel codes
# - added automatic of.sh modifications (to run on altix)
# - this file is used to generate cases for pseudo-homogenous walls
# - included cycle to generate all the inlet velocities at once
# - included a possibility to load different set of interpolation points
#   -> breaks backwards compatibility !!!
# 2018-06-06 (MI)
# current state:
# - initial version of the script is complete and seems to produce viable
#   cases (I ran a few tests on my laptop - small and coarse cases)
# to do:
# - the on-top coating structure is different than assumed in the
#   topoSetDict generation (there is only a single layer of "tiles")
#   => the porosity zones generation needs to be changed accordingly
#   => the new porosity zones structure should be more stable than the
#      current one (and it should be more permeable through "cracks")
# - include automatic generation of "porosityProperties" file
#   => I need only one type of entry, head and footer. I will put this
#      in the caseConstructorAuxFuncs file
#   => the question is, do I want it to be a class or just three funcs
#   => NOTE: all the openFoam file footers are actually the same
#   => NOTE: in the header, I should put "location"

#LICENSE================================================================
#  caseConstructor.py
#
#  Copyright 2018 Martin Isoz & Tomas Hlavaty
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

#########DO NOT EDIT####################################################
import os
import math
import shutil as sh
import math

#IMPORT BLOCK-CUSTOM====================================================
from blockMeshDictClassV5 import *
from topoSetDictClassV1 import *
from caseConstructorAuxFuncsV4Spaci import *


#=======================================================================
#							EDITABLE
#=======================================================================

# -- dictionary with different media info===============================
# -- name = list with used porous zones
# -- perm = list with permeabilities of the media
# -- inPt = list of interpolation points files 
#       -- NOTE: this version of the blockMeshDict is so far prepared to deal with uncurved geometry
# -- notUse = list of lists which blocks not use (plugs)
# -- zone = list of lists, where the zones are active ([[inlet],[outlet]])
# -- set used mediaType as mediaName
mediaInfoDict = {
    'CF1':          {'name':['wall'],
                     'perm':[12.8612e-12]},
    'CF2':          {'name':['wall','coat'],
                     'perm':[4.76617e-12,4.76617e-12]},
    'CF3':          {'name':['wall','coat'],
                     'perm':[0.623756e-12,4.76617e-12]},
    'uncurved1layer': {'name':['wall','coat'],
                     'perm':[5.0e-13,2.7600003091200346e-15],
                     'inPt':['SpaciPt1'],
                     'notUse': [[[0,0],[4,4]],[[0,4],[4,0]]],
                     'zone':[[[2,0],[2,1],[0,2],[1,2],[2,2],[3,2],[4,2],[2,3],[3,2],[4,2],[2,4]],
                             [[1,0],[3,0],[0,1],[1,1],[3,1],[4,1],[0,3],[1,3],[3,3],[4,3],[1,4],[3,4]]],
                     },
    'uncurved2layer': {'name':['wall','coat','lOncoat'],
                     'perm':[5.0e-13,2.7600003091200346e-15,2.7600003091200346e-15],
                     'inPt':['SpaciPt1','SpaciPt2'],
                     'notUse': [[[0,0],[6,6]],[[0,6],[6,0]]],
                     'zone':[[[3,0],[3,1],[3,2],[0,3],[1,3],[2,3],[3,3],[4,3],[5,3],[6,3],[3,4],[3,5],[3,6]],
                             [[2,0],[4,0],[2,1],[4,1],[0,2],[1,2],[2,2],[4,2],[5,2],[6,2],[0,4],[1,4],[2,4],[4,4],[5,4],[6,4],[2,4],[4,4],[2,5],[4,5],[2,6],[4,6]],
                             [[1,0],[5,0],[0,1],[1,1],[5,1],[6,1],[0,5],[1,5],[5,5],[6,5],[1,6],[5,6]]]
                    },
}

#-----------------------------------------------------------------------
# I/O DATA
#-----------------------------------------------------------------------

# -- media type 
mediaName = 'uncurved1layer'

# -- base openfoam folder for simulations
baseCase = '../10_baseCaseVSpaciV4/'

# -- out openfoam folder
outFolder= '../testMatyas/'

# -- name of the simulation folder
name = 'test_%s_'%mediaName

#-----------------------------------------------------------------------
# CASE/RUN SETTINGS
#-----------------------------------------------------------------------
nCores      = 6                                                         #number of cores to run the case on
startTime   = 0                                                         #simulation startTime
endTime     = 2000                                                      #simulation endTime
wrInt       = 50                                                        #simulation write interval

#-----------------------------------------------------------------------
# GEOMETRY DATA
#-----------------------------------------------------------------------
# -- width dimensions
WWl = 137e-6                                                            #wall thickness
WCh = 1.13e-3*0.5                                                       #channel width (square)                                                            
nChX= 2  # 1                                                            #number of channels in hor. dir.
nChY= 2  # 1                                                            #number of channels in ver. dir.

# --length dimensions
LBf1 = 1.0e-3   # 7.0e-3                                                #length of buffer in front of channels
LBf2 = 5.0e-3   #                                                       #length of buffer behind channels
LPl  = 1.0e-3   # 1.0e-3                                                #plug length
LCh  = 10.0e-3  # 76.2e-3, 30.2e-3                                      #channel length (total, including the plugs, but WITHOUT the buffers)

LCh = LCh - 2*LPl
print('Geometry info: WWl = %g, WCh = %g, LBf1 = %g, LBf2 = %g, LPl = %g, LCh = %g'%(WWl,WCh,LBf1,LBf2,LPl,LCh))

# -- list with length dimensions
LLst = [LBf1,LPl,LCh,LPl,LBf2]

#-----------------------------------------------------------------------
# MESH DATA
#-----------------------------------------------------------------------
dX,dY,dZ   = 80e-6,80e-6,120e-6

# -- multiplication factors for number of cells in coat and wall
nTimesCoat = 2
nTimesWall = 2

#-----------------------------------------------------------------------
# CASE PARAMETERS
#-----------------------------------------------------------------------
spVel    = 4000                                                         # space velocity

print('Mesh discretization: %s'%(str([dX,dY,dZ])))

#-----------------------------------------------------------------------
# POROSITY DATA
#-----------------------------------------------------------------------
porNames = mediaInfoDict[mediaName]['name']
perms = mediaInfoDict[mediaName]['perm']
print('Creating geometry with %d zones: %s, with permeabilities: %s'%(len(porNames),str(porNames),str(perms)))

dVals  = [1.0/val for val in perms]                                     #go from porosity do d coefficient
fVals  = [0.0 for val in dVals]                                         # f values (not used at the time)

#########PREFERABLY DO NOT EDIT#########################################
# -- create outFolder
os.system('mkdir %s'%outFolder)

# -- folder with source data
pyFolder = '../00_pyCodes/'                                             #source directory for python codes
pyList   = [                                                            #python files
    'caseConstructorV4Spaci',
    'blockMeshDictClassV5',
    'topoSetDictClassV1',
    'caseConstructorAuxFuncsV4Spaci',
    # edgePtsFl,
]
# NOTE: these files are copied into folder caseName/ZZ_genScripts/ for
#       future reference
# NOTE: update these files if needed - in case of new code versions


#########DO NOT EDIT####################################################
# -- for all capilary places (axial)
# for fraction in fracLst:
#     # -- for all capilatry inlet flowrate
#     for Vcap in VcapLst:
        # -- calculate capillary inlet flowrate in m3/s

# -- calculate inlet velocity acording to space velocity
uIn = LCh/((1.0/spVel)*3600)
# edgS = importInterpPoints(edgePtsFl)

#-------------------------------------------------------------------
# AUXILIARY COMPUTATIONS
#-------------------------------------------------------------------

# -- compute estimates for turbulence variables -- not used at the moment
# dHChan = WCh*2                                                      #channel hydraulic diameter (including symmetry)
# ReFlow  = uIn*dHChan/nuG
# ITurbAm = 0.227*ReFlow**-0.100                                      #turbulence intensity over the pipe area, Russo,F. and Basse, N.T. (2016)
# k0Turb  = 1.5*(uIn*ITurbAm)**2.0
# eps0Turb= 0.1643*k0Turb**1.5/(0.07*dHChan)
# om0Turb = eps0Turb/k0Turb
#print ReFlow,ITurbAm,k0Turb,om0Turb

#-------------------------------------------------------------------
# DATA GATHERING
#-------------------------------------------------------------------
# -- copy the data from the baseCase
caseName = '%s/'%(name)
caseDir = outFolder + caseName

if os.path.isdir(caseDir):                                          #ensure, that the caseDir is clear
    sh.rmtree(caseDir)

sh.copytree(baseCase,caseDir)

# -- copy the data from the pyFolder
for pyCode in pyList:
    sh.copyfile(pyFolder + pyCode + '.py',caseDir + 'ZZ_genScripts/' + pyCode + '.py')
    
#-------------------------------------------------------------------
# BLOCKMESHDICT PREPARATION
#-------------------------------------------------------------------

# -- Load interpolation points
edgS = []
for edgsPts in mediaInfoDict[mediaName]['inPt']:
    edgS.append(importInterpPoints(edgsPts))
print('Loaded %d points files: %s'%(len(edgS),edgS))

# -- create width lists
# -- if there is coating, substract its latest verPoint, then add layers
if len(edgS) >= 1:
    WChX = WCh-edgS[-1].verPoint2[0]+WWl/2       
    WChY = WCh-edgS[-1].verPoint2[1]+WWl/2
else:
    WChX, WChY = WCh, WCh

# -- old implementation 
WLst = WCh, WWl

# -- list with widths
# -- free channel
WLstX, WLstY = [WChX], [WChY]
infoBLock = ['free']
# -- layers between layers
for wInd in range(len(edgS)-1):
    WLstX.append(edgS[-1-wInd].verPoint2[0]-edgS[-2-wInd].verPoint2[0])
    WLstY.append(edgS[-1-wInd].verPoint2[1]-edgS[-2-wInd].verPoint2[1])
    infoBLock.append(mediaInfoDict[mediaName]['name'][len(edgS)-wInd])
# -- last layer on the wall
try:
    WLstX.append(edgS[0].verPoint2[0]-WWl/2)
    WLstY.append(edgS[0].verPoint2[1]-WWl/2)
    infoBLock.append(mediaInfoDict[mediaName]['name'][1])
except:
    pass
# -- Wall
WLstX.append(WWl)
WLstY.append(WWl)
infoBLock.append(mediaInfoDict[mediaName]['name'][0])
# -- last layer on the wall
try:
    WLstX.append(edgS[0].verPoint2[0]-WWl/2)
    WLstY.append(edgS[0].verPoint2[1]-WWl/2)
    infoBLock.append(mediaInfoDict[mediaName]['name'][1])
except:
    pass
# -- layers between layers
for wInd in range(len(edgS)-1):
    WLstX.append(edgS[1+wInd].verPoint2[0]-edgS[wInd].verPoint2[0])
    WLstY.append(edgS[1+wInd].verPoint2[1]-edgS[wInd].verPoint2[1])
    infoBLock.append(mediaInfoDict[mediaName]['name'][wInd+2])
# -- free channel
WLstX.append(WChX)
WLstY.append(WChY)
infoBLock.append('free')

print('Created lists WLstX %s and WLstY %s'%(str(WLstX),str(WLstY)))
print('Infoblock: ',infoBLock)

# list of real width and height of channel without coating and wall
# WSCh = WCh + WWl/2 + edgS.verPoint2[0],WCh + WWl/2 - edgS.verPoint2[1]

#list of widths of each block in x direction
# WLstX = [WSCh[0],-edgS.verPoint2[0]-WWl/2,WWl,-edgS.verPoint2[0]-WWl/2,WSCh[0]]
# WLstY = [WSCh[1],edgS.verPoint2[1]-WWl/2,WWl,edgS.verPoint2[1]-WWl/2,WSCh[1]]

# pt1 = [3*WCh+2*WWl+radius[1]+-edgS[-1].verPoint2[0]-WWl/2+offset[0],3*WCh+2*WWl+radius[1]+edgS[-1].verPoint2[1]-WWl/2+offset[1],LBf+fraction*LCh]
# pt1 = [4*WCh+2*WWl,4*WCh+2*WWl,LBf+fraction*LCh]
# pt2 = [pt1[0],pt1[1],LBf*2+LCh+0.1]


# for axial coating distribution, not used at the moment
# ~ LLst = [LBf,LPl,LCh,LPl,LBf]
# LLst = [LBf,LPl]
# for ind in range(len(LChP)-1):
#     LLst.append(float(LChP[ind+1]-LChP[ind]))
# LLst.append(LPl)
# LLst.append(LBf)

# ~ LLst = [LBf,LPl,LCh,LPl,LBf]
# print(LLst)

# bChP = []
# bChP.append([True, False])
# for ind in range(len(sCLstIO)-1):
#     if not (sCLstIO[ind+1] == sCLstIO[ind]):
#         mid = (sCLstIO[ind+1]+sCLstIO[ind])/2
#         print(mid)
#         bIO = []
#         if mid >= bCI and mid <= eCI:
#             bIO.append(True)
#         else:
#             bIO.append(False)
#         if mid >= bCO and mid <= eCO:
#             bIO.append(True)
#         else:
#             bIO.append(False)
#         bChP.append(bIO)
# bChP.append([False,True])
# print (bChP)

#-------------------------------------------------------------------
# -- prepare the inputs (mesh)
# ~ # x direction discretization

# -- number of blocks at X and Y direction
nBlocksX = len(WLstX)
nBlocksY = len(WLstY)

# -- number of cells in X and Y direction
nCXLst = []
for i in range(len((WLstX))):
    try:
        if infoBLock[i%nBlocksX] == 'wall':
            nCXLst.append(int(round(WLstX[i]/dX))*nTimesWall)
        elif 'coat' in infoBLock[i%nBlocksX]:
            nCXLst.append(int(round(WLstX[i]/dX))*nTimesCoat)
        else:
            nCXLst.append(int(round(WLstX[i]/dX)))
    except:
        nCXLst.append(int(round(WLstX[i]/dX)))
    if nCXLst[-1] == 0:
        nCXLst[-1] = 1
    # if mediaInfoDict[mediaName] i%nBlocksX:
    #     nCXLst.append(int(round(WLstX[i]/dX))*nTimesCoat)
    # elif i%nBlocksX == 2:
    #     nCXLst.append(int(round(WLstX[i]/dX))*nTimesWall)
    # elif i%nBlocksX == 3:
    #     nCXLst.append(int(round(WLstX[i]/dX))*nTimesCoat)
    # else:
    #     nCXLst.append(int(round(WLstX[i]/dX)))
        
nCYLst = []
for i in range(len((WLstY))):
    try:
        if infoBLock[i%nBlocksX] == 'wall':
            nCYLst.append(int(round(WLstY[i]/dY))*nTimesWall)
        elif 'coat' in infoBLock[i%nBlocksX]:
            nCYLst.append(int(round(WLstY[i]/dY))*nTimesCoat)
        else:
            nCYLst.append(int(round(WLstY[i]/dY)))
    except:
        nCYLst.append(int(round(WLstY[i]/dY)))
    if nCYLst[-1] == 0:
        nCYLst[-1] = 1
    # if i%nBlocksY == 1:
    #     nCYLst.append(int(round(WLstY[i]/dY))*nTimesCoat )
    # elif i%nBlocksY == 2:
    #     nCYLst.append(int(round(WLstY[i]/dY))*nTimesWall )
    # elif i%nBlocksY == 3:
    #     nCYLst.append(int(round(WLstY[i]/dY))*nTimesCoat )
    # else:
    #     nCYLst.append(int(round(WLstY[i]/dY)))

# print(nCXLst)
# ~ # y direction discretization
# ~ nCYLst = [int(round(length/dY)) for length in WLstX]

# ~ # z direction discretization
nCZLst = [int(round(length/dZ)) for length in LLst]

print('Created disretization arrays: X, Y, Z',(nCXLst,nCYLst,nCZLst))



# mesh grading (basic)
grX, grY, grZ = 1, 1, 1

# mesh scale
mScale  = 1

#-------------------------------------------------------------------
# -- recompute dX,dY,dZ to be of the actually used size
dXLst = [WLstX[ind]/float(nCXLst[ind]) for ind in range(len(WLstX))]
dYLst = [WLstY[ind]/float(nCYLst[ind]) for ind in range(len(WLstY))]
dZLst = [LLst[ind]/float(nCZLst[ind]) for ind in range(len(LLst))]

# width of big channel
WBCh = WWl + WCh*2

# NOTE: out of these, I actually use only dZLst

#-------------------------------------------------------------------
# -- define which blocks are to be used
l = 0
boolBlockUse = []
boolBlockUse.append([[True for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])#use the whole first layer
l+=1
boolBlockUse.append([[True for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])#plugs
for j in range(nChY*len(WLstY)):
    for k in range (nChX*len(WLstX)):
        nBChx = int(k/nBlocksX)                                                    #number of channel in X direction
        nBChy = int(j/nBlocksY)                                                    #number of channel in Y direction
        if [k%nBlocksX,j%nBlocksY] in mediaInfoDict[mediaName]['notUse'][(nBChx + nBChy)%2]:
            boolBlockUse[l][j][k] = False
l+=1
boolBlockUse.append([[True for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])#channel
l+=1
boolBlockUse.append([[True for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])#plugs
for j in range(nChY*len(WLstY)):
    for k in range (nChX*len(WLstX)):
        nBChx = int(k/nBlocksX)                                                    #number of channel in X direction
        nBChy = int(j/nBlocksY)                                                    #number of channel in Y direction
        if [k%nBlocksX,j%nBlocksY] in mediaInfoDict[mediaName]['notUse'][(nBChx + nBChy +1)%2]:
            boolBlockUse[l][j][k] = False
l+=1
boolBlockUse.append([[True for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])#use the whole last layer
l+=1

zones = []
l=0
zones.append([[None for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])
l+=1
zones.append([[None for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])
for zInd in range(len(mediaInfoDict[mediaName]['name'])):
    zoneName = mediaInfoDict[mediaName]['name'][zInd]
    for j in range(nChY*len(WLstY)):
        for k in range (nChX*len(WLstX)):
            if [k%nBlocksX,j%nBlocksY] in mediaInfoDict[mediaName]['zone'][zInd]:
                zones[l][j][k] = zoneName
l+=1
zones.append([[None for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])
for zInd in range(len(mediaInfoDict[mediaName]['name'])):
    zoneName = mediaInfoDict[mediaName]['name'][zInd]
    for j in range(nChY*len(WLstY)):
        for k in range (nChX*len(WLstX)):
            if [k%nBlocksX,j%nBlocksY] in mediaInfoDict[mediaName]['zone'][zInd]:
                zones[l][j][k] = zoneName
l+=1
zones.append([[None for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])
for zInd in range(len(mediaInfoDict[mediaName]['name'])):
    zoneName = mediaInfoDict[mediaName]['name'][zInd]
    for j in range(nChY*len(WLstY)):
        for k in range (nChX*len(WLstX)):
            if [k%nBlocksX,j%nBlocksY] in mediaInfoDict[mediaName]['zone'][zInd]:
                zones[l][j][k] = zoneName
l+=1
zones.append([[None for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])
l+=1
    # zones.append(zonesHere)

print('Prepared zones: ', zones)

#-------------------------------------------------------------------
# -- generate the block objects
blockCount = 0
blockLst   = [[[[] for k in range(nChX*len(WLstX))] for j in range(nChY*len(WLstY))] for i in range(len(LLst))]
x0,y0,z0 = 0.0,0.0,0.0

zC = z0
for i in range(len(LLst)):                                          #layers along z
    yC = y0
    auxLstY = []
    for j in range(nChY*len(WLstY)):                                 #y-dir
        xC = x0
        auxLstX = []
        for k in range(nChX*len(WLstX)):                             #x-dir
            if boolBlockUse[i][j][k]:
                cZone = zones[i][j][k]
                LcLst = xC,yC,zC
                grLst = grX,grY,grZ
                dmLst = WLstX[k%nBlocksX], WLstY[j%nBlocksY], LLst[i]
                nCLst = nCXLst, nCYLst, nCZLst
                blockPars = funBlockPars(blockCount, LcLst, dmLst, nCLst, grLst, cZone,i, j, k,WLst,WBCh,nBlocksX,nBlocksY)      #imported function from caseConstructorAuxFuncs
                blockLst[i][j][k] = hexBlockClass(blockPars) 
                blockCount+=1
            xC+=WLstX[k%nBlocksX]
        yC+=WLstY[j%nBlocksY]
    zC+=LLst[i]

#-------------------------------------------------------------------
# -- prepare boundaries
boundStrLst = []

# - symmetries
symLeftBoundary = []
for i in range(len(LLst)):
    for j in range(nChY*len(WLstY)):
        if isinstance(blockLst[i][j][-1],hexBlockClass):
            symLeftBoundary.append(blockLst[i][j][-1].retFYZE())
        
symLeftBoundaryStr = retBoundString('symmetryLeft','symmetryPlane',symLeftBoundary)

boundStrLst.append(symLeftBoundaryStr)

symRightBoundary = []
for i in range(len(LLst)):
    for j in range(nChY*len(WLstY)):
        if isinstance(blockLst[i][j][0],hexBlockClass):             #I need to actually check - there might be void (plug block)
            symRightBoundary.append(blockLst[i][j][0].retFYZ0())

symRightBoundaryStr = retBoundString('symmetryRight','symmetryPlane',symRightBoundary)

boundStrLst.append(symRightBoundaryStr)

symTopBoundary = []
for i in range(len(LLst)):
    for k in range(len(WLstX)*nChX):
        if isinstance(blockLst[i][-1][k],hexBlockClass): 
            symTopBoundary.append(blockLst[i][-1][k].retFXZE())
        
symTopBoundaryStr = retBoundString('symmetryTop','symmetryPlane',symTopBoundary)

boundStrLst.append(symTopBoundaryStr)

symBottomBoundary = []
for i in range(len(LLst)):
    for k in range(len(WLstX)*nChX):
        if isinstance(blockLst[i][0][k],hexBlockClass):
            symBottomBoundary.append(blockLst[i][0][k].retFXZ0())
        
symBottomBoundaryStr = retBoundString('symmetryBottom','symmetryPlane',symBottomBoundary)

boundStrLst.append(symBottomBoundaryStr)

# - wall boundary
# -- surrounding walls
wallBoundary = []
for i in range(len(LLst)):
    for j in range(len(WLstY)*nChY):
        for k in range(len(WLstX)*nChX):
            if boolBlockUse[i][j][k]:
                if k > 0 and not boolBlockUse[i][j][k-1]:
                        wallBoundary.append(blockLst[i][j][k].retFYZ0())
                if k < len(WLstX)*nChX - 1:
                    if not boolBlockUse[i][j][k+1]:
                        wallBoundary.append(blockLst[i][j][k].retFYZE())
                if j > 0 and not boolBlockUse[i][j-1][k]:
                        wallBoundary.append(blockLst[i][j][k].retFXZ0())
                if j < len(WLstX)*nChY -1:
                    if not boolBlockUse[i][j+1][k]:
                        wallBoundary.append(blockLst[i][j][k].retFXZE())
                if i > 0 and not boolBlockUse[i-1][j][k]:
                    wallBoundary.append(blockLst[i][j][k].retFXY0())
                if i < len(LLst)-1:
                    if not boolBlockUse[i+1][j][k]:
                        wallBoundary.append(blockLst[i][j][k].retFXYE())

# NOTE: this is if hell, but it is self-explanatory

wallBoundaryStr = retBoundString('walls','wall',wallBoundary)

boundStrLst.append(wallBoundaryStr)

# - inlet boundary
inletBoundary = []
i = 0
for j in range(len(WLstY)*nChY):
    for k in range(len(WLstX)*nChX):
        if boolBlockUse[i][j][k]:
            inletBoundary.append(blockLst[i][j][k].retFXY0())

inletBoundaryStr = retBoundString('inlet','patch',inletBoundary)

boundStrLst.append(inletBoundaryStr)

# - outlet boundary
outletBoundary = []
i = -1
for j in range(len(WLstY)*nChY):
    for k in range(len(WLstX)*nChX):
        if boolBlockUse[i][j][k]:
            outletBoundary.append(blockLst[i][j][k].retFXYE())

outletBoundaryStr = retBoundString('outlet','patch',outletBoundary)

boundStrLst.append(outletBoundaryStr)

# - internal faces to be removed
masterPerpBoundary = []
for i in range(len(LLst)-1):
    for j in range(len(WLstY)*nChY):
        for k in range(len(WLstX)*nChX):
            if boolBlockUse[i][j][k] and boolBlockUse[i+1][j][k]:
                masterPerpBoundary.append(blockLst[i][j][k].retFXYE())
        
masterPerpBoundaryStr = retBoundString('masterPerp','patch',masterPerpBoundary)

boundStrLst.append(masterPerpBoundaryStr)
        
slavePerpBoundary = []
for i in range(1,len(LLst)):
    for j in range(len(WLstY)*nChY):
        for k in range(len(WLstX)*nChX):
            if boolBlockUse[i][j][k] and boolBlockUse[i-1][j][k]:
                slavePerpBoundary.append(blockLst[i][j][k].retFXY0())
        
slavePerpBoundaryStr = retBoundString('slavePerp','patch',slavePerpBoundary)

boundStrLst.append(slavePerpBoundaryStr)

masterHorBoundary = []
slaveHorBoundary = []
masterVerBoundary = []
slaveVerBoundary = []
for i in range(len(LLst)):
    for j in range(0,len(WLstY)*nChY):
        for k in range(0,len(WLstX)*nChX):
            if boolBlockUse[i][j][k]:
                if k < len(WLstX)*nChX-1:
                    if boolBlockUse[i][j][k+1]:
                        masterVerBoundary.append(blockLst[i][j][k].retFYZE())
                        slaveVerBoundary.append(blockLst[i][j][k+1].retFYZ0())
                if j < len(WLstY)*nChY-1:
                    if boolBlockUse[i][j+1][k]:
                        masterHorBoundary.append(blockLst[i][j][k].retFXZE())
                        slaveHorBoundary.append(blockLst[i][j+1][k].retFXZ0())

masterHorBoundaryStr = retBoundString('masterHor','patch',masterHorBoundary)

boundStrLst.append(masterHorBoundaryStr)

slaveHorBoundaryStr = retBoundString('slaveHor','patch',slaveHorBoundary)

boundStrLst.append(slaveHorBoundaryStr)

masterVerBoundaryStr = retBoundString('masterVer','patch',masterVerBoundary)

boundStrLst.append(masterVerBoundaryStr)

slaveVerBoundaryStr = retBoundString('slaveVer','patch',slaveVerBoundary)

boundStrLst.append(slaveVerBoundaryStr)

mergePairs = [
]

stitchPairs = [
    ['masterPerp','slavePerp'],
    ['masterHor','slaveHor'],
    ['masterVer','slaveVer'],
]

writeStitching('stitchMeshSc.sh',stitchPairs)

#-------------------------------------------------------------------
# -- file generation
bMD = open(caseDir+'system/blockMeshDict','w')		                #open file for writing

#-------------------------------------------------------------------
# write the file header
headStr = blockLst[0][0][0].retFileHeaderStr()
for line in headStr:
    bMD.write(line)
    
#-------------------------------------------------------------------
# convert to metres
bMD.write('convertToMeters \t' + repr(mScale) + '; \n\n')

#-------------------------------------------------------------------
# write vertices
bMD.write('vertices \n( \n')
l = 0
for i in range(len(LLst)):
    for j in range(nChY*len(WLstY)):
        for k in range(nChX*len(WLstX)):
            if boolBlockUse[i][j][k]:
                for vert in blockLst[i][j][k].vCoords:
                    bMD.write('\t ( ' + ' '.join(str(e) for e in vert) + ' )\t//' + ' %03d'%l + '\n')
                    l += 1
bMD.write('); \n\n')

#-------------------------------------------------------------------
# write edges
bMD.write('edges \n( \n')
for i in range(len(LLst)):                                                                                              
    for j in range(nChY*len(WLstY)):
        for k in range(nChX*len(WLstX)):
            if boolBlockUse[i][j][k]:     
                for line in blockLst[i][j][k].retPolyLine():		#NOTE: can be used .retPolyLine() or .retSpline()
                    bMD.write(str(line))
bMD.write('\n); \n\n')


#-------------------------------------------------------------------
# write blocks
bMD.write('blocks \n( \n')

for i in range(len(LLst)):
    for j in range(nChY*len(WLstY)):
        for k in range(nChX*len(WLstX)):
            if boolBlockUse[i][j][k]:
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
footStr = blockLst[0][2][2].retFileFooterStr()
for line in footStr:
    bMD.write(line)

#-------------------------------------------------------------------
# close file
bMD.close()

#-------------------------------------------------------------------
# INITIAL GUESS UPDATE
#-------------------------------------------------------------------

# # U
# pVals   = ['value           uniform (0 0 %f)'%uIn]  

# idStr   = ['value           uniform (0 0 UIN)']

# # write everything to the file
# with open(caseDir + './0.org/U', 'r') as file:
#     # read a list of lines into data
#     data = file.readlines()
    
# for j in range(len(idStr)):
#     for i in range(len(data)):
#         fInd = data[i].find(idStr[j])
#         if fInd>-1:
#             data[i] = data[i][:fInd] + pVals[j] + ';\n'

# with open(caseDir + './0.org/U', 'w') as file:
#     file.writelines( data )

# NOTE: nothing else is changed at the moment

#-------------------------------------------------------------------
# FLUID PROPERTIES UPDATE
#-------------------------------------------------------------------

# -- NOTE: maybe I can put here some transport properties changes
    
#-------------------------------------------------------------------
# POROSITY PROPERTIES UPDATE
#-------------------------------------------------------------------
# NOTE: the code structure here is slightly different than in orher parts
#       of the properties modification codes (I do not modify an existing
#       file, I generate a new one)
#       => this is close to the blockMeshDict and topoSetDict generation

# -- open the file and write the header
pPD = open(caseDir +'constant/porosityProperties','w')

porClassMember = porosityPropsWriterClass()                         #prepare class member for easy writing

headStr = porClassMember.retFileHeaderStr()
for line in headStr:
    pPD.write(line)
    
# -- write the properties of the base porosity class
# currD = [dVals[0] for dInd in range(3)]
# currF = [fVals[0] for dInd in range(3)]

# currPorName = basePorName
# currActive  = 'yes'
# currZoneName= basePorName

# porProps    = [currPorName,currActive,currZoneName,currD,currF]

# porStr = porClassMember.retPorEntryStr(porProps)
# for line in porStr:
#     pPD.write(line)
    
# -- write the rest of the porosity classes
for ind in range(len(porNames)):
    currD = [dVals[ind] for dInd in range(3)]
    currF = [fVals[ind] for dInd in range(3)]
        
    currPorName = porNames[ind]
    currActive  = 'yes'
    currZoneName= porNames[ind]
    
    porProps    = [currPorName,currActive,currZoneName,currD,currF]
    
    porStr = porClassMember.retPorEntryStr(porProps)
    for line in porStr:
        pPD.write(line)
                
# -- write the file footer
footStr = porClassMember.retFileFooterStr()
for line in footStr:
    pPD.write(line)
    
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
    file.writelines( data )
    
#-------------------------------------------------------------------
# controlDict
idStr = ['startTime ','endTime ','writeInterval ']

pVals = [repr(startTime),repr(endTime),repr(wrInt)]

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
    file.writelines( data )

# # -- snappyHexMesh for capillary
# with open(caseDir + '/system/snappyHexMeshDict', 'r') as file:
#     # read a list of lines into data
#     data = file.readlines()

# for ind in range(len(data)):
#     if data[ind].find('point1') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\tpoint1\t(%g %g %g);\n'%(pt1[0],pt1[1],pt1[2])
#     if data[ind].find('point2') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\tpoint2\t(%g %g %g);\n'%(pt2[0],pt2[1],pt2[2])
#     if data[ind].find('chInt1MinMoje') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\tmin (-0.1 -0.1 %g);	//chInt1MinMoje\n'%(pt1[2]-1e-3)
#     if data[ind].find('chInt1MaxMoje') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\tmax (1 1 %g); //chInt1MinMoje\n'%(pt1[2]+1e-3)
#     if data[ind].find('chInt2MinMoje') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\tmin (-0.1 -0.1 %g);	//chInt2MinMoje\n'%(pt1[2]-100e-6)
#     if data[ind].find('chInt2MaxMoje') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\tmax (1 1 %g); //chInt2MinMoje\n'%(pt1[2]+100e-6)
#     if data[ind].find('radius') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\tradius\t%g;\n'%(radius[1])

# with open(caseDir + '/system/snappyHexMeshDict', 'w') as file:
#     file.writelines(data)

# # -- topoSetDict for capillary
# with open(caseDir + '/system/topoSetDict', 'r') as file:
#     # read a list of lines into data
#     data = file.readlines()

# for ind in range(len(data)):
#     if data[ind].find('p1') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\t\tp1\t(%g %g %g);\n'%(pt1[0],pt1[1],-0.1)
#     if data[ind].find('p2') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\t\tp2\t(%g %g %g);\n'%(pt2[0],pt2[1],pt2[2])
#     if data[ind].find('outerRadius') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\t\touterRadius\t%g;\n'%(radius[1]+50e-6)
#     if data[ind].find('innerRadius') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\t\tinnerRadius\t%g;\n'%(radius[0])

# with open(caseDir + '/system/topoSetDict', 'w') as file:
#     file.writelines(data)


# -- boundary velocity update
with open(caseDir + '/0.org/U', 'r') as file:
    # read a list of lines into data
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


# # -- sampling at different locations
# with open(caseDir + '/system/sample', 'r') as file:
#     # read a list of lines into data
#     data = file.readlines()

# for ind in range(len(data)):
#     if data[ind].find('z_StartCap') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\t\tpoint\t(0.002185 0.002185 %g);\n'%(pt1[2])

# with open(caseDir + '/system/sample', 'w') as file:
#     file.writelines(data)

# with open(caseDir + '/system/topoSetDict.3', 'r') as file:
#     # read a list of lines into data
#     data = file.readlines()

# x = [pt2[0]-WWl-WCh*2,pt2[0]+WWl+WCh*2,pt2[0]+2*(WWl+WCh*2),pt2[0]-WWl-WCh*2]
# y = [pt2[1]-WWl-WCh*2,pt2[1],pt2[1],pt2[1]]

# ch1 = WCh + WWl/2
# ch2 = 3*WCh + 1.5*WWl

# for ind in range(len(data)):
#     if data[ind].find('x1') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\t\tp1 (%g %g 0.0);\n'%(x[0],y[0])
#         data[ind+1] = '\t\t\tp2 (%g %g 1.0);\n'%(x[0],y[0])
    
#     if data[ind].find('boxCH1') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\t\tbox   (%g %g 0) (%g %g 1);\n'%(ch1,ch1,ch2,ch2)
#         # data[ind+1] = '\t\t\tp2 (%g %g 1.0);\n'%(xch1,ych1)

#     if data[ind].find('x2') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\t\tp1 (%g %g 0.0);\n'%(x[1],y[1])
#         data[ind+1] = '\t\t\tp2 (%g %g 1.0);\n'%(x[1],y[1])
    
#     if data[ind].find('x3') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\t\tp1 (%g %g 0.0);\n'%(x[2],y[2])
#         data[ind+1] = '\t\t\tp2 (%g %g 1.0);\n'%(x[2],y[2])
    
#     if data[ind].find('x4') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\t\tp1 (%g %g 0.0);\n'%(x[3],y[3])
#         data[ind+1] = '\t\t\tp2 (%g %g 1.0);\n'%(x[3],y[3])

#     # if data[ind].find('radiusCH') >= 0:#!!do not put comment on executed temperature
#     #     data[ind] = '\t\t\tradius\t%g;\n'%(WCh)
    
#     if data[ind].find('radiusCAP') >= 0:#!!do not put comment on executed temperature
#         data[ind] = '\t\t\tradius\t%g;\n'%(radius[0])

# with open(caseDir + '/system/topoSetDict.3', 'w') as file:
#     file.writelines(data)

