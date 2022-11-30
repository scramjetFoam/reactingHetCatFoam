#!/usr/bin/python

#FILE DESCRIPTION=======================================================
# a function to generate blockMeshDict for simple ejector nozzle-diffuser
# geometry
#
# The output is a simple blockMeshDict to generate the mesh
#
# axial symmetry - wedge
# diffuser control points - spline edge


#LICENSE================================================================
#  blockMeshDictGen.py
#  
#  Copyright 2018-2019 Martin Isoz & Lucie Kubickova
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

#IMPORT BLOCK=======================================================
import os
import copy
import math
import re
import sys

# from auxiliarFuncs import *

blockMeshClass = 'blockMeshDictClassV5'
command = 'from ' + blockMeshClass + ' import *'
exec(command)

# functions to return the boundary strings
def retBoundString(bName,bType,bFaces):
    nameStr     = '\t' + bName + '\n\t{\n'
    typeStr     = '\t\ttype ' + bType + ';\n'
    faceStr0    = '\t\tfaces\n\t\t(\n'
    faceStrE    = '\t\t);\n\t}\n\n'
    outStr      = [nameStr,typeStr,faceStr0]
    for bFace in bFaces:
        outStr.append('\t\t\t(' + ' '.join(' %d'%fc for fc in bFace) + ')\n')
    outStr.append(faceStrE)
    return outStr
    
#===================================================================
#							EDITABLE
#===================================================================
def genBlockMeshDict(
    geomPars,
    geomOrig,
    cellSize,
    caseDir,
):
    #GET THE DATA=======================================================
    #-------------------------------------------------------------------
    # GEOMETRY DATA
    #-------------------------------------------------------------------
    
    # -- unpack the input parameters
    DIn, DOut, L, wAng   = geomPars
    x0,y0,z0 = geomOrig
    
    # auxiliar computation
    RIn = DIn*0.5
    ROut = DOut*0.5

    yMax = ROut
    dZ = math.tan(wAng*0.5/180*math.pi)*yMax

    #-------------------------------------------------------------------
    # MESH DATA
    #-------------------------------------------------------------------
    
    # size of a single cell
    dX,dY = cellSize
    
    # grading
    grX = "1.0"
    grYB = "1.0"
    # grYT = "1.00001"
    grYT = "1.00"
    grZ = "1.0"

    # number of cells
    nCX  = int(round(L/dX))
    nCYB = int(round(RIn/dY))
    # nCYT = int(round(2*(ROut-RIn)/dY/(float(grYT) - 1)))
    nCYT = int(round((ROut-RIn)/dY))

    # mesh scale
    mScale  = 1

    #===================================================================
    #							DO NOT EDIT
    #===================================================================
    
    #-------------------------------------------------------------------
    # GENERATE THE BLOCK OBJECTS
    #-------------------------------------------------------------------
    
    blockCount = 0
    wBlockCount= 0
    blockLst   = []
    edgeLst = []
    
    #-----------------------------------------------------------------------
    # bottom block
    xC,yC,zC = x0,y0,z0
    xE,yE = xC+L,yC+RIn

    bottomBlockPars = [
        blockCount*8+wBlockCount*6,                                         #initial index
        [                                                                   #list of vertex coordinates
            [xC,yC,zC-yC/yMax*dZ],
            [xE,yC,zC-yC/yMax*dZ],
            [xE,yE,zC-yE/yMax*dZ],
            [xC,yE,zC-yE/yMax*dZ],
            [xE,yE,zC+yE/yMax*dZ],
            [xC,yE,zC+yE/yMax*dZ],
        ],
        [nCX,nCYB,1],                                                  #number of cells in the block
        [grX,grYB,grZ],                                                      #block grading
        [[4,0],[5,1]],                                                      #doubled vertices (replaced & actual)
        "bottomBlock"
    ]
    bottomBlock = wedgeBlockClass(bottomBlockPars)
    
    blockCount += 0
    wBlockCount += 1
    blockLst.append(bottomBlock)
    
    #-----------------------------------------------------------------------
    # top block
    xC,yC,zC = xC,yE,z0
    xE,yE = xE,y0+ROut

    topBlockPars = [
        blockCount*8+wBlockCount*6,                                         #initial index
        [                                                                   #list of vertex coordinates
            [xC,yC,zC-yC/yMax*dZ],
            [xE,yC,zC-yC/yMax*dZ],
            [xE,yE,zC-yE/yMax*dZ],
            [xC,yE,zC-yE/yMax*dZ],
            [xC,yC,zC+yC/yMax*dZ],
            [xE,yC,zC+yC/yMax*dZ],
            [xE,yE,zC+yE/yMax*dZ],
            [xC,yE,zC+yE/yMax*dZ],
        ],
        [nCX,nCYT,1],                                                  #number of cells in the block
        [grX,grYT,grZ],                                                      #block grading
        "topBlock"
    ]
    topBlock = hexBlockClass(topBlockPars)
    
    blockCount += 1
    wBlockCount += 0
    blockLst.append(topBlock)

    #-----------------------------------------------------------------------
    # prepare boundaries

    boundStrLst = []
    
    # -- wedge - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    wedgeZ0Boundary = []
    for block in blockLst:
        wedgeZ0Boundary.append(block.retFXY0())
    
    wedgeZ0BoundaryStr = retBoundString('wedgeZ0','wedge',wedgeZ0Boundary)
    
    boundStrLst.append(wedgeZ0BoundaryStr)
    
    wedgeZEBoundary = []
    for block in blockLst:
        wedgeZEBoundary.append(block.retFXYE())
    
    wedgeZEBoundaryStr = retBoundString('wedgeZE','wedge',wedgeZEBoundary)
    
    boundStrLst.append(wedgeZEBoundaryStr)
    
    # -- walls - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    walls = []
    walls.append(bottomBlock.retFXZ0())
    walls.append(topBlock.retFYZ0())
    
    boundStrLst.append(
        retBoundString('walls','wall',walls)
    )

    # -- atmosphere - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    atmo = []
    atmo.append(topBlock.retFXZE())
    
    boundStrLst.append(
        retBoundString('atmosphere','wall',atmo)
    )

    # -- in/outs - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # - inlet boundary
    inlet = []
    inlet.append(bottomBlock.retFYZ0())
    
    inletStr = retBoundString('inlet','patch',inlet)
    
    boundStrLst.append(inletStr)
    
    # - outlet boundary
    outlet = []
    outlet.append(bottomBlock.retFYZE())
    outlet.append(topBlock.retFYZE())
    
    outletStr = retBoundString('outlet','patch',outlet)
    
    boundStrLst.append(outletStr)
    
    # -- horizontal interfaces - - - - - - - - - - - - - - - - - - - - - - - 
    # - master boundary
    master = []
    master.append(bottomBlock.retFXZE())
    
    boundStrLst.append(
        retBoundString('master','patch',master)
    )
    
    # - slave boundary
    slave = []
    slave.append(topBlock.retFXZ0())
    
    boundStrLst.append(
        retBoundString('slave','patch',slave)
    )
    
    mergePairs = [
        '(master slave)',
    ]
    
    #-------------------------------------------------------------------
    # FILE GENERATION
    #-------------------------------------------------------------------
    bMD = open(caseDir + 'system/blockMeshDict','w')		            #open file for writing
    
    #-------------------------------------------------------------------
    # write the headline
    bMD.write('/*--------------------------------*- C++ -*----------------------------------*\ \n')
    bMD.write('| ========                 |                                                 | \n')
    bMD.write('| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n')
    bMD.write('|  \\    /   O peration     | Version:  4.1                                   | \n')
    bMD.write('|   \\  /    A nd           | Web:      www.OpenFOAM.org                      | \n')
    bMD.write('|    \\/     M anipulation  |                                                 | \n')
    bMD.write('\*---------------------------------------------------------------------------*/ \n')
    
    # write file description
    bMD.write('FoamFile \n')
    bMD.write('{ \n \t version \t 2.0; \n \t format \t ascii; \n')
    bMD.write(' \t class \t\t dictionary; \n \t object \t blockMeshDict; \n} \n')
    bMD.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n\n')
     
    #-------------------------------------------------------------------
    # convert to metres
    bMD.write('convertToMeters \t' + repr(mScale) + '; \n\n')
    
    #-------------------------------------------------------------------
    # write vertices
    bMD.write('vertices \n( \n')
    k = 0
    for block in blockLst:
        for vert in block.vCoords:
            bMD.write('\t ( ' + ' '.join(str(e) for e in vert) + ' )\t//' + ' %03d'%k + '\n')
            k = k+1
    bMD.write('); \n\n')
    
    #-----------------------------------------------------------------------
    # write edges
    bMD.write('edges \n( \n')
    for edge in edgeLst:
        if isinstance(edge[0],str):
            splineType = edge[0]
            del edge[0]

        else:
            splineType = 'polyLine'

        bMD.write('\t '+splineType+' ' + str(edge[0]) + ' ' + str(edge[1]))

        del edge[0]
        del edge[0]

        if splineType == 'polyLine':
            bMD.write('  (\n')

        for point in edge:
            bMD.write('\t( ' + ' '.join(str(e) for e in point) + ' )\n')

        if splineType == 'polyLine':
            bMD.write(')\n')

    bMD.write('); \n\n')
    
    #-------------------------------------------------------------------
    # write blocks
    bMD.write('blocks \n( \n')
    
    for block in blockLst:
        for line in block.retBlockString():
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
    # close file
    bMD.close()

    return blockMeshClass
