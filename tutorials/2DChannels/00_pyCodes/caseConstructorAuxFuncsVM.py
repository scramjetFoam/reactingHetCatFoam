#!/usr/bin/python

#FILE DESCRIPTION=======================================================
# File with auxiliary functions for caseConstructor script

import math

#============import interpolation points================================
#~ import edgS

#============Function to set location of points in channel==============

# -- import interpolation points
def importInterpPoints(edgePointsFile):
    """ ugly way of importing different interpolation points for
        different medias """

    global edgS                                                         # define global variable (kind of ugly)
    edgS = __import__(edgePointsFile)
    return edgS

def funSetPtsVer(parLst, nPM):
    for point in nPM:
        if 'RW' in point[1]:                                            # moving point right (x vertice of point - parameter of move in x axe)
            parLst[point[0]][0] -= -edgS.verPoint2[0]--edgS.inTogether[-1][0]
        if 'LW' in point[1]:
            parLst[point[0]][0] += -edgS.verPoint2[0]--edgS.inTogether[-1][0]
        if 'UW' in point[1]:
            parLst[point[0]][1] += edgS.verPoint2[1] - edgS.inTogether[0][1]
        if 'DW' in point[1]:
            parLst[point[0]][1] -= edgS.verPoint2[1] - edgS.inTogether[0][1]
    return parLst

#============Function to transfer point vertices to new coordinate======
#============system with center in the center of channel================
def funTransVerON(verIn,WLst):                                                      
    verOut=verIn-WLst[0]-WLst[1]/2
    return verOut
    
#============Function to transfer point vertices to old coordinate======
#=======system with center in the right down corner of channel==========    
def funTransVerNO(verIn,WLst):
    verOut = verIn+WLst[0]+WLst[1]/2
    return verOut
    
#============Function to rotate vertices of face of block===============
def funRotateRight(nPointsBefBlock,oldV,hMT):                           # hMT - how many times, oldV - old 4 verticles of the face
    newV = []
    for vert in range(len(oldV)):
        newV.append(nPointsBefBlock+(oldV[vert]+hMT)%4)                 # rotation
    return newV

#============Function to define if block is in coating zone=============
def funCoatZone(j,k):
    nBChx = k//5                                                        # number of channel in X direction
    nBChy = j//5                                                        # number of channel in Y direction
    # if  (nBChx + nBChy)%2==0:
    if (k%5==3 and j%5==0) or (k%5==4 and j%5==1) or (k%5==0 and j%5==3) or (k%5==1 and j%5==4) or (k%5==3 and j%5==1) or (k%5==1 and j%5==3) or (k%5 == 1 and j%5==0) or (k%5==0 and j%5==1) or (k%5==4 and j%5==3) or (k%5==3 and j%5==4) or (k%5==1 and j%5==1)or (k%5==3 and j%5==3):
        return True
    else:
        return False
    # else:
    #     if (k%5 == 1 and j%5==0) or (k%5==0 and j%5==1) or (k%5==4 and j%5==3) or (k%5==3 and j%5==4) or (k%5==1 and j%5==1)or (k%5==3 and j%5==3):
    #         return True
    #     else:
    #         return False
        
#=============Function to define which block to use=====================
def funBlockNotUsed(j,k,place):
    nBChx = k//5                                                        # number of channel in X direction
    nBChy = j//5                                                        # number of channel in Y direction
    if  (nBChx + nBChy)%2==0:
        if (((k%5==0 and j%5==0) or (k%5==4 and j%5==4) or (k%5==3 and j%5==4) or (k%5==3 and j%5==3) or (k%5==4 and j%5==3) or (k%5==1 and j%5==0) or (k%5==1 and j%5==1) or (k%5==0 and j%5==1)) and 'in' in place) or (((k%5 == 4 and j%5==0) or (k%5==0 and j%5==4) or (k%5==0 and j%5==3) or (k%5==1 and j%5==3) or (k%5==1 and j%5==4) or (k%5==4 and j%5==1) or (k%5==3 and j%5==0) or (k%5==3 and j%5==1)) and 'out' in place):
            return True
        else:
            return False
    else:
        if (((k%5==0 and j%5==0) or (k%5==4 and j%5==4) or (k%5==3 and j%5==4) or (k%5==3 and j%5==3) or (k%5==4 and j%5==3) or (k%5==1 and j%5==0) or (k%5==1 and j%5==1) or (k%5==0 and j%5==1)) and 'out' in place) or (((k%5 == 4 and j%5==0) or (k%5==0 and j%5==4) or (k%5==0 and j%5==3) or (k%5==1 and j%5==3) or (k%5==1 and j%5==4) or (k%5==4 and j%5==1) or (k%5==3 and j%5==0) or (k%5==3 and j%5==1)) and 'in' in place):
            return True
        else:
            return False
            
#===========Function to define which block is in wall zone==============
def funWallZone(j,k):
    if (k%5 == 2 and j%5==0) or (k%5==2 and j%5==1) or (k%5==0 and j%5==2) or (k%5==1 and j%5==2) or (k%5==2 and j%5==2) or (k%5==3 and j%5==2) or (k%5==4 and j%5==2) or (k%5==2 and j%5==3) or (k%5==2 and j%5==4):
        return True
    else:
        return False
            
#====Function to make list of edges in block, based on default block in=
#====right up corner of the channel===================================== 
def funEdgesRU(nPointsBefBlock,way,dmLst,lcLst, WLst,nBChLst,WBCh):     # variable way - rotation - 'D'-default, 'R' - right, 'L' - left, 'O' - opposite, nBChLst - numbers of big channel
    edgesLs=[]                                                          # list with edges in block
    wC, hC, lC = dmLst                                                  # dimension of the block (used just lC)
    xC,yC,zC = lcLst                                                    # location of point 0 in the block (used just zC)
    edg01F = []
    edg12F = []
    edg01B = []
    edg12B = []
    rotVert =[]
#===============points for front face==================================#
    if 'D' in way:
        rotVert = funRotateRight(nPointsBefBlock,[0,1,2,3],0)
        for i in range(len(edgS.inPB01)):
            edg01F.append([funTransVerNO(edgS.inPB01[i][0],WLst)+nBChLst[0]*WBCh,funTransVerNO(edgS.inPB01[i][1],WLst)+nBChLst[1]*WBCh,zC])
        for i in range(len(edgS.inPB12)):
            edg12F.append([funTransVerNO(edgS.inPB12[i][0],WLst)+nBChLst[0]*WBCh,funTransVerNO(edgS.inPB12[i][1],WLst)+nBChLst[1]*WBCh,zC])
    elif 'R' in way:
        rotVert = funRotateRight(nPointsBefBlock,[0,1,2,3],1)
        for i in range(len(edgS.inPB12)-1,-1,-1):
            edg01F.append([funTransVerNO(edgS.inPB12[i][0],WLst)+nBChLst[0]*WBCh,funTransVerNO(-float(edgS.inPB12[i][1]),WLst)+nBChLst[1]*WBCh,+zC])
        for i in range(len(edgS.inPB01)-1,-1,-1):
            edg12F.append([funTransVerNO(edgS.inPB01[i][0],WLst)+nBChLst[0]*WBCh,funTransVerNO(-float(edgS.inPB01[i][1]),WLst)+nBChLst[1]*WBCh,+zC])
    elif 'L' in way:
        rotVert = funRotateRight(nPointsBefBlock,[0,1,2,3],3)
        for i in range(len(edgS.inPB12)-1,-1,-1):
            edg01F.append([funTransVerNO(-float(edgS.inPB12[i][0]),WLst)+nBChLst[0]*WBCh,funTransVerNO(edgS.inPB12[i][1],WLst)+nBChLst[1]*WBCh,zC])
        for i in range(len(edgS.inPB01)-1,-1,-1):
            edg12F.append([funTransVerNO(-float(edgS.inPB01[i][0]),WLst)+nBChLst[0]*WBCh,funTransVerNO(edgS.inPB01[i][1],WLst)+nBChLst[1]*WBCh,zC])
    elif 'O' in way:
        rotVert = funRotateRight(nPointsBefBlock,[0,1,2,3],2)
        for i in range(len(edgS.inPB01)):
            edg01F.append([funTransVerNO(-edgS.inPB01[i][0],WLst)+nBChLst[0]*WBCh,funTransVerNO(-edgS.inPB01[i][1],WLst)+nBChLst[1]*WBCh,zC])
        for i in range(len(edgS.inPB12)):
            edg12F.append([funTransVerNO(-edgS.inPB12[i][0],WLst)+nBChLst[0]*WBCh,funTransVerNO(-edgS.inPB12[i][1],WLst)+nBChLst[1]*WBCh,zC])
#================points for back face==================================#
    for i in range(len(edg01F)):
        edg01B.append([edg01F[i][0],edg01F[i][1],float(edg01F[i][2])+lC])
    for i in range(len(edg12F)):
        edg12B.append([edg12F[i][0],edg12F[i][1],float(edg12F[i][2])+lC])
#================final list of edges===================================#
    edgesLs.append ([rotVert[0],rotVert[1],edg01F])
    edgesLs.append ([rotVert[1],rotVert[2],edg12F])
    edgesLs.append([edgesLs[0][0]+4,edgesLs[0][1]+4,edg01B])
    edgesLs.append([edgesLs[1][0]+4,edgesLs[1][1]+4,edg12B])
    return edgesLs

def funEdgesCU(nPointsBefBlock,way,dmLst, lcLst, WLst, points,nBChLst,WBCh):  # variable way - rotation - 'D'-default, 'R' - right, 'L' - left, 'O' - opposite, points - which edges to use
    edgesLs=[]                                                          # list with edges in block
    wC, hC, lC = dmLst                                                  # dimension of the block (used just lC)
    xC,yC,zC = lcLst                                                    # location of point 0 in the block (used just zC)
    edg03F = []
    edg12F = []
    edg03B = []
    edg12B = []
    rotVert =[]
#===============points for front face==================================#
    if 'D' in way:
        rotVert = funRotateRight(nPointsBefBlock,[0,1,2,3],0)
        for i in range(len(edgS.inPB12)):
            edg03F.append([funTransVerNO(edgS.inPB12[i][0],WLst)+nBChLst[0]*WBCh,funTransVerNO(edgS.inPB12[i][1],WLst)+nBChLst[1]*WBCh,zC])
        for i in range(len(edgS.inPB12)):
            edg12F.append([funTransVerNO(-float(edgS.inPB12[i][0]),WLst)+nBChLst[0]*WBCh,funTransVerNO(edgS.inPB12[i][1],WLst)+nBChLst[1]*WBCh,zC])
    elif 'L' in way:
        rotVert = funRotateRight(nPointsBefBlock,[0,1,2,3],3)
        for i in range(len(edgS.inPB01)-1,-1,-1):
            edg03F.append([funTransVerNO(-float(edgS.inPB01[i][0]),WLst)+nBChLst[0]*WBCh,funTransVerNO(edgS.inPB01[i][1],WLst)+nBChLst[1]*WBCh,zC])
        for i in range(len(edgS.inPB01)-1,-1,-1):
            edg12F.append([funTransVerNO(-edgS.inPB01[i][0],WLst)+nBChLst[0]*WBCh,funTransVerNO(-edgS.inPB01[i][1],WLst)+nBChLst[1]*WBCh,zC])
    elif 'R' in way:
        rotVert = funRotateRight(nPointsBefBlock,[0,1,2,3],1)
        for i in range(len(edgS.inPB01)-1,-1,-1):
            edg12F.append([funTransVerNO(edgS.inPB01[i][0],WLst)+nBChLst[0]*WBCh,funTransVerNO(edgS.inPB01[i][1],WLst)+nBChLst[1]*WBCh,zC])
        for i in range(len(edgS.inPB01)-1,-1,-1):
            edg03F.append([funTransVerNO(edgS.inPB01[i][0],WLst)+nBChLst[0]*WBCh,funTransVerNO(-float(edgS.inPB01[i][1]),WLst)+nBChLst[1]*WBCh,+zC])
    elif 'O' in way:
        rotVert = funRotateRight(nPointsBefBlock,[0,1,2,3],2)
        for i in range(len(edgS.inPB12)):
            edg12F.append([funTransVerNO(edgS.inPB12[i][0],WLst)+nBChLst[0]*WBCh,funTransVerNO(-float(edgS.inPB12[i][1]),WLst)+nBChLst[1]*WBCh,+zC])
        for i in range(len(edgS.inPB12)):
            edg03F.append([funTransVerNO(-edgS.inPB12[i][0],WLst)+nBChLst[0]*WBCh,funTransVerNO(-edgS.inPB12[i][1],WLst)+nBChLst[1]*WBCh,zC])
#================points for back face==================================#
    for i in range(len(edg03F)):
        edg03B.append([edg03F[i][0],edg03F[i][1],float(edg03F[i][2])+lC])
    for i in range(len(edg12F)):
        edg12B.append([edg12F[i][0],edg12F[i][1],float(edg12F[i][2])+lC])
#================final list of edges===================================#
    if '03' in points:
        edgesLs.append ([rotVert[0],rotVert[3],edg03F])
        edgesLs.append([edgesLs[0][0]+4,edgesLs[0][1]+4,edg03B])
    elif '12' in points:
        edgesLs.append ([rotVert[1],rotVert[2],edg12F])
        edgesLs.append([edgesLs[0][0]+4,edgesLs[0][1]+4,edg12B])
    else:
        edgesLs.append ([rotVert[0],rotVert[3],edg03F])
        edgesLs.append ([rotVert[1],rotVert[2],edg12F])
        edgesLs.append([edgesLs[0][0]+4,edgesLs[0][1]+4,edg03B])
        edgesLs.append([edgesLs[1][0]+4,edgesLs[1][1]+4,edg12B])
    return edgesLs

def funBlockPars(blockCount, lcLst, dmLst, nCLst, grLst, cZone, i, j, k, WLst, WBCh, nBlocksX, nBlocksY):  # parameters named like in the whole code
    #============Declaration of variables===============================
    nPointsBefBlock = 8 * blockCount                                    # number of points used before block
    xC, yC, zC = lcLst                                                  # (lc -> location)
    wC, hC, lC = dmLst                                                  # (dm -> dimension)
    nCXLst, nCYLst, nCZLst = nCLst
    grX, grY, grZ = grLst                                               # (gr -> grading)
    nPM = []                                                            # list with moving of points [[number of point, how to move(R -right,L-left, D -down, U -up],[...]]
    edgesLs = []                                                        # list with edges (polylines, splines, arcs) in current block
    nBChx = k//5                                                        # number of big channelX
    nBChy = j//5                                                        # number of big channelY
    nBChLst = [nBChx,nBChy]                                             # list of numbers of big channel
    parLst = [[0,0] for nV in range(8)]                                 # list with parameters of moving points in the block [[x0,y0],[x1,y1],...,[x7,y7]]   

    #=====Setting nPM variable for funSetPtsVer=========================
    #=====Using funMakeBlockRD in py====================================
    #even blocks
    # if (nBChx + nBChy)%2==0:
    # if j%5 == 0 and k%5 == 3:
    #     nPM = [[1,'RW'],[5,'RW']]
    #     edgesLs = funEdgesCU(nPointsBefBlock,'O',dmLst, lcLst,WLst,'03',nBChLst,WBCh)
    # elif j%5 == 0 and k%5 == 4:
    #     nPM = [[0,'RW'],[4,'RW'],[2,'UW'],[6,'UW']]
    #     edgesLs = funEdgesRU(nPointsBefBlock,'O',dmLst, lcLst,WLst,nBChLst,WBCh)
    # elif j%5 == 1 and k%5 == 4:
    #     nPM = [[1,'UW'],[5,'UW']]
    #     edgesLs = funEdgesCU(nPointsBefBlock,'L',dmLst, lcLst,WLst,'12',nBChLst,WBCh)
    # elif j%5 == 3 and k%5 == 0:
    #     nPM = [[3,'DW'],[7,'DW']]
    #     edgesLs = funEdgesCU(nPointsBefBlock,'R',dmLst, lcLst,WLst,'12',nBChLst,WBCh)
    # elif j%5 == 4 and k%5 == 0:
    #     nPM = [[0,'DW'],[4,'DW'],[2,'LW'],[6,'LW']]
    #     edgesLs = funEdgesRU(nPointsBefBlock,'D',dmLst, lcLst,WLst,nBChLst,WBCh)
    # elif j%5 == 4 and k%5 == 1:
    #     nPM = [[3,'LW'],[7,'LW']]
    #     edgesLs = funEdgesCU(nPointsBefBlock,'D',dmLst, lcLst,WLst,'03',nBChLst,WBCh)
    #odd blocks
    # else:
    # elif j%5 == 0 and k%5==0:
    #     nPM = [[1,'LW'],[5,'LW'],[3,'UW'],[7,'UW']]
    #     edgesLs = funEdgesRU(nPointsBefBlock,'R',dmLst, lcLst,WLst,nBChLst,WBCh)
    # elif j%5 == 0 and k%5==1:
    #     nPM = [[0,'LW'],[4,'LW']]
    #     edgesLs = funEdgesCU(nPointsBefBlock,'O',dmLst, lcLst,WLst,'12',nBChLst,WBCh)
    # elif j%5==1 and k%5==0:
    #     nPM = [[0,'UW'],[4,'UW']]
    #     edgesLs = funEdgesCU(nPointsBefBlock,'R',dmLst, lcLst,WLst,'03',nBChLst,WBCh)
    # elif j%5==3 and k%5==4:
    #     nPM = [[2,'DW'],[6,'DW']]
    #     edgesLs = funEdgesCU(nPointsBefBlock,'L',dmLst, lcLst,WLst,'03',nBChLst,WBCh)
    # elif j%5==4 and k%5==3:
    #     nPM = [[2,'RW'],[6,'RW']]
    #     edgesLs = funEdgesCU(nPointsBefBlock,'D',dmLst, lcLst,WLst,'12',nBChLst,WBCh)
    # elif j%5==4 and k%5==4:
    #     nPM = [[1,'DW'],[5,'DW'],[3,'RW'],[7,'RW']]
    #     edgesLs = funEdgesRU(nPointsBefBlock,'L',dmLst, lcLst,WLst,nBChLst,WBCh)
    
    #==============Setting verticles of points in nPM ==================
    parLst = funSetPtsVer(parLst, nPM)
    
    #===blockPars with implemented parameters of trasfer and list for edges
    blockPars = [                                                      
        nPointsBefBlock,
        [
            [xC +      parLst[0][0],yC +      parLst[0][1],zC],
            [xC + wC + parLst[1][0],yC +      parLst[1][1],zC],
            [xC + wC + parLst[2][0],yC + hC + parLst[2][1],zC],
            [xC +      parLst[3][0],yC + hC + parLst[3][1],zC],
            [xC +      parLst[4][0],yC +      parLst[4][1],zC+lC],
            [xC + wC + parLst[5][0],yC +      parLst[5][1],zC+lC],
            [xC + wC + parLst[6][0],yC + hC + parLst[6][1],zC+lC],
            [xC +      parLst[7][0],yC + hC + parLst[7][1],zC+lC],
        ],
        [nCXLst[k%nBlocksX],nCYLst[j%nBlocksY],nCZLst[i]],
        [grX,grY,grZ],
        cZone,
        edgesLs,                                                        # edges (polyLine, spline, arcs) list [[numPtFrom,numPtTo,[list with interpolation points]],...]                                           
    ]
    return blockPars 

# -- blockMeshDictGen functions-----------------------------------------
def retBoundString(b_name, b_type, b_faces):
    name_str = '\t' + b_name + '\n\t{\n'
    type_str = '\t\ttype ' + b_type + ';\n'
    face_str0 = '\t\tfaces\n\t\t(\n'
    face_strE = '\t\t);\n\t}\n\n'
    out_str = [name_str, type_str, face_str0]
    for face in b_faces:
        out_str.append('\t\t\t(' + ' '.join('%d'%fc for fc in face) + ')\n')
    out_str.append(face_strE)
    return out_str
    
def writeStitching(fileName, stitchPairs):
    """ function to automatically write script for mesh stitching """
    stitchSc = open(fileName,'w')                                       # open file for writing
    
    stitchSc.write('#!/bin/sh\n\n')
    
    k = 0
    for pair in stitchPairs:
        stitchSc.write('stitchMesh -perfect -overwrite ' + pair[0] + ' ' + pair[1] + ' >> log.stitchMesh_%d\n'%k)
        k +=1
        
    stitchSc.close()


class porosityPropsWriterClass:
    """ Python class to return strings to write porosity properties file"""
    
    def retFileHeaderStr(self):
        """ returns the string with file header"""
        retStr = []
        retStr.append('/*--------------------------------*- C++ -*----------------------------------*\ \n')
        retStr.append('  =========                 |                                                   \n')
        retStr.append('  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox             \n')
        retStr.append('   \\\\    /   O peration     | Website:  https://openfoam.org                    \n')
        retStr.append('    \\\\  /    A nd           | Version:  10                                      \n')
        retStr.append('     \\\\/     M anipulation  |                                                   \n')
        retStr.append('\*---------------------------------------------------------------------------*/ \n')
        
        retStr.append('FoamFile \n')
        retStr.append('{ \n \t version \t 2.0; \n \t format \t ascii; \n')
        retStr.append(' \t class \t\t dictionary; \n \t location \t "constant"; \n \t object \t porosityProperties; \n} \n')
        retStr.append('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n\n')
        
        return retStr
    
    def retFileFooterStr(self):
        """ returns the string with file footer"""
        retStr = []
        retStr.append('// ************************************************************************* //\n\n')
        return retStr
        
    def retPorEntryStr(self, porProps):
        """ returns string with porosity entry according to the porosity
            properties"""
            
        # NOTE: at the time, always the same coordinate system is used
        # NOTE: at the time, we always used DarcyForchheimer porosity
        porName,active,zoneName,dVec,fVec = porProps
        
        retStr = []
        retStr.append('%s\n{\n'%porName)
        retStr.append('\ttype\t\tDarcyForchheimer;\n')
        retStr.append('\tactive\t\t%s;\n'%active)
        retStr.append('\tcellZone\t%s;\n\n'%zoneName)
        
        retStr.append('\tDarcyForchheimerCoeffs\n\t{\n')
        retStr.append('\t\td\t(%e %e %e);\n'%tuple(dVec))
        retStr.append('\t\tdPerm\tdPerm\t[0 -2 0 0 0 0 0]\t(%e %e %e);\n'%tuple(dVec))
        retStr.append('\t\tf\t(%e %e %e);\n\n'%tuple(fVec))
        
        retStr.append('\t\tcoordinateSystem\n\t\t{\n')
        retStr.append('\t\t\ttype\t\tcartesian;\n')
        retStr.append('\t\t\torigin\t\t(0 0 0);\n')
        retStr.append('\t\t\tcoordinateRotation\n\t\t\t{\n')
        retStr.append('\t\t\t\ttype\taxesRotation;\n')
        retStr.append('\t\t\t\te1\t\t(1 0 0);\n')
        retStr.append('\t\t\t\te2\t\t(0 0 1);\n')
        retStr.append('\t\t\t}\n')
        retStr.append('\t\t}\n')
        retStr.append('\t}\n')
        retStr.append('}\n')
        
        return retStr
        

