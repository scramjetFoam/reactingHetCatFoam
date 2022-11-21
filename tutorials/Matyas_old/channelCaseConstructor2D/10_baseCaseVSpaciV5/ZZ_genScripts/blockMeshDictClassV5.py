#!/usr/bin/python

#FILE DESCRIPTION=======================================================

# Python class containing the functions to ease up generation of
# blockMeshDict files

#LICENSE================================================================
#  blockMeshDictClass.py
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

#IMPORT BLOCK===========================================================
import math

class blockClass:
    """ Base python class containing functions to generate blockMeshDict
        blocks"""
        
    def __init__(self,blockPars):
        self.v0     = blockPars[0]                                      #initial vertex
        self.vCoords= blockPars[1]                                      #vertex coordinates
        self.nCells = blockPars[2]                                      #number of cells for discretization
        self.grading= blockPars[3]                                      #grading
        self.cZone  = blockPars[4]                                      #does the block belong to any cellZone?
        self.edges	= blockPars[5]										#edges
        
    def retVIndices(self):
        """ returns indices of the block vertices"""
        raise NotImplementedError()
        
    # functions to return block faces
    def retFXY0(self):
        """ return vert. indices of face in XY plane with 
            lower Z coordinate"""
        v0 = self.v0
        ind=self.retVIndices()
        return [ind[0],ind[1],ind[2],ind[3]]
    
    def retFXYE(self):
        """ return vert. indices of face in XY plane with 
            higher Z coordinate"""
        v0 = self.v0
        ind=self.retVIndices()
        return [ind[4],ind[5],ind[6],ind[7]]
        
    def retFXZ0(self):
        """ return vert. indices of face in XZ plane with 
            lower Y coordinate"""
        v0 = self.v0
        ind=self.retVIndices()
        return [ind[0],ind[1],ind[5],ind[4]]
    
    def retFXZE(self):
        """ return vert. indices of face in XZ plane with 
            higher Y coordinate"""
        v0 = self.v0
        ind=self.retVIndices()
        return [ind[3],ind[2],ind[6],ind[7]]
        
    def retFYZ0(self):
        """ return vert. indices of face in YZ plane with 
            lower X coordinate"""
        v0 = self.v0
        ind=self.retVIndices()
        return [ind[0],ind[3],ind[7],ind[4]]
    
    def retFYZE(self):
        """ return vert. indices of face in YZ plane with 
            higher X coordinate"""
        v0 = self.v0
        ind=self.retVIndices()
        return [ind[1],ind[2],ind[6],ind[5]]
   
    # function to return splines in the block	
    def retSpline(self):
        """ function that returns list of defined edges"""
        verStr 	 	 = ''
        for i in range(len(self.edges)):
            verStr += '\tspline ' + str(self.edges[i][0]) + ' ' + str(self.edges[i][1]) + ' (\n'
            for points in self.edges[i][2]:
                verStr 	 += '\t\t ('
                for vert in points:
                    verStr += ' ' + str(vert)
                verStr += ')\n'
            verStr 		 += '\t\t)\n'
        return [verStr]
    
    #function to return polyLines in the block
    def retPolyLine(self):
        """ function that returns list of defined edges"""
        verStr 	 	 = ''
        for i in range(len(self.edges)):
            verStr += '\tpolyLine ' + str(self.edges[i][0]) + ' ' + str(self.edges[i][1]) + ' (\n'
            for points in self.edges[i][2]:
                verStr 	 += '\t\t ('
                for vert in points:
                    verStr += ' ' + str(vert)
                verStr += ')\n'
            verStr 		 += '\t\t)\n'
        return [verStr]
        # Note: this is rather high level - the actual edges are defined
        #       out of the blockClass.
        #   (-) more complicated to define edges
        #   (+) it is available to handle arbitrary edges
        # Note: in here, I am not sure if it wouldn't be better to return
        #       directly arc string or to write another function, which
        #       will return arc string using retEdges
        
    # functions to return the block for blockMeshDict
    def retBlockString(self):
        introStr    = '\thex\n'
        if self.cZone is not None:
            vertStr     = '\t\t(' + ' '.join(' %d'%vert for vert in self.retVIndices()) + ')\t' + self.cZone + '\n'
        else:
            vertStr     = '\t\t(' + ' '.join(' %d'%vert for vert in self.retVIndices()) + ')\n'
        discStr     = '\t\t(' + ' '.join(' %d'%nC for nC in self.nCells) + ')'
        grTpStr     = '\tsimpleGrading '
        gradStr     = '\t(' + ' '.join(' %f'%gr for gr in self.grading) + ')\n\n'
        return [introStr,vertStr,discStr,grTpStr,gradStr]
        
    def retFileHeaderStr(self):
        """ returns the string with file header"""
        retStr = []
        retStr.append('/*--------------------------------*- C++ -*----------------------------------*\ \n')
        retStr.append('| ========                 |                                                 | \n')
        retStr.append('| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n')
        retStr.append('|  \\    /   O peration     | Version:  4.1                                   | \n')
        retStr.append('|   \\  /    A nd           | Web:      www.OpenFOAM.org                      | \n')
        retStr.append('|    \\/     M anipulation  |                                                 | \n')
        retStr.append('\*---------------------------------------------------------------------------*/ \n')
        
        retStr.append('FoamFile \n')
        retStr.append('{ \n \t version \t 2.0; \n \t format \t ascii; \n')
        retStr.append(' \t class \t\t dictionary; \n \t object \t blockMeshDict; \n} \n')
        retStr.append('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n\n')
        
        return retStr
    
    def retFileFooterStr(self):
        """ returns the string with file footer"""
        retStr = []
        retStr.append('// ************************************************************************* //\n\n')
        return retStr
        
    

class hexBlockClass(blockClass):
    """ Python class containing functions to generate HEX blockMeshDict
        blocks"""
        
    def __init__(self,blockPars):
        self.v0     = blockPars[0]                                      #initial vertex
        self.vCoords= blockPars[1]                                      #vertex coordinates
        self.nCells = blockPars[2]                                      #number of cells for discretization
        self.grading= blockPars[3]                                      #grading
        self.cZone  = blockPars[4]                                      #does the block belong to any cellZone?
        self.edges 	= blockPars[5]
        
    def retVIndices(self):
        """ returns indices of the block vertices"""
        return [self.v0+ind for ind in range(8)]
        
class wedgeBlockClass(blockClass):
    """ Python class containing functions to generate blockMeshDict
        wedge blocks"""
        
    def __init__(self,blockPars):
        self.v0     = blockPars[0]                                      #initial vertex
        self.vCoords= blockPars[1]                                      #vertex coordinates
        self.nCells = blockPars[2]                                      #number of cells for discretization
        self.grading= blockPars[3]                                      #grading
        self.cZone  = blockPars[4]                                      #does the block belong to any cellZone?
        self.vDbl   = blockPars[5]                                      #doubled edges (two pairs)
        
    def retVIndices(self):
        """ returns indices of the block vertices"""
        virtInd = [self.v0+ind for ind in range(8)] 
        for i in range(len(virtInd)):
            for vInd in self.vDbl:
                if i == vInd[0]:
                    virtInd[i] = virtInd[vInd[1]]
                    for j in range(i+1,len(virtInd)):
                        virtInd[j] = virtInd[j]-1
        return virtInd                                                  #wedge has 6 vertices - doubled
