#!/usr/bin/python

#FILE DESCRIPTION=======================================================

# Python class containing the functions to ease up generation of
# blockMeshDict files

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
        self.name   = blockPars[4]
        
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
        
    # functions to return the block for blockMeshDict
    def retBlockString(self):
        introStr    = '\thex\n'
        vertStr     = '\t\t(' + ' '.join(' %d'%vert for vert in self.retVIndices()) + ')'
        if not self.name == None:
            nameStr     = '\t' + self.name
        discStr     = '\t(' + ' '.join(' %d'%nC for nC in self.nCells) + ')'
        grTpStr     = '\tsimpleGrading '
        gradStr     = '\t(' + ' '.join(' %s'%gr for gr in self.grading) + ')\n\n'
        if not self.name == None:
            return [introStr,vertStr,nameStr,discStr,grTpStr,gradStr]
        else:
            return [introStr,vertStr,discStr,grTpStr,gradStr]
        
class hexBlockClass(blockClass):
    """ Python class containing functions to generate HEX blockMeshDict
        blocks"""
        
    def __init__(self,blockPars):
        self.v0     = blockPars[0]                                      #initial vertex
        self.vCoords= blockPars[1]                                      #vertex coordinates
        self.nCells = blockPars[2]                                      #number of cells for discretization
        self.grading= blockPars[3]                                      #grading
        self.name   = blockPars[4]
        
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
        self.vDbl   = blockPars[4]                                      #doubled edges (two pairs)
        self.name   = blockPars[5]
        
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
