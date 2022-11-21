#!/usr/bin/python

# FILE DESCRIPTION======================================================

# Python class containing the functions to ease up generation of
# blockMeshDict files

# IMPORT BLOCK==========================================================
import math


# ======================================================================

class blockClass:
    """Class for generating blockMeshDict blocks."""

    def __init__(self, blockPars):
        self.v0 = blockPars[0]  # initial vertex
        self.vCoords = blockPars[1]  # vertex coords
        self.nCells = blockPars[2]  # nOf cells
        self.grading = blockPars[3]  # grading
        self.name = blockPars[4]

    def retVIndices(self):
        """Return indices of the block vertices."""
        raise NotImplementedError()

    # methods returning block fases
    def retFXY0(self):
        """Return vertex indices of XYplane face with lower Zcoord."""
        ind = self.retVIndices()
        return [ind[0], ind[1], ind[2], ind[3]]

    def retFXYE(self):
        """Return vertex indices of XYplane face with higher Zcoord."""
        ind = self.retVIndices()
        return [ind[4], ind[5], ind[6], ind[7]]

    def retFXZ0(self):
        """Return vertex indices of XZplane face with lower Ycoord."""
        ind = self.retVIndices()
        return [ind[0], ind[1], ind[5], ind[4]]

    def retFXZE(self):
        """Return vertex indices of XZplane face with higher Ycoord."""
        ind = self.retVIndices()
        return [ind[3], ind[2], ind[6], ind[7]]

    def retFYZ0(self):
        """Return vertex indices of YZplane face with lower Xcoord."""
        ind = self.retVIndices()
        return [ind[0], ind[3], ind[7], ind[4]]

    def retFYZE(self):
        """Return vertex indices of YZplane face with higher Xcoord."""
        ind = self.retVIndices()
        return [ind[1], ind[2], ind[6], ind[5]]

    # methods returning block for blockMeshDict
    def retBlockString(self):
        """Return list of strings for blockMeshDict."""
        intro = '\thex\n'
        vert = '\t\t(' + ' '.join('%d' % vi for vi in self.retVIndices()) + ')'
        disc = '\t(' + ' '.join('%d' % nc for nc in self.nCells) + ')'
        grading_type = '\tSimpleGrading'
        grading = '\t(' + ' '.join('%s'%gr for gr in self.grading) + ')\n\n'
        if self.name is not None:
            name = '\t' + self.name
            return [intro, vert, name, disc, grading_type, grading]
        else:
            return [intro, vert, disc, grading_type, grading]

class hexBlockClass(blockClass):
    """Class for generating hex blockMeshDict blocks."""
    def __init__(self, blockPars):
        super().__init__(blockPars)
    
    def retVIndices(self):
        """Return indices of the block vertices."""
        return [self.v0 + i for i in range(8)]

class wedgeBlockClass(blockClass):
    """Class for generating wedge blockMeshDict blocks."""
    def __init__(self, blockPars):
        super().__init__(blockPars)
        self.vDbl = blockPars[4]  # doubled edges (two pairs)
        self.name = blockPars[5]
    
    def retVIndices(self):
        """Return indices of the block vertices."""
        ind = [self.v0 + i for i in range(8)]  # virtual vertex indices
        for i in range(8):
            for e in self.vDbl:
                if i == e[0]:
                    ind[i] = ind[e[1]]
                    for j in range(i + 1, 8):
                        ind[j] -= 1
        return ind

