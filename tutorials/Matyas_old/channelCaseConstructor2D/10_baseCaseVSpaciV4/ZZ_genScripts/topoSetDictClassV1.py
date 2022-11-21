#!/usr/bin/python

#FILE DESCRIPTION=======================================================

# Python class containing the functions to ease up generation of
# topoSetDic files

#IMPORT BLOCK===========================================================
import math

class topoSetDictClass:
    """ Base python class containing functions to generate topoSetDict"""
    
    def retClearSetStr(self,setName,setType):
        """ return string for clearing the set of name and setType """
        retStr = []
        retStr.append('\t{\n')
        retStr.append('\t\tname\t%s;\n'%setName)
        retStr.append('\t\ttype\t%s;\n'%setType)
        retStr.append('\t\taction\tclear;\n')
        retStr.append('\t}\n')
        
        return retStr
        
    def retBoxToCellStr(self,setName,setType,action,boxPars):
        """ return string to create/add boxToCell """
        xMin,yMin,zMin,xMax,yMax,zMax = boxPars
        retStr = []
        retStr.append('\t{\n')
        retStr.append('\t\tname\t%s;\n'%setName)
        retStr.append('\t\ttype\t%s;\n'%setType)
        retStr.append('\t\taction\t%s;\n'%action)
        retStr.append('\t\tsource\tboxToCell;\n')
        retStr.append('\t\tsourceInfo\n')
        retStr.append('\t\t{\n')
        retStr.append('\t\t\tbox (%f %f %f) (%f %f %f);\n'%(xMin,yMin,zMin,xMax,yMax,zMax))
        retStr.append('\t\t}\n')
        retStr.append('\t}\n')
        
        return retStr
    
    def retZoneToCellStr(self,setName,setType,action,zonePars):
        """ return string to create/add boxToCell """
        zoneName = zonePars
        retStr = []
        retStr.append('\t{\n')
        retStr.append('\t\tname\t%s;\n'%setName)
        retStr.append('\t\ttype\t%s;\n'%setType)
        retStr.append('\t\taction\t%s;\n'%action)
        retStr.append('\t\tsource\tzoneToCell;\n')
        retStr.append('\t\tsourceInfo\n')
        retStr.append('\t\t{\n')
        retStr.append('\t\t\tname\t%s;\n'%(zoneName))
        retStr.append('\t\t}\n')
        retStr.append('\t}\n')
        
        return retStr
        
    def retPerfSourceActionStr(self,masterSetName,slaveSetName,setType,source,action):
        """ return string to perform new/add/remove/subset action """
        retStr = []
        retStr.append('\t{\n')
        retStr.append('\t\tname\t%s;\n'%masterSetName)
        retStr.append('\t\ttype\t%s;\n'%setType)
        retStr.append('\t\taction\t%s;\n'%action)
        retStr.append('\t\tsource\t%s;\n'%source)
        retStr.append('\t\tsourceInfo\n')
        retStr.append('\t\t{\n')
        retStr.append('\t\t\tset\t%s;\n'%slaveSetName)
        retStr.append('\t\t}\n')
        retStr.append('\t}\n')
        
        return retStr
        
    def retCellToZoneStr(self,setName,setType,action,cellPars):
        """ return string to create cellZone from cellSet """
        cellName = cellPars
        retStr = []
        retStr.append('\t{\n')
        retStr.append('\t\tname\t%s;\n'%setName)
        retStr.append('\t\ttype\t%s;\n'%setType)
        retStr.append('\t\taction\t%s;\n'%action)
        retStr.append('\t\tsource\tsetToCellZone;\n')
        retStr.append('\t\tsourceInfo\n')
        retStr.append('\t\t{\n')
        retStr.append('\t\t\tset\t%s;\n'%(cellName))
        retStr.append('\t\t}\n')
        retStr.append('\t}\n')
        
        return retStr
        
    def retFileOpenStr(self):
        """ returns the string to open file"""
        retStr = []
        retStr.append('actions\n(\n')
        return retStr
        
    def retFileCloseStr(self):
        """ returns the string to close file"""
        retStr = []
        retStr.append(');\n')
        return retStr
        
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
        retStr.append(' \t class \t\t dictionary; \n \t object \t topoSetDict; \n} \n')
        retStr.append('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n\n')
        
        return retStr
    
    def retFileFooterStr(self):
        """ returns the string with file footer"""
        retStr = []
        retStr.append('// ************************************************************************* //\n\n')
        return retStr
        
