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
# NOTES:
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
#   => Note: all the openFoam file footers are actually the same
#   => Note: in the header, I should put "location"

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
import io
import sys
import shutil as sh
import random
import math

#IMPORT BLOCK-CUSTOM====================================================
from blockMeshDictClassV5 import *
from topoSetDictClassV1 import *
from caseConstructorAuxFuncsVSpaci import *

#CUSTOM FUNCTIONS=======================================================
def dValsDict(medName):
    """ a dictionary to hold the data on different medias
    
        permeabilities """

    # Note: this should be easily extandable for "tiled" walls (just
    #       return more values in a list)

    return {
        'CF1' : 12.8612,
        'CF2' : 4.76617,
        'CF3' : 0.623756,
    }.get(medName,'CF1')

# -- I have put everything into function to call it from genContClass

#=======================================================================
#							EDITABLE
#=======================================================================

#-----------------------------------------------------------------------
# I/O DATA
#-----------------------------------------------------------------------
baseCase = '../10_baseCaseVSpaciV2/'
outFolder= '../20_toAltix/'
# name = 'TestSpaci'
name = 'unC_V16Spaci'

#-----------------------------------------------------------------------
# CASE/RUN SETTINGS
#-----------------------------------------------------------------------
# -- of.sh settings
queName     = 'batch'                                                   #que name for altix
wallTime    = '300:00:00'                                               #walltime in hours for altix
user        = 'hlavatyo'                                                #user name on altix
runScript   = 'Allrun-parallel'                                         #name of the runscript for the simulation
nNodes      = 1                                                         #number of nodes to use
nCoresPNode = 32                                                         #number of cores per node


nCores      = nNodes*nCoresPNode                                        #number of cores to run the case on
startTime   = 0                                                         #simulation startTime
endTime     = 1000                                                      #simulation endTime
wrInt       = 300                                                       #simulation write interval

#-----------------------------------------------------------------------
# GEOMETRY DATA
#-----------------------------------------------------------------------
# width dimensions
# width dimensions

    # name=dirTu
    # baseCase = dirBC

WWl = 137e-6                                                            #wall thickness
WCh = 1e-3*0.5                                                       #channel width (square)
# nChX= 3                                                                 #number of channels in hor. dir.
# nChY= 2                                                                 #number of channels in hor. dir.
nChX= 2
nChY= 2

# length dimensions
# LBf  = 7.0e-3                                                           #length of buffer in front of channels
# LPl  = 7.0e-3                                                           #plug length
# LCh  = 114.5e-3                                                     	#channel length (total, including the plugs, but WITHOUT the buffers)

# length dimensions
LBf  = 7.0e-3                                                           #length of buffer in front of channels
LPl  = 1.0e-3                                                           #plug length
LCh  = 76.2e-3                                                     	#channel length (total, including the plugs, but WITHOUT the buffers)
# LCh  = 7.02e-3                                                     	#channel length (total, including the plugs, but WITHOUT the buffers)





# relative parameters defining the position of the coating in the inlet and outlet channel (related to the length of the channel without plugs)
# NOTE: at the moment there is a problem with coating under plugs
# use topoSet script to cut it
# 0 -- end of the outlet plug
# 1 -- start of the inlet plug
bCI = 0.0								                                # begining of the coating in inlet channel (not working at the moment, use topoSet)
eCI = 1.0								                                # end of the coating inlet channel

bCO = 0.0								                                # begining of the coating in outlet channel
eCO = 1.0								                                # end of the coating outlet channel (not working at the moment, use topoSet)

# ~ LCh  = 34e-3							#reduced channel length for testing
# LCh  = LCh - LPl*2.0                                                    #plugs are inside the channel -> update LCh (remove the plugs)

LLst = [LBf,LCh,LBf]
#~ LLst = [LBf,LPl]
    

#-----------------------------------------------------------------------
# MESH DATA
#-----------------------------------------------------------------------
# ~ dX,dY,dZ   = 20e-6,20e-6,20e-6					#domain discretization
# dX,dY,dZ   = 80e-6,80e-6,80e-6
# dX,dY,dZ   = 100e-6,100e-6,100e-6
# dX,dY,dZ   = 60e-6,60e-6,60e-6
# dX,dY,dZ   = 24e-6,24e-6,24e-6
# dX,dY,dZ   = 50e-6,50e-6,50e-6
# dX,dY,dZ   = 30e-6,30e-6,30e-6
# dX,dY,dZ   = 40e-6,40e-6,40e-6
dX,dY,dZ   = 30e-6,30e-6,90e-6

# Note: significantly more should be available for changes, but this is
#       not implemented at the moment

#-----------------------------------------------------------------------
# CASE PARAMETERS
#-----------------------------------------------------------------------
spVel    = 40000
VcapLst = [20.0,40.0,60.0,80.0]
VcapLst = [20.0]
# VcapLst = [10e-3]
# VcapLst = [10.0,30.0,33.0]
# VcapLst = [80.0]
fracLst = [1./16.,1./8,2./8,3./8,4./8,5./8,6./8,7./8]
# fracLst = [30.0]
# fracLst = [0.5]
# fracLst = [40.0]
# fracLst = [60.0]
# kapilara Data
radius = [320e-6/2,450e-6/2]
# radius = [75e-6/2,200e-6/2]
# radius = [250e-6/2,350e-6/2]
offset = [0e-6,0e-6]
# fraction = 0.5

# ~ uIn = 3.0
# ~ parLstB = [								#list with geometry files
# ~'B05',
# ~'B1',
# ~ 'edgSAutoV1',
# ~'B2',
# ~'B25',
# ~'B3',
# ~'B35',
# ~'B4',
# ~'B45',
# ~'B5'
# ~ ]

nuG       = 1.5e-5                                                      #gas kinematic viscosity, m2s-1

#~ mediaName = 'CF2'
mediaName = 'CF3'

# edgePtsFl = 'testPtsSpaciMS'
edgePtsFl = 'testPts'
# ~ edgePtsFl = 'edgSAutoV1'                                                #name of the file with edge interpolation points file
# ~ edgePtsFl = 'A05'

#-----------------------------------------------------------------------
# POROSITY DATA
#-----------------------------------------------------------------------
nPorZones = 1 
# Note: number of different porous zones - with different porosities
# Note: the idea is to put together the zones as tiles to construct
#       the channel walls
# Note: this should help with determination of correct weights for
#       filters and so on

zoneL = 11.0
# Note: length of the zone. The zone is always taken to be as wide as
#       half of the channel and as deep as half of the wall. the only
#       remaining parameter is the length of one zone "tile"
# Note: this length is expressed in the number of cells (coding purposes)


porNames = ['porosity%d'%ind for ind in range(nPorZones)]

basePorName = 'inWallPorosity'                                          #porosity to be used in the walls

porNames[0] = 'coating'                                                 #porosity to be used in the coating

dVal   = dValsDict(mediaName)

dVals  = [5.0e-13,2.7600003091200346e-15,2.7600003091200346e-15]	    		        #porosity values
# ~dVals  = [val*1e-15 for val in dVals]                                #correct scaling
dVals  = [1.0/val for val in dVals]                                     #go from porosity do d coefficient
# Note: we use the same porosity in all the spatial directions.
# Note: the base porosity (the mean one) has to be alway at the 0th
#       position
# Note: shouldn't I use the same inWallPorosity for all the three media?
#       (the cross in the middle should not be affected by the coating)

fVals  = [0.0 for val in dVals]                                         # f values (not used at the time)

# Note: at the moment, we use always the same coordinateSystem for all
#       the porosities, so this is not changed

#########PREFERABLY DO NOT EDIT#########################################
pyFolder = '../00_pyCodes/'                                             #source directory for python codes
pyList   = [                                                            #python files
    'caseConstructorV3AllVels',
    'blockMeshDictClassV5',
    'topoSetDictClassV1',
    'caseConstructorAuxFuncsVSpaci',
    # ~ edgePtsFl,
]
# Note: these files are copied into folder caseName/ZZ_genScripts/ for
#       future reference
# Note: update these files if needed - in case of new code versions

# load file with interpolation points
# ~ edgS = importInterpPoints(edgePtsFl)

#########DO NOT EDIT####################################################
for fraction in fracLst:
    for Vcap in VcapLst:
        VcapTu = Vcap / 60. * 1e-6 
        uIn = LCh/((1.0/spVel)*3600)
        edgS = importInterpPoints(edgePtsFl)
        #-------------------------------------------------------------------
        # AUXILIARY COMPUTATIONS
        #-------------------------------------------------------------------
        # -- compute estimates for turbulence variables
        dHChan = WCh*2                                                      #channel hydraulic diameter (including symmetry)
        ReFlow  = uIn*dHChan/nuG
        ITurbAm = 0.227*ReFlow**-0.100                                      #turbulence intensity over the pipe area, Russo,F. and Basse, N.T. (2016)
        k0Turb  = 1.5*(uIn*ITurbAm)**2.0
        eps0Turb= 0.1643*k0Turb**1.5/(0.07*dHChan)
        om0Turb = eps0Turb/k0Turb
        
        #print ReFlow,ITurbAm,k0Turb,om0Turb
        
        #-------------------------------------------------------------------
        # DATA GATHERING
        #-------------------------------------------------------------------
        # -- copy the data from the baseCase
        caseName = '%s_%.2e_%.2f_%.2f/'%(name,radius[0]*2,Vcap,fraction)
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
        
        # -- prepare the inputs (geometry)    

        LChP = []								#length of different parts of the coating
        
        cLstIO = [0.0,LBf,eCI,bCO,eCO, 1.0]
        sCLstIO = sorted(cLstIO)
        for minInd in range(len(sCLstIO)):
            if minInd >=1:
                if not (sCLstIO[minInd] == sCLstIO[minInd-1]): 
                    LChP.append(LCh*sCLstIO[minInd])
            else:
                    LChP.append(LCh*sCLstIO[minInd])         
        print(LChP)
        
        WLst = WCh,WWl
        
        # list of real width and height of channel without coating and wall
        WSCh = WCh + WWl/2 + edgS.verPoint2[0],WCh + WWl/2 - edgS.verPoint2[1]
        
        #list of widths of each block in x direction
        WLstX = [WSCh[0],-edgS.verPoint2[0]-WWl/2,WWl,-edgS.verPoint2[0]-WWl/2,WSCh[0]]
        WLstY = [WSCh[1],edgS.verPoint2[1]-WWl/2,WWl,edgS.verPoint2[1]-WWl/2,WSCh[1]]

        pt1 = [1*WCh+WWl+radius[1]+-edgS.verPoint2[0]-WWl/2+offset[0],1*WCh+WWl+radius[1]+edgS.verPoint2[1]-WWl/2+offset[1],LBf+fraction*LCh]
        pt2 = [pt1[0],pt1[1],LBf*2+LCh+0.1]
        # ~ LLst = [LBf,LPl,LCh,LPl,LBf]
        # LLst = [LBf,LPl]
        # for ind in range(len(LChP)-1):
        #     LLst.append(float(LChP[ind+1]-LChP[ind]))
        # LLst.append(LPl)
        # LLst.append(LBf)

        # ~ LLst = [LBf,LPl,LCh,LPl,LBf]
        print(LLst)

        bChP = []
        bChP.append([True, False])
        for ind in range(len(sCLstIO)-1):
            if not (sCLstIO[ind+1] == sCLstIO[ind]):
                mid = (sCLstIO[ind+1]+sCLstIO[ind])/2
                print(mid)
                bIO = []
                if mid >= bCI and mid <= eCI:
                    bIO.append(True)
                else:
                    bIO.append(False)
                if mid >= bCO and mid <= eCO:
                    bIO.append(True)
                else:
                    bIO.append(False)
                bChP.append(bIO)
        bChP.append([False,True])
        print (bChP)
        
        #-------------------------------------------------------------------
        # -- prepare the inputs (mesh)
        # ~ # x direction discretization
        # ~ nCXLst = [int(round(length/dX)) for length in WLstX]
        nTimesCoat = 1
        nTimesWall = 1

        nCXLst = []
        for i in range(len((WLstX))):
            if i%5 == 1:
                nCXLst.append(int(round(WLstX[i]/dX))*nTimesCoat )
            elif i%5 == 2:
                nCXLst.append(int(round(WLstX[i]/dX))*nTimesWall )
            elif i%5 == 3:
                nCXLst.append(int(round(WLstX[i]/dX))*nTimesCoat )
            else:
                nCXLst.append(int(round(WLstX[i]/dX)))
                
        nCYLst = []
        for i in range(len((WLstY))):
            if i%5 == 1:
                nCYLst.append(int(round(WLstY[i]/dY))*nTimesCoat )
            elif i%5 == 2:
                nCYLst.append(int(round(WLstY[i]/dY))*nTimesWall )
            elif i%5 == 3:
                nCYLst.append(int(round(WLstY[i]/dY))*nTimesCoat )
            else:
                nCYLst.append(int(round(WLstY[i]/dY)) )
        
        print(nCXLst)
        # ~ # y direction discretization
        # ~ nCYLst = [int(round(length/dY)) for length in WLstX]
        
        # ~ # z direction discretization
        nCZLst = [int(round(length/dZ)) for length in LLst]


        
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
        
        # Note: out of these, I actually use only dZLst
        
        #-------------------------------------------------------------------
        # -- define which blocks are to be used
        l = 0
        boolBlockUse = []
        boolBlockUse.append([[True for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])#use the whole first layer
        l+=1
        boolBlockUse.append([[True for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])
        # for j in range(nChY*len(WLstY)):
        #     for k in range(nChX*len(WLstX)):
        #         if funBlockNotUsed(j,k,'out'):                                # using imported function (in - inlet plug, out - outlet plug) => return True if not used
        #             boolBlockUse[l][j][k] = False
        #         elif funBlockNotUsed(j,k,'in'):
        #             boolBlockUse[l][j][k] = False

        # for i in range(4):
        #     for j in range(4):
        #         boolBlockUse[l][8+i][8+j] = True

        # l+=1
        # for p in range(len(LChP)-1):
        #     boolBlockUse.append([[True for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])#use the middle layer
        #     l+=1
        # boolBlockUse.append([[True for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])

        # for j in range(nChY*len(WLstY)):
        #     for k in range(nChX*len(WLstX)):
                # if i == 9 and j == 9:
                # if funBlockNotUsed(j,k,'out'):
                    # boolBlockUse[l][j][k] = False
        # for i in range(4):
        #     for j in range(4):
        #         # print( boolBlockUse[l][8+i][8+i])
        #         boolBlockUse[l][8+i][8+j] = False
        l+=1
        boolBlockUse.append([[True for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])#use the whole last layer
        l+=1

        # for j in range(len(LLst)):
        #     for i in range(20):
        #         boolBlockUse[j][0][i] = False
        #         boolBlockUse[j][1][i] = False
        #         boolBlockUse[j][-1][i] = False
        #         boolBlockUse[j][-2][i] = False
        #         # boolBlockUse[j][-1][i] = False
        #         # boolBlockUse[j][-2][i] = False
        #         boolBlockUse[j][i][0] = False
        #         boolBlockUse[j][i][1] = False
        #         boolBlockUse[j][i][-1] = False
        #         boolBlockUse[j][i][-2] = False
                # boolBlockUse[-1][j][i] = False
                # boolBlockUse[-2][j][i] = False

        # define which blocks should be included into the porousZoneWall
        l=0
        boolPorousWall = []
        boolPorousWall.append([[False for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])
        l+=1
        # boolPorousWall.append([[False for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])
        # for j in range(nChY*len(WLstY)):
        #     for k in range (nChX*len(WLstX)):
        #         if funWallZone(j,k):                                        # using imported function => return True if included
        #             boolPorousWall[l][j][k] = True
        # l+=1
        # for p in range(len(LChP)-1):
        #     boolPorousWall.append([[False for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])
        #     for j in range(nChY*len(WLstY)):
        #         for k in range (nChX*len(WLstX)):
        #             if funWallZone(j,k):
        #                 boolPorousWall[l][j][k] = True
        #     l+=1
        boolPorousWall.append([[False for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])
        for j in range(nChY*len(WLstY)):
            for k in range (nChX*len(WLstX)):
                if funWallZone(j,k):
                    boolPorousWall[l][j][k] = True
        l+=1
        boolPorousWall.append([[False for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])
        l+=1
        
        # define which blocks should be included into the porousZoneCoating
        
        l=0
        boolPorousCouting = []
        boolPorousCouting.append([[False for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])
        l+=1
        boolPorousCouting.append([[False for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])
        for j in range(nChY*len(WLstY)):
            for k in range (nChX*len(WLstX)):
                if funCoatZone(j,k):                                        # using imported function => return True if includeds
                    boolPorousCouting[l][j][k] = True
        l+=1
        # for p in range(1,len(bChP)-1,1):
        #     print(p)
        #     boolPorousCouting.append([[False for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])
        #     for j in range(nChY*len(WLstY)):
        #         for k in range (nChX*len(WLstX)):
        #             if funCoatZone(j,k,p,bChP):
        #                 boolPorousCouting[l][j][k] = True
        #     l+=1
        # boolPorousCouting.append([[False for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])
        # for j in range(nChY*len(WLstY)):
        #     for k in range (nChX*len(WLstX)):
        #         if funCoatZone(j,k,-1,bChP):
        #             boolPorousCouting[l][j][k] = True
        # l+=1
        boolPorousCouting.append([[False for i in range(nChX*len(WLstX))] for val in range(nChY*len(WLstY))])
        l+=1
        
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
                        if boolPorousWall[i][j][k]:
                            cZone = basePorName
                        elif boolPorousCouting[i][j][k]:
                            cZone = porNames[0]
                        else:
                            cZone = None
                        LcLst = xC,yC,zC
                        grLst = grX,grY,grZ
                        dmLst = WLstX[k%5], WLstY[j%5], LLst[i]
                        nCLst = nCXLst, nCYLst, nCZLst
                        blockPars = funBlockPars(blockCount, LcLst, dmLst, nCLst, grLst, cZone,i, j, k,WLst,WBCh)      #imported function from caseConstructorAuxFuncs
                        blockLst[i][j][k] = hexBlockClass(blockPars) 
                        blockCount+=1
                    xC+=WLstX[k%5]
                yC+=WLstY[j%5]
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
        
        # Note: this is if hell, but it is self-explanatory
        
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
        
        # Note: nothing else is changed at the moment
        
        #-------------------------------------------------------------------
        # FLUID PROPERTIES UPDATE
        #-------------------------------------------------------------------
        
        # transportProperties - fluid kinematic viscosity
        idStr = [
            'nu              [0 2 -1 0 0 0 0]',
        ]         
        
        pVals = [nuG]
        
        # write everything to the file
        with open(caseDir + './constant/transportProperties', 'r') as file:
            # read a list of lines into data
            data = file.readlines()
            
        for j in range(len(idStr)):
            for i in range(len(data)):
                fInd = data[i].find(idStr[j])
                if fInd>-1:
                    data[i] = data[i][:fInd] + idStr[j] + '\t' + repr(pVals[j]) + ';\n'
        
        with open(caseDir + './constant/transportProperties', 'w') as file:
            file.writelines( data )
            
        #-------------------------------------------------------------------
        # POROSITY PROPERTIES UPDATE
        #-------------------------------------------------------------------
        # Note: the code structure here is slightly different than in orher parts
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
        currD = [dVals[0] for dInd in range(3)]
        currF = [fVals[0] for dInd in range(3)]
        
        currPorName = basePorName
        currActive  = 'yes'
        currZoneName= basePorName
        
        porProps    = [currPorName,currActive,currZoneName,currD,currF]
        
        porStr = porClassMember.retPorEntryStr(porProps)
        for line in porStr:
            pPD.write(line)
            
        # -- write the rest of the porosity classes
        for ind in range(nPorZones):
            currD = [dVals[ind+1] for dInd in range(3)]
            currF = [fVals[ind+1] for dInd in range(3)]
                
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
        
        #-------------------------------------------------------------------
        # RUN SCRIPTS UPDATE
        #-------------------------------------------------------------------
        
        # of.sh
        idStr = [
                    '#PBS -M ',
                    '#PBS -N ',
                    '#PBS -q ',
                    '#PBS -l nodes=',
                    '#PBS -l walltime=',
                    'caseDir=',
                    'bash '
                    ]
        
        caseName = caseDir.split('/')[-2]
        
        # -- auxiliary parameters
    

        # write everything to the file
        with open(caseDir + '/system/snappyHexMeshDict', 'r') as file:
            # read a list of lines into data
            data = file.readlines()

        for ind in range(len(data)):
            if data[ind].find('point1') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\tpoint1\t(%g %g %g);\n'%(pt1[0],pt1[1],pt1[2])
            if data[ind].find('point2') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\tpoint2\t(%g %g %g);\n'%(pt2[0],pt2[1],pt2[2])
            if data[ind].find('chInt1MinMoje') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\tmin (-0.1 -0.1 %g);	//chInt1MinMoje\n'%(pt1[2]-150e-6)
            if data[ind].find('chInt1MaxMoje') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\tmax (1 1 %g); //chInt1MinMoje\n'%(pt1[2]+150e-6)
            if data[ind].find('radius') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\tradius\t%g;\n'%(radius[1])

        with open(caseDir + '/system/snappyHexMeshDict', 'w') as file:
            file.writelines(data)

        # write everything to the file
        with open(caseDir + '/system/topoSetDict', 'r') as file:
            # read a list of lines into data
            data = file.readlines()

        for ind in range(len(data)):
            if data[ind].find('p1') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\t\tp1\t(%g %g %g);\n'%(pt1[0],pt1[1],-0.1)
            if data[ind].find('p2') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\t\tp2\t(%g %g %g);\n'%(pt2[0],pt2[1],pt2[2])
            if data[ind].find('outerRadius') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\t\touterRadius\t%g;\n'%(radius[1]+50e-6)
            if data[ind].find('innerRadius') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\t\tinnerRadius\t%g;\n'%(radius[0])

        with open(caseDir + '/system/topoSetDict', 'w') as file:
            file.writelines(data)


        with open(caseDir + '/0.org/U', 'r') as file:
            # read a list of lines into data
            data = file.readlines()

        uIn_cyl = Vcap/60*1e-6 / (math.pi*radius[0]**2)
        for ind in range(len(data)):
            if data[ind].find('UIN_intF') >= 0:#!!do not put comment on executed temperature
                data[ind] = 'internalField\tuniform\t(0 0 %g);\n'%(uIn)
            if data[ind].find('UIN_w') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\t\tvalue\tuniform (0 0 %g);\n'%(uIn)
            if data[ind].find('UIN_cyl') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\t\tvalue\tuniform (0 0 %g);\n'%(uIn_cyl)
        
        with open(caseDir + '/0.org/U', 'w') as file:
            file.writelines(data)

        with open(caseDir + '/system/sample', 'r') as file:
            # read a list of lines into data
            data = file.readlines()

        for ind in range(len(data)):
            if data[ind].find('z_StartCap') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\t\tpoint\t(0.002185 0.002185 %g);\n'%(pt1[2])
        
        with open(caseDir + '/system/sample', 'w') as file:
            file.writelines(data)

        with open(caseDir + '/system/topoSetDict.3', 'r') as file:
            # read a list of lines into data
            data = file.readlines()

        x = [pt2[0]-WWl-WCh*2,pt2[0]+WWl+WCh*2,pt2[0]+2*(WWl+WCh*2),pt2[0]-WWl-WCh*2]
        y = [pt2[1]-WWl-WCh*2,pt2[1],pt2[1],pt2[1]]

        ch1 = WCh + WWl/2
        ch2 = 3*WCh + 1.5*WWl

        for ind in range(len(data)):
            if data[ind].find('x1') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\t\tp1 (%g %g 0.0);\n'%(x[0],y[0])
                data[ind+1] = '\t\t\tp2 (%g %g 1.0);\n'%(x[0],y[0])
            
            if data[ind].find('boxCH1') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\t\tbox   (%g %g 0) (%g %g 1);\n'%(ch1,ch1,ch2,ch2)
                # data[ind+1] = '\t\t\tp2 (%g %g 1.0);\n'%(xch1,ych1)

            if data[ind].find('x2') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\t\tp1 (%g %g 0.0);\n'%(x[1],y[1])
                data[ind+1] = '\t\t\tp2 (%g %g 1.0);\n'%(x[1],y[1])
            
            if data[ind].find('x3') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\t\tp1 (%g %g 0.0);\n'%(x[2],y[2])
                data[ind+1] = '\t\t\tp2 (%g %g 1.0);\n'%(x[2],y[2])
            
            if data[ind].find('x4') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\t\tp1 (%g %g 0.0);\n'%(x[3],y[3])
                data[ind+1] = '\t\t\tp2 (%g %g 1.0);\n'%(x[3],y[3])

            # if data[ind].find('radiusCH') >= 0:#!!do not put comment on executed temperature
            #     data[ind] = '\t\t\tradius\t%g;\n'%(WCh)
            
            if data[ind].find('radiusCAP') >= 0:#!!do not put comment on executed temperature
                data[ind] = '\t\t\tradius\t%g;\n'%(radius[0])

        with open(caseDir + '/system/topoSetDict.3', 'w') as file:
            file.writelines(data)



