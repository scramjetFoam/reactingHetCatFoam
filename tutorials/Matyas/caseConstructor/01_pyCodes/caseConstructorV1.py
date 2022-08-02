#!//sr/bin/python

#FILE DESCRIPTION=======================================================
# Python script used for an automated case construction
#
# This is the variant of the script for the automatic creations of the
# steady state RANS simulations for the liquid-liquid flow in ejector
#
# used solver is simpleFoam
#
# Notes:
# - all the "constant" directory is directly copied from the base case
# - all the "system" directory is directly copied from the base case
#   (including fvSchemes, fvSolution and controlDict)
# - the only changed thing are the boundary conditions and geometry
#
#
   
#LICENSE================================================================
#  caseConstructor.py
#
#  Copyright 2015-2019 Martin Isoz <martin@Poctar>
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

#IMPORT BLOCK===========================================================
import os
import math
import io
import sys
import numpy as np
import shutil as sh
# from auxiliarFuncs import *

#IMPORT BLOCK-CUSTOM====================================================
# custom functions------------------------------------------------------
geomGenerator	= "fblockMeshDictGenV1"

importCommand = "from " + geomGenerator + " import *"
exec(importCommand)
geomGenerator   += ".py"

#########EDITABLE#######################################################

#INPUT PARAMETERS=======================================================
# -- case defining parameters
# flow rate at the inlet 
QInLst = [1.0e-3]

# outlet pressure
p0 = 101325

# TURBULENCE INTENSITY AT INLET
I0   = 0.02

# -- geometry and meshing parameters
# inlet diameter
DIn = 0.05

# outlet diameter
DOut = 0.06

# length
L = 0.3

# wedge angle in degrees
wAng = 5.0
    
# origin
geomOrig = [0.0,0.0,0.0]                                                 #geometry origin (middle of the inlet)

# cell size
cellSize = [1.0e-3, 1.0e-3]

# -- case run properties - of.sh file
nCores      = 2
startTime   = 0                                                         #simulation startTime
endTime     = 1
wrInt       = 0.01                                                      #simulation write interval
queName     = "Mshort"                                                   #que name for kraken
wallTime    = "2-00:00:00"                                               #walltime in hours for kraken
nameSpec    = "DIn"
specFolder  = "test/"

# -- case specification
genDir = "../"
baseCase = genDir + "10_baseCase/"
baseDir	= genDir + "ZZ_cases/" + specFolder

for QIn in QInLst:
    caseDir     = (baseDir + "sFLJ_QIn" + 
                repr(round(QIn*1e3,4)) + "_genCaseV1" + "_" + nameSpec + "_" + str(eval(nameSpec)) + "/"
                )

    geomPars = [
        DIn,
        DOut,
        L,
        wAng,
    ]   

    #########PREFERABLY DO NOT EDIT#########################################
    
    #COPY CASE BASICS FROM THE BASECASE=====================================
    if os.path.isdir(caseDir):                                              #ensure, that the caseDir is clear
        sh.rmtree(caseDir)
        
    sh.copytree(baseCase,caseDir)                                             #copy data to caseDir
    
    #GENERATE BLOCKMESHDICT (CASE GEOMETRY)=================================
    blockMeshClass = genBlockMeshDict(geomPars,geomOrig,cellSize,caseDir)
    
    #SPECIFY CURRENT SCRIPT VERSIONS========================================
    blockMeshClass  += ".py"
    caseConstructor = os.path.basename(__file__)
   
    #COPY CURRENT SCRIPT VERSIONS===========================================
    scFolder= genDir + "/01_pyCodes/"                                              #folder with Scripts sources
    scNames = [ 
                geomGenerator,
                blockMeshClass,
                caseConstructor,
    ]                                   
    for scName in scNames:
        sh.copyfile(scFolder + scName,caseDir + scName)     #copy current script version
    
    #CASE CONSTANTS AND CALCULATIONS========================================
    # input data------------------------------------------------------------
    # -- liquid properties
    rhoL,muL= 1000.0,1.0e-3                                             # water
    nu      = muL/rhoL                                                  # kinematic viscosity
    
    # ESTIMATE TURBULENCE VARIABLES
    uIn = QIn/(math.pi*(DIn*0.5)**2.0)

    # length scale
    l = 0.1*DIn

    # model coefficient
    Cmu = 0.09

    # turbulence variables
    k0 = 3/2*uIn**2*I0**2
    e0 = Cmu*k0**(3/2)/l
    w0 = e0/Cmu/k0
    nut0 = k0/w0
    
    #OPEN AUTOGENERATED README FILE=========================================
    README  = open(caseDir + "./README","a")                                #open file to append
    README.write("\ncaseDir:" + caseDir + "\n\n")
    # -- start by writing basic case info and geometry
    README.write("\ngeometrical parameters [m]\n")
    geomParsToSave = ["DIn", "DOut", "L", "wAng"]

    for i in geomParsToSave:
        README.write(i + " = " + str(eval(i)) + "\n")

    README.write("\nother ineteresting facts\n")
    README.write("cellSize  \t = \t " + repr(cellSize) + " m\n")
    README.write("rho       \t = \t " + repr(rhoL) + " kgm-3\n")
    README.write("nu        \t = \t " + repr(nu) + " m2s-1\n")
    README.write("k0       \t = \t %.4g"%(k0) + " Jkg-1\n")
    README.write("w0       \t = \t %.4g"%(w0) + " s-1\n")
    README.write("e0       \t = \t %.4g"%(e0) + " s-1\n")
    README.write("nut0     \t = \t %.4g"%(nut0) + " Jskg-1\n")
    README.write("nCores    \t = \t " + repr(nCores) + "\n")
    README.write("startTime \t = \t " + repr(startTime) + " s\n")
    README.write("endTime   \t = \t " + repr(endTime) + " s\n")
    
    #BC FILES MODIFICATION==================================================
    #-----------------------------------------------------------------------
    # FUNCTION CALL
    #-----------------------------------------------------------------------
    print ("ADJUSTING BC===============================\n\n")
    README.write("\n 0.org==============================================\n")
    #-----------------------------------------------------------------------
    # U
    #-----------------------------------------------------------------------
    #
    # Boundary conditions for the velocity field
    #
    README.write("\n U\n")
    
    pVals   = ["uniform (" + repr(uIn) + " 0 0)"]           #inlet liquid velocity speed
    
    idStr   = ["uniform (u0 0 0)"]
    
    # write everything to the file
    with open(caseDir + "./0.org/U", "r") as file:
        # read a list of lines into data
        data = file.readlines()
        
    for j in range(len(idStr)):
        for i in range(len(data)):
            fInd = data[i].find(idStr[j])
            if fInd>-1:
                data[i] = data[i][:fInd] + pVals[j] + ";\n"
    
    with open(caseDir + "./0.org/U", "w") as file:
        file.writelines( data )
        README.writelines( data )                                           #write to readme
        
    print ("DONE=======================================\n\n")
    
    #-----------------------------------------------------------------------
    # p
    #-----------------------------------------------------------------------
    #
    # Boundary conditions for the pressure field
    #
    README.write("\n p\n")
    
    pVals   = ["uniform " + repr(p0)]
   
    idStr   = ["uniform p0"]
    
    # write everything to the file
    with open(caseDir + "./0.org/p", "r") as file:
        # read a list of lines into data
        data = file.readlines()
        
    for j in range(len(idStr)):
        for i in range(len(data)):
            fInd = data[i].find(idStr[j])
            if fInd>-1:
                data[i] = data[i][:fInd] + pVals[j] + ";\n"
    
    with open(caseDir + "./0.org/p", "w") as file:
        file.writelines( data )
        README.writelines( data )                                           #write to readme
        
    print ("DONE=======================================\n\n")
    
    #-----------------------------------------------------------------------
    # k, omega and epsilon
    #-----------------------------------------------------------------------
    #
    # Boundary conditions for turbulent variables
    #
    varNames = ["k","omega","epsilon","nut"]
    varVals  = [
        k0,
        w0,
        e0,
        nut0,
    ]
    varNames = []
    
    for varInd in range(len(varNames)):
        varNm   = varNames[varInd]
        var0 = varVals[varInd]
        
        README.write("\n %s\n"%(varNm))
        
        pVals   = [
            "internalField   uniform " + repr(var0),
        ]
        
        idStr   = [
            "internalField   uniform var0",
        ]
        
        # write everything to the file
        with open(caseDir + "./0.org/%s"%(varNm), "r") as file:
            # read a list of lines into data
            data = file.readlines()
            
        for j in range(len(idStr)):
            for i in range(len(data)):
                fInd = data[i].find(idStr[j])
                if fInd>-1:
                    data[i] = data[i][:fInd] + pVals[j] + ";\n"
        
        with open(caseDir + "./0.org/%s"%(varNm), "w") as file:
            file.writelines( data )
            README.writelines( data )                                        #write to readme
            
    print ("DONE=======================================\n\n")
    
    #CONSTANTS DIRECTORY FILES MODIFICATIONS================================
    print ("ADJUSTING FILES IN ./CONSTANTS=============\n\n")
    README.write("\n CONSTANTS==========================================\n")
    README.write("\n transportProperties\n")
    
    idStr = [
        "nu0",
    ]         
    
    pVals = [[nu]]
    
    # write everything to the file
    with open(caseDir + "./constant/transportProperties", "r") as file:
        # read a list of lines into data
        data = file.readlines()
        
    for j in range(len(idStr)):
        k = 0
        for i in range(len(data)):
            fInd = data[i].find(idStr[j])
            if fInd>-1:
                data[i] = data[i][:fInd] + "\t" + repr(pVals[j][k]) + ";\n"
                k = k+1
    
    with open(caseDir + "./constant/transportProperties", "w") as file:
        file.writelines( data )
        README.writelines( data )                                       #write to readme
        
    print ("DONE=======================================\n\n")
    #SYSTEM DIRECTORY FILES MODIFICATIONS===================================
    print ("ADJUSTING FILES IN ./SYSTEM================\n\n")
    README.write("\n SYSTEM=============================================\n")
    #-----------------------------------------------------------------------
    # decomposeParDict
    #-----------------------------------------------------------------------
    #
    # decomposes the case for run on multiple cores
    #
    README.write("\n decomposeParDict\n")
    
    idStr = ["numberOfSubdomains "]
    
    pVals = [repr(nCores)]
    
    # write everything to the file
    with open(caseDir + "./system/decomposeParDict", "r") as file:
        # read a list of lines into data
        data = file.readlines()
        
    for j in range(len(idStr)):
        for i in range(len(data)):
            fInd = data[i].find(idStr[j])
            if fInd>-1:
                data[i] = data[i][:fInd] + idStr[j] + "\t" + pVals[j] + ";\n"
    
    with open(caseDir + "./system/decomposeParDict", "w") as file:
        file.writelines( data )
        README.writelines( data )                                           #write to readme

    #-----------------------------------------------------------------------
    # controlDict
    #-----------------------------------------------------------------------
    #
    # creates initial condition for the case (presence of the liquid)
    #
    README.write("\n controlDict\n")
    
    idStr = ["startTime ","endTime ","writeInterval "]
    
    pVals = [repr(startTime),repr(endTime),repr(wrInt)]
    
    # write everything to the file
    with open(caseDir + "./system/controlDict", "r") as file:
        # read a list of lines into data
        data = file.readlines()
        
    for j in range(len(idStr)):
        for i in range(len(data)):
            fInd = data[i].find(idStr[j])
            if fInd>-1:
                data[i] = data[i][:fInd] + idStr[j] + "\t" + pVals[j] + ";\n"
    
    with open(caseDir + "./system/controlDict", "w") as file:
        file.writelines( data )
        README.writelines( data )                                           #write to readme
        
        
    print ("DONE=======================================\n\n")
    #RUN SCRIPTS PREPARATION================================================
    print ("PREPARING RUN SCRIPTS======================\n\n")
    README.write("\n RUN SCRIPTS========================================\n")
        
    #-----------------------------------------------------------------------
    # ./Allrun-slurm
    #-----------------------------------------------------------------------
    #
    README.write("\n Allrun-slurm\n")
    
    idStr = [
                "#SBATCH -J ",
                "#SBATCH -p ",
                "#SBATCH -n",
                "#SBATCH --time=",
                "srun -n",
                ]
    
    caseName = caseDir.split("/")[-2]
    
    pVals = [caseName,queName,repr(nCores),wallTime,repr(nCores) + " $application -parallel >> log.$application"]
    
    # write everything to the file
    with open(caseDir + "./Allrun-slurm", "r") as file:
        # read a list of lines into data
        data = file.readlines()
        
    for j in range(len(idStr)):
        for i in range(len(data)):
            fInd = data[i].find(idStr[j])
            if fInd>-1:
                data[i] = data[i][:fInd] + idStr[j] + pVals[j] + "\n"
    
    with open(caseDir + "./Allrun-slurm", "w") as file:
        file.writelines( data )
        README.writelines( data )                                           #write to readme
    
    print ("DONE=======================================\n\n")
    #CLOSE THE AUTOGENERATED README FILE====================================
    README.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n\n")
    README.close()
