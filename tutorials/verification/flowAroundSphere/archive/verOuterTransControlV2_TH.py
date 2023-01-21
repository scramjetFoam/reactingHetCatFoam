# verOuterTransControlV2.py
# -- Script for verification simulation of conjugated mass transport
# -- NOTE: SCRIPT USES ARGUMENTS
#       -- "makeMesh": prepare mesh for each cellSize 
#           (this must be done manually if protoMesh doesn't exist yet)
#       -- "runSim": run the simulation
#       -- "showPlots": shows Plots
#       -- "getCsv": generate csv files for TeX plots 

# -- imports
import numpy as np
import os
import re
import shutil as sh
import matplotlib.pyplot as plt
import sys
from auxiliarFuncs import *

# -- set solver to be used
solverLst = ['reactingHetCatSimpleFoam','scalarTransportFoamCO']
solver = solverLst[0]

# -- set number of the enthalpy corrections
numOfTCorr = 0
isothermal = (True if numOfTCorr==0 else False)  # isothermal logic

# -- script arguments logic
args = sys.argv
runSim = (True if 'runSim' in args else False)
showPlots = (True if 'showPlots' in args else False)
makeMesh = (True if 'makeMesh' in args else False)
getCsv = (True if 'getCsv' in args else False)

# -- directory naming
baseCaseDir = 'baseCase'
outFolder = 'ZZ_cases'
ZZZ_path = 'ZZZ_res'
ZZZ_file = 'flow.csv'
ZZZ_filepath = ZZZ_path+'/'+ZZZ_file

# -- case parameters
yInf = 0.01         # molar fraction farfar from sphere
p = 101325          # presure
Runiv = 8.314       # universal gas constant
sHr = -283e3        # standard reaction enthalpy (physical 283e3)	
T = 500             # temperature
DFreeZ = 1e-5       # Diffusivity in fluid
kappaEff = 1        # mass transfer coefficient
nu = 5.58622522e-05 # kinematic viscosity !?
EA = 90e3           # activation energy

# -- geometry
R = 0.01            # sphere radius
length1 = 15*R      # inlet <-> sphere centre
length2 = 45*R      # sphere centre <-> outlet
width = 15*R        # top|bottom wall <-> sphere centre

# -- list parameters
invLst = [0.22]                # inlet velocity
k0Lst = [1e9]                       # reaction pre-exponential factor
cellSizeLst = [0.45*R]              # FV cell Size
tortLst = [50]                   # tortuosity

# == ARCHIVED SETTINGS: 
# -- 12/11/2022 khyrm@multipede (12 cases):
# invLst = [0.0275,0.11,0.22,0.44]
# k0Lst = [1e9]
# cellSizeLst = [0.35*R]
# tortLst = [0.5,5,50]

# -- prepare prototype mesh for each cellSize
if makeMesh:
    for cellSize in cellSizeLst:
        # NOTE MK: This could be an auxiliary function.
        if not os.path.isdir('%s'%outFolder): os.mkdir('%s'%outFolder)
        # if not os.path.isdir('%s/protoMesh'%outFolder): os.mkdir('%s/protoMesh'%outFolder)
        meshDir = '%s/protoMesh/%g'%(outFolder,cellSize)
        print('Preparing mesh %s',meshDir)
        # -- check that meshDir is clean
        if os.path.isdir(meshDir): sh.rmtree(meshDir)
        # -- copy files
        sh.copytree(baseCaseDir,meshDir)
        changeInCaseFolders(meshDir,'system/blockMeshDict',['length1', 'length2', 'width','nDiscX','nDiscYZ'],[str(length1),str(length2),str(width),str(int((length1+length2)/cellSize)),str(int(2*width/cellSize))])
        changeInCaseFolders(meshDir,'system/snappyHexMeshDict',['spR'],[str(R)])
        os.chdir(meshDir)
        os.system('chmod u=rwx All*') # NOTE: Just to make sure.
        os.system('./Allmesh')
        os.chdir('../../../')


# -- create cases for:
cases = [(inv,k0,cellSize,tort) for inv in invLst for k0 in k0Lst for cellSize in cellSizeLst for tort in tortLst]
for case in cases:
    # parameters
    inv,k0,cellSize,tort = case
    k = k0*np.exp(-EA/(Runiv*T))
    DFree = DFreeZ
    DEff = DFree/tort*0.5

    caseName = 'intraTrans_yInf_%g_R_%g_T_%g_cS_%g_k0_%g_tort_%g_inv_%g'%(yInf,R,T,cellSize,k0,tort,inv)
    caseDir = '%s/%s/'%(outFolder,caseName)
    meshDir = '%s/protoMesh/%g'%(outFolder,cellSize)
    
    if runSim:
        print('Preparing case %s'%caseName)
        # -- check that caseDir is clean
        if os.path.isdir(caseDir): sh.rmtree(caseDir)   # ensure the caseDir is clear
        # -- copy files
        sh.copytree(meshDir,caseDir)
        # -- write parameters
        changeInCaseFolders(caseDir,'0.org/T',['isoT'],[str(T)])
        changeInCaseFolders(caseDir,'0.org/CO',['yCOSet'],[str(yInf)])
        changeInCaseFolders(caseDir,'0.org/U', ['inv'],[str(inv)])
        changeInCaseFolders(caseDir,'system/controlDict',['customSolver'],[solver])
        changeInCaseFolders(caseDir,'system/fvSolution',['nTCorr'],[str(numOfTCorr)])
        changeInCaseFolders(caseDir,'constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(k0),str(EA),str(sHr)])
        changeInCaseFolders(caseDir,'constant/transportProperties',['kappaEffSet','tortSet','DSet'],[str(kappaEff),str(tort),str(DFreeZ)])
        # -- run simulation
        os.chdir(caseDir)
        os.system('chmod u=rwx Allrun-parallel')
        os.system('./Allrun-parallel')
    else: 
        os.chdir(caseDir)

    k0Art = k0*np.exp(-EA/(Runiv*T))                            # definition in 2020 Chandra Direct numerical simulation of a noniso...
    thiele = R*np.sqrt(k0Art/DEff)
    rSqIdeal = 4./3*np.pi*R**3*k0Art*yInf*p/Runiv/T             # ideal reaction source
    rS = read_real_source()                                     # simulation real source
    eta_sim = rS/rSqIdeal

    # 
    with open('log.integrace','r') as fl:
        lines = fl.readlines()  
    inds = [0,0,0,0]
    names = ['areaAverage(batt1) of cCOS',"areaAverage(batt1) of gradCCO",'areaAverage(batt1) of yCOS','areaAverage(batt1) of gradYCO']
    for lineInd in range(len(lines)):
        for j in range(len(names)):
            if names[j] in lines[lineInd]:
                inds[j] = lineInd
    cCO = float(lines[inds[0]].split('=')[1])
    gradCCO = float(lines[inds[1]].split('=')[1])
    yCO = float(lines[inds[2]].split('=')[1])
    gradYCO = float(lines[inds[3]].split('=')[1])
    
    j, jY = gradCCO, gradYCO #/(4*np.pi*R**2)
    km, kmY = j/(yInf*p/Runiv/T-cCO), jY/(yInf-yCO)
    Sh, ShY = km*(2*R)/DFree, kmY*(2*R)/DFree

    # eta_corr
    # correlation: Sherwood -> Biot -> eta_corr
    # NOTE MK: nu could be read from log if it is not set
    Re = inv * (R*2) / nu
    Sc = nu / DFree
    Sh_corr = 2 + .6 * Re**(1/2) * Sc**(1/3)
    BiM_corr = Sh_corr/2 * DFree/DEff
    eta_corr = 3/(thiele**2) * (thiele/np.tanh(thiele)-1)/(1+(thiele/np.tanh(thiele)-1)/BiM_corr)

    # eta_sim 
    eta_sim = rS/rSqIdeal

    # eta_anal [old]
    BiM = km*R/DEff
    eta_anal = 3/(thiele**2) * (thiele/np.tanh(thiele)-1)/(1+(thiele/np.tanh(thiele)-1)/BiM)


    log_report(thiele,DEff,DFree,Re,Sh_corr,Sc,Sh,eta_sim,eta_corr,eta_anal,gradCCO,cCO,j,km,gradYCO,yCO,jY,kmY,ShY)
    os.chdir('../../')
    flow_csv(ZZZ_path,ZZZ_filepath,tort,Re,eta_sim,eta_corr,eta_anal)

    

if showPlots:
    for tort in tortLst:
        eta_plt(ZZZ_filepath, tort)
    eta_err_plt(ZZZ_filepath, tortLst, invLst)

if getCsv:
    generate_eta_csvs(ZZZ_path,ZZZ_filepath,tortLst)