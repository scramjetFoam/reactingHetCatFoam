# verOuterTransControlV3.py
# -- Script for verification simulation of conjugated mass transport
# -- NOTE SCRIPT ARGUMENTS:
#       -- "makeMesh": prepare mesh for each cellSize 
#           (this must be done manually if protoMesh doesn't exist yet)
#       -- "runSim": run the simulation
#       -- "showPlots": shows Plots
#       -- "getCsv": generate csv files for TeX plots 
#       -- "errMesh": evaluate error mesh dependency 
#           (for otherwise identical parameters) 
# -- NOTE V3 changelog: 
#       -- thieleLst + tortLst --> k0
#       -- ReLst + R --> invLst
#       -- flow.csv stores Thiele modulus
#       -- removed eta_anal
#       -- uses auxiliarFuncsV3
#       -- added error mesh dependency tests

# -- imports
import numpy as np
import os
# import re
import shutil as sh
import matplotlib.pyplot as plt
import sys
from auxiliarFuncsV3 import *

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
errMesh = (True if 'errMesh' in args else False)

# -- directory naming
baseCaseDir = 'baseCase_flow'
outFolder = 'ZZ_cases_flow'
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
eps = 0.5           # material porosity
nu = 5.58622522e-05 # kinematic viscosity
EA = 90e3           # activation energy

# -- geometry
R = 0.01            # sphere radius
length1 = 15*R      # inlet <-> sphere centre
length2 = 45*R      # sphere centre <-> outlet
width = 15*R        # top|bottom wall <-> sphere centre

# -- list parameters [ORIGINAL]
ReLst = [10, 40, 80, 160]                       # Reynolds number
invLst = [round(Re*nu/2/R,4) for Re in ReLst]   # inlet velocity
thieleLst = [2, 6]                              # Thiele modulus
cellSizeLst = [0.4*R]                           # FV cell Size
tortLst = [0.5, 1.0, 2.5, 5.0]                  # tortuosity

# == ARCHIVED SETTINGS: 
# -- 17. 1. 2023: mesh independence tests
thieleLst = [2]
thieleLst = [6]
ReLst = [80]
invLst = [round(Re*nu/2/R,4) for Re in ReLst]
# cellSizeLst = [0.8*R, 0.4*R, 0.2*R]  # khyrm@WSL
cellSizeLst = [0.8*R, 0.7*R, 0.6*R, 0.5*R, 0.4*R, 0.3*R, 0.2*R]  # khyrm@multipede
tortLst = [1]

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
        changeInCaseFolders(meshDir,'system/snappyHexMeshDictIntraTrans',['spR'],[str(R)])
        os.chdir(meshDir)
        os.system('chmod u=rwx All*') # NOTE MK: Just to make sure.
        os.system('./Allmesh')
        os.chdir('../../../')
    if not runSim: sys.exit()

# -- numpy array for results:
if errMesh:
    emdNp = np.zeros((2,len(cellSizeLst)))

# -- create cases for:
cases = [(inv,cellSize,tort,thiele) for inv in invLst for cellSize in cellSizeLst for tort in tortLst for thiele in thieleLst]
for case in cases:
    # -- parameters
    inv,cellSize,tort,thiele = case
    Re = ReLst[invLst.index(inv)]
    DFree = DFreeZ
    DEff = DFree*eps/tort
    k0Art = DEff*(thiele/R)**2 
    k0 = k0Art*np.exp(EA/(Runiv*T))
    # thiele tort re cS
    caseName = 'flow_phi_%g_tort_%g_Re_%g_cS_%g'%(thiele,tort,Re,cellSize)
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

    if errMesh:
        etaErr = abs(eta_sim-eta_corr)
        emdNp[:,cellSizeLst.index(cellSize)] = np.array([cellSize, etaErr])

    log_report(thiele,DEff,DFree,Re,Sh_corr,Sc,Sh,eta_sim,eta_corr,gradCCO,cCO,j,km,gradYCO,yCO,jY,kmY,ShY)
    os.chdir('../../')
    flow_csv(ZZZ_path,ZZZ_filepath,thiele,tort,Re,eta_sim,eta_corr)


if getCsv:
    # -- create eta 
    generate_eta_csvs(ZZZ_path,ZZZ_filepath,thieleLst,tortLst)
    # -- create eta csv
    for thiele in thieleLst:
        for tort in tortLst:
            # generate arrays
            n = 60
            ReNp = np.array([ReLst[0]+i*(ReLst[-1]-ReLst[0])/n for i in range(n+1)])
            DFree,DEff = DFreeZ,DFree*eps/tort
            Sc = nu/DFree
            Sh_corrNp = 2+.6*Sc**(1/3)*ReNp**(1/2)
            BiM_corrNp = Sh_corrNp/2 * DFree/DEff
            k0Art = DEff*(thiele/R)**2
            eta_corrNp = np.array(3/(thiele**2) * (thiele/np.tanh(thiele)-1)/(1+(thiele/np.tanh(thiele)-1)/BiM_corrNp))
            # write to CSVs
            with open(ZZZ_path+'/etacsv/etaCorr_phi_%g_tort_%2.1f.csv'%(thiele,tort), 'w') as f1:
                f1.writelines(['x, y\n'])
                f1.writelines(['%g, %g\n'%(ReNp[i], eta_corrNp[i]) for i in range(len(ReNp))])

if showPlots:
    for tort in tortLst:
        eta_plt(ZZZ_filepath, thiele, tort)
        # eta_err_plt(ZZZ_filepath, tortLst, invLst)
if errMesh:
    print(emdNp)
    # NOTE MK: This is not controled by getCsv.
    if not os.path.exists(ZZZ_path+'/errMeshcsv'):
        os.makedirs(ZZZ_path+'/errMeshcsv')
    with open(ZZZ_path+'/errMeshcsv/errMesh_phi_%g_Re_%g'%(thiele,round(Re)), 'w') as f1:
        f1.writelines(['cS, \terr\n'])
        for i in range(len(emdNp[0])):
            f1.writelines(['%g,\t%g\n'%(emdNp[0,i], emdNp[1,i])])
    if showPlots:
        title = 'Dependence of error on the mesh for φ = %g, Re = %g.'%(thiele,Re)
        # -- centred slopes
        at = -1  # crosspoint at
        plt.plot(np.array(cellSizeLst), np.array(cellSizeLst)/cellSizeLst[at]*emdNp[1,at], label='slope = 1')
        plt.plot(np.array(cellSizeLst), np.array(cellSizeLst)**2/cellSizeLst[at]**2*emdNp[1,at], label='slope = 2')
        plt.plot(np.array(cellSizeLst), emdNp[1], marker='x', linestyle='--', label='absolute η error', color='black')
        
        # -- uncentered slopes
        # plt.plot(cellSizeLst, cellSizeLst, label='slope = 1')
        # plt.plot(cellSizeLst, np.array(cellSizeLst)**2, label='slope = 2')
        # plt.plot(cellSizeLst, emdNp[1], marker='x', linestyle='--', label='absolute η error', color='black')
        
        plt.yscale('log')
        plt.xscale('log')
        plt.title(title)
        plt.legend()
        plt.show()