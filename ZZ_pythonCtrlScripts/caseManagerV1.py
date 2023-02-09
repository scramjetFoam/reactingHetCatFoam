# -- Python script to create manage openFoam cases using python class OpenFOAMCase

# -- for the usage see please OF_caseClass.py

# -- imports 
from OF_caseClass import OpenFOAMCase
import sys
import numpy as np
import os
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

# -- baseCase directory 
verTestDir = "../tutorials/verification/flowAroundSphereMass"
# -- directory naming
baseCaseDir = '%s/baseCase'%verTestDir
outFolder = '%s/ZZ_cases'%verTestDir
ZZZ_path = '%s/ZZZ_res'%verTestDir
ZZZ_file = '/flow.csv'
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

# -- baseCase object
# bsCase = OpenFOAMCase()
# bsCase.loadOFCaseFromBaseCase(baseCaseDir)

# -- different cellSize meshes baseCase list
meshesCaseLst = []

# -- prepare prototype mesh for each cellSize
if makeMesh:
    for cellSize in cellSizeLst:
        meshCase = OpenFOAMCase()
        meshCase.loadOFCaseFromBaseCase(baseCaseDir)
        meshDir = '%s/protoMesh/%g'%(outFolder,cellSize)
        meshCase.changeOFCaseDir(meshDir)
        meshCase.copyBaseCase()
        replaceInBlockMesh = ["system/blockMeshDict",['length1', 'length2', 'width','nDiscX','nDiscYZ'],[str(length1),str(length2),str(width),str(int((length1+length2)/cellSize)),str(int(2*width/cellSize))]]
        replaceInSnappyHexMesh = ["system/snappyHexMeshDict", ['spR'], [str(R)]]             
        meshCase.replace([replaceInBlockMesh, replaceInSnappyHexMesh])
        meshCase.runCommands(['chmod u=rwx All*', './Allmesh'])
        meshesCaseLst.append(meshCase)
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
        caseHere = OpenFOAMCase()
        caseHere.loadOFCaseFromBaseCase(meshDir)
        caseHere.changeOFCaseDir(caseDir)
        meshCase.copyBaseCase()
        caseHere.replace(caseDir,'0.org/T',['isoT'],[str(T)])
        caseHere.replace(caseDir,'0.org/CO',['yCOSet'],[str(yInf)])
        caseHere.replace(caseDir,'0.org/U', ['inv'],[str(inv)])
        caseHere.replace(caseDir,'system/controlDict',['customSolver'],[solver])
        caseHere.replace(caseDir,'system/fvSolution',['nTCorr'],[str(numOfTCorr)])
        caseHere.replace(caseDir,'constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(k0),str(EA),str(sHr)])
        caseHere.replace(caseDir,'constant/transportProperties',['kappaEffSet','tortSet','DSet'],[str(kappaEff),str(tort),str(DFreeZ)])
        # -- run simulation
        caseHere.runCommands(['chmod u=rwx Allrun-parallel', './Allrun-parallel'])
    else: 
        caseHere = OpenFOAMCase()
        caseHere.loadOFCaseFromBaseCase(meshDir)
    
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
        title = 'Dependence of error on the mesh for φ = %g, Re = %g, Deff/DFree = %g.'%(thiele,Re, DEff/DFree)        # -- numerical slope fit
       # -- numerical slope fit
        try:
            fit, slope = logfit(emdNp)
        except NameError:
            fit = emdNp[1]
            slope = 0
            print('Warning: missing scipy')
        # -- centerinhg
        at = -1
        # at = False
        if at: b, c, d = emdNp[0,at], emdNp[1,at], fit[at]
        else:  b, c, d = 1, 1, 1
        # plt.plot(emdNp[0], emdNp[1]/c, marker='x', linewidth=0, label='absolute η error', color='black')
        plt.plot(emdNp[0], fit/d, linestyle='--', label='slope fit = %3.2f'%slope, color='black')
        plt.plot(emdNp[0], (emdNp[0])/b, label='slope = 1')
        plt.plot(emdNp[0], (emdNp[0]**2)/(b**2), label='slope = 2')
        
        plt.yscale('log')
        plt.xscale('log')
        plt.title(title)
        plt.legend()
        plt.show()
        
     