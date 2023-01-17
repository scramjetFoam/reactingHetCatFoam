# -- verInnerTransControlV2.py
# -- Script for verification cases with heat & mass transfer in a spherical particle
# -- NOTE: SCRIPT USES ARGUMENTS:
#       -- "makeMesh"
#       -- "runSim"
#       -- "showPlots"
#       -- "notParallel" (parallel by default)
# -- TODO: 
#       -- check order of the method

# -- imports
import numpy as np
import os
# import re
import shutil as sh
import matplotlib.pyplot as plt
import sys
from auxiliarFuncs import *
nonisoT_etaAnal = 42.094 # -- obtained from shootChandraVM.py

# -- set solver to be used
solverLst = ['reactingHetCatSimpleFoam']
solver = solverLst[0]

# -- set number of the enthalpy corrections
numOfTCorr = 1
isothermal = (True if numOfTCorr == 0 else False)  # isothermal logic

# -- script arguments logic
args = sys.argv
runSim = (True if 'runSim' in args else False)
showPlots = (True if 'showPlots' in args else False)
makeMesh = (True if 'makeMesh' in args else False)
getCsv = (True if 'getCsv' in args else False)
errMesh = (True if 'errMesh' in args else False)
# -- not sure why
parallel = (False if 'notParallel' in args else True)

# -- directory naming
baseCaseDir = 'baseCase'
outFolder = 'ZZ_cases'
if isothermal: outFolder += '_isoT'
ZZZ_path = 'ZZZ_res'
if isothermal: ZZZ_path += '_isoT'
# ZZZ_file = 'mult_steadySt'
# ZZZ_filepath = ZZZ_path+'/'+ZZZ_file

# -- set case parameters
yInf = 0.02      # molar fraction farfar from sphere
p = 101325      # presure	
Runiv = 8.314   # universal gas constant
R = 0.1           # sphere radius
T0 = 800     # used for multiple steady state, comment out or set to False if not used

# -- set geometry
domainSize = 1.1*R
tort = 2                # tortuosity
kappaEff = 2            # mass transfer coefficient
DFreeZ = 1e-5           # set diffusivity in fluid

# -- list parameters
thieleLst = [0.5]   # Thiele modulus
TLst = [300]
gammaLst = [20]
betaLst = [0.6]
# cellSizeLst = [0.4*R]  # NOTE: The mesh will be much more refined inside the sphere: (5 5)
# cellSizeLst = [0.2*R, 0.4*R, 0.8*R]
cellSizeLst = [0.4*R]

# -- prepare prototype mesh for each cellSize
if makeMesh:
    for cellSize in cellSizeLst:
        if not os.path.isdir('%s'%outFolder): os.mkdir('%s'%outFolder)
        meshDir = '%s/protoMesh/%g'%(outFolder,cellSize)
        print('Preparing mesh %s',meshDir)
        # -- check that meshDir is clean
        if os.path.isdir(meshDir): sh.rmtree(meshDir)
        # -- copy files
        sh.copytree(baseCaseDir,meshDir)
        changeInCaseFolders(meshDir,'system/blockMeshDict',['dSN', 'nDisc'],[str(domainSize),str(int(domainSize/cellSize*2))])
        changeInCaseFolders(meshDir,'system/snappyHexMeshDict',['spR'],[str(R)])
        changeInCaseFolders(meshDir,'system/snappyHexMeshDictIntraTrans',['spR'],[str(R)])
        os.chdir(meshDir)
        os.system('chmod u=rwx All*') # Just to make sure.
        os.system('./AllmeshIntraSphere')
        os.chdir('../../../')
    if not runSim:
        sys.exit()

# -- numpy array for results
if isothermal:
    resNp = np.zeros((len(cellSizeLst)))
    resNp2 = np.zeros((len(cellSizeLst)))
else:
    resNp = np.zeros((2,len(thieleLst)+1))
    if errMesh:
        emdNp = np.zeros((2,len(cellSizeLst)))

# -- create cases for
cases = [(T,thiele,cellSize,beta,gamma) for T in TLst for thiele in thieleLst for cellSize in cellSizeLst for beta in betaLst for gamma in gammaLst ]
for case in cases:
    # -- parameters
    T,thiele,cellSize,beta,gamma = case
    EA = gamma*Runiv*T
    DFree = DFreeZ
    DEff = DFree/tort*0.5
    sHr = -beta/yInf/p*Runiv*T/DEff*kappaEff*T
    k0 = (thiele/R)**2 * DEff/(np.exp(-gamma))
    if not T0: T0 = T

    caseName = 'intraTrans_phi_%g_beta_%g_cellSize_%g_T_%g'%(thiele,beta,cellSize,T)
    caseDir = '%s/%s/'%(outFolder,caseName)
    meshDir = '%s/protoMesh/%g'%(outFolder,cellSize)

    if runSim:
        print('Preparing case %s'%caseName)
        # -- check that caseDir is clean
        if os.path.isdir(caseDir): sh.rmtree(caseDir)
        # -- copy files
        sh.copytree(meshDir,caseDir)
        # -- change the case files
        # changeInCaseFolders(caseDir,'0.org/T',['isoT'],[str(T)]) # NOTE MK: modified to capture multiple steady states
        changeInCaseFolders(caseDir,'0.org/T',['initT','boundT'],[str(T0),str(T)])
        changeInCaseFolders(caseDir,'0.org/CO',['yCOSet'],[str(yInf)])
        changeInCaseFolders(caseDir,'system/controlDict',['customSolver'],[solver])
        changeInCaseFolders(caseDir,'system/fvSolution',['nTCorr'],[str(numOfTCorr)])
        changeInCaseFolders(caseDir,'constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(k0),str(EA),str(sHr)])
        changeInCaseFolders(caseDir,'constant/transportProperties',['kappaEffSet','tortSet','DSet'],[str(kappaEff),str(tort),str(DFreeZ)])
        # -- run simulation
        os.chdir(caseDir)
        os.system('chmod u=rwx All*')
        if parallel: 
            os.system('./AllrunIntraSphere-parallel')
        else:
            os.system('./AllrunIntraSphere')

    else:
        os.chdir(caseDir)
    
    k = k0*np.exp(-gamma)       # ???
    k0Art = k0*np.exp(-gamma)   # ???
    rSqIdeal = 4/3*np.pi*R**3*k0Art*yInf*p/Runiv/T     # ideal reaction source
    rS = read_real_source()                            # simulation reaction source
    etaSim = rS/rSqIdeal                               # simulation effectivness factor

    print('Case with thiele = %g'%thiele)
    
    if isothermal: 
        # -- analytical effectivness factor
        etaAnal = 3./(thiele**2)*(thiele*1./np.tanh(thiele)-1)
        print('Reaction source = %g, simulation effectivness factor is %g, relative error = %g'%(rS,etaSim,(etaAnal-etaSim)/etaAnal))
        print('\nThiele modulus = %g\nAnalytical effectivness factor is %g'%(thiele,etaAnal))   
        # -- load concetration profile to compare with analytical solution
        timeLst = os.listdir('./')
        timeLst = sorted([time for time in timeLst if isFloat(time)])
        for timeInd in range(len(timeLst)):
            with open('postProcessing/graphCellFace(start=(000),end=(100),fields=(CO))/%s/line.xy'%timeLst[timeInd],'r') as fl:
                lines = fl.readlines()
            simData = np.zeros((len(lines)-1,2))
            for lineInd in range(len(lines)-1):
                simData[lineInd,:] = np.array([num for num in np.array(lines[lineInd+1].replace('\n','').split(' ')) if isFloat(num)]).astype(float)
            plt.plot(simData[:,0],simData[:,1],label='sim, time = %s'%timeLst[timeInd])

        # -- analytical solution
        rAnal = simData[:,0]
        yAnal = yInf * R/rAnal * np.sinh(rAnal*(k/DEff)**0.5)/np.sinh(thiele)
        etaErr = abs(etaAnal-etaSim)
        etaErrRel = etaErr/etaAnal
        resNp[cellSizeLst.index(cellSize)] = etaErr
        resNp2[cellSizeLst.index(cellSize)] = np.linalg.norm(yAnal-simData[:,1])
        
        with open('anal.csv','w') as fl:
            fl.writelines('r,yinf\n')
            fl.writelines(['%g,%g\n'%(rAnal[i],yAnal[i]) for i in range(len(rAnal))])
        plt.plot(rAnal,yAnal,label='analytical')
        plt.legend()
        plt.title('Concentration profile comparison for k0=%d'%k0)
        plt.savefig('resHere.png')
        if showPlots:
            plt.show()
        if runSim or getCsv:
            id_parameters = (T,thiele)
            mesh_err_csv(ZZZ_path, id_parameters, cellSize, etaErrRel)

    else: # nonisothermal
        print('beta',beta,'thiele',thiele,'etaSim',etaSim,'R',R,'k0Art',k0Art,'Deff',DEff)
        resNp[:,thieleLst.index(thiele)] = np.array([thiele,etaSim])
        if errMesh:
            etaErr = abs(etaSim-nonisoT_etaAnal)
            emdNp[:,cellSizeLst.index(cellSize)] = np.array([cellSize, etaErr])

        # -- writing to a .csv
        if runSim or getCsv:
            # -- beta, thiele, etaSim
            csv_path = '../../%s/etaCsv'%ZZZ_path
            csv_file = '%s/simRes%g.csv'%(csv_path,beta)
            # -- ensure that the folder and the file (with a header) exist before writing
            if not os.path.exists(csv_path): 
                os.mkdir(csv_path)
            if not os.path.isfile(csv_file): 
                with open(csv_file, 'w') as f2: 
                    f2.writelines(['x,y\n'])
            # -- write
            # -- TODO: Only write new values
            with open(csv_file, 'a') as f3: f3.writelines(['%s,%s\n'%(str(thiele),str(etaSim))])


    os.chdir('../../')


# -- create plots
if showPlots:
    if isothermal:
        # -- error mesh dependence plot
        plt.plot(cellSizeLst,resNp,label='eta diff (my)')
        # plt.plot(cellSizeLst,resNp[1],label='eta diff (STF)')
        plt.plot(cellSizeLst,resNp2,label='whole sol diff (my)')
        # plt.plot(cellSizeLst,resNp2[1],label='whole sol diff (STF)')
        plt.plot(cellSizeLst,cellSizeLst,label='slope = 1')
        plt.plot(cellSizeLst,np.array(cellSizeLst)**2,label='slope = 2')
        title = 'Dependence of error on the mesh.'
        fileName = 'error_mesh.png'
        plt.title(title)
        plt.xlabel('FV cell size')
        plt.ylabel('error')
        plt.yscale('log')
        plt.xscale('log')
        plt.legend()
        # plt.savefig('%s/%s'%(ZZZ_path,fileName))
        plt.show()
    else:
        # -- eta-thiele diagram
        with open('baseCase/analRes06.csv', 'r') as fl:
        # with open('etaAnal_beta_%g_gamma_%g.csv'%(beta,gamma),'r') as fl:
            lines = fl.readlines()
        analRes = np.zeros((len(lines)-1,2))
        for i in range(1,len(lines)):
            analRes[i-1,:] = np.array(lines[i].split(',')) 
        plt.plot(analRes[:,0], analRes[:,1],'r',label='anal. res')
        plt.plot(resNp[0,:],resNp[1,:],'x',label='sim res.')
        plt.xlim((1e-1,1e1))
        plt.ylim((1e0,1e2))
        title = 'Multiple steady states.'
        fileName = 'mult_steadySt.png'
        plt.title(title)
        plt.yscale('log')
        plt.xscale('log')
        plt.legend()
        # plt.savefig('%s/%s'%(ZZZ_path,fileName))
        plt.show()
        # -- error mesh dependence plot
        print(emdNp)
        plt.plot(cellSizeLst, emdNp[1], marker='x', label='eta_err')
        plt.plot(cellSizeLst, cellSizeLst, label='slope = 1')
        plt.plot(cellSizeLst, np.array(cellSizeLst)**2, label='slope = 2')
        title = 'Dependence of error on the mesh.'
        plt.title(title)
        plt.yscale('log')
        plt.xscale('log')
        plt.legend()
        plt.show()
 
 