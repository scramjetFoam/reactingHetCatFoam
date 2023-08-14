# -- verInnerTransControlMassV1.py
# -- Script for verification cases with heat & mass transfer in a spherical particle
# -- NOTE SCRIPT ARGUMENTS:
#   -- "makeMesh"
#   -- "runSim"
#   -- "showPlots"
#   -- "errMesh"

# -- imports
import numpy as np
import os
import shutil as sh
import matplotlib.pyplot as plt
import sys
sys.path.append('../../../ZZ_pythonCtrlScripts')
from OF_caseClass import OpenFOAMCase
from auxiliarFuncsMassV1 import *

# -- obtained from shootChandraVM.py
# nonisoT_etaAnal = 42.094    # phi = 0.5
nonisoT_etaAnal =  7.516    # phi = 4.0

# -- set solver to be used
solverLst = ['reactingHetCatSimpleFoamM']
solver = solverLst[0]

# -- simulation parameters
endTime = 500
nProc = 12
numOfCCorr = 1
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

# -- directory naming
baseDir = '../baseCaseMass_local_intraTransp'
outFolder = 'ZZ_cases'
if isothermal: outFolder += '_isoT'
ZZZ_path = 'ZZZ_res'
if isothermal: ZZZ_path += '_isoT'
# ZZZ_file = 'mult_steadySt'
# ZZZ_filepath = ZZZ_path+'/'+ZZZ_file
if not os.path.exists(ZZZ_path): os.mkdir(ZZZ_path)
if not os.path.exists(outFolder): os.mkdir(outFolder)

# -- set case parameters
yInf = 0.02      # molar fraction farfar from sphere
p = 101325      # presure	
Runiv = 8.314   # universal gas constant
R = 0.1           # sphere radius
T0 = False     # used for multiple steady state (800), comment out or set to False if not used
# T0 = 800

# -- set geometry
domainSize = 1.1*R
tort = 2                # tortuosity
kappaEff = 2            # mass transfer coefficient
DFreeZ = 1e-5           # set diffusivity in fluid

# -- list parameters [ORIGINAL]
thieleLst = [0.2, 0.4, 0.5, 0.75, 1.0, 2.0, 4.0]
# thieleLst = [0.4, 0.5, 0.75, 1.0, 2.0, 4.0]
# thieleLst = [0.5]
TLst = [300]
gammaLst = [20]
betaLst = [0.6]
# cellSizeLst = [0.4*R]  # NOTE: The mesh will be much more refined inside the sphere: (5 5)
cellSizeLst = [0.2*R]

# -- chemical species
specieNames = np.array(["CO", "prod", "N2"])
species = ' '.join(specieNames) 
molMass = np.array([ 28e-3, 28e-3, 28e-3])
sigmaVs = np.array([ 18.0, 18.0, 18.0])
nuVec = np.array([-1,1,0])
alphaVec = np.array([1, 0, 0])
yIn = np.array([yInf, 0, 1-yInf])
MgIn = np.sum(yIn * molMass)
wIn = yIn * molMass / MgIn

# -- mesh tests
# thieleLst = [0.5]
# thieleLst = [4.0]
# cellSizeLst = [1.6*R, 0.8*R, 0.4*R, 0.2*R, 0.1*R]
# cellSizeLst = [0.4*R, 0.2*R, 0.1*R]
# cellSizeLst = [1.6*R, 0.8*R]

# -- prepare prototype mesh for each cellSize
meshesCaseLst = []

if makeMesh:
    for cellSize in cellSizeLst:
        meshCase = OpenFOAMCase()
        meshCase.loadOFCaseFromBaseCase(baseDir)
        meshDir = '%s/protoMesh/%g'%(outFolder,cellSize)
        meshCase.changeOFCaseDir(meshDir)
        meshCase.copyBaseCase()

        # -- create and update boundary conditions
        for chSpI in range(len(specieNames)):
            name = specieNames[chSpI]
            meshCase.runCommands(['cp 0.org/bsChemSp 0.org/%sMass' % name])
            meshCase.addToDictionary([['0.org/%sMass' % name, '\n\t"(sphere|infinity)"\n\t{\n\t\ttype\t\t\tfixedValue;\n\t\tvalue\t\t\tuniform wChSpSet;\n\t}\n', 'boundaryField']])
            meshCase.replace([["0.org/%sMass"% (name), ['wChSpSetInit','wChSpSet','nameSet'], ['%.5g'%(wIn[chSpI]), '%.5g'%(wIn[chSpI]), str(name)]]])
        meshCase.addToDictionary([['0.org/U', '\n\t"(sphere|infinity)"\n\t{\n\t\ttype\t\t\tfixedValue;\n\t\tvalue\t\t\tuniform (0 0 0.0);\n\t}\n', 'boundaryField']])
        meshCase.addToDictionary([['0.org/T', '\n\t"(sphere|infinity)"\n\t{\n\t\ttype\t\t\tfixedValue;\n\t\tvalue\t\t\tuniform boundT;\n\t}\n', 'boundaryField']])
        meshCase.addToDictionary([['0.org/p', '\n\tinfinity\n\t{\n\t\ttype\t\t\ttotalPressure;\n\t\tp0\t\t\t\tuniform 101325;\n\t\tvalue\t\t\tuniform 101325;\n\n\t\tU\t\t\t\tU;\n\t\tphi\t\t\t\tphi;\n\t}\n', 'boundaryField']])
        
            
        meshCase.runCommands(['rm 0.org/bsChemSp'])
        
        # -- update transport properties
        meshCase.addToDictionary([['constant/transportProperties', 'species (%s);\n' % species, '']])

        for nameInd in range(len(specieNames)):
            name = specieNames[nameInd]
            meshCase.addToDictionary([
                ['constant/transportProperties', '%s\n{\n}\n' % name, ''],
                ['constant/transportProperties', 'D  D\t[0 2 -1 0 0 0 0] DSet;\n', name ],
                ['constant/transportProperties', 'sigmaV\t%g;\n' % sigmaVs[nameInd], name ],
                ['constant/transportProperties', 'molM\t%g;\n' % molMass[nameInd], name ],
                ['constant/transportProperties', 'nuVec\t(0 %g);\n' % nuVec[nameInd], name ],
                ['constant/transportProperties', 'alphaVec\t(0 %g);\n' % alphaVec[nameInd], name ],
            ])

        # -- update blockMeshDict
        meshCase.replace([
            ['system/blockMeshDict', ['length1', 'length2', 'width', 'width'], [str(domainSize), str(domainSize), str(domainSize), str(domainSize)]],
            ['system/blockMeshDict', ['nDiscX', 'nDiscYZ'], [str(int(domainSize/cellSize*2)), str(int(domainSize/cellSize*2))]]
        ])
        # -- update others
        meshCase.replace([
            ["system/decomposeParDict",['np'],[str(nProc)]],
            ["system/snappyHexMeshDict",['spR'],[str(R)]]
        ])
        meshCase.replace([
            ['0.org/U',['inv'],[str(0)]]
        ])
        # -- setup verification case from the base case
        # -- creation and updates of the boundary conditions
        meshCase.setParameters([
            ['system/controlDict', 'endTime', str(endTime), ''],
            ['system/decomposeParDict', 'numberOfSubdomains', str(nProc), ''],
            ['system/fvSolution', 'nConcCorrectors', str(numOfCCorr), ''],
            ['system/fvSolution', 'nTempCorrectors', str(numOfTCorr), ''],
            ['system/fvSolution', 'p', str(0.18), 'fields'],
            ['constant/transportProperties', 'D', str(DFreeZ), 'CO']
        ])

        meshCase.runCommands([
            'rm -rf 0',
            'mkdir 0', 
            'cp -rf 0.org/* 0', 
            'mkdir dynamicCode',
            'blockMesh > log.blockMesh',
            'paraFoam -touch',
            'snappyHexMesh -overwrite > log.snappyHexMesh',
            'topoSet -dict "../../../baseCase/system/topoSetDict" > log.topoSet',
            'rm -rf */cellLevel',
            'rm -rf */pointLevel',
        ])

        meshesCaseLst.append(meshCase)

    if not runSim: sys.exit()


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
        caseHere = OpenFOAMCase()
        caseHere.loadOFCaseFromBaseCase(meshDir)
        caseHere.changeOFCaseDir(caseDir)
        caseHere.copyBaseCase()
        caseHere.replace([
            ['0.org/T',['initT','boundT'],[str(T0),str(T)]],
            # ['0.org/CO',['yCOSet'],[str(yInf)]],
            ['system/controlDict',['customSolver'],[solver]],
            ['system/fvSolution',['customSolver'],['customSolver|%s' %species.replace(' ','Mass|')]],
            ['constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(k0),str(EA),str(sHr)]],
            ['constant/transportProperties',['kappaEffSet','tortSet','DSet'],[str(kappaEff),str(tort),str(DFreeZ)]]
        ])
        caseHere.runCommands([
            'rm -rf 0',
            'mkdir 0',
            'cp -rf 0.org/* 0',
            'decomposePar -force > log2.decomposePar',
            'foamJob -parallel -screen renumberMesh -overwrite > log.renumberMesh', 
            'foamJob -parallel -screen %s > log1.%s' % (solver, solver),
        ])

        # -- rerun if the first run fails
        rerun = False
        with open('%slog1.%s'%(caseDir,solver), 'r') as f:
            if (f.readlines()[-3][:-1] != 'End'): rerun = True
        if rerun: caseHere.runCommands(['foamJob -parallel -screen %s > log2.%s' % (solver, solver)])

        caseHere.runCommands([
            'reconstructPar -latestTime > log2.reconstructPar',
            'intSrcSphereM > log.intSrcSphereM',
            # 'postProcess -func integrace -latestTime > log.integrace',
            "postProcess -func 'graphCellFace(start = (0 0 0), end = (1 0 0), fields=(COMass))' > log.postProcess",
            'rm -rf processor*'
        ])

    else: 
        caseHere = OpenFOAMCase()
        caseHere.loadOFCaseFromBaseCase(caseDir)

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
if errMesh:
    # -- error mesh dependence plot
    print(emdNp)
    if not os.path.exists(ZZZ_path+'/errMeshcsv'):
        os.makedirs(ZZZ_path+'/errMeshcsv')
    with open(ZZZ_path+'/errMeshcsv/errMesh_phi_%3.2f_beta_%g_gamma_%g'%(thiele,beta,gamma), 'w') as f1:
        f1.writelines(['cS, \terr\n'])
        for i in range(len(emdNp[0])):
            f1.writelines(['%g,\t%g\n'%(emdNp[0,i], emdNp[1,i])])
    if showPlots:
        title = 'Dependence of error on the mesh for φ = %g, β = %g, γ = %g.'%(thiele,beta,gamma)
        # -- numerical slope fit
        try:
            fit, slope = logfit(emdNp)
            print('slope =', slope)
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
