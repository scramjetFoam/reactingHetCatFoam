# -- Script to do verification simulation of isothermal reaction inside sphere

# -- NOTE:
# (1) isothermal study (numOfTCorr = 0)
#     -> effectivness comparison (analytical x simulation)
#     -> concentration profile (analytical x simulation)
#     -> mesh independence study (for various cell sizes)
#     -> dependence of the error on the mesh
# (2) non-isothermal study (numOfTCorr > 0)
#     -> multiple steady states
#     -> dependence of the error on the mesh
# (3) flow study (inv > 0)
#     -> comparison of effectivness for different Biot numbers (Ds/Df ratios)
#     -> dependence of the error on the mesh
# (4) flow study let to separate file
#     -> parameters set to get curve from Chandra


# -- imports
import numpy as np
import os
import re
import shutil as sh
import matplotlib.pyplot as plt
import sys

print(sys.argv)

# -- auxiliary function to determine if val is float
def isFloat(val):
    try:
        float(val)
        return True
    except:
        return False

# -- auxiliary function to change case params 
def changeInCaseFolders(file,whatLst,forWhatLst):
    """Set params in files copied from baseCase."""
    print('Changing file %s, %s --> %s'%(file,str(whatLst),str(forWhatLst)))
    with open(caseDir + '/%s'%file, 'r') as fl:
        data = fl.readlines()

    for whatInd in range(len(whatLst)):
        for ind in range(len(data)):
            data[ind] = data[ind].replace(whatLst[whatInd],forWhatLst[whatInd])

    with open(caseDir + '/%s'%file, 'w') as fl:
        fl.writelines(data)

# -- auxilary function to read parameters from log files
def pars_from_log(pars, solver):
    """Read params given by a dictionary from a temporary log file."""
    print('Reading parameters from log file.')
    par_vals = []
    par_keys = list(pars.keys())
    
    # create, read & remove temporary log file:
    with open('tmplog.%s'%solver, 'w') as tmplog: tmplog.write(os.popen('tail -n 100 log.%s'%solver).read())
    with open('tmplog.%s'%solver, 'r') as tmplog: lines = tmplog.readlines() 
    for par_key in par_keys:
        for lineInd in range(len(lines)-1,0,-1):
            if lines[lineInd].find(pars[par_key][0]) >= 0:
                num_lst = [t[0] for t in re.findall("(\d+[./]\d*e?-?\d*)|(\d+[./]\d*)", lines[lineInd])]
                val = num_lst[pars[par_key][1]]
                par_vals.append(float(val))
                print('%s = %s'%(par_key, val))
                break
            if lineInd == 0:
                par_vals.append("N/A")
                print('%s not found'%pars[par_key][0])
    os.remove('tmplog.%s'%solver)
    return par_vals

# -- number of the enthalpy corrections
numOfTCorr = 1

# isothermal logic:
isothermal = False
if numOfTCorr == 0: 
    isothermal = True

# -- swicher if the simulation should be run or to just check results
runSim = True  
# runSim = False

# -- case parameters
yInf = 0.02      # molar fraction farfar from sphere
p = 101325      # presure	
Runiv = 8.314   # universal gas constant
R = 0.1           # sphere radius

# -- tube dimensions:
length1 = 1.1*R       # inlet <-> sphere centre
length2 = length1       # sphere centre <-> outlet
width = length1         # top|bottom wall <-> sphere centre

DFreeZ = 1e-5
inv = 0

# -- setup study parameters here
baseCaseDir = 'baseCase'
outFolder = 'ZZ_cases'
  
# -- NOTE: k0Lst not used 
thieleLst = [0.5,0.75,1.,2,4]

TLst = [300]   # temperature (set according to gamma) 2020 Chandra Direct numerical simulation of a noniso...
gamma = 20      # see line above
beta = 0.6

tort = 2      # tortuosity
solverLst = ['reactingHetCatSimpleFoam','scalarTransportFoamCO'] # used solver
solver = solverLst[0]
kappaEff = 2

# -- blockMesh cell dimension -- NOTE: in snappy much more refined (4 4) on sphere
cellSizeLst = [0.5*R]  # FV cell Size

# MK: change naming for flow
if isothermal:
    outFolder += '_isoT'

# -- numpy array with results
if isothermal:
    resNp = np.zeros((len(cellSizeLst)))
    resNp2 = np.zeros((len(cellSizeLst)))
else:
    resNp = np.zeros((2,len(thieleLst)+1))

# -- create case for:

# for tortInd in range(len(tortLst)):
for TInd in range(len(TLst)):
    for thieleInd in range(len(thieleLst)):
        for cellSizeInd in range(len(cellSizeLst)):
            cellSize = cellSizeLst[cellSizeInd]
            T = TLst[TInd]
            EA = gamma*Runiv*T
            # tort = tortLst[k0Ind]
            
            DEff = DFreeZ/tort*0.5  

            # -- set parameters for kinetics to ensure thiele, gamma, beta
            k0 = (thieleLst[thieleInd]/R)**2*DEff/(np.exp(-EA/(Runiv*T)))
            sHr = -beta/yInf/p*Runiv*T/DEff*kappaEff*T
            

            # -- create caseFolder based on baseCase
            caseName = 'intraTrans_phi_%g_beta_%g_cellSize_%g'%(thieleLst[thieleInd],beta,cellSize)
            caseDir = '%s/%s/'%(outFolder,caseName)

            if runSim:
                print('Preparing case %s'%caseName)
                if os.path.isdir(caseDir):  # ensure the caseDir is clear
                    sh.rmtree(caseDir)

                sh.copytree(baseCaseDir,caseDir)

                # -- change the case files
                changeInCaseFolders('0.org/T',['isoT'],[str(T)])
                changeInCaseFolders('0.org/CO',['yCOSet'],[str(yInf)])
                changeInCaseFolders('0.org/U', ['inv'],[str(inv)])
                changeInCaseFolders('system/controlDict',['customSolver'],[solver])
                changeInCaseFolders('system/blockMeshDict',['dSN', 'nDisc'],[str(length1),str(int(length1/cellSize*2))])
                changeInCaseFolders('system/fvSolution',['nTCorr'],[str(numOfTCorr)])
                changeInCaseFolders('constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(k0),str(EA),str(sHr)])
                changeInCaseFolders('constant/transportProperties',['kappaEffSet','tortSet','DSet'],[str(kappaEff),str(tort),str(DFreeZ)])
                changeInCaseFolders('system/snappyHexMeshDict',['spR'],[str(R)])
                changeInCaseFolders('system/snappyHexMeshDictIntraTrans',['spR'],[str(R)])
                
                # -- run the simulation
                os.chdir(caseDir)
                os.system('./AllrunIntraSphere')

            else:
                os.chdir(caseDir)

            # -- load DFree, Deff, and k from log file
            pars = {'DFree':('DFree',-1), 'DEff':('DEff', 0), 'k':('max(k)', -1)}
            # par_vals = pars_from_log(pars, solver)
            # DFree, DEff, k = par_vals
            DFree, k = DFreeZ, k0*np.exp(-EA/(Runiv*T)) 

            # -- compute case parameters
            k0Art = k0*np.exp(-EA/(Runiv*T))                            # definition in 2020 Chandra Direct numerical simulation of a noniso...
            thiele = R*np.sqrt(k0Art/DEff)
            rSqIdeal = 4./3*np.pi*R**3*k0Art*yInf*p/Runiv/T             # ideal reaction source
            print('Case with thiele = %g'%thiele)        

            
            # -- compute simulation results
            # NOTE MK: use custom function
            # -- read real source
            with open('log.intSrcSphere', 'r') as fl:
                lines = fl.readlines()
            # reaction source
            for lineInd in range(len(lines)-1,0,-1):
                if lines[lineInd].find('reaction source') >= 0:
                    rS = float(lines[lineInd].split(' ')[-1].replace('\n',''))
                    # -- simulation effectivness factor
                    etaSim = rS/rSqIdeal  
                    break
                if lineInd == 0:
                    print('Reaction source not found.')
            
            # -- in the isothermal cases we can compare concentration profiles
            if isothermal: 
                # -- analytical effectivness factor
                etaAnal = 3./(thiele**2)*(thiele*1./np.tanh(thiele)-1) 
                print('Reaction source = %g, simulation effectivness factor is %g, relative error = %g'%(rS,etaSim,(etaAnal-etaSim)/etaAnal))

                print('\nThiele modulus = %g\nAnalytical effectivness factor is %g'%(thiele,etaAnal))   
                # -- load concetration profile to compare with analytical solution
                timeLst = os.listdir('./')
                timeLst = sorted([time for time in timeLst if isFloat(time)])
                print(timeLst)
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
                resNp[cellSizeInd] = abs(etaAnal-etaSim)
                resNp2[cellSizeInd] = np.linalg.norm(yAnal-simData[:,1])
                with open('anal.csv','w')as fl:
                    fl.writelines('r,yinf\n')
                    for i in range(len(rAnal)):
                        fl.writelines('%g,%g\n'%(rAnal[i],yAnal[i]))
                plt.plot(rAnal,yAnal,label='analytical')
                plt.legend()
                plt.title('Concentration profile comparison for k0=%d'%k0)
                plt.savefig('resHere.png')
                plt.show()
            
            # -- in the nonisothermal case we can compare eta-thiele diagram
            elif not isothermal:
                etaSim = rS/rSqIdeal  # simulation effectivness factor
                print('beta',beta,'thiele',thiele,'etaSim',etaSim,'R',R,'k0Art',k0Art,'Deff',DEff)
                resNp[:,thieleInd] = np.array([thiele,etaSim])
            os.chdir('../../')

# create &or name directory for res plot:
dirName = 'ZZZ_res'
if isothermal: dirName += '_isoT'
if not os.path.exists(dirName): os.mkdir(dirName)

# create plot:
if isothermal:
    plt.plot(cellSizeLst,resNp,label='eta diff (my)')
    # plt.plot(cellSizeLst,resNp[1],label='eta diff (STF)')
    plt.plot(cellSizeLst,resNp2,label='whole sol diff (my)')
    # plt.plot(cellSizeLst,resNp2[1],label='whole sol diff (STF)')
    plt.plot(cellSizeLst,cellSizeLst,label='slope = 1')
    plt.plot(cellSizeLst,np.array(cellSizeLst)**2,label='slope = 2')
    title = 'Dependence of the error on the mesh.'
    fileName = 'error_mesh.png'
else:
    with open('baseCase/analRes06.csv' ,'r') as fl:
        lines = fl.readlines()
    analRes = np.zeros((len(lines)-1,2))
    for i in range(1,len(lines)):
        analRes[i-1,:] = np.array(lines[i].split(',')) 
    # plt.plot(resNp[0,:],resNp[1,:],'x',label='sim res.')
    plt.plot(analRes[:,0], analRes[:,1],'r',label='anal. res')
    plt.plot(resNp[0,:],resNp[1,:],'x',label='sim res.')
    # plt.plot(analRes[:,0], analRes[:,1])#,resNp[0,:],resNp[1,:],'x',label='sim res.')
    title = 'Multiple steady states.'
    fileName = 'mult_steadySt.png'
plt.title(title)
plt.yscale('log')
plt.xscale('log')
plt.legend()
# plt.savefig('%s/%s'%(dirName,fileName))
plt.show()