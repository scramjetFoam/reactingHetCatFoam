# -- Script to do verification simulation of isothermal reaction inside sphere

# -- TODO:
# 1) try different parameters setups, check analytical profile in the isothermal cases, check order of the method (should be 2) in all cases
# 2) MK: script arguments for numOfTCorr, dynamic & runSim

# -- imports
import numpy as np
import os 
import re
import shutil as sh
import matplotlib.pyplot as plt

# -- auxiliary function to determine if val is float
def isFloat(val):
    try:
        float(val)
        return True
    except:
        return False

# -- auxiliary function to change case params 
# NOTETH can be done faster, but I did not want to think about it
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
    
    # prepare & open short temporary log file:
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
# NOTE TH: if 0 - isothermal study
numOfTCorr = 0
# isothermal logic:
isothermal = False
if numOfTCorr == 0: 
    isothermal = True

# -- swicher if the simulation should be run or to just check results
runSim = True  
# runSim = False

# -- case parameters
yInf = 0.1      # molar fraction farfar from sphere
p = 101325      # presure
sHr = -138725599.72220203    # standard reaction enthalpy (physical 283e3)	
Runiv = 8.314   # universal gas constant
R = 1           # sphere radius
Rinf = 1.1      # infinite radius
inv = 0         # inlet velocity
# dynamic logic:
dynamic = False
if inv != 0: 
    dynamic = True  
# tube dimensions:
length1 = 1.1        # inlet <-> sphere centre
length2 = 3*length1  # sphere centre <-> outlet
width = 1.1          # top|bottom wall <-> sphere centre

# -- setup study parameters here
baseCaseDir = 'baseCase'
outFolder = 'ZZ_cases'

# MK: change naming for dynamic
if dynamic:
    baseCaseDir += '_dynamic'
    outFolder += '_dyn'
if isothermal:
    outFolder += '_isoT'
  
cellSizeLst = [0.1]  # FV cell Size
# cellSizeLst = [0.5,0.25,0.125,0.0625,0.03125]  
# cellSizeLst = [0.5,0.25,0.125,0.0625]
k0Lst = [1e2, 5e2, 1e3, 5e3, 1e4, 1e5]      # reaction pre-exponential factor
# EA = 90e3     # reaction activation energy (set according to gamma) 2020 Chandra Direct numerical simulation of a noniso...
TLst = [500]   # temperature (set according to gamma) 2020 Chandra Direct numerical simulation of a noniso...
gamma = 20      # see line above
# betaLst = [0.4,0.6,0.8] set according to T) 2020 Chandra Direct numerical simulation of a noniso...
tort = 1      # tortuosity
solverLst = ['reactingHetCatSimpleFoam','scalarTransportFoamCO'] # used solver
solver = solverLst[0]
kappaEff = 1

# -- numpy array with results
if isothermal:
    resNp = np.zeros((len(cellSizeLst)))
    resNp2 = np.zeros((len(cellSizeLst)))
else:
    resNp = np.zeros((2,len(k0Lst)+1))

# -- create case for:
# -- 1) varying T
# -- 2) varying k0
# -- 3) varysing cS
for TInd in range(len(TLst)):
    for k0Ind in range(len(k0Lst)):
        for cellSizeInd in range(len(cellSizeLst)):
            cellSize = cellSizeLst[cellSizeInd]
            T = TLst[TInd]
            k0 = k0Lst[k0Ind]
            EA = gamma*Runiv*T

            # -- create caseFolder based on baseCase
            caseName = 'intraTrasn_yInf_%g_R_%g_T_%g_cS_%g_k0_%g_tort_%g'%(yInf,R,T,cellSize,k0,tort)
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
                if dynamic: changeInCaseFolders('system/blockMeshDict',['length1', 'length2', 'width','nDisc'],[str(length1),str(length2),str(width),str(int(Rinf/cellSize)*2)])
                else: changeInCaseFolders('system/blockMeshDict',['dSN', 'nDisc'],[str(Rinf),str(int(Rinf/cellSize)*2)])
                changeInCaseFolders('system/fvSolution',['nTCorr'],[str(numOfTCorr)])
                changeInCaseFolders('constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(k0),str(EA),str(sHr)])
                changeInCaseFolders('constant/transportProperties',['kappaEffSet','tortSet'],[str(kappaEff),str(tort)])
                
                # -- run the simulation
                os.chdir(caseDir)
                os.system('./AllrunIntraSphere')

            else:
                os.chdir(caseDir)

            # -- load DFree, Deff, and k from log file
            pars = {'DFree':('DFree',-1), 'DEff':('DEff', 0), 'k':('max(k)', -1)}
            par_vals = pars_from_log(pars, solver)
            DFree, DEff, k = par_vals

            # -- compute analytical results
            k0Art = k0*np.exp(-EA/(Runiv*T))                            # definition in 2020 Chandra Direct numerical simulation of a noniso...
            rSqIdeal = 4./3*np.pi*R**3*k0Art*yInf*p/Runiv/T             # ideal reaction source
            thiele = R*(k0Art/DEff)**(0.5)                              # thiele modulus
            
            # -- compute simulation results
            # NOTE MK: use my own function
            # -- read real source
            with open('log.intSrcSphere', 'r') as fl:
                lines = fl.readlines()
            # reaction source
            for lineInd in range(len(lines)-1,0,-1):
                if lines[lineInd].find('reaction source') >= 0:
                    rS = float(lines[lineInd].split(' ')[-1].replace('\n',''))
                    break
                if lineInd == 0:
                    print('Reaction source not found.')
            print('reaction source = %g'%rS)
            
            etaSim = rS/rSqIdeal                                        # simulation effectivness factor
            
            # -- in the isothermal cases we can compare concentration profiles
            if isothermal: 
                etaAnal = 3./(thiele**2)*(thiele*1./np.tanh(thiele)-1)  # analytical effectivness factor
                print('Simulation effectivness factor is %g, relative error = %g'%(etaSim,(etaAnal-etaSim)/etaAnal))
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
                plt.savefig('resHere.png')
                plt.show()
            
            # -- in the nonisothermal case we can compare eta-thiele diagram
            else:
                beta = yInf*p/Runiv/T*(-sHr)*DEff/kappaEff/T
                print('beta',beta,'thiele',thiele,'effect sim',etaSim)
                resNp[:,k0Ind] = np.array([thiele,etaSim])
                

            # NOTETH: inner + outer transport
            # BiM = kL*R/DEff             # Biot number
            # etaAnal = 3./(thiele**2) * (thiele*1./np.tanh(thiele)-1)/(1+(thiele*1./np.tanh(thiele)-1)/BiM)                      # analytical effectivness factor
            
            os.chdir('../../')


## -- create plot

# create &or name directory for res plot:
dirName = 'ZZZ_res'
if dynamic:    dirName += '_dyn'
if isothermal: dirName += '_isoT'
if not os.path.exists(dirName): os.mkdir(dirName)

# create plot:
if isothermal:  # for isothermal cases, dependence of error on mesh
    plt.plot(cellSizeLst,resNp,label='eta diff (my)')
    # plt.plot(cellSizeLst,resNp[1],label='eta diff (STF)')
    plt.plot(cellSizeLst,resNp2,label='whole sol diff (my)')
    # plt.plot(cellSizeLst,resNp2[1],label='whole sol diff (STF)')
    plt.plot(cellSizeLst,cellSizeLst,label='slope = 1')
    plt.plot(cellSizeLst,np.array(cellSizeLst)**2,label='slope = 2')
    title = 'Dependence of the error on the mesh.'
    fileName = 'error_mesh.png'
else:
    plt.plot(resNp[0,:],resNp[1,:],'x',label='sim res.')
    title = 'Dependence of the error on the mesh.'
    fileName = 'mult_steadySt.png'
plt.title(title)
plt.savefig('%s/%s'%(dirName,fileName))
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.show()