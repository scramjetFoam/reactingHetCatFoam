# -- Script to do verification simulation of isothermal reaction inside sphere

# -- TODO: try different parameters setups, check analytical profile in the isothermal cases, check order of the method (should be 2) in all cases

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
numOfTCorr = 0
# isothermal logic:
isothermal = False
if numOfTCorr == 0: 
    isothermal = True

# -- swicher if the simulation should be run or to just check results
runSim = True  
# runSim = False

# -- case parameters
yInf = 0.01      # molar fraction farfar from sphere
p = 101325      # presure
sHr = -283e3    # standard reaction enthalpy (physical 283e3)	
Runiv = 8.314   # universal gas constant
R = 0.01           # sphere radius
Rinf = 1.1      # infinite radius
# inv = 0.0558        # inlet velocity
inv = 0.1116        # inlet velocity
# flow logic:
flow = False
if inv != 0: 
    flow = True  
# tube dimensions:
length1 = 15*R       # inlet <-> sphere centre
length2 = 45*R       # sphere centre <-> outlet
width = 15*R         # top|bottom wall <-> sphere centre
DFreeZ = 1e-5

# -- setup study parameters here
baseCaseDir = 'baseCase'
outFolder = 'ZZ_cases'

# MK: change naming for flow
if flow:
    baseCaseDir += '_flow'
    outFolder += '_flow'
if isothermal:
    outFolder += '_isoT'
  
# cellSizeLst = [R/3]  # FV cell Size
cellSizeLst = [0.7*R]  # FV cell Size
# cellSizeLst = [0.5,0.25,0.125,0.0625,0.03125]  
# cellSizeLst = [0.5,0.25,0.125,0.0625]  # MK
#k0Lst = [1e2, 5e2, 1e3, 5e3, 1e4, 1e5]      # reaction pre-exponential factor
k0Lst = [1e9]
# EA = 90e3     # reaction activation energy (set according to gamma) 2020 Chandra Direct numerical simulation of a noniso...
TLst = [500]   # temperature (set according to gamma) 2020 Chandra Direct numerical simulation of a noniso...
gamma = 20      # see line above
# betaLst = [0.4,0.6,0.8] set according to T) 2020 Chandra Direct numerical simulation of a noniso...
tort = 5      # tortuosity
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

for TInd in range(len(TLst)):
    for k0Ind in range(len(k0Lst)):
        for cellSizeInd in range(len(cellSizeLst)):
            cellSize = cellSizeLst[cellSizeInd]
            T = TLst[TInd]
            k0 = k0Lst[k0Ind]
            EA = gamma*Runiv*T

            # -- create caseFolder based on baseCase
            caseName = 'intraTrans_yInf_%g_R_%g_T_%g_cS_%g_k0_%g_tort_%g_inv_%g'%(yInf,R,T,cellSize,k0,tort,inv)
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
                if flow: changeInCaseFolders('system/blockMeshDict',['length1', 'length2', 'width','nDiscX','nDiscYZ'],[str(length1),str(length2),str(width),str(int((length1+length2)/cellSize)),str(int(2*width/cellSize))])
                else: changeInCaseFolders('system/blockMeshDict',['dSN', 'nDisc'],[str(Rinf),str(int(Rinf/cellSize*2))])
                changeInCaseFolders('system/fvSolution',['nTCorr'],[str(numOfTCorr)])
                changeInCaseFolders('constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(k0),str(EA),str(sHr)])
                changeInCaseFolders('constant/transportProperties',['kappaEffSet','tortSet','DSet'],[str(kappaEff),str(tort),str(DFreeZ)])
                changeInCaseFolders('system/snappyHexMeshDict',['spR'],[str(R)])
                # -- run the simulation
                os.chdir(caseDir)
                # os.system('./AllrunIntraSphere') # --> only runs inside the sphere
                # os.system('./Allrun') # --> for flow cases
                os.system('./Allrun-parallel')

            else:
                os.chdir(caseDir)

            # -- load DFree, Deff, and k from log file
            pars = {'DFree':('DFree',-1), 'DEff':('DEff', 0), 'k':('max(k)', -1)}
            # par_vals = pars_from_log(pars, solver)
            # DFree, DEff, k = par_vals
            DFree, DEff, k = 1e-5, 1e-5, 2.06
            DEff = DFree/tort*0.5           

            # -- compute analytical results
            k0Art = k0*np.exp(-EA/(Runiv*T))                            # definition in 2020 Chandra Direct numerical simulation of a noniso...
            rSqIdeal = 4./3*np.pi*R**3*k0Art*yInf*p/Runiv/T             # ideal reaction source
            thiele = R*(k0Art/DEff)**(1/2)                              # thiele modulus
            
            # -- compute simulation results
            # NOTE MK: use custom function
            # -- read real source
            with open('log.intSrcSphere', 'r') as fl:
                lines = fl.readlines()
            # reaction source
            for lineInd in range(len(lines)-1,0,-1):
                if lines[lineInd].find('reaction source') >= 0:
                    rS = float(lines[lineInd].split(' ')[-1].replace('\n',''))
                    etaSim = rS/rSqIdeal  # simulation effectivness factor
                    print('reaction source = %g'%rS)
                    break
                if lineInd == 0:
                    print('Reaction source not found.')
            
            
            
            # -- in the isothermal cases we can compare concentration profiles
            if isothermal and not flow: 
                # -- analytical effectivness factor
                etaAnal = 3./(thiele**2)*(thiele*1./np.tanh(thiele)-1)  
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
                plt.title('Concentration profile comparison for k0=%d'%k0)
                plt.savefig('resHere.png')
                plt.show()
            
            # -- in the nonisothermal case we can compare eta-thiele diagram
            elif not isothermal and not flow:
                beta = yInf*p/Runiv/T*(-sHr)*DEff/kappaEff/T
                print('beta',beta,'thiele',thiele,'effect sim',etaSim)
                resNp[:,k0Ind] = np.array([thiele,etaSim])
                
            elif flow:
            # NOTETH: inner + outer transport
                nu = 5.5862252e-05
                Re = inv * (R*2) / nu
                Sc = nu / DFree
                ShC = 2 + 0.6 * Re**0.5 * Sc**(1./3)
                print('Re = %g, Sc = %g, ShC = %g'%(Re,Sc,ShC))

                with open('log.integrace','r') as fl:
                    lines = fl.readlines()

                inds = [0,0,0,0]
                names = ['areaAverage(batt1) of cCOS',"areaAverage(batt1) of gradCCO",'areaAverage(batt1) of yCOS','areaAverage(batt1) of gradYCO']
                for lineInd in range(len(lines)):
                    for j in range(len(names)):
                        if  names[j] in lines[lineInd]:
                            inds[j] = lineInd
                cCO = float(lines[inds[0]].split('=')[1])
                gradCCO = float(lines[inds[1]].split('=')[1])
                yCO = float(lines[inds[2]].split('=')[1])
                gradYCO = float(lines[inds[3]].split('=')[1])
                j = -gradCCO #/(4*np.pi*R**2)
                km = j/(yInf*p/Runiv/T-cCO)
                Sh = km*(2*R)/DFree

                jY = -gradYCO #/(4*np.pi*R**2)
                kmY = jY/(yInf-yCO)
                ShY = kmY*(2*R)/DFree

                print('correlation: Re = %g, Sc = %g, ShC = %g'%(Re,Sc,ShC))
                print('simulation: thiele = %g, gradCCO = %g, cCO = %g, j = %g, km = %g, Sh = %g\ngradYCO = %g, yCO = %g, jy = %g, kmy = %g, Shy = %g'%(thiele, gradCCO, cCO, j,km,Sh,gradYCO, yCO, jY,kmY,ShY))
            os.chdir('../../')


## -- create plot

# create &or name directory for res plot:
dirName = 'ZZZ_res'
if flow:    dirName += '_flow'
if isothermal: dirName += '_isoT'
if not os.path.exists(dirName): os.mkdir(dirName)

# create plot:
if not flow:
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
        plt.plot(resNp[0,:],resNp[1,:],'x',label='sim res.')
        title = 'Multiple steady states.'
        fileName = 'mult_steadySt.png'
    plt.title(title)
    plt.savefig('%s/%s'%(dirName,fileName))
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.show()
