# -- Script to do verification simulation of isothermal reaction inside sphere

# -- TODO:
# 1) try different parameters setups, check analytical profile in the isothermal cases, check order of the method (should be 2) in all cases

# -- imports
import numpy as np
import os 
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
    print('Changing file %s, %s --> %s'%(file,str(whatLst),str(forWhatLst)))
    with open(caseDir + '/%s'%file, 'r') as fl:
        data = fl.readlines()

    for whatInd in range(len(whatLst)):
        for ind in range(len(data)):
            data[ind] = data[ind].replace(whatLst[whatInd],forWhatLst[whatInd])

    with open(caseDir + '/%s'%file, 'w') as fl:
        fl.writelines(data)

# -- setup study parameters here
baseCaseDir = 'baseCase'
outFolder = 'ZZ_cases'

# -- number of the enthalpy corrections
# NOTETH: if 0 - isothermal study
numOfTCorr = 1  

# -- swicher if the simulation should be run or to just check results
runSim = False
# runSim = True

# -- case data
yInf = 0.1     # molar fraction farfar from sphere
p = 101325      # presure
sHr = -138725599.72220203    # standard reaction enthalpy (physical 283e3)	
Runiv = 8.314   # universal gas constant
R = 1           # sphere radius
Rinf = 1.1         # infinite radius 
# cellSizeLst = [0.5,0.25,0.125,0.0625,0.03125]  # cell Size
# cellSizeLst = [0.5,0.25,0.125,0.0625]  # cell Size
cellSizeLst = [0.1]  # FV cell Size
k0Lst = [1e2,5e2,1e3,5e3,1e4,1e5]      # reaction pre-exponential factor
# EA = 90e3     # reaction activation energy (set according to gamma) 2020 Chandra Direct numerical simulation of a noniso...
TLst = [500]   # temperature (set according to gamma) 2020 Chandra Direct numerical simulation of a noniso...
gamma = 20      # see line above
# betaLst = [0.4,0.6,0.8] set according to T) 2020 Chandra Direct numerical simulation of a noniso...
tort = 1      # tortuosity
# solverLst = ['reactingHetCatSimpleFoam','scalarTransportFoamCO'] # used solver
solver = 'reactingHetCatSimpleFoam' # used solver
kappaEff = 1

# -- numpy array with results
if numOfTCorr == 0:
    resNp = np.zeros((len(cellSizeLst)))
    resNp2 = np.zeros((len(cellSizeLst)))
else:
    resNp = np.zeros((2,len(k0Lst)+1))

# -- create case
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
                if os.path.isdir(caseDir):                                          #ensure, that the caseDir is clear
                    sh.rmtree(caseDir)

                sh.copytree(baseCaseDir,caseDir)

                # -- change the case files
                changeInCaseFolders('0.org/T',['isoT'],[str(T)])
                changeInCaseFolders('0.org/CO',['yCOSet'],[str(yInf)])
                changeInCaseFolders('system/controlDict',['customSolver'],[solver])
                changeInCaseFolders('system/blockMeshDict',['dSN','nDisc'],[str(Rinf),str(int(Rinf/cellSize)*2)])
                changeInCaseFolders('system/fvSolution',['nTCorr'],[str(numOfTCorr)])
                changeInCaseFolders('constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(k0),str(EA),str(sHr)])
                changeInCaseFolders('constant/transportProperties',['kappaEffSet','tortSet'],[str(kappaEff),str(tort)])
                
                # -- run the simulation
                os.chdir(caseDir)
                os.system('./AllrunIntraSphere')

            else:
                os.chdir(caseDir)

            # -- load Dfree, Deff, and k from log file
            # NOTETH: this can be slow (in that case maybe use something like tail?)
            # This can be nicely done in custom function
            print('Reading parameters from log file.')
            with open('log.%s'%solver, 'r') as fl:
                lines = fl.readlines()
            
            #DFree
            for lineInd in range(len(lines)-1,0,-1):    # -- go from back to get latest result
                if lines[lineInd].find('DFree') >= 0:
                    DFree = float(lines[lineInd].split(' ')[-1].replace('\n',''))
                    print('DFree = %g'%DFree)
                    break
                if lineInd == 0:
                    print('DFree not found.')
            #DEff
            for lineInd in range(len(lines)-1,0,-1):
                if lines[lineInd].find('DEff') >= 0:
                    DEff = float(lines[lineInd].split(',')[0].split(' ')[-1].replace('\n',''))
                    print('DEff = %g'%DEff)
                    break
                if lineInd == 0:
                    print('DEff not found.')
            #k
            for lineInd in range(len(lines)-1,0,-1):
                if lines[lineInd].find('max(k)') >= 0:
                    k = float(lines[lineInd].split(' ')[-1].replace(').\n',''))
                    print('k = %g'%k)
                    break
                if lineInd == 0:
                    print('k not found.')

            # -- compute analytical results
            k0Art = k0*np.exp(-EA/(Runiv*T))        # definition in 2020 Chandra Direct numerical simulation of a noniso...
            rSqIdeal = 4./3*np.pi*R**3*k0Art*yInf*p/Runiv/T
            thiele = R*(k0Art/DEff)**(0.5)  # thiele modulus
            
            # -- compute simulation results
            # -- read real source
            with open('log.intSrcSphere', 'r') as fl:
                lines = fl.readlines()
            #reaction source
            for lineInd in range(len(lines)-1,0,-1):
                if lines[lineInd].find('reaction source') >= 0:
                    rS = float(lines[lineInd].split(' ')[-1].replace('\n',''))
                    break
                if lineInd == 0:
                    print('Reaction source not found.')
            print('reaction source = %g'%rS)
            etaSim = rS/rSqIdeal
            
            # -- in the isothermal cases we can compare profiles
            if numOfTCorr == 0: 
                etaAnal = 3./(thiele**2) * (thiele*1./np.tanh(thiele) - 1)                      # analytical effectivness factor
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
if numOfTCorr == 0:
    plt.plot(cellSizeLst,resNp,label='eta diff (my)')
    # plt.plot(cellSizeLst,resNp[1],label='eta diff (STF)')
    plt.plot(cellSizeLst,resNp2,label='whole sol diff (my)')
    # plt.plot(cellSizeLst,resNp2[1],label='whole sol diff (STF)')
    plt.plot(cellSizeLst,cellSizeLst,label='slope = 1')
    plt.plot(cellSizeLst,np.array(cellSizeLst)**2,label='slope = 2')
    plt.title('Dependence of the error on the mesh.')
    plt.savefig('ZZZ_res/error_mesh.png')
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.show()
else:
    plt.plot(resNp[0,:],resNp[1,:],'x',label='sim res.')
    plt.title('Multiple steady states.')
    plt.savefig('ZZZ_res/mult_steadySt.png')
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.show()