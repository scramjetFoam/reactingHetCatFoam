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
    print('Reading parameters from log file:')
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
numOfTCorr = 0
# isothermal logic:
isothermal = False
if numOfTCorr == 0: 
    isothermal = True

# -- swicher if the simulation should be run or to just check results
runSim = True  
# runSim = False

# -- case parameters
yInf = 0.2      # molar fraction farfar from sphere
# yInfLst = [0.02]      # molar fraction farfar from sphere
p = 101325      # presure
# sHr = -283e3    # standard reaction enthalpy (physical 283e3)	
# sHrLst = [-12307920,-123079200,-1230792000]     # standard reaction enthalpy (physical 283e3)	
# sHrLst = [-300e3/10,-300e3/5,-300e3]     # standard reaction enthalpy (physical 283e3)	
Runiv = 8.314   # universal gas constant
R = 0.1           # sphere radius
# inv = 0.0558        # inlet velocity
# inv = 0        # inlet velocity
inv = 0.1116        # inlet velocity
# tube dimensions:
length1 = 1.1*R       # inlet <-> sphere centre
length2 = length1       # sphere centre <-> outlet
width = length1         # top|bottom wall <-> sphere centre

DFreeZ = 1e-5

# -- setup study parameters here
baseCaseDir = 'baseCase'
outFolder = 'ZZ_cases'
  
# cellSizeLst = [R/3]  # FV cell Size
# cellSizeLst = [0.5*R]  # FV cell Size
# cellSizeLst = [0.5,0.25,0.125,0.0625,0.03125]  
# cellSizeLst = [0.5,0.25,0.125,0.0625]  # MK
#k0Lst = [1e2, 5e2, 1e3, 5e3, 1e4, 1e5]      # reaction pre-exponential factor
# k0Lst = [1e9]
k0Lst = [1e9,1e9,1e9,1e9,1e9,1e9]
thieleLst = [0.2,0.5,0.75,1.,2,4]
# thieleLst = [1.]
# k0Lst = [1e5]
# EA = 90e3     # reaction activation energy (set according to gamma) 2020 Chandra Direct numerical simulation of a noniso...
TLst = [300]   # temperature (set according to gamma) 2020 Chandra Direct numerical simulation of a noniso...
gamma = 20      # see line above
beta = 0.6
# betaLst = [0.4,0.6,0.8] set according to T) 2020 Chandra Direct numerical simulation of a noniso...
tort = 2      # tortuosity
solverLst = ['reactingHetCatSimpleFoam','scalarTransportFoamCO'] # used solver
solver = solverLst[0]
kappaEff = 2
# kappaEff = 0.024375 

# flow logic:
flow = False
if inv != 0: 
    flow = True  
    R = 0.01
    length1 = 15*R
    length2 = 45*R
    width = 15*R
    yInf = 0.01
    kappaEff = 1
    sHr = -283e3	
    tort = 5
    TLst = [500]
    k0Lst = [1e11]

# cellSizeLst = [1*R,0.7*R,0.5*R,0.4*R,0.3*R]  # FV cell Size
# cellSizeLst = [R]  # FV cell Size
# cellSizeLst = [0.8*R,0.7*R,0.6*R,0.5*R]  # FV cell Size
cellSizeLst = [0.8*R]  # FV cell Size
nCellsSp = 30

# MK: change naming for flow
if flow:
    baseCaseDir += '_flow'
    outFolder += '_flow'
if isothermal:
    outFolder += '_isoT'

# -- numpy array with results
if isothermal:
    resNp = np.zeros((len(cellSizeLst)))
    resNp2 = np.zeros((len(cellSizeLst)))
else:
    resNp = np.zeros((2,len(k0Lst)+1))

# -- create case for:
if flow:
    tortLst = [
       0.5,
       0.5,
       0.5,
       0.5
    #    5,
#        50, 
#        5, 
#        5]
]
    # tortLst = [0.5]
    invLst = [0.11,0.05,0.22,0.44]
#    invLst =  [
    #    0.2234490088,
#        0.2234490088,
#        0.2234490088,
#        0.0279311261,
#        0.1117245044]
for flowInd in range(len(tortLst)):
    for TInd in range(len(TLst)):
        for k0Ind in range(len(k0Lst)):
            for cellSizeInd in range(len(cellSizeLst)):
                # -- set parameters
                cellSize = cellSizeLst[cellSizeInd]
                T = TLst[TInd]
                k0 = k0Lst[k0Ind]
                EA = gamma*Runiv*T
                inv = invLst[flowInd]
                tort = tortLst[flowInd]
                
                DEff = DFreeZ/tort*0.5  

                if not flow:
                    k0 = (thieleLst[k0Ind]/R)**2*DEff/(np.exp(-EA/(Runiv*T)))
                    sHr = -beta/yInf/p*Runiv*T/DEff*kappaEff*T
                

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
                    # if flow: changeInCaseFolders('system/blockMeshDict',['spR','inR','bMinVal','bMaxVal','xMinVal','xMaxVal','nCellsSp'],[str(R),str(R/2.),str(-width/2),str(width/2),str(-length1),str(length2),str(nCellsSp)])
                    if flow: changeInCaseFolders('system/blockMeshDict',['length1', 'length2', 'width','nDiscX','nDiscYZ'],[str(length1),str(length2),str(width),str(int((length1+length2)/cellSize)),str(int(2*width/cellSize))])
                    else: changeInCaseFolders('system/blockMeshDict',['dSN', 'nDisc'],[str(length1),str(int(length1/cellSize*2))])
                    changeInCaseFolders('system/fvSolution',['nTCorr'],[str(numOfTCorr)])
                    changeInCaseFolders('constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(k0),str(EA),str(sHr)])
                    changeInCaseFolders('constant/transportProperties',['kappaEffSet','tortSet','DSet'],[str(kappaEff),str(tort),str(DFreeZ)])
                    changeInCaseFolders('system/snappyHexMeshDict',['spR'],[str(R)])
                    changeInCaseFolders('system/snappyHexMeshDictIntraTrans',['spR'],[str(R)])
                    # -- run the simulation
                    os.chdir(caseDir)
                    if flow:
                        os.system('./Allrun-parallel')
                    else:
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
                if isothermal and not flow: 
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
                elif not isothermal and not flow:
                    etaSim = rS/rSqIdeal  # simulation effectivness factor
                    print('beta',beta,'thiele',thiele,'etaSim',etaSim,'R',R,'k0Art',k0Art,'Deff',DEff)
                    resNp[:,k0Ind] = np.array([thiele,etaSim])
                    
                elif flow:
                # NOTETH: inner + outer transport
                    # -- diffusivity ratio
                    # print("DEff = %g"%DEff)
                    # print("DFree = %g"%DFree)
                    print('DEff/DFree = %g'%(DEff/DFree))

                    # -- Sherwood
                    nu = 5.58622522e-05
                    Re = inv * (R*2) / nu
                    Sc = nu / DFree
                    ShC = 2 + 0.6 * Re**0.5 * Sc**(0.3333333)

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
                    j = +gradCCO #/(4*np.pi*R**2)
                    km = j/(yInf*p/Runiv/T-cCO)
                    Sh = km*(2*R)/DFree

                    jY = +gradYCO #/(4*np.pi*R**2)
                    kmY = jY/(yInf-yCO)
                    ShY = kmY*(2*R)/DFree
                    print('Re = %g'%Re)
                    print('correlation: ShC = %g, Sc = %g'%(ShC,Sc))
                    print('simulation:  Sh =  %g'%(Sh))
                    print('gradCCO = %g, cCO = %g, j = %g, km = %g, Sh = %g\ngradYCO = %g, yCO = %g, jy = %g, kmy = %g, Shy = %g'%(gradCCO, cCO, j,km,Sh,gradYCO, yCO, jY,kmY,ShY))

                    # -- effectivness
                    # cpmpute etaAnal
                    BiM = km*R/DEff
                    etaAnal = 3/(thiele**2) * (thiele/np.tanh(thiele)-1)/(1+(thiele/np.tanh(thiele)-1)/BiM)
                    
                    # read etaSim
                    with open('log.intSrcSphere', 'r') as fl:
                        lines = fl.readlines()
                    for lineInd in range(len(lines)-1,0,-1):
                        if lines[lineInd].find('reaction source') >= 0:
                            rS = float(lines[lineInd].split(' ')[-1].replace('\n',''))
                            etaSim = rS/rSqIdeal
                            # print('reaction source = %g'%rS)
                            print('etaSim =  %g\netaAnal = %g'%(etaSim, etaAnal))
                            break
                        if lineInd == 0:
                            print('Reaction source not found.')
                    
                    # write data for plotting
                    if not os.path.exists('../../ZZZ_res'): os.mkdir('../../ZZZ_res')
                    mode = 'r'
                    if not os.path.isfile('../../ZZZ_res/flow.csv'): mode = 'w'
                    with open('../../ZZZ_res/flow.csv', mode) as f:
                        if mode == 'w': pass
                        else:
                            lines = f.readlines()
                    wr = True
                    for lineInd in range(len(lines)):
                        test = ['%g'%tort, '%g'%Re]
                        if test==lines[lineInd].split(',')[0:2]: 
                            wr=False
                            break
                    if wr: 
                        with open('../../ZZZ_res/flow.csv','a') as f: 
                            f.writelines('%g,%g,%g,%g,%g,%g\n'%(tort,Re,Sh,ShC,etaSim,etaAnal))


                os.chdir('../../')

# flow = False
# isothermal=False

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
    plt.savefig('%s/%s'%(dirName,fileName))
    plt.show()
elif flow and not runSim:
    # Sh = f(Re)
    # eta = f(Re)
    # for various torts

    # 1. open cases run for flow and get the data
    # --> create output files: 
    #     
    #     ShC + Sh + Re + tort
    # --> 

    pass
