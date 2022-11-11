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
yInf = 0.01      # molar fraction farfar from sphere
p = 101325      # presure
Runiv = 8.314   # universal gas constant
R = 0.01           # sphere radius
invLst = [0.11,0.11,0.22,0.22]        # inlet velocity
sHr = -283e3    # standard reaction enthalpy (physical 283e3)	
tortLst = [0.5, 5, 0.5, 5]
TLst = [500]
k0Lst = [1e9]
EA = 90e3
# tube dimensions:
length1 = 15*R       # inlet <-> sphere centre
length2 = 45*R       # sphere centre <-> outlet
width = 15*R         # top|bottom wall <-> sphere centre

DFreeZ = 1e-5

# -- setup study parameters here
baseCaseDir = 'baseCase'
outFolder = 'ZZ_cases'
  
solverLst = ['reactingHetCatSimpleFoam','scalarTransportFoamCO'] # used solver
solver = solverLst[0]
kappaEff = 1

flow = True  

cellSizeLst = [0.7*R]  # FV cell Size

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

for flowInd in range(len(tortLst)):
    for TInd in range(len(TLst)):
        for k0Ind in range(len(k0Lst)):
            for cellSizeInd in range(len(cellSizeLst)):
                # -- set parameters
                cellSize = cellSizeLst[cellSizeInd]
                T = TLst[TInd]
                k0 = k0Lst[k0Ind]
                inv = invLst[flowInd]
                tort = tortLst[flowInd]
                
                DEff = DFreeZ/tort*0.5

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
                    else: changeInCaseFolders('system/blockMeshDict',['dSN', 'nDisc'],[str(length1),str(int(length1/cellSize*2))])
                    changeInCaseFolders('system/fvSolution',['nTCorr'],[str(numOfTCorr)])
                    changeInCaseFolders('constant/reactiveProperties',['k0Set','EASet','sHrSet'],[str(k0),str(EA),str(sHr)])
                    changeInCaseFolders('constant/transportProperties',['kappaEffSet','tortSet','DSet'],[str(kappaEff),str(tort),str(DFreeZ)])
                    changeInCaseFolders('system/snappyHexMeshDict',['spR'],[str(R)])
                    changeInCaseFolders('system/snappyHexMeshDictIntraTrans',['spR'],[str(R)])
                    
                    # -- run the simulation
                    os.chdir(caseDir)
                    os.system('./Allrun-parallel')

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
                kmC = ShC*DFree/(2*R)
                BiM = kmC*R/DEff
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
