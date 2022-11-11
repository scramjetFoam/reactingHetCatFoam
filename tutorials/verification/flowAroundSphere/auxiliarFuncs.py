import numpy as np
import os
import re
import shutil as sh
import matplotlib.pyplot as plt
import sys


def isFloat(val):
    """Determine if val is a float"""
    
    try:
        float(val)
        return True
    except:
        return False


def changeInCaseFolders(caseDir,file,whatLst,forWhatLst):
    """Write given parameters to files copied from baseCase."""

    print('Changing file %s, %s --> %s'%(file,str(whatLst),str(forWhatLst)))
    with open(caseDir + '/%s'%file, 'r') as fl:
        data = fl.readlines()

    for whatInd in range(len(whatLst)):
        for ind in range(len(data)):
            data[ind] = data[ind].replace(whatLst[whatInd],forWhatLst[whatInd])

    with open(caseDir + '/%s'%file, 'w') as fl:
        fl.writelines(data)


def read_real_source():
    """Read real reaction source from integration log."""
    with open('log.intSrcSphere', 'r') as fl:
        lines = fl.readlines()
    for lineInd in range(len(lines)-1,0,-1):
        if lines[lineInd].find('reaction source') >= 0:
            rS = float(lines[lineInd].split(' ')[-1].replace('\n',''))
            return rS
        if lineInd == 0:
            print('Reaction source not found.')
            return -1


def pars_from_log(pars, solver):
    """Read parameters given by a dictionary from a temporary log file."""

    print('Reading parameters from log file:')
    par_vals = []
    par_keys = list(pars.keys())
    
    # create, read & remove temporary log file:
    with open('tmplog.%s'%solver, 'w') as tmplog: tmplog.write(os.popen('tail -n 100 log.%s'%solver).read())
    with open('tmplog.%s'%solver, 'r') as tmplog: lines = tmplog.readlines() 
    for par_key in par_keys:
        for lineInd in range(len(lines)-1,0,-1):
            if lines[lineInd].find(pars[par_key][0]) >= 0:
                num_lst = [t[0] for t in re.findall("-?\d+.?\d*e?-?\d*", lines[lineInd])]
                val = num_lst[pars[par_key][1]]
                par_vals.append(float(val))
                print('%s = %s'%(par_key, val))
                break
            if lineInd == 0:
                par_vals.append("N/A")
                print('%s not found'%pars[par_key][0])
    os.remove('tmplog.%s'%solver)
    return par_vals


def log_report(thiele,DEff,DFree,Re,ShC,Sc,Sh,etaSim,etaAnal,gradCCO,cCO,j,km,gradYCO,yCO,jY,kmY,ShY):
    print('=============================')
    print('Case with thiele = %g'%thiele)
    print('DEff/DFree = %g'%(DEff/DFree))
    print('Re = %g'%Re)
    print('correlation: ShC = %g, Sc = %g'%(ShC,Sc))
    print('simulation:  Sh =  %g'%Sh)
    print('etaSim =  %g\netaAnal = %g'%(etaSim,etaAnal))
    print('gradCCO = %g, cCO = %g, j = %g, km = %g, Sh = %g'%(gradCCO,cCO,j,km,Sh))
    print('gradYCO = %g, yCO = %g, jy = %g, kmy = %g, Shy = %g'%(gradYCO,yCO,jY,kmY,ShY))


def write_to_csv(tort,Re,Sh,ShC,etaSim,etaAnal):
    """Write data for plotting to a csv file"""
    filepath = '../../ZZZ_res_flow/'
    filename = 'flow.csv'
    if not os.path.exists(filepath): os.mkdir(filepath)
    mode = ('r' if os.path.isfile(filepath+filename) else 'w')
    with open(filepath+filename, mode) as f:
        if mode != 'w': lines = f.readlines()
    # only write for new tort-Re combinantion
    wr = True
    for lineInd in range(len(lines)):
        test = ['%g'%tort, '%g'%Re]
        if test==lines[lineInd].split(',')[0:2]: 
            wr=False
            break
    if wr: 
        with open('../../ZZZ_res_flow/flow.csv','a') as f: 
            f.writelines('%g,%g,%g,%g,%g,%g\n'%(tort,Re,Sh,ShC,etaSim,etaAnal))


def frossling(Re, nu, DFree):
    """Compute ShC for givern parameters"""

    Sc = nu/DFree
    ShC = 2 + .6*Re**(1/2)*Sc**(1/3)
    return ShC

def Sh_plt(filepath, Re1, Re2, nu, DFree, tort):
    """Create Sh-Re plot for a given value of tort"""

    # -- Frossling correlation:
    Re_lst = np.linspace(Re1, Re2, 50)
    ShC_lst = frossling(Re_lst, nu, DFree)
    fig, ax = plt.subplots()
    ax.plot(Re_lst, ShC_lst, label='ShC')
    # -- Simulation results:
    with open(filepath) as f:
        # tort, Re, Sh, ShC, etaSim, etaAnal
        lines = f.readlines()
    flowRes = np.zeros((len(lines),len(lines[0].split(','))))
    for lineInd in range(len(lines)):
        flowRes[lineInd] = lines[lineInd].split(',')
    # filter by tort:
    filtered = flowRes[flowRes[:,0]==tort]
    plt.plot(filtered[:,1],filtered[:,2], 'x', label='Sh')

    plt.title('Sherwood Number for tort = %g'%tort)
    plt.xlabel('Re')
    plt.ylabel('Sh')
    plt.legend()
    plt.show()

 
def eta_plt(filepath, tort):
    """Create eta-Re plot for a given value of tort"""
    # -- Simulation results:
    with open(filepath) as f:
        # tort, Re, Sh, ShC, etaSim, etaAnal
        lines = f.readlines()
    flowRes = np.zeros((len(lines),len(lines[0].split(','))))
    for lineInd in range(len(lines)):
        flowRes[lineInd] = lines[lineInd].split(',')
    # filter by tort:
    filtered = flowRes[flowRes[:,0]==tort]
    plt.plot(filtered[:,1],filtered[:,4], 'x', label='eta_sim')
    plt.plot(filtered[:,1],filtered[:,5], 'x', label='eta_corr')
    
    plt.title('Effectivness Comparison for tort = %g'%tort)
    plt.xlabel('Re')
    plt.ylabel('eta')
    plt.legend()
    plt.show()