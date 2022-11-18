# -- auxiliarFuncsV3.py
# -- file storing auxiliary functions udes in verOuterTransControlV3.py
# -- NOTE V3 changelog:
#       -- csv files changed to include thiele
#       -- generate_eta_csvs changed accordingly
#       -- removed: eta_anal references, eta_err_plt

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

#======================================================================#
#   NOTE: following functions are used in verOuterTransControl


def frossling(Re, nu, DFree):
    """Compute ShC for givern parameters"""
    Sc = nu/DFree
    ShC = 2 + .6*Re**(1/2)*Sc**(1/3)
    return ShC


def log_report(thiele,DEff,DFree,Re,ShC,Sc,Sh,eta_sim,eta_corr,gradCCO,cCO,j,km,gradYCO,yCO,jY,kmY,ShY):
    print('==========================================================')
    print('Case with thiele = %g'%thiele)
    print('DEff/DFree = %g'%(DEff/DFree))
    print('Re = %g'%Re)
    print('correlation: ShC = %g, Sc = %g'%(ShC,Sc))
    print('simulation:  Sh =  %g'%Sh)
    print('eta_sim =  %g\neta_corr = %g'%(eta_sim,eta_corr))
    print('gradCCO = %g, cCO = %g, j = %g, km = %g, Sh = %g'%(gradCCO,cCO,j,km,Sh))
    print('gradYCO = %g, yCO = %g, jy = %g, kmy = %g, Shy = %g'%(gradYCO,yCO,jY,kmY,ShY))
    print('==========================================================')
 

def flow_csv(ZZZ_path,ZZZ_filepath,thiele,tort,Re,eta_sim,eta_corr):
    """Write data from 'flowAroundSphere' for plotting  to a csv file"""
    if not os.path.exists(ZZZ_path): os.mkdir(ZZZ_path)
    # NOTE MK: if the redundancy test below doesn't work, just delete it
    test_vars = [round(thiele,3),round(tort,3),round(Re,3)]
    if os.path.isfile(ZZZ_filepath):
        with open(ZZZ_filepath,'r') as f1:
            lines = f1.readlines()
        for line in lines:
            line_vars = [round(float(line.split(',')[i]),3) for i in range(3)]
            print('test', test_vars)
            print('line', line_vars)
            if test_vars == line_vars: 
                return -1
    
    with open(ZZZ_filepath,'a') as f2: 
        f2.writelines('%g,%g,%g,%g,%g\n'%(thiele,tort,Re,eta_sim,eta_corr))


def generate_eta_csvs(ZZZ_path, ZZZ_filepath, thieleLst, tortLst):
    """Create a separate csv file for each pair of etas"""
    ZZZ_path_csv = ZZZ_path+'/etacsv'
    if not os.path.isdir(ZZZ_path_csv): os.mkdir(ZZZ_path_csv)
    # read original csv:
    with open(ZZZ_filepath) as f1:
        lines = f1.readlines()
    # covert to a numpy array:
    flowRes = np.zeros((len(lines),len(lines[0].split(','))))
    for lineInd in range(len(lines)):
        flowRes[lineInd] = lines[lineInd].split(',')
    # filter and write:
    for thiele in thieleLst:
        for tort in tortLst:
            # -- filter by thiele & tort
            filtered = flowRes[flowRes[:,0]==thiele]
            filtered = filtered[filtered[:,1]==tort]
            Re_lst = filtered[:,2]
            eta_sim_lst = filtered[:,3]
            eta_corr_lst = filtered[:,4]
            # -- write to the file
            with open('%s/etaSim_phi_%g_tort_%2.1f.csv'%(ZZZ_path_csv,thiele,tort), 'w') as f2:
                f2.writelines(['x,y\n'])
                f2.writelines(['%g,%g\n'%(Re_lst[i],eta_corr_lst[i]) for i in range(len(Re_lst))])
            with open('%s/etaErr_phi_%g_tort_%2.1f.csv'%(ZZZ_path_csv,thiele,tort), 'w') as f3:
                f3.writelines(['x,y\n'])
                f3.writelines(['%g,%g\n'%(Re_lst[i],abs(eta_sim_lst[i]-eta_corr_lst[i])/eta_corr_lst[i]) for i in range(len(Re_lst))])


def eta_plt(ZZZ_filepath, thiele, tort):
    """Create Re-eta plots for a given value of tort and thiele"""
    # -- Simulation results:
    with open(ZZZ_filepath) as f1:
        # tort, Re, Sh, ShC, eta_sim, eta_anal
        lines = f1.readlines()
    flowRes = np.zeros((len(lines),len(lines[0].split(','))))
    for lineInd in range(len(lines)):
        flowRes[lineInd] = lines[lineInd].split(',')
    # filter by tort and thiele:
    print(thiele)
    print(tort)
    print(flowRes)
    filtered = flowRes[flowRes[:,0]==thiele]
    print(filtered)
    print(filtered[:,1])
    filtered = filtered[filtered[:,1]==tort]
    plt.plot(filtered[:,2],filtered[:,3], 'x', label='eta_sim')
    with open('ZZZ_res/etacsv/etaCorr_phi_%g_tort_%2.1f.csv'%(thiele,tort), 'r') as f1:    
        lines = f1.readlines()
    corrRes = [corr.replace('\n','').split(',') for corr in lines[1:]]
    Re_corrLst = [float(corr[0]) for corr in corrRes]
    eta_corrLst = [float(corr[1]) for corr in corrRes]
    plt.plot(Re_corrLst,eta_corrLst, '--', label='eta_corr')
    plt.title('Effectivness Comparison for tort = %g'%tort)
    plt.xlabel('Re')
    plt.ylabel('eta')
    plt.legend()
    plt.show()
