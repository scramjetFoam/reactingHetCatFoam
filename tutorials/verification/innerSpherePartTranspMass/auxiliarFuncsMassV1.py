import numpy as np
import os
import re
import shutil as sh
import matplotlib.pyplot as plt
import sys
try:
    from scipy.stats import linregress
except ModuleNotFoundError:
    print('Warning: missing scipy')
    noscipy = True
else:
    noscipy = False


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
    with open('log.intSrcSphereM', 'r') as fl:
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


def log_report(thiele,DEff,DFree,Re,ShC,Sc,Sh,eta_sim,eta_corr,eta_anal,gradCCO,cCO,j,km,gradYCO,yCO,jY,kmY,ShY):
    print('==========================================================')
    print('Case with thiele = %g'%thiele)
    print('DEff/DFree = %g'%(DEff/DFree))
    print('Re = %g'%Re)
    print('correlation: ShC = %g, Sc = %g'%(ShC,Sc))
    print('simulation:  Sh =  %g'%Sh)
    print('eta_sim =  %g\neta_corr = %g\neta_anal = %g'%(eta_sim,eta_corr,eta_anal))
    print('gradCCO = %g, cCO = %g, j = %g, km = %g, Sh = %g'%(gradCCO,cCO,j,km,Sh))
    print('gradYCO = %g, yCO = %g, jy = %g, kmy = %g, Shy = %g'%(gradYCO,yCO,jY,kmY,ShY))
    print('==========================================================')


# NOTE: NO LONGER IN USE
# def Sh_plt(ZZZ_filepath, Re1, Re2, nu, DFree, tort):
#     """Create Sh-Re plot for a given value of tort"""
#     # -- Frossling correlation:
#     Re_lst = np.linspace(Re1, Re2, 50)
#     ShC_lst = frossling(Re_lst, nu, DFree)
#     fig, ax = plt.subplots()
#     ax.plot(Re_lst, ShC_lst, label='ShC')
#     # -- Simulation results:
#     with open(ZZZ_filepath) as f:
#         # tort, Re, Sh, ShC, eta_sim, eta_anal
#         lines = f.readlines()
#     flowRes = np.zeros((len(lines),len(lines[0].split(','))))
#     for lineInd in range(len(lines)):
#         flowRes[lineInd] = lines[lineInd].split(',')
#     # filter by tort:
#     filtered = flowRes[flowRes[:,0]==tort]
#     plt.plot(filtered[:,1],filtered[:,2], 'x', label='Sh')

#     plt.title('Sherwood Number for tort = %g'%tort)
#     plt.xlabel('Re')
#     plt.ylabel('Sh')
#     plt.legend()
#     plt.show()
 

def flow_csv(ZZZ_path,ZZZ_filepath,tort,Re,eta_sim,eta_corr,eta_anal):
    """Write data from 'flowAroundSphere' for plotting  to a csv file"""
    if not os.path.exists(ZZZ_path): os.mkdir(ZZZ_path)
    # NOTE: if the redundancy test below doesn't work, just delete it
    test_tortRe = [round(tort,3),round(Re,3)]
    if os.path.isfile(ZZZ_filepath):
        with open(ZZZ_filepath,'r') as f1:
            lines = f1.readlines()
        for line in lines:
            line_tortRe = [round(float(line.split(',')[i]),3) for i in range(2)]
            if test_tortRe == line_tortRe: 
                return -1
    
    with open(ZZZ_filepath,'a') as f2: 
        f2.writelines('%g,%g,%g,%g,%g\n'%(tort,Re,eta_sim,eta_corr,eta_anal))


def generate_eta_csvs(ZZZ_path,ZZZ_filepath,tortLst):
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
    for tort in tortLst:
        # filter
        filtered = flowRes[flowRes[:,0]==tort]
        Re_lst = filtered[:,1]
        eta_sim_lst = filtered[:,2]
        eta_corr_lst = filtered[:,3]
        # write to the file
        with open('%s/etaCorr_tort%g.csv'%(ZZZ_path_csv,tort), 'w') as f2:
            f2.writelines(['x,y\n'])
            f2.writelines(['%g, %g\n'%(Re_lst[i],eta_sim_lst[i]) for i in range(len(Re_lst))])
        with open('%s/etaSim_tort%g.csv'%(ZZZ_path_csv,tort), 'w') as f3:
            f3.writelines(['x,y\n'])
            f3.writelines(['%g, %g\n'%(Re_lst[i],eta_corr_lst[i]) for i in range(len(Re_lst))])
        with open('%s/etaErr_tort%g.csv'%(ZZZ_path_csv,tort), 'w') as f4:
            f4.writelines(['x,y\n'])
            f4.writelines(['%g, %g\n'%(Re_lst[i],abs(eta_sim_lst[i]-eta_corr_lst[i])/eta_corr_lst[i]) for i in range(len(Re_lst))])


def eta_plt(ZZZ_filepath, tort):
    """Create Re-eta plots for a given value of tort"""
    # -- Simulation results:
    with open(ZZZ_filepath) as f:
        # tort, Re, Sh, ShC, eta_sim, eta_anal
        lines = f.readlines()
    flowRes = np.zeros((len(lines),len(lines[0].split(','))))
    for lineInd in range(len(lines)):
        flowRes[lineInd] = lines[lineInd].split(',')
    # filter by tort:
    filtered = flowRes[flowRes[:,0]==tort]
    # tort,Re,eta_sim,eta_corr,eta_anal
    plt.plot(filtered[:,1],filtered[:,2], 'x', label='eta_sim')
    plt.plot(filtered[:,1],filtered[:,3], 'x', label='eta_corr')
    plt.plot(filtered[:,1],filtered[:,4], 'x', label='eta_anal')
    plt.title('Effectivness Comparison for tort = %g'%tort)
    plt.xlabel('Re')
    plt.ylabel('eta')
    plt.legend()
    plt.show()


def eta_err_plt(ZZZ_filepath, tortLst, invLst):
    """Create Re-err(eta) plots for a given value of tort"""
    with open(ZZZ_filepath) as f:
        lines = f.readlines()
    flowRes = np.zeros((len(lines),len(lines[0].split(','))))
    for lineInd in range(len(lines)):
        flowRes[lineInd] = lines[lineInd].split(',')
    filtered = [0 for i in tortLst] # filter by tort
    for tort in tortLst:
        i = tortLst.index(tort)
        filtered[i] = flowRes[flowRes[:,0]==tort]
        plt.plot(filtered[i][:,1],abs(filtered[i][:,2]-filtered[i][:,3])/filtered[i][:,3], 'x', label='eta_err for tort %g'%tort)
    plt.title('Error comparrison for tort from %s'%str(tortLst))
    plt.xlabel('Re')
    plt.ylabel('eta')
    plt.legend()
    plt.show()

def mesh_err_csv(ZZZ_path, id_parameters, cell_size, err):
    # -- for given parameters, write meshSize x error --> must be in the main loop
    # -- make sure dir exists
    ZZZ_path = "../../" + ZZZ_path
    ZZZ_filepath = ZZZ_path + "/meshErr" + ".csv"
    if not os.path.isdir(ZZZ_path):
        print("creating directory") 
        os.mkdir(ZZZ_path)
    # -- make sure file with header exists
    if not os.path.isfile(ZZZ_filepath):
        print("creating file")
        with open(ZZZ_filepath, 'w') as f1: 
            f1.writelines(['cellSize, err\n'])
    # -- write data
    print("writing to file")
    with open(ZZZ_filepath, 'a') as f2:
        f2.writelines(['%g,%g\n'%(cell_size, err)])

# == regression for error mesh dependence
if not noscipy:
    def logfit(emdNp):
        lr = linregress(np.log(emdNp))
        l = lambda x : (x**lr.slope)*np.exp(lr.intercept)
        return l(emdNp[0]), lr.slope