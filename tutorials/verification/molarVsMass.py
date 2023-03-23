import matplotlib.pyplot as plt
import numpy as np

# etaODE
with open('innerSpherePartTransp/etaAnal_beta_0.6_gamma_20.csv', 'r') as f1:
    etaODELst = np.array([[float(data.split(',')[0]), float(data.split(',')[1][:-1])] for data in f1.readlines()[1:]])
    
# etaMolar
with open('innerSpherePartTransp/ZZZ_res/etaCsv/simRes0.6_cS_0.4.csv', 'r') as f2:
    etaMolar1Lst = np.array([[float(data.split(',')[0]), float(data.split(',')[1][:-1])] for data in f2.readlines()[1:]])  
with open('innerSpherePartTransp/ZZZ_res/etaCsv/simRes0.6.csv', 'r') as f3:
    etaMolar2Lst = np.array([[float(data.split(',')[0]), float(data.split(',')[1][:-1])] for data in f3.readlines()[1:]])  

# etaMass 1
with open('innerSpherePartTranspMass/ZZZ_res/etaCsv/simRes0.6_cS_0.4.csv', 'r') as f4:
    etaMass1Lst = np.array([[float(data.split(',')[0]), float(data.split(',')[1][:-1])] for data in f4.readlines()[1:]])
with open('innerSpherePartTranspMass/ZZZ_res/etaCsv/simRes0.6.csv', 'r') as f5:
    etaMass2Lst = np.array([[float(data.split(',')[0]), float(data.split(',')[1][:-1])] for data in f5.readlines()[1:]])



etaODEShortLst = etaODELst[0:10000,:]
f = 5
etaODEPltLst = np.array([etaODEShortLst[f*i,:] for i in range(int(len(etaODEShortLst)/f))])

plt.plot(etaODEPltLst[:,0], etaODEPltLst[:,1], linestyle='-', label='ODE1')
plt.plot(etaMolar1Lst[:,0], etaMolar1Lst[:,1], linestyle='none', marker='o', label='molar1')
plt.plot(etaMolar2Lst[:,0], etaMolar2Lst[:,1], linestyle='none', marker='s', label='molar2')
plt.plot(etaMass1Lst[:,0], etaMass1Lst[:,1], linestyle='none', marker='+', label='mass1')
plt.plot(etaMass2Lst[:,0], etaMass2Lst[:,1], linestyle='none', marker='x', label='mass2')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()