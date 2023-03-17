import matplotlib.pyplot as plt
import numpy as np

# etaODE
with open('innerSpherePartTransp/etaAnal_beta_0.6_gamma_20.csv', 'r') as f1:
    etaODELst = np.array([[float(data.split(',')[0]), float(data.split(',')[1][:-1])] for data in f1.readlines()[1:]])
    
# etaMolar
with open('innerSpherePartTransp/ZZZ_res/etaCsv/simRes0.6.csv', 'r') as f2:
    etaMolarLst = np.array([[float(data.split(',')[0]), float(data.split(',')[1][:-1])] for data in f2.readlines()])  

# etaMass
with open('innerSpherePartTranspMass/ZZZ_res/etaCsv/simRes0.6.csv', 'r') as f3:
    etaMassLst = np.array([[float(data.split(',')[0]), float(data.split(',')[1][:-1])] for data in f3.readlines()])

etaODEShortLst = etaODELst[0:-215000,:]
f = 5
etaODEPltLst = np.array([etaODEShortLst[f*i,:] for i in range(int(len(etaODEShortLst)/f))])

plt.plot(etaODEPltLst[:,0], etaODEPltLst[:,1],label='ODE')
plt.plot(etaMolarLst[:,0], etaMolarLst[:,1], linestyle='none', marker='o',label='molar')
plt.plot(etaMassLst[:,0], etaMassLst[:,1], linestyle='none', marker='x',label='mass')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()