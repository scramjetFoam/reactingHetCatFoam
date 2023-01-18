# Script to replicate curves from 1961_Weisz optimized for beta and gamma values

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import time

gamma = 20
beta = 0.6

def model2(x,y):
    y1, y2 = y[0], y[1]
    dy1dx = y2
    if x == 0:
        dy2dx = aphi**2. * y1 * np.exp(gamma*beta*(1.-y1)/(1.+beta*(1.-y1)))
    else:
        dy2dx = -2./x*y2 + aphi**2. * y1 * np.exp(gamma*beta*(1.-y1)/(1.+beta*(1.-y1)))
    return np.array([dy1dx, dy2dx])


def event(x,y):
    y1, y2 = y[0], y[1]
    return y1-1

event.direction = 1
event.terminal = True


aphi = 40
y10Lst0 = np.linspace(1.00e-12, 1.00e-09, 100)
y10Lst1 = np.linspace(1.00e-09, 4.00e-7, 100)
y10Lst2 = np.linspace(4.00e-07, 5.00e-7, 10)
y10Lst3 = np.linspace(5.00e-07, 1.00e-03, 100)
y10Lst4 = np.linspace(5.00e-03, 1.00e-02, 100)
y10Lst5 = np.linspace(1.00e-02, 0.990000, 100)
y10Lst6 = np.linspace(0.990000, 1.0, 10)
y10Lst = np.concatenate((
    y10Lst0, 
    y10Lst1, 
    y10Lst2, 
    y10Lst3, 
    y10Lst4, 
    y10Lst5,
    y10Lst6,
    ))

# -- NOTE MK:   The list below is for solution at phi=0.5 rather than the full curve, 
# --            comment out extrapolation to do that.
# y10Lst = np.linspace(4.88e-7, 4.89e-7, 100)  

start_time1 = time.perf_counter()
print('starting computation', end='')
thiele_eta = np.zeros((2, len(y10Lst)))
for i in range(len(y10Lst)):
    y10 = y10Lst[i]
    x = np.array([0.,10.])
    y0=np.array([y10,0.])

    res = integrate.solve_ivp(model2,x,y0,method='LSODA',events=event)

    for j in range(len(res.t)):
        if res.t[-(j+1)] != 0:
            a = 1/res.t[-(j+1)]
            break

    
    thiele_eta[0,i] = aphi/a
    thiele_eta[1,i] = 3./(a*(thiele_eta[0,i])**2.)*res.y[1,-1]

plt.plot(thiele_eta[0],thiele_eta[1], label='direct computation')

sort = np.argsort(thiele_eta[0])
thiele_eta_sorted = thiele_eta[:,sort]

end_time1 = time.perf_counter()
print('\tfinished in %6.5f s'%(end_time1-start_time1))

# -- extrapolation 

eta_ext = []
phi_ext = []

eta1 = thiele_eta_sorted[1,-2]
eta2 = thiele_eta_sorted[1,-1]
phi1 = thiele_eta_sorted[0,-2]
phi2 = thiele_eta_sorted[0,-1]
del_phi = 1e-3
phi = np.sqrt(phi1*phi2)

start_time2 = time.perf_counter()
print('starting extraplation', end='')

iter = 0
while phi < 1.5e2 and iter < 1e6:
    iter += 1
    eta = np.sqrt(eta1*eta2)
    phi = np.sqrt(phi1*phi2)
    phi_new = phi+del_phi

    dlnetadlnphi = np.log(eta1/eta2)/np.log(phi1/phi2)
    dyprime_1dphi = eta*phi*(dlnetadlnphi+2)/3
    yprime_1 = phi**2*eta/3
    eta_new = 3/(phi_new)**2*(yprime_1 + dyprime_1dphi * del_phi)
    
    phi_ext.append(phi_new)
    eta_ext.append(eta_new)
    eta1, eta2 = eta2, eta_new
    phi1, phi2 = phi2, phi_new

end_time2 = time.perf_counter()
print('\tfinished in %6.5f s'%(end_time2-start_time2))

plt.plot(phi_ext, eta_ext, label='extrapolation')

# -- writing to .csv file

filename = 'etaAnal_beta_%g_gamma_%g.csv'%(beta, gamma)
with open(filename, 'w') as f1:
    f1.writelines('phi,\t\teta\n')
    for i in range(np.shape(thiele_eta)[1]):
        f1.writelines(['%8.6f,\t%8.6f\n'%(thiele_eta[0,i],thiele_eta[1,i])])
    for i in range(len(phi_ext)):
        f1.writelines(['%8.6f,\t%8.6f\n'%(phi_ext[i],eta_ext[i])])

print('data written to\t\t%s'%filename)

# -- plotting

plt.yscale('log')
plt.xscale('log')
plt.xlim((1e-1,1e1))
plt.ylim((1e0,1e2))
plt.xlabel('φ')
plt.ylabel('η')
plt.title('Solution to ODE from 1961_Weisz for β=%g, γ=%g'%(beta, gamma))
plt.legend()
plt.show()