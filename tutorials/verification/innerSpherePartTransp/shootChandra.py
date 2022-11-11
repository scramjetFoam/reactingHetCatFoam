# Script to replicate curves from 2020_Chandra and 1961_Weisz

import numpy as np
from scipy import integrate
from scipy import optimize
import matplotlib.pyplot as plt

gamma = 20
beta = 0.6

# as boundary value problem -- not working
def model(x,y):
    y1, y2 = y[0], y[1]
    dy1dx = y2
    if x == 0:
        dy2dx = phi02 * y1 * np.exp(gamma*beta*(1.-y1)/(1.+beta*(1.-y1)))
    else:
        dy2dx = -2./x*y2 + phi02 * y1 * np.exp(gamma*beta*(1.-y1)/(1.+beta*(1.-y1)))
    return np.array([dy1dx, dy2dx])

def event(x,y):
    y1, y2 = y[0], y[1]
    return y1-1

event.direction = 1
event.terminal = True

def model2(x,y):
    y1, y2 = y[0], y[1]
    dy1dx = y2
    if x == 0:
        dy2dx = aphi**2. * y1 * np.exp(gamma*beta*(1.-y1)/(1.+beta*(1.-y1)))
    else:
        dy2dx = -2./x*y2 + aphi**2. * y1 * np.exp(gamma*beta*(1.-y1)/(1.+beta*(1.-y1)))
    return np.array([dy1dx, dy2dx])

def fUcel(y10):
    y20 = 0
    y0 = np.array([y10,y20])
    res = integrate.solve_ivp(model,x,y0,method='LSODA')
    y1v1 = res.y[0,-1]
    print(y1v1)
    return 1.0 - y1v1 

# x = np.array([0.,1.])
# phi02 = 1.
# y10Init = 0.0001
# res = optimize.root(fUcel,y10Init)

a = 10
phi = 4
aphi = a*phi
y10Lst = np.linspace(1e-9,1.,40)
# y10Lst = np.array([0.1])
eta = np.zeros(len(y10Lst))
thiele = np.zeros(len(y10Lst))

for i in range(len(y10Lst)):
    y10 = y10Lst[i]
    x = np.array([0.,5.])
    y0=np.array([y10,0.])

    res = integrate.solve_ivp(model2,x,y0,method='LSODA',events=event)

    a = 1./res.t[-1]
    thiele[i] = aphi/a
    eta[i] = 3./(a*(thiele[i])**2.)*res.y[1,-1]
    print(a,thiele[i],eta[i])

plt.plot(thiele,eta)
plt.yscale('log')
plt.xscale('log')
plt.xlim((1e-1,1e1))
plt.ylim((1e0,1e2))
plt.show()
print(res)

# plt.plot(res.t,res.y[0,:])
# plt.show()



# res = integrate.solve_ivp(model,x,[res.x[0],0],method='LSODA')
# plt.plot(res.t,res.y[0])
# plt.show()

# xPlot = np.linspace(-5,5,50)
# yPlot = np.zeros(len(xPlot))
# for i in range(len(xPlot)):
#     yPlot[i] = fUcel(xPlot[i])

# plt.plot(xPlot,yPlot)
# plt.show()




