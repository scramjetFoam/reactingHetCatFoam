# -- Script to calculate theoretical numerical solution of the Weisz and Hicks case (1961 Weisz The behavior of por cat...)
# -- problem solved acording to section IV in article

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import root

# # -- function with equation
# def equation(y,x,*pars):
#     thiele,beta,gamma = pars
#     # print(np.array([y[1],thiele**2*y[0]*np.exp(gamma*beta*(1-y[0])/(1+beta*(1-y[0])))]))
#     return np.array([y[1],thiele**2*y[0]*np.exp(gamma*beta*(1-y[0])/(1+beta*(1-y[0])))-2/x*y[1]])

# # -- function to compute difference for shooting method
# def toRoot(y00,*pars):
#     y0 = np.array([float(y00[0]),0])
#     res = odeint(equation,y0,np.linspace(1e-7,1,100),pars)
#     return res[-1,0]-1

# # -- function to solve shooting method
# def shootingMethod(pars):
#     res = root(toRoot,np.array([0.]),pars)
#     resFin = odeint(equation,np.array([float(res.x[0]),0]),np.linspace(1e-7,1,100),pars)
#     return resFin[-1,1]

# # -- function to control solve ODE for each thiele
# def analSol(beta, gamma, thieleLst):
#     etaLst = np.zeros(len(thieleLst))
#     for thieleInd in range(len(thieleLst)):
#         thiele = thieleLst[thieleInd]
#         pars = thiele,beta,gamma
#         etaLst[thieleInd] = 3./thiele**2*shootingMethod(pars)
#     return etaLst

# -- function to solve ODE for each thiele
def analSol(beta, gamma, thieleLst):
    etaLst = np.zeros(len(thieleLst))
    for thieleInd in range(len(thieleLst)):
        thiele = thieleLst[thieleInd]
        
    return etaLst
