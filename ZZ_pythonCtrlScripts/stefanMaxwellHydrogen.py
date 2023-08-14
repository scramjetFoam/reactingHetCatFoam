# python script to test Stefan Maxwell diffusion in Stefan tube agains analytical solution

import numpy as np  
from OF_caseClass import OpenFOAMCase
from sympy import symbols, Function, Eq, dsolve
import matplotlib.pyplot as plt
from scipy.optimize import root
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

# specify the species
# test case 1: experimental case with methanol, acetone and N2


# parameters
N3 = 0

cT = 36.39291352

# y1I = 0.3173
# y2I = 0.5601
# y1O = 0.1
# y2O = 0.7774

y1I = 0.7
y2I = 0.0
y1O = 0.0
y2O = 0.7

yIn = np.array([y1I,y2I])
yO = np.array([y1O,y2O])
yIn = np.r_[yIn, 1 - np.sum(yIn)]
yO = np.r_[yO, 1 - np.sum(yO)]
# Mg = np.array([58.08e-3, 32.04e-3, 28.0e-3])
Mg = np.array([2e-3, 17e-3, 28.0e-3])
MgIn = np.sum( yIn * Mg )
MgO = np.sum( yO * Mg )
wIn = yIn * Mg / MgIn

wO = yO * Mg / MgO

p = 101308
R = 8.314
T = 303
cT = p/R/T

# D12 = 8.48e-6
# D13 = 13.72e-6
# D23 = 19.91e-6

D12 = 0.806e-4
D13 = 0.805e-4
D23 = 0.210e-4

l = 0.238 

# l = 0.2843


def returns_dydt(z,x,*par):
    N1,N2 = par
    x1,x2 = x
    dx1dt =  ((x1*N2-x2*N1)/D12+(x1*N3-(1-x1-x2)*N1)/D13)/cT
    dx2dt =  ((x2*N1-x1*N2)/D12+(x2*N3-(1-x1-x2)*N2)/D23)/cT
    return [dx1dt,dx2dt]

# numerical solution
def newton(x):
    N1, N2 = x

    # z = symbols('z')
    # x1 = symbols('x1', cls=Function)
    # x2 = symbols('x2', cls=Function)

    # eq1 = Eq(x1(z).diff(z) - ((x1(z)*N2-x2(z)*N1)/D12+(x1(z)*N3-(1-x1(z)-x2(z))*N1)/D13)/cT, 0)
    # eq2 = Eq(x2(z).diff(z) - ((x2(z)*N1-x1(z)*N2)/D12+(x2(z)*N3-(1-x1(z)-x2(z))*N2)/D23)/cT, 0)

    # sol = dsolve((eq1, eq2),[x1(z),x2(z)],ics={x1(0):0.319, x2(0):0.528})

    sol = solve_ivp(returns_dydt, [0.0,l],  [y1I,y2I], args=(N1,N2,),method='DOP853',max_step=1e-3)

    x1_sol = sol.y[0,-1]
    x2_sol = sol.y[1,-1]

    # print(sol)
    # return [x1_sol.subs(z,0.238),x2_sol.subs(z,0.238)]
    return [x1_sol-y1O,x2_sol-y2O]

sol = root(newton, np.array([1.7e-3,3.1e-3]))
print(sol)

N1, N2 = sol.x

# z = symbols('z')
# x1 = symbols('x1', cls=Function)
# x2 = symbols('x2', cls=Function)

# eq1 = Eq(x1(z).diff(z) - ((x1(z)*N2-x2(z)*N1)/D12+(x1(z)*N3-(1-x1(z)-x2(z))*N1)/D13)/cT, 0)
# eq2 = Eq(x2(z).diff(z) - ((x2(z)*N1-x1(z)*N2)/D12+(x2(z)*N3-(1-x1(z)-x2(z))*N2)/D23)/cT, 0)

# sol = dsolve((eq1, eq2),[x1(z),x2(z)],ics={x1(0):0.319, x2(0):0.528})

z_vals = np.linspace(0, l, 100)

# x1_sol = sol[0].rhs
# x2_sol = sol[1].rhs

# x_vals = [x1_sol.subs(z, t_val) for t_val in t_vals]
# y_vals = [x2_sol.subs(z, t_val) for t_val in t_vals]
# z_vals = [1-x1_sol.subs(z, t_val)-x2_sol.subs(z, t_val) for t_val in t_vals]

sol = solve_ivp(returns_dydt, [0.0,l], [y1I,y2I],args=(N1,N2,),method='DOP853',max_step=1e-3)
# print(sol.y)

# print(x_vals)

figure, axis = plt.subplots(1, 2, figsize=(30, 15))
axis[1].plot(sol.t, sol.y[0,:], 'b', label='yAc(z)')
axis[1].plot(sol.t, sol.y[1,:], 'g', label='yMeth(z)')
axis[1].plot(sol.t, 1-sol.y[0,:]-sol.y[1,:], 'r', label='yN2(z)')
axis[1].set_xlabel('z')

# with open('/home/tomas/01_materials/00_myMats/2301_STC/vzor_LaTeX/02_graphs/stefMax/anal.dat','w') as f:
#     f.writelines('z\tac\tmet\tN2\n')
#     for i in range(len(sol.t)):
#         f.writelines('%g\t%g\t%g\t%g\n'%(sol.t[i],sol.y[0,i],sol.y[1,i],1-sol.y[0,i]-sol.y[1,i]))


# # openfoam reseni
baseCaseDir = '../ZZ_cases/SMTestH2_2/'
outFolder = '../ZZ_cases/SMTestH2_2run'

case = OpenFOAMCase()
case.loadOFCaseFromBaseCase(baseCaseDir)
case.changeOFCaseDir(outFolder)
case.copyBaseCase()

case.runCommands(
    [
        # './Allclean',
        # './Allrun',
        'reactingHetCatSimpleFoamMSM',
        'postProcess -func graphUniform',
    ]
)

# dataOF = np.genfromtxt(case.dir + '/postProcessing/graphUniform/800/line.xy')
dataOF = np.genfromtxt(case.dir + '/postProcessing/graphUniform/400/line.xy')
# dataOF = np.genfromtxt(case.dir + 'postProcessing/graphUniform/400/line.xy', delimiter='\t')


axis[1].plot(dataOF[:,0],dataOF[:,2],"g--",label='yMethOFSM')
axis[1].plot(dataOF[:,0],dataOF[:,1],"b--",label='yAcOFSM')
axis[1].plot(dataOF[:,0],dataOF[:,3],"r--",label='yN2OFSM')
print(wIn)
print(wO)
print(MgIn)

# # openfoam reseni
baseCaseDir = '../ZZ_cases/SMTestH2_2/'
outFolder = '../ZZ_cases/SMTestH2_2runF'

case = OpenFOAMCase()
case.loadOFCaseFromBaseCase(baseCaseDir)
case.changeOFCaseDir(outFolder)
case.copyBaseCase()

case.runCommands(
    [
        # './Allclean',
        # './Allrun',
        'reactingHetCatSimpleFoamM',
        'postProcess -func graphUniform',
    ]
)

# dataOF = np.genfromtxt(case.dir + '/postProcessing/graphUniform/800/line.xy')
dataOF = np.genfromtxt(case.dir + '/postProcessing/graphUniform/400/line.xy')
# dataOF = np.genfromtxt(case.dir + 'postProcessing/graphUniform/400/line.xy', delimiter='\t')


axis[1].plot(dataOF[:,0],dataOF[:,2],"g:",label='yMethOFF')
axis[1].plot(dataOF[:,0],dataOF[:,1],"b:",label='yAcOFF')
axis[1].plot(dataOF[:,0],dataOF[:,3],"r:",label='yN2OFF')
print(wIn)
print(wO)
print(MgIn)

# with open('/home/tomas/01_materials/00_myMats/2301_STC/vzor_LaTeX/02_graphs/stefMax/OFSM.dat','w') as f:
#     f.writelines('z\tac\tmet\tN2\n')
#     for i in range(len(dataOF[:,0])):
#         f.writelines('%g\t%g\t%g\t%g\n'%(dataOF[i,0],dataOF[i,1],dataOF[i,2],dataOF[i,3]))

# # openfoam reseni
# baseCaseDir = '../tutorials/untested/massStefTubeV1/'
# outFolder = '../ZZ_cases/FTest'

# case = OpenFOAMCase()
# case.loadOFCaseFromBaseCase(baseCaseDir)
# case.changeOFCaseDir(outFolder)
# case.copyBaseCase()

# case.runCommands(
#     [
#         # './Allclean',
#         # './Allrun',
#         'reactingHetCatSimpleFoamM',
#         'postProcess -func graphUniform',
#     ]
# )

# # dataOF = np.genfromtxt(case.dir + '/postProcessing/graphUniform/400/line.xy')
# dataOF = np.genfromtxt(case.dir + '/postProcessing/graphUniform/400/line.xy')
# # dataOF = np.genfromtxt('../ZZ_cases/data/lineSMMatMult.xy')
# # dataOF = np.genfromtxt(case.dir + 'postProcessing/graphUniform/400/line.xy', delimiter='\t')

# axis[0].plot(dataOF[:,0],dataOF[:,2],"g:",label='yMethOFF')
# axis[0].plot(dataOF[:,0],dataOF[:,1],"b:",label='yAcOFF')
# axis[0].plot(dataOF[:,0],dataOF[:,3],"r:",label='yN2OFF')

# # experimental data
# # data = np.loadtxt("../ZZ_cases/data/myData.dat", delimiter="\t", skiprows=1)
# data = np.genfromtxt('../ZZ_cases/data/myData.dat', delimiter='\t')
# axis[0].plot(data[:,7],data[:,8],"or",label='yN2Exp')
# axis[0].plot(data[:,9],data[:,10],"og",label='yMethExp')
# axis[0].plot(data[:,11],data[:,12],"ob",label='yAcExp')

# axis[0].set_xlim((0,0.238))
# axis[0].set_ylim((0,1))
# axis[0].grid()
# axis[0].legend(loc='best')

# with open('/home/tomas/01_materials/00_myMats/2301_STC/vzor_LaTeX/02_graphs/stefMax/OFF.dat','w') as f:
#     f.writelines('z\tac\tmet\tN2\n')
#     for i in range(len(dataOF[:,0])):
#         f.writelines('%g\t%g\t%g\t%g\n'%(dataOF[i,0],dataOF[i,1],dataOF[i,2],dataOF[i,3]))



# # openfoam reseni
# # baseCaseDir = '../ZZ_cases/SM1Test'
# # outFolder = '../ZZ_cases/SMTest2'

# # case = OpenFOAMCase()
# # case.loadOFCaseFromBaseCase(baseCaseDir)
# # case.changeOFCaseDir(outFolder)
# # case.copyBaseCase()

# # case.runCommands(
# #     [
# #         # './Allclean',
# #         # './Allrun',
# #         'reactingHetCatSimpleFoamMSM > log.reac1',
# #         'postProcess -func graphUniform',
# #     ]
# # )

# # # dataOF = np.genfromtxt(case.dir + '/postProcessing/graphUniform/800/line.xy')
# # dataOF = np.genfromtxt(case.dir + '/postProcessing/graphUniform/400/line.xy')
# # # dataOF = np.genfromtxt(case.dir + 'postProcessing/graphUniform/400/line.xy', delimiter='\t')

# # axis[1].plot(dataOF[:,0],dataOF[:,2],"g--",label='yH2OFSM')
# # axis[1].plot(dataOF[:,0],dataOF[:,1],"b--",label='yCOOFSM')
# # axis[1].plot(dataOF[:,0],dataOF[:,3],"r--",label='yN2OFSM')

# # # openfoam reseni
# # baseCaseDir = '../ZZ_cases/F1Test/'
# # outFolder = '../ZZ_cases/F1Test2'

# # case = OpenFOAMCase()
# # case.loadOFCaseFromBaseCase(baseCaseDir)
# # case.changeOFCaseDir(outFolder)
# # case.copyBaseCase()

# # case.runCommands(
# #     [
# #         # './Allclean',
# #         # './Allrun',
# #         'reactingHetCatSimpleFoamM > log.reac1',
# #         'postProcess -func graphUniform',
# #     ]
# # )

# # # dataOF = np.genfromtxt(case.dir + '/postProcessing/graphUniform/400/line.xy')
# # dataOF = np.genfromtxt(case.dir + '/postProcessing/graphUniform/400/line.xy')
# # # dataOF = np.genfromtxt('../ZZ_cases/data/lineSMMatMult.xy')
# # # dataOF = np.genfromtxt(case.dir + 'postProcessing/graphUniform/400/line.xy', delimiter='\t')

# # axis[1].plot(dataOF[:,0],dataOF[:,2],"g:",label='yH2OFF')
# # axis[1].plot(dataOF[:,0],dataOF[:,1],"b:",label='yCOOFF')
# # axis[1].plot(dataOF[:,0],dataOF[:,3],"r:",label='yN2OFF')

# # # experimental data
# # # data = np.loadtxt("../ZZ_cases/data/myData.dat", delimiter="\t", skiprows=1)
# # # data = np.genfromtxt('../ZZ_cases/data/myData.dat', delimiter='\t')
# # # axis[1].plot(data[:,7],data[:,8],"or",label='yN2Exp')
# # # axis[1].plot(data[:,9],data[:,10],"og",label='yMethExp')
# # # axis[1].plot(data[:,11],data[:,12],"ob",label='yAcExp')

axis[1].set_xlim((0,0.238))
axis[1].set_ylim((0,1))
axis[1].grid()
axis[1].legend(loc='best')



plt.show()