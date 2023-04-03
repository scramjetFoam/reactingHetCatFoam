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
yIn = np.array([0.319,0.528])
yIn = np.r_[yIn, 1 - np.sum(yIn)]
Mg = np.array([58.08e-3, 32.04e-3, 28.0e-3])
MgIn = np.sum( yIn * Mg )
wIn = yIn * Mg / MgIn

# parameters
N3 = 0

cT = 36.39291352

D12 = 8.48e-6
D13 = 13.72e-6
D23 = 19.91e-6

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

    sol = solve_ivp(returns_dydt, [0.0,0.238],  [0.319,0.528], args=(N1,N2,),method='DOP853',max_step=1e-3)

    x1_sol = sol.y[0,-1]
    x2_sol = sol.y[1,-1]

    # print(sol)
    # return [x1_sol.subs(z,0.238),x2_sol.subs(z,0.238)]
    return [x1_sol,x2_sol]

sol = root(newton, np.array([1.7e-3,3.1e-3]))
print(sol)

N1, N2 = sol.x

# z = symbols('z')
# x1 = symbols('x1', cls=Function)
# x2 = symbols('x2', cls=Function)

# eq1 = Eq(x1(z).diff(z) - ((x1(z)*N2-x2(z)*N1)/D12+(x1(z)*N3-(1-x1(z)-x2(z))*N1)/D13)/cT, 0)
# eq2 = Eq(x2(z).diff(z) - ((x2(z)*N1-x1(z)*N2)/D12+(x2(z)*N3-(1-x1(z)-x2(z))*N2)/D23)/cT, 0)

# sol = dsolve((eq1, eq2),[x1(z),x2(z)],ics={x1(0):0.319, x2(0):0.528})

z_vals = np.linspace(0, 0.238, 100)

# x1_sol = sol[0].rhs
# x2_sol = sol[1].rhs

# x_vals = [x1_sol.subs(z, t_val) for t_val in t_vals]
# y_vals = [x2_sol.subs(z, t_val) for t_val in t_vals]
# z_vals = [1-x1_sol.subs(z, t_val)-x2_sol.subs(z, t_val) for t_val in t_vals]

sol = solve_ivp(returns_dydt, [0.0,0.238], [0.319,0.528],args=(N1,N2,),method='DOP853',max_step=1e-3)
# print(sol.y)

# print(x_vals)

figure, axis = plt.subplots(1, 2, figsize=(30, 15))
axis[0].plot(sol.t, sol.y[0,:], 'b', label='yAc(z)')
axis[0].plot(sol.t, sol.y[1,:], 'g', label='yMeth(z)')
axis[0].plot(sol.t, 1-sol.y[0,:]-sol.y[1,:], 'r', label='yN2(z)')
axis[0].set_xlabel('z')


# # openfoam reseni
baseCaseDir = '../tutorials/untested/massStefTubeVSM/'
outFolder = '../ZZ_cases/SMTest'

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


axis[0].plot(dataOF[:,0],dataOF[:,2],"g--",label='yMethOFSM')
axis[0].plot(dataOF[:,0],dataOF[:,1],"b--",label='yAcOFSM')
axis[0].plot(dataOF[:,0],dataOF[:,3],"r--",label='yN2OFSM')
print(wIn)
print(MgIn)

# openfoam reseni
baseCaseDir = '../tutorials/untested/massStefTubeV1/'
outFolder = '../ZZ_cases/FTest'

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

# dataOF = np.genfromtxt(case.dir + '/postProcessing/graphUniform/400/line.xy')
dataOF = np.genfromtxt(case.dir + '/postProcessing/graphUniform/400/line.xy')
# dataOF = np.genfromtxt('../ZZ_cases/data/lineSMMatMult.xy')
# dataOF = np.genfromtxt(case.dir + 'postProcessing/graphUniform/400/line.xy', delimiter='\t')

axis[0].plot(dataOF[:,0],dataOF[:,2],"g:",label='yMethOFF')
axis[0].plot(dataOF[:,0],dataOF[:,1],"b:",label='yAcOFF')
axis[0].plot(dataOF[:,0],dataOF[:,3],"r:",label='yN2OFF')

# experimental data
# data = np.loadtxt("../ZZ_cases/data/myData.dat", delimiter="\t", skiprows=1)
data = np.genfromtxt('../ZZ_cases/data/myData.dat', delimiter='\t')
axis[0].plot(data[:,7],data[:,8],"or",label='yN2Exp')
axis[0].plot(data[:,9],data[:,10],"og",label='yMethExp')
axis[0].plot(data[:,11],data[:,12],"ob",label='yAcExp')

axis[0].set_xlim((0,0.238))
axis[0].set_ylim((0,1))
axis[0].grid()
axis[0].legend(loc='best')

# openfoam reseni
baseCaseDir = '../ZZ_cases/SM1Test'
outFolder = '../ZZ_cases/SMTest2'

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

axis[1].plot(dataOF[:,0],dataOF[:,2],"g--",label='yH2OFSM')
axis[1].plot(dataOF[:,0],dataOF[:,1],"b--",label='yCOOFSM')
axis[1].plot(dataOF[:,0],dataOF[:,3],"r--",label='yN2OFSM')

# openfoam reseni
baseCaseDir = '../ZZ_cases/F1Test/'
outFolder = '../ZZ_cases/F1Test2'

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

# dataOF = np.genfromtxt(case.dir + '/postProcessing/graphUniform/400/line.xy')
dataOF = np.genfromtxt(case.dir + '/postProcessing/graphUniform/400/line.xy')
# dataOF = np.genfromtxt('../ZZ_cases/data/lineSMMatMult.xy')
# dataOF = np.genfromtxt(case.dir + 'postProcessing/graphUniform/400/line.xy', delimiter='\t')

axis[1].plot(dataOF[:,0],dataOF[:,2],"g:",label='yH2OFF')
axis[1].plot(dataOF[:,0],dataOF[:,1],"b:",label='yCOOFF')
axis[1].plot(dataOF[:,0],dataOF[:,3],"r:",label='yN2OFF')

# experimental data
# data = np.loadtxt("../ZZ_cases/data/myData.dat", delimiter="\t", skiprows=1)
# data = np.genfromtxt('../ZZ_cases/data/myData.dat', delimiter='\t')
# axis[1].plot(data[:,7],data[:,8],"or",label='yN2Exp')
# axis[1].plot(data[:,9],data[:,10],"og",label='yMethExp')
# axis[1].plot(data[:,11],data[:,12],"ob",label='yAcExp')

axis[1].set_xlim((0,0.238))
axis[1].set_ylim((0,1))
axis[1].grid()
axis[1].legend(loc='best')



plt.show()