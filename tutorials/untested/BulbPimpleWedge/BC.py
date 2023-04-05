import numpy as np


# H2 CO N2
yIn1 = np.array([0,0.49914])
yIn1 = np.r_[yIn1, 1 - np.sum(yIn1)]
Mg1 = np.array([2.04e-3, 44.01e-3, 28.0e-3])
MgIn1 = np.sum( yIn1 * Mg1 )
wIn1 = yIn1 * Mg1 / MgIn1
print(wIn1)

# H2 CO N2
yIn1 = np.array([0.50121,0])
yIn1 = np.r_[yIn1, 1 - np.sum(yIn1)]
Mg1 = np.array([2.04e-3, 44.01e-3, 28.0e-3])
MgIn1 = np.sum( yIn1 * Mg1 )
wIn1 = yIn1 * Mg1 / MgIn1
print(wIn1)