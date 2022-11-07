import numpy as np

# -- H2, CO, acetone
x = np.array([0.25,0.25,0.5])
Dii = np.array([[0, 9.477e-5, 5.082e-5],
                [9.477e-5,0, 1.26758420701e-05],
                [5.082e-5,1.26758420701e-05,0]])

# -- acetone, CO, H2
Dii = np.array([[0, 1.26758420701e-05, 5.082e-5],
                [1.26758420701e-05,0, 9.477e-5],
                [5.082e-5,9.477e-5,0]])

# # -- methanol, acetone, N2
# Dii = np.array([[0, 8.48, 13.72],
#                 [8.48,0, 19.91],
#                 [13.72,19.91,0]])

B = np.zeros((len(x)-1,len(x)-1))
n = len(x)
for i in range(len(x)-1):
    for k in range(len(x)-1):
        if i == k:
            sum = 0
            for j in range(n):
                if j != i:
                    sum += x[j] / Dii[i,j]
            B[i,k] = x[i]/Dii[i,n-1] + sum

        else:
            B[i,k] = -x[i]*(1./Dii[i,k]-1./Dii[i,n-1])
print(B)
