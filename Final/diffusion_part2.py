# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 11:06:48 2022

@author: zacht
"""

import numpy as np
import matplotlib.pyplot as plt
#define input parameters
# sigma_a = float(input('enter absorption cross section for the material   '))
sigma_t = 1
# Q = float(input('enter your value of source for the material   '))
# M = int(input('enter length of the slab  ')) # length of slab
# n = int(input('enter the number of cells  ')) # number of cells

Z = [5,10,15,20,25,30,35,40,45,50]
# flux = np.zeros((len(Z),2*n))
# x_matrix = np.zeros((len(Z),2*n))
for m in Z:
    #define input parameters
    n = 50*int(m/5)
    sigma_a = 5.0/m**2
    Q = 10/m**2
    M = 2*m
    D = 1/(3*sigma_t)    # Diffussion coefficient for material
    l_t = m +2*D        # extrapolated width
    h = l_t/n            # mesh spacing
    c = (2*D/h**2) + sigma_a
    f = D/ h**2
    A = np.zeros([n-1, n-1])
    B = np.zeros(n-1)
    x = np.linspace(-(l_t), (l_t), 2*n)
    for i in range(1,n-2):
        # print(i)
        A[i,i-1] = -f
        A[i,i] = c
        A[i,i+1] = -f
        B[i] = Q
    
    # building matrices A and B
    A[0,0] = c
    A[0,1] = -2*f
    A[n-2,n-3] = -f
    A[n-2,n-2] = c
    B[0] = Q 
    B[n-2] = Q 
    # print (A)
    # print(B)
    
    
    guess = np.ones(n-1)
    
    # Let's define our convergence criteria
    tol = 10**(-8)
    #   Initial guess
    mat_x_estimate = guess
    
    # Error
    difference = 1.0
    #   number of iterations
    k = 0
    
    while (difference > tol):
        k = k + 1 # update Number of iteration
        mat_x_ini = 1.0 * mat_x_estimate
        if (k % 1000) == 0:
            print('iterations (M = ' + str(m) + '): ' + str(k))
        for i in range (len(A)):
             mat_x_estimate[i] = (B[i]- np.dot(A[i], mat_x_estimate)+ A[i,i]*mat_x_estimate[i])/A[i,i]
        # Determine 'Error'
        difference = np.linalg.norm(mat_x_estimate - mat_x_ini,ord=2)/np.linalg.norm(mat_x_estimate,ord=2)
        
    print (mat_x_estimate)
    solution = np.zeros(n)
    solution[0:n-1] = mat_x_estimate
    print(solution)
    
    mat_y_estimate = np.zeros(2*(n))
    for i in range(0,n):
        mat_y_estimate[i]= solution[n-1-i]
        mat_y_estimate[n+i]=solution[i]
    # flux[int((m/5)-1),:] = mat_y_estimate
    # x_matrix[int((m/5)-1),:] = x
    np.savetxt("Part2_flux_from_diffusion_norm2_convergence_M"+str(m)+".csv",mat_y_estimate,delimiter=",")
    np.savetxt("x_from_diffusin_M"+str(m)+".csv",x,delimiter=",")

# print(mat_y_estimate)
# print(x)
plt.figure(1)
plt.plot(x,mat_y_estimate)
plt.show()
# np.savetxt("Part2_flux_from_diffusion.csv",flux,delimiter=",")
# np.savetxt("x_from_diffusion.csv",x_matrix,delimiter=",")