# Code for the Neutrons Slowing Down Transport Project

# In this version, the parameters are:
    # Vacuum boundaries on both sides
    # Homogeneous material
    # Spatial discretations are constant
    # 2 angular discretezation
    # All spatial discretizations equal
    
import numpy as np
from gS import gaussSeidel
import matplotlib.pyplot as plt

Sigma_t = 1
Sigma_s = 0.9
S2_quad = np.polynomial.legendre.leggauss(2)
S2_quad = np.array(S2_quad).tolist()

mu = S2_quad[0]
mu.reverse()
w = S2_quad[1]
w.reverse()

h = .1
Q = 1

I = 10
A = np.zeros((2*(I-1),2*(I-1)))
for i in range(I-1):
    A[i,i] = abs(mu[0])/h + Sigma_t - Sigma_s*w[0]
    A[2*(I-1)-(1+i),2*(I-1)-(1+i)] = abs(mu[1])/h + Sigma_t - Sigma_s*w[1]
for i in range(I-2):
    A[i+1,i] = -abs(mu[0])/h
    A[i,2*(I-1)-(2+i)] = -Sigma_s*w[1]
    A[2*(I-1)-(i+1),2*(I-1)-(2+i)] = -abs(mu[1])/h
    A[I+i-1,I-(3+i)] = -Sigma_s*w[0]
    
print(A)
psi = np.ones(2*(I-1))
b = Q/2*np.ones(2*(I-1))
convergence_check = [1]
num_iterations = 0
A_size = len(A)
while max(convergence_check)>.0001:
        #while num_iterations < 1000:
            num_iterations = num_iterations + 1
            old_psi = []
            convergence_check = []
            for i in range(2*(I-1)):
                old_psi.append(psi[i])
                psi[i] = gaussSeidel(psi,A,b,i,A_size)
                convergence_check.append(abs(old_psi[i]-psi[i]))
print(psi)

print(f'Number of iterations: {num_iterations}. Max diff between iterations: {max(convergence_check)}')
phi = np.zeros(I)
phi[0] = psi[2*(I-1)-1]
phi[I-1] = psi[I-2]
for i in range(I-2):
    phi[i+1] = psi[i] + psi[2*(I-1)-(i+2)]
print(phi)

x = np.zeros(I)
for i in range(I):
    x[i] = i
plt.plot(x,phi)
