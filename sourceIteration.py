# Source iteration transport
import numpy as np
import matplotlib.pyplot as plt


S2_quad = np.polynomial.legendre.leggauss(2)
S2_quad = np.array(S2_quad).tolist()
mu = S2_quad[0]
mu.reverse()
w = S2_quad[1]
w.reverse()

# Material one
Sigma_t_m1 = 1
Sigma_s_m1 = 0.1
Q_m1 = .1
width_m1 = 10
I_m1 = 1000
h_m1 = width_m1/I_m1
# Material two
Sigma_t_m2 = 1
Sigma_s_m2 = 0.5
Q_m2 = .1
width_m2 = 10
I_m2 = 100
h_m2 = width_m2/I_m2


psi_guess = np.zeros((I_m1+I_m2,2))
flux_guess = np.sum(psi_guess,axis=1).tolist()

iterations = 0
convergence_check = [1]
while (max(convergence_check)>.0001):
    iterations += 1
    psi_new = np.zeros((I_m1+I_m2,2))
    for i in range(I_m1+I_m2-1):
        # Sweep from left to right (from vacuum boundary to vacuum boundary)
        if (i<I_m1):
            q_left = .5*Sigma_s_m1*flux_guess[i+1] + 0.5*Q_m1
            psi_new[i+1,0] = (1+(Sigma_t_m1*h_m1)/(2*mu[0]))**(-1)*(psi_new[i,0]*(1-(h_m1*Sigma_t_m1)/(2*mu[0])) + h_m1*q_left/mu[0])
        else:
            q_left = .5*Sigma_s_m2*flux_guess[i+1] + 0.5*Q_m2
            psi_new[i+1,0] = (1+(Sigma_t_m2*h_m2)/(2*mu[0]))**(-1)*(psi_new[i,0]*(1-(h_m2*Sigma_t_m2)/(2*mu[0])) + h_m2*q_left/mu[0])
        # Sweep from right to left (from vacuum boundary to vacuum boundary)
        if (i<I_m2):
            q_right = .5*Sigma_s_m2*flux_guess[-(i+2)] + 0.5*Q_m2
            psi_new[-(i+2),1] = ((Sigma_t_m2*h_m2)/(2*(mu[1]))-1)**(-1)*(h_m2*q_right/mu[1] - psi_new[-(i+1),1]*((h_m2*Sigma_t_m2)/(2*(mu[1]))+1))
        else:
            q_right = .5*Sigma_s_m1*flux_guess[-(i+2)] + 0.5*Q_m1
            psi_new[-(i+2),1] = ((Sigma_t_m1*h_m1)/(2*(mu[1]))-1)**(-1)*(h_m1*q_right/mu[1] - psi_new[-(i+1),1]*((h_m1*Sigma_t_m1)/(2*(mu[1]))+1))
    flux_new = np.zeros((I_m1+I_m2,2))
    # print(flux_new)
    for i in range(len(psi_new[0,:])):
        flux_new[:,i] = w[i]*psi_new[:,i]
    flux_new = np.sum(flux_new,axis=1).tolist()
    #put convergence check here
    convergence_check = []
    for i in range(len(flux_new)):
        convergence_check.append(abs((flux_new[i]-flux_guess[i])/flux_new[i]))
    psi_guess = np.zeros((I_m1+I_m2,2))
    flux_guess = np.zeros((I_m1+I_m2,1))
    for i in range(I_m1+I_m2):
        for j in range(2):
            psi_guess[i,j] = psi_new[i,j]
        flux_guess[i] = flux_new[i]
            
print(iterations)
print('Max flux is limited in material 1 by: ' + str(Q_m1/(Sigma_t_m1-Sigma_s_m1)))
print('Max flux is limited in material 2 by: ' + str(Q_m2/(Sigma_t_m2-Sigma_s_m2)))
x = np.zeros(I_m1+I_m2)
for i in range(I_m1+I_m2):
    if (i<I_m1):
        x[i] = i*h_m1
    else:
        x[i] = width_m1 + (i-I_m1)*h_m2
plt.plot(x,flux_new)