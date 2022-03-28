# Source iteration transport
import numpy as np
import matplotlib.pyplot as plt


S2_quad = np.polynomial.legendre.leggauss(2)
S2_quad = np.array(S2_quad).tolist()
mu = S2_quad[0]
mu.reverse()
w = S2_quad[1]
w.reverse()
reflective = True

# Material one, fast group
Sigma_t_g1_m1 = 1               # Total cross-section
Sigma_s_g1_m1 = 0.1             # Scatter within group cross-section
Sigma_s_g1_m1_12 = 0.0          # Scatter from group 1 to 2
Q_g1_m1 = .1                    # Group 1 source
# Material two, fast group
Sigma_t_g1_m2 = 1               # Total cross-section
Sigma_s_g1_m2 = 0.5             # Scatter within group cross-section
Sigma_s_g1_m2_12 = 0.0          # Scatter from group 1 to 2
Q_g1_m2 = .1                    # Group 1 source
# Material one, slow group
Sigma_t_g2_m1 = 1               # Total cross-section
Sigma_s_g2_m1 = 0.1             # Scatter within group cross-section
Q_g2_m1 = .1                    # Group 2 source
# Material two, fast group
Sigma_t_g2_m2 = 1               # Total cross-section
Sigma_s_g2_m2 = 0.0             # Scatter within group cross-section
Q_g2_m2 = .1                    # Group 2 source

# Material dimensions
width_m1 = 10
I_m1 = 1000
h_m1 = width_m1/I_m1
width_m2 = 10
I_m2 = 100
h_m2 = width_m2/I_m2

psi_guess_g1 = np.zeros((I_m1+I_m2,2))
flux_guess_g1 = np.sum(psi_guess_g1,axis=1).tolist()
psi_guess_g2 = np.zeros((I_m1+I_m2,2))
flux_guess_g2 = np.sum(psi_guess_g2,axis=1).tolist()

iterations = 0
convergence_check = [1]
while (max(convergence_check)>.0001):
    iterations += 1
    psi_new_g1 = np.zeros((I_m1+I_m2,2))
    psi_new_g2 = np.zeros((I_m1+I_m2,2))
    for i in range(I_m1+I_m2-1):
        # Sweep from left to right (from vacuum boundary to vacuum boundary)
        if (i<I_m1):
            q_left_g1 = .5*Sigma_s_g1_m1*flux_guess_g1[i+1] + 0.5*Q_g1_m1
            psi_new_g1[i+1,0] = (1+(Sigma_t_g1_m1*h_m1)/(2*mu[0]))**(-1)*(psi_new_g1[i,0]*(1-(h_m1*Sigma_t_g1_m1)/(2*mu[0])) + h_m1*q_left_g1/mu[0])
            q_left_g2 = .5*Sigma_s_g1_m1_12*flux_guess_g1[i+1] + .5*Sigma_s_g2_m1*flux_guess_g2[i+1] + 0.5*Q_g2_m1
            psi_new_g2[i+1,0] = (1+(Sigma_t_g2_m1*h_m1)/(2*mu[0]))**(-1)*(psi_new_g2[i,0]*(1-(h_m1*Sigma_t_g2_m1)/(2*mu[0])) + h_m1*q_left_g2/mu[0])
        else:
            q_left_g1 = .5*Sigma_s_g1_m2*flux_guess_g1[i+1] + 0.5*Q_g1_m2
            psi_new_g1[i+1,0] = (1+(Sigma_t_g1_m2*h_m2)/(2*mu[0]))**(-1)*(psi_new_g1[i,0]*(1-(h_m2*Sigma_t_g1_m2)/(2*mu[0])) + h_m2*q_left_g1/mu[0])
            q_left_g2 = .5*Sigma_s_g1_m2_12*flux_guess_g2[i+1] + .5*Sigma_s_g2_m2*flux_guess_g2[i+1] + 0.5*Q_g2_m2
            psi_new_g2[i+1,0] = (1+(Sigma_t_g2_m2*h_m2)/(2*mu[0]))**(-1)*(psi_new_g2[i,0]*(1-(h_m2*Sigma_t_g2_m2)/(2*mu[0])) + h_m2*q_left_g2/mu[0])
    if reflective:
        psi_new_g1[-1,1] = psi_new_g1[-1,0]
        psi_new_g2[-1,1] = psi_new_g2[-1,0]
    for i in range(I_m1+I_m2-1):
        # Sweep from right to left (from vacuum boundary to vacuum boundary)
        if (i<I_m2):
            q_right_g1 = .5*Sigma_s_g1_m2*flux_guess_g1[-(i+2)] + 0.5*Q_g1_m2
            psi_new_g1[-(i+2),1] = ((Sigma_t_g1_m2*h_m2)/(2*(mu[1]))-1)**(-1)*(h_m2*q_right_g1/mu[1] - psi_new_g1[-(i+1),1]*((h_m2*Sigma_t_g1_m2)/(2*(mu[1]))+1))
            q_right_g2 = .5*Sigma_s_g1_m2_12*flux_guess_g1[-(i+2)] + .5*Sigma_s_g2_m2*flux_guess_g2[-(i+2)] + 0.5*Q_g2_m2
            psi_new_g2[-(i+2),1] = ((Sigma_t_g2_m2*h_m2)/(2*(mu[1]))-1)**(-1)*(h_m2*q_right_g2/mu[1] - psi_new_g2[-(i+1),1]*((h_m2*Sigma_t_g2_m2)/(2*(mu[1]))+1))
        else:
            q_right_g1 = .5*Sigma_s_g1_m1*flux_guess_g1[-(i+2)] + 0.5*Q_g1_m1
            psi_new_g1[-(i+2),1] = ((Sigma_t_g1_m1*h_m1)/(2*(mu[1]))-1)**(-1)*(h_m1*q_right_g1/mu[1] - psi_new_g1[-(i+1),1]*((h_m1*Sigma_t_g1_m1)/(2*(mu[1]))+1))
            q_right_g2 = .5*Sigma_s_g1_m1_12*flux_guess_g1[-(i+2)] + .5*Sigma_s_g2_m1*flux_guess_g2[-(i+2)] + 0.5*Q_g2_m1
            psi_new_g2[-(i+2),1] = ((Sigma_t_g2_m1*h_m1)/(2*(mu[1]))-1)**(-1)*(h_m1*q_right_g2/mu[1] - psi_new_g2[-(i+1),1]*((h_m1*Sigma_t_g2_m1)/(2*(mu[1]))+1))
    flux_new_g1 = np.zeros((I_m1+I_m2,2))
    flux_new_g2 = np.zeros((I_m1+I_m2,2))
    # print(flux_new)
    for i in range(len(psi_new_g1[0,:])):
        flux_new_g1[:,i] = w[i]*psi_new_g1[:,i]
        flux_new_g2[:,i] = w[i]*psi_new_g2[:,i]
    flux_new_g1 = np.sum(flux_new_g1,axis=1).tolist()
    flux_new_g2 = np.sum(flux_new_g2,axis=1).tolist()
    # Convergence check
    convergence_check = []
    for i in range(len(flux_new_g1)):
        flux_new = flux_new_g1[i] + flux_new_g2[i]
        flux_guess = flux_guess_g1[i] + flux_guess_g2[i]
        convergence_check.append(abs((flux_new-flux_guess)/flux_new))
    psi_guess_g1 = np.zeros((I_m1+I_m2,2))
    psi_guess_g2 = np.zeros((I_m1+I_m2,2))
    flux_guess_g1 = np.zeros((I_m1+I_m2,1))
    flux_guess_g2 = np.zeros((I_m1+I_m2,1))
    for i in range(I_m1+I_m2):
        for j in range(2):
            psi_guess_g1[i,j] = psi_new_g1[i,j]
            psi_guess_g2[i,j] = psi_new_g2[i,j]
        flux_guess_g1[i] = flux_new_g1[i]
        flux_guess_g2[i] = flux_new_g2[i]
            
print("Number of iterations: " + str(iterations))
#print('Max flux is limited in material 1 by: ' + str(Q_g1_m1/(Sigma_t_g1_m1-Sigma_s_g1_m1)))
#print('Max flux is limited in material 2 by: ' + str(Q_g1_m2/(Sigma_t_g1_m2-Sigma_s_g1_m2)))
x = np.zeros(I_m1+I_m2)
for i in range(I_m1+I_m2):
    if (i<I_m1):
        x[i] = i*h_m1
    else:
        x[i] = width_m1 + (i-I_m1)*h_m2
plt.plot(x,flux_new_g1,x,flux_new_g2)
plt.legend(['Group 1 (fast)', 'Group 2 (slow)'])