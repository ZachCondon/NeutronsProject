# Source iteration transport
import numpy as np
import matplotlib.pyplot as plt

def initialize_problem_values():
    # This is a simple initialization function to set certain values and 
    # establish the quadrature angles and weights.
    M = 2           # number of materials
    G = 3           # number of groups
    N = 16           # number of angles (quadrature)
    # The following lines establish the quadrature and determine the angles (mu)
    # and weights (w). The reverse commands are used to follow the same order
    # of angles and weights later in the code.
    SN_quad = np.polynomial.legendre.leggauss(N)
    SN_quad = np.array(SN_quad).tolist()
    mu = SN_quad[0]
    mu.reverse()
    w = SN_quad[1]
    w.reverse()
    # The line below sets whether the right side of the simulation should have
    # a reflective boundary. Setting this equality to any other variable (except
    # for 1, since Python interprets that as True) will make the right side a vacuum.
    reflective = 0
    return M,G,N,mu,w,reflective

def define_cross_sections(numMaterials, numGroups):
    # If G is the number of groups and M is the number of materials, Sigma_t is
    # a two-dimensional matrix of size MxG. Sigma_s is a MxGxG matrix. The
    # first index refers to the material, the second index refers to the group,
    # and the following index points to which group it is scattering to. For 
    # this project, the upscattering cross-section is always zero.
    
    # Any value greater than zero. The size of the simulation will scale based
    # on the largest mean free path within each material. 
    # Theoretically, zero should be an acceptable input, but I haven't gotten
    # it to work. The function that determines the size of the mfp will ignore
    # the zero value 
    Sigma_t = np.zeros((numMaterials,numGroups))
    # Sigma_t[material index, group index]
    Sigma_t[0,0],Sigma_t[0,1],Sigma_t[0,2] = 1,1,1  # Material 1 total cross-sections
    Sigma_t[1,0],Sigma_t[1,1],Sigma_t[1,2] = 0,0,0  # Material 2 total cross-sections

    
    Sigma_s = np.zeros((numMaterials,numGroups,numGroups))
    # Sigma_s[material index, group index, scatter-to-group index]
    # Material 1 scattering cross-sections
    Sigma_s[0,0,0] = 0.4    # G1->G1
    Sigma_s[0,0,1] = 0.3    # G1->G2
    Sigma_s[0,0,2] = 0.0    # G1->G3
    Sigma_s[0,1,0] = 0.0    # G2->G1
    Sigma_s[0,1,1] = 0.2    # G2->G2
    Sigma_s[0,1,2] = 0.1    # G2->G3
    Sigma_s[0,2,0] = 0.0    # G3->G1
    Sigma_s[0,2,1] = 0.0    # G3->G2
    Sigma_s[0,2,2] = 0.1    # G3->G3
    # Material 2 scattering cross-sections
    Sigma_s[1,0,0] = 0.0    # G1->G1
    Sigma_s[1,0,1] = 0.0    # G1->G2
    Sigma_s[1,0,2] = 0.0    # G1->G3
    Sigma_s[1,1,0] = 0.0    # G2->G1
    Sigma_s[1,1,1] = 0.0    # G2->G2
    Sigma_s[1,1,2] = 0.0    # G2->G3
    Sigma_s[1,2,0] = 0.0    # G3->G1
    Sigma_s[1,2,1] = 0.0    # G3->G2
    Sigma_s[1,2,2] = 0.0    # G3->G3
    return Sigma_t, Sigma_s

def define_source(numMaterials,numGroups):
    # This places all the source values into a two-dimensional matrix
    Q = np.zeros((numMaterials,numGroups))
    # Q[material index, group index]
    # Material 1 source strength
    Q[0,0] = 0.1    # G1
    Q[0,1] = 0.1    # G2
    Q[0,2] = 0.1    # G3
    # Material 2 source strength
    Q[1,0] = 0.0    # G1
    Q[1,1] = 0.0    # G2
    Q[1,2] = 0.0    # G3
    return Q

def define_material(Sigma_t,M):
    # This function defines the width of each material (in mean free paths using
    # the minimum total cross-section). This ensures that the group with the 
    # largest mean free path (since mfp=1/Sigma_t) will be able to be accurately
    # represented by the code. 
    max_mfp = np.zeros((M))
    for m in range(M):
        mat_Sigma_t = Sigma_t[m]
        # The following 4 lines handle the error that would arise if there is
        # a void. If Sigma_t was set to be zero, the nonzero function would
        # cause an error. The exception sets the max_mfp to be one for each group
        # and in the width definition below, the number 10 can be modified to 
        # change the width of the void.
        try:
            max_mfp[m] = 1/np.min(mat_Sigma_t[np.nonzero(mat_Sigma_t)])
        except:
            max_mfp[m] = 1
        # if np.min(mat_Sigma_t[np.nonzero(mat_Sigma_t)])==0:
        #     max_mfp[m] = 1
        # else:
        #     max_mfp[m] = 1/np.min(mat_Sigma_t[np.nonzero(mat_Sigma_t)])
    width = [10*max_mfp[0],10*max_mfp[1]] # width array, in order of materials from left to right
    I = [1000,1000]                             # number of mesh points, same order as above
    h = [width[0]/I[0], width[1]/I[1]]          # mesh size, same order as above
    return h,I,width

def calculate_q(material_index, group, position_index, flux_guess, source, Sigma_s, G, direction):
    # This function calculates q, the source value for the transport equation.
    # The needed physical variables are:
        # the previously calculated flux (flux_guess),
        # the source (source), 
        # the scattering cross-section matrix (Sigma_s),
        # the total number of groups (G).
    # The variables that will be iterated upon when using this function are:
        # material_index,   defines which material when looking up values
        # group,            defines the current group when looking up values
        # position_index,   denotes the current position for the equations below
        # direction,        used to signify whether the calculations are performed left-to-right or vice-versa.
    if (direction == 0):    # direction == 0 signifies marching from left-to-right
        in_scatter = 0
        for scatter_from_group in range(G):
            in_scatter += Sigma_s[material,scatter_from_group,group]*flux_guess[scatter_from_group,position_index+1]
        q = 0.5*(in_scatter + Q[material,group])
    elif (direction == 1):
        in_scatter = 0
        for scatter_from_group in range(G):
            in_scatter += Sigma_s[material,scatter_from_group,group]*flux_guess[scatter_from_group,-(position_index+2)]
        q = 0.5*(in_scatter + Q[material,group])
    return q

def march_psi(material_index, group, position_index, angle_index, psi_new, Sigma_t, h, mu, q, direction):
    # This function takes the angular flux variable at the ith position and
    # uses it to calculate the angular flux at the (i+1)th position. It is based
    # on the transport equation for a slab, using the diamond difference method.
    if (direction==0):
        # The following if statement handles a void material. If the total
        # cross-section is defined as 0, then it will keep the flux at the same
        # value throughout.
        if Sigma_t[material,group]==0:
            next_psi = psi_new[group,i,n]
        else:
            next_psi = (1+(Sigma_t[material,group]*h[material])/(2*mu[n]))**(-1)*(psi_new[group,i,n]*(1-(h[material]*Sigma_t[material,group])/(2*mu[n])) + h[material]*q/mu[n])
    else:
        # The following if statement handles a void material. If the total
        # cross-section is defined as 0, then it will keep the flux at the same
        # value throughout.
        if Sigma_t[material,group]==0:
            next_psi = psi_new[group,-(i+1),-(n+1)]
        else:
            next_psi = ((Sigma_t[material,group]*h[material])/(2*(mu[-(n+1)]))-1)**(-1)*(h[material]*q/mu[-(n+1)] - psi_new[group,-(i+1),-(n+1)]*((h[material]*Sigma_t[material,group])/(2*(mu[-(n+1)]))+1))
    return next_psi

def calculate_flux(G,I,N,w,psi_new,):
    # The scalar flux is calculated directly from the angular flux. First a new
    # matrix (weighted_psi) is defined with the exact same dimensions as the
    # psi_new matrix. Then, this matrix is filled in for every group and position
    # with the appropriate weights for each angle. Finally, at each position and
    # for each group,the weighted angular flux is summed up to provide the 
    # scalar flux for each group at each position.
    weighted_psi = np.zeros((G,sum(I),N))
    for n in range(N):
        weighted_psi[:,:,n] = w[n]*psi_new[:,:,n]
    flux_new = np.zeros((G,sum(I)))
    for group in range(G):
        for i in range(sum(I)):
            flux_new[group,i] = sum(weighted_psi[group,i])
    return flux_new

def update_guesses(G,I,N,psi_new,flux_new):
    # This function simply updates the variables psi_guess and flux_guess with
    # the newly calculated psi and flux after a whole transport sweep is
    # completed. I have this function because I've had issues with how python
    # calculates and overrides variables that point to matrices. This is my
    # best method for ensuring the "guess" matrices are completely separate from 
    # the "new" matrices.
    psi_guess = np.zeros((G,sum(I),N))
    flux_guess = np.zeros((G,sum(I)))
    for i in range(sum(I)):
        for group in range(G):
            flux_guess[group,i] = flux_new[group,i]
            for n in range(int(N/2)):
                psi_guess[group,i,n] = psi_new[group,i,n]
    return psi_guess,flux_guess

def determine_material(i,I,direction):
    # This function simply determines the material at the current position, i.
    if (direction == 0):
        if (i<I[0]):
            material = 0
        else:
            material = 1
    elif (direction == 1):
        if (i<I[1]):
            material = 1
        else:
            material = 0
    return material

# -----Initialize all global variables-----
M,G,N,mu,w,reflective = initialize_problem_values()
Sigma_t, Sigma_s = define_cross_sections(M, G)
Q = define_source(M,G)
h,I,width = define_material(Sigma_t,M)
# -----Initialize the guesses for psi and flux-----
psi_guess = np.zeros((G,sum(I),N))
flux_guess = calculate_flux(G,I,N,w,psi_guess)

iterations = 0                          # This is used to show how many iterations the code took.
# After each sweep, convergence check is an array of size sum(I) containing the
# change in the flux at each position. It is calculated using (new-old)/new
# If any position has convergence check value greater than .0001, it will run again.
convergence_check = [1]
while (max(convergence_check)>.0001):
    iterations += 1
    psi_new = np.zeros((G,sum(I),N))    # Reset psi_new
    for n in range(int(N/2)):           # Sweep through all angles
        # the range is N/2 because the loop will iterate for the positive angle
        # followed by the corresponding negative angle. In effect, it sweeps from
        # left-to-right then from right-to-left.
        # -----Left-to-right sweep-----
        for i in range(sum(I)-1):
            direction = 0               # the number 0 signifies a left-to-right sweep
            material = determine_material(i, I, direction)
            for group in range(G):
                q = calculate_q(material,group,i,flux_guess,Q,Sigma_s,G,direction)
                psi_new[group,i+1,n] = march_psi(material,group,i,n,psi_new,Sigma_t,h,mu,q,direction)
        # -----Reflective boundary determination-----
        if reflective:
            # If reflective is true, the next two lines set the reflective
            # boundary condition for the sweep from right-to-left.
            for group in range(G):
                psi_new[group,-1,-(n+1)] = psi_new[group,-1,n]
        # -----Right-to-left sweep-----
        for i in range(sum(I)-1):
            direction = 1               # the number 1 signifies a right-to-left sweep
            material = determine_material(i, I, direction)
            for group in range(G):
                q = calculate_q(material,group,i,flux_guess,Q,Sigma_s,G,direction)
                psi_new[group,-(i+2),-(n+1)] = march_psi(material,group,i,n,psi_new,Sigma_t,h,mu,q,direction)
    
    flux_new = calculate_flux(G,I,N,w,psi_new)  # calculates the new flux

    # -----Convergence check-----
    # The following lines establish an array with all the convergence check 
    # values. If any value in the array is greater than .0001, the transport
    # sweep will run again.
    convergence_check = []
    for i in range(sum(I)):
        convergence_check.append(abs((sum(flux_new[:,i]) - sum(flux_guess[:,i]))/sum(flux_new[:,i])))
    
    # Reset psi and flux matrices.
    psi_guess,flux_guess = update_guesses(G,I,N,psi_new,flux_new)
            
print("Number of iterations: " + str(iterations))

# Graph the results
x = np.zeros(sum(I))
for i in range(sum(I)):
    if (i<I[0]):
        x[i] = i*h[0]
    else:
        x[i] = width[0] + (i-I[0])*h[1]
plt.plot(x,flux_new[0],x,flux_new[1],x,flux_new[2])
plt.legend(['Group 1 (fast)', 'Group 2 (medium)','Group 3 (slow)'])
plt.title('Flux using S' + str(N) + ' method')
plt.xlabel('Position (units of mean free path)')