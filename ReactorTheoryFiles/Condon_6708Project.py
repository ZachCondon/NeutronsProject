# Project for NE 6708
# Zachary Condon

#----Project prompt----#
# Write a multigroup 1D diffusion solver with vacuum boundary on the left and 
# reflecting boundary on the right. It must account for heterogeneous media
# with at least two materials. Optional features (if you want to include them):
# eigenvalue calculations, general boundary conditions, adaptive mesh 
# refinement, parallel computing, ... (?)

import numpy as np
import matplotlib.pyplot as plt
from gS import gaussSeidel
import time
from scenario1 import scenario1
from scenario2 import scenario2
from scenario3 import scenario3
from scenario4 import scenario4
from scenario5 import scenario5

print('-------Welcome to Zach Condon\'s NE6708 project--------')
print('This project is my 1D, multigroup, heterogeneous diffusion solver.')
print('This solver simulates two materials and solves for the flux using finite volumes.')
print('------------------------------------------------------')
print('There are three test cases to verify that my solver is working correctly')
print('*For more information regarding any case, see the corresponding file')
print('------------------------------------------------------')
print('------------------------------------------------------')
print('Scenario 1 is a homogeneous material, with one neutron group, and vacuum boundaries on both sides.')
print('Scenario 2 is a homogeneous material, with one neutron group, and a vacuum on left and reflective on right.')
print('Scenario 3 is a heterogeneous material, with one neutron group, and a vacuum on left and reflective on right.')
print('Scenario 4 is the heterogeneous, two-group, vacuum-on-left, reflective-on-right solution.')
print('Scenario 5 is identical to Scenario 4, except that the sources are defined per cell.')
print('------------------------------------------------------')
print()

# Set the Scenario number here:
scenario = 5
start = time.time()

if scenario == 1:   # Homogenous, 1 group, with vacuum boundaries on either side
    n = 60              # the number of finite volumes (cells)
    w = 10              # [cm], width of the material, 10 cm on either side of x-axis
    Sig_ai_frac = 0.1   # The ratio of absorption cross section to total cross section
    Sig_ti = 1          # [1/cm], macroscopic total cross section
    Si = 1              # [n/cm^3], source strength per cell
    scenario1(n,w,Sig_ai_frac,Sig_ti,Si)
    
elif scenario == 2: # Homogenous, 1 group, vacuum on left, reflective on right
    n = 60              # number of finite volumes (cells)
    w = 10              # [cm], width of the block, starting from x=0 to x=-w
    Sig_ai_frac = 0.1   # The ratio of absorption cross section to total cross section
    Sig_ti = 1          # [1/cm], macroscopic total cross section
    Si = 1              # [n/cm^3], source strength per cell
    scenario2(n,w,Sig_ai_frac,Sig_ti,Si)

elif scenario == 3: # Heterogenous, 1 group, vacuum on left, reflective on right

    #----- Material 1 -----
    n1 = 20             # number of cells in material 1
    w1 = 5              # [cm], width of material one, located from x=-w2 to x=-w2-w1
    Sig_ai_frac1 = 0.01  # The ratio of absorption cross section to total cross section
    Sig_ti_m1 = 1       # [1/cm], macroscopic total cross section
    Si_m1 = 1           # [n/cm^3], source strength per cell
    
    #----- Material 2 -----
    n2 = 20             # number of cells in material 2
    w2 = 5              # [cm], width of material 2, located from x=0 to x=-w2
    Sig_ai_frac2 = 0.2  # The ratio of absorption cross section to total cross section
    Sig_ti_m2 = 3       # [1/cm], macroscopic total cross section
    Si_m2 = 2           # [n/cm^3], source strength per cell
    
    scenario3(n1,w1,Sig_ai_frac1,Sig_ti_m1,Si_m1,n2,w2,Sig_ai_frac2,Sig_ti_m2,Si_m2)    
    
elif scenario == 4: # Heterogenous, 2 group, vacuum on left, reflective on right
    #------------------------- DEFINE THE PARAMETERS -------------------------
    #----- Geometry -----
    n1 = 20                         # number of cells in material 1
    n2 = 20                         # number of cells in material 2
    w1 = 5                          # [cm], width of material 1
    w2 = 5                          # [cm], width of material 2
    
    
    #----- Material 1 -- Group 1 -----
    Sig_ti_m1_g1 = 1                # [1/cm], macroscopic total cross section
    Sig_ai_m1_g1_frac = 0.01         # [1/cm], macroscopic absorption cross section
    Sig_s12_m1_frac = 0.2           # [1/cm], macroscopic scatter cross section from group 1->2
    Si_m1_g1 = 2                    # [n/cm^3], source strength per cell
    
    #----- Material 2 -- Group 1 -----
    Sig_ti_m2_g1 = 1                # [1/cm], macroscopic total cross section
    Sig_ai_m2_g1_frac = 0.002         # [1/cm], macroscopic absorption cross section
    Sig_s12_m2_frac = 0.4           # [1/cm], macroscopic scatter cross section from group 1->2
    Si_m2_g1 = 1                    # [n/cm^3], source strength per cell
    
    #----- Material 1 -- Group 2 -----
    Sig_ti_m1_g2 = 1                # [1/cm], macroscopic total cross section
    Sig_ai_m1_g2_frac = 0.2 # [1/cm], macroscopic absorption cross section
    Si_m1_g2 = 0                    # [n/cm^3], source strength per cell
    
    #----- Material 2 -- Group 2 -----
    Sig_ti_m2_g2 = 4                # [1/cm], macroscopic total cross section
    Sig_ai_m2_g2_frac = 0.1         # [1/cm], macroscopic absorption cross section
    Si_m2_g2 = 1                    # [n/cm^3], source strength per cell
    
    scenario4(n1,n2,
              w1,w2,
              Sig_ti_m1_g1,Sig_ai_m1_g1_frac,Sig_s12_m1_frac,Si_m1_g1,
              Sig_ti_m2_g1,Sig_ai_m2_g1_frac,Sig_s12_m2_frac,Si_m2_g1,
              Sig_ti_m1_g2,Sig_ai_m1_g2_frac,Si_m1_g2,
              Sig_ti_m2_g2,Sig_ai_m2_g2_frac,Si_m2_g2)
    
elif scenario == 5: # Same as scenario 4 with the ability to place point sources
    #------------------------- DEFINE THE PARAMETERS -------------------------
    #----- Geometry -----
    n1 = 20                         # number of cells in material 1
    n2 = 20                         # number of cells in material 2
    w1 = 5                          # [cm], width of material 1
    w2 = 5                          # [cm], width of material 2
    
    #----- Sources -----
    S_m1_g1 = np.zeros(n1)
    S_m2_g1 = np.zeros(n2)
    S_m1_g2 = np.zeros(n1)
    S_m2_g2 = np.zeros(n2)
    # Define a source anywhere in the array from index 0 to n-1
    # For example: S_m2_g1[12] = 1 puts a group 1 source of strength 1 in the 
    # 12th cell in material 2. this can be done for any up to all cells.
    S_m2_g1[12] = 1
    # S_m1_g1[14] = 1
    S_m1_g2[4] = 1
    
    S_g1 = np.concatenate((S_m1_g1,S_m2_g1))
    S_g2 = np.concatenate((S_m1_g2,S_m2_g2))
    
    #----- Material 1 -- Group 1 -----
    Sig_ti_m1_g1 = 1                # [1/cm], macroscopic total cross section
    Sig_ai_m1_g1_frac = 0.001         # [1/cm], macroscopic absorption cross section
    Sig_s12_m1_frac = 0.1           # [1/cm], macroscopic scatter cross section from group 1->2
    
    #----- Material 2 -- Group 1 -----
    Sig_ti_m2_g1 = 1                # [1/cm], macroscopic total cross section
    Sig_ai_m2_g1_frac = 0.002         # [1/cm], macroscopic absorption cross section
    Sig_s12_m2_frac = 0.1           # [1/cm], macroscopic scatter cross section from group 1->2
    
    #----- Material 1 -- Group 2 -----
    Sig_ti_m1_g2 = 1                # [1/cm], macroscopic total cross section
    Sig_ai_m1_g2_frac = 0.2         # [1/cm], macroscopic absorption cross section
    
    #----- Material 2 -- Group 2 -----
    Sig_ti_m2_g2 = 1                # [1/cm], macroscopic total cross section
    Sig_ai_m2_g2_frac = 0.1         # [1/cm], macroscopic absorption cross section
    
    scenario5(n1,n2,
              w1,w2,
              S_g1,S_g2,
              Sig_ti_m1_g1,Sig_ai_m1_g1_frac,Sig_s12_m1_frac,
              Sig_ti_m2_g1,Sig_ai_m2_g1_frac,Sig_s12_m2_frac,
              Sig_ti_m1_g2,Sig_ai_m1_g2_frac,
              Sig_ti_m2_g2,Sig_ai_m2_g2_frac,)
    
end = time.time()
print(f'This took {end-start} seconds')
