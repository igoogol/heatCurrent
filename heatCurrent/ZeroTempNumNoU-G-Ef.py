'''
Created on Aug 27, 2014
This code is used to calculate the heat conductance.
Plot: Gh v.s Ef

@author: googol
'''
__author__ = 'googol'

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

# Parameters
Gamma1 = 0.5
Gamma2 = 0.5
Gamma = Gamma1 + Gamma2
SigmaR = -1j*Gamma/2.0
SigmaA = +1j*Gamma/2.0
Omega = 1.0     # response frequency
#kT = 1.0        # 25.6  # room temperature
E0 = 0.0        # single energy level of the QD
NE = 2000       # number of energy points
Ueq = 0.0       # equilibrium potential
 
# Energy grid
FermiEnergy = sp.linspace(-10,10,200) # 200
Gh = []
for Ef in FermiEnergy: 
    Emin = Ef - Omega
    Emax = Ef   
    E = np.linspace(Emin,Emax,NE)
    dE = E[1] - E[0]
#     # one way to declare 2-D matrix
#     Gr = [[0 for col in range(1)]for row in range(NE)]
#     Gr[2][0] = 3
#     print Gr
#     # another way to declare 2-D matrix
#     Gr = [[]for row in range(NE)]
#     Gr[1].append([1,2])
#     print Gr
    Gr = [[]for row in range(NE)]
    Ga = [[]for row in range(NE)]
    GrBar = [[]for row in range(NE)]
    GaBar = [[]for row in range(NE)]
    g1_12 = [[]for row in range(NE)]
    g2_12 = [[]for row in range(NE)]
    gh1_12 = [[]for row in range(NE)]
    gh2_12 = [[]for row in range(NE)]
    for ii in np.arange(NE):
        Gr[ii] = 1.0/(E[ii] - E0 - SigmaR -Ueq)
        Ga[ii] = 1.0/(E[ii] - E0 - SigmaA -Ueq)
        GrBar[ii] = 1.0/(E[ii] + Omega - E0 - SigmaR -Ueq)
        GaBar[ii] = 1.0/(E[ii] + Omega - E0 - SigmaA -Ueq)
        # g1_12[ii] = - GrBar[ii]*Gamma2*Ga[ii]*Gamma1 + Omega*GrBar[ii]*(Gamma2/Gamma)*Ga[ii]*1j*Gamma1
        g1_12[ii] = - GrBar[ii]*Gamma2*Ga[ii]*Gamma1 
        # g2_12[ii] = (Ga[ii] + GrBar[ii])*1j*Gamma1*(0-Gamma2/Gamma)
        g2_12[ii] = 0
        gh1_12[ii] = 1/Omega*(E[ii]+Omega/2.0-Ef)*g1_12[ii]
        gh2_12[ii] = 1/Omega*(Omega/2.0)*g2_12[ii]  
    Gh1 = -(dE/(2.0*sp.pi))*sp.sum(gh1_12)
    Gh2 = +(dE/(2.0*sp.pi))*sp.sum(gh2_12)
    Gh.append(Gh1+Gh2)
    
# save data
np.save("Gh-Ef-noU.npy",Gh)
np.save("Ef.npy",FermiEnergy)
# np.save("Gh-Ef.npy",Gh)
# np.save("Ef.npy",FermiEnergy)

plt.plot(FermiEnergy,sp.real(Gh),label='real',color='blue')
plt.plot(FermiEnergy,sp.imag(Gh),label='imag',color='red')
plt.xlabel("$E_f/\Gamma$")
plt.ylabel("$G^h$")
plt.title("$\Omega = 1 \Gamma$, $\Gamma_L = \Gamma_R = 0.5 \Gamma$")
plt.legend()
plt.show()






