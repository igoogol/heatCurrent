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
FermiEnergy = sp.linspace(-10,10,200)
Gh = []
for Ef in FermiEnergy: 
#     realGh1 = Gamma1*Gamma2/(8*sp.pi*Gamma*Omega)*(-4*Omega*sp.arctan(2*(E0-Ef)/Gamma) \
#                                                +(4*(Ef-E0)+2*Omega)*sp.arctan(2*(E0-Ef-Omega)/Gamma) \
#                                                +(4*(E0-Ef)+2*Omega)*sp.arctan(2*(E0-Ef+Omega)/Gamma) \
#                                                -Gamma*sp.log(4*(E0-Ef)**2 + Gamma**2 + 8*(E0-Ef)*Omega + 4*Omega**2) \
#                                                +Gamma*sp.log(4*(E0-Ef)**2 + Gamma**2 - 8*(E0-Ef)*Omega + 4*Omega**2) )
    realGh1 = Gamma1*Gamma2/(8*sp.pi*(Gamma**2+Omega**2)*Omega)*( \
                                               -(4*Gamma*(E0-Ef))*sp.arctan(2*(E0-Ef-Omega)/Gamma) \
                                               +(4*Gamma*(E0-Ef))*sp.arctan(2*(E0-Ef+Omega)/Gamma) \
                                               +(4*Omega*(E0-Ef))*sp.log(4*(E0-Ef)**2 + Gamma**2) \
                                               +(-2*Omega*(E0-Ef)-Omega**2-Gamma**2)  \
                                                 *sp.log(4*(E0-Ef)**2 + Gamma**2 + 8*(E0-Ef)*Omega + 4*Omega**2) \
                                               +(-2*Omega*(E0-Ef)+Omega**2+Gamma**2)  \
                                                 *sp.log(4*(E0-Ef)**2 + Gamma**2 - 8*(E0-Ef)*Omega + 4*Omega**2) )
    imagGh1 = Gamma1*Gamma2/(4*sp.pi*(Gamma**2+Omega**2)*Omega)*( \
                                                   -2*(Gamma**2+Omega**2)*sp.arctan(2*(E0-Ef)/Gamma) \
                                                   +(-2*Omega*(E0-Ef)+Omega**2+Gamma**2)*sp.arctan(2*(E0-Ef-Omega)/Gamma) \
                                                   +(+2*Omega*(E0-Ef)+Omega**2+Gamma**2)*sp.arctan(2*(E0-Ef+Omega)/Gamma) \
                                                   -2*Gamma*(E0-Ef)*sp.log(4*(E0-Ef)**2 + Gamma**2) \
                                                   +Gamma*(E0-Ef)*sp.log(4*(E0-Ef)**2 + Gamma**2 + 8*(E0-Ef)*Omega + 4*Omega**2) \
                                                   +Gamma*(E0-Ef)*sp.log(4*(E0-Ef)**2 + Gamma**2 - 8*(E0-Ef)*Omega + 4*Omega**2) )
#     realGh2 = -Gamma1*Gamma2/(4*Gamma*sp.pi)*(2*sp.arctan(2*(E0-Ef)/Gamma) \
#                                               -sp.arctan(2*(E0-Ef-Omega)/Gamma) \
#                                               -sp.arctan(2*(E0-Ef+Omega)/Gamma))
#     imagGh2 = Gamma1*Gamma2/(8*Gamma*sp.pi)*(sp.log(4*(E0-Ef)**2 + Gamma**2 + 8*(E0-Ef)*Omega + 4*Omega**2) \
#                                              -sp.log(4*(E0-Ef)**2 + Gamma**2 - 8*(E0-Ef)*Omega + 4*Omega**2))
    Gh.append(realGh1 + imagGh1*1j)
    
# print Gh

plt.plot(FermiEnergy,sp.real(Gh),label='real',color='blue')
plt.plot(FermiEnergy,sp.imag(Gh),label='imag',color='red')
plt.xlabel("$E_f/\Gamma$")
plt.ylabel("$G^h$")
plt.title("$\Omega = 1 \Gamma$, $\Gamma_L = \Gamma_R = 0.5 \Gamma$")
plt.legend()
plt.show()






