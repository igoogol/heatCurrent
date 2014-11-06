import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# plot Gh-Ef
GhNoU = np.load("Gh-Ef-noU.npy")
Gh = np.load("Gh-Ef.npy")
FermiEnergy = np.load("Ef.npy")
  
plt.plot(FermiEnergy,sp.real(Gh),label='real',color='blue')
plt.plot(FermiEnergy,sp.imag(Gh),label='imag',color='red')
plt.plot(FermiEnergy,sp.real(GhNoU),'.',markersize = 3.0,label='real',color='blue')
plt.plot(FermiEnergy,sp.imag(GhNoU),'.',markersize = 3.0,label='imag',color='red')
plt.xlabel("$(E_f-E_0)/\Gamma$")
plt.ylabel("$G^h$")
plt.title("")
plt.legend()
plt.show()

# # plot Gh-Omega
# GhNoU = np.load("Gh-Omega-NoU.npy")
# Gh = np.load("Gh-Omega.npy")
# Frequency = np.load("Omega.npy")
#  
# plt.plot(Frequency,sp.real(Gh),label='real',color='blue')
# plt.plot(Frequency,sp.imag(Gh),label='imag',color='red')  
# plt.plot(Frequency,sp.real(GhNoU),'.',markersize = 3.0,label='real',color='blue')
# plt.plot(Frequency,sp.imag(GhNoU),'.',markersize = 3.0,label='imag',color='red') 
# plt.xlabel("$\Omega/\Gamma$")
# plt.ylabel("$G^h$")
# plt.title("")
# plt.legend()
# plt.show()