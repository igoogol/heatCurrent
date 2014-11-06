import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

GhNoU = np.load("Gh-Ef-noU.npy")
Gh = np.load("Gh-Ef.npy")
FermiEnergy = np.load("Ef.npy")

plt.plot(FermiEnergy,sp.real(Gh),label='real',color='blue')
plt.plot(FermiEnergy,sp.imag(Gh),label='imag',color='red')
plt.plot(FermiEnergy,sp.real(GhNoU),'--',label='real',color='blue')
plt.plot(FermiEnergy,sp.imag(GhNoU),'--',label='imag',color='red')
plt.xlabel("$(E_f-E_0)/\Gamma$")
plt.ylabel("$G^h$")
plt.title("")
plt.legend()
plt.show()