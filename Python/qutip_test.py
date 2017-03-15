import numpy as np
import time
from qutip import *


start = time.time()
def bra(spin):
	return basis(2, spin).dag()

def ket(spin):
	return basis(2, spin)

def identity(size):
	return(qeye(size))

def spintotal(z_or_x):
	if z_or_x == "z":
		S = (tensor(sigmaz(), identity(2), identity(2)) +
			tensor(identity(2), sigmaz(), identity(2)) +
			tensor(identity(2), identity(2), sigmaz()))
		return S
	elif z_or_x == "x":
		S = (tensor(sigmax(), identity(2), identity(2)) +
			tensor(identity(2), sigmax(), identity(2)) +
			tensor(identity(2), identity(2), sigmax()))
		return S

def E(k, dt):
	E2 = np.sqrt(dt) * spintotal("z")
	E1 = np.sqrt(dt) * spintotal("x")
	E0 = tensor(identity(2), identity(2), identity(2)) - ((E1 ** 2) / 2) - ((E2 ** 2) / 2)
	
	if k == 2:
		return E2 
	if k == 1:
		return E1 
	if k == 0:
		return E0 

# print(E(1, 0.1))
# print(E(1, 0.1)**2)

def rho_0(spin):
	return .5 * (tensor((ket(spin) * bra(spin)), (ket(spin) * bra(spin)), identity(2)))

def epsilon(spin, dt):
	epsilon = (E(0, dt) * rho_0(spin) * E(0, dt).dag() +
			   E(1, dt) * rho_0(spin) * E(1, dt).dag() +
			   E(2, dt) * rho_0(spin) * E(2, dt).dag())		
	return epsilon

PartTrace = (epsilon(0, 0.4) * rho_0(0)).ptrace(2)
rho_1 = ket(0) * bra(0)

print(rho_0(0))
print((rho_1 * PartTrace).tr())



# eind = time.time()
# print(eind-start)