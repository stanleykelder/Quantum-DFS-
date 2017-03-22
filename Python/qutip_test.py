import numpy as np
import time
from qutip import *
import scipy.sparse as sp

start = time.time()

# Takes 0 or 1 as input, returns corresponding bra
def bra(spin):
	return basis(2, spin).dag()

# Takes 0 or 1 as input, returns corresponding ket
def ket(spin):
	return basis(2, spin)

# Takes number of dimensions as input, returns corresponding identity matrix
# def identity(size):
	# return(qeye(size))

# Takes 0 or 1 as input. Returns density matrix
def rho_1(spin):
	return ket(spin) * bra(spin)

# Takes z or x as input, returns corresponding total Spin matrix
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

# Takes as input a value for k (0, 1 or 2) and a time t. Returns corresponding Error
def E(k, t):
	E2 = np.sqrt(t) * spintotal("z")
	E1 = np.sqrt(t) * spintotal("x")
	E0 = tensor(identity(2), identity(2), identity(2)) - ((E1 ** 2) / 2) - ((E2 ** 2) / 2)
	
	if k == 2:
		return E2 
	if k == 1:
		return E1 
	if k == 0:
		return E0 

# Takes 0 or 1 as input, returns corresponding initial rho encoded
def rho_0(spin):
	return .5 * (tensor((ket(spin) * bra(spin)), identity(2)))

# 
def epsilon(spin, dt):
	epsilon = (E(0, dt) * direct_sum(rho_0(spin), 4) * E(0, dt).dag() +
			   E(1, dt) * direct_sum(rho_0(spin), 4) * E(1, dt).dag() +
			   E(2, dt) * direct_sum(rho_0(spin), 4) * E(2, dt).dag())		
	return epsilon

##########################################################################
# Use numpy to calculate direct sum. input qutip matrix and n3. Returns 
def direct_sum(rho, n3):
	zero = np.zeros([n3, n3])
	rho = rho.full()

	rho_encoded = np.zeros(np.add(rho.shape,zero.shape))
	rho_encoded[:rho.shape[0], :rho.shape[1]] = rho
	rho_encoded[rho.shape[0]:, rho.shape[1]:] = zero
	return Qobj(rho_encoded, [[2, 2, 2], [2, 2, 2]])
###########################################################################

def projection(M):
	return Qobj(M[:4, :4], [[2, 2], [2, 2]])

# print(direct_sum(rho_0(0), 4))
# print(projection(direct_sum(rho_0(0), 4)).ptrace(0))
# print((projection(direct_sum(rho_0(0), 4)).ptrace(0) * rho_1(0)).tr())


Epsilon = epsilon(0, 0.04)
Projection = projection(Epsilon)
Ptrace = Projection.ptrace(0)
Fidelity = (Projection.ptrace(0) * rho_1(0)).tr()

# print(Epsilon)
# print(Projection.ptrace(0))
# print(Fidelity)


eind = time.time()
print(eind-start)