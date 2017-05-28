from __future__ import division
import numpy as np
import time
from qutip import *
from definitions_numpy import *
from numpy.core.umath_tests import inner1d


start = time.time()
print tensor(sigmaz(), identity(2), identity(2)) - tensor(identity(2), sigmaz(), identity(2)) +  tensor(identity(2), identity(2), sigmaz())


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
def ancillary():
	return tensor(rho_1(0), rho_1(0))

# EDIT
def error(M, dt):
	epsilon = (E(0, dt) * M * E(0, dt).dag() +
			   E(1, dt) * M * E(1, dt).dag() +
			   E(2, dt) * M * E(2, dt).dag())		
	return epsilon

# Use numpy to calculate direct sum. input qutip matrix and n3. Returns 
def rho_system(rho):
	ancillary = tensor(rho_1(0), rho_1(0))
	return tensor(rho, ancillary)

# Use U to encode a matrix
def encode(U, M):
	return(U * M * U.dag())

# Use (same) U to decode a matrix
def decode(U, M):
	return(U.dag() * M * U)

def fidelity(rho, M):
	return((rho * M).tr())

# Calculate the Gell-Mann matrices for n dimensions. Returns a lis of all matrices
def Gell_Mann(n):
	gellmann = []
	gellmann.append(Qobj(1 / (np.sqrt(n)) * identity(n), [[2, 2], [2, 2]]))
	for k in range (1, n):
		A = np.zeros((n, n), dtype=np.complex)
		for j in range (0, n):
			M = np.zeros((n, n), dtype=np.complex)
			if j < k:	
				M[k][j] = 1 / (np.sqrt(2))
				M[j][k] = 1 / (np.sqrt(2))
				gellmann.append(Qobj(M, [[2, 2], [2, 2]]))
				M[k][j] = 1j / (np.sqrt(2))
				M[j][k] = -1j / (np.sqrt(2))
				gellmann.append(Qobj(M, [[2, 2], [2, 2]]))
				A[j][j] =  np.sqrt(1 / (k * (k + 1)))
		A[k][k] =  -k * np.sqrt(1 / (k * (k + 1)))
		gellmann.append(Qobj(A, [[2, 2], [2, 2]]))
	return gellmann

def E_U(k, t):	
	return U_().dag() * E(k, t) * U_()

E_U_array = []
for k in range(0, 3):
	E_U_array.append(E_U(k, 1/1048))

anc = ancillary()


Gell_Mann_array = Gell_Mann(4)

Gell_Mann_array8 = Gell_Mann(4)
for x in range(len(Gell_Mann_array)):
	Gell_Mann_array8[x] = tensor(Gell_Mann_array8[x], identity(2))



start = time.time()

# slow original method to check outcome of other methods
totalsum = 0
for i in range(0,16):
	for j in range(0, 16):
		A = Gell_Mann_array[i] * Gell_Mann_array[j]
		if A[0][0][0] != 0:
			# print i, j
			for k in range(0, 3):
				trace_i = (tensor(identity(2), Gell_Mann(4)[i]) * E_U_array[k]).tr()
				trace_j = np.conj((tensor(identity(2), Gell_Mann(4)[j]) * E_U_array[k]).tr())
				totalsum += trace_i * trace_j * A[0][0][0]




eind = time.time()

# Mind the /4!
print totalsum 

######################################### TEST AREA ########################
print "TEST AREA"


print (E_U(0, 0.1) * E_U(0, 0.1).dag() + E_U(01, 0.1)* E_U(1, 0.1).dag() + E_U(2, 0.1)* E_U(2, 0.1).dag()).tr() / 8



############################################################################
