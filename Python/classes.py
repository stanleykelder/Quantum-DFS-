import numpy as np
from scipy import sparse
import math
# from qutip import *


def bra(spin):
	return np.matrix([[1 - spin, 0 + spin]])

def ket(spin):
	return np.matrix([[1 - spin], [0 + spin]])

def Id(size):
	return np.identity(size)

# Pauli matrices. type defines X or Z
def Pauli(type):
	if type == "x":
		return np.matrix(([[0, 1], [1, 0]]))
	if type == "z":
		return np.matrix(([[1, 0], [0, -1]]))

# Total spin. type defines X or Z
def SpinTotal(type):
	X1 = sparse.kron((sparse.kron(Pauli(type), Id(2)).toarray()), Id(2)).toarray()
	X2 = sparse.kron((sparse.kron(Id(2), Pauli(type)).toarray()), Id(2)).toarray()
	X3 = sparse.kron((sparse.kron(Id(2), Id(2)).toarray()), Pauli(type)).toarray()
	return X1 + X2 + X3

def E(k, dt):
	E2 = math.sqrt(dt) * SpinTotal("z")
	E1 = math.sqrt(dt) * SpinTotal("x")
	E0 = Id(8) - ((E1 ** 2) / 2) - ((E2 ** 2) / 2)
	
	if k == 2:
		return E2 
	if k == 1:
		return E1 
	if k == 0:
		return E0 

print(E(1, 0.1))
print(E(1, 0.1)**2)

# print(E(1, 0.1) * E(1, 0.1))

# print((E(1, 0.1) ** 2) / 2)
# print((E(2, 0.1) ** 2) / 2)

# print(E(0, 0.1))





