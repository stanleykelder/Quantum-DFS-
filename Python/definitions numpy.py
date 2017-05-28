from __future__ import division
import numpy as np
import time
from Uencode import *
from numpy.core.umath_tests import inner1d
from qutip import Qobj, identity, tensor


start = time.time()

# Tensor product np: np.tensordot(A, B, 0) 0 is for specifying np to use tensor

sigmaz = np.array([[ 1,  0], [ 0, -1]])
sigmax = np.array([[ 0,  1], [ 1, 0]])

spintotalz = np.array([[ 3.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
				       [ 0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.],
					   [ 0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.],
					   [ 0.,  0.,  0., -1.,  0.,  0.,  0.,  0.],
					   [ 0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.], 
					   [ 0.,  0.,  0.,  0.,  0., -1.,  0.,  0.], 
					   [ 0.,  0.,  0.,  0.,  0.,  0., -1.,  0.], 
					   [ 0.,  0.,  0.,  0.,  0.,  0.,  0., -3.]])

spintotalx = np.array([[ 0.,  1.,  1.,  0.,  1.,  0.,  0.,  0.],
					   [ 1.,  0.,  0.,  1.,  0.,  1.,  0.,  0.],
					   [ 1.,  0.,  0.,  1.,  0.,  0.,  1.,  0.],
					   [ 0.,  1.,  1.,  0.,  0.,  0.,  0.,  1.],
					   [ 1.,  0.,  0.,  0.,  0.,  1.,  1.,  0.],
					   [ 0.,  1.,  0.,  0.,  1.,  0.,  0.,  1.],
					   [ 0.,  0.,  1.,  0.,  1.,  0.,  0.,  1.],
					   [ 0.,  0.,  0.,  1.,  0.,  1.,  1.,  0.]])



# Takes as input a value for k (0, 1 or 2) and a time t. Returns corresponding Error
def E(k, t):
	E2 = np.sqrt(t) * spintotalz
	E1 = np.sqrt(t) * spintotalx
	E0 = np.identity(8) - (np.dot(E1, E1) / 2) - (np.dot(E2, E2) / 2)
	
	if k == 2:
		return E2 
	if k == 1:
		return E1 
	if k == 0:
		return E0 

# 4D ancillary 
def ancillary():
	zeros = np.zeros((4, 4))
	zeros[0, 0] = 1.
	return zeros

# Calculate the Gell-Mann matrices for n dimensions. Returns a lis of all matrices
def Gell_Mann_np(n):
	gellmann = []
	gellmann.append((1 / (np.sqrt(n)) * np.identity(n)))
	for k in range (1, n):
		A = np.zeros((n, n), dtype=np.complex)
		for j in range (0, n):
			M = np.zeros((n, n), dtype=np.complex)
			if j < k:	
				M[k, j] = 1 / (np.sqrt(2))
				M[j, k] = 1 / (np.sqrt(2))
				gellmann.append(M)
				M[k, j] = 1j / (np.sqrt(2))
				M[j, k] = -1j / (np.sqrt(2))
				gellmann.append(M)
				A[j, j] =  np.sqrt(1 / (k * (k + 1)))
		A[k][k] =  -k * np.sqrt(1 / (k * (k + 1)))
		gellmann.append(A)
	return gellmann

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


# Create lists of Gell-Mann matrices (all in numpy)
# List of 4D Gell-Mann matrix
Gell_Mann_array4 = Gell_Mann_np(4)
# List of 8D Gell-Mann matrix
Gell_Mann_array8 = Gell_Mann(4)
for x in range(len(Gell_Mann_array8)):
	Gell_Mann_array8[x] = tensor(Gell_Mann_array8[x], identity(2)).full()



def E_U(k, t):	
	return U * E(k, t) * U_dag

A = np.dot(E_U(0, 0.1), Gell_Mann_array8[4])
for i in range(0, 8):
	print A[i, i]
print type(E_U(0, .1)), type(Gell_Mann_array8[4].T)

# print np.sum(inner1d(E_U(0, 0.1), Gell_Mann_array8[4].T))
print np.einsum('ij,ji->', E_U(0, 0.1), Gell_Mann_array8[4])



####################################################################################
#############          Calculations for p_0         ################################


p0 = 0

# list for all j=j'!=0
list_even = [0, 1, 2, 3, 4, 5, 8, 9, 10, 15]
list_j = [[0, 3], [0, 8], [0, 15], [1, 2], [3, 8], [3, 15], [4, 5], [8, 15], [9, 10]]

for k in range(0, 3):
	# start with j=j'=0 (sigmaj * sigmaj' [0,0] = 0.25)
	# E_k real and E_k' real are equal and right part is 0.25 + 0j -> tr(E_k)^2 * 8
	print "a"
	# For 


# start adding for j=j'
for i in list_even:
	print i, (Gell_Mann_array[i] * Gell_Mann_array[i])[0,0], (Gell_Mann_array[i] * Gell_Mann_array[i]).tr()