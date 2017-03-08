import numpy as np
import math

# Define bra and ket
# for spin-down
Bra = np.matrix([[1, 0]])
Ket = np.matrix([[1], [0]])

# Rho is |KetxBra|
Rho = np.dot(Ket, Bra)

# initiate dt at 1
dt = 0.05

# Define Sx and Sz
Sx = np.matrix([[0, 1], [1, 0]])
Sz = np.matrix([[1, 0], [0, -1]])

# Identity matrix
Id = np.identity(2)

# Define E_k as function of dt
E1 = math.sqrt(dt) * Sx
E2 = math.sqrt(dt) * Sz
E0 = Id - ((E1 ** 2) / 2) - ((E2 ** 2) / 2)  

# Define conjugate transpose E_k
E1_T = (math.sqrt(dt) * Sx).getH()
E2_T = math.sqrt(dt) * Sz.getH()
E0_T = Id - ((E1 ** 2) / 2) - ((E2 ** 2) / 2).getH()


ERho = np.sum([(E0 * Rho * E0_T), (E1 * Rho * E1_T), (E2 * Rho * E2_T)], 0) 

# print(Rho * E1_T)
# print(E1)
# print(E2)
# print(E0)
# print(((E1 ** 2) / 2))
# # print((E0 * Rho * E0_T))
# # print((E1 * Rho * E1_T))
# # print((E2 * Rho * E2_T))
# print(ERho)



for x in range(0, 10):
	# Define E_k as function of dt
	E1 = math.sqrt(dt) * Sx
	E2 = math.sqrt(dt) * Sz
	E0 = Id - ((E1 ** 2) / 2) - ((E2 ** 2) / 2)  

	# Define conjugate transpose E_k
	E1_T = (math.sqrt(dt) * Sx).getH()
	E2_T = math.sqrt(dt) * Sz.getH()
	E0_T = Id - ((E1 ** 2) / 2) - ((E2 ** 2) / 2).getH()
	

	ERho = np.sum([(E0 * Rho * E0_T), (E1 * Rho * E1_T), (E2 * Rho * E2_T)], 0) 

	F = Rho * ERho
	print(dt)
	print(F)
	dt = dt + 0.05
