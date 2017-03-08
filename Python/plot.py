import numpy as np
import math
import matplotlib.pyplot as plt
import test


# Define bra and ket
# for spin-down
Bra = np.matrix([[1, 0]])
Ket = np.matrix([[1], [0]])

# Rho is |KetxBra|
Rho = np.dot(Ket, Bra)

# dt = 0.05

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

# red dashes, blue squares and green triangles
plt.plot(dt, E1, 'r--')
plt.show()	