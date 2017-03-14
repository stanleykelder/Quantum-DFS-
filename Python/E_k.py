import numpy as np
import math
from sympy.physics.quantum import *

Ket(0) = np.matrx([[1, 0]])
bra = Bra(0)
print(Ket(0))
print(bra)

psi = ket * bra
print(psi)


# Define bra and ket
# for spin-down
# Bra = np.matrix([[1, 0]])
# Ket = np.matrix([[1], [0]])

# Rho is |KetxBra|
# Rho = np.dot(Ket, Bra)

