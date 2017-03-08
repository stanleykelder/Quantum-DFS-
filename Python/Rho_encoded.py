import numpy as np
from scipy import sparse

# Define bra and ket
# for spin-down 
Bra = np.matrix([[1, 0]])
Ket = np.matrix([[1], [0]])

# Rho is |KetxBra|
Rho = np.dot(Ket, Bra)
Id = np.identity(2)

O = np.zeros((4, 4))

Kronprod = sparse.kron(Rho, Id).toarray()

print(Kronprod)

Encoded = np.zeros( np.add(Kronprod.shape,O.shape) )

Encoded[:Kronprod.shape[0],:Kronprod.shape[1]]=Kronprod
print(Encoded)

trace = Encoded.trace()
print(trace)