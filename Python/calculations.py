 from definitions import *
import numpy as np
import time
from qutip import *
from Uencode import *

# Just to determine running time
start = time.time()


# Variables to be defined by user
spin = 0


# Calculations. Uncomment prints to see partial results

# Rho. logic qbit to be encoded. 2x2 matrix
# rho_logic = rho_1(spin)
rho_logic = Qobj([[.5,.5], [.5,.5]])
# print(rho_logic)

# Rho within system of 3 qbits. 8x8 matrix
rho_system = rho_system(rho_logic)
# print(rho_system)

# Rho encoded. Rho (of 3 qbit system) encoded by encoding matrix U. 8x8 matrix
rho_encoded = encode(U(), rho_system)
# print(rho_encoded)

# Error on rho_encoded
error = error(rho_encoded, 0.01)
# print(error)

# Rho decoded. Rho_system after error. 8x8 matrix
rho_decoded = decode(U(), error)
# print(rho_decoded)

# Rho after time t, after encoding and decoding. 2x2 matrix 
rho_t = rho_decoded.ptrace(0)
# print(rho_t)

fidelity = fidelity(rho_logic, rho_t)



print("fid:", fidelity)





# Just to determine running time
eind = time.time()