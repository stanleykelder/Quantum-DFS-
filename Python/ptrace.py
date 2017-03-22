import numpy as np
import time
from qutip import *

# def bra(spin):
# 	return basis(2, spin).dag()

# def ket(spin):
# 	return basis(2, spin)

def identity(size):
	return(qeye(size))

# ket0 = tensor(ket("0"), ket("0"))
# bra0 = tensor(bra("0"), bra("0"))

# ket1 = tensor(ket("1"), ket(1))
# bra1 = tensor(bra(1), bra(1))

rho = (ket("00") + ket("11")) * (bra("00") + bra("11"))
PartialTrace = rho.ptrace(1)

O = create(4)


# test = to_super((ket(0) * bra(0)))

print(O * rho)








# print(PartialTrace)
