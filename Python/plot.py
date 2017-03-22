# import numpy as np
import matplotlib.pyplot as plt
from qutip_test import *

# Set variables 
spin = 0
n3 = 4

# Define time-range
t_0 = 0
t_eind = 0.1
dt = 0.001

t_list = np.arange(t_0, t_eind, dt)

# dt = np.arange(0, 0.1, 0.001)
F_list = []
for t in t_list:
	F_list.append(((projection(epsilon(spin, t)).ptrace(0)) * rho_1(spin)).tr())
# print(F_list)


plt.plot(t_list, F_list)
plt.show()

# Gaat epsilon goed?
# imaginair

