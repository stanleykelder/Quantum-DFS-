from __future__ import division
import numpy as np
# import random
import time
from definitions_numpy_second import *
from qutip import Qobj, identity, tensor


start = time.time()

# List of 8D Gell-Mann matrix
Gell_Mann_np_array8 = Gell_Mann(4)
for x in range(len(Gell_Mann_np_array8)):
	Gell_Mann_np_array8[x] = tensor(identity(2), Gell_Mann_np_array8[x]).full()

# List for errors
E_np_list = []
for k in range(0, 5):
	E_np_list.append(E_np(k, 1/1024))


stepsize = [1, -1]
best = 0


phi_best = [0, 2.9952672283097694, 2.980003648994243, 2.3097590034962843, 2.392971986986662, 0.9330021545405287, 4.779137754867795, 1.4203132342064115, 5.208909766000058, 5.71410517126595, 5.976718126168045, 2.204476335740357, 2.883792584365478, 5.883824609655843, 5.790628080647954, 3.1499051384967265, 3.12327881236131, 1.2983757141528631, 0.8189793546357893, 4.322751096103236, 0.007213084717866636, 1.5812004971827343, 1.503237769312593, 2.0637866824710467, 0.2995173621126904, 1.2553284798475253, 0.4655126526499396, 1.618645763920969, 2.7937404565538535, 0.8250188171194267, 0.8154405568634413, 4.460732421732778, 5.41797945697566, 1.1150261884452337, 5.469159717922341, 0.20080388384594894, 5.872722096772457]
Cos_best = [0, -0.93417739336217809, 0.9254350258761962, -0.9995902730489987, -0.14225535582500634, 0.62738115517081106, -0.99999937971134834, -0.024551501640582639, 0.9768416434287011, 0.97585965390252938, -0.78446479632891253, 0.81353434887623932, 0.84121553902515755, -0.5510349792568418, 0.58232034721926973, -0.48323960759099138, -0.9989052543909569, 0.65729188585910936, -0.50655345897338677, -0.9865818284185698, -0.79621540882644393, -0.26231427548759551, 0.0011818457451394787, -0.36144706447741359, -0.9993133594940041, -0.99984229644703437, 0.4906171265974994, -0.89684455496887883, -0.32040223339557206]
Sin_best = [0, 0.35680890926523451, -0.3789063378745782, -0.02862317288191865, 0.98982999234166513, 0.77871232566111259, 0.0011138118865750019, 0.99969856645250443, -0.21396355685843887, 0.21839857116115841,  0.6201733494117897, 0.58151686407060754, 0.54069993240485392, -0.8344821457858832, 0.81295941670813376, 0.87548813907174361, 0.046779191422445063, 0.75363610368915823, 0.86220855551432418, -0.16326755903813733, 0.60501324179503602, 0.96498249770419042, 0.99999930162007344, 0.9323926316636465, -0.037051444463169624, 0.017759004350488664, 0.87137525503609137, 0.44234584232551982, 0.94728158898773562]



# loop several times to find maximum p0
for y in range(500):
	
	prev_best = best
	phi = list(phi_best)	
	Cos = list(Cos_best)
	Sin = list(Sin_best)


	# stepsize eiher posiive or negative
	for x in range(2):
		step = stepsize[x]	

		# loop through all variables to find biggest increase
		for n in range(1, len(phi) + len(Cos)):
			phi_check = list(phi)
			Cos_check = list(Cos)
			Sin_check = list(Sin)

			# first through phases
			if n < len(phi):
				
				phi_check[n] = phi_check[n] + step
				if phi_check[n] < 0:
					phi_check[n] = phi_check[n] + 2 * np.pi
				elif phi_check[n] > 2 * np.pi:
					phi_check[n] = phi_check[n] - 2 * np.pi
			
			# then through angles
			else:
				i = n - len(phi)
				if i == 0:
					continue
				new = np.arccos(Cos_check[i]) + step
				Cos_check[i] = np.cos(new)
				Sin_check[i] = np.sin(new)

			# Make U
			U = np.array(U_(Cos_check, Sin_check, phi_check))
			U_dag = np.conj(U.transpose())

			# Create list E'_k
			E_U_k = []
			for k in range(0, 5):
				E_U_k.append(np.dot(U_dag, np.dot(E_np_list[k], U)))


			# calculate p0
			p_0 = np.real(p0(E_U_k, Gell_Mann_np_array8) / 4)
			
			# find best
			if p_0 > best:
				best = p_0
				phi_best = list(phi_check)
				Cos_best = list(Cos_check)
				Sin_best = list(Sin_check)

	print best 

	if prev_best == best:
		if stepsize[0] < 1/(2**20):
			break
		else:
			print stepsize[0]
			stepsize[0] = stepsize[0] / 2
			stepsize[1] = stepsize[1] / 2



end = time.time()

print best
print phi_best
print Cos_best
print Sin_best


print end - start




