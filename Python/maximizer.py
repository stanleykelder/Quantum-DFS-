from __future__ import division
import numpy as np
# import random
import time
from definitions_numpy import *
from qutip import Qobj, identity, tensor


start = time.time()

# List of 8D Gell-Mann matrix
Gell_Mann_np_array8 = Gell_Mann(4)
for x in range(len(Gell_Mann_np_array8)):
	Gell_Mann_np_array8[x] = tensor(identity(2), Gell_Mann_np_array8[x]).full()

# List for errors
E_np_list = []
for k in range(0, 6):
	E_np_list.append(E_np(k, 1/1024))


stepsize = [1, -1]
best = 0


phi_best = [0, 1.9952672283097694, 2.792503648994243, 2.2472590034962843, 2.330471986986662, 1.6205021545405287, 4.841637754867795, 1.5453132342064115, 4.896409766000058, 1.8059198640863636, 6.101718126168045, 1.079476335740357, 3.8837925843654784, 5.758824609655843, 5.540628080647954, 2.8999051384967265, 1.6232788123613102, 1.2983757141528631, 0.8814793546357893, 4.260251096103236, 6.165398391897453, 1.6437004971827343, 1.503237769312593, 3.0012866824710467, 0.5495173621126904, 2.0678284798475253, 5.998697959829526, 1.618645763920969, 6.13942576373344, 0.3250188171194267, 1.2529405568634413, 0.1150471145531915, 4.54297945697566, 1.3650261884452337, 5.344159717922341, 0.7633038838459489, 5.372722096772457]
Cos_best = [0, -0.9341773933621781, 0.9254350258761962, -0.9995902730489987, -0.080153497663315615, 0.6273811551708111, -0.53936473033970311, -0.02455150164058264, 0.34774566750987923, 0.89148986931161012, -0.88631808435502479, 0.90767244764772925, 0.8412155390251576, -0.99716571734965187, 0.67913238862222913, -0.4832396075909914, -0.99694364640839672, 0.70308012742287196, -0.85790728337690059, -0.94408164592650923, -0.66948392779656374, -0.015419069368624073, 0.0011818457451394787, -0.30250473521656784, -0.97741380125142263, -0.92385624211677397, 0.37815069091268583, -0.46849251196761754, 0.41126909997268174]
Sin_best = [0, 0.3568089092652345, -0.3789063378745782, -0.02862317288191865, 0.99678253235715208, 0.7787123256611126, 0.84207225798358853, 0.9996985664525044, 0.93758890283967655, 0.45304063053413679, 0.46307694106405162, 0.41967931540782455, 0.5406999324048539, 0.075236508043331082, 0.73401580277557077, 0.8754881390717436, -0.078124041664072247, 0.7111106344465945, 0.51380452813192146, -0.32971176173241573, 0.74282654127466752, 0.99988111908356658, 0.9999993016200734, 0.95314788211040691, -0.21133447689206894, 0.3827396555125594, 0.92574405478094124, 0.88346746755626027, 0.91151397543189672]



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
			for k in range(len(E_np_list)):
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




