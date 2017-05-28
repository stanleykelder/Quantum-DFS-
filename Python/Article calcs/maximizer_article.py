from __future__ import division
import numpy as np
# import random
import time
from defenitions_article import *
from qutip import Qobj, identity, tensor


start = time.time()

# List of 8D Gell-Mann matrix
Gell_Mann_article8 = Gell_Mann(4)
for x in range(len(Gell_Mann_article8)):
	Gell_Mann_article8[x] = tensor(identity(2), Gell_Mann_article8[x]).full()

E_np_list = []
for k in range(0, 4):
	E_np_list.append(E_np(k, 1/1024))


stepsize = [1/(2**3), -1/(2**3)]
best = 0


phi_best = [0, 3.0030797283097694, 2.3002983311738294, 1.8107355659962843, 2.234768861986662, 0.8880802795405287, 0.6990774476882087, 0.9505866717064115, 5.210862891000058, 5.71410517126595, 6.275546251168045, 1.817757585740357, 2.498050396865478, 1.219779927476257, 0.6324427734683677, 3.6499051384967265, 2.61058349986131, 1.6538444641528631, 2.1549168546357893, 2.0336885961032363, 3.4212577668974524, 1.8585442471827343, 1.114565894312593, 2.0696460574710467, 0.3034236121126904, 1.2436097298475253, 0.5993017151499396, 2.731927013920969, 2.9451076440538535, 0.8240422546194267, 2.0214952443634413, 0.06914867705319128, 5.06348726947566, 1.4392449384452337, 5.401776905422341, 6.146098566025535, 4.821940846772457]
Cos_best = [0, -0.91087513741979531, 0.92543502587619619, -0.99641897531666757, 0.034011679452355705, 0.83425450305902804, -0.64484188116668495, -0.077211419965364583, 0.70783391731522127, 0.93868191469102702, -0.96890788184776599, 0.79615105290837529, 0.8641569763404835, -0.9999179829255279, 0.59180701302514238, 0.38396459958663937, -0.98119132735087355, 0.8385211872808132, -0.8833399549809868, -0.9994827247946706, -0.5440955044613552, 0.19044121720729648, -0.16890301915595332, -0.14663255743231604, -0.99132872147910178, -0.9384362546056485, 0.30753606841442876, -0.79096733286634735, 0.56272125802419792]
Sin_best = [0, 0.41268206167759341, 0.37890633787457828, 0.084553093550042632, 0.99942143546195272, 0.55137956447961878, 0.76431600028614544, 0.99701474243209265, 0.7063788965549499, 0.34478437179197585, 0.24742173811788551, 0.60509792674647744, 0.50322233678769979, -0.012807319084950477, 0.80607968547424569, 0.9233478143496479, 0.19303776607552964, 0.5448689920349435, 0.4687328918842673, -0.0321602679873909, 0.8390232905140019, 0.98169860078773841, 0.98563267504684704, 0.98919102963020167, 0.13140534985535923, 0.3454524512023652, 0.95153642422357942, 0.6118583809496253, 0.82664671158098979]


# loop several times to find maximum p0
for y in range(1000):
	
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

			# calculate p0
			p_0 = np.real(p0(E_np_list, U, U_dag, Gell_Mann_article8))
			
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




