import numpy as np
from scipy.spatial.transform import Rotation as R






def rotateData(phi_theta_psi, particle_data, PBC, L=None):

	phi, theta, psi = phi_theta_psi

	num = [i for i in particle_data]
	x = [particle_data[i][0] for i in particle_data]
	y = [particle_data[i][1] for i in particle_data]
	z = [particle_data[i][2] for i in particle_data]

	comx = sum(x)/len(x)
	comy = sum(y)/len(y)
	comz = sum(z)/len(z)


	if PBC==False:
		# Shift the COM to (0,0,0)
		particle_data = {i:[particle_data[i][0]-comx, particle_data[i][1]-comy, particle_data[i][2]-comz, particle_data[i][3]] for i in particle_data}

		x = [particle_data[i][0] for i in particle_data]
		y = [particle_data[i][1] for i in particle_data]
		z = [particle_data[i][2] for i in particle_data]

	else:
		# Repeat the packing all around the main box
		# rotate it all, then enforce the original PBC at the end
		particle_data_ext = {i:[] for i in particle_data}
		for i in particle_data:
			for shiftx in [-L,0,L]:
				for shifty in [-L,0,L]:
					for shiftz in [-L,0,L]:
						particle_data_ext[i].append([particle_data[i][0]+shiftx, particle_data[i][1]+shifty, particle_data[i][2]+shiftz, particle_data[i][3]])			
						

		x = [particle_data_ext[i][j][0] for i in particle_data_ext for j in range(len(particle_data_ext[i]))]
		y = [particle_data_ext[i][j][1] for i in particle_data_ext for j in range(len(particle_data_ext[i]))]
		z = [particle_data_ext[i][j][2] for i in particle_data_ext for j in range(len(particle_data_ext[i]))]

		comx = sum(x)/len(x)
		comy = sum(y)/len(y)
		comz = sum(z)/len(z)


		# Shift the COM to (0,0,0)
		particle_data_ext = {i:[[particle_data_ext[i][j][0]-comx, particle_data_ext[i][j][1]-comy, particle_data_ext[i][j][2]-comz, particle_data_ext[i][j][3]] for j in range(len(particle_data_ext[i]))] for i in particle_data_ext}

		x = [particle_data_ext[i][j][0] for i in particle_data_ext for j in range(len(particle_data_ext[i]))]
		y = [particle_data_ext[i][j][1] for i in particle_data_ext for j in range(len(particle_data_ext[i]))]
		z = [particle_data_ext[i][j][2] for i in particle_data_ext for j in range(len(particle_data_ext[i]))]





	rot = R.from_euler('ZYZ', [-psi,-theta,-phi], degrees=False)




	# Rotate all particles
	if PBC==False:
		particle_data = {i:list(rot.apply(particle_data[i][:3]))+[particle_data[i][3]] for i in particle_data}
	else:
		particle_data_ext = {i:[list(rot.apply(particle_data_ext[i][j][:3]))+[particle_data_ext[i][j][3]] for j in range(len(particle_data_ext[i]))] for i in particle_data_ext}


	if PBC==True:
		# pick which realization of each particle is inside the box
		particle_data = {}
		for i in particle_data_ext:
			for j in range(len(particle_data_ext[i])):
				if np.abs(particle_data_ext[i][j][0]) < L/2 and np.abs(particle_data_ext[i][j][1]) < L/2 and np.abs(particle_data_ext[i][j][2]) < L/2:
					particle_data[i] = [particle_data_ext[i][j][0]+comx, particle_data_ext[i][j][1]+comy, particle_data_ext[i][j][2]+comz, particle_data_ext[i][j][3]]
					break
	else:
		particle_data = {i:[particle_data[i][0]+comx, particle_data[i][1]+comy, particle_data[i][2]+comz, particle_data[i][3]] for i in particle_data}


	x = [particle_data[i][0] for i in particle_data]
	y = [particle_data[i][1] for i in particle_data]
	z = [particle_data[i][2] for i in particle_data]

	comx = sum(x)/len(x)
	comy = sum(y)/len(y)
	comz = sum(z)/len(z)

	return particle_data





