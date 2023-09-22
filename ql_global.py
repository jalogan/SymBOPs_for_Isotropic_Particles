import numpy as np




def ql_global(self, l, particles):

	# Keep only particles that have neighbors (this was changed 5/23/2020)
	particles = [i for i in particles if len(self.centers[i].neighs)>0]

	neigh_total = sum([len(self.centers[i].neighs) for i in particles])


	if isinstance(l, int):

		if len(particles)!=0:
	
			# average slmbar weighted by the number of neighbors
			Qlmbar = list(sum([np.array(self.centers[p].qlmbar[l], dtype=complex)*len(self.centers[p].neighs)/neigh_total for p in particles]))
			Qlmtilde = list(sum([np.array(self.centers[p].qlmtilde[l], dtype=complex)*len(self.centers[p].neighs)/neigh_total for p in particles]))
				
			Qlmbar_mag_sq = np.abs(np.vdot(np.array(Qlmbar, dtype=complex), np.array(Qlmbar, dtype=complex)))
			Ql = np.abs(np.sqrt((4*np.pi/(2*l+1))*Qlmbar_mag_sq))

		else:
			Qlmbar = [0]*(2*l+1)
			Qlmtilde = [0]*(2*l+1)
			Ql = 0.0


		return [Ql, Qlmbar, Qlmtilde]






