import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from copy import deepcopy
from matplotlib import pyplot
from pathlib import Path
from tqdm import tqdm
from computeBondTypeRadii import computeBondTypeRadii

plt.rc('text', usetex=True)





def getqlm(self, lval, init_phi_theta_psi, include_all_bonds=False):

	print("\nCollecting bonds for potential ordered domains using l="+str(lval)+".")

	phi, theta, psi = init_phi_theta_psi
	
	# Get magnitude of ideal qlm for each unique bond
	# This sets the radii of the rings for ideal bonds for the lattice type in the qlm complex plane plot
	qlm_bond_radii = computeBondTypeRadii(lval, self.lattice_type)


	used_bonds = []

	chi_l = {4:np.sqrt(7/12), 6:-np.sqrt(2)/4}
	ql4 = []

	for i in tqdm([part for part in self.centers if part not in self.excluded_parts], leave=False):

		# Get all ql4 values 
		# Write file of all ql4 components for each bond
		for j in self.centers[i].neighs:
			bond_num = [b for b in self.bonds if (self.bonds[b][0]==i and self.bonds[b][1]==j) or (self.bonds[b][0]==j and self.bonds[b][1]==i)]
			if len(bond_num)>0:
				if bond_num[0] not in used_bonds and bond_num[0] not in self.excluded_bonds:
					used_bonds += bond_num
					if lval==4:
						qlm_val = self.centers[i].neighs_qlm[lval][j][-1]
						ang = np.arctan2(np.imag(qlm_val), np.real(qlm_val))
						if ang < 0.0:
							ang += 2*np.pi


						ql4.append([bond_num[0], qlm_val, ang])
						if len(bond_num)>1:
							for bb in range(1,len(bond_num)):
								ql4.append([bond_num[bb], qlm_val, ang])

						# Get full symbop
						# Cubic
						ql4_val = self.centers[i].neighs_qlm[lval][j][-1]
						ql0_val = self.centers[i].neighs_qlm[lval][j][lval]
						ql4_mag = np.absolute(ql4_val)
						symbop_val = np.absolute(ql4_mag*np.sqrt(2.0*(1.0 - chi_l[lval]*chi_l[lval])) + chi_l[lval]*ql0_val)
					

					elif lval==6:
						qlm_val = self.centers[i].neighs_qlm[lval][j][-3]
						ang = np.arctan2(np.imag(qlm_val), np.real(qlm_val))

						ql4.append([bond_num[0], qlm_val, ang])
						if len(bond_num)>1:
							for bb in range(1,len(bond_num)):
								ql4.append([bond_num[bb], qlm_val, ang])

						# Get full symbop
						# Cubic
						ql4_val = self.centers[i].neighs_qlm[lval][j][-3]
						ql0_val = self.centers[i].neighs_qlm[lval][j][lval]
						ql4_mag = np.absolute(ql4_val)
						symbop_val = np.absolute(ql4_mag*np.sqrt(2.0*(1.0 - chi_l[lval]*chi_l[lval])) + chi_l[lval]*ql0_val)



	##########################################################
	## Method using the fact that we can find the orientation of all

	if lval==4:
		#min_max_angles = [[np.pi-self.half_angle, np.pi+self.half_angle]]
		min_max_angles = [np.pi-self.half_angle, np.pi+self.half_angle]
	elif lval==6:
		#min_max_angles = [[-self.half_angle, self.half_angle]]
		min_max_angles = [-self.half_angle, self.half_angle]
			

	min_angle = min_max_angles[0]#min_max_angles[0][0]
	max_angle = min_max_angles[1]#min_max_angles[0][1]



	ql4_count = 0
	ql4_sum = 0
	if lval==4:
		equatorial_pts = []
		for i in ql4:

			# Take points that you believe belong in a domain
			# constrain the radial position in the ql4-complex plane
			if np.abs(i[1]) >= max(qlm_bond_radii)-self.ql4_ring_width:

				equatorial_pts.append([i[0],i[1],i[2]])

				ql4_count += 1
				ql4_sum += i[1]


		self.equatorial_pts_doms = []
		pts_to_plot_dom = []
		bonds_to_plot_dom = []
		for i in equatorial_pts:
			if min_angle <= i[2] <= max_angle:
				self.equatorial_pts_doms.append([i[0], i[1], i[2]])
				pts_to_plot_dom.append(i[1])
				bonds_to_plot_dom.append(i[0])
			else:
				pass


	elif lval==6:

		maxql4 = max([np.abs(i[1]) for i in ql4])	

		self.non_equatorial_pts_doms = []
		pts_to_plot_dom = []
		bonds_to_plot_dom = []
		for i in ql4:
			if np.abs(i[1]) >= maxql4 - self.ql4_ring_width and (min_angle <= i[2] <= max_angle):
				self.non_equatorial_pts_doms.append([i[0], i[1], i[2]])
				pts_to_plot_dom.append(i[1])
				bonds_to_plot_dom.append(i[0])				
				
				ql4_count += 1
				ql4_sum += i[1]




	##########################################################
	# PLOTS

	#aspect = 160/256
	fig, ax = plt.subplots(figsize=(15,15))



	
	###############
	# Plot all bonds	
	ax.scatter([np.real(i) for i in [k[1] for k in ql4]], [np.imag(i) for i in [k[1] for k in ql4]], s=3, c="silver")
	###############



	###############
	# Plot bonds that are in the potential domains
	if lval==4:
		pts_real = [np.real(i) for i in pts_to_plot_dom]
		pts_imag = [np.imag(i) for i in pts_to_plot_dom]
		col = [[0.004,0,0.5725,1] for cc in range(len(pts_real))]
		ax.scatter(pts_real, pts_imag, s=3, c=col)# c="r")
	elif lval==6:
		pts_real = [np.real(i) for i in pts_to_plot_dom]
		pts_imag = [np.imag(i) for i in pts_to_plot_dom]
		col = [[0.6,0.6,0.184,1] for cc in range(len(pts_real))]
		ax.scatter(pts_real, pts_imag, s=3, c=col)# c="r")
	###############



	## Plot the circles representing the ideal FCC values
	for br in qlm_bond_radii:
		angles = np.linspace(0.0, 2*np.pi, 3000)
		xvals_r1 = [br*np.cos(i) for i in angles]
		yvals_r1 = [br*np.sin(i) for i in angles]

		ax.scatter(xvals_r1, yvals_r1, s=0.1, c="k", zorder=3)




	##########
	# Draw vectors from the origin marking the region of interest
	if lval==4 and min_angle!=None and max_angle!=None:
		max_qlm_rad = max(qlm_bond_radii)

		xs_min = [0.0, max_qlm_rad*np.cos(min_angle)]
		ys_min = [0.0, max_qlm_rad*np.sin(min_angle)]

		xs_max = [0.0, max_qlm_rad*np.cos(max_angle)]
		ys_max = [0.0, max_qlm_rad*np.sin(max_angle)]

		ax.plot(xs_min, ys_min, "k-")
		ax.plot(xs_max, ys_max, "k-")

	elif lval==6 and min_angle!=None and max_angle!=None:

		xs_min = [0.0, maxql4*np.cos(min_angle)]
		ys_min = [0.0, maxql4*np.sin(min_angle)]

		xs_max = [0.0, maxql4*np.cos(max_angle)]
		ys_max = [0.0, maxql4*np.sin(max_angle)]

		ax.plot(xs_min, ys_min, "k-")
		ax.plot(xs_max, ys_max, "k-")
	#########



	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	ax.set_xlabel(r"$\mathrm{Re}(q_{"+str(lval)+r"}^4)$", fontsize=90)
	ax.set_ylabel(r"$\mathrm{Im}(q_{"+str(lval)+r"}^4)$", fontsize=90)








	if lval==4: 
		ax.set_xticks(np.arange(-0.5,0.75,0.25))
		ax.set_yticks(np.arange(-0.5,0.75,0.25))
		ax.set_xticklabels(["-0.5", "-0.25", "0.0", "0.25", "0.5"])
		ax.set_yticklabels(["-0.5", "-0.25", "0.0", "0.25", "0.5"])

	elif lval==6:
		ax.set_xticks(np.arange(-0.4, 0.6, 0.2))
		ax.set_yticks(np.arange(-0.4, 0.6, 0.2))
		ax.set_xticklabels(["-0.4", "-0.2", "0.0", "0.2", "0.4"])
		ax.set_yticklabels(["-0.4", "-0.2", "0.0", "0.2", "0.4"])

		

	ax.axis("equal")

	ax.tick_params(labelsize=80, width=5, length=15)

	plt.tight_layout()


	plt.savefig(Path(str(self.OUT[lval])+"q"+str(lval)+"4_scatter_in_complex_plane_phi="+str(round(phi,3))+"_theta="+str(round(theta,3))+"_psi="+str(round(psi,3))+".pdf"), dpi=300, transparent=True)
	plt.savefig(Path(str(self.OUT[lval])+"q"+str(lval)+"4_scatter_in_complex_plane_phi="+str(round(phi,3))+"_theta="+str(round(theta,3))+"_psi="+str(round(psi,3))+".png"), dpi=300, transparent=True)



	#print("\nShowing the good plot...\n")
	#plt.show()
	plt.close()



	
	# Write domain points to file
	min_angle = min_max_angles[0]#min_max_angles[0][0]
	max_angle = min_max_angles[1]#min_max_angles[0][1]





	return min_max_angles












