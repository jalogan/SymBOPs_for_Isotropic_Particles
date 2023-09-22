import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc,colors
rc('text', usetex=True)
from wigner_D import Wigner_D
import copy
from pathlib import Path



def Rz(a):
	return np.asarray([[np.cos(a),-np.sin(a),0],[np.sin(a),np.cos(a),0],[0,0,1]])

def Ry(a):
	return np.asarray([[np.cos(a),0,np.sin(a)],[0,1,0],[-np.sin(a),0,np.cos(a)]])





def maxQl4(self, l):

	idealcubic2 = [-0.2082449933578368-0.35492033774429954j, 0, -0.8329799734313472, 0, -0.1041224966789184+0.35492033774429954j]
	idealcubic4 = [-0.4564354645876385, 0, 0, 0, -0.7637626158259734, 0, 0, 0, -0.4564354645876387]
	idealcubic6 = [0, 0, 0.6614378277661475, 0, 0, 0, -0.3535533905932746, 0, 0, 0, 0.6614378277661475, 0, 0]

	idealcubic = {2:idealcubic2, 4:idealcubic4, 6:idealcubic6}

	hot_spot_signals = []
	'''
	with open(str(self.OUT[l])+"hot_spot_signals.txt","w") as bestf:
			bestf.write('signal_strength,\tphi,\ttheta,\tpsi\tin radians\n')
	'''

	#with open(str(self.OUT[l])+"hot_spots.txt","w") as bestf:
	#		bestf.write('Pl value,\tphi,\ttheta,\tpsi\tin radians\n')



	I = 1j


	start_phi = 0.0
	end_phi = 2*np.pi
	step_phi = 0.02
	start_theta = 0.0
	end_theta = np.pi
	step_theta = 0.02

	phi_range = np.arange(start_phi,end_phi+step_phi,step_phi)
	theta_range = np.arange(start_theta,end_theta+step_theta,step_theta)

	z = np.zeros((len(theta_range),len(phi_range)))

	best_list = []
	best_list.append([self.Qlmtilde[l][-1], 0.0, 0.0, 0.0])

	maxQl4_mag = 0.0
	best = 0.0 + 0.0j
	Ql4 = 0.0 + 0.0j

	psi = 0.0
	j=-1
	for phi in phi_range:
		j+=1
		i=-1
		for theta in theta_range:
			i+=1
			local_max = -10.0
				
			D = Wigner_D(l,phi,theta,psi)
			first = np.dot(D,np.conj(self.Qlmtilde[l]).T)

			sec = np.dot(idealcubic[l], first)

			psi0_first = first
			psi0_sec = sec
			psi0_Ql4 = psi0_first[3-l]
			psi0_Ql4_phase = np.angle(psi0_Ql4)
			psi0_Ql4_mag = np.absolute(psi0_Ql4)
			if psi0_Ql4_phase < 0.0:
				psi0_Ql4_phase += 2*np.pi

				
			# Because for real data the domain will not have perfect symmetry,
			# check all of the rotations by np.pi/2 to see which gives the highest signal
			for psi in [-0.25*psi0_Ql4_phase + nn*0.5*np.pi for nn in range(-4, 5)]:

				if l==6:
					sig = 0.5*np.sqrt(7)*psi0_Ql4_mag*np.cos(4*psi + psi0_Ql4_phase) - psi0_first[6]/np.sqrt(8)
				elif l==4:
					sig = 0.5*np.sqrt(7)*psi0_Ql4_mag*np.cos(4*psi + psi0_Ql4_phase) - psi0_first[4]/np.sqrt(8)
			
				sec_val = np.absolute(sig)
			
				if sec_val > local_max:
					psi_local_max = psi
					local_max = sec_val
					local_sec_max = sec
					local_first_max = first
					if l==6:
						Ql4_mag = np.absolute(first[-3])
						Ql0 = first[6]
					elif l==4:
						Ql4_mag = np.absolute(first[-1])
						Ql0 = first[4]


			# Check if this phi, theta, psi local_max is actually a global max
			if local_max > maxQl4_mag:
				maxQl4_mag = local_max
				best_phi = phi
				best_theta = theta
				best_psi = psi_local_max
				bestQlm = local_first_max


			# Record all values for all angles
			# and for the best psi
			z[i][j] = local_max
		
			# Record hot spot
			hot_spot_signals.append([local_max, phi, theta, psi_local_max])

			
	z_vals = [j for i in z for j in i]
	max_z_val = max(z_vals)






	################
	# Plot heatmap


	aspect = 8/6
	fig = plt.figure(figsize=(20,20/aspect))
	ax = fig.add_subplot(1,1,1)

	heatmap = ax.pcolor(phi_range,theta_range,z,vmin=0,vmax=max_z_val,cmap="rainbow",rasterized=True)

	cbar = fig.colorbar(heatmap)

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	ax.set_xlabel(r"$\phi$",fontsize=65)
	ax.set_ylabel(r"$\cos\theta$",fontsize=65)

	xticks = np.arange(start_phi,end_phi+step_phi,np.pi/3)
	yticks = np.arange(start_theta,end_theta+step_theta,np.pi/3)

	ax.set_xticks(xticks)
	ax.set_yticks(yticks)

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	ax.set_xticklabels([r"$0$",r"$\frac{\pi}{3}$",r"$\frac{2\pi}{3}$",r"$\pi$",r"$\frac{4\pi}{3}$",r"$\frac{5\pi}{3}$",r"$2\pi$"])#,rotation=0,fontsize=2)

	ax.set_yticklabels([str(round(np.cos(i),2)) for i in yticks])

	plt.tick_params(axis='both', which='major', labelsize=65)#, fontweight='bold')
	plt.tick_params(axis='x', pad=10)#, fontweight='bold')
	for m in cbar.ax.yaxis.get_ticklabels():
		m.set_size(65)


	plt.tight_layout()


	#plt.savefig(str(self.OUT[l])+"q"+str(l)+"/Max_Q"+str(l)+"4_heatmap_theta_0_2pi_vs_phi_0_2pi_phi="+str(round(best_phi,2))+"_theta="+str(round(best_theta,2))+"_psi="+str(round(best_psi,2))+"_step_"+str(step_phi)+".png",dpi=300)
	plt.savefig(Path(self.OUT[l]+"Max_Q"+str(l)+"4_heatmap_theta_0_2pi_vs_phi_0_2pi_phi="+str(round(best_phi,2))+"_theta="+str(round(best_theta,2))+"_psi="+str(round(best_psi,2))+"_step_"+str(step_phi)+".png"),dpi=300)


	#plt.savefig(str(self.OUT[l])+"q"+str(l)+"/Max_Q"+str(l)+"4_heatmap_theta_0_2pi_vs_phi_0_2pi_phi="+str(round(best_phi,2))+"_theta="+str(round(best_theta,2))+"_psi="+str(round(best_psi,2))+"_step_"+str(step_phi)+".pdf",dpi=300)
	plt.savefig(Path(self.OUT[l]+"Max_Q"+str(l)+"4_heatmap_theta_0_2pi_vs_phi_0_2pi_phi="+str(round(best_phi,2))+"_theta="+str(round(best_theta,2))+"_psi="+str(round(best_psi,2))+"_step_"+str(step_phi)+".pdf"),dpi=300)


	#plt.show()
	plt.close()




	



	# Write all hot spots to file
	# Sort hot_spots by signal strength from largest to smallest
	hot_spot_signals.sort(key = lambda x: x[0], reverse=True)

	'''
	#with open(str(self.OUT[l])+"q"+str(l)+"/hot_spot_signals.txt","w") as bestf:
	with open(str(self.OUT[l])+"hot_spot_signals.txt","w") as bestf:
		#bestf.write('Pl value,\tphi,\ttheta,\tpsi\tin radians\n')
		for hot in hot_spot_signals:
			bestf.write(str(hot[0])+"\t"+str(hot[1])+"\t"+str(hot[2])+"\t"+str(hot[3])+"\n")
	'''


	# Write data to file
	# make signals list of all data to use below

	norm = colors.Normalize(vmin=0.0, vmax=max_z_val)

	signals = []
	#with open(self.OUT+"q"+str(l)+"/sphere_colormap_phi="+str(round(best_phi,2))+"_theta="+str(round(best_theta,2))+"_psi="+str(round(best_psi,2))+"_step_"+str(step_phi)+"_data.txt", "w") as f:
	'''
	with open(self.OUT[l]+"sphere_colormap_phi="+str(round(best_phi,2))+"_theta="+str(round(best_theta,2))+"_psi="+str(round(best_psi,2))+"_step_"+str(step_phi)+"_data.txt", "w") as f:


		f.write("phi_index, theta_index, signal, phi, theta\n")
		f.write("use spherical_map.py in .../SymBOPs/extra/ to plot.\n")
		vals = []
		pts = []
		for ph in range(len(phi_range)):
			for th in range(len(theta_range)):
				#for ps in range(len(psi_range)):

				theta = th*step_theta
				phi = ph*step_phi

				val = z[th, ph]
		
				vals.append(norm(val))
				signals.append([phi, theta, val])

				f.write(str(ph)+"\t"+str(th)+"\t"+str(val)+"\t"+str(phi)+"\t"+str(theta)+"\n")			
	'''

	#vals = []
	pts = []
	for ph in range(len(phi_range)):
		for th in range(len(theta_range)):
			#for ps in range(len(psi_range)):
			theta = th*step_theta
			phi = ph*step_phi
			val = z[th, ph]
			#vals.append(norm(val))
			signals.append([phi, theta, val])



	# Find best coordinates based on hot spots
	# Find the orthogonal coordiantes
	# Take the strongest signal as the first axis

	orig_best_psi = copy.deepcopy(best_psi)

	# Sort signals largest to smallest
	signals.sort(key = lambda x: x[2], reverse=True)
	#print("\n\nsignals: ", signals, "\n\n")

	best_axes = [[signals[0][0], signals[0][1]]]
	v = best_axes[0]
	best_axes_vecs = [[np.sin(v[1])*np.cos(v[0]), np.sin(v[1])*np.sin(v[0]), np.cos(v[1])]]


	####################
	# Original method

	perp_thresh = 0.01 #0.02
	while len(best_axes)<self.dim:
		perp_thresh += 0.01

		if signals[0][1] > np.pi:
			signals[0][1] -= 2*np.pi
		if signals[0][1] < 0.0:
			signals[0][1]*=-1
			signals[0][0] += np.pi
		if signals[0][0] < 0.0:
			signals[0][0] += 2*np.pi
		elif signals[0][0] > 2*np.pi:
			signals[0][0] -= 2*np.pi


		best_axes = [[signals[0][0], signals[0][1]]]
		v = best_axes[0]
		best_axes_vecs = [np.array([np.sin(v[1])*np.cos(v[0]), np.sin(v[1])*np.sin(v[0]), np.cos(v[1])])]


		for s in signals:

			phi_val = s[0]
			theta_val = s[1]

			if theta_val > np.pi:
				theta_val -= 2*np.pi
			if theta_val < 0.0:
				theta_val*=-1
				phi_val += np.pi
			if phi_val < 0.0:
				phi_val += 2*np.pi
			elif phi_val > 2*np.pi:
				phi_val -= 2*np.pi


			vec = np.array([np.sin(theta_val)*np.cos(phi_val), np.sin(theta_val)*np.sin(phi_val), np.cos(theta_val)])

			# Check if this is close to being orthogonal to the currently found best_axes (Looking for orthogonal coordinates)
			keep = True
			for v in best_axes_vecs:
				#axis_vec = [np.sin(v[1])*np.cos(v[0]), np.sin(v[1])*np.sin(v[0]), np.cos(v[1])]
				cosangle = np.dot(v, vec)
				#print("cos(angle): ", np.abs(cosangle))
				if np.abs(cosangle) > perp_thresh:#0.03:
					keep = False
			if keep==True:
				#print("phi, theta: ", phi_val, theta_val)
				#print("vec: ", vec)				
	
				# Use Gram-Schmidt process to make sure the vector is orthogonal to others
				new_vec = np.array(copy.deepcopy(vec))
				
				for axis in best_axes_vecs:
					proj = np.array(np.dot(np.array(new_vec), np.array(axis))*np.array(axis))
					new_vec -= proj
				
				#print("new_vec: ", new_vec)

				# Make sure we have right-handed coords
				if len(best_axes)==2:
					cross = best_axes_vecs[1][0]*new_vec[1] - best_axes_vecs[1][1]*new_vec[0]
					if cross < 0:
						new_vec *= -1

				#print("new_vec: ", new_vec)

				# Make sure the vectors remain unit
				mag = np.sqrt(sum([c**2 for c in new_vec]))
				new_vec /= mag

				#print("new_vec unit: ", new_vec)

				# get new angles
				if new_vec[2]!=0.0:		
					theta_val = np.arctan(np.sqrt(new_vec[0]**2 + new_vec[1]**2)/(new_vec[2]))
				else:
					theta_val = np.pi/2
				phi_val = np.arctan2(new_vec[1],new_vec[0])


				if theta_val > np.pi:
					theta_val -= 2*np.pi
				if theta_val < 0.0:
					theta_val*=-1
					phi_val += np.pi
				if phi_val < 0.0:
					phi_val += 2*np.pi
				elif phi_val > 2*np.pi:
					phi_val -= 2*np.pi

				vec_check = np.array([np.sin(theta_val)*np.cos(phi_val), np.sin(theta_val)*np.sin(phi_val), np.cos(theta_val)])

				for goodvec in best_axes_vecs:
					inner = np.dot(vec_check, goodvec)
					#print("vec_check: ", vec_check)
					#print("goodvec: ", goodvec)
					#print("inner: ", inner)
					if np.abs(inner) > perp_thresh:
						keep = False



				if keep==True and sum([i**2 for i in new_vec]) > 0.0:				
					best_axes.append([phi_val, theta_val])
					best_axes_vecs.append(vec_check)#new_vec)
					

				#if sum([i**2 for i in vec]) > 0.0:				
				#	best_axes.append([phi_val, theta_val])
				#	best_axes_vecs.append(vec)
				#print("best_axes_vecs: ", best_axes_vecs)

			if len(best_axes)==self.dim:
				break
			#print("\n")

	###print("best_axes [phi, theta]: ", best_axes, len(best_axes))
	###print("best_axes_vecs: ", best_axes_vecs)


	#with open(self.OUT+"q"+str(l)+"/best_axes_phi_theta.txt", "w") as g:
	with open(self.OUT[l]+"best_axes_phi_theta.txt", "w") as g:
		g.write("best axes spherical phi, theta for x,y,z\nfirst vector becomes z-axis, other two x, and y (right-handed coords)\nLast lines are the unit vectors representing these axes\n\n")
		for ba in best_axes:
			g.write(str(ba[0])+"\t"+str(ba[1])+"\n")
		for ba in best_axes_vecs:
			g.write(str(ba[0])+"\t"+str(ba[1])+"\t"+str(ba[2])+"\n")



	###################
	## Original method to find Euler angles

	# Rotate best axis to global z-axis
	rotation_angles_zyz = [-best_axes[0][0], -best_axes[0][1]]
	best_axes_vecs_RyRz = []
	best_axes_RyRz = []
	for v in best_axes_vecs:
		####R = np.dot(Ry(-best_axes[0][1]),Rz(-best_axes[0][0]))
		####unit = np.dot(R, v)

		v = np.dot(Rz(-best_axes[0][0]), v)
		unit = np.dot(Ry(-best_axes[0][1]), v)

		# New vecs
		best_axes_vecs_RyRz.append(unit)

		#print("(maxQl4.py) unit: ", unit)

		if unit[2]!=0.0:		
			theta = np.arctan(np.sqrt(unit[0]**2 + unit[1]**2)/(unit[2]))
		else:
			theta = np.pi/2
		phi = np.arctan2(unit[1],unit[0])

		#########
		if theta > np.pi:
			theta -= 2*np.pi
		if theta < 0.0:
			theta*=-1
			phi += np.pi
		if phi < 0.0:
			phi += 2*np.pi
		elif phi > 2*np.pi:
			phi -= 2*np.pi
		#########


		# New angles
		best_axes_RyRz.append([phi, theta])
		#print("new angles: ", best_axes_RyRz[-1])


	


	# Rotate to align other vecs to x- and y-axes
	cross_prod = np.cross(best_axes_vecs_RyRz[1], best_axes_vecs_RyRz[2])	
	#print("(should only have z-comp) cross_prod: ", cross_prod)
	#print("phi_rot to bring the X- and Y- axes to the global X-axis: ", best_axes_RyRz[1][0], best_axes_RyRz[2][0])

	if cross_prod[2] > 0:
		phi_rot = best_axes_RyRz[1][0] 
	else:
		phi_rot = best_axes_RyRz[2][0]

	rotation_angles_zyz.append(-phi_rot)

	best_axes_vecs_RzRyRz = []
	best_axes_RzRyRz = []
	for v in best_axes_vecs_RyRz:
		R = Rz(-phi_rot)
		unit = np.dot(R, v)

		# New vecs
		best_axes_vecs_RzRyRz.append(unit)
		
		if unit[2]!=0.0:		
			theta = np.arctan(np.sqrt(unit[0]**2 + unit[1]**2)/(unit[2]))
		else:
			theta = np.pi/2
		phi = np.arctan2(unit[1],unit[0])


		#########
		if theta > np.pi:
			theta -= 2*np.pi
		if theta < 0.0:
			theta*=-1
			phi += np.pi
		if phi < 0.0:
			phi += 2*np.pi
		elif phi > 2*np.pi:
			phi -= 2*np.pi
		#########



		# New angles
		best_axes_RzRyRz.append([phi, theta])



	#print("(maxQl4.py) best_axes_RzRyRz: ", best_axes_RzRyRz)
	#print("(maxQl4.py) best_axes_vecs_RzRyRz: ", best_axes_vecs_RzRyRz)
	#print("(maxQl4.py) rotation_angles_zyz: ", rotation_angles_zyz)
	#print("psi rotation: ", phi_rot)
	#print("orig_best_psi: ", orig_best_psi)

	# Update best euler angles
	# The extra np.pi/4 on psi depends on the ideal cubic system
	best_phi = best_axes[0][0]
	best_theta = best_axes[0][1]
	best_psi = phi_rot


	################################################################



	best_list.append([maxQl4_mag, best_phi, best_theta, best_psi])

	best_list.sort(reverse=True)

	print("l={}: best_phi: {}, best_theta: {}, best_psi: {} ".format(l, best_phi, best_theta, best_psi))



	return [best_phi, best_theta, best_psi]














