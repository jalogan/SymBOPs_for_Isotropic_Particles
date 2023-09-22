from Sample import Sample
from rotateData import rotateData
from makeSettingsFile import makeSettingsFile
from plotBondSymBOP import plotBondSymBOP
from makeVisFile import makeVisFile
from pathlib import Path


"""
This is the brains of the analysis. From here all of the functions are called.
"""


def main(settingsIN):


	# Make instance of Sample
	init_sam = Sample(iteration=0, settingsIN=settingsIN)
	objects = {0:[init_sam]}


	print("\nSettings Used:\n")
	print("boundaries: ", init_sam.boundaries)
	print("rmin, rmax: ", init_sam.rmin, init_sam.rmax)
	print("dom_min: ", init_sam.dom_min)
	print("l: ", init_sam.l)
	print("ql4_ring_width: ", init_sam.ql4_ring_width)
	print("half_angle: ", init_sam.half_angle)
	print("lattice: ", init_sam.lattice_type)
	print("dataIN: ", init_sam.IN)
	print("dataOUT: ", init_sam.OUT)
	print("dataFILE: ", init_sam.FILE)
	print("#------------------------------------------#\n\n")

	print("iteration: 0")



	# Make directories to store output files
	iteration = -1
	new_doms_found = 1
	while new_doms_found and len(Sample.excluded_parts) < len(objects[0][0].centers):
		iteration += 1
		
		if iteration > 0:
			# Make new Settings.txt file to initialize new instance of Sample
			makeSettingsFile(objects[0][0].OUT[0], "data_remaining.txt", objects[0][0].OUT[-1], objects[0][0].l, objects[0][0].half_angle, objects[0][0].ql4_ring_width, objects[0][0].boundaries, objects[0][0].rmin, objects[0][0].rmax, objects[0][0].dom_min, objects[0][0].lattice_type)

			print("iteration: ", iteration)

			# Make instance of Sample
			init_sam = Sample(iteration, objects[0][0].OUT[0])
			objects[iteration] = [init_sam]



		## Get all attributes of the remaining particles in the sample ##

		# Get neighbors for all particles
		init_sam.getAllNeighbors()

		# Compute local ql order parameters for each particle
		# This includes computing the single-bond ql for each
		# neighbor as well as the standard BOP ql for the particle
		for lval in init_sam.l:
			init_sam.getAllLocalql(lval, init_sam)

		# Compute global Ql order parameter for whole sample or region
		for lval in init_sam.l:
			init_sam.getGlobalQl(lval)

		# Find Euler angles that align reference tensor to highest signal domains
		for lval in init_sam.l:
			init_sam.maximizeQl4(lval)	





		# Rotate data to align highest signal domains with reference tensor	
		data_to_rot = {p:[init_sam.centers[p].x, init_sam.centers[p].y, init_sam.centers[p].z, init_sam.centers[p].center_type] for p in init_sam.centers}
		rotated_data = rotateData(init_sam.phi_theta_psi[init_sam.analysis_lval], data_to_rot, init_sam.pbc, init_sam.boxL)

		
		# Write rotated data to file to be read back in as new analysis
		with open(str(init_sam.OUT[0])+"rotated_data_phi="+str(round(init_sam.phi_theta_psi[init_sam.analysis_lval][0],2))+"_theta="+str(round(init_sam.phi_theta_psi[init_sam.analysis_lval][1],2))+"_psi="+str(round(init_sam.phi_theta_psi[init_sam.analysis_lval][2],2))+".txt", "w") as f:
			f.write("\n\n")
			for i in rotated_data:
				f.write(str(i)+"\t"+str(rotated_data[i][0])+"\t"+str(rotated_data[i][1])+"\t"+str(rotated_data[i][2])+"\t"+str(rotated_data[i][3])+"\n")
		

		dataIN = init_sam.OUT[0]#+"iter"+str(iteration)+"/"
		dataFILE = "rotated_data_phi="+str(round(init_sam.phi_theta_psi[init_sam.analysis_lval][0],2))+"_theta="+str(round(init_sam.phi_theta_psi[init_sam.analysis_lval][1],2))+"_psi="+str(round(init_sam.phi_theta_psi[init_sam.analysis_lval][2],2))+".txt"
		

		# Make new Settings.txt file to initialize new instance of Sample
		makeSettingsFile(dataIN, dataFILE, init_sam.OUT[-1], init_sam.l, init_sam.half_angle, init_sam.ql4_ring_width, init_sam.boundaries, init_sam.rmin, init_sam.rmax, init_sam.dom_min, init_sam.lattice_type)




		#########################
		# Now working with rotated data


		# Make instance of Sample
		rot_sam = Sample(iteration, dataIN)
		objects[iteration].append(rot_sam)


		## Get all attributes of the rotated remaining particles in the sample ##

		# Get neighbors for all particles
		rot_sam.getAllNeighbors()

		# Compute local ql order parameters for each particle
		# This includes computing the single-bond ql for each
		# neighbor as well as the standard BOP ql for the particle
		for lval in rot_sam.l:
			rot_sam.getAllLocalql(lval, rot_sam)

		# Compute global Ql order parameter for whole sample or region
		for lval in rot_sam.l:
			rot_sam.getGlobalQl(lval)

		# Compute 
		for lval in rot_sam.l:
			rot_sam.getPotentialDomainParticles(lval, init_sam.phi_theta_psi)


		# Get domains
		new_doms_found = rot_sam.getDomains(objects[iteration][0], iteration)


		# Plot Bond SymBOPs for this orientation
		rot_sam.plotBondSymBOP()




	# Make file for visualizing domains in OVITO
	makeVisFile(objects)





