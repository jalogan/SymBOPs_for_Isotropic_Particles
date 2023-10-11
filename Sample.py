from Particle import Particle
from importSettings import importSettings
from ql_global import ql_global
from maxQl4 import maxQl4
from getqlm import getqlm
from bonds_to_parts import bonds_to_parts
from plotBondSymBOP import plotBondSymBOP
from pathlib import Path
import numpy as np
from distPBC import distPBC
from scipy.spatial import distance_matrix
import sys


"""
The main class for the analysis. It has functions to read the data, find the neighbors, compute BOPs,
find the orientation of the highest SymBOP signal, find domains, etc. 
"""


class Sample:

	# global variables
	parts_used_in_doms = []
	bonds_used_in_doms = []
	excluded_bonds = []
	excluded_parts = []
	domain_summary = {}
	bonds_for_each_domain = {}
	# The l-value used for the SymBOP Analysis
	analysis_lval = 6


	def __init__(self, iteration, settingsIN):

		# After a domain is found, the bonds that were used 
		# to make the domain are removed to look for new 
		# domains. These are the excluded bonds.
		self.bonds = {}
		self.Ql = {}
		self.Qlmbar = {}
		self.Qlmtilde = {}
		self.phi_theta_psi = {}
		self.equatorial_pts_doms = []
		self.non_equatorial_pts_doms = []
		self.min_max_angles = {}
		self.domains = []
		self.bonds_in_domains = []

		self.readSettings(settingsIN)

		# Rename the output directories to be iter0, iter1, etc.
		out = self.OUT[:]
		self.OUT = {-1:out}
		self.OUT[0] = out+"iter"+str(iteration)+"/"


		# Make the directory if it doesn't exist
		Path(self.OUT[0]).mkdir(parents=False, exist_ok=True)

		# Make ql directory if it doesn't exist
		for lval in self.l:
			self.OUT[lval] = self.OUT[0]+"q"+str(lval)+"/"
			Path(self.OUT[lval]).mkdir(parents=False, exist_ok=True)

		# Read particle data
		self.readData()



	def readSettings(self, settingsIN):
		importSettings(self, settingsIN)


	def setSampleVariables(self, boundaries, rmin, rmax, dom_min, l, ql4_ring_width, half_angle, lattice, IN, FILE, OUT):


		self.l = l
		self.boundaries = boundaries
		self.rmin, self.rmax = rmin, rmax
		self.lattice_type = lattice
		self.dom_min = dom_min
		self.half_angle = half_angle
		self.ql4_ring_width = ql4_ring_width
		self.IN = IN
		self.FILE = FILE
		self.OUT = OUT
	

	def readData(self):

		inpath = Path(self.IN + self.FILE)

		inp = []
		with open(inpath,"r") as f:
			fileread = f.readlines()[self.startline:]
			for line in fileread:
				inp.append(line.strip("\n").split(self.delim))
			
		x = [float(i[1]) for i in inp]
		y = [float(i[2]) for i in inp]
		if self.dim==3:	
			z = [float(i[3]) for i in inp]

		try:
			num = [int(i[0]) for i in inp]
		except:
			#if len(list(set(num))) < len(num):
			num = range(len(x))
		
		COM = [sum(x)/len(x), sum(y)/len(y), sum(z)/len(z)]

		x = [i-COM[0] for i in x]
		y = [i-COM[1] for i in y]
		if self.dim==3:
			z = [i-COM[2] for i in z]
		

		if self.boundaries==[] or self.boundaries==None:
			self.boxL = 2.05*max(max(x), max(y), max(z), np.abs(min(x)), np.abs(min(y)), np.abs(min(z)))	
			self.boundaries = [[-0.5*self.boxL, 0.5*self.boxL], [-0.5*self.boxL, 0.5*self.boxL], [-0.5*self.boxL, 0.5*self.boxL]]
		else:
			self.boxL = max(self.boundaries[0][1] - self.boundaries[0][0], self.boundaries[0][1] - self.boundaries[0][0], self.boundaries[0][1] - self.boundaries[0][0])


		self.centers = {}
		if self.dim==3:
			for i in range(len(num)):
				self.centers[num[i]] = Particle(x[i], y[i], z[i], 1)

		elif self.dim==2:
			for i in range(len(num)):
				self.centers[num[i]] = Particle(x[i], y[i], 1)





	
	def getAllNeighbors(self, rot=0):

		if not self.pbc:
			# Reset bonds
			self.bonds = {}

			nppos = np.array([self.centers[i].pos for i in self.centers], dtype=float)
			dist = distance_matrix(nppos, nppos)
		
			numbonds = 0
			cent_keys = list(self.centers.keys())
			for piind in range(0, nppos.shape[0]-1):
				for pjind in range(piind+1, nppos.shape[0]):
					pi = cent_keys[piind]
					pj = cent_keys[pjind]
					if self.rmin <= dist[piind][pjind] <= self.rmax:
						self.centers[pi].neighs.append(pj)
						self.centers[pj].neighs.append(pi)
						self.bonds[numbonds] = [pi, pj]
						self.bonds[numbonds+1] = [pj, pi]
						numbonds += 2

			neighsum = sum([len(self.centers[p].neighs) for p in self.centers])
			avgnumneighs = neighsum/len(self.centers)

			if not rot:
				print("\nNumber of particles: ", len(self.centers))
				print("Number of bonds: ", len(self.bonds))
				print("Avg Neighs/particle: ", round(avgnumneighs, 2))

		else:
			center_keys = list(self.centers.keys())
			total_neighs = 0
			for i in range(len(center_keys)):
				#excluded_i_pairs = [pp[0] for pp in excluded_part_pairs if pp[1]==center_keys[i]] + [pp[1] for pp in excluded_part_pairs if pp[0]==center_keys[i]] 	
				pi = center_keys[i]
				pj = center_keys[j]
				for j in range(i+1, len(center_keys)):
					#if center_keys[j] not in excluded_i_pairs:
					dist = distPBC(self.centers[pi].pos(), self.centers[pj].pos(), self.boundaries)
					if self.rmin <= dist <= self.rmax:
						self.centers[pi].neighs.append(pj)
						self.centers[pj].neighs.append(pi)
						total_neighs += 1
					
						self.bonds[len(self.bonds)] = [pi, pj]
						self.bonds[len(self.bonds)] = [pj, pi]

			neighsum = sum([len(self.centers[p].neighs) for p in self.centers])
			avgnumneighs = neighsum/len(self.centers)

		if neighsum==0:
			print("\nNo Particles have any neighbors. Try a different rmin and rmax.\n")
			sys.exit(0)

	










	def getAllLocalql(self, lval, obj):
		for part in self.centers:
			self.centers[part].getLocalql(lval, obj)

	def getGlobalQl(self, lval, particles=None):
		if particles==None:
			particles = [part for part in self.centers if part not in Sample.excluded_parts]
		self.Ql[lval], self.Qlmbar[lval], self.Qlmtilde[lval] = ql_global(self, lval, particles) 


	def plotBondSymBOP(self):
		plotBondSymBOP(self)

	def maximizeQl4(self, lval):
		# Find the Euler angles that maximize |Ql4|
		phi, theta, psi = maxQl4(self, lval)
		self.phi_theta_psi[lval] = [phi, theta, psi]	

	def getPotentialDomainParticles(self, lval, init_phi_theta_psi):
		self.min_max_angles[lval] = getqlm(self, lval, init_phi_theta_psi[self.analysis_lval], include_all_bonds=False)
		# Erase last line printed to screen
		sys.stdout.write('\x1b[1A')
		sys.stdout.write('\x1b[2K')
		sys.stdout.write('\x1b[1A')
		sys.stdout.write('\x1b[2K')


	def getDomains(self, unrot_obj, iteration):
		numdoms_found = 0
		if iteration > 0:
			numdoms_found = len([1 for it in Sample.domain_summary for d in Sample.domain_summary[it]])

		print("\nFinding domains for the current orientation.")

		#for dom in range(len(self.min_max_angles[4])):
		doms_parts_, bonds_used_in_doms_, parts_used_in_doms_, new_domains, num_rem_parts = bonds_to_parts(self, unrot_obj)#, dom, self.min_max_angles[4])

		# Erase last line printed to screen
		sys.stdout.write('\x1b[1A')
		sys.stdout.write('\x1b[2K')


		
		for dom_num, dom_parts in enumerate(new_domains):#bonds_used_in_doms_:

			dom_num_consec = dom_num + numdoms_found

			'''
			# write bonds_used_in_current_dom to file to recreate the percolation in a visualization
			with open(self.OUT[0]+"bonds_used_in_dom_"+str(dom_num_consec)+".txt", "w") as f:
				f.write("bond\tparticle_i\tparticle_j\n\n")
				for bb in bonds_used_in_doms_[dom_num]:
					f.write(str(bb)+"\t"+str(self.bonds[bb][0])+"\t"+str(self.bonds[bb][1])+"\n")
			'''
			# write particle positions in domain to file
			with open(str(self.OUT[0])+"domain_"+str(dom_num_consec)+".txt", "w") as f:
				f.write("particle_number\tx\ty\tz\n"+str(len(new_domains[dom_num]))+"\n")
			with open(str(self.OUT[0])+"domain_"+str(dom_num_consec)+".xyz", "w") as f:
				f.write(str(len(new_domains[dom_num]))+"\n\n")

			with open(str(self.OUT[0])+"domain_"+str(dom_num_consec)+".txt", "a") as f:
				for i in new_domains[dom_num]:
					f.write(str(i)+"\t"+str(unrot_obj.centers[i].x)+"\t"+str(unrot_obj.centers[i].y)+"\t"+str(unrot_obj.centers[i].z)+"\n")
			with open(str(self.OUT[0])+"domain_"+str(dom_num_consec)+".xyz", "a") as f:
				for i in new_domains[dom_num]:
					f.write("c\t"+str(unrot_obj.centers[i].x)+"\t"+str(unrot_obj.centers[i].y)+"\t"+str(unrot_obj.centers[i].z)+"\t"+str(dom_num_consec)+"\n")




		# Collect all bonds and particles found in some domain up to this point
		Sample.bonds_used_in_doms += [b for d in bonds_used_in_doms_ for b in d]
		Sample.parts_used_in_doms += parts_used_in_doms_

		Sample.bonds_used_in_doms = list(set(Sample.bonds_used_in_doms))
		Sample.parts_used_in_doms = list(set(Sample.parts_used_in_doms))

		# Update the bonds and particles that should be excluded from future domain searches
		Sample.excluded_bonds = Sample.bonds_used_in_doms #list(set(Sample.bonds_used_in_doms))
		Sample.excluded_parts = Sample.parts_used_in_doms #list(set(Sample.parts_used_in_doms))

		self.domains = new_domains
		self.bonds_in_domains = bonds_used_in_doms_
		Sample.domain_summary[iteration] = new_domains
		Sample.bonds_for_each_domain[iteration] = bonds_used_in_doms_

		domain_sizes = [len(j) for i in Sample.domain_summary for j in Sample.domain_summary[i]]

		total_parts = num_rem_parts + len(Sample.excluded_parts)

		#print("\n#------------------------------------------#\nEnd of Iteration")
		print("\n#------------------------------------------#")
		print("Summary of Domains: ")
		for ind, num in enumerate(domain_sizes):
			print("domain {}: {}".format(ind, num))
		print("{} particles found across {} ordered domains.".format(len(Sample.parts_used_in_doms), len(domain_sizes)))
		print("{}% particles found in some domain.".format(round(100*len(Sample.parts_used_in_doms)/(total_parts),1)))
		print("End of Iteration\n#------------------------------------------#\n\n")

		if len(parts_used_in_doms_) > 0:
			return 1
		return 0






