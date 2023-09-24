import numpy as np







def bonds_to_parts(self, unrot_obj, dom_num, min_max_angles):



	# Get all of the bonds found from l=4 and l=6
	bs = [int(i[0]) for i in self.equatorial_pts_doms]	
	bs += [int(i[0]) for i in self.non_equatorial_pts_doms]	

	#print("# bs: ", len(bs))
	with open(str(self.OUT[0])+"bonds_used_in_doms.txt", "w") as f:
		pass


	# Make partbonds, dict of all bonds each particle is part of
	partbonds = {i:[] for i in self.centers}

	for i in bs:
		parta = self.bonds[i][0]
		partb = self.bonds[i][1]
		partbonds[parta].append(i)
		partbonds[partb].append(i)


	# Rotate the data back to its original orientation to compare with full data set
	# This is now fixed in setup.py by recording the original positions as orig_data
	bonds_used = []
	domcount = -1
	bonds_used_in_doms = []
	parts_used_in_doms = []
	domain_summary = []
	while len(bonds_used) < len(bs):
		b = np.random.choice([i for i in bs if i not in bonds_used])
		parta = self.bonds[b][0]
		partb = self.bonds[b][1]
		dom = [parta, partb]
		parts_used = []
		bonds_used.append(b)
		bonds_used_in_current_dom = [b]

		# Find all bonds that branch off of all particles found to be in domain
		while len(parts_used) < len(dom):
			for p in [i for i in dom if i not in parts_used]:
				parts_used.append(p)
				for b in partbonds[p]:
					if b in bs and b not in bonds_used:
						partb = [i for i in self.bonds[b] if i!=p]
						bonds_used.append(b)
						if partb[0] not in dom:
							dom.append(partb[0])	
							bonds_used_in_current_dom.append(b)


		# This allows for there to be one empty bond between domains and to call them the same domain
		# it looks at second nearest neighbors
		'''
		while len(parts_used) < len(dom):
			for p in [i for i in dom if i not in parts_used]:
				parts_used.append(p)
				for b in partbonds[p]:#+[bds for bb in [partbonds[pp] for pp in Particle.data[p].neighs] for bds in bb]:
					#if b in bs:
					if b in bs and b not in bonds_used:
						partb = [i for i in Particle.bonds[b] if i!=p]
						bonds_used.append(b)
						if partb[0] not in dom:
							dom.append(partb[0])	
							bonds_used_in_current_dom.append(b)
				# Look at second neighbors
				for pp in Particle.data[p].neighs:
					#parts_used.append(pp)
					for b in partbonds[pp]:
						if b in bs and b not in bonds_used:
							partb = [i for i in Particle.bonds[b] if i!=pp]
							bonds_used.append(b)
							if partb[0] not in dom:
								dom.append(partb[0])	
								bonds_used_in_current_dom.append(b)
		'''






		with open(str(self.OUT[0])+"bonds_used_in_doms.txt", "a") as f:
			for bb in bonds_used_in_current_dom:
				f.write(str(bb)+"\t"+str(self.bonds[bb][0])+"\t"+str(self.bonds[bb][1])+"\n")


		if len(dom) > self.dom_min-1:
			domcount += 1

			bonds_used_in_doms.append(bonds_used_in_current_dom)
			parts_used_in_doms += dom
			domain_summary.append(dom)

			'''
			# write ql values for each particle in a domain to a common file
			with open(str(self.OUT[0])+"ql_values_for_all_particles_in_some_domain.txt", "a") as f: 
				for i in dom:
					f.write(str(i)+"\t")
					for lval in range(len(self.l)):
						if lval != len(self.l)-1:
							f.write(str(self.centers[i].ql[self.l[lval]])+"\t")
						else:
							f.write(str(self.centers[i].ql[self.l[lval]])+"\n")
			'''



	# write the particels that are NOT in domains to a file to resubmit as a "new" sample
	remaining_data = [k for k in self.centers if k not in parts_used_in_doms]
	with open(str(self.OUT[0])+"data_remaining.txt", "w") as f:
		f.write(str(len(remaining_data))+"\n\n")

	with open(str(self.OUT[0])+"data_remaining.txt", "a") as f:
		for i in remaining_data:
			f.write(str(i)+"\t"+str(unrot_obj.centers[i].x)+"\t"+str(unrot_obj.centers[i].y)+"\t"+str(unrot_obj.centers[i].z)+"\n")

	with open(str(self.OUT[0])+"data_remaining.xyz", "w") as f:
		f.write(str(len(remaining_data))+"\n\n")

	with open(str(self.OUT[0])+"data_remaining.xyz", "a") as f:
		for i in remaining_data:
			f.write("c\t"+str(unrot_obj.centers[i].x)+"\t"+str(unrot_obj.centers[i].y)+"\t"+str(unrot_obj.centers[i].z)+"\n")




	return [dom, bonds_used_in_doms, parts_used_in_doms, domain_summary]



