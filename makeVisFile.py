from numpy.random import uniform





def getColor():#num, total):
	#scale = num/(total-1)
	r, g, b = uniform(0.0, 1.0, 3)
	return [r, g, b]
	#return [scale*1.0, 1.0-scale, (1.0-scale)*scale]





def makeVisFile(objs):

	colors = [[1,0,0],[0,1,0],[0,0,1],[1,0,127/255],[204/255,0,204/255],[153/255,153/255,0],[0,204/255,204/255],[0, 153/255, 0],[1,153/255,153/255],[139/255,69/255,19/255]]

	first_obj = objs[list(objs.keys())[0]][0]
	last_obj = objs[list(objs.keys())[-1]][1]

	parts_in_doms = []
	for it in last_obj.domain_summary:
		parts_in_doms += [part for d in last_obj.domain_summary[it] for part in d]

	parts_not_in_doms = set(first_obj.centers) - set(parts_in_doms)

	total_num_doms = len([1 for it in last_obj.domain_summary for d in last_obj.domain_summary[it]])




	# write OVITO extended XYZ file to visualize domains along with the rest of the sample

	radius = 0.25*first_obj.rmin
	domnum = -1
	partcount = 0
	dompartcount = 0
	with open(str(first_obj.OUT[0])+"visualize_sample.xyz", "w") as f:
		f.write(str(len(first_obj.centers))+"\n")
		f.write("Lattice=\"1 0 0 0 1 0 0 0 1\" pbc=\"F F F\" Properties=\"type:I:1:pos:R:3:color:R:3:radius:R:1:\"\n")
		
		for orient in objs:
			obj = objs[orient][1]
			for domind, dom in enumerate(obj.domains):
				domnum += 1
				if domnum < len(colors):
					col = colors[domnum]
				else:
					col = getColor()#domnum, total_num_doms)
				for part in dom:				
					f.write(str(domnum)+"\t"+str(first_obj.centers[part].x)+"\t"+str(first_obj.centers[part].y)+"\t"+str(first_obj.centers[part].z)+"\t"+str(col[0])+"\t"+str(col[1])+"\t"+str(col[2])+"\t"+str(radius)+"\n")
					dompartcount += 1


		# include the particles that are not found in any domains
		domnum += 1
		col = [0.4,0.4,0.4]
		for part in parts_not_in_doms:
			f.write(str(domnum)+"\t"+str(first_obj.centers[part].x)+"\t"+str(first_obj.centers[part].y)+"\t"+str(first_obj.centers[part].z)+"\t"+str(col[0])+"\t"+str(col[1])+"\t"+str(col[2])+"\t"+str(radius)+"\n")
			partcount += 1







	# write OVITO extended XYZ file to visualize domains only

	with open(str(first_obj.OUT[0])+"visualize_domains_only.xyz", "w") as f:
		f.write(str(dompartcount)+"\n")
		f.write("Lattice=\"1 0 0 0 1 0 0 0 1\" pbc=\"F F F\" Properties=\"type:I:1:pos:R:3:color:R:3:radius:R:1:\"\n")
	
		domnum = -1
		partcount = 0
		dompartcount = 0
		
		for orient in objs:
			obj = objs[orient][1]
			for domind, dom in enumerate(obj.domains):
				domnum += 1
				if domnum < len(colors):
					col = colors[domnum]
				else:
					col = getColor()#domnum, total_num_doms)
				for part in dom:				
					f.write(str(domnum)+"\t"+str(first_obj.centers[part].x)+"\t"+str(first_obj.centers[part].y)+"\t"+str(first_obj.centers[part].z)+"\t"+str(col[0])+"\t"+str(col[1])+"\t"+str(col[2])+"\t"+str(radius)+"\n")
					dompartcount += 1




