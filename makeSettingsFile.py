

"""
This module makes a new Settings.txt file to run the next iteration of the analysis.
The new file is stored in the directory .../iter<iteration number>/ for the current iteration.
"""


def makeSettingsFile(dataIN, dataFILE, dataOUT, l, half_angle, ql4_ring_width, boundaries, rmin, rmax, dom_min, lattice_type):



	with open(str(dataIN)+"Settings.txt", "w") as f:
		f.write("dataIN\n"+str(dataIN)+"\n")	
		f.write("dataFILE\n"+str(dataFILE)+"\n")
		f.write("dataOUT\n"+str(dataOUT)+"\n")

		f.write("l\n")
		for lval in l[:-1]:
			f.write(str(lval)+",")
		f.write(str(str(l[-1])+"\n"))

		f.write("half_angle\n"+str(half_angle)+"\n")
		f.write("ql4_ring_width\n"+str(ql4_ring_width)+"\n")

		f.write("Box\n")
		if boundaries == [] or boundaries==None:
			f.write("None\n")
		else:
			for b in boundaries:
				f.write(str(b[0])+","+str(b[1])+"\n")

		f.write("rmin_rmax\n"+str(rmin)+","+str(rmax)+"\n")
		f.write("dom_min\n"+str(dom_min)+"\n")
		f.write("lattice_type\n"+str(lattice_type)+"\n")

