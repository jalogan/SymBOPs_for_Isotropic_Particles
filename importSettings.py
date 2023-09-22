from pathlib import Path



def importSettings(self, settingsIN):

	"""Import the parameters for the analysis."""



	######################
	## General Settings ##
	######################

	# Set parameters
	# Changing these parameters may lead to 
	# undefined paths through the code, but you're free to try
	self.startline = 2
	self.delim = "\t"	
	self.dim = 3
	self.pbc = False




	##########################################
	## Read Settings from file Settings.txt ##
	##########################################

	OUT = None
	if settingsIN!=None and Path.is_file(Path(str(settingsIN)+"Settings.txt")):
		bound = []
		boundaries = []
		with open(str(settingsIN)+"Settings.txt", "r") as f:
			lines = f.readlines()
			for i in range(len(lines)):
				#if lines[i].strip("\n")==str(ratio):
				#	line = lines[i+1].strip("\n").split("\t")
				#	rmin = float(line[0])
				#	rmax = float(line[1])
				if lines[i].strip("\n")=="Box":
					if lines[i+1].strip("\n")!="None":
						for j in [1,2,3]:
							bound.clear()
							line = lines[i+j].strip("\n").split("\t")[0].split(",")
							bound.append(float(line[0]))
							bound.append(float(line[1]))
							boundaries.append(bound)
					else:
						boundaries = None
				elif (lines[i].strip("\n")).lower()=="lattice_type":
					line = lines[i+1].strip("\n").split("\t")
					lattice_type = str(line[0]).lower()

				elif (lines[i].strip("\n")).lower()=="l":
					line = lines[i+1].strip("\n").split(",")
					l = [int(i) for i in line]
				elif (lines[i].strip("\n")).lower()=="dom_min":
					line = lines[i+1].strip("\n")
					dom_min = int(line)
				#elif (lines[i].strip("\n")).lower()=="ideal_qs":
				#	line = lines[i+1].strip("\n").split(",")
				#	line_ = []
				#	ideal_qs = {}
				#	for li in line:
				#		line_ = li.split(":")
				#		#print(line_)
				#		ideal_qs[int(line_[0])] = [float(line_[1])]
				elif (lines[i].strip("\n")).lower()=="half_angle":
					line = lines[i+1].strip("\n")
					half_angle = float(line)
				elif (lines[i].strip("\n")).lower()=="ql4_ring_width":
					line = lines[i+1].strip("\n")
					ql4_ring_width = float(line)
				elif (lines[i].strip("\n")).upper()=="DATAIN":
					line = lines[i+1].strip("\n")
					IN = str(line)
				elif (lines[i].strip("\n")).upper()=="DATAOUT":
					line = lines[i+1].strip("\n")
					OUT = str(line)
				elif (lines[i].strip("\n")).upper()=="DATAFILE":
					line = lines[i+1].strip("\n")
					FILE = str(line)
				elif (lines[i].strip("\n")).lower()=="rmin_rmax":
					line = lines[i+1].strip("\n").split(",")
					rmin,rmax = [float(i) for i in line]


		if OUT == None:
			OUT = IN

	##########################################
	##########################################


	## make output directory if it doesn't exist ##
	Path(OUT).mkdir(parents=True, exist_ok=True)	


	# Set class variables for Sample
	#self.setSampleVariables(boundaries, rmin, rmax, dom_min, l, ideal_qs, ql4_ring_width, half_angle, lattice_type, IN, FILE, OUT)
	self.setSampleVariables(boundaries, rmin, rmax, dom_min, l, ql4_ring_width, half_angle, lattice_type, IN, FILE, OUT)







