import numpy as np




def vecPBC(a, b, boundaries):

	delta = [a[0] - b[0], a[1] - b[1], a[2] - b[2]]

	Lx = boundaries[0][1] - boundaries[0][0]
	Ly = boundaries[1][1] - boundaries[1][0]
	Lz = boundaries[2][1] - boundaries[2][0]


	# X
	if np.abs(delta[0]) > 0.5*Lx:
	
		if np.abs(delta[0] - Lx) > 0.5*Lx:			
			x = delta[0] + Lx
		
		else:			
			x = delta[0] - Lx					
	else:
		x = delta[0]

	# Y
	if np.abs(delta[1]) > 0.5*Ly:
	
		if np.abs(delta[1] - Ly) > 0.5*Ly:			
			y = delta[1] + Ly

		else:
			y = delta[1] - Ly
	else:
		y = delta[1]
	
	# Z
	if Lz!=None:
		if np.abs(delta[2]) > 0.5*Lz:
		
			if np.abs(delta[2] - Lz) > 0.5*Lz:			
				z = delta[2] + Lz

			else:
				z = delta[2] - Lz
		else:
			z = delta[2]

	return [x,y,z]







