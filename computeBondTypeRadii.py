from scipy.special import sph_harm
import numpy as np
from latticeBondTypes import latticeBondTypes




def computeBondTypeRadii(lval, lattice):

	
	non_trivial_mvals, unique_bond_angles = latticeBondTypes(lattice)
	#print("mvals: ", non_trivial_mvals)	
	#print("bond angles: ", unique_bond_angles)

	radii = []
	for mval in non_trivial_mvals:
		for phi,theta in unique_bond_angles:
			radii.append(np.absolute(np.sqrt(4.0*np.pi/(2.0*lval + 1.0))*sph_harm(mval, lval, phi, theta)))

	#print("radii: ", radii)


	return radii

