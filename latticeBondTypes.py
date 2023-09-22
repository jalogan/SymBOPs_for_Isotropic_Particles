import numpy as np


"""
 This file is for listing the unique bond types for
 different lattices, with respect to l-fold symmetry. 
 For instance, FCC has two unique bond types with l=4,6
 which can be represented in spherical coords 
 by (phi, theta) = (0, pi/2) and (phi, theat) = (pi/4, pi/4)
 The non-trivial m-values for the spherical harmonic vector |q)
 must be known. 
 This returns the non-trivial m-values and the phi, theta angles for each
 [[m-values], [[angles]]]
"""

def latticeBondTypes(lattice):


	if lattice=="fcc":
		return [[4], [[0.0, np.pi/2.0],[np.pi/4.0, np.pi/4.0]]]

	elif lattice=="bcc":
		return [[4], [[np.pi/4.0, np.pi/4.0]]]

	elif lattice=="sc":
		return [[4], []]


