import numpy as np
import math


# Convention: phi,theta, and psi correspond to euler angles alpha, beta, gamma in the active zyz rotation convention. The rotation matrices rotate a vector V to V' as		
# V'= Rz(phi)Ry(theta)Rz(psi)V, with psi acting first. 





# Rotation matrices for rotations around the z- and y- axes.
#def Rz(a):
#	return np.asarray([[np.cos(a),-np.sin(a),0],[np.sin(a),np.cos(a),0],[0,0,1]])

#def Ry(a):
#	return np.asarray([[np.cos(a),0,np.sin(a)],[0,1,0],[-np.sin(a),0,np.cos(a)]])





# Wigner D matrix elements
# the matrix has m' as rows and m as columns -- Wigner_D(l,phi,theta,psi)[m+l][mp+l] should be equal to sympy's Rotation.D(l,m,mp,phi,theta,psi).doit().evalf().  


def Wigner_D(l,phi,theta,psi):

	wigner_D = np.zeros((2*l+1,2*l+1),dtype=complex)
	check = np.zeros((2*l+1,2*l+1))

	a=-1
	for mp in range(-l,l+1):
		a+=1
		b=-1
		for m in range(-l,l+1):
			b+=1
			I = 1j
			wigner_D[a][b] = np.exp(-I*(mp*phi+m*psi))*Wigner_d(l,mp,m,theta)

	wigner_D_active = np.conj(wigner_D.T)

	return wigner_D_active




# Wigner d matrix elements
# Wigner_d(l,m,mp,theta) should be equal to sympy's Rotation.d(l,m,mp,theta).doit().evalf(). 
# m and mp can be flipped in this function iff wigner_d is multiplied by ((-1)**(m-mp)) at the end. 
# In this case the call to Wigner_d above in the Wigner_D function should have m and mp flipped also.

def Wigner_d(l,mp,m,theta):

	n_min = max(0,mp-m)
	n_max = min(l-m,l+mp)

	I = 1j

	if n_min<n_max:
		wigner_d = 0
		for n in np.arange(n_min,n_max+1):

			w = np.sqrt(float(math.factorial(l+m)*math.factorial(l-m)*math.factorial(l+mp)*math.factorial(l-mp)))/(math.factorial(l-m-n)*math.factorial(l+mp-n)*math.factorial(n+m-mp)*math.factorial(n))
			Wn = w*(np.cos(theta/2)**(2*l+mp-m-2*n))*(np.sin(theta/2)**(m-mp+2*n))		
			wigner_d += ((-1)**n)*Wn

	elif n_min==n_max:

		n = n_min

		w = np.sqrt(float(math.factorial(l+m)*math.factorial(l-m)*math.factorial(l+mp)*math.factorial(l-mp)))/(math.factorial(l-m-n)*math.factorial(l+mp-n)*math.factorial(n+m-mp)*math.factorial(n))

		Wn = w*(np.cos(theta/2)**(2*l+mp-m-2*n))*(np.sin(theta/2)**(m-mp+2*n))

		wigner_d = ((-1)**(n))*Wn

	return wigner_d





