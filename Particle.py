from ql_local import ql_local






class Particle:


	def __init__(self, position, cent_type):
		self.pos = position
		self.x = position[0]
		self.y = position[1]
		if position[2]:
			self.z = position[2]
		else:
			self.z = 0.0
		self.center_type = cent_type
		self.neighs = []
		self.ql = {}
		self.qlmbar = {}
		self.qlmtilde = {}
		self.neighs_qlm = {}
		self.bond_symbops = {}

	def __init__(self, x, y, z=0.0, cent_type=1):
		self.pos = [x,y,z]
		self.x = x
		self.y = y
		self.z = z
		self.center_type = cent_type
		self.neighs = []
		self.ql = {}
		self.qlmbar = {}
		self.qlmtilde = {}
		self.neighs_qlm = {}
		self.bond_symbops = {}



	def getLocalql(self, lval, sample_obj):
		self.ql[lval], self.qlmbar[lval], self.qlmtilde[lval], self.neighs_qlm[lval] = ql_local(self, lval, sample_obj)








