import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc, colors, colorbar, cm, tri
rc('text', usetex=True)
from pathlib import Path




def plotBondSymBOP(self):


	idealcubic = {4:[-0.4564354645876385, 0, 0, 0, -0.7637626158259734, 0, 0, 0, -0.4564354645876387], 6:[0, 0, 0.6614378277661475, 0, 0, 0, -0.3535533905932746, 0, 0, 0, 0.6614378277661475, 0, 0]}




	bond_symbops = {lval:[] for lval in self.l}
	for bs in self.bonds:
		i = self.bonds[bs][0]
		j = self.bonds[bs][1]

		for lval in self.l:
			qlm_val = self.centers[i].neighs_qlm[lval][j]
			symbop_val = np.vdot(idealcubic[lval], qlm_val)

			bond_symbops[lval].append([bs, symbop_val])





	# PLOTS

	aspect = 4.8/6.4 

	for l in self.l: 

		# Plot bond symbops for each lval
		fig, ax = plt.subplots(figsize=(12,12*aspect))

		symbopdf = pd.DataFrame({"bond_num":[i[0] for i in bond_symbops[l]], "bond_symbop":[np.real(i[1]) for i in bond_symbops[l]], "color":[1.0]*len(bond_symbops[l])})

		hist = sns.histplot(symbopdf["bond_symbop"], ax=ax, kde=False, bins=25, stat="density")
		
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')

		plt.xlabel(r"$\left| b^{*}_{ij} \right|_{"+str(l)+r"}$", fontsize=22)
		plt.ylabel(r"$\mathcal{P}\left( \left| b^{*}_{ij} \right|_{"+str(l)+r"} \right)$", fontsize=22)

		ax.tick_params(labelsize=18)

		ax.set_xticks(np.arange(-1.0, 1.2, 0.2))

		plt.savefig(Path(str(self.OUT[l])+"bond_symbops_l"+str(l)+".png"), dpi=300, transparent=True)
		plt.savefig(Path(str(self.OUT[l])+"bond_symbops_l"+str(l)+".pdf"), dpi=300)

		plt.close()




		aspect = 1.0

		# bond symbop scatter plot l1 vs l2
		for l2 in [k for k in self.l if k<l]:

			fig, ax = plt.subplots(figsize=(12,12*aspect))

			ax.scatter([np.real(i[1]) for i in bond_symbops[l]], [np.real(i[1]) for i in bond_symbops[l2]], c="k", marker=".", s=2, alpha=0.3)

			plt.rc('text', usetex=True)
			plt.rc('font', family='serif')

			plt.xlabel(r"$\left| b^{*}_{ij} \right|_{"+str(l)+r"}$", fontsize=24)
			plt.ylabel(r"$\left| b^{*}_{ij} \right|_{"+str(l2)+r"}$", fontsize=24)

			ax.tick_params(labelsize=22)

			ax.set_xticks(np.arange(-1.0, 1.2, 0.2))
			ax.set_yticks(np.arange(-1.0, 1.2, 0.2))

			plt.savefig(Path(str(self.OUT[l])+"bond_symbops_l"+str(l)+"_bond_symbops_l"+str(l2)+"_scatter.png"), dpi=300, transparent=True)
			plt.savefig(Path(str(self.OUT[l])+"bond_symbops_l"+str(l)+"_bond_symbops_l"+str(l2)+"_scatter.pdf"), dpi=300)

			plt.close()
		






