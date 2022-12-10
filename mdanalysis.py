# Project Title: "Automated computational workflow to prioritize potential resistance variants identified in HIV Integrase Subtype C and CRF02_AG" 
# This script is developed for the fufuillment for Masters at the South African National Bioinformatics Institute at the University of the Western Cape
# The project is funded by the Poliomyelitis Research Foundation and the UWC Ada & Bertie Levenstein Bursary Programme
# Currently any licensing and usage of this software is governed under the regulations of the afore mentioned parties

#Author:	Keaghan Brown (3687524) - MSc Bioinformatics Candidate (3687524@myuwc.ac.za)
#Author:	Ruben Cloete (Supervisor) - Lecturer at South African National Bioinformatics Institute (ruben@sanbi.ac.za)

# This script completes the final objective for The Automated Computational Workflow to Prioritize Potential Resistance Variants Identified in HIV Integrase Subtype C and CRF02_AG

import MDAnalysis as mda
from MDAnalysis.analysis import rms
import MDAnalysis.transformations as trans
from MDAnalysis.tests.datafiles import PDB, GRO, XTC
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import datetime
import os

class TrajAn:

	def __init__(self):
		return
	
	def traj_repair(self, systems_folder):
		print(datetime.datetime.now())
		self.systems_folder = systems_folder
		for system in os.listdir(self.systems_folder):
			traj_file = ''
			tpr_file = ''
			os.chdir(self.systems_folder + '/' + system)
			print(system)
			for file in os.listdir():
				if file.endswith('.xtc') and not 'step' in file and system.lower() in file:
					traj_file = str(file)
					print(traj_file)
				elif file.endswith('.tpr') and not 'step' in file and system.lower() in file:
					tpr_file = str(file)
					print(tpr_file)
			subprocess.call(['/home/akuma/Desktop/analysis.sh', traj_file, tpr_file])
				
			
p = TrajAn()
p.traj_repair('/home/akuma/Desktop/MD_results')
"""			
subprocess.call(['/home/akuma/Desktop/analysis.sh', '/home/akuma/Desktop/MD_results/S39C/md_p1.xtc', '/home/akuma/Desktop/MD_results/S39C/md_p1.tpr'])

print(datetime.datetime.now())

u = mda.Universe('/home/akuma/Desktop/MD_results/S39C/md_p1.gro', '/home/akuma/Desktop/MD_results/S39C/md_p1_nojump.xtc')

R = rms.RMSD(u,  # universe to align
             u,  # reference universe or atomgroup
             select='backbone',  # group to superimpose and calculate RMSD
             ref_frame=0)  # frame index of the reference
             
R.run()

print(datetime.datetime.now())

df = pd.DataFrame(R.rmsd,
                  columns=['Frame', 'Time (ps)',
                           'Backbone'])
                           
df.plot(x="Time (ps)", y=["Backbone"], linewidth='0.25')

print(datetime.datetime.now())
plt.margins(0)
plt.show()
"""
