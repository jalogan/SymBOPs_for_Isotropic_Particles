from pathlib import Path
from main import main
from makeSettingsFile import makeSettingsFile


"""
This module is where you specify the input path of
the Settings.txt file, OR manually specify each 
parameter below without using a settings.txt input file.

If Settings.txt has been made and is located in the 
directory settingsIN defined here, then it will be used.
Remove any Settings.txt file from the settingsIN directory
(or rename it) to manually put the parameters in below.

If Settings.txt has not been made, the variable settingsIN
will be ignored, and only the manually input settings will
be used. Do not comment out settingsIN.
"""

###############################################
## Path to Settings.txt file (if being used) ##
###############################################
settingsIN = ""



#-------------------------------------------------------
#NOTE
## Only need to fill out parameters below this line if 
## you do not want to use a "Settings.txt" file to define the parameters


#############################
## Manually input Settings ##
#############################

# Paths to where to find data and save output files
IN = ""
OUT = ""
FILE = ""

# l-values to use for spherical harmonics (set as list)
l = [4,6]

# PBC box boundaries
# This can be left blank
boundaries = []

# Min and Max radius to look for neighbors
rmin = 67.9
rmax = 89.6

# Expected lattice
lattice_type = "fcc"

# parameters for choosing bonds in potential domains
half_angle = 0.5
qlm_ring_width = 0.04

# Minimum domain size to keep (default is 100)
dom_min = 35








#--------------------------------------------------------------------
# No need to change anything below this line.
#--------------------------------------------------------------------

if settingsIN!=None and Path.is_file(Path(settingsIN+"Settings.txt")):
	main(settingsIN)
else:
	makeSettingsFile(IN, FILE, OUT, l, half_angle, qlm_ring_width, boundaries, rmin, rmax, dom_min, lattice_type)
	main(IN)


print("\n\n--------------Thanks! Have a great day!--------------\n\n")



