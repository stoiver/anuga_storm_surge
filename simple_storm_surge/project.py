""" Common filenames and locations for topographic data, meshes and outputs.
    This file defines the parameters of the scenario you wish to run.
"""

import anuga


#------------------------------------------------------------------------------
# Runtime parameters
#------------------------------------------------------------------------------
cache = False
verbose = False

#------------------------------------------------------------------------------
# Define scenario as pressure_cell or wind_stress or both
#------------------------------------------------------------------------------
scenario = 'pressure_cell'  # Low pressure applied over sea
#scenario = 'wind_stress'    # On shore wind
#scenario = 'pressure_and_wind'  # Pressure_cell and On shore wind


#------------------------------------------------------------------------------
# Filenames
#------------------------------------------------------------------------------
name_stem = 'cairns'
meshname = name_stem + '.msh'

#------------------------------------------------------------------------------
# domain parameters
#------------------------------------------------------------------------------
rf = 20
len1 = 400_000
len2 = 200_000

#------------------------------------------------------------------------------
# Data for pressure cell
#------------------------------------------------------------------------------
pressure_cell_center = (150_000.0, 100_000.0)   # Assume to be on continental shelf
pressure_cell_radius = 20_000 # 20 kilometers
pressure_cell_min = 0.9 # of 1 atm
pressure_cell_max = 1.0 # 0f 1 atm

#------------------------------------------------------------------------------
# Data for Tides
#------------------------------------------------------------------------------
tide = 0.0
