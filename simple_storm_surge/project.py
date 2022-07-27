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

data_dir = 'data'

# Filename for locations where timeseries are to be produced
#gauge_filename = anuga.join(data_dir, 'gauges.csv')
gauge_filename = anuga.join(data_dir, 'low_cell.csv')

#------------------------------------------------------------------------------
# Domain definitions
#------------------------------------------------------------------------------
# bounding polygon for study area
bounding_polygon = anuga.read_polygon(anuga.join(data_dir, 'extent.csv'))

A = anuga.polygon_area(bounding_polygon) / 1000000.0
print ('Area of bounding polygon = %.2f km^2' % A)

#------------------------------------------------------------------------------
# Interior region definitions
#------------------------------------------------------------------------------
# Read interior polygons
poly_cairns = anuga.read_polygon(anuga.join(data_dir, 'cairns.csv'))
poly_island0 = anuga.read_polygon(anuga.join(data_dir, 'islands.csv'))
poly_island1 = anuga.read_polygon(anuga.join(data_dir, 'islands1.csv'))
poly_island2 = anuga.read_polygon(anuga.join(data_dir, 'islands2.csv'))
poly_island3 = anuga.read_polygon(anuga.join(data_dir, 'islands3.csv'))
poly_shallow = anuga.read_polygon(anuga.join(data_dir, 'shallow.csv'))

# Optionally plot points making up these polygons
#plot_polygons([bounding_polygon, poly_cairns, poly_island0, poly_island1,
#               poly_island2, poly_island3, poly_shallow],
#               style='boundingpoly', verbose=False)

# Define resolutions (max area per triangle) for each polygon
# Make these numbers larger to reduce the number of triangles in the model,
# and hence speed up the simulation

# bigger base_scale == less triangles
just_fitting = False
#base_scale = 25000 # 635763 triangles
#base_scale = 50000 # 321403 triangles
#base_scale = 100000 # 162170 triangles
base_scale = 400000 # 42093 triangles
base_scale = 1000000 # 17831 triangles
default_res = 100 * base_scale   # Background resolution as area of triangles
islands_res = base_scale
cairns_res = base_scale
shallow_res = 5 * base_scale

# Define list of interior regions with associated resolutions
interior_regions = [[poly_cairns,  cairns_res],
                    [poly_island0, islands_res],
                    [poly_island1, islands_res],
                    [poly_island2, islands_res],
                    [poly_island3, islands_res],
                    [poly_shallow, shallow_res]]

#------------------------------------------------------------------------------
# Data for exporting ascii grid
#------------------------------------------------------------------------------
eastingmin = 363000
eastingmax = 418000
northingmin = 8026600
northingmax = 8145700
gauge_filename = anuga.join(data_dir, 'low_cell.csv')
#------------------------------------------------------------------------------
# Data for pressure cell
#------------------------------------------------------------------------------
#pressure_cell_center = (451_871, 8_128_376)   # Assume to be on continental shelf
easting=301278.510000
northing=7961069.140000
pressure_cell_center = (200_000.0, 100_000.0)   # Assume to be on continental shelf

pressure_cell_radius = 20_000 # 20 kilometers
pressure_cell_min = 0.9 # of 1 atm
pressure_cell_max = 1.0 # 0f 1 atm

#------------------------------------------------------------------------------
# Data for Tides
#------------------------------------------------------------------------------
tide = 0.0
