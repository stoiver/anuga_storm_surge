"""Script for running a storm surge inundation scenario for Cairns, QLD Australia.

Source data such as elevation and boundary data is assumed to be available in
directories specified by project.py
The output sww file is stored in directory named after the scenario, i.e
pressure_cell or wind_shear or pressure_and_wind

The scenario is defined by a triangular mesh created from project.polygon,
and the elevation data.

Geoscience Australia, 2004-present

Modified:
Stephen Roberts, 2022-present
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
# Standard modules
import os
import time
import sys

# Related major packages
import anuga
import numpy as np


# Application specific imports
import project                 # Definition of file names and polygons

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

#------------------------------------------------------------------------------
# Do the domain creation on processor 0
#------------------------------------------------------------------------------
if anuga.myid == 0:

    #------------------------------------------------------------------------------
    # Create a simple rectangular domain
    #------------------------------------------------------------------------------
    domain = anuga.rectangular_cross_domain(m=2*project.rf, n=4*project.rf, len1=project.len1, len2=project.len2)

    # Print some stats about mesh and domain
    print ('Number of triangles = ', len(domain))
    print ('The extent is ', domain.get_extent())
    print (domain.geo_reference)
    print (domain.statistics())

    #------------------------------------------------------------------------------
    # Setup parameters of computational domain
    #------------------------------------------------------------------------------
    domain.set_name('simple_' + project.scenario) # Name of sww file
    domain.set_datadir('.')                       # Store sww output here
    domain.set_minimum_storable_height(0.01)      # Store only depth > 1cm
    domain.set_flow_algorithm(alg)


    #------------------------------------------------------------------------------
    # Setup initial conditions
    #------------------------------------------------------------------------------
    tide = project.tide
    domain.set_quantity('stage', project.tide)
    domain.set_quantity('friction', 0.0)

    #------------------------------------------------------------------------------
    # Setup elevation
    #------------------------------------------------------------------------------
    def topography(x,y):

        km = 1000.0
        continential_shelf = -10.0
        ocean = -1000.0

        z = (x-100*km)/(100*km)*ocean

        z = np.minimum(z, continential_shelf)
        z = np.maximum(z, ocean)

        return z

    domain.set_quantity('elevation', topography, verbose=project.verbose)
else:
    domain = None



#------------------------------------------------------------------------------
# Now produce parallel domain
#------------------------------------------------------------------------------
domain = anuga.distribute(domain,verbose=project.verbose)

domain.set_store_vertices_uniquely(False)

#------------------------------------------------------------------------------
# Setup information for pressure_cell scenario
#------------------------------------------------------------------------------
if project.scenario == 'pressure_cell' or project.scenario == 'pressure_and_wind':
    # Function for pressure cell
    
    def pressure_cell(t,x,y):
        import numpy as np

        x0, y0 = project.pressure_cell_center
        r = project.pressure_cell_radius
        maxp = project.pressure_cell_max
        minp = project.pressure_cell_min

        x = np.array(x)
        y = np.array(y)

          # pressure cell in atms
        p = maxp - (maxp-minp)*np.exp( -((x-x0)**2 + (y-y0)**2)/r**2 ) 

        # convert to pascals, 1 atm = 101325 pascals
        p = p*101325

        p = p.reshape((1,-1))

        return p

    def exact_stage(x,y):

        from anuga.config import rho_w
        from anuga import g
        maxp = project.pressure_cell_max*101325
        minp = project.pressure_cell_min*101325


        p = pressure_cell(0.0, x,y)

        p = p.flatten()

        print(x.shape, p.shape)

        dp = maxp - p

        print(np.max(dp))
        print(np.min(dp))

        # g h grad(w)  =  - 1/ rho_w  h grad(p)

        w = 1/(rho_w*g) * dp



        return w

    domain.set_quantity('stage', exact_stage)
    
    from anuga.operators.barometric_pressure import Barometric_pressure_operator
    Barometric_pressure_operator(domain, pressure_cell, use_coordinates=True)
    
    from anuga.shallow_water.forcing import Barometric_pressure
    #P = Barometric_pressure(pressure_cell, use_coordinates=True)
    #domain.forcing_terms.append(P)


#------------------------------------------------------------------------------
# Setup information for wind stress scenario
#------------------------------------------------------------------------------
if project.scenario == 'wind_stress' or project.scenario == 'pressure_and_wind':
    # Function for wind stress
    
    def wind_speed(t,x,y):
        import numpy as np

        x = np.array(x)
        y = np.array(y)

        # magnitude of wind in km/hr
        s = 100*np.ones_like(x)

        # change to m/s
        s = s * 1000/3600
        return s

    def wind_phi(t,x,y):
        import numpy as np

        x = np.array(x)
        y = np.array(y)

        # direction of wind in degrees
        phi = 180*np.ones_like(x)
        return phi

    from anuga.shallow_water.forcing import Wind_stress

    W = Wind_stress(wind_speed, wind_phi, use_coordinates=True)

    domain.forcing_terms.append(W)    

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
print ('Available boundary tags', domain.get_boundary_tags())

Br = anuga.Reflective_boundary(domain)

domain.set_boundary({'top': Br,
                     'bottom': Br,
                     'left': Br,
                     'right': Br})

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
import time
t0 = time.time()

for t in domain.evolve(yieldstep=100, finaltime=5000):
    if anuga.myid == 0:
        domain.print_timestepping_statistics()


domain.sww_merge(delete_old=True)

if anuga.myid == 0:
    print ('That took %.2f seconds' %(time.time()-t0))

    
anuga.finalize()
