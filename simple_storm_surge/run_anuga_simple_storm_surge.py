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
import numpy


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
    # Preparation of topographic data, extract from zip file
    #------------------------------------------------------------------------------
    # Unzip asc from zip file
    import zipfile as zf
    if project.verbose: print ('Reading ASC from cairns.zip')
    zf.ZipFile(anuga.join(project.data_dir,project.name_stem+'.zip')).extract(project.name_stem+'.asc')

    #------------------------------------------------------------------------------
    # Create a simple rectangular domain
    #------------------------------------------------------------------------------
    domain = anuga.rectangular_cross_domain(m=100, n=400, len1=400_000, len2=200_000)

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

        z = -100.0*numpy.ones_like(x)

        channel = np.logical_and(y>5,y<15)

        z = np.where(np.logical_and(channel,x<10), x/300, z)
        z = np.where(np.logical_and(channel,x>20), x/300, z)

        z = numpy.where()


    domain.set_quantity('elevation',
                        filename=project.name_stem + '.asc',
                        use_cache=project.cache,
                        verbose=project.verbose,
                        )



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

        #print('minx',np.min(x))
        #print('maxx',np.max(x))
        #print('miny',np.min(y))
        #print('maxy',np.max(y))

        # pressure cell in atms
        p = maxp - (maxp-minp)*np.exp( -((x-x0)**2 + (y-y0)**2)/r**2 ) 

        # convert to pascals 1 atm = 101325 pascals
        p = p*101325

        # print('x', x)
        # print('y', y)
        # print('p', p)
        # print('len(p)', len(p))

        p = p.reshape((1,-1))

        #print('p shape', p.shape)

        return p

    
    
    
    from anuga.operators.barometric_pressure import Barometric_pressure_operator
    Barometric_pressure_operator(domain, pressure_cell, use_coordinates=True)
    
    #from anuga.shallow_water.forcing import Barometric_pressure
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
        s = 300*np.ones_like(x)

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

Bd = anuga.Dirichlet_boundary([project.tide, 0, 0]) # Mean water level
Bs = anuga.Transmissive_stage_zero_momentum_boundary(domain) # Neutral boundary

# Set to tide
Bw = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
                    domain=domain, 
                    function=lambda t: [project.tide, 0, 0])

domain.set_boundary({'ocean_east': Bw,
                        'bottom': Bw,
                        'onshore': Bw,
                        'top': Bw})


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
