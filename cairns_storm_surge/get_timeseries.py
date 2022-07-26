"""
Generate time series of nominated "gauges" read from project.gauge_filename. This 
is done by first running sww2csv_gauges on two different directories to make 
'csv' files. Then running csv2timeseries_graphs detailing the two directories 
containing the csv file and produces one set of graphs in the 'output_dir' containing
the details at the gauges for both these sww files.

Note, this script will only work if pylab is installed on the platform
"""

from os import sep
import os
import shutil
import glob

import project
import anuga

verbose = True

possible_scenerios = ['cairns_pressure_and_wind', 'cairns_wind_shear', 'cairns_pressure_cell']

for scenerio in possible_scenerios:

    process_scenerio = True

    try:
        anuga.sww2csv_gauges(scenerio+'.sww',
                    project.gauge_filename,
                    quantities=['stage','speed','depth','elevation'],
                    use_cache=False,
                    verbose=verbose)

    except:
        print ('Failed to process %s' % scenerio )
        process_scenerio = False

    if process_scenerio:
        # check if directory exists or not yet
        if not os.path.exists(scenerio):
            os.makedirs(scenerio)
        
        # move files into created directory
        pattern = 'gauge_*'

        files = glob.glob(pattern)

        for file in files:
            print('Moving ',file,' to directory ',scenerio)
            shutil.copy(file, scenerio)
            os.remove(file)

        import pylab
        anuga.csv2timeseries_graphs(directories_dic={scenerio+sep: [scenerio,0,0]},
                            output_dir=scenerio+sep,
                            base_name='gauge_',
                            plot_numbers='',
                            quantities=['stage','speed','depth'],
                            extra_plot_name='',
                            assess_all_csv_files=True,                            
                            create_latex=False,
                            verbose=verbose)
    # except ImportError:
    #     #ANUGA does not rely on pylab to work 
    #     print ('must have pylab installed to generate plots')
    # except FileNotFoundError:
    #     print ('data for scenerio '+ scenerio + ' not found')

    

