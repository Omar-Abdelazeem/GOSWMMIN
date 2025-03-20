#%%
import pyswmm
import re
import pandas as pd
import numpy as np
import time
import pathlib
import os

# Location of input files
network = 'Modena'
directory = pathlib.Path('./Networks/'+network+'/Continuity and Computational Cost Files')
filepath = directory.resolve()
print(filepath)

# List of delta x max values (spatial discretization)
delta_xs=["10000m", "250m", "100m", "50m", "25m"]

continuities = []
times_exec = []
# Number of times to run the simulation to get an average execution time
timing_count = 5
# Loops over each line in the input file 
for delta_x in delta_xs:
    file= directory / ('Modena_12hr_SWMMIN_'+delta_x+'_VTS.inp')
    print(file.resolve())
    print(os.path.isfile(file.resolve()))
    print("\nDelta x=",delta_x)

    cont_error=0

    # Start timing
    start_time=time.perf_counter()
    # Run the simulation and get the continuity error for timing_count number of times
    for i in np.arange(timing_count):
        print("Timing count: ",i+1,"/",str(timing_count),"for  Delta x=",delta_x," with variable timesteps between 1 and 0.1 s ")
        with pyswmm.Simulation(str(file)) as sim:
            for step in sim:
                pass
            sim._model.swmm_end()
            cont_error=sim.flow_routing_error
    end_time=time.perf_counter()
    # Calculate the average execution time
    elapsed_time=(end_time-start_time)/timing_count

    continuities.append(cont_error)
    times_exec.append(elapsed_time)
print(zip(delta_xs,continuities,times_exec))
results_df = pd.DataFrame(zip(delta_xs,continuities,times_exec),columns=['Delta x','Continuity Error','Execution Time'])
results_df.to_csv(directory/'Variable Time Step Results.csv')