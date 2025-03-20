#%%
import pyswmm
import re
import pandas as pd
import numpy as np
import time
import pathlib
import multiprocessing
import os
from functools import partial

global timesteps

# Location of input files
dir = pathlib.Path('../Networks/Modena/Continuity and Computational Cost Files')
# Make a directory to store the temporary input files

if not os.path.exists(dir/'Temp'):
    os.makedirs(dir/'Temp')

# List of delta x max values (spatial discretization)
delta_xs=["10000m"]
# Log spaced array of solution speeds to sweep over
timesteps = np.logspace(-2,1.8,77)
timesteps[-1] = 60
print(timesteps)


def time_simulation(results_dict,delta_x, timestep):
    file=pathlib.Path('../Networks/Modena/Continuity and Computational Cost Files/Modena_12hr_SWMMIN_'+delta_x+'.inp')
    dir = file.parent
    # Calculate the time step based on the solution speed and spatial discretization
    print("\nDelta x=",delta_x," Delta t =",round(timestep,2))
    linecount=0
    file_read=open(file,'r')
    
    # Read the input file to find the line numbers for the routing and reporting steps
    for line in file_read:
        if re.search("ROUTING_STEP",line):
            step=linecount
        if re.search("REPORT_STEP",line):
            report=linecount
        linecount+=1

    file_read.close()
    file_read=open(file,'r')
    lines=file_read.readlines()
    file_read.close()

    # Update the routing and reporting steps in the input file
    # Reporting step is set to 1 hour to minimize storage requirements     
    lines[step]="ROUTING_STEP         "+str(timestep)+"\n"
    lines[report]="REPORT_STEP          01:00:00\n"

    # Save temp file
    file_upd=file.parent/pathlib.Path("Temp")/(file.stem+str(timestep)+".inp")
    file_o=open(file_upd,"w")
    file_o.writelines(lines)
    file_o.close()

    cont_error=0

    # Start timing
    start_time=time.perf_counter()
    # Run the simulation and get the continuity error for timing_count number of times
    if timestep <= 0.1:
        timing_count=1
    elif timestep <= 1: timing_count=1
    else: timing_count=1

    for i in np.arange(timing_count):
        print("Timing count: ",i,"/",str(timing_count),"for  Delta x=",delta_x," at time step Delta t =",round(timestep,2))
        with pyswmm.Simulation(str(file_upd)) as sim:
            for step in sim:
                pass
            sim._model.swmm_end()
            cont_error=sim.flow_routing_error
    end_time=time.perf_counter()
    # Calculate the average execution time
    elapsed_time=(end_time-start_time)/timing_count

    # Remove the temporary files to save space
    os.remove(file_upd)
    os.remove(file_upd.with_suffix(".rpt"))
    os.remove(file_upd.with_suffix(".out"))

    # Save the intermediate results to a dictionary and then to a csv file
    dict_temp = results_dict[delta_x]
    dict_temp[timestep]=(cont_error,elapsed_time)
    results_dict[delta_x]=dict_temp
    results_df = pd.DataFrame.from_dict(dict(results_dict))
    results_df.to_csv(dir/'Test_DF_TS.csv')

    return cont_error, elapsed_time


# Parallelize the simulation runs
if __name__=='__main__':
    # Assign the results_dict as manager.dict() to allow for shared memory
    with multiprocessing.Manager() as manager:
        results_dict = dict()
        for dx in delta_xs:
            results_dict[dx]={}
        # Create a list of all run combinations to parallelize
        run_combos = [(results_dict,delta_x, solspeed) for delta_x in delta_xs for solspeed in timesteps]
        # Run the simulations in parallel. Adjust the number of processes to the number of cores available
        with multiprocessing.Pool() as pool:
                pool=multiprocessing.Pool(processes=multiprocessing.cpu_count()-6)
                results = pool.starmap(time_simulation, run_combos)
    pool.close()
    pool.join()

    # Unpack the results and save them to csv files
    conts_df = pd.DataFrame( columns=delta_xs, index=timesteps)
    comps_df = pd.DataFrame( columns=delta_xs, index=timesteps)
    for i in range(len(run_combos)):
        conts_df.loc[run_combos[i][2],run_combos[i][1]]=results[i][0]
        comps_df.loc[run_combos[i][2],run_combos[i][1]]=results[i][1]
        
    comps_df.to_csv(dir/"Comp_times_Timesteps.csv")
    conts_df.to_csv(dir/"Continuity_Errors_Timesteps.csv")