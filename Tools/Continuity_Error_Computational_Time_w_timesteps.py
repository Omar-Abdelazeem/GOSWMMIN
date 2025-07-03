
import pyswmm
import re
import pandas as pd
import numpy as np
import time
import pathlib
import multiprocessing
import os
import swmmio

global timesteps

# Loops over each line in the input file 
def time_simulation(results_dict,delta_x, timestep, filepath, network):
    file = filepath/(network+"_12hr_SWMMIN_"+delta_x+".inp")
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
    lines[step]="ROUTING_STEP         "+str(round(timestep,4))+"\n"
    lines[report]="REPORT_STEP          01:00:00\n"

    # Save temp file
    file_upd=filepath/pathlib.Path("Temp")/(file.stem+str(round(timestep,5))+".inp")
    file_o=open(file_upd,"w")
    file_o.writelines(lines)
    file_o.close()

    cont_error=0

    # Start timing
    # Run the simulation and get the continuity error for timing_count number of times
    try:
        print("Simulating Delta x=",delta_x," at time step Delta t =",round(timestep,3))
        with pyswmm.Simulation(str(file_upd)) as sim:
            for step in sim:
                pass
            sim._model.swmm_end()
            cont_error=sim.flow_routing_error
    except: 
        cont_error = 9999

    # Remove the temporary files to save space
    os.remove(file_upd)
    os.remove(file_upd.with_suffix(".rpt"))
    os.remove(file_upd.with_suffix(".out"))

    # Save the intermediate results to a dictionary and then to a csv file
    dict_temp = results_dict[delta_x]
    dict_temp[timestep]=cont_error
    results_dict[delta_x]=dict_temp
    results_df = pd.DataFrame.from_dict(dict(results_dict))
    results_df.to_csv(filepath/'Test_DF_SS.csv')

    return cont_error



# Parallelize the simulation runs
if __name__=='__main__':
    # Location of input files
    network = 'Farina et al (2014)'
    directory = pathlib.Path('./Networks/'+network+'/Continuity and Computational Cost Files')
    filepath = directory.resolve()
    print(filepath)
    # Make a directory to store the temporary input files
    if not os.path.exists(directory/'Temp'):
        os.makedirs(directory/'Temp')

    # List of delta x max values (spatial discretization)
    delta_xs=["10000m", "250m", "100m", "50m", "25m"]
    # Log spaced array of timesteps to sweep over
    timesteps = np.logspace(np.log10(0.01), np.log10(60), 76)
    print(timesteps)

    # Assign the results_dict as manager.dict() to allow for shared memory
    with multiprocessing.Manager() as manager:
        results_dict = dict()
        for dx in delta_xs:
            results_dict[dx]={}
        # Create a list of all run combinations to parallelize
        run_combos = [(results_dict, delta_x, timestep, filepath, network)  for timestep in timesteps for delta_x in delta_xs]
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
        conts_df.loc[run_combos[i][2],run_combos[i][1]]=results[i]

    try:      
        conts_df.to_csv(directory/"Continuity_Errors_timestep.csv")
    except:
        conts_df.to_csv(directory/"Continuity_Errors_timestep_Failsafe.csv")

