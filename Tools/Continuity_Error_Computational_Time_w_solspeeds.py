import pyswmm
import re
import pandas as pd
import numpy as np
import time
import pathlib
import multiprocessing
import os
import swmmio


# Loops over each line in the input file 
def time_simulation(results_dict,delta_x, solspeed, mean_xs, timing_count, filepath, network):
    file = filepath/(network+"_12hr_SWMMIN_"+delta_x+".inp")
    # Calculate the time step based on the solution speed and spatial discretization
    timestep = mean_xs[delta_x] / solspeed
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
    file_upd=filepath/pathlib.Path("Temp")/(file.stem+str(round(solspeed,5))+".inp")
    file_o=open(file_upd,"w")
    file_o.writelines(lines)
    file_o.close()

    cont_error=0

    # Start timing
    start_time=time.perf_counter()
    # Run the simulation and get the continuity error for timing_count number of times
    try:
        for i in np.arange(timing_count):
            print("Timing count: ",i+1,"/",str(timing_count),"for  Delta x=",delta_x," at time step Delta t =",round(timestep,3))
            with pyswmm.Simulation(str(file_upd)) as sim:
                for step in sim:
                    pass
                sim._model.swmm_end()
                cont_error=sim.flow_routing_error
    except: 
        cont_error = 9999
    end_time=time.perf_counter()
    # Calculate the average execution time
    elapsed_time=(end_time-start_time)/timing_count

    # Remove the temporary files to save space
    os.remove(file_upd)
    os.remove(file_upd.with_suffix(".rpt"))
    os.remove(file_upd.with_suffix(".out"))

    # Save the intermediate results to a dictionary and then to a csv file
    dict_temp = results_dict[delta_x]
    dict_temp[solspeed]=(cont_error,elapsed_time)
    results_dict[delta_x]=dict_temp
    results_df = pd.DataFrame.from_dict(dict(results_dict))
    results_df.to_csv(filepath/'Test_DF_SS.csv')

    return cont_error, elapsed_time



# Parallelize the simulation runs
if __name__=='__main__':
    # Location of input files
    network = 'Pescara'
    directory = pathlib.Path('./Networks/'+network+'/Continuity and Computational Cost Files')
    filepath = directory.resolve()
    print(filepath)
    # Make a directory to store the temporary input files
    if not os.path.exists(directory/'Temp'):
        os.makedirs(directory/'Temp')

    # List of delta x max values (spatial discretization)
    delta_xs=["10000m", "250m", "100m", "50m", "25m"]
    # Log spaced array of solution speeds to sweep over
    sol_speeds = np.logspace(3, 0, 75)
    # sol_speeds = np.append(sol_speeds, np.logspace(3,0,75))
    # sol_speeds = np.unique(sol_speeds)
    # sol_speeds = sol_speeds[::-1]
    print(sol_speeds)

    mean_xs = {}

    # Find mean conduit length for each delta x
    for x in delta_xs:
        filename = network+ "_12hr_SWMMIN_"+x+".inp"
        # Read the SWMM input file and get conduit data
        model = swmmio.Model(str(filepath / filename))
        links = model.links()
        # Calculate mean length of conduits
        mean_xs[x]=links['Length'].mean()

        # Export mean_xs to a csv file
        mean_xs_df = pd.DataFrame(list(mean_xs.items()), columns=['Delta_x', 'Mean_Conduit_Length'])
        # mean_xs_df.to_csv(directory/'mean_delta_xs.csv', index=False)

    timing_count=1

    # Assign the results_dict as manager.dict() to allow for shared memory
    with multiprocessing.Manager() as manager:
        results_dict = dict()
        for dx in delta_xs:
            results_dict[dx]={}
        # Create a list of all run combinations to parallelize
        run_combos = [(results_dict,delta_x, solspeed, mean_xs, timing_count,filepath, network)  for solspeed in sol_speeds for delta_x in delta_xs]
        # Run the simulations in parallel. Adjust the number of processes to the number of cores available
        with multiprocessing.Pool() as pool:
                pool=multiprocessing.Pool(processes=multiprocessing.cpu_count()-6)
                results = pool.starmap(time_simulation, run_combos)
    pool.close()
    pool.join()

    # Unpack the results and save them to csv files
    conts_df = pd.DataFrame( columns=delta_xs, index=sol_speeds)
    comps_df = pd.DataFrame( columns=delta_xs, index=sol_speeds)
    for i in range(len(run_combos)):
        conts_df.loc[run_combos[i][2],run_combos[i][1]]=results[i][0]
        comps_df.loc[run_combos[i][2],run_combos[i][1]]=results[i][1]

    try:      
        comps_df.to_csv(directory/"Comp_times_solspeed.csv")
        conts_df.to_csv(directory/"Continuity_Errors_solspeed.csv")
    except:
        comps_df.to_csv(directory/"Comp_times_solspeed_Failsafe.csv")
        conts_df.to_csv(directory/"Continuity_Errors_solspeed_Failsafe.csv")