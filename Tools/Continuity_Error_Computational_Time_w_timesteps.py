import pyswmm
import re
from swmm.toolkit.shared_enum import LinkAttribute,NodeAttribute
import pandas as pd
import numpy as np
import datetime
import matplotlib as mpl
import matplotlib.pyplot as plt
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
import time
import pathlib
import multiprocessing
import sys,os
import tempfile
import random
from functools import partial
import json

global timesteps

dir = pathlib.Path(r'C:\IWS_Modelling\Github\ISWMM_Development\Network-Files\Network 3\\')
if not os.path.exists(dir/'Temp'):
    os.makedirs(dir/'Temp')
delta_xs=["1000m", "250m", "100m","50m","25m"]
timesteps=np.arange(10,60.1,1)
timesteps=np.arange(0.01,0.101,0.01)
timesteps=np.append(timesteps,np.arange(5,10.01,0.25))
timesteps=np.append(timesteps,np.arange(0.1,1.01,0.05))

results_dict = dict()
for dx in delta_xs:
    results_dict[dx]={}
run_combos = [(results_dict,delta_x, timestep) for delta_x in delta_xs for timestep in timesteps]

timing_count = 1
# Loops over each line in the input file 
def time_simulation(results_dict,delta_x, timestep):
    file=pathlib.Path(r'C:\IWS_Modelling\Github\ISWMM_Development\Network-Files\Network 3\Network3_12hr_Float-Storage_'+delta_x+'.inp')
    dir = file.parent
    print("\nDelta x=",delta_x," Delta t =",timestep)
    linecount=0
    file_read=open(file,'r')
    
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

    lines[step]="ROUTING_STEP         "+str(timestep)+"\n"
    lines[report]="REPORT_STEP          01:00:00\n"
    file_upd=file.parent/pathlib.Path("Temp")/(file.stem+str(timestep)+".inp")
    file_o=open(file_upd,"w")
    file_o.writelines(lines)
    file_o.close()
    cont_error=0
    start_time=time.perf_counter()
    for i in np.arange(timing_count):
        print("Timing count: ",i,"/",str(timing_count),"for  Delta x=",delta_x," at time step Delta t =",timestep)
        with pyswmm.Simulation(str(file_upd)) as sim:
            for step in sim:
                pass
            sim._model.swmm_end()
            cont_error=sim.flow_routing_error
    end_time=time.perf_counter()
    elapsed_time=(end_time-start_time)/timing_count
    os.remove(file_upd)
    os.remove(file_upd.with_suffix(".rpt"))
    os.remove(file_upd.with_suffix(".out"))
    dict_temp = results_dict[delta_x]
    dict_temp[timestep]=(cont_error,elapsed_time)
    results_dict[delta_x]=dict_temp
    results_df = pd.DataFrame.from_dict(dict(results_dict))
    results_df.to_csv(r'C:\IWS_Modelling\Github\ISWMM_Development\Network-Files\Network 3\Test_DF.csv')

    return cont_error, elapsed_time



if __name__=='__main__':
    with multiprocessing.Manager() as manager:
        
        # print(results_dict.keys(),results_dict.values())
        with multiprocessing.Pool() as pool:
    
                pool=multiprocessing.Pool(processes=multiprocessing.cpu_count()-6)
                results = pool.starmap(time_simulation, run_combos)
    pool.close()
    pool.join()

    conts_df = pd.DataFrame( columns=delta_xs, index=timesteps)
    comps_df = pd.DataFrame( columns=delta_xs, index=timesteps)
    for i in range(len(run_combos)):
        conts_df.loc[run_combos[i][2],run_combos[i][1]]=results[i][0]
        comps_df.loc[run_combos[i][2],run_combos[i][1]]=results[i][1]
        
    comps_df.to_csv(r"C:\IWS_Modelling\Github\ISWMM_Development\Network-Files\Network 3\Comp_times_2.csv")
    conts_df.to_csv(r"C:\IWS_Modelling\Github\ISWMM_Development\Network-Files\Network 3\Continuity_Errors_2.csv")