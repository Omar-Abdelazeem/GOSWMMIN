import pyswmm
import re
import pandas as pd
import numpy as np
import time
import pathlib
import multiprocessing
import os

global timesteps

dir = pathlib.Path(r'C:\IWS_Modelling\Github\ISWMM_Development\Network-Files\Network 3\\')
if not os.path.exists(dir/'Temp'):
    os.makedirs(dir/'Temp')
delta_xs=["1000m","250m","100m","50m","25m"]
# sol_speeds = np.logspace(0,3,50)
sol_speeds = [1,1.132541315,1.282649831,1.452653926,1.645190588,1.863246312,2.110203429,2.389892566,2.70665207,3.065395295,3.471686819,3.931828756,4.45295851,
5.043159487,5.711586478,6.468607662,7.325965428,8.296958521,9.396648315,10.64209244,12.05260937,13.65007807,15.45927736,17.50827032,19.82883949,22.45697996,
25.43345761, 28.80444153, 32.6222201, 36.94601205, 41.84288508, 47.3887961, 53.66976946, 60.78323128, 68.8395207, 77.9636013, 88.29699955, 100,]

# sol_speeds = [
#     1, 1.097843767, 1.205260937, 1.323188207, 1.452653926, 1.594787058, 1.750827032, 
#     1.922134544, 2.110203429, 2.316673681, 2.543345761, 2.792196292, 3.065395295, 
#     3.365325118, 3.694601205, 4.056094905, 4.45295851, 4.888652745, 5.366976946, 
#     5.892102188, 6.468607662, 7.101520603, 7.79636013, 8.559185375, 9.396648315, 
#     10.31605178, 11.32541315, 12.43353424, 13.65007807, 14.98565312, 16.45190588, 
#     18.06162232]

# sol_speeds = [1,1.132541315,1.282649831,1.452653926,1.645190588,1.863246312,2.110203429,2.389892566,2.70665207,3.065395295,3.471686819,3.931828756,4.45295851,
# 5.043159487,5.711586478,6.468607662,7.325965428,8.296958521,9.396648315,10.64209244,12.05260937,13.65007807,15.45927736,17.50827032,19.82883949,22.45697996,
# 25.43345761, 28.80444153, 32.6222201, 36.94601205, 41.84288508, 47.3887961, 53.66976946, 60.78323128, 68.8395207, 77.9636013, 88.29699955, 100, 113.2541315,
# 128.2649831, 145.2653926, 164.5190588, 186.3246312, 211.0203429, 238.9892566, 270.665207, 306.5395295, 347.1686819, 393.1828756, 445.295851, 504.3159487,
# 571.1586478, 646.8607662, 732.5965428, 829.6958521, 939.6648315, 1064.209244]


timing_count = 5
# Loops over each line in the input file 
def time_simulation(results_dict,delta_x, solspeed):
    file=pathlib.Path(r'C:\IWS_Modelling\Github\ISWMM_Development\Network-Files\Network 3\Network3_12hr_Float-Storage_'+delta_x+'.inp')
    dir = file.parent
    timestep = float(delta_x[:-1]) / solspeed
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
    dict_temp[solspeed]=(cont_error,elapsed_time)
    results_dict[delta_x]=dict_temp
    results_df = pd.DataFrame.from_dict(dict(results_dict))
    results_df.to_csv(r'C:\IWS_Modelling\Github\ISWMM_Development\Network-Files\Network 3\Test_DF.csv')

    return cont_error, elapsed_time



if __name__=='__main__':
    with multiprocessing.Manager() as manager:
        results_dict = dict()
        for dx in delta_xs:
            results_dict[dx]={}
        run_combos = [(results_dict,delta_x, solspeed) for delta_x in delta_xs for solspeed in sol_speeds]
        with multiprocessing.Pool() as pool:
                pool=multiprocessing.Pool(processes=multiprocessing.cpu_count()-6)
                results = pool.starmap(time_simulation, run_combos)
    pool.close()
    pool.join()

    conts_df = pd.DataFrame( columns=delta_xs, index=sol_speeds)
    comps_df = pd.DataFrame( columns=delta_xs, index=sol_speeds)
    for i in range(len(run_combos)):
        conts_df.loc[run_combos[i][2],run_combos[i][1]]=results[i][0]
        comps_df.loc[run_combos[i][2],run_combos[i][1]]=results[i][1]
        
    comps_df.to_csv(r"C:\IWS_Modelling\Github\ISWMM_Development\Network-Files\Network 3\Comp_times_solspeed.csv")
    conts_df.to_csv(r"C:\IWS_Modelling\Github\ISWMM_Development\Network-Files\Network 3\Continuity_Errors_solspeed.csv")