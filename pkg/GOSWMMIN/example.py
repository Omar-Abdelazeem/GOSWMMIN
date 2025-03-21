#%%
from SWMMIN_sim import SWMMIN_sim
import pathlib

# %%
filepath = pathlib.Path('../Resources/Modena_12hr.inp')
sim = SWMMIN_sim(filepath)
min_pressure = 0
des_pressure = 10
pdw_exponent = 0.5

# Create the SWMMIN simulation
# Creates an 8-hr supply duration with a non-adaptive
# sim.Convert_to_SWMMIN(supply_duration= 8.0, 
#                   minimum_pressure=min_pressure, desired_pressure = des_pressure, pdw_exponent=pdw_exponent,
#                   n_days= 1,length_to_diameter=30, solution_speed=100, q_des=q_des, tank_heights=tank_heights, 
#                   tank_areas=tank_areas, consum_pattern=consum_pattern)

sim.Convert_to_SWMMIN(supply_duration= 8.0, 
                  minimum_pressure=min_pressure, desired_pressure = des_pressure, pdw_exponent=pdw_exponent,
                  n_days= 1,maximum_xdelta=100, timestep=0.5)

#%%
# Running the simulation
sim.Run_SWMMIN()

#%% 
sim.track_water_balances()
sim.mass_balance()
sim.effective_supply_duration()
#%%
# Creating a SWMMIN simulation object with optional CSV inputs for PDW variables
# tank sizes and consumption patterns
filepath = pathlib.Path('../Resources/Linear_Network.inp')
sim = SWMMIN_sim(filepath)
min_pressure = '../Resources/min_pressure.csv'
des_pressure = '../Resources/des_pressure.csv'
pdw_exponent = '../Resources/pdw_exponent.csv'
q_des = '../Resources/q_des.csv'
tank_areas = '../Resources/tank_areas.csv'
tank_heights = '../Resources/tank_heights.csv'
consum_pattern = '../Resources/consum_pattern.csv'

# Create the SWMMIN simulation
# Creates an 8-hr supply duration with a non-adaptive
# sim.Convert_to_SWMMIN(supply_duration= 8.0, 
#                   minimum_pressure=min_pressure, desired_pressure = des_pressure, pdw_exponent=pdw_exponent,
#                   n_days= 1,length_to_diameter=30, solution_speed=100, q_des=q_des, tank_heights=tank_heights, 
#                   tank_areas=tank_areas, consum_pattern=consum_pattern)

sim.Convert_to_SWMMIN(supply_duration= 8.0, 
                  minimum_pressure=min_pressure, desired_pressure = des_pressure, pdw_exponent=pdw_exponent,
                  n_days= 1,maximum_xdelta=100, timestep=0.5, q_des=q_des, tank_heights=tank_heights, 
                  tank_areas=tank_areas, consum_pattern=consum_pattern)

#%%
# Running the simulation
sim.Run_SWMMIN()

#%%
# Getting Pressure Results
pressures = sim.get_pressures(specific_nodes=['DN1','DN2'])
tank_vols, tank_heights = sim.get_tank_vols_heights(specific_nodes=['DN1','DN2'])
